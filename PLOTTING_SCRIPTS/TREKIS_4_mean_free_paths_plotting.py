#!/usr/bin/env python3
"""
Standalone script to batch parse TREKIS mean free path (IMFP/MFP) files with varied naming conventions, 
group them by particle type, prioritize total files, exclude data > 1e20, and generate summary plots.
Scans folders matching 'MFPs_and_Ranges_in_*' if present.

Usage:
    python TREKIS_5_mfp_plotting.py
    python TREKIS_5_mfp_plotting.py --dir /path/to/data
    python TREKIS_5_mfp_plotting.py -o ./plots --show
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np

# Updated legend decoding mapping as specified
PROCESS_DECODE = {
    "imfps": "inelastic",
    "emfps": "elastic",
    "brems": "bremsstrahlung",
    "brems_mfps": "bremsstrahlung",
    "pair": "pair creation",
}

KNOWN_PARTICLES = {"electron", "positron", "photon", "proton", "ion", "neutron", "muon", "hole"}
MODEL_KEYWORDS = {"cdf", "no", "bhw", "penelope", "dos", "delta", "sp", "total"}


def parse_filename(filename):
    """
    Parses varied filenames matching OUTPUT_[process]_[material]_[particle]_[model].dat
    Returns a dictionary with process, material, particle, model, and total flag.
    """
    base = os.path.splitext(filename)[0]
    if not base.startswith("OUTPUT_"):
        return None
    
    remainder = base[len("OUTPUT_"):]
    tokens = remainder.split("_")
    # Filter out empty tokens from double underscores
    tokens = [t for t in tokens if t]
    
    if len(tokens) < 3:
        return None

    # Identify particle
    particle = None
    particle_idx = -1
    for i, token in enumerate(tokens):
        if token.lower() in KNOWN_PARTICLES:
            particle = token.lower()
            particle_idx = i
            break
            
    if particle_idx == -1:
        # Fallback: find token right before a known model keyword
        for i, token in enumerate(tokens):
            if token.lower() in MODEL_KEYWORDS and i > 0:
                particle = tokens[i-1].lower()
                particle_idx = i - 1
                break
                
    if particle_idx == -1:
        particle = "unknown"
        particle_idx = len(tokens) - 2 if len(tokens) >= 2 else 0

    # Determine process name (handles multi-word processes like Annihilation_MFPs or Brems_MFPs)
    if len(tokens) > 1 and tokens[1].lower() == "mfps":
        process_name = "_".join(tokens[:2])
        material_start = 2
    else:
        process_name = tokens[0]
        material_start = 1

    material_element = "_".join(tokens[material_start:particle_idx])
    model = "_".join(tokens[particle_idx + 1:])
    is_total = "total" in model.lower()
    
    return {
        "filename": filename,
        "process": process_name,
        "material": material_element,
        "particle": particle,
        "model": model,
        "is_total": is_total
    }


def read_mfp_file(filepath):
    """Reads energy (col 0) and MFP value (last column) from a .dat file, excluding values > 1e20."""
    energies = []
    values = []
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                val_row = [float(p) for p in parts]
                if len(val_row) >= 2:
                    val = val_row[-1]
                    # Exclude data > 1e20
                    if val <= 1e20:
                        energies.append(val_row[0])
                        values.append(val)
            except ValueError:
                continue
                
    return np.array(energies), np.array(values)


def process_directory(source_dir, output_dir, show_plots):
    """Processes MFP files inside a specific directory and saves summary plots."""
    all_files = [f for f in os.listdir(source_dir) if f.lower().endswith(".dat") and f.startswith("OUTPUT_")]
    
    parsed_files = []
    for f in all_files:
        info = parse_filename(f)
        if info:
            # Exclude stopping process as requested
            if "stopping" in info["process"].lower():
                continue
            parsed_files.append(info)

    if not parsed_files:
        print(f"No valid MFP data files found in: {source_dir}")
        return

    particles = set(item["particle"] for item in parsed_files)

    for particle in particles:
        particle_files = [item for item in parsed_files if item["particle"] == particle]
        
        process_groups = {}
        for item in particle_files:
            proc = item["process"].lower()
            if proc not in process_groups:
                process_groups[proc] = []
            process_groups[proc].append(item)

        selected_files = []
        for proc, group in process_groups.items():
            total_matches = [g for g in group if g["is_total"]]
            if total_matches:
                selected_files.append(total_matches[0])
            else:
                selected_files.append(group[0])

        # Generate plot for this particle with figsize=(5, 4) and no grid
        fig, ax = plt.subplots(figsize=(5, 4))

        material_name = selected_files[0]["material"] if selected_files else ""
        
        for item in selected_files:
            filepath = os.path.join(source_dir, item["filename"])
            energies, values = read_mfp_file(filepath)
            
            if len(energies) == 0:
                continue

            raw_proc = item["process"]
            label_name = PROCESS_DECODE.get(raw_proc.lower(), raw_proc)

            ax.plot(energies, values, label=label_name, lw=1.5)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Energy (eV)", fontsize=10)
        ax.set_ylabel("Mean free path (A)", fontsize=10)
        ax.set_title(f"Mean Free Paths for {particle.capitalize()} ({material_name})", fontsize=10, pad=10)
        
        ax.grid(False)
        ax.legend(fontsize=8, frameon=True, framealpha=0.85)
        
        plt.tight_layout()

        save_filename = f"Mean_free_paths_{particle}_{material_name}.png"
        out_path = os.path.join(output_dir, save_filename)
        plt.savefig(out_path, dpi=300)
        print(f"  -> Saved MFP plot: {out_path}")

        if show_plots:
            plt.show()
        else:
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Batch parse and plot Mean Free Path (IMFP) data grouped by particle type."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Directory to scan for data files",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default=None,
        help="Directory to save generated plots",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display interactive plot windows",
    )

    args = parser.parse_args()

    target_dir = os.path.abspath(args.dir)
    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    matching_subdirs = [
        d for d in glob.glob(os.path.join(target_dir, "MFPs_and_Ranges_in_*")) 
        if os.path.isdir(d)
    ]

    if matching_subdirs:
        print(f"Found {len(matching_subdirs)} matching directory(ies) matching 'MFPs_and_Ranges_in_*'.")
        for sub_d in matching_subdirs:
            print(f"\nProcessing directory: {sub_d}")
            out_d = os.path.abspath(args.output_dir) if args.output_dir else sub_d
            os.makedirs(out_d, exist_ok=True)
            process_directory(sub_d, out_d, args.show)
    else:
        print(f"No directories matching 'MFPs_and_Ranges_in_*' found inside {target_dir}. Scanning target directory directly.")
        out_dir = os.path.abspath(args.output_dir) if args.output_dir else target_dir
        os.makedirs(out_dir, exist_ok=True)
        process_directory(target_dir, out_dir, args.show)

    print("\nMFP plotting complete.")


if __name__ == "__main__":
    main()