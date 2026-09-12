#!/usr/bin/env python3
"""
Standalone script to batch parse TREKIS stopping files, group them by particle type, 
prioritize total files, and generate 3 summary plots per particle:
1. Stopping power vs Energy (Col 1 vs Col 0)
2. Range vs Energy (Col 2 vs Col 0)
3. Stopping power vs Range (Col 1 vs Col 2)
Scans folders matching 'MFPs_and_Ranges_in_*' if present.

Usage:
    python TREKIS_6_stopping_plotting.py
    python TREKIS_6_stopping_plotting.py --dir /path/to/data
    python TREKIS_6_stopping_plotting.py -o ./plots --show
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np

KNOWN_PARTICLES = {"electron", "positron", "photon", "proton", "ion", "neutron", "muon", "hole"}
MODEL_KEYWORDS = {"cdf", "no", "bhw", "penelope", "dos", "delta", "sp", "total"}


def parse_filename(filename):
    """
    Parses varied filenames matching OUTPUT_[process]_[material]_[particle]_[model].dat
    Specifically looking for stopping-related processes.
    """
    base = os.path.splitext(filename)[0]
    if not base.startswith("OUTPUT_"):
        return None
    
    remainder = base[len("OUTPUT_"):]
    tokens = remainder.split("_")
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
        for i, token in enumerate(tokens):
            if token.lower() in MODEL_KEYWORDS and i > 0:
                particle = tokens[i-1].lower()
                particle_idx = i - 1
                break
                
    if particle_idx == -1:
        particle = "unknown"
        particle_idx = len(tokens) - 2 if len(tokens) >= 2 else 0

    process_name = tokens[0]
    material_element = "_".join(tokens[1:particle_idx])
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


def read_stopping_file(filepath):
    """
    Reads Energy (col 0), Stopping Power (col 1), and Range (col 2) from a .dat file.
    Excludes data > 1e20 or invalid entries.
    """
    energies = []
    stopping_powers = []
    ranges = []
    
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                val_row = [float(p) for p in parts]
                if len(val_row) >= 3:
                    e = val_row[0]
                    sp = val_row[1]
                    rng = val_row[2]
                    
                    if e <= 1e20 and sp <= 1e20 and rng <= 1e20:
                        energies.append(e)
                        stopping_powers.append(sp)
                        ranges.append(rng)
            except ValueError:
                continue
                
    return np.array(energies), np.array(stopping_powers), np.array(ranges)


def process_directory(source_dir, output_dir, show_plots):
    """Processes stopping files inside a specific directory and saves the 3 required plots per particle."""
    all_files = [f for f in os.listdir(source_dir) if f.lower().endswith(".dat") and f.startswith("OUTPUT_")]
    
    parsed_files = []
    for f in all_files:
        info = parse_filename(f)
        if info:
            # Filter strictly for stopping files
            if "stop" in info["process"].lower():
                parsed_files.append(info)

    if not parsed_files:
        print(f"No valid stopping data files found in: {source_dir}")
        return

    particles = set(item["particle"] for item in parsed_files)

    for particle in particles:
        particle_files = [item for item in parsed_files if item["particle"] == particle]
        
        # Prioritize 'total' file if available
        total_matches = [g for g in particle_files if g["is_total"]]
        selected_file = total_matches[0] if total_matches else particle_files[0]

        filepath = os.path.join(source_dir, selected_file["filename"])
        energies, stopping_powers, ranges = read_stopping_file(filepath)

        if len(energies) == 0:
            continue

        material_name = selected_file["material"]

        # -------------------------------------------------------------
        # Plot 1: Stopping Power vs Energy (Col 1 vs Col 0)
        # -------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(energies, stopping_powers, lw=1.5, color="tab:blue")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Energy (eV)", fontsize=10)
        ax.set_ylabel("Stopping Power (eV/A)", fontsize=10)
        ax.set_title(f"Stopping Power vs Energy — {particle.capitalize()} ({material_name})", fontsize=10, pad=10)
        ax.grid(False)
        plt.tight_layout()

        save_path_1 = os.path.join(output_dir, f"Stopping_power_vs_energy_{particle}_{material_name}.png")
        plt.savefig(save_path_1, dpi=300)
        print(f"  -> Saved plot: {save_path_1}")
        if show_plots:
            plt.show()
        plt.close(fig)

        # -------------------------------------------------------------
        # Plot 2: Range vs Energy (Col 2 vs Col 0)
        # -------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(energies, ranges, lw=1.5, color="tab:orange")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Energy (eV)", fontsize=10)
        ax.set_ylabel("Range (A)", fontsize=10)
        ax.set_title(f"Range vs Energy — {particle.capitalize()} ({material_name})", fontsize=10, pad=10)
        ax.grid(False)
        plt.tight_layout()

        save_path_2 = os.path.join(output_dir, f"Range_vs_energy_{particle}_{material_name}.png")
        plt.savefig(save_path_2, dpi=300)
        print(f"  -> Saved plot: {save_path_2}")
        if show_plots:
            plt.show()
        plt.close(fig)

        # -------------------------------------------------------------
        # Plot 3: Stopping Power vs Range (Col 1 vs Col 2)
        # -------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(ranges, stopping_powers, lw=1.5, color="tab:green")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Range (A)", fontsize=10)
        ax.set_ylabel("Stopping Power (eV/A)", fontsize=10)
        ax.set_title(f"Stopping Power vs Range — {particle.capitalize()} ({material_name})", fontsize=10, pad=10)
        ax.grid(False)
        plt.tight_layout()

        save_path_3 = os.path.join(output_dir, f"Stopping_power_vs_range_{particle}_{material_name}.png")
        plt.savefig(save_path_3, dpi=300)
        print(f"  -> Saved plot: {save_path_3}")
        if show_plots:
            plt.show()
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Batch parse and plot stopping power and range data grouped by particle type."
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
        out_d = os.path.abspath(args.output_dir) if args.output_dir else target_dir
        os.makedirs(out_d, exist_ok=True)
        process_directory(target_dir, out_d, args.show)

    print("\nStopping power & range plotting complete.")


if __name__ == "__main__":
    main()