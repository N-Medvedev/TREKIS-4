#!/usr/bin/env python3
"""
Standalone script to batch parse TREKIS 1D distribution files in a directory.
Scans exclusively for .dat files containing 'spectrum_1d_[axis]' or 'velocity_theta_distr_1d_[axis]'.

Usage:
    python TREKIS_4_1d_spectra_plotting.py
    python TREKIS_4_1d_spectra_plotting.py --dir /path/to/data --xscale log --yscale linear
    python TREKIS_4_1d_spectra_plotting.py -o ./plots --show
"""

import argparse
import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# File matching pattern targeting ONLY .dat files
FILE_PATTERN = re.compile(r"(spectrum_1d_\w+|velocity_theta_distr_1d_\w+).*\.dat$", re.IGNORECASE)


def parse_and_plot_file(filepath, output_dir=None, show_plots=False, xscale="auto", yscale="auto"):
    """Parses a single 1D distribution file and generates plots for all spatial bins with configurable scales."""
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # 1. Extract metadata from header comments
    bins = []
    col_names = []
    col_units = []

    for line in lines:
        if line.startswith("#"):
            clean = line.strip("#").strip()
            tokens = clean.split()
            if not tokens:
                continue

            if not col_names and "Time" in line:
                col_names = tokens
            elif not col_units and any(u in line for u in ["fs", "eV", "rad", "deg"]):
                col_units = tokens
            else:
                try:
                    num_tokens = [float(x) for x in tokens]
                    if len(num_tokens) >= 2:
                        bins = num_tokens
                except ValueError:
                    pass

    # Fallbacks
    if not col_names:
        col_names = ["Time", "Variable", "Distribution"]

    # Detect axis name (e.g., 'z', 'x', 'r') from filename
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    axis_match = re.search(r"1d_([A-Za-z0-9]+)", base_name)
    axis_name = axis_match.group(1).lower() if axis_match else "z"

    # Extract clean base label for titles
    clean_label_match = re.search(r"OUTPUT_(.*?)_1d", base_name)
    if clean_label_match:
        clean_title_base = clean_label_match.group(1).replace("_", " ")
    else:
        clean_title_base = base_name.replace("_", " ")

    # 2. Parse data blocks
    blocks = []
    current_block = []

    for line in lines:
        clean_line = line.strip()
        if clean_line.startswith("#") or not clean_line:
            if current_block:
                blocks.append(current_block)
                current_block = []
        else:
            try:
                row_vals = [float(val) for val in clean_line.split()]
                current_block.append(row_vals)
            except ValueError:
                continue

    if current_block:
        blocks.append(current_block)

    if not blocks:
        print(f"Skipping {os.path.basename(filepath)}: No numerical data found.")
        return

    # 3. Restructure data into DataFrame
    full_data = np.vstack([np.array(b) for b in blocks])
    num_cols = full_data.shape[1]
    num_bins = num_cols - 2

    df_cols = [col_names[0], col_names[1]] + [f"Bin_{i+1}" for i in range(num_bins)]
    df = pd.DataFrame(full_data, columns=df_cols[:num_cols])

    # Dynamic label formatting
    var_label = (
        f"{col_names[1]} ({col_units[1]})"
        if len(col_units) > 1
        else col_names[1]
    ).replace("_", " ")

    base_dist_name = "Distribution"
    if len(col_names) > 2:
        raw_dist = col_names[2].replace("_", " ")
        base_dist_name = re.sub(r"\b(space\s+)?[xyzr]\b", f"space {axis_name.upper()}", raw_dist, flags=re.IGNORECASE)
    
    dist_label = f"{base_dist_name} ({axis_name.upper()})" if "space" not in base_dist_name.lower() else base_dist_name
    time_unit = f" {col_units[0]}" if len(col_units) > 0 else ""

    # Generate title area labels
    bin_titles = []
    if len(bins) >= num_bins + 1:
        for i in range(num_bins):
            bin_titles.append(
                f"{clean_title_base} (area {axis_name}: {bins[i]:.2e}—{bins[i+1]:.2e})"
            )
    else:
        for i in range(num_bins):
            bin_titles.append(f"{clean_title_base} (bin {i+1})")

    # Determine global ranges for auto-scaling evaluation
    x_min, x_max = df[col_names[1]].min(), df[col_names[1]].max()
    bin_cols = [f"Bin_{i+1}" for i in range(num_bins)]
    y_min = df[bin_cols].min().min()
    y_max = df[bin_cols].max().max()

    # Evaluate scale choices based on user input and auto rule (> 1e4 range)
    use_xlog = False
    if xscale == "log":
        use_xlog = True
    elif xscale == "auto":
        if (x_max - x_min > 1e4) and (x_min > 0):
            use_xlog = True

    use_ylog = False
    if yscale == "log":
        use_ylog = True
    elif yscale == "auto":
        if (y_max - y_min > 1e4) and (y_min > 0):
            use_ylog = True

    # 4. Generate plots per spatial bin
    unique_times = df[col_names[0]].unique()
    num_times = len(unique_times)
    ncol_legend = 2 if num_times > 4 else 1

    for bin_idx in range(num_bins):
        bin_col = f"Bin_{bin_idx+1}"

        fig, ax = plt.subplots(figsize=(6, 4.5))

        for t in unique_times:
            sub_df = df[df[col_names[0]] == t]
            ax.plot(
                sub_df[col_names[1]],
                sub_df[bin_col],
                label=f"t = {t}{time_unit}",
            )

        if use_xlog:
            ax.set_xscale("log")
        if use_ylog:
            ax.set_yscale("log")

        ax.set_xlabel(var_label, fontsize=10)
        ax.set_ylabel(dist_label, fontsize=10)
        ax.set_title(bin_titles[bin_idx], fontsize=10, pad=10)
        
        ax.grid(False)
        
        ax.legend(
            loc="upper right",
            ncol=ncol_legend,
            fontsize=8,
            frameon=True,
            framealpha=0.85,
            labelspacing=0.3,
            handletextpad=0.4,
            columnspacing=0.8
        )
        
        plt.tight_layout(pad=1.2)

        save_filename = f"{base_name}_bin_{bin_idx+1}.png"
        out_path = (
            os.path.join(output_dir, save_filename)
            if output_dir
            else save_filename
        )

        plt.savefig(out_path, dpi=300)
        print(f"  -> Saved: {out_path}")

        if show_plots:
            plt.show()
        else:
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Automatically find and plot 1D spectra/velocity data files in a folder."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Directory to scan for data files (default: directory containing this script)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default=None,
        help="Directory to save generated plots (default: same directory as data files)",
    )
    parser.add_argument(
        "--xscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="X-axis scale ('auto' switches to log if xmax - xmin > 1e4)",
    )
    parser.add_argument(
        "--yscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="Y-axis scale ('auto' switches to log if ymax - ymin > 1e4)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display interactive plot windows in addition to saving",
    )

    args = parser.parse_args()

    target_dir = os.path.abspath(args.dir)
    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    out_dir = os.path.abspath(args.output_dir) if args.output_dir else target_dir
    os.makedirs(out_dir, exist_ok=True)

    matched_files = [
        f for f in os.listdir(target_dir) 
        if f.lower().endswith(".dat") and FILE_PATTERN.search(f)
    ]

    if not matched_files:
        print(f"No matching .dat files found in directory: {target_dir}")
        print("Expected files matching: 'spectrum_1d_[axis].dat' or 'velocity_theta_distr_1d_[axis].dat'.")
        return

    print(f"Found {len(matched_files)} matching .dat file(s) in {target_dir}:")
    for fname in matched_files:
        print(f"\nProcessing: {fname}")
        file_path = os.path.join(target_dir, fname)
        parse_and_plot_file(file_path, output_dir=out_dir, show_plots=args.show, xscale=args.xscale, yscale=args.yscale)

    print("\nBatch processing complete.")


if __name__ == "__main__":
    main()