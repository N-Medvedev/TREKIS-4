#!/usr/bin/env python3
"""
Standalone script to batch parse and plot TREKIS 1D density/energy histogram files.
Supports configurable linear or logarithmic scales for X and Y axes with automatic detection.

Usage:
    python script_name.py
    python script_name.py --dir /path/to/data --xscale log --yscale linear
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import pandas as pd


def parse_1d_block_file(file_path):
    """
    Parses a 1D dat file containing:
    Col 0: Time
    Col 1: Bin End / Boundary Coordinates (Edges)
    Col 2: Function Values BEFORE this boundary (Heights)
    """
    col_names = []
    col_units = []
    
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]

    comment_lines = [line[1:].strip() for line in lines if line.startswith("#")]
    
    if len(comment_lines) >= 2:
        col_names = comment_lines[0].split()
        col_units = comment_lines[1].split()
    elif len(comment_lines) == 1:
        col_names = comment_lines[0].split()
        col_units = [""] * len(col_names)

    try:
        df = pd.read_csv(file_path, comment='#', sep=r'\s+', header=None)
    except Exception as e:
        print(f"Error reading data in {file_path}: {e}")
        return None, [], []

    if len(col_names) >= df.shape[1]:
        df.columns = col_names[:df.shape[1]]
    else:
        df.columns = [f"Col_{i}" for i in range(df.shape[1])]
        col_names = list(df.columns)

    if len(col_units) < len(col_names):
        col_units.extend([""] * (len(col_names) - len(col_units)))

    return df, col_names, col_units


def plot_1d_histogram_file(out_folder, fname, verbose=False, xscale="auto", yscale="auto"):
    file_path = os.path.join(out_folder, fname)
    df, col_names, col_units = parse_1d_block_file(file_path)

    if df is None or df.empty or df.shape[1] < 3:
        if verbose:
            print(f"Skipping {fname}: Needs at least 3 columns (Time, Grid, Values).")
        return

    time_col = col_names[0]
    grid_col = col_names[1]
    func_col = col_names[2]

    time_unit = col_units[0] if col_units[0] else "fs"
    grid_unit = f" ({col_units[1]})" if col_units[1] else ""
    func_unit = f" ({col_units[2]})" if col_units[2] else ""

    base_name = os.path.splitext(fname)[0]
    title_label = base_name.replace("OUTPUT_", "").replace("_", " ")

    # Determine global ranges for auto-scaling evaluation
    x_min, x_max = df[grid_col].min(), df[grid_col].max()
    y_min, y_max = df[func_col].min(), df[func_col].max()

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

    fig, ax = plt.subplots(figsize=(5, 4))

    if use_xlog:
        ax.set_xscale("log")
    if use_ylog:
        ax.set_yscale("log")

    grouped = df.groupby(time_col, sort=False)
    cmap = plt.cm.plasma
    num_groups = len(grouped)

    custom_xmin = None
    custom_xmax = None

    for i, (time_val, group) in enumerate(grouped):
        color = cmap(i / max(1, num_groups - 1))
        
        edges = group[grid_col].values
        # Drop the first dummy 0 value to match N-1 bins between N edges
        y_vals = group[func_col].values[1:]

        if len(edges) < 2:
            continue

        x_0 = edges[0]
        x_last = edges[-1]
        x_prev = edges[-2]

        # Check grid limit conditions (> 1e5 / < -1e5)
        if x_last > 1e5:
            calculated_xmax = x_prev * 2.0
            custom_xmax = calculated_xmax if custom_xmax is None else max(custom_xmax, calculated_xmax)

            if x_0 < -1e5:
                calculated_xmin = -x_prev * 2.0
                custom_xmin = calculated_xmin if custom_xmin is None else min(custom_xmin, calculated_xmin)

        # Draw histogram stairs: N edges and N-1 bin heights
        ax.stairs(
            y_vals, 
            edges=edges, 
            label=f"t = {time_val:g} {time_unit}", 
            linewidth=1.2, 
            color=color
        )

    # Apply window bounds if thresholds were exceeded
    if custom_xmin is not None or custom_xmax is not None:
        current_xmin, current_xmax = ax.get_xlim()
        new_xmin = custom_xmin if custom_xmin is not None else current_xmin
        new_xmax = custom_xmax if custom_xmax is not None else current_xmax
        ax.set_xlim(new_xmin, new_xmax)

    ax.set_xlabel(f"{grid_col}{grid_unit}")
    ax.set_ylabel(f"{func_col}{func_unit}")
    ax.grid(False)

    if num_groups > 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
    else:
        ax.legend(loc="best", fontsize=8)

    plt.title(title_label, wrap=True)
    plt.tight_layout()

    out_path = os.path.join(out_folder, f"{base_name}_histogram_plot.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    if verbose:
        print(f"Saved 1D histogram plot to: {out_path}")


def process_1d_directory(directory_path='.', xscale="auto", yscale="auto"):
    all_dat_files = glob.glob(os.path.join(directory_path, "*.dat"))
    target_files = []

    for f_path in all_dat_files:
        fname = os.path.basename(f_path)
        
        if not re.search(r'(?:^|_|\b)1d(?:_|\b)', fname, re.IGNORECASE):
            continue

        if not re.search(r'density|energy', fname, re.IGNORECASE):
            continue

        target_files.append(f_path)

    target_files = sorted(list(set(target_files)))

    if not target_files:
        print(f"No valid 1D density/energy .dat files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} 1D density/energy .dat file(s) to process as histograms:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        plot_1d_histogram_file(directory_path, fname, verbose=True, xscale=xscale, yscale=yscale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and plot 1D histogram files with configurable axis scales."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=".",
        help="Directory containing 1D .dat files (default: current directory)",
    )
    parser.add_argument(
        "--xscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="X-axis scale ('auto' switches to log if xmax - xmin > 1e4 and values are positive)",
    )
    parser.add_argument(
        "--yscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="Y-axis scale ('auto' switches to log if ymax - ymin > 1e4 and values are positive)",
    )

    args = parser.parse_args()
    target_dir = os.path.abspath(args.dir)

    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    process_1d_directory(target_dir, xscale=args.xscale, yscale=args.yscale)