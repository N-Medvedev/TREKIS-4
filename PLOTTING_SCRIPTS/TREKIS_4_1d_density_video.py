#!/usr/bin/env python3
"""
Standalone script to batch animate 1D histogram files with configurable linear or logarithmic axis scales.
Supports automatic or manual specification of X and Y axis scaling.

Usage:
    python script_name.py
    python script_name.py --dir /path/to/data --fps 15 --xscale log --yscale linear
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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


def animate_1d_histogram_file(out_folder, fname, fps=10, verbose=False, xscale="auto", yscale="auto"):
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

    grouped = [group for _, group in df.groupby(time_col, sort=False)]
    time_values = list(df.groupby(time_col, sort=False).groups.keys())

    if not grouped:
        return

    # Determine global axes limits
    all_edges = df[grid_col].values
    all_y = df[func_col].values

    x_min, x_max = all_edges.min(), all_edges.max()
    y_min, y_max = all_y.min(), all_y.max()

    # Check grid boundary conditions (> 1e5 / < -1e5) across the dataset
    for group in grouped:
        edges = group[grid_col].values
        if len(edges) >= 2:
            x_0 = edges[0]
            x_last = edges[-1]
            x_prev = edges[-2]

            if x_last > 1e5:
                x_max = x_prev * 2.0
                if x_0 < -1e5:
                    x_min = -x_prev * 2.0

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

    ax.set_xlim(x_min, x_max)

    if use_ylog:
        pos_y = [v for v in all_y if v > 0]
        min_pos = min(pos_y) if pos_y else 1e-30
        ax.set_ylim(max(1e-30, min_pos * 0.9), y_max * 1.1)
    else:
        y_padding = (y_max - y_min) * 0.05 if y_max != y_min else 1.0
        ax.set_ylim(y_min - y_padding, y_max + y_padding)

    ax.set_xlabel(f"{grid_col}{grid_unit}")
    ax.set_ylabel(f"{func_col}{func_unit}")
    ax.grid(False)

    time_text = ax.text(0.05, 0.88, '', transform=ax.transAxes, fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Padding margins to ensure titles and y-axis labels never clip
    fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.12)

    # Global container to clear and redraw step stairs on frame update
    stairs_container = []

    def init():
        time_text.set_text('')
        return time_text,

    def update(frame_idx):
        group = grouped[frame_idx]
        t_val = time_values[frame_idx]

        edges = group[grid_col].values
        y_vals = group[func_col].values[1:]  # Drop dummy index 0

        # Remove previous frame's stairs line
        while stairs_container:
            stairs_container.pop().remove()

        # Draw updated step histogram
        stairs_poly = ax.stairs(y_vals, edges=edges, color="tab:blue", linewidth=1.5)
        stairs_container.append(stairs_poly)

        time_text.set_text(f"t = {t_val:g} {time_unit}")
        ax.set_title(f"{title_label}", pad=12, wrap=True)
        return stairs_poly, time_text

    anim = animation.FuncAnimation(
        fig, update, init_func=init, 
        frames=len(grouped), interval=1000//fps, blit=False
    )

    out_path = os.path.join(out_folder, f"{base_name}_histogram_anim.mp4")

    try:
        writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='Matplotlib'), bitrate=3000)
        anim.save(out_path, writer=writer, dpi=200)
        if verbose:
            print(f"Saved HD 1D histogram animation to: {out_path}")
    except Exception as e:
        print(f"Error saving MP4 for {fname} (Ensure ffmpeg is installed): {e}")
    finally:
        plt.close()


def process_1d_directory_animations(directory_path='.', fps=10, xscale="auto", yscale="auto"):
    all_dat_files = glob.glob(os.path.join(directory_path, "*.dat"))
    target_files = []

    for f_path in all_dat_files:
        fname = os.path.basename(f_path)
        
        # 1. Require "1d" in filename
        if not re.search(r'(?:^|_|\b)1d(?:_|\b)', fname, re.IGNORECASE):
            continue

        # 2. Require "density" or "energy" in filename
        if not re.search(r'density|energy', fname, re.IGNORECASE):
            continue

        target_files.append(f_path)

    target_files = sorted(list(set(target_files)))

    if not target_files:
        print(f"No valid 1D density/energy .dat files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} 1D density/energy .dat file(s) to animate:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        animate_1d_histogram_file(directory_path, fname, fps=fps, verbose=True, xscale=xscale, yscale=yscale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and batch animate 1D histogram files with configurable linear or logarithmic axis scales."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=".",
        help="Directory containing 1D .dat files (default: current directory)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=10,
        help="Frames per second for output animation (default: 10)",
    )
    parser.add_argument(
        "--xscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="X-axis scale ('auto' switches to log if range > 1e4 and values are positive)",
    )
    parser.add_argument(
        "--yscale",
        type=str,
        choices=["auto", "linear", "log"],
        default="auto",
        help="Y-axis scale ('auto' switches to log if range > 1e4 and values are positive)",
    )

    args = parser.parse_args()
    target_dir = os.path.abspath(args.dir)

    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    process_1d_directory_animations(target_dir, fps=args.fps, xscale=args.xscale, yscale=args.yscale)