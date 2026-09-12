#!/usr/bin/env python3
"""
Standalone script to batch animate 0D distribution files with configurable linear or logarithmic axis scales.
Supports automatic or manual specification of X and Y axis scaling.

Usage:
    python script_name.py
    python script_name.py --dir /path/to/data --fps 20 --xscale log --yscale linear
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd


def parse_0d_block_file(file_path):
    """
    Parses a 0D dat file containing:
    Col 0: Time
    Col 1: Grid Axis (e.g., Theta Distribution)
    Col 2: Function Values
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


def animate_0d_file(out_folder, fname, fps=10, verbose=False, xscale="auto", yscale="auto"):
    file_path = os.path.join(out_folder, fname)
    df, col_names, col_units = parse_0d_block_file(file_path)

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

    grouped = [group.sort_values(by=grid_col) for _, group in df.groupby(time_col)]
    time_values = list(df.groupby(time_col).groups.keys())

    if not grouped:
        return

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

    # Figure setup
    fig, ax = plt.subplots(figsize=(5, 4))
    
    if use_xlog:
        ax.set_xscale("log")
    if use_ylog:
        ax.set_yscale("log")

    if use_xlog:
        pos_x = [v for v in df[grid_col] if v > 0]
        min_x_pos = min(pos_x) if pos_x else 1e-30
        ax.set_xlim(max(1e-30, min_x_pos * 0.9), x_max * 1.1)
    else:
        ax.set_xlim(x_min, x_max)

    if use_ylog:
        pos_y = [v for v in df[func_col] if v > 0]
        min_y_pos = min(pos_y) if pos_y else 1e-30
        ax.set_ylim(max(1e-30, min_y_pos * 0.9), y_max * 1.1)
    else:
        y_padding = (y_max - y_min) * 0.05 if y_max != y_min else 1.0
        ax.set_ylim(y_min - y_padding, y_max + y_padding)

    ax.set_xlabel(f"{grid_col}{grid_unit}")
    ax.set_ylabel(f"{func_col}{func_unit}")
    ax.grid(False)

    line, = ax.plot([], [], lw=1.5, color="tab:blue")
    time_text = ax.text(0.05, 0.88, '', transform=ax.transAxes, fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Expanded margins:
    # left=0.18 reserves space for Y-axis units
    # top=0.85 reserves extra space for long or multi-line titles
    fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.12)

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def update(frame_idx):
        group = grouped[frame_idx]
        t_val = time_values[frame_idx]
        
        line.set_data(group[grid_col], group[func_col])
        time_text.set_text(f"t = {t_val:g} {time_unit}")
        
        # Added pad=12 and wrap=True to prevent top clipping
        ax.set_title(f"{title_label}", pad=12, wrap=True)
        return line, time_text

    anim = animation.FuncAnimation(
        fig, update, init_func=init, 
        frames=len(grouped), interval=1000//fps, blit=True
    )

    out_path = os.path.join(out_folder, f"{base_name}_anim.mp4")
    
    try:
        writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='Matplotlib'), bitrate=3000)
        anim.save(out_path, writer=writer, dpi=200)
        
        if verbose:
            print(f"Saved HD 0D animation to: {out_path}")
    except Exception as e:
        print(f"Error saving MP4 for {fname} (Ensure ffmpeg is installed): {e}")
    finally:
        plt.close()


def process_0d_animations(directory_path='.', fps=20, xscale="auto", yscale="auto"):
    all_dat_files = glob.glob(os.path.join(directory_path, "*.dat"))
    target_files = []

    for f_path in all_dat_files:
        fname = os.path.basename(f_path)
        
        # 1. Ignore DOS and Surface files
        if re.search(r'dos|surface|total', fname, re.IGNORECASE):
            continue

        # 2. Ignore 1D, 2D, or 3D dimensional files
        if re.search(r'(?:^|_|\b)[1-3]d(?:_|\b)', fname, re.IGNORECASE):
            continue

        target_files.append(f_path)

    target_files = sorted(list(set(target_files)))

    if not target_files:
        print(f"No valid 0D .dat files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} 0D .dat file(s) to animate:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        animate_0d_file(directory_path, fname, fps=fps, verbose=True, xscale=xscale, yscale=yscale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and batch animate 0D distribution files with configurable linear or logarithmic axis scales."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=".",
        help="Directory containing 0D .dat files (default: current directory)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=20,
        help="Frames per second for output animation (default: 20)",
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

    process_0d_animations(target_dir, fps=args.fps, xscale=args.xscale, yscale=args.yscale)