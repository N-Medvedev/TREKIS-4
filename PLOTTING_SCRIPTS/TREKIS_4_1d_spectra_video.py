#!/usr/bin/env python3
"""
Standalone script to batch parse TREKIS 1D distribution files in a directory 
and generate MP4 animations showing the evolution across time steps for each spatial bin.
Scans exclusively for .dat files containing 'spectrum_1d_[axis]' or 'velocity_theta_distr_1d_[axis]'.
Supports configurable linear or logarithmic scales for X and Y axes with automatic detection.

Usage:
    python TREKIS_4_1d_video_animation.py
    python TREKIS_4_1d_video_animation.py --dir /path/to/data --fps 10 --xscale log --yscale linear
"""

import argparse
import os
import re
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd

# File matching pattern targeting ONLY .dat files
FILE_PATTERN = re.compile(r"(spectrum_1d_\w+|velocity_theta_distr_1d_\w+).*\.dat$", re.IGNORECASE)


def parse_and_animate_file(filepath, output_dir=None, fps=5, xscale="auto", yscale="auto"):
    """Parses a single 1D distribution file and generates an MP4 animation for each spatial bin."""
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
    axis_match = re.search(r"1d_([A-Za-z0-9])", base_name, re.IGNORECASE)
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

    # Determine global axis limits for consistent animation scaling across time steps
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

    unique_times = sorted(df[col_names[0]].unique())

    # 4. Generate MP4 Animation per spatial bin
    for bin_idx in range(num_bins):
        bin_col = f"Bin_{bin_idx+1}"

        fig, ax = plt.subplots(figsize=(6, 4.5))

        if use_xlog:
            ax.set_xscale("log")
        if use_ylog:
            ax.set_yscale("log")

        if use_xlog:
            pos_x = [v for v in df[col_names[1]] if v > 0]
            min_x_pos = min(pos_x) if pos_x else 1e-30
            ax.set_xlim(max(1e-30, min_x_pos * 0.9), x_max * 1.1)
        else:
            ax.set_xlim(x_min, x_max)

        if use_ylog:
            pos_y = [v for v in df[bin_cols].values.flatten() if v > 0]
            min_y_pos = min(pos_y) if pos_y else 1e-30
            ax.set_ylim(max(1e-30, min_y_pos * 0.9), y_max * 1.1)
        else:
            y_range = y_max - y_min if y_max != y_min else 1.0
            y_min_lim = y_min - 0.05 * y_range
            y_max_lim = y_max + 0.1 * y_range
            ax.set_ylim(y_min_lim, y_max_lim)
        
        # Initialize line object for animation
        (line,) = ax.plot([], [], lw=2, color='tab:blue')

        ax.set_xlabel(var_label, fontsize=10)
        ax.set_ylabel(dist_label, fontsize=10)
        ax.grid(False)

        title_text = ax.set_title("", fontsize=10, pad=10)

        def init():
            line.set_data([], [])
            title_text.set_text("")
            return line, title_text

        def update(frame_idx):
            t = unique_times[frame_idx]
            sub_df = df[df[col_names[0]] == t]
            x_data = sub_df[col_names[1]]
            y_data = sub_df[bin_col]
            
            line.set_data(x_data, y_data)
            title_text.set_text(f"{bin_titles[bin_idx]}\nt = {t}{time_unit}")
            return line, title_text

        ani = animation.FuncAnimation(
            fig, update, init_func=init, frames=len(unique_times), interval=1000/fps, blit=True
        )

        save_filename = f"{base_name}_bin_{bin_idx+1}.mp4"
        out_path = (
            os.path.join(output_dir, save_filename)
            if output_dir
            else save_filename
        )

        try:
            # Requires ffmpeg installed on system
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
            ani.save(out_path, writer=writer, dpi=200)  # Re-add dpi=200

            print(f"  -> Saved video: {out_path}")
        except Exception as e:
            print(f"  -> Error saving video {save_filename}: {e}")
            print("     (Make sure ffmpeg is installed and accessible in your system PATH)")

        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Automatically find and create MP4 animations for 1D spectra/velocity data files in a folder."
    )
    parser.add_argument(
        "-d",
        "--dir",
        type=str,
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Directory to scan for data files (default: directory containing this script)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save generated MP4 videos (default: same directory as data files)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=5,
        help="Frames per second for the output MP4 videos (default: 5)",
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

    print(f"Found {len(matched_files)} matching .dat file(s) for animation in {target_dir}:")
    for fname in matched_files:
        print(f"\nProcessing animation for: {fname}")
        file_path = os.path.join(target_dir, fname)
        parse_and_animate_file(file_path, output_dir=out_dir, fps=args.fps, xscale=args.xscale, yscale=args.yscale)

    print("\nBatch video creation complete.")


if __name__ == "__main__":
    main()