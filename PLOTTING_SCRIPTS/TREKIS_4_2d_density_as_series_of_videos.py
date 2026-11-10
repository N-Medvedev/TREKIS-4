#!/usr/bin/env python3
"""
0000000000000000000000000000000000000000000000000000000000000
This file is part of TREKIS-4
available at: https://github.com/N-Medvedev/TREKIS-4
1111111111111111111111111111111111111111111111111111111111111
This module is written by N. Medvedev
in 2026
-----------------------------------
Standalone script to batch parse TREKIS 2D spatial distribution files 
and compile them into animated video files (MP4 or GIF) showing 1D profiles 
along one coordinate across time, split per bin of the second coordinate.

Usage:
    python plot_trekis_1d_profiles_video.py --ext mp4
    python plot_trekis_1d_profiles_video.py --keep-original-axes --ext gif
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd


def parse2dblockfile(filepath):
    """
    Parses a 2D multi-block data file containing timestamp sections:
    Header: Column names and units.
    Blocks: Separated by timestamp lines (e.g., '#       0.000000 fs')
    Columns: [Coord 1 (e.g., L), Coord 2 (e.g., R), Value]
    """
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        return None, None, [], []

    col_names = []
    col_units = []
    
    if lines[0].startswith("#"):
        col_names = lines[0][1:].split()
    
    unit_line = lines[1].lstrip("#").strip()
    if "[" in unit_line or "]" in unit_line or any(u in unit_line for u in ["A", "eV", "fs", "cm"]):
        col_units = unit_line.split()
        start_idx = 2
    else:
        col_units = [""] * len(col_names)
        start_idx = 1

    if len(col_units) < len(col_names):
        col_units.extend([""] * (len(col_names) - len(col_units)))

    # Strip brackets "[" and "]" from unit names
    col_units = [u.strip("[]") for u in col_units]

    blocks = []
    current_time = None
    current_time_unit = "fs"
    current_data = []

    for line in lines[start_idx:]:
        if line.startswith("#"):
            if current_time is not None and current_data:
                blocks.append((current_time, current_time_unit, pd.DataFrame(current_data, columns=col_names[:3])))
                current_data = []
            
            parts = line.lstrip("#").split()
            if len(parts) >= 1:
                try:
                    current_time = float(parts[0])
                    if len(parts) >= 2:
                        current_time_unit = parts[1].strip("[]")
                except ValueError:
                    current_time = 0.0
        else:
            row_vals = line.split()
            if len(row_vals) >= 3:
                try:
                    current_data.append([float(val) for val in row_vals[:3]])
                except ValueError:
                    continue

    if current_time is not None and current_data:
        blocks.append((current_time, current_time_unit, pd.DataFrame(current_data, columns=col_names[:3])))

    return col_names, col_units, blocks


def create_animations_from_file(outfolder, fname, verbose=False, ext="mp4", fps=10, keep_original_axes=True):
    filepath = os.path.join(outfolder, fname)
    col_names, col_units, blocks = parse2dblockfile(filepath)

    if not blocks or not col_names or len(col_names) < 3:
        if verbose:
            print(f"Skipping {fname}: Insufficient data or columns.")
        return

    # Determine X-axis and Y-binning coordinate configuration
    if keep_original_axes:
        # Col 0 is Coord 1 (e.g., L), Col 1 is Coord 2 (e.g., R)
        x_raw, bin_raw = col_names[0], col_names[1]
        x_unit = f" ({col_units[0]})" if col_units[0] else ""
        bin_unit = f" ({col_units[1]})" if col_units[1] else ""
        bin_axis_name = col_names[1]
    else:
        # Col 1 is Coord 2 (e.g., R), Col 0 is Coord 1 (e.g., L)
        x_raw, bin_raw = col_names[1], col_names[0]
        x_unit = f" ({col_units[1]})" if col_units[1] else ""
        bin_unit = f" ({col_units[0]})" if col_units[0] else ""
        bin_axis_name = col_names[0]

    x_col = x_raw.replace("_", " ")
    bin_col = bin_raw.replace("_", " ")
    val_col = col_names[2].replace("_", " ")
    val_unit = f" (${col_units[2]}$)" if col_units[2] else ""

    base_name = os.path.splitext(fname)[0]
    title_prefix = base_name.replace("OUTPUT_", "").replace("_", " ")

    # Pre-process all blocks to build uniform pivot tables across time steps
    processed_blocks = []
    global_ymin = float('inf')
    global_ymax = float('-inf')

    for time_val, time_unit, df_block in blocks:
        if df_block.empty:
            continue
        try:
            # Pivot with bin_raw as index (rows) and x_raw as columns
            pivot_df = df_block.pivot(index=bin_raw, columns=x_raw, values=col_names[2])
            
            X_full = pivot_df.columns.values
            Y_bins_full = pivot_df.index.values
            Z_full = pivot_df.values

            # Drop dummy index 0 row/column matching the parser standard
            X = X_full[1:] if len(X_full) > 1 else X_full
            Y_bins = Y_bins_full[1:] if len(Y_bins_full) > 1 else Y_bins_full
            Z_vals = Z_full[1:, 1:] if (len(Y_bins_full) > 1 and len(X_full) > 1) else Z_full

            global_ymin = min(global_ymin, np.nanmin(Z_vals))
            global_ymax = max(global_ymax, np.nanmax(Z_vals))
            
            processed_blocks.append((time_val, time_unit, X, Y_bins, Z_vals))
        except Exception:
            continue

    if not processed_blocks:
        if verbose:
            print(f"No valid data blocks found for {fname}.")
        return

    # Extract unique bin coordinates from the first valid block to generate individual bin video files
    _, _, X_sample, Y_bins_sample, _ = processed_blocks[0]
    num_bins = len(Y_bins_sample)

    # Compute bin boundaries if possible (similar to 1D script area spacing)
    # If Y_bins has N elements, we can construct N intervals: [Y_bins[i] — Y_bins[i+1]]
    full_bin_edges = Y_bins_full if len(Y_bins_full) > num_bins else Y_bins_sample

    for bin_idx in range(num_bins):
        # Format bin range text for title
        if len(full_bin_edges) >= num_bins + 1:
            bin_label_str = f"({bin_axis_name}: {full_bin_edges[bin_idx]:.2e}—{full_bin_edges[bin_idx+1]:.2e}{bin_unit})"
        else:
            bin_label_str = f"({bin_axis_name} bin {bin_idx+1}: {Y_bins_sample[bin_idx]:g}{bin_unit})"

        fig, ax = plt.subplots(figsize=(6, 4.5))

        # Initialize plot with the first time frame curve
        t0, tu0, X_curr, _, Z_curr = processed_blocks[0]
        y_curve = Z_curr[bin_idx, :]
        
        line, = ax.plot(X_curr, y_curve, color='tab:blue', lw=1.5)

        ax.set_xlabel(f"{x_col}{x_unit}")
        ax.set_ylabel(f"{val_col}{val_unit}")
        
        # Add a buffer margin to Y limits for stable clean animation
        y_range = global_ymax - global_ymin if global_ymax > global_ymin else 1.0
        ax.set_ylim(global_ymin - 0.05 * y_range, global_ymax + 0.05 * y_range)
        ax.set_xlim(X_curr[0], X_curr[-1])

        title_text = ax.set_title(f"{title_prefix}\nBin {bin_idx+1} {bin_label_str}\nt = {t0:g} {tu0}", fontsize=9)
        plt.tight_layout()

        def update(frame_idx):
            tv, tu, X_f, _, Z_f = processed_blocks[frame_idx]
            line.set_ydata(Z_f[bin_idx, :])
            title_text.set_text(f"{title_prefix}\nBin {bin_idx+1} {bin_label_str}\nt = {tv:g} {tu}")
            return line, title_text

        anim = animation.FuncAnimation(fig, update, frames=len(processed_blocks), interval=1000/fps, blit=False)
        
        out_filename = f"{base_name}_bin_{bin_idx+1}_anim.{ext.lower()}"
        out_path = os.path.join(outfolder, out_filename)

        try:
            if ext.lower() == "gif":
                writer = animation.PillowWriter(fps=fps)
                anim.save(out_path, writer=writer, dpi=200)
            else:
                writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
                anim.save(out_path, writer=writer, dpi=200)
            
            if verbose:
                print(f"  -> Saved animation for bin {bin_idx+1}: {out_path}")
        except Exception as e:
            print(f"Error saving video for {fname} bin {bin_idx+1}: {e}")
        
        plt.close(fig)


def process_directory(directorypath='.', ext="mp4", fps=10, keep_original_axes=True):
    all_dat_files = glob.glob(os.path.join(directorypath, "*.dat"))
    target_files = []

    for f_path in all_dat_files:
        fname = os.path.basename(f_path)
        if not re.search(r'(?:^|_|\b)2d(?:_|\b)', fname, re.IGNORECASE):
            continue
        if not re.search(r'density|energy', fname, re.IGNORECASE):
            continue
        target_files.append(f_path)

    target_files = sorted(list(set(target_files)))

    if not target_files:
        print(f"No valid 2D density/energy .dat files found in '{directorypath}'.")
        return

    print(f"Found {len(target_files)} 2D density/energy file(s) to convert to 1D profile videos:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        print(f"\nProcessing file: {fname}")
        create_animations_from_file(directorypath, fname, verbose=True, ext=ext, fps=fps, keep_original_axes=keep_original_axes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert TREKIS 2D distribution files into sets of animated 1D curve videos per spatial bin."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=".",
        help="Directory containing 2D .dat files (default: current directory)",
    )
    parser.add_argument(
        "-ext", "--ext", "--extension",
        type=str,
        default="mp4",
        choices=["mp4", "gif"],
        help="Video file format/extension: 'mp4' (requires ffmpeg) or 'gif' (default: mp4)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=10,
        help="Frames per second for the output video (default: 10)",
    )
    parser.add_argument(
        "--keep-original-axes",
        action="store_true",
        default=True,
        help="Keep original axis arrangement (Coord 1 on X-axis, Coord 2 as bins).",
    )
    
    args = parser.parse_args()
    target_dir = os.path.abspath(args.dir)

    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    process_directory(target_dir, ext=args.ext, fps=args.fps, keep_original_axes=args.keep_original_axes)