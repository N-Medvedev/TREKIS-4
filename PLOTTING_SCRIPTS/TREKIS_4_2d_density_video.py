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
and compile them into animated video files (MP4 or GIF) with high-contrast options.

Usage:
    python plot_trekis_2d_video.py --cmap turbo --log-colors --ext mp4
    python plot_trekis_2d_video.py --cmap inferno --ext gif
"""

import argparse
import glob
import os
import re
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
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
        return None, None, []

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


def create_animation_from_file(outfolder, fname, verbose=False, cmap="turbo", ext="mp4", fps=10, keep_original_axes=True, log_colors=False):
    filepath = os.path.join(outfolder, fname)
    col_names, col_units, blocks = parse2dblockfile(filepath)

    if not blocks or not col_names or len(col_names) < 3:
        if verbose:
            print(f"Skipping {fname}: Insufficient data or columns for 2D map.")
        return

    if keep_original_axes:
        x_raw, y_raw = col_names[0], col_names[1]
        x_unit, y_unit = f" ({col_units[0]})" if col_units[0] else "", f" ({col_units[1]})" if col_units[1] else ""
    else:
        x_raw, y_raw = col_names[1], col_names[0]
        x_unit, y_unit = f" ({col_units[1]})" if col_units[1] else "", f" ({col_units[0]})" if col_units[0] else ""

    # Strip underscores for display labels
    x_col = x_raw.replace("_", " ")
    y_col = y_raw.replace("_", " ")
    val_col = col_names[2].replace("_", " ")
    val_unit = f" (${col_units[2]}$)" if col_units[2] else ""

    base_name = os.path.splitext(fname)[0]
    title_prefix = base_name.replace("OUTPUT_", "").replace("_", " ")

    global_min = float('inf')
    global_max = float('-inf')
    valid_blocks = []

    custom_xmin = None
    custom_xmax = None
    custom_ymin = None
    custom_ymax = None

    for time_val, time_unit, df_block in blocks:
        if df_block.empty:
            continue
        try:
            pivot_col = x_raw
            pivot_idx = y_raw
            pivot_df = df_block.pivot(index=pivot_idx, columns=pivot_col, values=col_names[2])

            X_full = pivot_df.columns.values
            Y_full = pivot_df.index.values
            Z_full = pivot_df.values

            X = X_full
            Y = Y_full
            Z_vals = Z_full[1:, 1:] if (len(Y_full) > 1 and len(X_full) > 1) else Z_full

            # Evaluate grid limit conditions (> 1e5 / < -1e5) on X boundaries
            if len(X) >= 2:
                x_0, x_1, x_last, x_prev = X[0], X[1], X[-1], X[-2]
                if x_last > 1e5:
                    calculated_xmax = x_prev * 2.0
                    custom_xmax = calculated_xmax if custom_xmax is None else max(custom_xmax, calculated_xmax)
                if x_0 < -1e5:
                    calculated_xmin = max(min(-x_prev, x_1), x_0)
                    custom_xmin = calculated_xmin if custom_xmin is None else min(custom_xmin, calculated_xmin)

            # Evaluate grid limit conditions (> 1e5 / < -1e5) on Y boundaries
            if len(Y) >= 2:
                y_0, y_1, y_last, y_prev = Y[0], Y[1], Y[-1], Y[-2]
                if y_last > 1e5:
                    calculated_ymax = y_prev * 2.0
                    custom_ymax = calculated_ymax if custom_ymax is None else max(custom_ymax, calculated_ymax)
                if y_0 < -1e5:
                    calculated_ymin = max(min(-y_prev, y_1), y_0)
                    custom_ymin = calculated_ymin if custom_ymin is None else min(custom_ymin, calculated_ymin)

            global_min = min(global_min, np.nanmin(Z_vals))
            global_max = max(global_max, np.nanmax(Z_vals))
            valid_blocks.append((time_val, time_unit, X, Y, Z_vals))
        except Exception:
            continue

    if not valid_blocks:
        if verbose:
            print(f"No valid pivotable data blocks found for {fname}.")
        return

    fig, ax = plt.subplots(figsize=(5, 4))
    t_val, t_unit, X, Y, Z = valid_blocks[0]

    # Handle Logarithmic vs Linear Color Normalization for higher contrast
    if log_colors:
        # Prevent log(0) or negative errors by clipping positive bounds safely
        pos_vals = [v[v > 0] for _, _, _, _, v in valid_blocks if np.any(v > 0)]
        vmin_log = min([np.min(pv) for pv in pos_vals]) if pos_vals else 1e-6
        vmax_log = max([np.max(v) for _, _, _, _, v in valid_blocks]) if global_max > 0 else 1.0
        norm = mcolors.LogNorm(vmin=max(vmin_log, 1e-16), vmax=vmax_log)
        mesh = ax.pcolormesh(X, Y, Z, shading='flat', cmap=cmap, norm=norm)
    else:
        mesh = ax.pcolormesh(X, Y, Z, shading='flat', cmap=cmap, vmin=global_min, vmax=global_max)

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(f"{val_col}{val_unit}")

    if custom_xmin is not None or custom_xmax is not None:
        current_xmin, current_xmax = ax.get_xlim()
        ax.set_xlim(
            custom_xmin if custom_xmin is not None else current_xmin,
            custom_xmax if custom_xmax is not None else current_xmax
        )

    if custom_ymin is not None or custom_ymax is not None:
        current_ymin, current_ymax = ax.get_ylim()
        ax.set_ylim(
            custom_ymin if custom_ymin is not None else current_ymin,
            custom_ymax if custom_ymax is not None else current_ymax
        )

    ax.set_xlabel(f"{x_col}{x_unit}")
    ax.set_ylabel(f"{y_col}{y_unit}")
    title_text = ax.set_title(f"{title_prefix}\nt = {t_val:g} {t_unit}", fontsize=10)
    plt.tight_layout()

    def update(frame_idx):
        tv, tu, X_f, Y_f, Z_f = valid_blocks[frame_idx]
        mesh.set_array(Z_f.ravel())
        title_text.set_text(f"{title_prefix}\nt = {tv:g} {tu}")
        return mesh, title_text

    anim = animation.FuncAnimation(fig, update, frames=len(valid_blocks), interval=1000/fps, blit=False)
    out_path = os.path.join(outfolder, f"{base_name}_anim.{ext.lower()}")

    try:
        if ext.lower() == "gif":
            writer = animation.PillowWriter(fps=fps)
            anim.save(out_path, writer=writer, dpi=200)
        else:
            writer = animation.FFMpegWriter(fps=fps, bitrate=1800)
            anim.save(out_path, writer=writer, dpi=200)
        
        if verbose:
            print(f"Saved animation video to: {out_path}")
    except Exception as e:
        print(f"Error saving video for {fname}: {e}")
        if ext.lower() == "mp4":
            print("Tip: MP4 output requires 'ffmpeg' installed on your system. Try running with '--ext gif' if ffmpeg is unavailable.")
    
    plt.close()


def process2ddirectory(directorypath='.', cmap="turbo", ext="mp4", fps=10, keep_original_axes=True, log_colors=False):
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

    print(f"Found {len(target_files)} 2D density/energy .dat file(s) to convert to video:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        create_animation_from_file(directorypath, fname, verbose=True, cmap=cmap, ext=ext, fps=fps, keep_original_axes=keep_original_axes, log_colors=log_colors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process 2D spatial distribution files and compile them into video animations with high-contrast options."
    )
    parser.add_argument(
        "-d", "--dir",
        type=str,
        default=".",
        help="Directory containing 2D .dat files (default: current directory)",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="turbo",
        help="Matplotlib colormap (default: 'turbo', options: 'inferno', 'plasma', 'magma', 'hot', 'jet')",
    )
    parser.add_argument(
        "--log-colors",
        action="store_true",
        help="Use logarithmic color scaling (LogNorm) for vastly improved contrast across multiple orders of magnitude.",
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
        help="Revert axes back to original file order (L on X-axis, R on Y-axis).",
    )
    
    args = parser.parse_args()
    target_dir = os.path.abspath(args.dir)

    if not os.path.isdir(target_dir):
        print(f"Error: Directory '{target_dir}' does not exist.")
        sys.exit(1)

    process2ddirectory(target_dir, cmap=args.cmap, ext=args.ext, fps=args.fps, keep_original_axes=args.keep_original_axes, log_colors=args.log_colors)