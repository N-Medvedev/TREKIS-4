#!/usr/bin/env python3
"""
Standalone script to batch animate 1D restored files with simulation box scaling and configurable axis scales.
Supports configurable linear or logarithmic scales for X and Y axes with automatic detection.

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
import numpy as np


def parse_simulation_box(out_folder):
    """
    Extracts xmin, xmax, ymin, ymax, zmin, zmax using fixed line indexing:
    Line 11: X bounds (index 10)
    Line 12: Y bounds (index 11)
    Line 13: Z bounds (index 12)
    """
    candidate_files = ["INPUT.txt", "NUMERICAL_PARAMETERS.txt"]
    param_file = None

    for fname in candidate_files:
        fpath = os.path.join(out_folder, fname)
        if os.path.isfile(fpath):
            param_file = fpath
            break

    if not param_file:
        print("Warning: Neither INPUT.txt nor NUMERICAL_PARAMETERS.txt was found.")
        return None

    with open(param_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # Find block start
    block_start = -1
    for i, line in enumerate(lines):
        if "::: NUMERICAL PARAMETERS :::" in line:
            block_start = i
            break

    if block_start == -1:
        print(f"Warning: '::: NUMERICAL PARAMETERS :::' block not found in {param_file}.")
        return None

    block_lines = lines[block_start + 1:]

    def parse_two_floats(line_str):
        clean_line = line_str.split('!')[0].strip()
        tokens = clean_line.split()
        
        floats = []
        for token in tokens:
            token_clean = token.lower().strip()
            if token_clean in ['t', 'f', 'true', 'false']:
                continue
            
            token_clean = token_clean.replace('d', 'e')
            try:
                floats.append(float(token_clean))
            except ValueError:
                continue

        if len(floats) < 2:
            raise ValueError(f"Could not extract 2 boundary floats from line: '{line_str.strip()}'")

        return floats[0], floats[1]

    try:
        xmin, xmax = parse_two_floats(block_lines[10])
        ymin, ymax = parse_two_floats(block_lines[11])
        zmin, zmax = parse_two_floats(block_lines[12])

        return {
            'dx_trans': abs(xmax - xmin),
            'dy_trans': abs(ymax - ymin),
            'dz_trans': abs(zmax - zmin)
        }
    except Exception as e:
        print(f"Error parsing simulation box boundaries from {param_file}: {e}")
        return None


def parse_1d_block_file(file_path):
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


def get_axis_and_cross_section(fname, box_dims):
    fname_lower = fname.lower()
    
    if re.search(r'1d[_\s]*z', fname_lower):
        return 'z', box_dims['dx_trans'] * box_dims['dy_trans']
    elif re.search(r'1d[_\s]*x', fname_lower):
        return 'x', box_dims['dy_trans'] * box_dims['dz_trans']
    elif re.search(r'1d[_\s]*y', fname_lower):
        return 'y', box_dims['dx_trans'] * box_dims['dz_trans']
    else:
        return None, None


def animate_1d_histogram_file(out_folder, fname, box_dims, fps=10, verbose=False, xscale="auto", yscale="auto"):
    axis, cross_section = get_axis_and_cross_section(fname, box_dims)
    if not axis:
        if verbose:
            print(f"Skipping {fname}: Axis (X, Y, or Z) could not be identified.")
        return

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
    
    # Extract original Y unit from column headers (e.g., "1/A^3" or "eV/A^3")
    raw_y_unit = col_units[2] if len(col_units) > 2 else ""

    # Adjust Y-label unit: replace A^3 with incident
    fname_lower = fname.lower()
    if 'energy' in fname_lower:
        val_type = "Total Energy"
        if "A^3" in raw_y_unit or "A3" in raw_y_unit:
            adjusted_unit = raw_y_unit.replace("A^3", "incident").replace("A3", "incident")
        else:
            adjusted_unit = "eV/incident"
        y_label = f"Total Energy per Bin ({adjusted_unit})"
    else:
        val_type = "Particle Count"
        if "A^3" in raw_y_unit or "A3" in raw_y_unit:
            adjusted_unit = raw_y_unit.replace("A^3", "incident").replace("A3", "incident")
        else:
            adjusted_unit = "1/incident"
        y_label = f"Particle Count per Bin ({adjusted_unit})"

    base_name = os.path.splitext(fname)[0]
    title_label = f"{base_name.replace('OUTPUT_', '').replace('_', ' ')} ({val_type})"

    grouped = [group for _, group in df.groupby(time_col, sort=False)]
    time_values = list(df.groupby(time_col, sort=False).groups.keys())

    if not grouped:
        return

    # Pre-calculate integrated bin values and frame-by-frame global limits
    processed_frames = []
    all_y_restored = []
    x_min, x_max = float('inf'), float('-inf')

    for group in grouped:
        edges = group[grid_col].values
        dens_vals = group[func_col].values[1:]  # Drop dummy index 0

        if len(edges) < 2:
            continue

        bin_widths = np.diff(edges)
        bin_volumes = cross_section * bin_widths
        restored_vals = dens_vals * bin_volumes

        processed_frames.append((edges, restored_vals))
        all_y_restored.extend(restored_vals)

        # Handle boundary padding limits
        x_0 = edges[0]
        x_last = edges[-1]
        x_prev = edges[-2]

        curr_xmax = x_prev * 2.0 if x_last > 1e5 else x_last
        curr_xmin = -x_prev * 2.0 if (x_last > 1e5 and x_0 < -1e5) else x_0

        x_min = min(x_min, curr_xmin)
        x_max = max(x_max, curr_xmax)

    if not processed_frames:
        return

    y_min, y_max = min(all_y_restored), max(all_y_restored)

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
        pos_y = [v for v in all_y_restored if v > 0]
        min_pos = min(pos_y) if pos_y else 1e-30
        ax.set_ylim(max(1e-30, min_pos * 0.9), y_max * 1.1)
    else:
        y_padding = (y_max - y_min) * 0.05 if y_max != y_min else 1.0
        ax.set_ylim(max(0, y_min - y_padding), y_max + y_padding)

    ax.set_xlabel(f"{grid_col}{grid_unit}")
    ax.set_ylabel(y_label)
    ax.grid(False)

    time_text = ax.text(0.05, 0.88, '', transform=ax.transAxes, fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.12)
    stairs_container = []

    def init():
        time_text.set_text('')
        return time_text,

    def update(frame_idx):
        edges, restored_vals = processed_frames[frame_idx]
        t_val = time_values[frame_idx]

        while stairs_container:
            stairs_container.pop().remove()

        stairs_poly = ax.stairs(restored_vals, edges=edges, color="tab:blue", linewidth=1.5)
        stairs_container.append(stairs_poly)

        time_text.set_text(f"t = {t_val:g} {time_unit}")
        ax.set_title(f"{title_label}", pad=12, wrap=True)
        return stairs_poly, time_text

    anim = animation.FuncAnimation(
        fig, update, init_func=init,
        frames=len(processed_frames), interval=1000 // fps, blit=False
    )

    out_path = os.path.join(out_folder, f"{base_name}_histogram_anim.mp4")

    try:
        writer = animation.FFMpegWriter(fps=fps, metadata=dict(artist='Matplotlib'), bitrate=3000)
        anim.save(out_path, writer=writer, dpi=200)
        if verbose:
            print(f"Saved HD 1D animation to: {out_path}")
    except Exception as e:
        print(f"Error saving MP4 for {fname}: {e}")
    finally:
        plt.close()


def process_1d_directory_animations(directory_path='.', fps=10, xscale="auto", yscale="auto"):
    box_dims = parse_simulation_box(directory_path)
    if not box_dims:
        print("Aborting animation creation due to missing simulation box parameters.")
        return

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

    print(f"Found {len(target_files)} 1D density/energy .dat file(s) to animate:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        animate_1d_histogram_file(directory_path, fname, box_dims, fps=fps, verbose=True, xscale=xscale, yscale=yscale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and animate 1D restored files with simulation box scaling and configurable axis scales."
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

    process_1d_directory_animations(target_dir, fps=args.fps, xscale=args.xscale, yscale=args.yscale)