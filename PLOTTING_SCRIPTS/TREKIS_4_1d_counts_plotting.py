#!/usr/bin/env python3
"""
Standalone script to batch parse and plot TREKIS 1D restored files with simulation box scaling.
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
import numpy as np


def parse_simulation_box(out_folder):
    """
    Looks for INPUT.txt or NUMERICAL_PARAMETERS.txt in out_folder.
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
        # Line 11 of block -> index 10 (X limits)
        xmin, xmax = parse_two_floats(block_lines[10])

        # Line 12 of block -> index 11 (Y limits)
        ymin, ymax = parse_two_floats(block_lines[11])

        # Line 13 of block -> index 12 (Z limits)
        zmin, zmax = parse_two_floats(block_lines[12])

        return {
            'dx_trans': abs(xmax - xmin),
            'dy_trans': abs(ymax - ymin),
            'dz_trans': abs(zmax - zmin)
        }
    except IndexError:
        print(f"Error: '::: NUMERICAL PARAMETERS :::' block in {param_file} does not contain enough lines.")
        return None
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


def plot_1d_restored_file(out_folder, fname, box_dims, verbose=False, xscale="auto", yscale="auto"):
    axis, cross_section = get_axis_and_cross_section(fname, box_dims)
    if not axis:
        if verbose:
            print(f"Skipping {fname}: Axis (X, Y, or Z) could not be identified from filename.")
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
    
    # Determine type (Density or Energy)
    fname_lower = fname.lower()
    if 'energy' in fname_lower:
        val_type = "Total Energy"
        y_unit = " (eV)"
        y_label = f"Total Energy per Bin{y_unit}"
        suffix = "total_energy_histogram_plot.png"
    else:
        val_type = "Particle Count"
        y_label = "Particle Count per Bin"
        suffix = "particle_count_histogram_plot.png"

    base_name = os.path.splitext(fname)[0]
    title_label = f"{base_name.replace('OUTPUT_', '').replace('_', ' ')} ({val_type})"

    # Determine global ranges for auto-scaling evaluation
    x_min, x_max = df[grid_col].min(), df[grid_col].max()
    
    # Evaluate scale choices based on user input and auto rule (> 1e4 range)
    use_xlog = False
    if xscale == "log":
        use_xlog = True
    elif xscale == "auto":
        if (x_max - x_min > 1e4) and (x_min > 0):
            use_xlog = True

    # Pre-calculate restored values to find true global y range for yscale auto check
    all_restored_max = 0
    all_restored_min = float('inf')
    grouped = df.groupby(time_col, sort=False)
    for _, group in grouped:
        edges = group[grid_col].values
        dens_vals = group[func_col].values[1:]
        if len(edges) >= 2:
            bin_widths = np.diff(edges)
            bin_volumes = cross_section * bin_widths
            restored_values = dens_vals * bin_volumes
            if len(restored_values) > 0:
                all_restored_max = max(all_restored_max, restored_values.max())
                all_restored_min = min(all_restored_min, restored_values.min())

    use_ylog = False
    if yscale == "log":
        use_ylog = True
    elif yscale == "auto":
        if (all_restored_max - all_restored_min > 1e4) and (all_restored_min > 0):
            use_ylog = True

    fig, ax = plt.subplots(figsize=(5, 4))

    if use_xlog:
        ax.set_xscale("log")
    if use_ylog:
        ax.set_yscale("log")

    cmap = plt.cm.plasma
    num_groups = len(grouped)

    custom_xmin = None
    custom_xmax = None

    for i, (time_val, group) in enumerate(grouped):
        color = cmap(i / max(1, num_groups - 1))
        
        edges = group[grid_col].values
        dens_vals = group[func_col].values[1:]  # Drop dummy index 0

        if len(edges) < 2:
            continue

        x_0 = edges[0]
        x_last = edges[-1]
        x_prev = edges[-2]

        if x_last > 1e5:
            calculated_xmax = x_prev * 2.0
            custom_xmax = calculated_xmax if custom_xmax is None else max(custom_xmax, calculated_xmax)

            if x_0 < -1e5:
                calculated_xmin = -x_prev * 2.0
                custom_xmin = calculated_xmin if custom_xmin is None else min(custom_xmin, calculated_xmin)

        bin_widths = np.diff(edges)
        bin_volumes = cross_section * bin_widths
        restored_values = dens_vals * bin_volumes

        ax.stairs(
            restored_values, 
            edges=edges, 
            label=f"t = {time_val:g} {time_unit}", 
            linewidth=1.2, 
            color=color
        )

    if custom_xmin is not None or custom_xmax is not None:
        current_xmin, current_xmax = ax.get_xlim()
        new_xmin = custom_xmin if custom_xmin is not None else current_xmin
        new_xmax = custom_xmax if custom_xmax is not None else current_xmax
        ax.set_xlim(new_xmin, new_xmax)

    ax.set_xlabel(f"{grid_col}{grid_unit}")
    ax.set_ylabel(y_label)
    ax.grid(False)

    if num_groups > 10:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
    else:
        ax.legend(loc="best", fontsize=8)

    plt.title(title_label, wrap=True)
    plt.tight_layout()

    out_path = os.path.join(out_folder, f"{base_name}_{suffix}")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    if verbose:
        print(f"Saved plot to: {out_path}")


def process_1d_restored_files(directory_path='.', xscale="auto", yscale="auto"):
    box_dims = parse_simulation_box(directory_path)
    if not box_dims:
        print("Aborting calculations due to missing box parameters.")
        return

    all_dat_files = glob.glob(os.path.join(directory_path, "*.dat"))
    target_files = []

    for f_path in all_dat_files:
        fname = os.path.basename(f_path)
        
        # Must be 1D
        if not re.search(r'(?:^|_|\b)1d(?:_|\b)', fname, re.IGNORECASE):
            continue

        # Must contain density OR energy
        if not (re.search(r'density', fname, re.IGNORECASE) or re.search(r'energy', fname, re.IGNORECASE)):
            continue

        target_files.append(f_path)

    target_files = sorted(list(set(target_files)))

    if not target_files:
        print(f"No valid 1D density or energy .dat files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} 1D file(s) to process:")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        plot_1d_restored_file(directory_path, fname, box_dims, verbose=True, xscale=xscale, yscale=yscale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and plot 1D restored files with simulation box scaling and configurable axis scales."
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

    process_1d_restored_files(target_dir, xscale=args.xscale, yscale=args.yscale)