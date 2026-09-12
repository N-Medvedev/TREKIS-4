import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors

def parse_header(filename):
    """
    Reads header lines to extract dimension labels and units.
    """
    with open(filename, "r") as f:
        header1 = f.readline().strip()
        header2 = f.readline().strip()

    h1 = header1[1:].split()
    h2 = header2[1:].split()

    x_label = f"{h1[0]} ({h2[0]})"
    y_label = f"{h1[1]} ({h2[1]})"
    front_label = f"{h1[2]} ({h2[2]})"
    back_label  = f"{h1[3]} ({h2[3]})"

    return x_label, y_label, front_label, back_label

def parse_block_data(filename):
    """
    Parses block data across time steps from the data file.
    """
    blocks = []
    current_block = []
    current_time = None

    with open(filename, "r") as f:
        lines = f.readlines()[2:]  # Skip header lines

    for line in lines:
        line = line.strip()

        if line.startswith("# Time"):
            parts = line.split()
            current_time = float(parts[-1])
            continue

        if line == "":
            if current_block:
                arr = np.array(current_block, dtype=float)
                blocks.append((current_time, arr))
                current_block = []
            continue

        parts = line.split()
        if len(parts) == 4:
            current_block.append(parts)

    if current_block:
        arr = np.array(current_block, dtype=float)
        blocks.append((current_time, arr))

    return blocks

def block_to_grid(block):
    """
    Converts 1D block data into structured 2D coordinate mesh and field grids.
    """
    block = block[np.lexsort((block[:, 0], block[:, 1]))]

    x = block[:, 0]
    y = block[:, 1]
    front = block[:, 2]
    back = block[:, 3]

    nx = len(np.unique(x))
    ny = len(np.unique(y))

    X = x.reshape(ny, nx)
    Y = y.reshape(ny, nx)
    FRONT = front.reshape(ny, nx)
    BACK = back.reshape(ny, nx)

    return X, Y, FRONT, BACK

def make_animation(data_index, cmap, norm, label, extent, frames, x_label, y_label, outfile):
    """
    Generates a non-blinking FuncAnimation MP4 video.
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    t0, FRONT0, BACK0 = frames[0]
    Z0 = FRONT0 if data_index == 0 else BACK0

    im = ax.imshow(
        Z0, origin="lower", extent=extent, cmap=cmap, norm=norm, aspect="auto"
    )

    cb = fig.colorbar(im, ax=ax)
    cb.set_label(label)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    def update(i):
        t, FRONT, BACK = frames[i]
        Z = FRONT if data_index == 0 else BACK
        im.set_data(Z)
        ax.set_title(f"{label}, t = {t:.3f}")
        return [im]

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), blit=True
    )

    ani.save(outfile, fps=10, dpi=150)
    plt.close(fig)

def plot_surface_file(filename, output_dir=None):
    """
    Processes a single Surface file and outputs two MP4 videos (front and back).
    """
    if output_dir is None:
        output_dir = os.path.dirname(filename) or '.'

    base_name = os.path.basename(filename).replace('.dat', '')

    try:
        x_label, y_label, front_label, back_label = parse_header(filename)
        blocks = parse_block_data(filename)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return

    if not blocks:
        print(f"No block data found in {filename}.")
        return

    # Precompute frames and grid spatial dimensions
    frames = []
    X = Y = None
    for t, block in blocks:
        X, Y, FRONT, BACK = block_to_grid(block)
        frames.append((t, FRONT, BACK))

    extent = [X.min(), X.max(), Y.min(), Y.max()]

    # Dynamic auto-vmax calculation across all frames for normalized coloring
    all_front_max = max(np.max(FRONT) for _, FRONT, _ in frames)
    all_back_max = max(max(np.max(BACK) for _, _, BACK in frames), 1e-5)

    # Use auto-scaling normalization (or set fixed limits if preferred)
    norm_front = colors.Normalize(vmin=0, vmax=all_front_max)
    norm_back = colors.Normalize(vmin=0, vmax=all_back_max)

    # Output paths for MP4s
    front_outfile = os.path.join(output_dir, f"{base_name}_front.mp4")
    back_outfile = os.path.join(output_dir, f"{base_name}_back.mp4")

    # Front surface video
    make_animation(
        0, "inferno", norm_front, front_label, extent, frames, x_label, y_label, front_outfile
    )
    print(f"Saved: {front_outfile}")

    # Back surface video
    make_animation(
        1, "viridis", norm_back, back_label, extent, frames, x_label, y_label, back_outfile
    )
    print(f"Saved: {back_outfile}")

def process_directory(directory_path='.'):
    """
    Finds all 'OUTPUT_*Surface*.dat' files in the given directory and plots them.
    """
    search_pattern = os.path.join(directory_path, "OUTPUT_*Surface*.dat")
    target_files = glob.glob(search_pattern)

    if not target_files:
        print(f"No matching surface data files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} file(s) matching 'OUTPUT_*Surface*.dat':")
    for file_path in target_files:
        print(f"\nProcessing: {os.path.basename(file_path)}")
        plot_surface_file(file_path, output_dir=directory_path)

if __name__ == "__main__":
    import sys

    # Allow passing directory via command line argument, defaulting to current folder
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    process_directory(target_dir)