import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D

def parse_header(filename):
    """
    Reads header lines to extract dimension labels and units.
    """
    with open(filename, "r") as f:
        header1 = f.readline().strip()
        header2 = f.readline().strip()

    h1 = header1[1:].split()
    h2 = header2[1:].split()

    x_label = f"{h1[0]} [{h2[0]}]"
    y_label = f"{h1[1]} [{h2[1]}]"
    front_label = f"{h1[2]} [{h2[2]}]"
    back_label  = f"{h1[3]} [{h2[3]}]"

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

def animate_3d_surface(blocks, quantity="front", filename="output.mp4",
                       cmap="inferno", zlabel="Value", x_label="X", y_label="Y"):
    """
    Generates a 3D animated surface video with a rotating viewpoint.
    """
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')

    writer = FFMpegWriter(fps=10, bitrate=2000)

    with writer.saving(fig, filename, dpi=150):
        for i, (t, block) in enumerate(blocks):

            X, Y, FRONT, BACK = block_to_grid(block)
            Z = FRONT if quantity == "front" else BACK

            # Flatten for triangulated surface plotting
            x = X.ravel()
            y = Y.ravel()
            z = Z.ravel()

            ax.cla()

            ax.plot_trisurf(
                x, y, z,
                cmap=cmap,
                linewidth=0.2,
                antialiased=True,
                edgecolor='none'
            )

            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_zlabel(zlabel)
            ax.set_title(f"{quantity.capitalize()} surface, t = {t:.3f}")

            # Dynamic camera rotation
            ax.view_init(elev=30, azim=40 + i * 2)

            writer.grab_frame()

    plt.close(fig)

def plot_surface_file_3d(filename, output_dir=None):
    """
    Processes a single Surface file and outputs two 3D MP4 videos (front and back).
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

    front_outfile = os.path.join(output_dir, f"{base_name}_front_3d.mp4")
    back_outfile = os.path.join(output_dir, f"{base_name}_back_3d.mp4")

    # Generate 3D animation for Front surface
    animate_3d_surface(
        blocks,
        quantity="front",
        filename=front_outfile,
        cmap="inferno",
        zlabel=front_label,
        x_label=x_label,
        y_label=y_label
    )
    print(f"Saved: {front_outfile}")

    # Generate 3D animation for Back surface
    animate_3d_surface(
        blocks,
        quantity="back",
        filename=back_outfile,
        cmap="viridis",
        zlabel=back_label,
        x_label=x_label,
        y_label=y_label
    )
    print(f"Saved: {back_outfile}")

def process_directory(directory_path='.'):
    """
    Finds all 'OUTPUT_*Surface*.dat' files in the given directory and plots them in 3D.
    """
    search_pattern = os.path.join(directory_path, "OUTPUT_*Surface*.dat")
    target_files = glob.glob(search_pattern)

    if not target_files:
        print(f"No matching surface data files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} file(s) matching 'OUTPUT_*Surface*.dat':")
    for file_path in target_files:
        print(f"\nProcessing 3D surface: {os.path.basename(file_path)}")
        plot_surface_file_3d(file_path, output_dir=directory_path)

if __name__ == "__main__":
    import sys

    # Allow passing directory via command line argument, defaulting to current folder
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    process_directory(target_dir)