import os
import glob
import re
import matplotlib.pyplot as plt

def parse_dos_line(line_str):
    """
    Parses a line to extract E (col 0), DOS (col 1), and me_eff (col 3).
    Fallback to regex matching if standard whitespace split fails.
    """
    parts = line_str.split()
    if len(parts) >= 4:
        try:
            return float(parts[0]), float(parts[1]), float(parts[3])
        except ValueError:
            pass

    # Regex fallback for scientific notation numbers without clear spacing
    float_pattern = r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?'
    matches = re.findall(float_pattern, line_str)
    
    if len(matches) >= 4:
        return float(matches[0]), float(matches[1]), float(matches[3])
    elif len(matches) == 3:
        return float(matches[0]), float(matches[1]), float(matches[2])
    
    return None, None, None

def parse_dos_file(file_path):
    """
    Parses a DOS data file containing columns:
    E, DOS, k, me_eff
    
    Excludes effective mass values greater than 100 (or less than -100).
    """
    val_E, val_DOS, cond_E, cond_DOS = [], [], [], []
    val_effm_E, val_effm_vals, cond_effm_E, cond_effm_vals = [], [], [], []
    
    mode = None
    has_explicit_modes = False

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # Check if file has explicit mode markers
    for line in lines:
        if line.startswith("# --- Valence") or line.startswith("# --- Conduction"):
            has_explicit_modes = True
            break

    for line in lines:
        line_str = line.strip()
        
        if line_str.startswith("# --- Valence"):
            mode = "valence"
            continue
        elif line_str.startswith("# --- Conduction"):
            mode = "conduction"
            continue
        elif line_str.startswith("#") or not line_str:
            continue

        e, dos, effm = parse_dos_line(line_str)
        if e is None or dos is None:
            continue

        # Determine mode based on header or energy sign
        current_mode = mode
        if not has_explicit_modes:
            current_mode = "valence" if e < 0 else "conduction"

        if current_mode == "valence":
            val_E.append(e)
            val_DOS.append(dos)
            # Filter condition: Exclude effective mass values > 100 (or < -100)
            if effm is not None and abs(effm) <= 100:
                val_effm_E.append(e)
                val_effm_vals.append(effm)
        elif current_mode == "conduction":
            cond_E.append(e)
            cond_DOS.append(dos)
            # Filter condition: Exclude effective mass values > 100 (or < -100)
            if effm is not None and abs(effm) <= 100:
                cond_effm_E.append(e)
                cond_effm_vals.append(effm)

    return (val_E, val_DOS, cond_E, cond_DOS, 
            val_effm_E, val_effm_vals, cond_effm_E, cond_effm_vals)

def plot_processed_dos(out_folder, fname, file_label, out_name, out_extension=".png", verbose=False):
    """
    Load a DOS file and generate a plot of valence and conduction bands,
    together with the effective mass on a secondary axis.
    """
    file_path = os.path.join(out_folder, fname)
    if verbose:
        print("DOS plotting:", out_folder, fname)

    (val_E, val_DOS, cond_E, cond_DOS, 
     val_effm_E, val_effm_vals, cond_effm_E, cond_effm_vals) = parse_dos_file(file_path)

    if not val_E and not cond_E:
        print(f"Warning: No valid numeric data found in {file_path}")
        return

    fig, ax1 = plt.subplots(figsize=(5, 4))

    # Left axis: DOS
    if val_E:
        ax1.plot(val_E, val_DOS, color="skyblue", lw=1.5, label="Valence DOS")
    if cond_E:
        ax1.plot(cond_E, cond_DOS, color="orange", lw=1.5, label="Conduction DOS")
    
    ax1.set_xlabel("Energy (eV)")
    ax1.set_ylabel("DOS (states/eV)")
    ax1.grid(False)
    ax1.set_ylim(0.0, None)

    # Right axis: Effective mass
    if val_effm_vals or cond_effm_vals:
        ax2 = ax1.twinx()
        if val_effm_vals:
            ax2.plot(val_effm_E, val_effm_vals, color="blue", ls="--", lw=1.2, label="Valence m*")
        if cond_effm_vals:
            ax2.plot(cond_effm_E, cond_effm_vals, color="red", ls="--", lw=1.2, label="Conduction m*")
        
        ax2.set_ylabel("Effective mass ($m_e$)")

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8)
    else:
        ax1.legend(loc="best", fontsize=8)

    plt.title(f"{file_label}")
    plt.tight_layout()

    out_path = os.path.join(out_folder, out_name + out_extension)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    if verbose:
        print(f" Saved DOS plot to: {out_path}")

def process_directory(directory_path='.'):
    """
    Finds all matching DOS text/dat files in the directory and plots them.
    Excludes image and video files.
    """
    patterns = ["*DOS*.dat", "*DOS*.txt", "OUTPUT_DOS_*"]
    target_files = []
    
    for pat in patterns:
        target_files.extend(glob.glob(os.path.join(directory_path, pat)))
    
    # Filter out generated images/videos from matching list
    ignore_extensions = ('.png', '.jpg', '.jpeg', '.gif', '.mp4', '.avi')
    target_files = [f for f in set(target_files) if not f.lower().endswith(ignore_extensions)]
    target_files = sorted(target_files)

    if not target_files:
        print(f"No matching DOS data files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} DOS file(s):")
    for file_path in target_files:
        fname = os.path.basename(file_path)
        base_name = os.path.splitext(fname)[0]
        
        # Clean title string: "OUTPUT_DOS_of_Al" -> "DOS of Al"
        clean_label = base_name.replace("OUTPUT_", "").replace("_", " ")
        if clean_label.startswith("DOS DOS"):
            clean_label = clean_label.replace("DOS DOS", "DOS")
            
        plot_processed_dos(
            out_folder=directory_path,
            fname=fname,
            file_label=clean_label,
            out_name=f"{base_name}_plot",
            out_extension=".png",
            verbose=True
        )

if __name__ == "__main__":
    import sys
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    process_directory(target_dir)