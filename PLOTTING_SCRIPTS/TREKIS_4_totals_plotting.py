import os
import glob
import re
from itertools import cycle
import pandas as pd
import matplotlib.pyplot as plt

# Mapping dictionary for translating column headers to full descriptive names in legends
COLUMN_LABELS = {
    # Particle numbers
    'Nph': 'Photons',
    'Ne': 'Electrons',
    'Nh': 'Holes',
    'Np': 'Positrons',
    'NSHI': 'Swift Heavy Ions',
    'Nion': 'Ions',
    'Nmu': 'Muons',
    
    # Energies
    'Eph': 'Photon',
    'Ee': 'Electron',
    'Eh_kin': 'Hole (kinetic) ',
    'Eh_pot': 'Hole (potential)',
    'Ep': 'Positron',
    'Eat': 'Atomic',
    'Eion': 'Ion',
    'ESHI': 'SHI',
    'Etot': 'Total',
}

def get_total_modifier(file_path):
    """
    Extracts the modifier string after 'OUTPUT_total_' and before 'all' or file end.
    Example: 'OUTPUT_total_above_cutoff_all.dat' -> 'above cutoff'
             'OUTPUT_total_all.dat' -> ''
    """
    base = os.path.basename(file_path)
    # Match pattern between 'OUTPUT_total_' and '_all' or '.dat'
    match = re.search(r'OUTPUT_total_(.*?)(?:_all)?\.dat$', base)
    if match:
        modifier = match.group(1).replace('_', ' ').strip()
        # Ignore if it's just 'all' or empty
        if modifier and modifier != 'all':
            return modifier
    return ""

def parse_total_file(file_path):
    """
    Parses a 'total' data file to extract header column names, units, and data.
    """
    column_names = []
    units = []
    data_start_line = 0

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        line_str = line.strip()
        if line_str.startswith('#'):
            clean_line = line_str.lstrip('#').strip()
            parts = clean_line.split()
            
            if 'Time' in parts or 'Nph' in parts or 'Ne' in parts:
                column_names = parts
            elif any(u in parts for u in ['fs', 'a.u.', 'eV']):
                units = parts
            data_start_line = idx + 1
        else:
            break

    df = pd.read_csv(
        file_path,
        skiprows=data_start_line,
        sep=r'\s+',
        names=column_names,
        engine='python'
    )
    
    unit_dict = {}
    if len(units) == len(column_names):
        unit_dict = dict(zip(column_names, units))
    
    return df, unit_dict

def plot_total_file(file_path, output_dir=None):
    """
    Generates two plots (Numbers and Energies) for a given total data file.
    """
    if output_dir is None:
        output_dir = os.path.dirname(file_path) or '.'

    base_name = os.path.basename(file_path).replace('.dat', '')
    
    # Extract modifier (e.g., "above cutoff") for title display
    modifier = get_total_modifier(file_path)
    title_suffix = f" ({modifier})" if modifier else ""

    try:
        df, units = parse_total_file(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    time_col = df.columns[0]  # Usually 'Time'
    time_unit = units.get(time_col, 'fs')
    
    # Separate columns into particle counts (N*) and energies (E*)
    number_cols = [col for col in df.columns if col.startswith('N')]
    energy_cols = [col for col in df.columns if col.startswith('E')]

    linestyles = ['-', '--', '-.', ':']

    # 1. Plot Particle Numbers
    if number_cols:
        plt.figure(figsize=(5, 4))
        style_cycler = cycle(linestyles)
        
        for col in number_cols:
            unit_str = f" ({units[col]})" if col in units else ""
            label_name = COLUMN_LABELS.get(col, col)
            #plt.plot(df[time_col], df[col], label=f"{col}{unit_str}", linestyle=next(style_cycler))
            plt.plot(df[time_col], df[col], label=f"{label_name}", linestyle=next(style_cycler))
        
        plt.xlabel(f"{time_col} ({time_unit})")
        plt.ylabel("Particle (1/incident)")
        plt.title(f"Particle Numbers vs Time{title_suffix}")
        
        plt.legend(loc='best', frameon=True)
        plt.grid(False)
        plt.tight_layout()
        
        num_out_path = os.path.join(output_dir, f"{base_name}_numbers.png")
        plt.savefig(num_out_path, dpi=300)
        plt.close()
        print(f"Saved: {num_out_path}")

    # 2. Plot Energies
    if energy_cols:
        plt.figure(figsize=(5, 4))
        style_cycler = cycle(linestyles)
        
        for col in energy_cols:
            unit_str = f" ({units[col]})" if col in units else ""
            label_name = COLUMN_LABELS.get(col, col)
            #plt.plot(df[time_col], df[col], label=f"{col}{unit_str}", linestyle=next(style_cycler))
            plt.plot(df[time_col], df[col], label=f"{label_name}", linestyle=next(style_cycler))
        
        plt.xlabel(f"{time_col} ({time_unit})")
        plt.ylabel("Energy (eV)")
        plt.title(f"Energies vs Time{title_suffix}")
        
        plt.legend(loc='best', frameon=True)
        plt.grid(False)
        plt.tight_layout()
        
        energy_out_path = os.path.join(output_dir, f"{base_name}_energies.png")
        plt.savefig(energy_out_path, dpi=300)
        plt.close()
        print(f"Saved: {energy_out_path}")

def process_directory(directory_path='.'):
    """
    Finds all 'OUTPUT_total*.dat' files in the given directory and plots them.
    """
    search_pattern = os.path.join(directory_path, "OUTPUT_total*.dat")
    target_files = glob.glob(search_pattern)

    if not target_files:
        print(f"No matching total data files found in '{directory_path}'.")
        return

    print(f"Found {len(target_files)} file(s) matching 'OUTPUT_total*.dat':")
    for file_path in target_files:
        print(f"\nProcessing: {os.path.basename(file_path)}")
        plot_total_file(file_path, output_dir=directory_path)

if __name__ == "__main__":
    import sys
    
    target_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
    process_directory(target_dir)