import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import shutil
from import_path import path

module_path = os.path.abspath("../../source")
if module_path not in sys.path:
    sys.path.append(module_path)
    
from read_density_series import read_drugs_coordinates_density, read_last_series_data
from compute_partitioning import partitioning

# Constants
if len(sys.argv) != 2:
    print("Usage: python find_matching_parameters.py <DRUG>")
    sys.exit(1)

DRUG = sys.argv[1]
PROTEINS = ["YYYYYY"]
BASE_PATH = path()
DATA_FILE = '../2.1.-Get_AA_partitioning/data/6Y_partitioning_bar_plot_data.txt'
# OUTPUT_FILE = f'data/CG_best_fit/best_multipliers_{DRUG}.txt'
TEMPERATURES_TO_PLOT = [300]
CONCENTRATIONS_TO_PLOT = [350]
MARKERS = ['o', 'v', 's', 'v', 'X', '*', '>', '.']
COLORS = ['#ca0020', '#f4a582', '#92c5de', '#0571b0', '#9A88B3', '#7FB069', '#F9C80E', '#F86624']
MPIPI_RESIDUES = {
    "M": 1, "G": 2, "K": 3, "T": 4, "R": 5,
    "A": 6, "D": 7, "E": 8, "Y": 9, "V": 10,
    "L": 11, "Q": 12, "W": 13, "F": 14, "S": 15,
    "H": 16, "N": 17, "P": 18, "C": 19, "I": 20
}

# Update matplotlib parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'lines.markersize': 18,
    "font.weight": "bold"
})

# Read protein fractions from file
protein_fractions = {}
with open(DATA_FILE, 'r') as file:
    next(file)  # Skip the header
    for line in file:
        parts = line.strip().split(',')
        if len(parts) >= 4:
            protein_drug = parts[0].strip()
            drug_name = parts[1].strip()
            fraction = float(parts[2].strip())
            if protein_drug.upper() in PROTEINS and drug_name == DRUG:
                protein_fractions[protein_drug.upper()] = fraction

partitioned_N_per_protein = [protein_fractions.get(protein, 0.0) for protein in PROTEINS]

# Initialize dictionaries and lists
best_multipliers_dict = {}
N_partitioned_match_dict = {}
N_partitioned_match_list = []
protein_bar_list = []

# Process each protein
for ind_protein, protein in enumerate(PROTEINS):
    fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
    partitioned_N_per_protein_first_script = partitioned_N_per_protein[ind_protein]
    config_path = f"{BASE_PATH}/{protein}/runs/config1_{DRUG}"

    runs = []
    for conc in CONCENTRATIONS_TO_PLOT:
        mol_path = f"{config_path}/mol_{conc}"
        if os.path.exists(mol_path):
            runs += [d for d in os.listdir(mol_path) if d.startswith('multiplier_')]

    N_partitioned = []
    for run in runs:
        for conc in CONCENTRATIONS_TO_PLOT:
            try:
                file_path = f'{config_path}/mol_{conc}/{run}/densities_drug_chunked2_{TEMPERATURES_TO_PLOT[0]}.dat'
                file_path_protein = f'{config_path}/mol_{conc}/{run}/densities_cond_chunked2_{TEMPERATURES_TO_PLOT[0]}.dat'

                coordinate_data, density_data_last_series, _, _ = read_drugs_coordinates_density(file_path)
                coordinate_data_protein, density_data_last_series_protein, _, _ = read_last_series_data(file_path_protein)

                ax.plot(coordinate_data, density_data_last_series, label=f'{run}', markersize=8, alpha=0.8, linewidth=1)

                N = partitioning(file_path)
                N_partitioned.append(N[0] / conc)
                ax.plot(N[3], N[2], linestyle='--', markersize=8, color='grey', alpha=0.6, linewidth=5)

            except Exception as e:
                print(f"Error processing {run} and concentration {conc}: {e}")

    if N_partitioned:
        best_match_index = np.argmin(np.abs(np.array(N_partitioned) - partitioned_N_per_protein_first_script))
        best_match_config = runs[best_match_index]
        best_multipliers_dict[protein] = best_match_config.split('_')[1]
        N_partitioned_match_dict[protein] = N_partitioned[best_match_index]
        N_partitioned_match_list.append(N_partitioned[best_match_index])
        protein_bar_list.append(protein)
    else:
        print(f"No partitioned data found for protein {protein}")

    best_match_path_drug = f'{config_path}/mol_{CONCENTRATIONS_TO_PLOT[0]}/{best_match_config}/densities_drug_chunked2_{TEMPERATURES_TO_PLOT[0]}.dat'
    coordinate_data, density_data_last_series, _, _ = read_drugs_coordinates_density(best_match_path_drug)
    ax.plot(coordinate_data, density_data_last_series, markersize=8, color='r', label=f'best match {best_multipliers_dict[protein]}', linestyle='dotted', alpha=0.5, linewidth=4)
    ax.legend(prop={'size': 8}, loc='upper right', ncols=2)

    ax.set_title(protein.upper() + ' + ' + DRUG)
    ax.set_xlabel('Box Length')
    ax.set_ylabel('Density ($N$/chunk)')
    fig.tight_layout()
    # fig.savefig(f'figures/multipliers_{protein}_{DRUG}.png', dpi=300)
    plt.show()

    # Copy sm_ff.dat file
    best_run_path = f'{config_path}/mol_{CONCENTRATIONS_TO_PLOT[0]}/{best_match_config}/sm_ff.dat'
    if os.path.exists(best_run_path):
        shutil.copy(best_run_path, f'data/CG_best_fit/sm_ff_{DRUG}.dat')

# # Save best multipliers to a file
# with open(OUTPUT_FILE, 'w') as file:
#     for protein, multiplier in best_multipliers_dict.items():
#         file.write(f"{protein}: {multiplier}\n")

# Calculate well depth multipliers and partitioned bar
well_depth_multipliers = {}
N_partitioned_bar = {}
for protein, multiplier in best_multipliers_dict.items():
    for char in protein.lower():
        if char != 'y':
            atom_type = MPIPI_RESIDUES.get(char.upper())
            if atom_type:
                well_depth_multipliers[atom_type] = float(multiplier)
        if char == 'y':
            atom_type = 9

for protein, partitioned in N_partitioned_match_dict.items():
    for char in protein.lower():
        if char != 'y':
            atom_type = MPIPI_RESIDUES.get(char.upper())
            if atom_type:
                N_partitioned_bar[atom_type] = partitioned
