import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import math
import sys
import argparse

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'lines.markersize': 6, 
    "font.weight": "bold"  
})

fig2 = plt.figure(figsize=(10, 8), dpi=80)
ax2 = plt.axes()
colors = ['#543005','#8c510a','#bf812d','#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695']

protein_list = ['YYYYYY']

# Set up argument parser
parser = argparse.ArgumentParser(description='Process some drugs.')
parser.add_argument('drugs', metavar='D', type=str, nargs='+', help='a list of drugs')

# Parse arguments
args = parser.parse_args()

# Get the list of drugs from the command line
drugs = args.drugs

print(f"Processing the following drugs: {drugs}")

def skip_comments_and_metadata(f):
    for line in f:
        if not line.startswith(('@', '#')):
            yield line

def get_partitioning_aa(protein, drug):
    def load_data(file_path):
        with open(file_path) as f:
            return np.loadtxt(skip_comments_and_metadata(f))

    data_set1 = load_data(f"{protein}/{drug}/protein.xvg")
    data_set2 = load_data(f"{protein}/{drug}/drug.xvg")
    error_data = load_data(f"{protein}/{drug}/density_error_drug.xvg")

    y_values_set1 = data_set1[:, 1]
    epsilon = 5.0  
    min_samples = 25   
    dbscan = DBSCAN(eps=epsilon, min_samples=min_samples).fit(y_values_set1.reshape(-1, 1))
    labels = dbscan.labels_

    cluster2_label_set1 = np.unique(labels)[0]  
    y_values_set2_cluster2 = error_data[labels == cluster2_label_set1, 1]
    errors_set2_cluster2 = error_data[labels == cluster2_label_set1, 2]
    
    sum_cluster2_set2 = np.sum(y_values_set2_cluster2)
    total_sum_set2 = np.sum(error_data[:, 1])

    total_sum_set2_error_density = np.sqrt(np.sum(error_data[:, 2] ** 2))
    total_sum_set2_error_volume = 0.05 * total_sum_set2
    total_sum_set2_error = np.sqrt(total_sum_set2_error_density**2 + total_sum_set2_error_volume**2)

    sum_cluster2_set2_error_density = np.sqrt(np.sum(errors_set2_cluster2**2))
    sum_cluster2_set2_error_volume = 0.05 * sum_cluster2_set2
    sum_cluster2_set2_error = np.sqrt(sum_cluster2_set2_error_density**2 + sum_cluster2_set2_error_volume**2)

    fraction_cluster2_set2 = sum_cluster2_set2 / (total_sum_set2) # - sum_cluster2_set2)  # measured here as N_droplet/N_total
    fraction_cluster2_set2_error = fraction_cluster2_set2 * np.sqrt(
        (sum_cluster2_set2_error / sum_cluster2_set2) ** 2 +
        (total_sum_set2_error / total_sum_set2) ** 2
    )
    
    # Normalize data for plotting
    data_set1[:, 1] /= np.max(data_set1[:, 1])
    data_set2[:, 1] /= np.max(data_set2[:, 1])
    error_data[:, 1] /= np.max(data_set2[:, 1])  # normalize error data relative to drug data

    return fraction_cluster2_set2, fraction_cluster2_set2_error, labels, data_set1[:, 0], data_set1[:, 1], data_set2[:, 0], data_set2[:, 1], error_data[:, 1], y_values_set1

partitioned_N_per_protein = {drug: [] for drug in drugs}
partitioned_N_per_protein_errors = {drug: [] for drug in drugs}
cluster_data = {drug: [] for drug in drugs}
protein_profiles = {drug: [] for drug in drugs}
drug_profiles = {drug: [] for drug in drugs}
drug_errors = {drug: [] for drug in drugs}
protein_profiles_raw = {drug: [] for drug in drugs}

for drug in drugs:
    for protein in protein_list:
        fraction_molecules_partitioned, fraction_error, labels, x_coords, y_coords, drug_x_coords, drug_y_coords, drug_errors_vals, y_values_set1 = get_partitioning_aa(protein, drug)
        partitioned_N_per_protein[drug].append(fraction_molecules_partitioned)
        partitioned_N_per_protein_errors[drug].append(fraction_error)
        cluster_data[drug].append(labels)
        protein_profiles_raw[drug].append(y_values_set1)
        protein_profiles[drug].append((x_coords, y_coords))
        drug_profiles[drug].append((drug_x_coords, drug_y_coords))
        drug_errors[drug].append((drug_x_coords, drug_errors_vals))

all_x_coords = [x for drug in drugs for x, _ in protein_profiles[drug]]
all_y_coords = [y for drug in drugs for _, y in protein_profiles[drug]]

x_min, x_max = np.min(all_x_coords), np.max(all_x_coords)
y_min, y_max = 0, 1

with open("data/6Y_partitioning_bar_plot_data.txt", "w") as f:
    f.write("Protein, Drug, Partitioned N, Partitioning Error\n")
    x_labels = [i.split('/')[0] for i in protein_list]
    bar_width = 0.2
    num_proteins = len(protein_list)
    num_drugs = len(drugs)
    positions = np.arange(num_proteins) * (bar_width * num_drugs + 0.5)

    for i, drug in enumerate(drugs):
        pos = positions + i * bar_width
        ax2.bar(pos, partitioned_N_per_protein[drug], bar_width, label=drug, color=colors[i % len(colors)], edgecolor='black', yerr=partitioned_N_per_protein_errors[drug], capsize=5)
        for j, protein in enumerate(protein_list):
            f.write(f"{protein.split('/')[0]}, {drug}, {partitioned_N_per_protein[drug][j]}, {partitioned_N_per_protein_errors[drug][j]}\n")

    ax2.set_xticks(positions + bar_width * (num_drugs - 1) / 2)
    ax2.set_xticklabels(x_labels)
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.legend(loc='best', fontsize=12, edgecolor='black', framealpha=0.9, ncols=2)
    ax2.set_ylabel('Partitioning')
    ax2.set_title(r'Partitioning of Drugs ($N_{\rm droplet}/N_{\rm total}$) in AA sim.')
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
    # fig2.savefig('figures/partitioning_all.png', dpi=300, bbox_inches='tight')
    plt.show()

with open("density_data_AA.txt", "w") as f:
    for i, protein in enumerate(protein_list):
        num_drugs = len(drugs)
        cols = min(num_drugs, 7)
        rows = math.ceil(num_drugs / cols)
        
        fig, axs = plt.subplots(rows, cols, figsize=(20, 6), dpi=80, squeeze=False)
        
        for j, drug in enumerate(drugs):
            row, col = divmod(j, cols)
            ax = axs[row, col]
            
            labels = cluster_data[drug][i]
            x_coords, y_coords = protein_profiles[drug][i]
            protein_density_raw = protein_profiles_raw[drug][i]
            drug_x_coords, drug_y_coords = drug_profiles[drug][i]
            drug_errors_vals = drug_errors[drug][i][1]
            unique_labels = np.unique(labels)
            cluster_colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
            
            ax.set_title(f'{drug}')
            ax.plot(x_coords, y_coords, 'k-', lw=2, label='Protein')
            cluster_protein_densities_array = []

            for k, col in zip(unique_labels, cluster_colors):
                class_member_mask = (labels == k)
                cluster_protein_densities_array.append(protein_density_raw[class_member_mask])
                ax.plot(x_coords[class_member_mask], y_coords[class_member_mask], 'o', markerfacecolor=col, markeredgecolor='k', markersize=10)

            ax.plot(drug_x_coords, drug_y_coords, color='b', label='Small Molecule')
            
            f.write(f"Protein: {protein.split('/')[0]}, Drug: {drug}\n")
            f.write("x_coords: " + ','.join(map(str, x_coords)) + "\n")
            f.write("y_coords: " + ','.join(map(str, y_coords)) + "\n")
            f.write("Drug x_coords: " + ','.join(map(str, drug_x_coords)) + "\n")
            f.write("Drug y_coords: " + ','.join(map(str, drug_y_coords)) + "\n")
            f.write("Drug errors: " + ','.join(map(str, drug_errors_vals)) + "\n")
            f.write("Average Density Condensate: " + str(np.average(np.concatenate(cluster_protein_densities_array))) + "\n\n")
            
            # ax.set_xlim(x_min, x_max)
            # ax.set_ylim(y_min, y_max)

        for j in range(num_drugs, rows * cols):
            fig.delaxes(axs[j // cols, j % cols])

        fig.text(0.5, 0.04, 'Box long axis, nm', ha='center', va='center')
        fig.text(0.04, 0.5, 'Normalized Density', ha='center', va='center', rotation='vertical')
        plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
        # fig.savefig(f'figures/clusters_with_profiles_{protein.split("/")[0]}.pdf', dpi=300, bbox_inches='tight')
        plt.show()
