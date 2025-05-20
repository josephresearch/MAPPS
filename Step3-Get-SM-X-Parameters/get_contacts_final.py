import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Process drug, parameter value, and atom1 type sigma.')
parser.add_argument('--drug', type=str, default='Aceticacid', help='Name of the drug')
parser.add_argument('--eps_0', type=float, required=True, help='Scaling epsilon parameter value from PMFs')
parser.add_argument('--sigma', type=float, required=True, help='CG sigma value for the drug bead')
args = parser.parse_args()

drug = args.drug
atom1_type_drug_sigma = args.sigma
parameter_value = args.eps_0

chain_systems = [2, 3, 4, 5, 6, 7, 8]
sets = [1, 2, 3]

density_scaled = {
    1: 0.804, 2: 0.703, 3: 0.692, 4: 0.709, 5: 0.991, 
    6: 0.723, 7: 0.710, 8: 0.737, 9: 1.0, 10: 0.747,
    11: 0.770, 12: 0.882, 13: 1.043, 14: 0.967, 15: 0.843, 
    16: 0.909, 17: 0.859, 18: 0.781, 19: 0.784, 20: 0.739
}

mpipi_res = {
    "MET": 1, "GLY": 2, "LYS": 3, "THR": 4, "ARG": 5, "ALA": 6, "ASP": 7, "GLU": 8, "TYR": 9, "VAL": 10,
    "LEU": 11, "GLN": 12, "TRP": 13, "PHE": 14, "SER": 15, "HIS": 16, "ASN": 17, "PRO": 18, "CYS": 19, "ILE": 20
}

mpipi_sigma = {
    "MET": 6.467950871, "GLY": 4.69511024, "LYS": 6.671340524, "THR": 5.889061575, "ARG": 6.839051042,
    "ALA": 5.27007355, "ASP": 5.823521721, "GLU": 6.177670316, "TYR": 6.733634028, "VAL": 6.265999286,
    "LEU": 6.534072422, "GLN": 6.277850874, "TRP": 7.066547962, "PHE": 6.629549396, "SER": 5.412669704,
    "HIS": 6.337782404, "ASN": 5.923349124, "PRO": 5.806155151, "CYS": 5.724362788, "ILE": 6.921680089
}

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "font.size": 20,
    "font.weight": "heavy",
    "axes.labelweight": "heavy"
})

def find_profile_xvg_files(directory, prefix='numcont'):
    return [os.path.join(root, file) for root, _, files in os.walk(directory) for file in files if file.startswith(prefix) and file.endswith('.xvg')]

def process_xvg_files(xvg_files):
    residue_contacts = {}
    times_non_zero = []
    res_names = []

    for file in xvg_files:
        with open(file, 'r') as f:
            data = [line.split() for line in f if not line.startswith(('#', '@'))]
            data = list(zip(*data))
            time, contacts = list(map(float, data[0])), list(map(float, data[1]))
            contacts = [1 if i > 0 else 0 for i in contacts]
            total_contacts = np.sum(contacts)

            res_name = file.split('/')[6].split('_')[1].split('.')[0]
            res_names.append(res_name)
            residue_contacts[res_name] = total_contacts / len(contacts)
            time_intervals = [(time[i+1] - time[i]) / 1000 for i in range(len(time)-1) if contacts[i] != 0 and contacts[i+1] != 0]
            times_non_zero.append(np.sum(time_intervals))
    return residue_contacts, times_non_zero, res_names

def average_results(results):
    num_elements = len(results[0])
    averages = [np.mean([results[j][i] for j in range(len(results))]) for i in range(num_elements)]
    std_devs = [stats.sem([results[j][i] for j in range(len(results))], axis=None, ddof=0) for i in range(num_elements)]
    return averages, std_devs

average_res_contacts_over_sets_AA = []
res_names = []

for chain_system in chain_systems:
    for set_number in sets:
        set_dir = f"single_chain/{drug}/gmx/seq{chain_system}/seq{chain_system}_{set_number}"
        if not os.path.exists(set_dir):
            continue

        xvg_files = find_profile_xvg_files(set_dir)
        if not xvg_files:
            continue

        residue_contacts, times_non_zero, res_names_temp = process_xvg_files(xvg_files)
        res_names.extend(res_names_temp)
        average_res_contacts_over_sets_AA.append(list(residue_contacts.values()))

average_of_res_contacts_AA, std_dev_of_res_contacts_AA = average_results(average_res_contacts_over_sets_AA)

residue_contact_pairs = list(zip(res_names, average_of_res_contacts_AA, std_dev_of_res_contacts_AA))
sorted_residue_contact_pairs = sorted(residue_contact_pairs, key=lambda x: x[1], reverse=True)
sorted_res_names_contacts, sorted_average_of_res_contacts, sorted_std_dev_of_res_contacts = zip(*sorted_residue_contact_pairs)

custom_colors = ['#ffffff', '#f4a582']
cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)

fig6 = plt.figure(figsize=(10, 8))
ax6 = plt.axes()

for i, (res_name, avg_contact) in enumerate(zip(sorted_res_names_contacts, sorted_average_of_res_contacts)):
    color = cmap(i / len(sorted_res_names_contacts))
    ax6.bar(res_name, avg_contact, yerr=sorted_std_dev_of_res_contacts[i], color=color, edgecolor='black')

ax6.tick_params(axis='x', rotation=45)
ax6.set_ylabel('Number of Contacts')
ax6.set_title(f'Average Contacts Number for {drug}, AA')
# fig6.savefig(f'figures/AA_average_contacts_sets_{drug}.png', dpi=300, bbox_inches='tight')

if "TYR" in sorted_res_names_contacts:
    tyr_contact_value = sorted_average_of_res_contacts[sorted_res_names_contacts.index("TYR")]
else:
    raise ValueError("TYR not found in sorted_res_names_contacts")

eps_from_contacts = []
scaled_parameters = []
contacts = []
residues = []

for res, avg_contact in zip(sorted_res_names_contacts, sorted_average_of_res_contacts):
    residue_number = mpipi_res.get(res)
    scaled_parameters.append(avg_contact / tyr_contact_value * density_scaled[residue_number])
    contacts.append(avg_contact)
    residues.append(res)

residue_data = [(res, scaled_param, mpipi_res.get(res)) for res, scaled_param in zip(sorted_res_names_contacts, scaled_parameters)]
sorted_residue_data = sorted(residue_data, key=lambda x: x[2])

for res, scaled_param, residue_number in sorted_residue_data:
    eps_from_contacts_per_res = parameter_value * scaled_param
    sigma_residue = mpipi_sigma.get(res)
    sigma_average = (sigma_residue + atom1_type_drug_sigma) / 2
    print(f"pair_coeff\t{residue_number}\t41\twf/cut\t{eps_from_contacts_per_res:.4f}\t{sigma_average:.4f}\t1\t3\t{3.0*sigma_average:.4f}")
    eps_from_contacts.append(eps_from_contacts_per_res)
