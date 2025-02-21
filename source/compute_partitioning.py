import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from read_density_series import find_flat_regions_drug, find_flat_regions_drug_full_profile, read_last_series_data, read_drugs_coordinates_density

def partitioning(file_path):

    path_protein = file_path.replace('drug', 'cond')
    protein_coordinate_data, protein_density_data_last_series, protein_result, protein_midpoint_index = read_last_series_data(path_protein)

    # Find flat regions in the high-density region (condensate)
    if 'YYYYYY' in file_path:
        droplet_start, droplet_end = find_flat_regions_drug_full_profile(protein_density_data_last_series)
    else:
        droplet_start, droplet_end = find_flat_regions_drug(protein_density_data_last_series)

    # Read drug data, including the error values
    updated_coord_drug, drug_average_density, result, average_densities_error = read_drugs_coordinates_density(file_path)

    # Total number of drug molecules in the droplet
    total_number_of_drug_molecules = np.sum(drug_average_density[droplet_start:droplet_end + 1])

    # Error in the total number of drug molecules (propagated from density errors)
    total_number_of_drug_molecules_error = np.sqrt(np.sum(np.array(average_densities_error[droplet_start:droplet_end + 1])**2))

    # Return both the total number of drug molecules and the associated error
    return total_number_of_drug_molecules, total_number_of_drug_molecules_error, drug_average_density[droplet_start:droplet_end + 1], updated_coord_drug[droplet_start:droplet_end + 1], protein_density_data_last_series[droplet_start:droplet_end + 1]
