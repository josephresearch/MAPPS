import numpy as np

def read_drugs_coordinates_density(file_path):
    """
    Read the data of the last series from the file.

    Parameters:
        file_path (str): The path to the file.

    Returns:
        List of density data for the last series.
    """
    path_protein = file_path.replace('drug', 'cond')
    shifted_coordinates_protein, average_densities_protein, result_protein, midpoint_index = read_last_series_data(path_protein)
    
    start_line = 5
    coordinate_data = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        coord_lines = lines[-50:]
        for line in coord_lines:
            if len(line.split()) >= 4:
                coordinate_data.append(float(line.split()[1]))
                
    densities = []
    for i in range(start_line - 1, len(lines)):
        line_s = lines[i].split()
        if len(line_s) >= 4:
            densities.append(float(line_s[2]))

    result = [densities[i:i + 50] for i in range(0, len(densities), 50)]
    result = result[len(result) // 2:]

    densities_collected = list(zip(*result))
    densities_collected = [list(t) for t in densities_collected]
    average_densities = [np.average(dens_series) for dens_series in densities_collected]
    average_densities_error = [np.std(dens_series) / np.sqrt(len(dens_series)) for dens_series in densities_collected]

    x_values = np.array(coordinate_data)
    midpoint_value = x_values[midpoint_index]
    centered_x_values = x_values - midpoint_value

    xnew, ynew = reverse_profiles(list(centered_x_values), list(average_densities))
    updated_coord_drug = [x + 1 if x < -0.5 else x - 1 if x > 0.5 else x for x in xnew]

    sorted_x_values, sorted_y_values = reorder_data(updated_coord_drug, ynew)
    
    return sorted_x_values, sorted_y_values, result, average_densities_error

def read_last_series_data(file_path):
    """
    Read the data of the last series from the file.

    Parameters:
        file_path (str): The path to the file.

    Returns:
        List of density data for the last series.
    """
    start_line = 5
    coordinate_data = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        coord_lines = lines[-50:]
        for line in coord_lines:
            if len(line.split()) >= 4:
                coordinate_data.append(float(line.split()[1]))
                
    densities = []
    for i in range(start_line - 1, len(lines)):
        line_s = lines[i].split()
        if len(line_s) >= 4:
            densities.append(float(line_s[3]))

    result = [densities[i:i + 50] for i in range(0, len(densities), 50)]
    result = result[len(result) // 2:]

    densities_collected = list(zip(*result))
    densities_collected = [list(t) for t in densities_collected]
    average_densities = [np.average(dens_series) for dens_series in densities_collected]

    midpoint_value, shifted_coordinates, shifted_densities, midpoint_index = find_flat_midpoint_and_shift(coordinate_data, average_densities)
    
    series_shifted = [find_flat_midpoint_and_shift(coordinate_data, set)[2] for set in result]
    
    sorted_x_values, sorted_y_values = reorder_data(shifted_coordinates, shifted_densities)
    return sorted_x_values, sorted_y_values, series_shifted, midpoint_index

def reorder_data(x_values, y_values):
    sorted_indices = np.argsort(x_values)
    return np.array(x_values)[sorted_indices], np.array(y_values)[sorted_indices]

def find_flat_midpoint_and_shift(x_values, y_values, threshold=0.2):
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    
    flat_region_indices = np.where(y_values > threshold)[0]
    if len(flat_region_indices) == 0:
        raise ValueError("No flat region with values above the threshold.")
    
    max_length = 0
    longest_midpoint_index = None
    for i in range(len(flat_region_indices)):
        start_index = flat_region_indices[i]
        end_index = start_index + 1
        while end_index < len(y_values) and y_values[end_index] > threshold:
            end_index += 1
        length = end_index - start_index
        if length > max_length:
            max_length = length
            longest_midpoint_index = (start_index + end_index) // 2
    
    if longest_midpoint_index is None:
        raise ValueError("No flat region with values above the threshold.")
    
    longest_midpoint_value = x_values[longest_midpoint_index]
    centered_x_values = x_values - longest_midpoint_value

    xnew, updated_density = reverse_profiles(list(centered_x_values), list(y_values))
    updated_coordinate = [x + 1 if x < -0.5 else x - 1 if x > 0.5 else x for x in xnew]

    return [longest_midpoint_value, updated_coordinate, updated_density, longest_midpoint_index]

def find_flat_regions(density, threshold=0.1):
    flat_regions = []
    is_flat = False
    start_idx = 0
    for i in range(1, len(density)):
        relative_change = abs((density[i] - density[i - 1]) / density[i - 1]) if density[i - 1] != 0 else 0
        if relative_change <= threshold:
            if not is_flat:
                start_idx = i - 1
                is_flat = True
        else:
            if is_flat:
                flat_regions.append((start_idx, i - 1))
                is_flat = False
    if is_flat:
        flat_regions.append((start_idx, len(density) - 1))
    return flat_regions

def find_flat_regions_drug(density):
    density = np.array(density)
    gradients = np.diff(density)
    left_boundary = np.argmax(gradients)
    right_boundary = np.argmin(gradients)
    closest_index = right_boundary if density[left_boundary] >= density[right_boundary] else np.abs(density[right_boundary:] - density[left_boundary]).argmin() + right_boundary
    return left_boundary, closest_index

def find_flat_regions_drug_full_profile(density, flat_threshold=0.2, slope_threshold=0.01):
    gradients = np.diff(density)
    left_boundary = next((i for i in range(1, len(gradients)) if abs(gradients[i]) > slope_threshold and all(abs(gradients[i-j]) < flat_threshold for j in range(1, i))), None)
    right_boundary = next((i + 1 for i in range(len(gradients) - 1, left_boundary, -1) if abs(gradients[i]) > slope_threshold and all(abs(gradients[i+j]) < flat_threshold for j in range(1, len(gradients) - i))), None)
    return left_boundary, right_boundary

def reverse_profiles(x, y):
    indices_to_remove = [i for i, val in enumerate(x) if val < -0.5 or val > 0.5]
    removed_values = [x.pop(i) for i in reversed(indices_to_remove)]
    removed_y_values = [y.pop(i) for i in reversed(indices_to_remove)]
    x.extend(removed_values)
    y.extend(removed_y_values)
    return x, y
