import os
import shutil

# Define the function to create directories and files with modified values
def modify_and_copy_files(base_dir, input_filename, start_multiplier, end_multiplier, num_folders):
    # Create a list of multipliers
    multipliers = [start_multiplier + x * (end_multiplier - start_multiplier) / (num_folders - 1) for x in range(num_folders)]
    
    input_filepath = os.path.join(base_dir, input_filename)
    
    # Read the input file
    with open(input_filepath, 'r') as file:
        lines = file.readlines()

    # Get the list of files in the base directory
    files_in_dir = [f for f in os.listdir(base_dir) if os.path.isfile(os.path.join(base_dir, f))]

    # Process each multiplier
    for multiplier in multipliers:
        # Create a new directory for the current multiplier
        dir_name = os.path.join(base_dir, f"multiplier_{multiplier:.2f}")
        os.makedirs(dir_name, exist_ok=True)
        
        # Modify the values in the 4th column and save to a new file in the new directory
        modified_lines = []
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 5 and parts[0] == "pair_coeff":
                try:
                    parts[4] = f"{float(parts[4]) * multiplier:.4f}"
                except ValueError:
                    pass  # Skip lines where conversion to float fails
            modified_lines.append("  ".join(parts) + "\n")
        
        new_filename = os.path.join(dir_name, input_filename)
        with open(new_filename, 'w') as new_file:
            new_file.writelines(modified_lines)
        
        # Copy other files into the new directory
        for file in files_in_dir:
            if file != input_filename:
                shutil.copy(os.path.join(base_dir, file), dir_name)

# Specify the base directory and parameters
base_dir = 'YYYYYY/runs/config1_2-Methylpropane/mol_350'  # Replace with your actual base directory name
input_filename = 'sm_ff.dat'  # Replace with your actual input filename
start_multiplier = 0.5
end_multiplier = 2.0
num_folders =20

# Run the function
modify_and_copy_files(base_dir, input_filename, start_multiplier, end_multiplier, num_folders)
