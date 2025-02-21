import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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


def read_xvg(file_path):
    """
    Reads a .xvg file and extracts the data (skipping comments and metadata).
    If the file has three columns, it returns x, y, and errors.
    """
    x, y, errors = [], [], []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):  # Skip GROMACS comments and metadata
                continue
            data = line.split()
            x.append(float(data[0]))
            y.append(float(data[1]))
            if len(data) > 2:  # Check if there's a third column for errors
                errors.append(float(data[2]))
    return np.array(x), np.array(y), np.array(errors) if errors else None

def get_max_density(file_pattern):
    """
    Finds the maximum density value across all drug density files matching the pattern.
    """
    file_paths = glob.glob(file_pattern)
    max_density = 0
    for file_path in file_paths:
        _, y, _ = read_xvg(file_path)
        max_density = max(max_density, np.max(y))
    return max_density

def plot_combined_with_scaled_protein(file_pattern, average_file, protein_file, output_file, current_directory):
    """
    Plots all density_segment_*.xvg files together, overlays the average with error bars,
    and includes scaled protein.xvg data on the same plot with dual y-axes.
    """
    # Collect file paths matching the given pattern
    file_paths = glob.glob(file_pattern)
    if not file_paths:
        print("No files found matching the pattern:", file_pattern)
        return

    # Determine the maximum drug density for scaling protein data
    max_drug_density = get_max_density(file_pattern)  # Ensure this function is implemented

    # Initialize the plot and the first axis
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    # Parse time intervals from filenames
    time_intervals = []
    for file_path in file_paths:
        file_name = os.path.basename(file_path)
        # Extract start and end times in nanoseconds
        start_time_ns = float(file_name.split("_")[2]) * 1e-3  # Convert to ns
        end_time_ns = float(file_name.split("_")[3].split(".")[0]) * 1e-3  # Convert to ns
        time_intervals.append((file_path, start_time_ns, end_time_ns))

    # Sort by start time
    time_intervals.sort(key=lambda x: x[1])

    # Create a colormap for time intervals
    norm = mcolors.Normalize(vmin=min(t[1] for t in time_intervals), vmax=max(t[1] for t in time_intervals))
    cmap = cm.viridis

    # Plot individual density_segment_*.xvg files
    for file_path, start_time_ns, end_time_ns in time_intervals:
        x, y, _ = read_xvg(file_path)  # Read data from the file
        time_interval_midpoint = (start_time_ns + end_time_ns) / 2  # Use midpoint for color mapping
        color = cmap(norm(time_interval_midpoint))  # Map midpoint to a color
        ax1.plot(x, y, label=f"{start_time_ns:.2f}-{end_time_ns:.2f} ns", color=color, alpha=0.7)
    
    ax1.set_ylabel("Protein Number Density [$N$/chunk]")
    ax1.set_xlabel("Box long axis [nm]")
    # Title and plot configuration
    ax1.set_title(f"YYYYYY + {current_directory}")
    # ax1.grid(True)

    # Display the legend
    ax1.legend(loc='upper left', bbox_to_anchor=(0.8, 1), ncol=2, fontsize=6)
    plt.tight_layout()      
    # Save the plot with folder name included in the filename
    plt.savefig(output_file)
    plt.show()
    
    print(f"Combined plot with dual y-axes saved as: {output_file}")


current_directory = os.path.basename(os.getcwd())

if __name__ == "__main__":
    file_pattern = "density_segment_*_protein.xvg"  
    average_file = "density_error_drug.xvg"   
    protein_file = "protein.xvg"            
    output_file = f"protein_vs_time.pdf"
    
    plot_combined_with_scaled_protein(file_pattern, average_file, protein_file, output_file, current_directory)
