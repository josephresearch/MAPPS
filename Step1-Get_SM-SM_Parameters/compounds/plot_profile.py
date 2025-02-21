import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams

rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    'font.size': 22,
    'axes.labelsize': 22,
    'axes.titlesize': 22,
    'lines.markersize': 4
})

colors = [plt.cm.tab20(i) for i in range(14)]
markers = ["o", 'X', '^', 'P', '>', 'd', 'v', 's', '*', '.', 'x', '>'] * 5

def find_profile_xvg_files(directory):
    xvg_files = []
    for root, _, files in os.walk(directory):
        xvg_files.extend(os.path.join(root, file) for file in files if file == 'profile.xvg')
    return xvg_files

def plot_profiles(directory):
    fig, main_ax = plt.subplots(figsize=(10, 8))
    xvg_files = find_profile_xvg_files(directory)

    for ind, file in enumerate(xvg_files):
        with open(file, 'r') as f:
            data = [line.split() for line in f if not line.startswith(('#', '@'))]
        
        x, y = zip(*((float(d[0]), float(d[1])) for d in data))
        closest_x = min(x, key=lambda v: abs(v - 1.4))
        shift_value = y[x.index(closest_x)]
        y = [val - shift_value for val in y]

        main_ax.plot(x, y, color=colors[ind], marker=markers[ind], markersize=4, linewidth=1.5, label=os.path.basename(os.path.dirname(file)))
        print(file)
        print(max(y))

    main_ax.set_xlabel('$d$, nm')
    main_ax.set_ylabel('PMF, kcal/mol')
    main_ax.legend(loc='upper right', fontsize=12, ncols=2)
    main_ax.set_xlim([0.25, 1.0])
    main_ax.set_ylim([-4.0, 3.0])

    fig.savefig('pmfs_all_molecules.png', transparent=False, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_profiles('./')
