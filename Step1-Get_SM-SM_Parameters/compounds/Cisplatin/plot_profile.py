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

# cmap = plt.cm.tab20c
# colors = [cmap(i) for i in np.linspace(0, 1, 20)]

markers = ["o", 'X', '^', 'P', '>', 'd', 'v', 's', '*', '.', 'x', '>'] * 5

fig, main_ax = plt.subplots()
fig.set_size_inches(10, 8)

def find_profile_xvg_files(directory):
    xvg_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == 'profile.xvg':
                xvg_files.append(os.path.join(root, file))
    return xvg_files

for dir in ['./']:
    xvg_files = find_profile_xvg_files(dir)
    shift_value = 0.0
    for ind, file in enumerate(xvg_files):
        with open(file, 'r') as f:
            data = [line.split() for line in f if not line.startswith(('#', '@'))]
            data = list(zip(*data))
            x, y = list(map(float, data[0])), list(map(float, data[1]))
            closest_x = min(x, key=lambda v: abs(v - 1.4))
            index_closest_x = x.index(closest_x)
            shift_value = y[index_closest_x]
            y = [val - shift_value for val in y]
            plt.plot(x, y,   linewidth=1.5, label=os.path.basename(os.path.dirname(file)))
            print(file)
            print(max(y))
            # im1 = plt.imread('./2-Methylpropane_pair_2A.png')
            # im2 = plt.imread('./2-Methylpropane_pair_5A.png')
            # im3 = plt.imread('./2-Methylpropane_pair_10A.png')

            # ax1 = fig.add_axes([0.1, 0.5, 0.15, 0.15], anchor='NE', zorder=0)
            # ax2 = fig.add_axes([0.25, 0.3, 0.15, 0.15], anchor='NE', zorder=0)
            # ax3 = fig.add_axes([0.6, 0.21, 0.25, 0.25], anchor='NE', zorder=0)

            # ax1.imshow(im1)
            # ax2.imshow(im2)
            # ax3.imshow(im3)
            # ax1.axis('off')
            # ax2.axis('off')
            # ax3.axis('off')

            main_ax.set_xlabel('$d$, nm')
            main_ax.set_ylabel('PMF, kcal/mol')

plt.legend(loc='upper right', fontsize=12, ncols = 2)
plt.xlim([0.25, 1.4])
fig.savefig('pmfs_all_molecules.png', transparent=False, bbox_inches='tight')
plt.show()
