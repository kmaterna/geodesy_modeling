
import numpy as np
import matplotlib.pyplot as plt

"""
Assumes a 1d profile of slip in the following format:

# Top_depth(km) bottom_depth(km) slip(mm).
2.50 3.50 0
3.50 4.00 50
4.00 5.00 0
"""


def read_1d_profile(filename):
    print("Reading filename %s" % filename)
    top_depth, bottom_depth, slip_mm = np.loadtxt(filename, skiprows=1, unpack=True)
    return top_depth, bottom_depth, slip_mm


def plot_1d_profile(txt_filename, plot_filename):
    print("Plotting slip profile in %s " % plot_filename)
    tops, bottoms, slips = read_1d_profile(txt_filename)
    plt.figure(figsize=(4, 5), dpi=300)
    x_vals, y_vals = [], []
    for x1, x2, slip in zip(tops, bottoms, slips):
        x_vals.append(slip)
        y_vals.append(x1)
        x_vals.append(slip)
        y_vals.append(x2)
    plt.plot(x_vals, y_vals, linewidth=3, color='black')
    plt.plot([0, np.max(x_vals)], [0, 0], linewidth=0.3, color='black')
    plt.ylim([0, 4.3])
    plt.xlim([-1, 38])
    plt.xlabel('Slip (mm)', fontsize=14)
    plt.ylabel('Depth (km)', fontsize=14)
    plt.gca().invert_yaxis()
    plt.gca().tick_params(axis='both', labelsize=14)
    plt.tight_layout()
    plt.savefig(plot_filename)
    plt.close()
    return
