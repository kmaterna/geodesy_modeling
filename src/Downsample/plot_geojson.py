# Takes a geojson file from downsampled pixels (see Tectonic_Utils.geodesy)
# and makes a patches plot.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Polygon


def plot_downsampled_InSAR(pixel_list, plot_name, vmin=None, vmax=None):
    plt.figure(figsize=(16, 8), dpi=300);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # Getting lat and lon information for plotting.
    xdim = [pixel.BL_corner[0] for pixel in pixel_list];
    ydim = [pixel.BL_corner[1] for pixel in pixel_list];

    # Plotting each pixel
    for pixel in pixel_list:
        x_total = [pixel.BL_corner[0], pixel.BL_corner[0], pixel.TR_corner[0], pixel.TR_corner[0]]
        y_total = [pixel.BL_corner[1], pixel.TR_corner[1], pixel.TR_corner[1], pixel.BL_corner[1]]

        fault_vertices = np.column_stack((x_total, y_total));
        patch_color = custom_cmap.to_rgba(pixel.mean * 1000);  # in mm
        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        ax = plt.gca();
        ax.add_patch(mypolygon);

    # Plot formatting
    plt.title('Downsampled %d pixels ' % (len(pixel_list)));
    plt.ylim([np.min(ydim) - 0.05, np.max(ydim) + 0.05]);
    plt.xlim([np.min(xdim) - 0.05, np.max(xdim) + 0.05]);
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('mm', fontsize=22);
    for tick in cb.ax.yaxis.get_ticklabels():
        tick.set_size(18)
    plt.xlabel("Longitude", fontsize=18);
    plt.ylabel("Latitude", fontsize=18);
    plt.savefig(plot_name);
    plt.close();

    # New figure for standard deviation of pixels
    plt.figure(figsize=(16, 8), dpi=300);
    color_boundary_object = matplotlib.colors.Normalize(vmin=0.0, vmax=40, clip=True);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');

    # Getting lat and lon information for plotting.
    xdim = [pixel.BL_corner[0] for pixel in pixel_list];
    ydim = [pixel.BL_corner[1] for pixel in pixel_list];

    # Plotting each pixel
    for pixel in pixel_list:
        x_total = [pixel.BL_corner[0], pixel.BL_corner[0], pixel.TR_corner[0], pixel.TR_corner[0]]
        y_total = [pixel.BL_corner[1], pixel.TR_corner[1], pixel.TR_corner[1], pixel.BL_corner[1]]

        fault_vertices = np.column_stack((x_total, y_total));
        patch_color = custom_cmap.to_rgba(pixel.std * 1000);  # in mm
        mypolygon = Polygon(fault_vertices, color=patch_color, alpha=1.0);
        ax = plt.gca();
        ax.add_patch(mypolygon);

    # Plot formatting
    plt.title('Downsampled %d pixels ' % (len(pixel_list)));
    plt.ylim([np.min(ydim) - 0.05, np.max(ydim) + 0.05]);
    plt.xlim([np.min(xdim) - 0.05, np.max(xdim) + 0.05]);
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    custom_cmap.set_array(np.arange(vmin, vmax, 100));
    cb = plt.colorbar(custom_cmap);
    cb.set_label('mm', fontsize=22);
    for tick in cb.ax.yaxis.get_ticklabels():
        tick.set_size(18)
    plt.xlabel("Longitude", fontsize=18);
    plt.ylabel("Latitude", fontsize=18);
    len_of_ext = len(plot_name.split('.')[-1]) + 1  # moving three spaces away from the end to add the file descriptor
    plt.savefig(plot_name[0:-len_of_ext] + "_std" + plot_name[-len_of_ext:]);
    plt.close();
    return;
