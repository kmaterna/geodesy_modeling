#!/usr/bin/env python3

"""
March 2021
A driver for visualizing the result of a time series calculation
Not really well-tested after a refactor, but it shouldn't be too hard to get working.
"""

import numpy as np
import sys
from . import gridts_object
from matplotlib import pyplot as plt, cm as cm
import datetime as dt
import os
import matplotlib
from geodesy_modeling import general_utils


def visualize_timedep_timeseries(config_filename, gps_filename, outdir):
    file_dict = get_file_dictionary(config_filename)
    uav_los, uav_lon, uav_lat = file_dict["uavsar_file"], file_dict["uavsar_lon"], file_dict["uavsar_lat"]
    myUAVSAR_TS = gridts_object.inputs_TS_grd(uav_los, uav_lon, uav_lat)
    [gps_lons, gps_lats, gps_names] = np.loadtxt(gps_filename, unpack=True, dtype={"formats": (float, float, 'U4'),
                                                                                   "names": ("lons", "lats", "name")})
    selected_epochs = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # allows you to combine intervals if necessary
    # selected currently breaks because it's a list not slices. But this is not a big deal.
    total_ts_visualizing(myUAVSAR_TS, gps_lons, gps_lats, gps_names, selected_epochs, outdir)
    return


def get_file_dictionary(config_filename):
    """GET FILE NAMES"""
    this_dict = {}
    print("Reading file %s " % config_filename)
    ifile = open(config_filename)
    for line in ifile:
        data_type = line.split(':')[0]
        total_data_files = line.split()[1]  # assuming one file per list entry
        this_dict[data_type] = total_data_files
    ifile.close()
    return this_dict


def plot_grid_TS_redblue(TS_GRD_OBJ, outfile, vmin=-50, vmax=200, aspect=1, incremental=False, gps_i=None,
                         gps_j=None, selected=None):
    """Make a nice time series plot.
    With incremental or cumulative displacement data."""
    num_rows_plots = 3
    num_cols_plots = 4

    # With selected, you can select certain intervals and combine others
    xdates = TS_GRD_OBJ.dtarray
    if selected is None:
        selected = np.arange(len(xdates))
    TS_array = TS_GRD_OBJ.TS[selected, :, :]
    xdates = [xdates[i] for i in range(max(selected) + 1) if i in selected]

    f, axarr = plt.subplots(num_rows_plots, num_cols_plots, figsize=(16, 10), dpi=300)

    for i in range(0, len(xdates)):
        # Collect the data for each plot.
        if incremental:
            data = np.subtract(TS_array[i, :, :], TS_array[i - 1, :, :])
        else:
            data = TS_array[i, :, :]
        data = data.T  # for making a nice map with west approximately to the left.
        data = np.fliplr(data)  # for making a nice map with west approximately to the left.

        gps_j_flipped = gps_j
        gps_i_flipped = [np.shape(data)[1] - x for x in gps_i]

        # Plotting now
        rownum, colnum = get_axarr_numbers(num_cols_plots, i)
        if i == 0 and incremental is True:
            axarr[rownum][colnum].set_visible(False)
            continue
        else:
            axarr[rownum][colnum].imshow(data, aspect=aspect, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
            titlestr = dt.datetime.strftime(xdates[i], "%Y-%m-%d")
            axarr[rownum][colnum].plot(gps_i_flipped, gps_j_flipped, 'v', markersize=6, color='black')
            axarr[rownum][colnum].get_xaxis().set_visible(False)
            axarr[rownum][colnum].set_title(titlestr, fontsize=20)

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical')
    cb.set_label('Displacement (mm)', fontsize=18)
    cb.ax.tick_params(labelsize=12)
    plt.savefig(outfile)
    return


def get_axarr_numbers(cols, idx):
    """Given an incrementally counting idx number and a subplot dimension, where is our plot?
    total_plots = rows * cols """
    col_num = np.mod(idx, cols)
    row_num = int(np.floor(idx / cols))
    return row_num, col_num


def plot_pixel_ts(TS, dtarray, i, j, name, outdir):
    pixel_value, pixel_value2 = [], []
    width_pixels = 80
    for date in range(len(dtarray)):
        pixel_value.append(np.nanmean(TS[date, i - width_pixels:i + width_pixels, j - width_pixels:j + width_pixels]))
        pixel_value2.append(TS[date, i, j])
    plt.figure(figsize=(8, 8))
    plt.plot(dtarray, pixel_value, '.--', markersize=12)
    plt.plot(dtarray, pixel_value2, '.--', color='red', markersize=12)
    plt.xlabel("Time")
    plt.ylabel("Displacement (mm)")
    plt.savefig(os.path.join(outdir, name + "_onepixel.png"))
    return


def total_ts_visualizing(myTS, gps_lon, gps_lat, gps_names, selected, outdir):
    """Vizualizations for grid time series format
    GPS points can be added by index. For future work. """
    plot_grid_TS_redblue(myTS, os.path.join(outdir, "increments.png"), vmin=-100, vmax=100, aspect=4,
                         incremental=True, gps_i=None, gps_j=None, selected=selected)
    plot_grid_TS_redblue(myTS, os.path.join(outdir, "full_TS.png"), vmin=-160, vmax=160, aspect=4,
                         incremental=False, gps_i=None, gps_j=None, selected=selected)
    # Comparing InSAR TS with GPS
    for i in range(len(gps_lon)):
        ipt, jpt = general_utils.get_nearest_pixel_in_raster(myTS.lon, myTS.lat, gps_lon[i], gps_lat[i])
        plot_pixel_ts(myTS.TS, myTS.dtarray, ipt, jpt, gps_names[i], outdir)
    return


if __name__ == "__main__":
    config_file, my_gps_filename, my_outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    visualize_timedep_timeseries(config_file, my_gps_filename, my_outdir)
