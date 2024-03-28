#!/usr/bin/env python3

"""
March 2021
A driver for visualizing the result of a time series calculation
Not really well-tested after a refactor, but it shouldn't be too hard to get working.
"""

import numpy as np
import sys
from . import uavsar_readwrite


def visualize_timedep_timeseries(config_filename, gps_filename, outdir):
    file_dict = get_file_dictionary(config_filename)
    uav_los, uav_lon, uav_lat = file_dict["uavsar_file"], file_dict["uavsar_lon"], file_dict["uavsar_lat"]
    myUAVSAR_TS = uavsar_readwrite.inputs_TS_grd(uav_los, uav_lon, uav_lat)
    [gps_lons, gps_lats, gps_names] = np.loadtxt(gps_filename, unpack=True, dtype={"formats": (float, float, 'U4'),
                                                                                   "names": ("lons", "lats", "name")})
    selected_epochs = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])  # allows you to combine intervals if necessary
    # selected currently breaks because it's a list not slices. But this is not a big deal.
    uavsar_readwrite.total_ts_visualizing(myUAVSAR_TS, gps_lons, gps_lats, gps_names, selected_epochs, outdir)
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


if __name__ == "__main__":
    config_file, gps_filename, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    visualize_timedep_timeseries(config_file, gps_filename, outdir)
