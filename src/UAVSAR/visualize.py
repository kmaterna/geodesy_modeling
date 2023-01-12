#!/usr/bin/env python3

"""
March 2021
A driver for visualizing the result of a time series calculation
Not really well-tested after a refactor, but it shouldn't be too hard to get working.
"""

import numpy as np
import sys
from . import uavsar_readwrite
from .. import general_utils


def visualize_timedep_timeseries(config_filename, gps_filename, outdir):
    file_dict = general_utils.get_file_dictionary(config_filename);
    uav_los, uav_lon, uav_lat = file_dict["uavsar_file"], file_dict["uavsar_lon"], file_dict["uavsar_lat"];
    myUAVSAR_TS = uavsar_readwrite.inputs_TS_grd(uav_los, uav_lon, uav_lat);
    [gps_lons, gps_lats, gps_names] = np.loadtxt(gps_filename, unpack=True, dtype={"formats": (float, float, 'U4'),
                                                                                   "names": ("lons", "lats", "name")});
    selected_epochs = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);  # allows you to combine intervals if necessary
    # selected currently breaks because it's a list not slices. But this is not a big deal.
    uavsar_readwrite.total_ts_visualizing(myUAVSAR_TS, gps_lons, gps_lats, gps_names, selected_epochs, outdir);
    return;


if __name__ == "__main__":
    config_file = sys.argv[1]
    gps_filename = sys.argv[2]
    outdir = sys.argv[3]
    visualize_timedep_timeseries(config_file, gps_filename, outdir);
