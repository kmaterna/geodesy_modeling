"""
Functions to get GPS TS into downsampled format
"""

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from ..GNSS_Object import gnss_object


def get_displacements_show_ts(stations, starttime, endtime, gps_sigma, prep_dir):
    """Get the values of TS at starttime and endtime"""
    startlim = starttime - dt.timedelta(days=1065);
    endlim = endtime + dt.timedelta(days=1065);
    gps_displacements_object = [];

    for station in stations:
        [start_pos, end_pos] = subsample_ts_start_end(station, starttime, endtime);
        E0, N0, U0 = start_pos[0], start_pos[1], start_pos[2];
        E1, N1, U1 = end_pos[0], end_pos[1], end_pos[2];
        one_object = gnss_object.Timeseries(name=station.name, coords=station.coords, dtarray=[starttime, endtime],
                                            dN=[0, N1 - N0], dE=[0, E1 - E0], dU=[0, U1 - U0],
                                            Sn=[gps_sigma, gps_sigma], Se=[gps_sigma, gps_sigma],
                                            Su=[3 * gps_sigma, 3 * gps_sigma], EQtimes=station.EQtimes);
        gps_displacements_object.append(one_object);

        f, axarr = plt.subplots(3, 1, figsize=(12, 8), dpi=300);
        axarr[0].plot(station.dtarray, station.dE, '.');
        axarr[0].set_xlim([startlim, endlim]);
        axarr[0].plot(starttime, E0, '.', color='red', markersize=15);
        axarr[0].plot(endtime, E1, '.', color='red', markersize=15);
        axarr[0].plot([starttime, endtime], [E0, E1], color='red');
        axarr[0].set_ylabel('East (mm)')
        axarr[1].plot(station.dtarray, station.dN, '.');
        axarr[1].set_xlim([startlim, endlim]);
        axarr[1].plot(starttime, N0, '.', color='red', markersize=15);
        axarr[1].plot(endtime, N1, '.', color='red', markersize=15);
        axarr[1].set_ylabel('North (mm)');
        axarr[2].plot(station.dtarray, station.dU, '.');
        axarr[2].set_xlim([startlim, endlim]);
        axarr[2].plot(starttime, U0, '.', color='red', markersize=15);
        axarr[2].plot(endtime, U1, '.', color='red', markersize=15);
        axarr[2].set_ylabel('Up (mm)');
        plt.savefig(prep_dir + "gps_" + station.name + "_ts.png");
        plt.close();

    return gps_displacements_object;


def subsample_in_time(station, starttime, window_days=30):
    """
    Downsample TS: return position corresponding a given time by averaging over a month around target date.
    return E0, N0, U0
    """
    dE_start, dN_start, dU_start = [], [], [];
    for i in range(len(station.dtarray)):
        if abs((station.dtarray[i] - starttime).days) < window_days:
            dE_start.append(station.dE[i]);
            dN_start.append(station.dN[i]);
            dU_start.append(station.dU[i]);
    if len(dE_start) > 2:
        E0 = np.nanmean(dE_start);
        N0 = np.nanmean(dN_start);
        U0 = np.nanmean(dU_start);
    else:
        E0, N0, U0 = np.nan, np.nan, np.nan;
    return [E0, N0, U0];


def subsample_ts_start_end(station, starttime, endtime, window_days=30):
    start_pos = subsample_in_time(station, starttime, window_days);
    end_pos = subsample_in_time(station, endtime, window_days);
    return [start_pos, end_pos];
