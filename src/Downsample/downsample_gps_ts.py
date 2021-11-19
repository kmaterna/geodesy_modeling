# June 2020
# A series of functions to help get GPS TS into downsampled format

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
        E0, N0, U0, E1, N1, U1 = subsample_in_time(station, starttime, endtime);
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


def subsample_in_time(station, starttime, endtime):
    """
    Downsample a TS: give us the data points corresponding to the starttime and endtime
    By averaging over a month around each target date.
    return E0, N0, U0, E1, N1, U1
    """
    dE_start, dN_start, dU_start = [], [], [];
    dE_end, dN_end, dU_end = [], [], [];
    for i in range(len(station.dtarray)):
        if abs((station.dtarray[i] - starttime).days) < 30:
            dE_start.append(station.dE[i]);
            dN_start.append(station.dN[i]);
            dU_start.append(station.dU[i]);
        if abs((station.dtarray[i] - endtime).days) < 30:
            dE_end.append(station.dE[i]);
            dN_end.append(station.dN[i]);
            dU_end.append(station.dU[i]);
    if len(dE_start) > 2:
        E0 = np.nanmean(dE_start);
        N0 = np.nanmean(dN_start);
        U0 = np.nanmean(dU_start);
    else:
        E0, N0, U0 = np.nan, np.nan, np.nan;
    if len(dE_end) > 2:
        E1 = np.nanmean(dE_end);
        N1 = np.nanmean(dN_end);
        U1 = np.nanmean(dU_end);
    else:
        E1, N1, U1 = np.nan, np.nan, np.nan;
    return E0, N0, U0, E1, N1, U1;
