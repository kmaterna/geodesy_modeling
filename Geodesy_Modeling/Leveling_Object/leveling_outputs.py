# A series of functions for io of leveling data
# LEVELING INPUT FUNCTIONS FOR CEC SALTON TROUGH LEVELING DATA
import datetime as dt

import numpy as np
from matplotlib import pyplot as plt


# -------------- WRITE FUNCTIONS ------------- #

def plot_leveling_displacements(leveling_stations, outfile):
    """Plot a leveling object"""
    disp = [item.leveling[-1] for item in leveling_stations];
    lon = [item.lon for item in leveling_stations];
    lat = [item.lat for item in leveling_stations];
    plt.figure();
    plt.scatter(lon, lat, c=disp, s=40, cmap='rainbow');
    plt.plot(leveling_stations[0].reflon, leveling_stations[0].reflat, marker='*');
    plt.colorbar();
    plt.savefig(outfile);
    return;


def write_leveling_invertible_format(myLev, idx1, idx2, unc, filename):
    """One header line
    One datum line (automatically first in the leveling array anyway)
    Lon, lat, disp, sigma, 0, 0, 1 (in m)"""
    print("Writing leveling to file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# Displacement for %s to %s: Lon, Lat, disp(m), sigma, 0, 0, 1 \n" %
                (dt.datetime.strftime(myLev[0].dtarray[idx1], "%Y-%m-%d"),
                 dt.datetime.strftime(myLev[0].dtarray[idx2], "%Y-%m-%d")))
    for station in myLev:
        data = station.leveling[idx2] - station.leveling[idx1];
        if ~np.isnan(data):
            ofile.write("%f %f %f %f 0 0 1\n" % (station.lon, station.lat, data, unc))
    ofile.close();
    return;


def plot_simple_leveling(txtfile, plotname):
    print("Plotting leveling in file %s " % plotname);
    [lon, lat, disp] = np.loadtxt(txtfile, unpack=True, skiprows=1, usecols=(0, 1, 2));
    plt.figure(dpi=300);
    plt.scatter(lon, lat, c=disp, s=40, cmap='rainbow')
    plt.colorbar();
    plt.savefig(plotname);
    return;
