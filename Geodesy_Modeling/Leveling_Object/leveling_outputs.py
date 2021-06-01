# A series of functions for io of leveling data
# LEVELING OUTPUT FUNCTIONS

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt


# -------------- WRITE FUNCTIONS ------------- #

def plot_leveling_displacements(leveling_stations, outfile):
    """
    Plot a leveling object at the last time-step
    Benchmarks without data at the last time-step are not shown.
    """
    disp = [item.leveling[-1] for item in leveling_stations];
    lon = [item.lon for item in leveling_stations];
    lat = [item.lat for item in leveling_stations];
    plt.figure();
    plt.scatter(lon, lat, c=disp, s=40, cmap='rainbow');
    plt.plot(leveling_stations[0].reflon, leveling_stations[0].reflat, marker='*');
    plt.colorbar(label="Leveling displacements (m)");
    plt.savefig(outfile);
    return;


def write_leveling_invertible_format(myLev, idx1, idx2, unc, filename):
    """
    One header line
    One datum line (assuming the whole network is with respect to the same reference datum)
    Lon, lat, disp, sigma, 0, 0, 1 (in m)
    """
    print("Writing leveling to file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# Displacement for %s to %s: Lon, Lat, disp(m), sigma, 0, 0, 1 \n" %
                (dt.datetime.strftime(myLev[0].dtarray[idx1], "%Y-%m-%d"),
                 dt.datetime.strftime(myLev[0].dtarray[idx2], "%Y-%m-%d")))
    # write reference line, hard coded to be 0
    ofile.write("%f %f 0.0 %f 0 0 1\n" % (myLev[0].reflon, myLev[0].reflat, unc) );
    # write all other lines
    for station in myLev:
        if station.lon == station.reflon and station.lat == station.reflat:
            continue;
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


def basic_ts_plot(leveling_object, plotname):
    plt.figure(figsize=(10, 7), dpi=300);
    for item in leveling_object:
        plt.plot(item.dtarray, item.leveling);
    plt.xlabel('Time');
    plt.ylabel('Leveling (m)');
    plt.savefig(plotname);
    return;


def plot_leveling_slopes(leveling_object, slopes, description, plotname):
    """
    Plot a color map of leveling benchmark slopes
    """
    lon = [item.lon for item in leveling_object];
    lat = [item.lat for item in leveling_object];
    plt.figure();
    plt.scatter(lon, lat, c=slopes, s=40, vmin=-0.020, vmax=0.020, cmap='RdBu');
    plt.plot(leveling_object[0].reflon, leveling_object[0].reflat, marker='*');
    plt.title(description, fontsize=20);
    plt.colorbar(label='Displacement Rate (m/yr)');
    plt.savefig(plotname);
    return;
