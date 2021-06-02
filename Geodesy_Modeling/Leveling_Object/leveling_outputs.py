# A series of functions for io of leveling data
# LEVELING OUTPUT FUNCTIONS

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm


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


def multi_panel_increment_plot2(leveling_object_list, plotname):
    """
    Incremental leveling plot
    Currently used for Brawley
    """
    f, axarr = plt.subplots(3, 3, figsize=(30, 18));

    idx1, idx2 = 0, 0;
    num_plots = np.shape(leveling_object_list[0].leveling);
    lons = [x.lon for x in leveling_object_list];
    lats = [x.lat for x in leveling_object_list];
    vmin = -6.0;
    vmax = 6.0;
    annotations = ['T1', 'T2', 'T3', 'T4A', 'T4B/C', 'T4D', 'T4E', 'T4F', 'T5'];

    for i in range(1, num_plots[0]):
        later_data = [x.leveling[i] for x in leveling_object_list];
        earlier_data = [x.leveling[i-1] for x in leveling_object_list];
        data = [];
        for j in range(len(earlier_data)):
            data.append((later_data[j]-earlier_data[j]) * 100);  # units of centimeters
        str1 = dt.datetime.strftime(leveling_object_list[0].dtarray[i-1], "%Y-%b");
        str2 = dt.datetime.strftime(leveling_object_list[0].dtarray[i], "%Y-%b");
        label = str1 + " to " + str2;
        print(label)

        single_panel_plot(axarr[idx2][idx1], lons, lats, data, vmin, vmax, annotations[i-1]+': '+label, 30);

        idx1 = idx1+1;
        if idx1 == 3 or idx1 == 6:
            idx2 = idx2+1;
            idx1 = 0;

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));

    f.add_axes([0.9, 0.05, 0.05, 0.9], visible=False);
    cb = plt.colorbar(custom_cmap, aspect=8, fraction=0.5);
    cb.set_label('Incremental Disp. (cm)', fontsize=40);
    cb.ax.tick_params(labelsize=32);

    plt.savefig(plotname);
    return;


def single_panel_plot(ax, lons, lats, data, vmin, vmax, label, plotting_multiplier=200):
    #  A plot within a multi panel plot
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');

    for i in range(len(data)):
        plotting_value = data[i];
        dot_color = custom_cmap.to_rgba(plotting_value);
        markersize = np.sqrt(abs(plotting_value))*plotting_multiplier;
        # if markersize<4:
        # 	markersize=4;
        if markersize > 35:
            markersize = 35;
        ax.plot(lons[i], lats[i], marker='o', markersize=markersize, color=dot_color, fillstyle="full");
        ax.plot(lons[0], lats[0], marker='*', markersize=10, color='black', fillstyle="full");

    ax.set_xticklabels([]);
    ax.set_yticklabels([]);
    ax.set_title(label, fontsize=35);
    return;
