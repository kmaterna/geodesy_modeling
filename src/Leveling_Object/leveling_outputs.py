"""
Functions for plotting/output of leveling data
"""

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm


# -------------- WRITE FUNCTIONS ------------- #

def plot_leveling_displacements(LevList, outfile, vmin=-0.25, vmax=0.15, scale=False, title=None):
    """
    Map the displacements of a leveling object at the last time-step
    Benchmarks without data at last time-step are not shown.
    Scaling the symbols by displacement size is possible.
    """
    print("Plotting leveling in file %s " % outfile)
    disp = [item.leveling[-1] for item in LevList]
    lon = [item.lon for item in LevList]
    lat = [item.lat for item in LevList]
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow')
    fig = plt.figure(figsize=(8, 6), dpi=300)
    if scale:
        for i, item in enumerate(LevList):
            dot_color = custom_cmap.to_rgba(disp[i])
            plt.plot(item.lon, item.lat, marker='o', markersize=abs(disp[i]) * 60, color=dot_color, fillstyle="full")
    else:
        plt.scatter(lon, lat, c=disp, s=60, cmap='rainbow')
    plt.plot(LevList[0].reflon, LevList[0].reflat, marker='*', markersize=20, color='black')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = fig.colorbar(mappable=custom_cmap, ax=plt.gca())
    cb.set_label("Leveling displacements (m)")
    if title:
        plt.title(title, fontsize=20)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(outfile)
    return


def write_leveling_invertible_format(LevList, idx1, idx2, unc, filename):
    """
    One header line.
    One datum line (assuming whole network is with respect to same reference datum).
    Lon, lat, disp, sigma, 0, 0, 1 (in m).

    :param LevList: list of Leveling objects
    :param idx1: int
    :param idx2: int
    :param unc: float, in meters
    :param filename: string
    """
    print("Writing leveling to file %s " % filename)
    ofile = open(filename, 'w')
    ofile.write("# Displacement for %s to %s: Lon, Lat, disp(m), sigma, 0, 0, 1 \n" %
                (dt.datetime.strftime(LevList[0].dtarray[idx1], "%Y-%m-%d"),
                 dt.datetime.strftime(LevList[0].dtarray[idx2], "%Y-%m-%d")))
    # write reference line, hard coded to be 0
    ofile.write("%f %f 0.0 %f 0 0 1\n" % (LevList[0].reflon, LevList[0].reflat, unc))
    # write all other lines
    for station in LevList:
        if station.lon == station.reflon and station.lat == station.reflat:
            continue
        data = station.leveling[idx2] - station.leveling[idx1]
        if ~np.isnan(data):
            ofile.write("%f %f %f %f 0 0 1\n" % (station.lon, station.lat, data, unc))
    ofile.close()
    return


def write_leveling_slopes(LevList, slopes, filename, unc=0, write_nans=True):
    """
    :param LevList: list of leveling objects
    :param slopes: list of floats
    :param filename: string
    :param unc: float
    :param write_nans: bool, optional. Default True.
    """
    print("Writing leveling to file %s " % filename)
    ofile = open(filename, 'w')
    ofile.write("# Slopes for leveling stations: lon, lat, slope[m/yr], uncertainty[m/yr], lkvE, lkvN, lkvU\n")
    ofile.write("%f %f 0.0 %f 0 0 1\n" % (LevList[0].reflon, LevList[0].reflat, unc))
    for station, slope in zip(LevList, slopes):
        if station.lon == station.reflon and station.lat == station.reflat:
            continue
        if np.isnan(slope) and not write_nans:   # if we ignore nans.
            continue
        ofile.write("%f %f %f %f 0 0 1\n" % (station.lon, station.lat, slope, unc))
    ofile.close()
    return


def plot_simple_leveling(txtfile, plotname):
    """Read basic leveling displacements (one epoch) from a simple text file. Map them with colors."""
    print("Plotting leveling in file %s " % plotname)
    [lon, lat, disp] = np.loadtxt(txtfile, unpack=True, skiprows=1, usecols=(0, 1, 2))
    fig = plt.figure(dpi=300)
    plt.scatter(lon, lat, c=disp, s=40, cmap='rainbow')
    _cb = fig.colorbar(ax=plt.gca())
    plt.savefig(plotname)
    return


def basic_ts_plot(LevList, plotname, title=None):
    """Plot the traces of time-dependent position at each leveling benchmark.

    :param LevList: list of Leveling objects
    :param plotname: string
    :param title: string, optional
    """
    plt.figure(figsize=(10, 7), dpi=300)
    for item in LevList:
        plt.plot(item.dtarray, item.leveling)
    plt.xlabel('Time')
    plt.ylabel('Leveling (m)')
    if title:
        plt.title(title)
    else:
        plt.title("Timeseries from %d leveling stations " % (len(LevList)) )
    plt.savefig(plotname)
    return


def plot_leveling_slopes(LevList, slopes, description, plotname):
    """
    Plot a color map of leveling benchmark slopes (slopes calculated outside this function)
    """
    print("Making plot %s " % plotname)
    lon = [item.lon for item in LevList]
    lat = [item.lat for item in LevList]
    fig = plt.figure()
    sc = plt.scatter(lon, lat, c=slopes, s=40, vmin=-0.020, vmax=0.020, cmap='RdBu')
    plt.plot(LevList[0].reflon, LevList[0].reflat, marker='*')
    plt.title(description, fontsize=20)
    _cb = fig.colorbar(mappable=sc, label='Displacement Rate (m/yr)', ax=plt.gca())
    plt.savefig(plotname)
    return


def multi_panel_increment_plot_brawley(LevList, plotname, fields_lon, fields_lat):
    """
    Incremental leveling plot for Brawley.  The annotations are specific to Brawley.
    """
    print("Plotting leveling in file %s " % plotname)
    f, axarr = plt.subplots(3, 3, figsize=(30, 18))

    idx1, idx2 = 0, 0
    num_plots = np.shape(LevList[0].leveling)
    lons = [x.lon for x in LevList]
    lats = [x.lat for x in LevList]
    vmin = -6.0
    vmax = 6.0
    annotations = ['T1', 'T2', 'T3', 'T4A', 'T4B/C', 'T4D', 'T4E', 'T4F', 'T5']

    for i in range(1, num_plots[0]):
        later_data = [x.leveling[i] for x in LevList]
        earlier_data = [x.leveling[i-1] for x in LevList]
        data = []
        for j in range(len(earlier_data)):
            data.append((later_data[j]-earlier_data[j]) * 100)  # units of centimeters
        str1 = dt.datetime.strftime(LevList[0].dtarray[i - 1], "%Y-%b")
        str2 = dt.datetime.strftime(LevList[0].dtarray[i], "%Y-%b")
        label = str1 + " to " + str2
        print(label)

        single_panel_plot(axarr[idx2][idx1], lons, lats, data, vmin, vmax, annotations[i-1]+': '+label, 17)
        axarr[idx2][idx1].plot(fields_lon, fields_lat, linewidth=2, color='darkred')

        if idx2 < 2:
            axarr[idx2][idx1].set_xticklabels([])
        if idx1 > 0:
            axarr[idx2][idx1].set_yticklabels([])

        idx1 = idx1+1
        if idx1 == 3 or idx1 == 6:
            idx2 = idx2+1
            idx1 = 0

    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r')
    custom_cmap.set_array(np.arange(vmin, vmax))

    f.add_axes([0.9, 0.05, 0.05, 0.9], visible=False)
    cb = plt.colorbar(custom_cmap, aspect=8, fraction=0.5)
    cb.set_label('Incremental Disp. (cm)', fontsize=40)
    cb.ax.tick_params(labelsize=32)

    plt.savefig(plotname)
    return


def multi_panel_increment_plot(LevList, outname, vmin=-0.07, vmax=0.07):
    """Incremental plot for leveling"""
    print("Plotting leveling in file %s " % outname)

    f, axarr = plt.subplots(2, 5, figsize=(30, 18))

    idx1, idx2 = 0, 0
    num_plots = len(LevList[0].leveling)
    lons = [x.lon for x in LevList]
    lats = [x.lat for x in LevList]

    for i in range(1, num_plots):
        later_data = [x.leveling[i] for x in LevList]
        earlier_data = [x.leveling[i-1] for x in LevList]
        data = []
        for j in range(len(earlier_data)):
            data.append(later_data[j] - earlier_data[j])
        str1 = dt.datetime.strftime(LevList[0].dtarray[i - 1], "%Y-%m")
        str2 = dt.datetime.strftime(LevList[0].dtarray[i], "%Y-%m")
        label = str1 + " to " + str2
        print(label)

        single_panel_plot(axarr[idx2][idx1], lons, lats, data, vmin, vmax, label, 150)
        idx1 = idx1 + 1
        if idx1 == 5:
            idx1 = 0
            idx2 = 1

    axarr[1][4].set_visible(False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=8, fraction=0.5)
    cb.set_label('Incremental Displacement (m)', fontsize=25)
    cb.ax.tick_params(labelsize=20)

    plt.savefig(outname)
    return


def cumulative_multi_panel_plot(LevList, outname, vmin=-0.30, vmax=0.20):
    """Generalized cumulative multi-panel plot for leveling data"""
    print("Plotting leveling in file %s " % outname)

    f, axarr = plt.subplots(2, 5, figsize=(30, 18))  # currently hard-coded to have these plot dimensions

    idx1, idx2 = 0, 0
    num_plots = len(LevList[0].leveling)
    lons = [x.lon for x in LevList]
    lats = [x.lat for x in LevList]

    for i in range(1, num_plots):
        data = [x.leveling[i] for x in LevList]
        label = dt.datetime.strftime(LevList[0].dtarray[i], "%Y-%m")
        single_panel_plot(axarr[idx2][idx1], lons, lats, data, vmin, vmax, label, 80)
        idx1 = idx1 + 1
        if idx1 == 5:   # 5 rows
            idx1 = 0
            idx2 = 1

    axarr[1][4].set_visible(False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='jet')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=8, fraction=0.5)
    cb.set_label('Cumulative Displacement (m)', fontsize=25)
    cb.ax.tick_params(labelsize=20)
    plt.savefig(outname)
    return


def single_panel_plot(ax, lons, lats, data, vmin, vmax, label, plotting_multiplier=200):
    """A plot within a multi-panel plot"""
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r')

    for i in range(len(data)):
        plotting_value = data[i]
        dot_color = custom_cmap.to_rgba(plotting_value)
        markersize = np.sqrt(abs(plotting_value))*plotting_multiplier
        # if markersize<4:
        # 	markersize=4
        if markersize > 35:
            markersize = 35
        ax.plot(lons[i], lats[i], marker='o', markersize=markersize, color=dot_color, fillstyle="full")
        ax.plot(lons[0], lats[0], marker='*', markersize=25, color='black', fillstyle="full")

    ax.xaxis.set_ticks(np.arange(-115.59, -115.48, 0.01))
    ax.tick_params(axis='x', labelsize=16, rotation=45)
    ax.yaxis.set_ticks(np.arange(33.00, 33.05, 0.01))
    ax.tick_params(axis='y', labelsize=16)
    ax.set_title(label, fontsize=35)
    return
