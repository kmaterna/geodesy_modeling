# THIS FILE IS ON THE CHOPPING BLOCK.

# April 1, 2020
# Tasks: 
# Read in leveling and UAVSAR datasets
# Associate long/lat with UAVSAR
# Find a pixel that represents the leveling datum
# Find pixels that represent the other benchmarks too
# Plot a few one-to-one comparisons: 
#    - 2009-2011 is close enough
#    - 2011 to 2011 is good
#    - Coseismic is good
# Plot the GPS on top of UAVSAR, and pull out those pixels. 
# 
# Bottom line: 
# It looks like our noise level is about 2.5 cm
# Above that, we can see many consistent features between the two data sets. 
# Subsidence and uplift lobes look largely consistent. 
# In 2009-2011, there is difference in where the subsidence happens and magnitude of subsidence
# I think this might be related to unwrapping issues near the edge of the UAVSAR scene? 
# It looks like there's localized subsidence between June 2014 and October 2014. 
# There is a difference in the 2012 Brawley Earthquakes for track 26509

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import Geodesy_Modeling.UAVSAR.uavsar_readwrite


def avg_uavsar_disp(TS, slicenum, row, col):
    """Average around a few pixels"""
    width_pixels = 10;
    return np.nanmean(TS[slicenum, row - width_pixels:row + width_pixels, col - width_pixels:col + width_pixels]);


def one_to_one_comparison(myLev, myUAVSAR, row, col, lev1, lev2, uav1, uav2, outdir):
    filename = outdir + "/one_to_one_" + str(lev1) + str(lev2) + "_" + str(uav1) + str(uav2);

    oto_lev, oto_uavsar = [], [];
    lon_plotting, lat_plotting = [], [];

    # With respect to the datum
    UAVSAR_general_scene = np.subtract(myUAVSAR.TS[uav2, :, :], myUAVSAR.TS[uav1, :, :]);
    UAVSAR_general_scene = np.subtract(UAVSAR_general_scene, UAVSAR_general_scene[row[0]][col[0]]);

    # How much the reference leveling benchmark moved.
    UAVSAR_ref_offset = avg_uavsar_disp(myUAVSAR.TS, uav2, row[0], col[0]) - avg_uavsar_disp(myUAVSAR.TS, uav1, row[0],
                                                                                             col[0])

    for i in range(len(myLev.lon)):
        if row[i] == np.nan:
            continue;
        else:
            leveling_offset = 1000 * (myLev.leveling[i][lev2] - myLev.leveling[i][lev1])  # negative sign convention
            uavsar_offset = avg_uavsar_disp(myUAVSAR.TS, uav2, row[i], col[i]) - \
                            avg_uavsar_disp(myUAVSAR.TS, uav1, row[i], col[i]) - UAVSAR_ref_offset;
            uavsar_offset = -1 * uavsar_offset;
            if ~np.isnan(leveling_offset) and ~np.isnan(uavsar_offset):
                oto_lev.append(leveling_offset);
                oto_uavsar.append(uavsar_offset);
                lon_plotting.append(myLev.lon[i]);
                lat_plotting.append(myLev.lat[i]);

    print("one-to-one: ", np.shape(oto_lev), np.shape(oto_uavsar));
    vmin = -50;
    vmax = 50;

    fig, axarr = plt.subplots(2, 2, figsize=(14, 10));

    # Individual plot of leveling or UAVSAR
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    for i in range(len(oto_lev)):
        dot_color = custom_cmap.to_rgba(oto_lev[i]);
        dot_color_uavsar = custom_cmap.to_rgba(oto_uavsar[i]);
        axarr[0][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color,
                         fillstyle="full");
        axarr[1][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color_uavsar,
                         fillstyle="full");

    axarr[0][0].set_title(
        "Leveling: " + dt.datetime.strftime(myLev.dtarray[lev1], "%m-%Y") + " to " + dt.datetime.strftime(
            myLev.dtarray[lev2], "%m-%Y"), fontsize=15);
    axarr[0][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12, color='black');
    axarr[1][0].set_title(
        "UAVSAR: " + dt.datetime.strftime(myUAVSAR.dtarray[uav1], "%m-%Y") + " to " + dt.datetime.strftime(
            myUAVSAR.dtarray[uav2], "%m-%Y"), fontsize=15);
    axarr[1][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12, color='black');

    # The one-to-one plot
    axarr[0][1].plot([-80, 80], [-80, 80], linestyle='--', color='gray')
    axarr[0][1].plot(oto_lev, oto_uavsar, markersize=10, marker='v', linewidth=0);
    axarr[0][1].set_xlabel('Leveling offset (mm)', fontsize=20);
    axarr[0][1].set_ylabel('UAVSAR offset (mm)', fontsize=20);
    axarr[0][1].tick_params(axis='both', which='major', labelsize=16)
    axarr[0][1].set_xlim([-80, 80])
    axarr[0][1].set_ylim([-80, 80])
    axarr[0][1].grid(True)

    axarr[1][1].scatter(myUAVSAR.lon, myUAVSAR.lat, c=UAVSAR_general_scene, s=8,
                        marker='o', cmap='RdYlBu', vmin=vmin, vmax=vmax);
    axarr[1][1].plot(myLev.lon, myLev.lat, '*', color='black');
    axarr[1][1].plot(myLev.lon[0], myLev.lat[0], '*', color='red');
    axarr[1][1].plot(-115.510, 33.081, 'v', markersize=10, color='black');  # P506
    axarr[1][1].plot(-115.628392, 33.044960, 'v', markersize=10, color='black');  # P495
    axarr[1][1].plot(-115.581895, 33.038325, 'v', markersize=10, color='black');  # WMDG
    axarr[1][1].plot(-115.735041, 33.069808, 'v', markersize=10, color='black');  # CRRS
    axarr[1][1].plot(-115.613, 33.072, 'v', markersize=10, color='black');  # WMCA

    cbarax = fig.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(filename + ".png");

    return;


def leveling_uavsar_ts_comparison_driver():
    """
    If you want to use the UAVSAR time series format (e.g., a track processed in ISCE)
    Use this to compare with leveling
    """
    outdir = "UAVSAR_Apr29";
    myUAVSAR = Geodesy_Modeling.UAVSAR.uavsar_readwrite.inputs_TS_grd(file_dict["uavsar_file"], file_dict["uavsar_lon"],
                                                                      file_dict["uavsar_lat"]);
    # Comparing UAVSAR with leveling.
    find_leveling_in_uavsar(myLev, myUAVSAR);  # only have to do the first time.
    row, col = np.loadtxt("leveling_index_cache.txt", usecols=(3, 4), skiprows=1, unpack=True,
                          dtype={'names': ('row', 'col'), 'formats': (int, int)});
    one_to_one_comparison(myLev, myUAVSAR, row, col, 0, 1, 1, 4, outdir);  # 2009 to 2011
    one_to_one_comparison(myLev, myUAVSAR, row, col, 1, 2, 4, 6, outdir);  # 2011 to 2011
    one_to_one_comparison(myLev, myUAVSAR, row, col, 2, 3, 6, 7, outdir);  # Brawley coseismic
    one_to_one_comparison(myLev, myUAVSAR, row, col, 3, 4, 7, 8, outdir);  # 2012 to 2013
    one_to_one_comparison(myLev, myUAVSAR, row, col, 4, 5, 8, 9, outdir);  # 2013 to 2014
    one_to_one_comparison(myLev, myUAVSAR, row, col, 3, 5, 7, 9, outdir);  # 2012 to 2014
    one_to_one_comparison(myLev, myUAVSAR, row, col, 5, 8, 9, 10, outdir);  # 2014 to 2017
    one_to_one_comparison(myLev, myUAVSAR, row, col, 3, 8, 7, 10, outdir);  # 2012 to 2017
    return;
