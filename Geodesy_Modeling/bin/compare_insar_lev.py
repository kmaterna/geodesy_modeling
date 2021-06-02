#!/usr/bin/env python3

# Compare InSAR Displacements with Leveling:
# TSX from 2012-2013
# S1 from 2014-2018
# S1 from 2018-2019
# S1 from OU/Cornell
# Individual UAVSAR intfs
# UAVSAR time series slices
# Holy Cow I reproduced Mariana's results.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import sys
from Geodesy_Modeling import multiSAR_utilities
from Geodesy_Modeling import InSAR_1D_Object
from Geodesy_Modeling import Leveling_Object
from Geodesy_Modeling import UAVSAR

def one_to_one_comparison(myLev, InSAR_Data, sat, filename, vmin=-50, vmax=50, gps_lon=None, gps_lat=None):
    """
    myLev : list of leveling objects (displacements) with two time intervals, a start and an end
    InSAR_Data : InSAR object (displacements) with a start and end time
    sat : string
    filename : string
    vmin/vmax : float, in mm
    right now, compares vertical leveling with LOS displacements (a little messy)
    """
    oto_lev, oto_tsx = [], [];
    lon_plotting, lat_plotting = [], [];

    lon_leveling_list = [item.lon for item in myLev];
    lat_leveling_list = [item.lat for item in myLev];
    vector_index, close_pixels = multiSAR_utilities.find_pixels_idxs_in_InSAR_Obj(InSAR_Data, lon_leveling_list,
                                                                                  lat_leveling_list);

    reference_insar_los = np.nanmean(np.array(InSAR_Data.LOS)[close_pixels[0]]);  # InSAR disp near leveling refpixel.
    # the first element of leveling is the datum Y-1225, so it should be used as reference for InSAR

    # Get the one-to-one pixels
    for i in range(len(lon_leveling_list)):
        if np.isnan(vector_index[i]):
            continue;
        else:
            leveling_disp = 1000 * (myLev[i].leveling[1] - myLev[i].leveling[0])  # negative sign convention
            insar_disp = np.nanmean(np.array(InSAR_Data.LOS)[close_pixels[i]]) - reference_insar_los;
            if ~np.isnan(leveling_disp) and ~np.isnan(insar_disp):
                oto_lev.append(leveling_disp);
                oto_tsx.append(insar_disp);
                lon_plotting.append(lon_leveling_list[i]);
                lat_plotting.append(lat_leveling_list[i]);

    fig, axarr = plt.subplots(2, 2, figsize=(14, 10));
    # Individual plot of leveling or InSAR
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    for i in range(len(oto_lev)):
        dot_color = custom_cmap.to_rgba(oto_lev[i]);
        dot_color_tsx = custom_cmap.to_rgba(oto_tsx[i]);
        axarr[0][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color,
                         fillstyle="full");
        axarr[1][0].plot(lon_plotting[i], lat_plotting[i], marker='o', markersize=10, color=dot_color_tsx,
                         fillstyle="full");

    axarr[0][0].set_title(
        "Leveling: " + dt.datetime.strftime(myLev[0].dtarray[0], "%m-%Y") + " to " + dt.datetime.strftime(
            myLev[0].dtarray[1], "%m-%Y"), fontsize=15);
    axarr[0][0].plot(myLev[0].lon, myLev[0].lat, '*', markersize=12, color='black');
    axarr[1][0].set_title(
        sat + ": " + dt.datetime.strftime(InSAR_Data.starttime, "%m-%Y") + " to " + dt.datetime.strftime(
            InSAR_Data.endtime, "%m-%Y"), fontsize=15);
    axarr[1][0].plot(myLev[0].lon, myLev[0].lat, '*', markersize=12, color='black');

    # The one-to-one plot
    axarr[0][1].plot([-80, 80], [-80, 80], linestyle='--', color='gray')
    axarr[0][1].plot(oto_lev, oto_tsx, markersize=10, marker='v', linewidth=0);
    axarr[0][1].set_xlabel('Leveling offset (mm)', fontsize=20);
    axarr[0][1].set_ylabel(sat + ' offset (mm)', fontsize=20);
    axarr[0][1].tick_params(axis='both', which='major', labelsize=16)
    axarr[0][1].set_xlim([-80, 80])
    axarr[0][1].set_ylim([-80, 80])
    axarr[0][1].grid(True)

    # Plotting the InSAR data in the bottom panel, as used in other panels
    plotting_data = np.subtract(InSAR_Data.LOS, reference_insar_los);
    axarr[1][1].scatter(InSAR_Data.lon, InSAR_Data.lat, c=plotting_data, s=8,
                        marker='o', cmap='RdYlBu_r', vmin=vmin, vmax=vmax);
    axarr[1][1].plot(lon_leveling_list, lat_leveling_list, '*', color='black');
    axarr[1][1].plot(myLev[0].lon, myLev[0].lat, '*', color='red');
    axarr[1][1].plot(-115.510, 33.081, 'v', markersize=10, color='black');
    if gps_lon:
        for i in range(len(gps_lon)):
            axarr[1][1].plot(gps_lon[i], gps_lat[i], 'v', markersize=10, color='black');

    _ = fig.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);
    plt.savefig(filename);
    print("Saving %s " % filename);
    return;


def drive_ou_cornell_comparison(myLev, data_file, los_file, s1_slice, lev_slice, outfile):
    InSAR_Data = InSAR_1D_Object.inputs.inputs_cornell_ou_velocities_hdf5(data_file, los_file, s1_slice);
    InSAR_Data = InSAR_1D_Object.utilities.remove_nans(InSAR_Data);
    myLev = Leveling_Object.utilities.get_onetime_displacements(myLev, lev_slice[0], lev_slice[1]);  # one lev slice
    one_to_one_comparison(myLev, InSAR_Data, "S1", outfile);
    return;


def drive_single_uavsar_intf_comparison(myLev, bounds, uavsar_filename, los_filename, lev_slice, outfile):
    """Read the UAVSAR Data"""
    InSAR_Data = InSAR_1D_Object.inputs.inputs_isce_unw_geo_losrdr(uavsar_filename, los_filename, bounds);
    InSAR_Data = InSAR_1D_Object.utilities.remove_nans(InSAR_Data);
    InSAR_Data = InSAR_1D_Object.utilities.flip_los_sign(InSAR_Data);
    myLev = Leveling_Object.utilities.get_onetime_displacements(myLev, lev_slice[0], lev_slice[1]);  # one lev slice
    one_to_one_comparison(myLev, InSAR_Data, "UAV", outfile);
    return;


def drive_tre_comparison(myLev, los_filename, lev_slice, outfile):
    """Read TRE data and compare with leveling"""
    VertTSXData, EastTSXData = InSAR_1D_Object.inputs.inputs_TRE_vert_east(los_filename);
    myLev = Leveling_Object.utilities.get_onetime_displacements(myLev, lev_slice[0], lev_slice[1]);  # one lev slice
    one_to_one_comparison(myLev, VertTSXData, "TSX", outfile);
    InSAR_1D_Object.outputs.plot_insar(VertTSXData, outfile + "InSAR_velo.png");
    return;


def drive_uavsar_ts_comparison(myLev, losfile, lonfile, latfile, gps_lon, gps_lat, lev_slice, uav_slice, outfile):
    """A different type of UAVSAR format, the SBAS time series analysis"""
    myUAVSAR_TS = UAVSAR.uavsar_readwrite.inputs_TS_grd(losfile, lonfile, latfile);
    myLev = Leveling_Object.utilities.get_onetime_displacements(myLev, lev_slice[0], lev_slice[1]);  # one lev slice
    myUAVSAR_insarobj = UAVSAR.utilities.get_onetime_displacements(myUAVSAR_TS, uav_slice[0], uav_slice[1]);
    myUAVSAR_insarobj = InSAR_1D_Object.utilities.flip_los_sign(myUAVSAR_insarobj);
    one_to_one_comparison(myLev, myUAVSAR_insarobj, "UAVSAR", outfile, gps_lon=gps_lon, gps_lat=gps_lat);
    return;


if __name__ == "__main__":
    # Opening stuff
    config_filename = sys.argv[1];
    file_dict = multiSAR_utilities.get_file_dictionary(config_filename);
    myLev = Leveling_Object.leveling_inputs.inputs_brawley_leveling(file_dict["leveling"].split()[0],
                                                                    file_dict["leveling"].split()[1]);
    myLev = Leveling_Object.leveling_inputs.compute_rel_to_datum_nov_2009(myLev);

    # # # S1_Cornell experiment (2015 data), leveling slice 5-6
    # output_dir = "S1_OU/T4D/";
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice=0, lev_slice=[5, 6], outfile=output_dir+'ascending_56.png');
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice=0, lev_slice=[5, 6], outfile=output_dir+'descending_56.png');
    #
    # # # S1_Cornell experiment (2016 data), leveling slice 6-7
    # output_dir = "S1_OU/T4E/";
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice=1, lev_slice=[6, 7], outfile=output_dir+'ascending_67.png');
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice=1, lev_slice=[6, 7], outfile=output_dir+'descending_67.png');
    #
    # # # S1_Cornell experiment (2017 data), leveling slice 7-8
    # output_dir = "S1_OU/T4F/";
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice=2, lev_slice=[7, 8], outfile=output_dir+'ascending_78.png');
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice=2, lev_slice=[7, 8], outfile=output_dir+'descending_78.png');
    #
    # # S1_Cornell experiment (2018 data), leveling slice 8-9
    # output_dir = "S1_OU/T5/";
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice=3, lev_slice=[8, 9], outfile=output_dir+'ascending_89.png');
    # drive_ou_cornell_comparison(myLev, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice=3, lev_slice=[8, 9], outfile=output_dir+'descending_89.png');

    # TSX experiment: 2012-2013, leveling slice 3-4
    output_dir = "TSX/"
    drive_tre_comparison(myLev, file_dict["tsx"], lev_slice=[3, 4], outfile=output_dir+"one_to_one_34.png");

    # # S1 experiment: 2014-2018, leveling slice 5-8
    # output_dir = "SNT1/"
    # drive_tre_comparison(myLev, file_dict["snt1"], lev_slice=[5, 8], outfile=output_dir+"one_to_one_58.png");
    #
    # # S1 experiment: 2018-2019, leveling slice 8-9
    # output_dir = "SNT2/"
    # drive_tre_comparison(myLev, file_dict["snt2"], lev_slice=[8, 9], outfile=output_dir+"one_to_one_89.png");

    # # Individual UAVSAR experiments, starting with 2011-2012, leveling slice 2-3
    # # Should set vmin/vmax to -150/150 for this one.
    # output_dir = "UAVSAR_intfs/"
    # bounds = (dt.datetime.strptime("20111110", "%Y%m%d"), dt.datetime.strptime("20120926", "%Y%m%d"));
    # drive_single_uavsar_intf_comparison(myLev, bounds, file_dict["uavsar_08508_2011_2012_unw"],
    #                                     file_dict["uavsar_08508_2011_2012_los"], lev_slice=[2, 3],
    #                                     outfile=output_dir+"one_to_one_23.png");
    #
    # # Individual UAVSAR experiments, 2010-2011, leveling slice 1-2
    # output_dir = "UAVSAR_intfs/"
    # bounds = (dt.datetime.strptime("20101215", "%Y%m%d"), dt.datetime.strptime("20111110", "%Y%m%d"));
    # drive_single_uavsar_intf_comparison(myLev, bounds, file_dict["uavsar_08508_2010_2011_unw"],
    #                                     file_dict["uavsar_08508_2010_2011_los"], lev_slice=[1, 2],
    #                                     outfile=output_dir+"one_to_one_12.png");
    #
    # # Individual UAVSAR experiments, 2009-2010, leveling slice 0-1
    # output_dir = "UAVSAR_intfs/"
    # bounds = (dt.datetime.strptime("20091015", "%Y%m%d"), dt.datetime.strptime("20101215", "%Y%m%d"));
    # drive_single_uavsar_intf_comparison(myLev, bounds, file_dict["uavsar_08508_2009_2010_unw"],
    #                                     file_dict["uavsar_08508_2009_2010_los"], lev_slice=[0, 1],
    #                                     outfile=output_dir+"one_to_one_01.png");
    #
    # outdir = "UAVSAR_Apr29/";
    # gps_lons = [-115.510, -115.628392, -115.581895, -115.613, -115.735];  # Stations P506, P495, WMDG, WMCA and CRRS
    # gps_lats = [33.081, 33.044960, 33.038325, 33.072, 33.070];
    # uav_los, uav_lon, uav_lat = file_dict["uavsar_file"], file_dict["uavsar_lon"], file_dict["uavsar_lat"];
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[0, 1], uav_slice=[1, 4],
    #                            outfile=outdir+"one_01_14.png");  # 2009-2011
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[1, 2], uav_slice=[4, 6],
    #                            outfile=outdir+"one_12_46.png");  # 2011 to 2011
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[2, 3], uav_slice=[6, 7],
    #                            outfile=outdir+"one_23_67.png");  # 2011 to 2012
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[3, 4], uav_slice=[7, 8],
    #                            outfile=outdir+"one_34_78.png");  # 2012 to 2013
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[4, 5], uav_slice=[8, 9],
    #                            outfile=outdir+"one_45_89.png");  # 2013 to 2014
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[3, 5], uav_slice=[7, 9],
    #                            outfile=outdir+"one_35_79.png");  # 2012 to 2014
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[5, 8],
    #                            uav_slice=[9, 10], outfile=outdir+"one_58_910.png");  # 2014 to 2017
    # drive_uavsar_ts_comparison(myLev, uav_los, uav_lon, uav_lat, gps_lons, gps_lats, lev_slice=[3, 8],
    #                            uav_slice=[7, 10], outfile=outdir+"one_38_710.png");  # 2012 to 2014
