#!/usr/bin/env python3

"""
Driver for Brawley Project, 2020-2021
Compare InSAR Displacements with Leveling:
TSX from 2012-2013
S1 from 2014-2018
S1 from 2018-2019
S1 from OU/Cornell
Individual UAVSAR intfs
UAVSAR time series slices
Holy Cow I reproduced Mariana's results.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import sys
from Geodesy_Modeling.src import general_utils, InSAR_1D_Object, Leveling_Object, UAVSAR
from Tectonic_Utils.read_write import general_io


def compute_difference_metrics_on_same_pixels(list1, list2):
    """
    :param list1: a 1d array of LOS data from platform 1 (like Leveling)
    :param list2: a matching 1d array of LOS data from platform 2 (like UAVSAR)
    :returns: average misfit value, and r^2 coefficient.
    """
    misfit_metric = np.nanmean(np.abs(np.subtract(list1, list2)));  # average deviation from 1-to-1
    corr_matrix = np.corrcoef(list1, list2)
    corr = corr_matrix[0, 1]
    r2 = corr ** 2
    return misfit_metric, r2;


def get_nearest_in_pixel_list(tuple_coord_list, target_lons, target_lats):
    """
    Get the nearest index (and neighbors) for each given coordinate.

    :param tuple_coord_list: list of matching tuples of (lon, lat)
    :param target_lons: 1d array of longitudes to be found
    :param target_lats: 1d array of latitudes to be found
    :returns closest_index: a list that matches the length of target_lons
    :returns close_indices: a list of lists (might return the 100 closest pixels for a given coordinate)
    """
    print("Finding target leveling pixels in vector of data");
    closest_index, close_indices = [], [];
    for tlon, tlat in zip(target_lons, target_lats):
        i_found, _, closer_pts = general_utils.get_nearest_pixels_in_list(tuple_coord_list, tlon, tlat);
        if i_found != -1:
            closest_index.append(i_found);
            close_indices.append(closer_pts);
        else:
            closest_index.append(np.nan);
            close_indices.append(np.nan);
    return closest_index, close_indices;


def one_to_one_comparison(myLev, InSAR_Data, sat, filename, vmin=-50, vmax=50, gps_lon=None, gps_lat=None,
                          gps_names=None, graph_scale=80, label="LOS", proj_vertical=0, lkv=(0, 0, 1)):
    """
    myLev : list of leveling objects (displacements) with two time intervals, a start and an end
    InSAR_Data : InSAR object (displacements) with a start and end time
    sat : string
    filename : string
    vmin/vmax : float, in mm
    Can do LOS or projected pseudo-vertical, depending on experimental step
    """
    if proj_vertical:  # experimental step
        proj_InSAR_Data = InSAR_1D_Object.utilities.proj_los_into_vertical_no_horiz(InSAR_Data, const_lkv=lkv);
    else:
        proj_InSAR_Data = InSAR_Data;

    lon_plotting, lat_plotting = [], [];
    lon_leveling_list = [item.lon for item in myLev];
    lat_leveling_list = [item.lat for item in myLev];
    vector_index, close_pixels = get_nearest_in_pixel_list(InSAR_Data.get_coordinate_tuples(),
                                                           lon_leveling_list, lat_leveling_list);

    reference_insar_los = np.nanmean(np.array(proj_InSAR_Data.LOS)[close_pixels[0]]);  # InSAR disp near lev. refpixel.
    # the first element of leveling is the datum Y-1225, so it should be used as reference for InSAR

    # Get the one-to-one pixels
    oto_lev, oto_tsx = [], [];
    for i in range(len(lon_leveling_list)):
        if np.isnan(vector_index[i]):
            continue;
        else:
            leveling_disp = 1000 * (myLev[i].leveling[1] - myLev[i].leveling[0])  # negative sign convention
            insar_disp = np.nanmean(np.array(proj_InSAR_Data.LOS)[close_pixels[i]]) - reference_insar_los;
            if ~np.isnan(leveling_disp) and ~np.isnan(insar_disp):
                oto_lev.append(leveling_disp);
                oto_tsx.append(insar_disp);
                lon_plotting.append(lon_leveling_list[i]);
                lat_plotting.append(lat_leveling_list[i]);
    # Computing misfit
    misfit_metric, r2 = compute_difference_metrics_on_same_pixels(oto_lev, oto_tsx);

    # Comparison plot between leveling and InSAR
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

    leveling_title = "Leveling: " + dt.datetime.strftime(myLev[0].dtarray[0], "%m-%Y") + " to " + dt.datetime.strftime(
            myLev[0].dtarray[1], "%m-%Y");
    axarr[0][0].set_title(leveling_title, fontsize=15);
    axarr[0][0].plot(myLev[0].lon, myLev[0].lat, '*', markersize=12, color='black');
    axarr[0][0].set_xticks([-115.59, -115.57, -115.55, -115.53, -115.51]);
    axarr[0][0].ticklabel_format(useOffset=False)

    satellite_title = sat + ": " + dt.datetime.strftime(InSAR_Data.starttime, "%m-%Y") + " to " + dt.datetime.strftime(
            InSAR_Data.endtime, "%m-%Y")
    axarr[1][0].set_title(satellite_title, fontsize=15);
    axarr[1][0].plot(myLev[0].lon, myLev[0].lat, '*', markersize=12, color='black');
    axarr[1][0].set_xticks([-115.59, -115.57, -115.55, -115.53, -115.51]);
    axarr[1][0].ticklabel_format(useOffset=False)

    # The one-to-one plot
    axarr[0][1].plot([-80, 80], [-80, 80], linestyle='--', color='gray')
    axarr[0][1].plot(oto_lev, oto_tsx, markersize=10, marker='v', linewidth=0);
    axarr[0][1].set_xlabel('Leveling offset (mm)', fontsize=20);
    axarr[0][1].set_ylabel(sat + ' offset (mm)', fontsize=20);
    axarr[0][1].tick_params(axis='both', which='major', labelsize=16)
    axarr[0][1].set_xlim([-graph_scale, graph_scale])
    axarr[0][1].set_ylim([-graph_scale, graph_scale])
    axarr[0][1].set_title("Avg misfit = %.2f mm, Rsq = %.2f" % (misfit_metric, r2), fontsize=15);
    axarr[0][1].grid(True)

    # Plotting the InSAR data in the bottom panel, as used in other panels
    plotting_data = np.subtract(InSAR_Data.LOS, reference_insar_los);
    axarr[1][1].scatter(InSAR_Data.lon, InSAR_Data.lat, c=plotting_data, s=8,
                        marker='o', cmap='RdYlBu_r', vmin=vmin, vmax=vmax);
    axarr[1][1].plot(lon_leveling_list, lat_leveling_list, '*', color='black');
    axarr[1][1].plot(myLev[0].lon, myLev[0].lat, '*', color='red');
    axarr[1][1].plot(-115.510, 33.081, 'v', markersize=10, color='black');
    axarr[1][1].text(-115.510, 33.081, '  P506', color='black');
    if gps_lon is not None:
        for i in range(len(gps_lon)):
            axarr[1][1].plot(gps_lon[i], gps_lat[i], 'v', markersize=10, color='black');
            axarr[1][1].text(gps_lon[i], gps_lat[i], '  ' + gps_names[i], color='black');

    _ = fig.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label(label + ' Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);
    plt.savefig(filename);
    print("Saving %s " % filename);

    # Writing file out
    filename_txt = filename.split('.')[0]+'.txt'
    print("Writing %s" % filename_txt);
    ofile = open(filename_txt, 'w');
    ofile.write("# leveling los / "+leveling_title + ' / ' + satellite_title + '\n');
    for i in range(len(oto_lev)):
        ofile.write("%f %f \n" % (oto_lev[i], oto_tsx[i]) );
    ofile.close();
    return;


def read_leveling_data(data_file, error_file):
    myLev = Leveling_Object.leveling_inputs.inputs_brawley_leveling(data_file, error_file);
    myLev = Leveling_Object.leveling_inputs.compute_rel_to_datum_nov_2009(myLev);
    return myLev;


def drive_ou_cornell_compare(file_dict, s1_slice, lev_slice, track, outfile):
    if track == 'ascending' or track == 'asc':
        data_file, los_file = file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"];
        lkv = np.array([float(x) for x in file_dict["s1_ascending_lkv"].split('/')]);
    else:
        data_file, los_file = file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"];
        lkv = np.array([float(x) for x in file_dict["s1_descending_lkv"].split('/')]);
    myLev = read_leveling_data(file_dict["leveling"], file_dict["lev_error"]);
    InSAR_Data = InSAR_1D_Object.inputs.inputs_cornell_ou_velocities_hdf5(data_file, los_file, s1_slice);
    InSAR_Data = InSAR_Data.remove_nans();
    myLev = Leveling_Object.LevStation.get_onetime_displacements_list(myLev, lev_slice[0], lev_slice[1]);  # one slice
    one_to_one_comparison(myLev, InSAR_Data, "S1", outfile, label="LOS", graph_scale=50, proj_vertical=1, lkv=lkv);
    return;


def drive_single_uavsar_intf_compare(file_dict, bounds, uavsar_filename, los_filename, lev_slice, outfile):
    """Read the UAVSAR Data"""
    myLev = read_leveling_data(file_dict["leveling"], file_dict["lev_error"]);
    lkv = np.array([float(x) for x in file_dict["uavsar_08508_lkv"].split('/')]);
    InSAR_Data = InSAR_1D_Object.inputs.inputs_isce_unw_geo_losrdr(uavsar_filename, los_filename, starttime=bounds[0],
                                                                   endtime=bounds[1]);
    InSAR_Data = InSAR_Data.remove_nans();
    InSAR_Data = InSAR_Data.flip_los_sign();
    InSAR_Data = InSAR_1D_Object.remove_ramp.remove_ramp(InSAR_Data);  # experimental step
    myLev = Leveling_Object.LevStation.get_onetime_displacements_list(myLev, lev_slice[0], lev_slice[1]);  # one slice
    one_to_one_comparison(myLev, InSAR_Data, "UAV", outfile, label="LOS", proj_vertical=1, lkv=lkv);
    return;


def drive_tre_compare(file_dict, insar_key, lev_slice, outdir, outfile):
    """Read TRE data and compare with leveling. Don't need to project into vertical."""
    insar_datafile = file_dict[insar_key];
    myLev = read_leveling_data(file_dict["leveling"], file_dict["lev_error"]);
    VertTSXData, EastTSXData = InSAR_1D_Object.inputs.inputs_TRE_vert_east(insar_datafile);
    myLev = Leveling_Object.LevStation.get_onetime_displacements_list(myLev, lev_slice[0], lev_slice[1]);  # one slice
    one_to_one_comparison(myLev, VertTSXData, insar_key, outdir + outfile, label="Vertical", graph_scale=50);
    fields = general_io.read_gmt_multisegment_latlon(file_dict["field_file"], split_delimiter=',');
    InSAR_1D_Object.outputs.plot_insar(VertTSXData, outdir+"InSAR_velo.png", lons_annot=fields[0][0],
                                       lats_annot=fields[1][0]);
    return;


def drive_uavsar_ts_compare(file_dict, lev_slice, uav_slice, outfile):
    """A different type of UAVSAR format, the SBAS time series analysis"""
    losfile, lonfile, latfile = file_dict["uavsar_file"], file_dict["uavsar_lon"], file_dict["uavsar_lat"];
    myUAVSAR_TS = UAVSAR.uavsar_readwrite.inputs_TS_grd(losfile, lonfile, latfile);
    myLev = read_leveling_data(file_dict["leveling"], file_dict["lev_error"]);
    lkv = np.array([float(x) for x in file_dict["uavsar_26509_lkv"].split('/')]);
    gps_lon, gps_lat, gps_names = np.loadtxt(file_dict["gnss_file"], dtype={'names': ('lon', 'lat', 'name'),
                                                                            'formats': (float, float, 'U4')},
                                             skiprows=1, unpack=True);
    myLev = Leveling_Object.LevStation.get_onetime_displacements_list(myLev, lev_slice[0], lev_slice[1]);  # one slice
    myUAVSAR_insarobj = UAVSAR.utilities.get_onetime_displacements(myUAVSAR_TS, uav_slice[0], uav_slice[1]);
    myUAVSAR_insarobj = myUAVSAR_insarobj.flip_los_sign();
    myUAVSAR_insarobj = InSAR_1D_Object.remove_ramp.remove_ramp(myUAVSAR_insarobj);  # experimental step
    one_to_one_comparison(myLev, myUAVSAR_insarobj, "UAVSAR", outfile, gps_lon=gps_lon, gps_lat=gps_lat,
                          gps_names=gps_names, label="LOS", proj_vertical=1, lkv=lkv);
    return;


if __name__ == "__main__":
    config_filename = sys.argv[1];
    file_dict = general_utils.get_file_dictionary(config_filename);

    # # TSX experiment: 2012-2013, leveling slice 3-4
    drive_tre_compare(file_dict, 'tsx', lev_slice=[3, 4], outdir="TSX/", outfile="one_to_one_34.png");
    # # S1 experiment: 2014-2018, leveling slice 5-8
    # drive_tre_compare(file_dict, "snt1", lev_slice=[5, 8], outdir="SNT1/", outfile="vert_58.png");
    # # S1 experiment: 2018-2019, leveling slice 8-9
    # drive_tre_compare(file_dict, "snt2", lev_slice=[8, 9], outdir="SNT2/", outfile="vert_89.png");

    # # # S1_Cornell experiment (2015 data), leveling slice 5-6
    drive_ou_cornell_compare(file_dict, s1_slice=0, lev_slice=[5, 6], track='asc', outfile='S1_OU/T4D/asc_56.png');
    drive_ou_cornell_compare(file_dict, s1_slice=0, lev_slice=[5, 6], track='desc', outfile='S1_OU/T4D/desc_56.png');

    # # # S1_Cornell experiment (2016 data), leveling slice 6-7
    drive_ou_cornell_compare(file_dict, s1_slice=1, lev_slice=[6, 7], track='asc', outfile='S1_OU/T4E/asc_67.png');
    drive_ou_cornell_compare(file_dict, s1_slice=1, lev_slice=[6, 7], track='desc', outfile='S1_OU/T4E/desc_67.png');
    #
    # # # S1_Cornell experiment (2017 data), leveling slice 7-8
    drive_ou_cornell_compare(file_dict, s1_slice=2, lev_slice=[7, 8], track='asc', outfile='S1_OU/T4F/asc_78.png');
    drive_ou_cornell_compare(file_dict, s1_slice=2, lev_slice=[7, 8], track='desc', outfile='S1_OU/T4F/desc_78.png');
    #
    # S1_Cornell experiment (2018 data), leveling slice 8-9
    drive_ou_cornell_compare(file_dict, s1_slice=3, lev_slice=[8, 9], track='asc', outfile='S1_OU/T5/asc_89.png');
    drive_ou_cornell_compare(file_dict, s1_slice=3, lev_slice=[8, 9], track='desc', outfile='S1_OU/T5/desc_89.png');

    # Individual UAVSAR experiments, starting with 2011-2012, leveling slice 2-3
    # Should set vmin/vmax to -150/150 for this one.
    output_dir = "UAVSAR_intfs/"
    bounds = (dt.datetime.strptime("20111110", "%Y%m%d"), dt.datetime.strptime("20120926", "%Y%m%d"));
    drive_single_uavsar_intf_compare(file_dict, bounds, file_dict["uavsar_08508_2011_2012_unw"],
                                     file_dict["uavsar_08508_2011_2012_los"], lev_slice=[2, 3],
                                     outfile=output_dir+"one_to_one_23.png");
    # Individual UAVSAR experiments, 2010-2011, leveling slice 1-2
    bounds = (dt.datetime.strptime("20101215", "%Y%m%d"), dt.datetime.strptime("20111110", "%Y%m%d"));
    drive_single_uavsar_intf_compare(file_dict, bounds, file_dict["uavsar_08508_2010_2011_unw"],
                                     file_dict["uavsar_08508_2010_2011_los"], lev_slice=[1, 2],
                                     outfile=output_dir+"one_to_one_12.png");
    # Individual UAVSAR experiments, 2009-2010, leveling slice 0-1
    bounds = (dt.datetime.strptime("20091015", "%Y%m%d"), dt.datetime.strptime("20101215", "%Y%m%d"));
    drive_single_uavsar_intf_compare(file_dict, bounds, file_dict["uavsar_08508_2009_2010_unw"],
                                     file_dict["uavsar_08508_2009_2010_los"], lev_slice=[0, 1],
                                     outfile=output_dir+"one_to_one_01.png");

    outdir = "UAVSAR_Apr29/";
    drive_uavsar_ts_compare(file_dict, lev_slice=[0, 1], uav_slice=[1, 4], outfile=outdir+"o_01_14.png");  # 2009-11
    drive_uavsar_ts_compare(file_dict, lev_slice=[1, 2], uav_slice=[4, 6], outfile=outdir+"o_12_46.png");  # 2011-11
    drive_uavsar_ts_compare(file_dict, lev_slice=[2, 3], uav_slice=[6, 7], outfile=outdir+"o_23_67.png");  # 2011-12
    drive_uavsar_ts_compare(file_dict, lev_slice=[3, 4], uav_slice=[7, 8], outfile=outdir+"o_34_78.png");  # 2012-13
    drive_uavsar_ts_compare(file_dict, lev_slice=[4, 5], uav_slice=[8, 9], outfile=outdir+"o_45_89.png");  # 2013-14
    drive_uavsar_ts_compare(file_dict, lev_slice=[3, 5], uav_slice=[7, 9], outfile=outdir+"o_35_79.png");  # 2012-14
    drive_uavsar_ts_compare(file_dict, lev_slice=[5, 8], uav_slice=[9, 10], outfile=outdir+"o_58_910.png");  # 2014-17
    drive_uavsar_ts_compare(file_dict, lev_slice=[3, 8], uav_slice=[7, 10], outfile=outdir+"o_38_710.png");  # 2012-14
