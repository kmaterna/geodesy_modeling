# Compare InSAR Displacements with Leveling:
# TSX from 2012-2013
# S1 from 2014-2018
# S1 from 2018-2019
# S1 from OU/Cornell
# Individual UAVSAR intfs
# Holy Cow I reproduced Mariana's results. 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import datetime as dt
import sys
import multiSAR_utilities
import multiSAR_input_functions
import InSAR_Object.inputs
import InSAR_Object.utilities


def plot_InSAR(InSAR_Data, output_dir):
    plt.figure(figsize=(10, 10));
    plt.scatter(InSAR_Data.lon, InSAR_Data.lat, c=InSAR_Data.LOS, s=10);
    plt.colorbar();
    plt.savefig(output_dir + "InSAR_vel.png");
    return;


def one_to_one_comparison(myLev, InSAR_Data, vector_index, lev1, lev2, sat, output_dir, flip_insar=False):
    filename = output_dir + "one_to_one_" + str(lev1) + str(lev2);

    oto_lev, oto_tsx = [], [];
    lon_plotting, lat_plotting = [], [];

    reference_insar_los = InSAR_Data.LOS[vector_index[0]];
    print(reference_insar_los);
    # the first line of leveling is the datum Y-1225, so it should be used as reference for InSAR

    # Get the one-to-one pixels
    for i in range(len(myLev.lon)):
        if vector_index[i] == np.nan:
            continue;
        else:
            leveling_disp = 1000 * (myLev.leveling[i][lev2] - myLev.leveling[i][lev1])  # negative sign convention
            insar_disp = InSAR_Data.LOS[vector_index[i]] - reference_insar_los;
            if flip_insar:
                insar_disp = insar_disp * -1;
            if ~np.isnan(leveling_disp) and ~np.isnan(insar_disp):
                oto_lev.append(leveling_disp);
                oto_tsx.append(insar_disp);
                lon_plotting.append(myLev.lon[i]);
                lat_plotting.append(myLev.lat[i]);

    vmin = -50;
    vmax = 50;
    fig, axarr = plt.subplots(2, 2, figsize=(14, 10));

    # Individual plot of leveling or TSX
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
        "Leveling: " + dt.datetime.strftime(myLev.dtarray[lev1], "%m-%Y") + " to " + dt.datetime.strftime(
            myLev.dtarray[lev2], "%m-%Y"), fontsize=15);
    axarr[0][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12, color='black');
    axarr[1][0].set_title(
        sat + ": " + dt.datetime.strftime(InSAR_Data.starttime, "%m-%Y") + " to " + dt.datetime.strftime(
            InSAR_Data.endtime, "%m-%Y"), fontsize=15);
    axarr[1][0].plot(myLev.lon[0], myLev.lat[0], '*', markersize=12, color='black');

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
    if flip_insar:
        plotting_data = plotting_data * -1;
    axarr[1][1].scatter(InSAR_Data.lon, InSAR_Data.lat, c=plotting_data, s=8,
                        marker='o', cmap='RdYlBu_r', vmin=vmin, vmax=vmax);
    axarr[1][1].plot(myLev.lon, myLev.lat, '*', color='black');
    axarr[1][1].plot(myLev.lon[0], myLev.lat[0], '*', color='red');
    axarr[1][1].plot(-115.510, 33.081, 'v', markersize=10, color='black');

    cbarax = fig.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(filename + ".png");

    return;


def drive_ou_cornell_comparison(myLev, leveling_slice, data_file, los_file, s1_slice, output_dir):
    InSAR_Data = InSAR_Object.inputs.inputs_cornell_ou_velocities_hdf5(data_file, los_file, s1_slice);
    InSAR_Data = InSAR_Object.utilities.remove_nans(InSAR_Data);
    vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, InSAR_Data);
    one_to_one_comparison(myLev, InSAR_Data, vector_index, leveling_slice[0], leveling_slice[1], "S1", output_dir);
    return;


def drive_single_uavsar_intf_comparison(myLev, leveling_slice, bounds, uavsar_filename, los_filename, output_dir):
    # Read the UAVSAR Data
    InSAR_Data = InSAR_Object.inputs.inputs_isce_unw_geo_losrdr(uavsar_filename, los_filename, bounds);
    InSAR_Data = InSAR_Object.utilities.remove_nans(InSAR_Data);
    vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, InSAR_Data);
    one_to_one_comparison(myLev, InSAR_Data, vector_index, leveling_slice[0], leveling_slice[1], "UAV", output_dir, flip_insar=True);
    return;


if __name__ == "__main__":
    # Opening stuff
    config_filename = "config_file.txt"
    file_dict = multiSAR_input_functions.get_file_dictionary(config_filename);
    myLev = multiSAR_input_functions.inputs_leveling(file_dict["leveling"].split()[0], file_dict["leveling"].split()[1]);
    myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);

    # # # S1_Cornell experiment (2015 data)
    # output_dir = "S1_OU/T4D/";
    # leveling_slice = [5, 6];
    # s1_slice = 0;
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice, output_dir + 'ascending_');
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice, output_dir + 'descending_');
    #
    # # # S1_Cornell experiment (2016 data)
    # output_dir = "S1_OU/T4E/";
    # leveling_slice = [6, 7];
    # s1_slice = 1;
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice, output_dir + 'ascending_');
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice, output_dir + 'descending_');
    #
    # # # S1_Cornell experiment (2017 data)
    # output_dir = "S1_OU/T4F/";
    # leveling_slice = [7, 8];
    # s1_slice = 2;
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice, output_dir + 'ascending_');
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice, output_dir + 'descending_');
    #
    # # # S1_Cornell experiment (2018 data)
    # output_dir = "S1_OU/T5/";
    # leveling_slice = [8, 9];
    # s1_slice = 3;
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_ascending"], file_dict["s1_ou_ascending_los"],
    #                             s1_slice, output_dir + 'ascending_');
    # drive_ou_cornell_comparison(myLev, leveling_slice, file_dict["s1_ou_descending"], file_dict["s1_ou_descending_los"],
    #                             s1_slice, output_dir + 'descending_');
    #

    # Individual UAVSAR experiments, starting with 2011-2012
    output_dir = "UAVSAR_intfs/"
    leveling_slice = [2, 3];
    bounds = (dt.datetime.strptime("20111110", "%Y%m%d"), dt.datetime.strptime("20120926", "%Y%m%d"));
    drive_single_uavsar_intf_comparison(myLev, leveling_slice, bounds,
                                        file_dict["uavsar_08508_2011_2012_unw"],
                                        file_dict["uavsar_08508_2011_2012_los"], output_dir);

    # #TSX experiment
    # output_dir = "TSX/"
    # VertTSXData, EastTSXData = InSAR_Object.inputs.inputs_TRE(file_dict["tsx"]);
    # vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, VertTSXData);
    # one_to_one_comparison(myLev, VertTSXData, vector_index, 3, 4, "TSX", output_dir);  # the 2012-2013 interval
    # plot_InSAR(VertTSXData,output_dir);

    # # S1 experiment
    # output_dir = "SNT1/"
    # VertS1Data, EastS1Data = InSAR_Object.inputs.inputs_TRE(file_dict["snt1"]);
    # vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, VertS1Data);
    # one_to_one_comparison(myLev, VertS1Data, vector_index, 5, 8, "S1", output_dir); # the 2014-2018 interval
    # plot_InSAR(VertS1Data,output_dir);

    # # S1 experiment
    # output_dir = "SNT2/"
    # VertS1Data, EastS1Data = InSAR_Object.inputs.inputs_TRE(file_dict["snt2"]);
    # vector_index = multiSAR_utilities.find_leveling_in_vector(myLev, VertS1Data);
    # one_to_one_comparison(myLev, VertS1Data, vector_index, 8, 9, "S1", output_dir); # the 2014-2018 interval
    # plot_InSAR(VertS1Data,output_dir);
