#!/usr/bin/env python

# Prepare data for a multi-datasource, multi-time-step inversion
# GNSS: make relative to benchmark
# INSAR: make relative to chosen benchmark; solve for ramps; downsample
# Leveling: write out

from Geodesy_Modeling import Leveling_Object
import numpy as np
import sys, json
import subprocess
import datetime as dt
from Geodesy_Modeling import Downsample
from Geodesy_Modeling import InSAR_1D_Object
from Geodesy_Modeling import GNSS_Object
from read_write_insar_utilities import isce_read_write


def welcome_and_parse(argv):
    print("Welcome to the multiSAR data input pipeline.")
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ");
        sys.exit(0);
    else:
        configname = argv[1];
    config_file = open(configname, 'r');
    config = json.load(config_file);
    return config;


def get_starttime_endtime(epochs_dict, select_interval_dict):
    """Given a dictionary of epochs and a list called span, extract the starttime and endtime for that interval"""
    start_time_candidates = [];
    end_time_candidates = [];
    span = select_interval_dict["span"];  # a list of intervals, such as ["A","B"];
    for item in span:
        for interval_key in epochs_dict:
            if epochs_dict[interval_key]["name"] in item:
                start_time_candidates.append(dt.datetime.strptime(epochs_dict[interval_key]["starttime"], "%Y-%m-%d"));
                end_time_candidates.append(dt.datetime.strptime(epochs_dict[interval_key]["endtime"], "%Y-%m-%d"));
    starttime = min(start_time_candidates);
    endtime = max(end_time_candidates);
    return starttime, endtime


def write_gps_displacements(config):
    """
    For GPS, we have to write the proper format text file, and sometimes we make other corrections.
    """
    # For each interval in gps:
    prep_dir = config["prep_inputs_dir"];
    if "gps_data" not in config.keys():
        print("\nNo GPS in this inversion");
        return;
    for interval_dict_key in config["gps_data"]:
        new_interval_dict = config["gps_data"][interval_dict_key];  # for each interval in GPS
        gps_sigma = new_interval_dict["gps_sigma"];
        starttime, endtime = get_starttime_endtime(config["epochs"], new_interval_dict);

        print("\nFor GPS %s, starting to extract GPS from %s to %s " % (interval_dict_key, starttime, endtime));
        stations = GNSS_Object.read_gnss.read_station_ts(new_interval_dict["gps_bbox"],
                                                         new_interval_dict["gps_reference"],
                                                         remove_coseismic=new_interval_dict["remove_coseismic"]);
        displacement_objects = Downsample.downsample_gps_ts.get_displacements_show_ts(stations, starttime, endtime,
                                                                                      gps_sigma,
                                                                                      prep_dir);
        if "gps_add_offset_mm" in new_interval_dict.keys():  # an option to add a constant (in enu) to the GNSS offsets
            displacement_objects = GNSS_Object.utilities.add_gps_constant_offset(displacement_objects,
                                                                                 new_interval_dict[
                                                                                     "gps_add_offset_mm"]);
        GNSS_Object.utilities.write_gps_invertible_format(displacement_objects, config["prep_inputs_dir"]
                                                          + new_interval_dict["gps_textfile"]);
    return;


def write_leveling_displacements(config):
    """
    For leveling, we only have to write the proper format text file.
    """
    if "leveling_data" not in config.keys():
        print("\nNo Leveling in this inversion");
        return;
    for interval_dict_key in config["leveling_data"]:
        new_interval_dict = config["leveling_data"][interval_dict_key];  # for each interval in Leveling
        print("\nPreparing leveling for file %s" % new_interval_dict["lev_outfile"])
        myLev = Leveling_Object.leveling_inputs.inputs_brawley_leveling(new_interval_dict["leveling_filename"],
                                                                        new_interval_dict["leveling_errors_filename"]);
        myLev = Leveling_Object.leveling_inputs.compute_rel_to_datum_nov_2009(myLev);  # Relative disp after 2009
        Leveling_Object.leveling_outputs.write_leveling_invertible_format(myLev, new_interval_dict["leveling_start"],
                                                                          new_interval_dict["leveling_end"],
                                                                          new_interval_dict["leveling_unc"],
                                                                          config["prep_inputs_dir"] +
                                                                          new_interval_dict["lev_outfile"]);
        Leveling_Object.leveling_outputs.plot_simple_leveling(config["prep_inputs_dir"] +
                                                              new_interval_dict["lev_outfile"],
                                                              config["prep_inputs_dir"] +
                                                              new_interval_dict["lev_plot"]);
    return;


def write_uavsar_displacements(config):
    """
    For UAVSAR, we quadtree downsample, multiply by -1, and chop.
    """
    if "uavsar_data" not in config.keys():
        print("\nNo UAVSAR in this inversion");
        return;
    for interval_dict_key in config["uavsar_data"]:
        print("\nStarting to prepare UAVSAR data for %s" % interval_dict_key);
        new_interval_dict = config["uavsar_data"][interval_dict_key];  # for each interval in UAVSAR

        # Get uavsar data
        if 'uav_sourcefile_begin' in new_interval_dict.keys():  # Using time series format
            source_xml_name = new_interval_dict["uav_sourcefile_begin"] + ".xml";
            scene0 = isce_read_write.read_scalar_data(new_interval_dict["uav_sourcefile_begin"], band=2);
            scene1 = isce_read_write.read_scalar_data(new_interval_dict["uav_sourcefile_end"], band=2);
            data = np.float32(np.subtract(scene1, scene0));  # must be 4-byte floats for quadtree
        else:  # using an individual interferogram format
            source_xml_name = new_interval_dict["uav_sourcefile_unw_geo"] + ".xml";
            data = isce_read_write.read_scalar_data(new_interval_dict["uav_sourcefile_unw_geo"], band=2);

        if new_interval_dict["flip_uavsar_los"] == 1:
            data = -1 * data;  # away from satellite = negative motion
            print("Multiplying UAVSAR by -1 for LOS motion sign convention. ");

        # Write the .unw.geo file and metadata
        uavsar_unw_file = config["prep_inputs_dir"] + new_interval_dict["uavsar_unw_file"];  # output file
        subprocess.call(['cp', source_xml_name, uavsar_unw_file + ".xml"]);  # write the xml out
        isce_read_write.data_to_file_2_bands(data, data, filename=uavsar_unw_file);  # write data bytes out.

        # Quadtree downsampling by Kite
        uav_textfile = config["prep_inputs_dir"] + new_interval_dict["uav_textfile"];  # .txt, invertible format
        drive_uavsar_kite_downsampling(new_interval_dict, config["prep_inputs_dir"]);

        # Now we optionally remove a ramp.
        if new_interval_dict["remove_ramp"] == 1:
            InSAR_1D_Object.remove_ramp.remove_ramp_filewise(uav_textfile, uav_textfile,
                                                             ref_coord=config['reference_ll']);

        # Now we optionally remove a constant
        if new_interval_dict["remove_constant"] == 1:
            InSAR_1D_Object.remove_ramp.remove_constant_filewise(uav_textfile, uav_textfile);

        # Now we make a plot
        InSAR_Obj = InSAR_1D_Object.inputs.inputs_txt(uav_textfile);
        InSAR_1D_Object.outputs.plot_insar(InSAR_Obj, config["prep_inputs_dir"] + new_interval_dict["uav_ending_plot"]);
    return;


def drive_uavsar_kite_downsampling(interval_dictionary, inputs_dir):
    """Setup Downsampling: rdrfile, xmlfile, datafile"""
    uavsar_unw_file = inputs_dir + interval_dictionary["uavsar_unw_file"];
    geojson_file = inputs_dir + interval_dictionary["geojson_file"];
    uav_plotfile = inputs_dir + interval_dictionary["uav_plotfile"];
    uav_textfile = inputs_dir + interval_dictionary["uav_textfile"];
    print("Copying %s into directory %s" % (interval_dictionary['rdrfile'], inputs_dir));
    subprocess.call(['cp', interval_dictionary["rdrfile"], inputs_dir], shell=False);

    # Downsample, bbox, and Print
    Downsample.quadtree_downsample_kite.kite_downsample_isce_unw(uavsar_unw_file, geojson_file,
                                                                 interval_dictionary["epsilon"],
                                                                 interval_dictionary["nan_allowed"],
                                                                 interval_dictionary["tile_size_min"],
                                                                 interval_dictionary["tile_size_max"]);
    Downsample.quadtree_downsample_kite.geojson_to_outputs(geojson_file, uav_plotfile, uav_textfile,
                                                           bbox=interval_dictionary["uavsar_bbox"],
                                                           std_min=interval_dictionary["uavsar_std_min"]);
    return;


def write_tsx_tre_displacements(config):
    """
    For TSX, we read the format, and write the insar file
    We also downsample and impose bounding box
    In the case of TSX, now trying both east and vertical in the inversion.
    """
    if "tsx_data" not in config.keys():
        print("\nNo TSX in this inversion");
        return;
    else:
        for interval_dict_key in config["tsx_data"]:
            new_interval_dict = config["tsx_data"][interval_dict_key];  # for each interval in TSX
            print("\nStarting to extract TSX TRE-format from %s " % (new_interval_dict["tsx_filename"]));

            # We can get both vertical and east from the TRE data.
            Vert_InSAR, East_InSAR = InSAR_1D_Object.inputs.inputs_TRE_vert_east(new_interval_dict["tsx_filename"]);
            Vert_InSAR = InSAR_1D_Object.utilities.impose_InSAR_bounding_box(Vert_InSAR, new_interval_dict[
                "tsx_bbox"]);  # bounding box vertical
            East_InSAR = InSAR_1D_Object.utilities.impose_InSAR_bounding_box(East_InSAR, new_interval_dict[
                "tsx_bbox"]);  # bounding box vertical
            Vert_InSAR = Downsample.uniform_downsample.uniform_downsampling(Vert_InSAR,
                                                                            new_interval_dict[
                                                                                "tsx_downsample_interval"],
                                                                            new_interval_dict["tsx_averaging_window"]);
            East_InSAR = Downsample.uniform_downsample.uniform_downsampling(East_InSAR,
                                                                            new_interval_dict[
                                                                                "tsx_downsample_interval"],
                                                                            new_interval_dict["tsx_averaging_window"]);

            Total_InSAR = InSAR_1D_Object.utilities.combine_objects(Vert_InSAR, East_InSAR);

            InSAR_1D_Object.outputs.write_insar_invertible_format(Total_InSAR, new_interval_dict["tsx_unc"],
                                                                  config["prep_inputs_dir"] +
                                                                  new_interval_dict["tsx_datafile"]);
            InSAR_obj = InSAR_1D_Object.inputs.inputs_simplest_txt(config["prep_inputs_dir"] +
                                                                   new_interval_dict["tsx_datafile"]);
            InSAR_1D_Object.outputs.plot_insar(InSAR_obj, config["prep_inputs_dir"] + new_interval_dict["tsx_plot"]);

    return;


def write_s1_displacements(config):
    """
    For S1, we read the format, and write the insar file
    We also downsample and impose bounding box
    """
    if "s1_data" not in config.keys():
        print("\nNo S1 in this inversion");
        return;
    else:
        for interval_dict_key in config["s1_data"]:
            new_interval_dict = config["s1_data"][interval_dict_key];  # for each interval in TSX
            print("\nStarting to extract S1 Cornell/OU-format from %s " % (new_interval_dict["s1_filename"]));
            InSAR_Data = InSAR_1D_Object.inputs.inputs_cornell_ou_velocities_hdf5(new_interval_dict["s1_filename"],
                                                                                  new_interval_dict["s1_lkv_filename"],
                                                                                  new_interval_dict["s1_slicenum"]);
            InSAR_Data = InSAR_1D_Object.utilities.impose_InSAR_bounding_box(InSAR_Data, new_interval_dict["s1_bbox"]);
            InSAR_Data = Downsample.uniform_downsample.uniform_downsampling(InSAR_Data,
                                                                            new_interval_dict["s1_downsample_interval"],
                                                                            new_interval_dict["s1_averaging_window"]);

            InSAR_1D_Object.outputs.write_insar_invertible_format(InSAR_Data, new_interval_dict["s1_unc"],
                                                                  config["prep_inputs_dir"] + new_interval_dict[
                                                                   "s1_datafile"]);
            InSAR_obj = InSAR_1D_Object.inputs.inputs_simplest_txt(config["prep_inputs_dir"] +
                                                                   new_interval_dict["s1_datafile"]);
            InSAR_1D_Object.outputs.plot_insar(InSAR_obj, config["prep_inputs_dir"] + new_interval_dict["s1_plot"]);

    return;


if __name__ == "__main__":
    config = welcome_and_parse(sys.argv);
    write_uavsar_displacements(config);
    write_leveling_displacements(config);
    write_gps_displacements(config);
    write_tsx_tre_displacements(config);
    write_s1_displacements(config);
