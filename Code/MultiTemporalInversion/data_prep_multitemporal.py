# June 2020
# A similar pipeline to the single-interval case, but with more config

import numpy as np 
import matplotlib.pyplot as plt 
import sys, json, os
import subprocess
import datetime as dt 
import quadtree_downsample_kite
import downsample_gps_ts
import multiSAR_input_functions
import remove_insar_ramp
import insar_LOS_tools

def welcome_and_parse(argv):
	print("Welcome to the multiSAR data input pipeline.")
	if len(argv)<2:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		configname = argv[1];
	config_file = open(configname,'r');
	config = json.load(config_file);
	return config;

def get_starttime_endtime(epochs_dict, select_interval_dict):
	# Given a dictionary of epochs and a list called span, extract the starttime and endtime for that interval
	start_time_candidates = [];
	end_time_candidates = [];
	span = select_interval_dict["span"];  # a list of intervals, such as ["A","B"];
	for item in span:
		for interval_key in epochs_dict:
			if epochs_dict[interval_key]["name"] in item:
				start_time_candidates.append(dt.datetime.strptime(epochs_dict[interval_key]["starttime"],"%Y-%m-%d"));
				end_time_candidates.append(dt.datetime.strptime(epochs_dict[interval_key]["endtime"],"%Y-%m-%d"));
	starttime = min(start_time_candidates);
	endtime = max(end_time_candidates);
	return starttime, endtime

def write_gps_displacements_multiple(config):
	"""
	For GPS, we only have to write the proper format text file, and sometimes we make other corrections. 
	"""
	# For each interval in gps: 
	prep_dir = config["prep_inputs_dir"];
	for interval_dict_key in config["gps_data"]:
		new_interval_dict = config["gps_data"][interval_dict_key];  # for each interval in GPS
		gps_sigma = new_interval_dict["gps_sigma"];
		starttime, endtime = get_starttime_endtime(config["epochs"],new_interval_dict);
		print("For GPS %s, starting to extract GPS from %s to %s " % (interval_dict_key, starttime, endtime) );
		stations = downsample_gps_ts.read_station_ts(new_interval_dict["gps_bbox"],new_interval_dict["gps_reference"]);
		displacement_objects = downsample_gps_ts.get_displacements_show_ts(stations, starttime, endtime, gps_sigma, prep_dir); 
		multiSAR_input_functions.write_gps_invertible_format(displacement_objects, config["prep_inputs_dir"]+new_interval_dict["gps_textfile"]);
	return;

def write_leveling_displacements_multiple(config):
	"""
	For leveling, we only have to write the proper format text file. 
	"""
	for interval_dict_key in config["leveling_data"]:
		new_interval_dict = config["leveling_data"][interval_dict_key];  # for each interval in Leveling
		myLev = multiSAR_input_functions.inputs_leveling(new_interval_dict["leveling_filename"], new_interval_dict["leveling_errors_filename"]);
		myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);  # Computing relative displacements from 2009 onward
		multiSAR_input_functions.write_leveling_invertible_format(myLev, new_interval_dict["leveling_start"], new_interval_dict["leveling_end"], new_interval_dict["leveling_unc"], config["prep_inputs_dir"]+new_interval_dict["lev_outfile"]);
		multiSAR_input_functions.plot_leveling(config["prep_inputs_dir"]+new_interval_dict["lev_outfile"], config["prep_inputs_dir"]+new_interval_dict["lev_plot"]);
	return;

def downsample_cut_write_uavsar_multiple(config):
	"""
	For UAVSAR, we quadtree downsample, multiply by -1, and chop. 
	"""
	for interval_dict_key in config["uavsar_data"]:
		new_interval_dict = config["uavsar_data"][interval_dict_key];  # for each interval in UAVSAR

		# Prepping UAVSAR: Input and Compute
		print("Starting to prepare UAVSAR data");
		scene0 = multiSAR_input_functions.inputs_uavsar_unw_geo(new_interval_dict["uav_datafile0"]);
		scene1 = multiSAR_input_functions.inputs_uavsar_unw_geo(new_interval_dict["uav_datafile1"]);
		data = np.float32(np.subtract(scene1, scene0));  # must be 4-byte floats for quadtree
		if new_interval_dict["flip_uavsar_los"]:
			data = -1*data;  # away from satellite = negative motion

		# Setup Downsampling: rdrfile, xmlfile, datafile
		uavsar_unw_file = config["prep_inputs_dir"]+new_interval_dict["uavsar_unw_file"];
		geojson_file = config["prep_inputs_dir"]+new_interval_dict["geojson_file"];
		uav_plotfile = config["prep_inputs_dir"]+new_interval_dict["uav_plotfile"];
		uav_textfile = config["prep_inputs_dir"]+new_interval_dict["uav_textfile"];
		subprocess.call(['cp',new_interval_dict["rdrfile"],config["prep_inputs_dir"]],shell=False);
		subprocess.call(['cp',new_interval_dict["uav_datafile0"]+'.xml',uavsar_unw_file+".xml"],shell=False);
		data.tofile(uavsar_unw_file) # write data bytes out

		# Downsample, bbox, and Print
		quadtree_downsample_kite.kite_downsample_isce_unw(uavsar_unw_file, geojson_file, 
			new_interval_dict["epsilon"], new_interval_dict["nan_allowed"], new_interval_dict["tile_size_min"], new_interval_dict["tile_size_max"]);
		quadtree_downsample_kite.geojson_to_outputs(geojson_file, uav_plotfile, uav_textfile, bbox=new_interval_dict["uavsar_bbox"], std_min=new_interval_dict["uavsar_std_min"]);

		# Now we remove a ramp. 
		if "remove_uavsar_ramp" in new_interval_dict.keys():
			ramp_adjusted_file = uav_textfile.split(".txt")[0]+"_ramp_removed.txt";  # no-ramps file
			remove_insar_ramp.remove_ramp(uav_textfile, ramp_adjusted_file);

	return;

def write_tsx_tre_displacements_multiple(config):
	for interval_dict_key in config["tsx_data"]:
		new_interval_dict = config["tsx_data"][interval_dict_key];  # for each interval in TSX
		print("Starting to extract TSX TRE-format from %s " % (new_interval_dict["tsx_filename"]) );
	
		TRE_TSX = multiSAR_input_functions.inputs_TRE(new_interval_dict["tsx_filename"]);
		Vert_InSAR, East_InSAR = insar_LOS_tools.TRE_to_InSAR_Obj(TRE_TSX);  # convert to displacements
		Vert_InSAR = insar_LOS_tools.impose_InSAR_bounding_box(Vert_InSAR, new_interval_dict["tsx_bbox"]);
		Vert_InSAR = insar_LOS_tools.uniform_downsampling(Vert_InSAR,new_interval_dict["tsx_downsample_interval"], new_interval_dict["tsx_averaging_window"]);  # uniform downsampling

		multiSAR_input_functions.write_insar_invertible_format(Vert_InSAR, new_interval_dict["tsx_unc"], config["prep_inputs_dir"]+new_interval_dict["tsx_datafile"]);
		multiSAR_input_functions.plot_insar(config["prep_inputs_dir"]+new_interval_dict["tsx_datafile"], config["prep_inputs_dir"]+new_interval_dict["tsx_plot"]);
	return;

if __name__=="__main__":
	config=welcome_and_parse(sys.argv);
	# downsample_cut_write_uavsar_multiple(config);
	# write_leveling_displacements_multiple(config);
	# write_gps_displacements_multiple(config);
	write_tsx_tre_displacements_multiple(config);
	# write_s1_tre_displacements(config);
