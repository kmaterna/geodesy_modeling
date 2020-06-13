# May 2020
# Here we are getting ready for a multi-datasource inversion
# STEPS: 
# 1. Get GPS; make relative to benchmark
# 2. Put InSAR relative to chosen GPS benchmark; solve for ramps etc; downsample
# 3. Get Leveling

import numpy as np 
import subprocess
import json
import sys
import multiSAR_input_functions
import quadtree_downsample_kite
import downsample_gps_ts

def welcome_and_parse(argv):
	print("Welcome to the multiSAR data input pipeline.")
	if len(argv)<2:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		configname = argv[1];
	config_file = open(configname,'r');
	config = json.load(config_file);
	return config;


def downsample_and_cut_uavsar(config):
	"""
	For UAVSAR, we quadtree downsample, multiply by -1, and chop. 
	"""
	# Prepping UAVSAR: Input and Compute
	print("Starting to prepare UAVSAR data");
	scene0 = multiSAR_input_functions.inputs_uavsar_unw_geo(config["uav_datafile0"]);
	scene1 = multiSAR_input_functions.inputs_uavsar_unw_geo(config["uav_datafile1"]);
	data = np.float32(np.subtract(scene1, scene0));  # must be 4-byte floats for quadtree
	if config["flip_uavsar_los"]:
		data = -1*data;  # away from satellite = negative motion

	# Setup Downsampling: rdrfile, xmlfile, datafile
	uavsar_unw_file = config["prep_inputs_dir"]+config["uavsar_unw_file"];
	geojson_file = config["prep_inputs_dir"]+config["geojson_file"];
	uav_plotfile = config["prep_inputs_dir"]+config["uav_plotfile"];
	uav_textfile = config["prep_inputs_dir"]+config["uav_textfile"];
	subprocess.call(['cp',config["rdrfile"],config["prep_inputs_dir"]],shell=False);
	subprocess.call(['cp',config["uav_datafile0"]+'.xml',uavsar_unw_file+".xml"],shell=False);
	data.tofile(uavsar_unw_file) # write data bytes out

	# # Downsample, bbox, and Print
	quadtree_downsample_kite.kite_downsample_isce_unw(uavsar_unw_file, geojson_file, 
		config["epsilon"], config["nan_allowed"], config["tile_size_min"], config["tile_size_max"]);
	quadtree_downsample_kite.geojson_to_outputs(geojson_file, uav_plotfile, uav_textfile, bbox=config["uavsar_bbox"], std_min=config["uavsar_std_min"]);
	return;

def get_leveling_displacements(config):
	"""
	For leveling, we only have to write the proper format text file. 
	"""
	myLev = multiSAR_input_functions.inputs_leveling(config["leveling_filename"], config["leveling_errors_filename"]);
	myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);  # Computing relative displacements from 2009 onward
	multiSAR_input_functions.write_leveling_invertible_format(myLev, config["leveling_start"], config["leveling_end"], config["leveling_unc"], config["prep_inputs_dir"]+config["lev_outfile"]);
	multiSAR_input_functions.plot_leveling(config["prep_inputs_dir"]+config["lev_outfile"], config["prep_inputs_dir"]+config["lev_plot"]);
	return;

def get_gps_displacements(config):
	"""
	For GPS, we only have to write the proper format text file, and sometimes we make other corrections. 
	"""
	print("Starting to extract GPS from %s to %s " % (config["starttime"], config["endtime"]) );
	stations = downsample_gps_ts.read_station_ts(config);
	displacement_objects = downsample_gps_ts.get_displacements_show_ts(config, stations);
	multiSAR_input_functions.write_gps_invertible_format(displacement_objects, config["prep_inputs_dir"]+config["gps_textfile"]);
	return;



if __name__=="__main__":
	config=welcome_and_parse(sys.argv);
	downsample_and_cut_uavsar(config);
	get_leveling_displacements(config);
	get_gps_displacements(config);


