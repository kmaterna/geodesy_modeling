# May 2020
# Here we are getting ready for a multi-datasource inversion
# STEPS: 
# 1. Get GPS; make relative to benchmark
# 2. Put InSAR relative to chosen GPS benchmark; solve for ramps etc; downsample
# 3. Get Leveling

import numpy as np 
import subprocess
import json
import sys, os
import datetime as dt
import multiSAR_input_functions
import quadtree_downsample_kite
import downsample_gps_ts
import remove_coseismic_model
import remove_insar_ramp

def welcome_and_parse(argv):
	print("Welcome to the multiSAR data input pipeline.")
	if len(argv)<2:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		configname = argv[1];
	config_file = open(configname,'r');
	config = json.load(config_file);
	return config;


def downsample_cut_write_uavsar(config):
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

	# Downsample, bbox, and Print
	quadtree_downsample_kite.kite_downsample_isce_unw(uavsar_unw_file, geojson_file, 
		config["epsilon"], config["nan_allowed"], config["tile_size_min"], config["tile_size_max"]);
	quadtree_downsample_kite.geojson_to_outputs(geojson_file, uav_plotfile, uav_textfile, bbox=config["uavsar_bbox"], std_min=config["uavsar_std_min"]);

	# Run the EMC model on the uavsar_los.txt and get the results back into this directory.
	# Subtract the model LOS and write it out. 
	if "adjust_EMC_uavsar" in config.keys():
		points_for_modeling_file = "Edit_for_EMC/EMC_coseismic/EMC_coseismic_model/Inputs/uavsar_los.txt";
		model_output_points = "Edit_for_EMC/EMC_coseismic/EMC_coseismic_model/Outputs/EMC/ll_disps.txt";
		subprocess.call(["cp",uav_textfile,points_for_modeling_file],shell=False);
		add_reference_pixel_los_to_end(points_for_modeling_file, config["reference_ll"]); # Adding the reference pixel for easy subtracting
		os.chdir("Edit_for_EMC/EMC_coseismic/EMC_coseismic_model");
		subprocess.call(["python","/Users/kmaterna/Documents/B_Research/Github_Tools/Elastic_stresses_py/Code/elastic_stresses_driver.py","config.txt"],shell=False);
		os.chdir("../../../");	
		subprocess.call(["cp",model_output_points,config["adjust_EMC_uavsar"]],shell=False);
		adjusted_file = uav_textfile.split(".txt")[0]+"_cos_corrected.txt";  # new file
		remove_coseismic_model.remove_model_los(uav_textfile, config["adjust_EMC_uavsar"], adjusted_file);
		uav_textfile = adjusted_file;  # get ready for the next step. 
	
	# A quick fix if you're skipping the EMC correction for speed. 
	# uav_textfile = uav_textfile.split(".txt")[0]+"_cos_corrected.txt";  # new file
	# ramp_adjusted_file = uav_textfile.split(".txt")[0]+"_cos_corrected_ramp_removed.txt";  # no-ramps file

	# Now we remove a ramp. 
	if "remove_uavsar_ramp" in config.keys():
		ramp_adjusted_file = uav_textfile.split(".txt")[0]+"_ramp_removed.txt";  # no-ramps file
		remove_insar_ramp.remove_ramp(uav_textfile, ramp_adjusted_file);

	return;


def add_reference_pixel_los_to_end(uav_textfile,reference_ll):
	# This is called before EMC correction, as we need the reference for applying the correction. 
	# We usually remove it later. 
	ofile=open(uav_textfile,'a');
	ofile.write("%f %f 0.0000 0.01 -1 -1 -1" % (reference_ll[0], reference_ll[1]) );
	ofile.close();
	return;


def write_leveling_displacements(config):
	"""
	For leveling, we only have to write the proper format text file. 
	"""
	myLev = multiSAR_input_functions.inputs_leveling(config["leveling_filename"], config["leveling_errors_filename"]);
	myLev = multiSAR_input_functions.compute_rel_to_datum_nov_2009(myLev);  # Computing relative displacements from 2009 onward
	multiSAR_input_functions.write_leveling_invertible_format(myLev, config["leveling_start"], config["leveling_end"], config["leveling_unc"], config["prep_inputs_dir"]+config["lev_outfile"]);
	multiSAR_input_functions.plot_leveling(config["prep_inputs_dir"]+config["lev_outfile"], config["prep_inputs_dir"]+config["lev_plot"]);
	return;

def write_gps_displacements(config):
	"""
	For GPS, we only have to write the proper format text file, and sometimes we make other corrections. 
	"""
	print("Starting to extract GPS from %s to %s " % (config["starttime"], config["endtime"]) );
	stations = downsample_gps_ts.read_station_ts(config);
	displacement_objects = downsample_gps_ts.get_displacements_show_ts(config, stations);
	multiSAR_input_functions.write_gps_invertible_format(displacement_objects, config["prep_inputs_dir"]+config["gps_textfile"]);
	
	if "adjust_EMC_gps" in config.keys():
		print("Adjusting GPS for EMC gradients");
		# Manually collect the modeled gps points here. 
		remove_coseismic_model.remove_model_gps(config["prep_inputs_dir"]+config["gps_textfile"], config["adjust_EMC_gps"]);	
	return;

def write_s1_tre_displacements(config):
	"""
	For S1, we read the format, and write the insar file 
	Might want to downsample in the future?  That doesn't have to be so fancy as quadtree. 
	"""
	if "s1_filename" not in config.keys():
		print("No S1 in this inversion");
		return;
	print("Starting to extract S1 TRE-format from %s " % (config["s1_filename"]) );
	TRE_S1 = multiSAR_input_functions.inputs_S1(config["s1_filename"]);
	multiSAR_input_functions.write_tre_invertible_format(TRE_S1, config["s1_bbox"], config["s1_unc"], config["prep_inputs_dir"]+config["s1_datafile"]);
	return;

def write_tsx_tre_displacements(config):
	"""
	For TSX, we read the format, and write the insar file 
	Might want to downsample in the future?  That doesn't have to be so fancy as quadtree. 
	"""	
	if "tsx_filename" not in config.keys():
		print("No TSX in this inversion");
		return;	
	print("Starting to extract TSX TRE-format from %s " % (config["tsx_filename"]) );
	TRE_TSX = multiSAR_input_functions.inputs_tsx(config["tsx_filename"]);
	multiSAR_input_functions.write_tre_invertible_format(TRE_TSX, config["tsx_bbox"], config["tsx_unc"], config["prep_inputs_dir"]+config["tsx_datafile"]);	
	return; 


if __name__=="__main__":
	config=welcome_and_parse(sys.argv);
	downsample_cut_write_uavsar(config);
	write_leveling_displacements(config);
	write_gps_displacements(config);
	write_s1_tre_displacements(config);
	write_tsx_tre_displacements(config);


