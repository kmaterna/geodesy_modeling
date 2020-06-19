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
import datetime as dt
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

	# For this to work, I've manually taken "uavsar_los.txt" and produced a E/N/U predictions at each pixel. 
	# It lives in the same directory as uav_textfile. 	
	if "adjust_EMC_uavsar" in config.keys():
		# add_reference_pixel_to_end(uav_textfile, config["reference_ll"]); # Adding the reference pixel for easy subtracting
		# HERE We need to run the EMC model on the uavsar_los.txt and get the results back into this directory. 
		# Manual right now, but could be automatic later. 
		# We need to subtract the model LOS and write it out. 
		remove_emc_model_los(uav_textfile, config["adjust_EMC_uavsar"]);
	return;


def remove_emc_model_los(los_file, model_disps_file):
	# Pseudocode: 
	# What is the modeled delta-e, delta-n, and delta-u of each point relative to the reference point? 
	# Then, how does that project into the local LOS? 
	# Finally, subtract that offset from each point
	# Then we write everything except the reference line
	# Into a file with "_updated" in its name 	
	[lon_pred, lat_pred, u_pred, v_pred, w_pred] = np.loadtxt(model_disps_file,unpack=True, skiprows=1);
	[lon_meas, lat_meas, disp, sig, unit_e, unit_n, unit_u] = np.loadtxt(los_file, unpack=True, skiprows=1);
	corrected_los = [];
	los_reference_pixel = [u_pred[-1], v_pred[-1], w_pred[-1]];
	for i in range(len(lon_meas)-1):
		model_prediction = [u_pred[i], v_pred[i], w_pred[i]];
		model_deltas = np.subtract(model_prediction, los_reference_pixel);
		
		los_unitv = [unit_e[i], unit_n[i], unit_u[i]];
		los_view = np.dot(model_deltas, los_unitv);

		corrected_los.append ( disp[i] - los_view );

	namestem = los_file.split(".txt")[0];
	ofile=open(namestem+"_emccorrected.txt",'w');
	ofile.write("# Header: lon, lat, disp(m), sig(m), unitE, unitN, unitU from ground to satellite\n");
	for i in range(len(lon_meas)-1):
		ofile.write("%f %f %f %f %f %f %f \n" % (lon_meas[i], lat_meas[i], corrected_los[i], sig[i], unit_e[i], unit_n[i], unit_u[i]) );
	ofile.close();
	return;

def remove_emc_model_gps(gps_file, model_disps_file):
	# Take the prediction at each station, 
	# subtract the reference prediction
	# and then adjust by that amount. 
	[lon_pred, lat_pred, u_pred, v_pred, w_pred] = np.loadtxt(model_disps_file,unpack=True, skiprows=1);
	[lon_meas, lat_meas, u_meas, v_meas, w_meas, sige_meas, sign_meas, sigu_meas] = np.loadtxt(gps_file, unpack=True, skiprows=1);
	disp_reference_pixel = [u_pred[0], v_pred[0], w_pred[0]];
	namestem = gps_file.split(".txt")[0];
	ofile=open(namestem+"_emccorrected.txt",'w');
	ofile.write("# Header: lon, lat, dE, dN, dU, Se, Sn, Su (m)\n");
	for i in range(len(lon_meas)):
		disp_at_station = [u_meas[i], v_meas[i], w_meas[i]];
		model_for_station = [u_pred[i], v_pred[i], w_pred[i]];
		model_delta = np.subtract(model_for_station, disp_reference_pixel);
		new_data = np.subtract(disp_at_station, model_delta);
		ofile.write("%f %f %f %f %f %f %f %f\n" % (lon_meas[i], lat_meas[i], new_data[0], new_data[1], new_data[2], sige_meas[i], sign_meas[i], sigu_meas[i]) );
	ofile.close();
	return;

def add_reference_pixel_to_end(uav_textfile,reference_ll):
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
		remove_emc_model_gps(config["prep_inputs_dir"]+config["gps_textfile"], config["adjust_EMC_gps"]);
	
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


