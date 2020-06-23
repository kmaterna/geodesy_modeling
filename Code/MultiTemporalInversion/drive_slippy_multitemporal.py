# June 2020
# Run Slippy with multiple time intervals in the input 
# (A really big G matrix)


import numpy as np 
import sys
import json
import buildG
import slippy.io



def welcome_and_parse(argv):
	print("Welcome to the MULTITEMPORAL INVERSION.");
	if len(argv)<3:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		fault_config = argv[1];
		data_config = argv[2];
	fault_config_file = open(fault_config,'r');
	config1 = json.load(fault_config_file);
	data_config_file = open(data_config,'r');
	config2 = json.load(data_config_file);
	returnval = {**config1, **config2};
	return returnval;

def read_in_data(config):
	if "gps_data" in config.keys():
		for interval in config["gps_data"].keys():
			filename = config["prep_inputs_dir"]+config["gps_data"][interval]["gps_textfile"];
			print("Reading %s " % filename);
			gps_input = slippy.io.read_gps_data(filename);
			config["gps_data"][interval]["data"] = gps_input;

	if "uavsar_data" in config.keys():
		for interval in config["uavsar_data"].keys():
			filename = config["prep_inputs_dir"]+config["uavsar_data"][interval]["uav_textfile"];
			print("Reading %s " % filename);
			uavsar_input = slippy.io.read_insar_data(filename);
			config["uavsar_data"][interval]["data"] = uavsar_input;

	if "leveling_data" in config.keys():
		for interval in config["leveling_data"].keys():
			filename = config["prep_inputs_dir"]+config["leveling_data"][interval]["lev_outfile"];
			print("Reading %s " % filename)
			leveling_input = slippy.io.read_insar_data(filename);
			config["leveling_data"][interval]["data"] = leveling_input;

	if "tsx_data" in config.keys():
		for interval in config["tsx_data"].keys():
			filename = config["prep_inputs_dir"]+config["tsx_data"][interval]["tsx_datafile"];
			print("Reading %s " % filename)
			tsx_input = slippy.io.read_insar_data(filename);
			config["tsx_data"][interval]["data"] = tsx_input;

	return config;


if __name__=="__main__":
	config=welcome_and_parse(sys.argv);
	# data_input = read_in_data(config);
	buildG.beginning_calc(config);
	# print(data_input);
	# FROM CONFIG: 
	# BUILD G
	# INVERT
	# PROCESS OUTPUTS