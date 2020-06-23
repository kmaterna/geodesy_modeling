# June 2020
# Run Slippy with multiple time intervals in the input 
# (A really big G matrix)


import numpy as np 
import sys
import json
import buildG

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



if __name__=="__main__":
	config=welcome_and_parse(sys.argv);
	buildG.beginning_calc(config);
