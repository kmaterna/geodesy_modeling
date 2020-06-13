# After you've done an inversion, what are the results? 
# How much moment is created, and how big is the misfit? 


import numpy as np 
import sys
import json
import moment_calculations
import slippy.io

def parse_json(argv):
	print("Metrics for this inversion:")
	if len(argv)<2:
		print("Error! Please provide the name of a config json. Exiting. "); sys.exit(0);
	else:
		configname = argv[1];
	config_file = open(configname,'r');
	config = json.load(config_file);
	for key in config:  # adding the output directory onto the output files, for ease of use. 
		if "output" in key and key != "output_dir":
			config[key]=config["output_dir"]+config[key];
	return config;

def get_slip_moment(slip_filename):
	# From the inversion results, what is the moment? 
	moment_total = 0;
	mu=30e9;  # Pa
	length, width, leftlat, thrust, tensile = np.loadtxt(slip_filename,skiprows=1,unpack=True,usecols=(5,6,7,8,9));
	for i in range(len(length)):
		slip = np.sqrt(leftlat[i]*leftlat[i] + thrust[i]*thrust[i]);
		area = length[i]*width[i]; # m^2
		momenti = moment_calculations.moment_from_muad(mu, area, slip);
		moment_total = moment_total+momenti;
	mw = moment_calculations.mw_from_moment(moment_total);
	scinot = "{:e}".format(moment_total);
	print("From file %s: Total Slip Moment is %s N-m, equivalent to mw=%f \n" % (slip_filename, scinot, mw) );
	return moment_total, mw;

def get_total_misfit(config):
	if "leveling_input_file" in config.keys():
		get_misfit_leveling(config["leveling_input_file"],config["leveling_output_file"]);
	if "insar_input_file" in config.keys():
		get_misfit_insar(config["insar_input_file"],config["insar_output_file"]);
	if "gps_input_file" in config.keys():
		get_misfit_gps(config["gps_input_file"],config["gps_output_file"]);
	return;

def get_misfit_gps(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	gps_input = slippy.io.read_gps_data(obs_file); 
	gps_pred = slippy.io.read_gps_data(pred_file); 
	abs_misfit = np.abs(gps_input[1]-gps_pred[1]);
	norm_misfit = np.divide(abs_misfit, gps_input[2]);  # divide by sigma
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	print("Average GPS misfit: %f mm" % (1000*mean_average_misfit) );
	print("Average normalized GPS misfit: %f sigma \n" % (mean_norm_average_misfit) );
	return;

def get_misfit_insar(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	insar_input = slippy.io.read_insar_data(obs_file)
	insar_pred = slippy.io.read_insar_data(pred_file)
	abs_misfit = np.abs(insar_input[1]-insar_pred[1]);
	norm_misfit = np.divide(abs_misfit, insar_input[2]);  # divide by sigma	
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	print("Average InSAR misfit: %f mm" % (1000*mean_average_misfit) );
	print("Average normalized InSAR misfit: %f sigma \n" % (mean_norm_average_misfit) );
	return;

def get_misfit_leveling(obs_file, pred_file):
	# Misfit from each data pair (GPS, UAVSAR, Leveling, S1, TSX)
	# Want in both absolute numbers and relative to the respective uncertainties. 
	leveling_input = slippy.io.read_insar_data(obs_file)
	leveling_pred = slippy.io.read_insar_data(pred_file)
	abs_misfit = np.abs(leveling_input[1]-leveling_pred[1]);
	norm_misfit = np.divide(abs_misfit, leveling_input[2]);  # divide by sigma	
	mean_average_misfit = np.nanmean(abs_misfit);
	mean_norm_average_misfit = np.nanmean(norm_misfit);
	print("Average Leveling misfit: %f mm" % (1000*mean_average_misfit) );
	print("Average normalized Leveling misfit: %f sigma \n" % (mean_norm_average_misfit) );	
	return;


def write_metrics_report(config):
	# Write the smoothing parameters, the misfits, and the slip moment into a file. 
	# This will eventually be useful for L-curves
	return;


if __name__=="__main__":
	config=parse_json(sys.argv);
	get_slip_moment(config["slip_output_file"]);
	get_total_misfit(config);




