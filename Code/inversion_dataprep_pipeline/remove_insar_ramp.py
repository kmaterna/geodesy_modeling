# June 2020
# Remove a ramp, either naturally or by GPS


import numpy as np 
import sys


def remove_ramp(insar_textfile, ramp_removed_file, data_segments = []):
	[lon_meas, lat_meas, disp, sig, unit_e, unit_n, unit_u] = np.loadtxt(insar_textfile, unpack=True, skiprows=1);
	# Here we will solve the least squares problem for the equation of a plane, and then remove it. 
	# Then write out the data again. 
	# LATER: Splitting up the data segments if the file contains more than one satellite. 

	# Plane equation: ax + by + c = z
	# Solving Ax = B
	Z = [];
	A = np.zeros((len(lon_meas),3));
	for i in range(len(lon_meas)):
		A[i,:] = [lon_meas[i], lat_meas[i], 1];
		Z.append(disp[i]);
	model = np.linalg.lstsq(A, Z);
	model = model[0];

	# Removing the planar model
	new_disp = [];
	for i in range(len(lon_meas)):
		ramp_solution = model[0]*lon_meas[i]+model[1]*lat_meas[i]+model[2];
		new_disp.append(disp[i] - ramp_solution);

	ofile=open(ramp_removed_file,'w');
	print("Writing ramp-removed data into file %s " % ramp_removed_file);
	for i in range(len(lon_meas)):
		ofile.write("%f %f %f %f %f %f %f\n" % (lon_meas[i], lat_meas[i], new_disp[i], sig[i], unit_e[i], unit_n[i], unit_u[i]) );
	ofile.close();

	return;

def remove_ramp_with_GPS(insar_textfile, gps_textfile):
	return;