#!/usr/bin/env python
# July 2020
# Manipulate a UAVSAR file from JPL into a properly unwrapped and geocoded data file. 
# This is a complicated multi-step process that needs both ISCE and GMTSAR functions. 
# It writes the filtered, unwrapped, masked, geocoded interferograms in isce format, 
# and writes los.rdr.geo in isce format too.  Useful for Kite downsampling next. 

import numpy as np 
import matplotlib.pyplot as plt
import subprocess
import haversine
import jpl_uav_read_write
import netcdf_read_write
import isce_read_write
import mask_and_interpolate

def step1(data_file_slant, corr_file_slant, ann_file, filetype):
	real, imag = jpl_uav_read_write.read_igram_data(data_file_slant, ann_file, igram_type=filetype, return_type='real_imag') 
	corr = jpl_uav_read_write.read_corr_data(corr_file_slant, ann_file, igram_type=filetype) 
	return real, imag, corr;

def step2(real, imag, corr, cut_rowcol):
	# Cut imag, real, and corr arrays; write them in lots of formats. 
	real = real[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]];
	imag = imag[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]];
	corr = corr[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]];
	phase = np.arctan2(imag,real);
	complex_numbers = np.float32(real) + 1j*np.float32(imag);
	cor32 = np.float32(corr);
	xdata = range(0, np.shape(phase)[1]);
	ydata = range(0, np.shape(phase)[0]);
	(ny, nx) = np.shape(phase);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, 'radians', 'phase.grd', dtype=float);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, corr, 'corr', 'corr.grd', dtype=float);
	netcdf_read_write.produce_output_plot('phase.grd', 'Phase', 'phase.png', 'phase', aspect=1.0,invert_yaxis=False);
	netcdf_read_write.produce_output_plot('corr.grd', 'Coherence', 'corr.png', 'corr', aspect=1.0,cmap='binary_r',invert_yaxis=False);
	isce_read_write.write_isce_data(complex_numbers, nx, ny, 'CFLOAT', 'cut_slc.int');  # THIS HAS CHANGED
	isce_read_write.write_isce_data(cor32, nx, ny, 'FLOAT', 'cut_cor.cor');
	return;

def plotting_filtering(before_file, after_file):
	before_phase = isce_read_write.read_phase_data(before_file);
	after_phase = isce_read_write.read_phase_data(after_file);
	f, axarr = plt.subplots(1,2,figsize=(10,8),dpi=300);
	axarr[0].imshow(before_phase, cmap='rainbow');
	axarr[0].set_title('Before')
	axarr[1].imshow(after_phase, cmap='rainbow');
	axarr[1].set_title('After');
	plt.savefig("before_after_filtering.png");
	return;

def step3(after_filtering, after_filtering_corr, cutoff, nx, ny):
	# Multiply by coherence mask
	phase = isce_read_write.read_phase_data(after_filtering);
	corr = isce_read_write.read_scalar_data(after_filtering_corr);
	coherence_mask = mask_and_interpolate.make_coherence_mask(corr, cutoff);
	masked_phase = mask_and_interpolate.apply_coherence_mask(phase, coherence_mask);
	xdata = range(0, np.shape(phase)[1]);
	ydata = range(0, np.shape(phase)[0]);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, 'radians','phase_filtered.grd');
	netcdf_read_write.produce_output_netcdf(xdata, ydata, masked_phase, 'radians', 'phase_masked.grd', dtype=float);
	netcdf_read_write.produce_output_plot('phase_filtered.grd', 'Phase', 'phase_filtered.png', 'phase', aspect=1.0,invert_yaxis=False);
	netcdf_read_write.produce_output_plot('phase_masked.grd', 'Phase', 'phase_masked.png', 'phase', aspect=1.0,invert_yaxis=False);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, corr, 'corr', 'corr.grd', dtype=float);
	netcdf_read_write.produce_output_plot('corr.grd', 'Coherence', 'corr.png', 'corr', aspect=1.0,cmap='binary_r',invert_yaxis=False);		
	return masked_phase, coherence_mask;

def step4(phase): # interpolate 
	# Perform phase interpolation
	interp_array = mask_and_interpolate.interpolate_2d(phase);
	xdata = range(0, np.shape(phase)[1]);
	ydata = range(0, np.shape(phase)[0]);	
	netcdf_read_write.produce_output_netcdf(xdata, ydata, interp_array, 'radians', 'phase_interp.grd', dtype=float);	
	netcdf_read_write.produce_output_plot('phase_interp.grd', 'Phase', 'phase_interp.png', 'phase', aspect=1.0,invert_yaxis=False);	
	return interp_array;  

def step5(mask):  # re-apply the mask
	unw_grd = netcdf_read_write.read_netcdf4("unwrap.grd");
	unw_grd = np.multiply(unw_grd, mask);
	xdata = range(0, np.shape(unw_grd)[1]);
	ydata = range(0, np.shape(unw_grd)[0]);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, unw_grd, 'radians','unwrap_masked.grd');
	netcdf_read_write.produce_output_plot('unwrap_masked.grd', 'Unwrapped Phase', 'unw_masked.png', 'phase', aspect=1.0,invert_yaxis=False);
	return;

def step6(ann_file, cut_rowcol, looks_x, looks_y):  # geocode the ground pixels and multiply by wavelength
	num_rows, num_cols = jpl_uav_read_write.get_rows_cols(ann_file, 'ground');
	start_lon, start_lat, lon_inc, lat_inc = jpl_uav_read_write.get_ground_range_corner_increment(ann_file);
	x_orig = [start_lon + i*lon_inc for i in range(0,num_cols)];
	y_orig = [start_lat + i*lat_inc for i in range(0,num_rows)];
	x_cut = x_orig[cut_rowcol[2]: cut_rowcol[3]];
	y_cut = y_orig[cut_rowcol[0]: cut_rowcol[1]];  # implement the grid cut
	# next, implement the multilooking
	x_filt = []; y_filt = [];
	counter=np.arange(0,len(x_cut),looks_x)
	for i in range(len(counter)):
		region = np.mean(x_cut[counter[i]:counter[i]+looks_x])
		x_filt.append(region);
	counter=np.arange(0,len(y_cut),looks_y);
	for i in range(len(counter)):
		region = np.mean(y_cut[counter[i]:counter[i]+looks_y])
		y_filt.append(region);

	unw = netcdf_read_write.read_netcdf4('unwrap_masked.grd');

	plt.figure(figsize=(11,7),dpi=300)
	X,Y = np.meshgrid(x_filt, y_filt);
	plt.pcolormesh(X, Y, unw,cmap='jet',vmin=0, vmax=20);
	plt.colorbar();
	plt.savefig('unwrapped_geocoded_phase.png');	
	
	
	# CONVERT TO MM using the wavelength of UAVSAR	
	unw = np.multiply(unw, 237.9/2);   
	(ny, nx) = np.shape(unw);	

	# ISCE UNW.GEO
	isce_read_write.write_isce_data(unw, nx, ny, "FLOAT", "uavsar.unw.geo", 
		firstLat=max(y_filt), firstLon=min(x_filt), deltaLon=lon_inc*looks_x, deltaLat=lat_inc*looks_y,
		Xmin=min(x_filt), Xmax=max(x_filt)); # 1 band, floats

	return x_filt, y_filt;

def cross_track_pos(target_lon, target_lat, nearrange_lon, nearrange_lat, heading_cartesian):
	# Just get the cross-track position of a coordinate system centered at (nearrange_lon, nearrange_lat) with given heading 
	distance = haversine.distance((target_lat, target_lon),(nearrange_lat, nearrange_lon));
	compass_bearing = haversine.calculate_initial_compass_bearing((nearrange_lat, nearrange_lon),(target_lat, target_lon));  # this comes as CW from north
	theta = bearing_to_cartesian(compass_bearing);  # the angle of the position vector in cartesian coords
	# heading_cartesian is the angle between the east unit vector and the flight direction
	x0 = distance * np.cos(np.deg2rad(theta));
	y0 = distance * np.sin(np.deg2rad(theta));  # in the east-north coordinate systeem
	x_prime, y_prime = rotation_matrix(x0, y0, heading_cartesian);
	return y_prime;

def rotation_matrix(x0, y0, heading):
	# heading in degrees CCW from East, like mathematical definition. 
	x_prime = x0*np.cos(np.deg2rad(-heading)) - y0*np.sin(np.deg2rad(-heading));
	y_prime = x0*np.sin(np.deg2rad(-heading)) + y0*np.cos(np.deg2rad(-heading));
	return x_prime, y_prime;

def bearing_to_cartesian(heading):
	return 90-heading; 

def complement_angle(angle):
	return 90-angle;

def cartesian_to_ccw_from_north(angle):
	return angle-90; 

def incidence_angle(xtp, cross_track_max, near_inc_angle, far_inc_angle):
	# Using the incidence angles (to the vertical) at the upper and lower corners of the track, 
	# what's the incidence angle at some location in between (xtp=cross-track-position)? 
	# near_angle is the incidence angle between the viewing geometry and the vertical. 
	# nearrange is the complement of that angle. 
	# This function is kind of like linear interpolation, but a little bit curved
	nearrange = np.deg2rad(complement_angle(near_inc_angle));
	farrange = np.deg2rad(complement_angle(far_inc_angle)); # angles measured from the ground to the satellite
	h = (np.tan(nearrange)*np.tan(farrange)*cross_track_max) / (np.tan(nearrange)-np.tan(farrange)) ;
	angle_to_horizontal = np.rad2deg(np.arctan(h / (xtp + (h/np.tan(nearrange))) ));
	return complement_angle(angle_to_horizontal);

def step7(ann_file, x_filt, y_filt):  # Make los.rdr.geo
	near_angle, far_angle, heading =jpl_uav_read_write.get_nearrange_farrange_heading_angles(ann_file);
	print("Heading is %f degrees CW from north" % heading);
	heading_cartesian = bearing_to_cartesian(heading);  # CCW from east
	print("Cartesian Heading is %f" % heading_cartesian)
	# Get the upper and lower left corners, so we can compute the length of the across-track extent in km
	ul_lon, ul_lat, ll_lon, ll_lat = jpl_uav_read_write.get_ground_range_left_corners(ann_file);
	cross_track_max = haversine.distance((ll_lat, ll_lon),(ul_lat,ul_lon));  # in km

	# Get the azimuth angle for the pixels looking up to the airplane
	# My own documentation says CCW from north, even though that's really strange. 
	azimuth = heading_cartesian-90; # 90 degrees to the right of the airplane heading (for the look vector from ground to plane)
	azimuth = cartesian_to_ccw_from_north(azimuth);   # degrees CCW from North
	print("azimuth from ground to plane is:", azimuth)

	[X, Y] = np.meshgrid(x_filt, y_filt);
	(ny, nx) = np.shape(X);
	grid_az = azimuth * np.ones(np.shape(X));
	grid_inc = np.zeros(np.shape(X));

	print("Computing incidence angles for all pixels")
	for i in range(ny):
		for j in range(nx):
			xtp = cross_track_pos(X[i,j],Y[i,j],ll_lon,ll_lat, heading_cartesian);  # THIS WILL HAVE TO CHANGE FOR ASCENDING AND DESCENDING
			inc = incidence_angle(xtp, cross_track_max, near_angle, far_angle);
			grid_inc[i,j] = inc;
	
	# Finally, write the thing
	isce_read_write.write_isce_unw(grid_inc, grid_az, nx, ny, "FLOAT", 'los.rdr.geo');
	
	return;



# CONFIGURE
filetype = 'ground';  # options: ground, slant
ann_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.ann";
data_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.int.grd";  # 1 GB
corr_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.cor.grd";  # 500 Mb
after_filtering = "cut_filtered_slc.int"
after_filtering_corr = "cut_filtered_cor.cor"
cut_rowcol = [2500, 5100, 7800, 13000];  # how to cut the frame for unwrapping and filtering. 
cor_cutoff = 0.21;
looks_y = 5;
looks_x = 5; # nx and ny are 1040, 520 for this scheme	
ny = int((cut_rowcol[1]-cut_rowcol[0]) / looks_y);
nx = int((cut_rowcol[3]-cut_rowcol[2]) / looks_x);

# TO START, RUN THESE. BEGIN THE PROCESS WITH CUTTING AND WRITING. 
real, imag, corr = step1(data_file, corr_file, ann_file, filetype);
step2(real, imag, corr, cut_rowcol); # cut and write in many formats

# # FILTER STEP: Filter using looks.py.
subprocess.call(['looks.py','-i','cut_slc.int','-o',after_filtering,'-r',str(looks_x),'-a',str(looks_y)],shell=False);
subprocess.call(['looks.py','-i','cut_cor.cor','-o',after_filtering_corr,'-r',str(looks_x),'-a',str(looks_y)],shell=False);
plotting_filtering('cut_slc.int',after_filtering);

# # # MASK, INTERPOLATE, AND UNWRAP STEP
phase, mask = step3(after_filtering, after_filtering_corr, cor_cutoff, nx, ny); # multiply by coherence mask
phase = step4(phase); # perform phase interpolation
subprocess.call(['/Users/kmaterna/Documents/B_Research/Salton/Brawley_multiSAR_project/Code/custom_unwrap.sh'],shell=True); # THEN UNWRAP
re_masked = step5(mask);  # re-apply the mask and write

# # THEN, GEOCODE BASED ON ANNOTATION FILE, RANGES, AND MULTILOOKING PARAMETERS
x_filt, y_filt = step6(ann_file, cut_rowcol, looks_x, looks_y);

step7(ann_file, x_filt, y_filt);  # produce los.rdr.geo


