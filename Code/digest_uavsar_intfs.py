#!/usr/bin/env python
# July 2020
# Manipulate a UAVSAR file from JPL into a properly unwrapped and geocoded data file. 

import numpy as np 
import matplotlib.pyplot as plt
import subprocess
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
	netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, 'radians', 'phase.grd', dtype=float);
	netcdf_read_write.produce_output_netcdf(xdata, ydata, corr, 'corr', 'corr.grd', dtype=float);
	netcdf_read_write.produce_output_plot('phase.grd', 'Phase', 'phase.png', 'phase', aspect=1.0,invert_yaxis=False);
	netcdf_read_write.produce_output_plot('corr.grd', 'Coherence', 'corr.png', 'corr', aspect=1.0,cmap='binary_r',invert_yaxis=False);
	isce_read_write.isce_igram_dump_without_xml(complex_numbers, "cut_slc.int");  # Can change on a computer with ISCE
	isce_read_write.isce_scalar_dump_without_xml(cor32, "cut_cor.cor");  # Can change on a computre with ISCE
	return;

def step3(after_filtering, after_filtering_corr, cutoff, nx, ny):
	# Multiply by coherence mask
	phase = isce_read_write.read_phase_data_no_isce(after_filtering, nx, ny);  # Can change this on a computer with ISCE
	corr = isce_read_write.read_scalar_data_no_isce(after_filtering_corr, nx, ny);  # Can change this on a computer with ISCE
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

def step4(phase):
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

def step6(ann_file, cut_rowcol, looks_x, looks_y):  # geocode the ground pixels
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

	unw = netcdf_read_write.read_netcdf4('geo_coords_results/unwrap.grd');
	netcdf_read_write.produce_output_netcdf(x_filt, y_filt, unw, "phase", "geo_coords_results/geo_unw.grd");
	plt.figure(figsize=(11,7),dpi=300)
	X,Y = np.meshgrid(x_filt, y_filt);
	plt.pcolormesh(X, Y, unw,cmap='jet',vmin=0, vmax=20);
	plt.colorbar();
	plt.savefig('unwrapped_geocoded.png');

	return x_filt, y_filt;


# CONFIGURE
filetype = 'ground';  # options: ground, slant
if filetype=='ground':
	ann_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.ann";
	data_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.int.grd";  # 1 GB
	corr_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.cor.grd";  # 500 Mb
	after_filtering = "Filtered/geo_coords/cut_filtered_slc.int"
	after_filtering_corr = "Filtered/geo_coords/cut_filtered_cor.cor"
	cut_rowcol = [2500, 5100, 7800, 13000];  # how to cut the frame for unwrapping and filtering. 
	cor_cutoff = 0.21;
	looks_y = 5;
	looks_x = 5; # nx and ny are 1040, 520 for this scheme	
	ny = int((cut_rowcol[1]-cut_rowcol[0]) / looks_y);
	nx = int((cut_rowcol[3]-cut_rowcol[2]) / looks_x);
else:  # slant coordinates
	ann_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.ann";
	data_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.int";  # 300 Mb
	corr_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.cor";
	after_filtering = "Filtered/radar_coords/cut_filtered_slc.int"
	after_filtering_corr = "Filtered/radar_coords/cut_filtered_cor.cor"
	cut_rowcol = [5500, 9000, 0, 1800];  # how to cut the frame for unwrapping and filtering.  Good for Brawley with 0.05 correlation  # cut_rowcol = [5500, 6000, 0, 500];  # test cut the frame for unwrapping and filtering. 
	cor_cutoff = 0.29;
	looks_azimuth = 7;
	looks_range = 5;  # nx and ny are 360, 500 for this scheme. 
	ny = int((cut_rowcol[1]-cut_rowcol[0]) / looks_azimuth);
	nx = int((cut_rowcol[3]-cut_rowcol[2]) / looks_range);


# TO START, RUN THESE. BEGIN THE PROCESS WITH CUTTING AND WRITING. 
# real, imag, corr = step1(data_file, corr_file, ann_file, filetype);
# step2(real, imag, corr, cut_rowcol); # cut and write in many formats

# THEN FILTER STEP: MUST MANUALLY TRANSFER .int and .cor to a computer with ISCE, and adjust the xml and vrt properly. 
# Filter using looks.py.
# Then return the files to the right machine. 

# # MASK, INTERPOLATE, AND UNWRAP STEP
# phase, mask = step3(after_filtering, after_filtering_corr, cor_cutoff, nx, ny); # multiply by coherence mask
# phase = step4(phase); # perform phase interpolation
# subprocess.call(['./custom_unwrap.sh'],shell=True); # THEN UNWRAP
# re_masked = step5(mask);  # re-apply the mask and write

# THEN, GEOCODE BASED ON ANNOTATION FILE, RANGES, AND MULTILOOKING PARAMETERS
x_filt, y_filt = step6(ann_file, cut_rowcol, looks_x, looks_y);




