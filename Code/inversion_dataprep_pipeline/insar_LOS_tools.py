# July 2020
# Perform uniform downsampling on an InSAR_Obj
# Impose Bounding Box
# Convert between InSAR data formats

import numpy as np 
import sys
import multiSAR_input_functions
import multiSAR_utilities


# InSAR Object is similar to the Hines format: 
# InSAR_Object = collections.namedtuple('InSAR_Object',['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime']);
# where LOS is in mm


def uniform_downsampling(InSAR_obj, sampling_interval, averaging_window=0):
	print("Uniform downsampling: Starting with %d points " % (len(InSAR_obj.lon)) );

	# Step 1: Create uniform downsampled arrays
	x_array = np.arange(np.min(InSAR_obj.lon), np.max(InSAR_obj.lon), sampling_interval);
	y_array = np.arange(np.min(InSAR_obj.lat), np.max(InSAR_obj.lat), sampling_interval);
	[X, Y] = np.meshgrid(x_array, y_array);
	new_obs_array = np.zeros(np.shape(X)); new_obs_unc = np.zeros(np.shape(X));
	new_e = np.zeros(np.shape(X)); new_n = np.zeros(np.shape(X)); new_u = np.zeros(np.shape(X));
	if len(x_array)*len(y_array)>len(InSAR_obj.lon):
		# Defensive programming
		print("ERROR!  Trying to uniformly downsample but the number of pixels actually increases.  Try again!");
		return InSAR_obj;

	# Step 2: Populate uniform arrays
	for i in range(len(y_array)):
		for j in range(len(x_array)):
			if averaging_window==0:  # If we just want to find THE nearest pixel
				idx, min_dist = multiSAR_utilities.get_nearest_pixel_in_vector(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i]);
				if min_dist < sampling_interval*110:  # rough degrees to km conversion
					new_obs_array[i][j] = InSAR_obj.LOS[idx];
					new_obs_unc[i][j] = InSAR_obj.LOS_unc[idx];
					new_e[i][j] = InSAR_obj.lkv_E[idx];
					new_n[i][j] = InSAR_obj.lkv_N[idx];
					new_u[i][j] = InSAR_obj.lkv_U[idx];
				else:  # the nearest pixel was too far away
					new_obs_array[i][j] = np.nan;
					new_obs_unc[i][j] = np.nan;
					new_e[i][j] = np.nan;
					new_n[i][j] = np.nan;
					new_u[i][j] = np.nan;
			else:  # If we want to average over a spatial window
				new_obs_array[i][j] = multiSAR_utilities.get_average_within_box(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i], averaging_window, InSAR_obj.LOS);
				new_obs_unc[i][j] = multiSAR_utilities.get_average_within_box(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i], averaging_window, InSAR_obj.LOS_unc);
				new_e[i][j] = multiSAR_utilities.get_average_within_box(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i], averaging_window, InSAR_obj.lkv_E);
				new_n[i][j] = multiSAR_utilities.get_average_within_box(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i], averaging_window, InSAR_obj.lkv_N);
				new_u[i][j] = multiSAR_utilities.get_average_within_box(InSAR_obj.lon, InSAR_obj.lat, x_array[j], y_array[i], averaging_window, InSAR_obj.lkv_U);

	ds_lon = np.reshape(X, (len(x_array)*len(y_array),));
	ds_lat = np.reshape(Y, (len(x_array)*len(y_array),));
	ds_LOS = np.reshape(new_obs_array, (len(x_array)*len(y_array),));
	ds_LOS_unc = np.reshape(new_obs_unc, (len(x_array)*len(y_array),));
	ds_lkv_e = np.reshape(new_e, (len(x_array)*len(y_array),));
	ds_lkv_n = np.reshape(new_n, (len(x_array)*len(y_array),));
	ds_lkv_u = np.reshape(new_u, (len(x_array)*len(y_array),));

	ds_InSAR_obj = multiSAR_input_functions.InSAR_Object(lon=ds_lon, lat=ds_lat, LOS=ds_LOS, LOS_unc=InSAR_obj.LOS_unc, 
		lkv_E=ds_lkv_e, lkv_N=ds_lkv_n, lkv_U=ds_lkv_u, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
	print("Done with downsampling: Ending with %d points " % (len(ds_lon)) );
	return ds_InSAR_obj;



def impose_InSAR_bounding_box(InSAR_obj, bbox=[-180, 180, -90, 90]):
	# Impose a bounding box on some InSAR data
	lon=[]; lat=[]; LOS=[]; LOS_unc=[]; unit_E=[]; unit_N=[]; unit_U=[];
	for i in range(len(InSAR_obj.lon)):
		if InSAR_obj.lon[i]>=bbox[0] and InSAR_obj.lon[i]<=bbox[1]:
			if InSAR_obj.lat[i]>=bbox[2] and InSAR_obj.lat[i]<=bbox[3]:
				if np.isnan(InSAR_obj.LOS[i]):
					continue;
				else:
					lon.append(InSAR_obj.lon[i]);
					lat.append(InSAR_obj.lat[i]);
					LOS.append(InSAR_obj.LOS[i]);
					LOS_unc.append(InSAR_obj.LOS_unc[i]);
					unit_E.append(InSAR_obj.lkv_E[i]);
					unit_N.append(InSAR_obj.lkv_N[i]);
					unit_U.append(InSAR_obj.lkv_U[i]);
	newInSAR_obj = multiSAR_input_functions.InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N, lkv_U=unit_U, 
		starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
	return newInSAR_obj;



def remove_nans(InSAR_Object):
	# Remove Nans from some InSAR object
	lon=[]; lat=[]; LOS=[]; LOS_unc=[]; unit_E=[]; unit_N=[]; unit_U=[];
	for i in range(len(InSAR_obj.lon)):
		if np.isnan(InSAR_obj.LOS[i]):
			continue;
		else:
			lon.append(InSAR_obj.lon[i]);
			lat.append(InSAR_obj.lat[i]);
			LOS.append(InSAR_obj.LOS[i]);
			LOS_unc.append(InSAR_obj.LOS_unc[i]);
			unit_E.append(InSAR_obj.lkv_E[i]);
			unit_N.append(InSAR_obj.lkv_N[i]);
			unit_U.append(InSAR_obj.lkv_U[i]);
	newInSAR_obj = multiSAR_input_functions.InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N, lkv_U=unit_U, 
		starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
	return newInSAR_obj;

