# Some useful functions for comparing pixels in leveling with 
# pixels in other things (i.e. uavsar, tsx, s1)
# This file could also be called pixel utilities

import numpy as np 
import haversine
import sys


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat):
	# Take a raster and find the grid location closest to the target location
	dist = np.zeros(np.shape(raster_lon));
	lon_shape = np.shape(raster_lon);
	for i in range(lon_shape[0]):
		for j in range(lon_shape[1]):
			mypt = [raster_lat[i][j], raster_lon[i][j]];
			dist[i][j]=haversine.distance((target_lat, target_lon), mypt);
	minimum_distance = np.nanmin(dist);	
	if minimum_distance<0.25:  # if we're inside the domain.
		idx = np.where(dist==np.nanmin(dist));
		i_found = idx[0][0];
		j_found = idx[1][0];
		print(raster_lon[i_found][j_found], raster_lat[i_found][j_found]);
	else:
		i_found = -1;
		j_found = -1;  # error codes
	return i_found, j_found, minimum_distance;

def get_nearest_pixel_in_vector(vector_lon, vector_lat, target_lon, target_lat):
	# Take a vector and find the location closest to the target location
	dist = np.zeros(np.shape(vector_lon));
	for i in range(len(vector_lon)):
		mypt = [vector_lat[i], vector_lon[i]];
		dist[i]=haversine.distance((target_lat, target_lon), mypt);
	minimum_distance = np.nanmin(dist);	
	if minimum_distance<0.25:  # if we're inside the domain.
		idx = np.where(dist==np.nanmin(dist));
		i_found = idx[0][0];
	else:
		i_found = -1;  # error codes
	return i_found, minimum_distance;	


def find_leveling_in_vector(myLev, vector_data):
	# Get the index for each leveling benchmark. 
	vector_index=[];
	for bm in range(len(myLev.lat)):
		name = myLev.name[bm].split()[0];
		i_found, mindist = get_nearest_pixel_in_vector(vector_data.lon, vector_data.lat, myLev.lon[bm], myLev.lat[bm]);
		vector_index.append(i_found);
	return vector_index;


def get_average_within_box(lonlist, latlist, target_lon, target_lat, averaging_window, data):
	# averaging window in degrees
	# We search the averaging window in both directions. 
	new_data = [];
	for i in range(len(lonlist)):
		if lonlist[i]>=target_lon-averaging_window and lonlist[i]<=target_lon+averaging_window:
			if latlist[i]>=target_lat-averaging_window and latlist[i]<=target_lat+averaging_window:
				new_data.append(data[i]);
	return np.nanmean(new_data); 


