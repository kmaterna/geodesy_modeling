# Some useful functions for comparing pixels in leveling with 
# pixels in other things (i.e. uavsar, tsx, s1)
# This file could also be called pixel utilities

import numpy as np
from Tectonic_Utils.geodesy import haversine


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat):
    # Take a raster and find the grid location closest to the target location
    dist = np.zeros(np.shape(raster_lon));
    lon_shape = np.shape(raster_lon);
    for i in range(lon_shape[0]):
        for j in range(lon_shape[1]):
            mypt = [raster_lat[i][j], raster_lon[i][j]];
            dist[i][j] = haversine.distance((target_lat, target_lon), mypt);
    minimum_distance = np.nanmin(dist);
    if minimum_distance < 0.25:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
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
        dist[i] = haversine.distance((target_lat, target_lon), mypt);
    minimum_distance = np.nanmin(dist);
    if minimum_distance < 0.25:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
        i_found = idx[0][0];
    else:
        i_found = -1;  # error codes
    return i_found, minimum_distance;


def find_leveling_in_vector(myLev, vector_data):
    # Get the index for each leveling benchmark.
    print("Finding target leveling pixels in vector of data");
    vector_index = [];
    for bm in range(len(myLev.lat)):
        # name = myLev.name[bm].split()[0];
        i_found, mindist = get_nearest_pixel_in_vector(vector_data.lon, vector_data.lat, myLev.lon[bm], myLev.lat[bm]);
        vector_index.append(i_found);
    return vector_index;


def get_file_dictionary(config_filename):
    # GET FILE NAMES
    this_dict = {};
    ifile = open(config_filename);
    for line in ifile:
        data_type = line.split(":")[0];
        total_data_files = line.split(":")[1][1:-1];
        this_dict[data_type] = total_data_files;
    ifile.close();
    return this_dict;


def write_gps_invertible_format(gps_object_list, filename):
    # One header line
    # GPS Data in meters
    print("Writing GPS displacements into file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# Header: lon, lat, dE, dN, dU, Se, Sn, Su (m)\n");
    for station in gps_object_list:
        if np.isnan(station.dE[1]):
            continue;
        else:
            ofile.write('%f %f ' % (station.coords[0], station.coords[1]));
            ofile.write("%f %f %f " % (0.001 * station.dE[1], 0.001 * station.dN[1], 0.001 * station.dU[1]));
            ofile.write("%f %f %f\n" % (station.Se[1], station.Sn[1], station.Su[1]));
    ofile.close();
    return;