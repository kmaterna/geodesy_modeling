# Some useful functions for comparing pixels in leveling with 
# pixels in other things (i.e. uavsar, tsx, s1)
# This file could also be called pixel utilities

import numpy as np
from Tectonic_Utils.geodesy import haversine


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat):
    """Take a raster (2d arrays with lat and lon)
    and find the grid location closest to the target location
    This could be optimized with a distance formula or walking algorithm in the future.
    """
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
        i_found, j_found = np.nan, np.nan;  # error codes
    return i_found, j_found, minimum_distance;


def get_nearest_pixel_in_vector(vector_lon, vector_lat, target_lon, target_lat):
    """Take a vector and find the location closest to the target location.
    Fast function because of numpy math"""
    dist = np.sqrt(np.power(np.subtract(vector_lon, target_lon), 2) + np.power(np.subtract(vector_lat, target_lat), 2));
    minimum_distance = np.nanmin(dist);
    close_pixels = np.where(dist<0.0009);  # the runner-up close pixels, about 100 of them
    if minimum_distance < 0.003:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
        i_found = idx[0][0];
    else:  # if we're outside the domain, return an error code.
        i_found = np.nan;
    return i_found, minimum_distance, close_pixels;


def find_pixels_idxs_in_InSAR_Obj(InSAR_Data, target_lons, target_lats):
    """
    Get the nearest index (and neighbors) for each given coordinate.
    InSAR_Data : insar object, or any object with 1D lists of lon/lat attributes
    closest_index : a list that matches the length of InSAR_Data.lon
    close_indices : a list of lists (might return the 100 closest pixels for a given coordinate)
    """
    print("Finding target leveling pixels in vector of data");
    closest_index, close_indices = [], [];
    for tlon, tlat in zip(target_lons, target_lats):
        i_found, mindist, closer_pts = get_nearest_pixel_in_vector(InSAR_Data.lon, InSAR_Data.lat, tlon, tlat);
        if i_found != -1:
            closest_index.append(i_found);
            close_indices.append(closer_pts);
        else:
            closest_index.append(np.nan);
            close_indices.append(np.nan);
    return closest_index, close_indices;


def get_file_dictionary(config_filename):
    """GET FILE NAMES"""
    this_dict = {};
    ifile = open(config_filename);
    for line in ifile:
        data_type = line.split(":")[0];
        total_data_files = line.split(":")[1][1:-1];
        this_dict[data_type] = total_data_files;
    ifile.close();
    return this_dict;


def write_gps_invertible_format(gps_object_list, filename):
    """One header line
    GPS Data in meters"""
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
