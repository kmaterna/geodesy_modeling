"""
Mathematical functions independent of object specifics
"""

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


def get_nearest_pixel_in_vector(vector_lon, vector_lat, target_lon, target_lat, tolerance=0.005):
    """
    A simpler implementation of finding the nearest pixel to a latitude/longitude point (doesn't do distance formula)
    Assumes a 2D grid of points
    Throws an error if the target point is outside the domain.

    :param vector_lon: 1d array, longitudes
    :param vector_lat: 1d array, latitudes
    :param target_lon: longitude of the point of interest
    :param target_lat: latitude of the point of interest
    :param tolerance: tolerance in degrees
    """
    idx_lon = np.abs(vector_lon - target_lon).argmin();
    idx_lat = np.abs(vector_lat - target_lat).argmin();
    # Throw an error if the recovered location is really far from the target (i.e., off the arrays)
    if np.abs(target_lon - vector_lon[idx_lon]) > tolerance:
        print("refidx[0]-InSAR_obj.lon[idx_lon]: %f degrees " % np.abs(target_lon - vector_lon[idx_lon]));
        raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (longitude)");
    if np.abs(target_lat - vector_lat[idx_lat]) > tolerance:
        print("refidx[1]-InSAR_obj.lon[idx_lat]: %f degrees " % np.abs(target_lat - vector_lat[idx_lat]));
        raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (latitude)");
    return idx_lon, idx_lat;


def get_many_nearest_pixels_in_vector(vector_lon, vector_lat, target_lon, target_lat):
    """
    Find the element closest to the target location in vector of lon and vector of lat.
    Fast function because of numpy math.

    This implementation seems unnecessarily complex and possibly bad error handling. Might want to replace it.
    """
    dist = np.sqrt(np.power(np.subtract(vector_lon, target_lon), 2) + np.power(np.subtract(vector_lat, target_lat), 2));
    minimum_distance = np.nanmin(dist);
    close_pixels = np.where(dist < 0.0009);  # the runner-up close pixels, about 100 of them
    if minimum_distance < 0.003:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
        i_found = idx[0][0];
    else:  # if we're outside the domain, return an error code.
        i_found = np.nan;
    return i_found, minimum_distance, close_pixels;


def find_pixels_idxs_in_ll_arrays(vector_lons, vector_lats, target_lons, target_lats):
    """
    Get the nearest index (and neighbors) for each given coordinate.

    :param vector_lons: 1D array of lon attributes
    :param vector_lats: 1D array of lat attributes
    :param target_lons: 1d array of longitudes to be found
    :param target_lats: 1d array of latitudes to be found
    :returns closest_index: a list that matches the length of InSAR_Data.lon
    :returns close_indices: a list of lists (might return the 100 closest pixels for a given coordinate)
    """
    print("Finding target leveling pixels in vector of data");
    closest_index, close_indices = [], [];
    for tlon, tlat in zip(target_lons, target_lats):
        i_found, mindist, closer_pts = get_many_nearest_pixels_in_vector(vector_lons, vector_lats, tlon, tlat);
        if i_found != -1:
            closest_index.append(i_found);
            close_indices.append(closer_pts);
        else:
            closest_index.append(np.nan);
            close_indices.append(np.nan);
    return closest_index, close_indices;


def compute_difference_metrics_on_same_pixels(list1, list2):
    """
    :param list1: a 1d array of LOS data from platform 1 (like Leveling)
    :param list2: a matching 1d array of LOS data from platform 2 (like UAVSAR)
    :returns: average misfit value, and r^2 coefficient.
    """
    misfit_metric = np.nanmean(np.abs(np.subtract(list1, list2)));  # average deviation from 1-to-1
    corr_matrix = np.corrcoef(list1, list2)
    corr = corr_matrix[0, 1]
    r2 = corr ** 2
    return misfit_metric, r2;


def convert_rates_to_disps(LOS_rates, starttime, endtime):
    """
    Compute displacement = rate * time

    :param LOS_rates: a vector
    :param starttime: a datetime object
    :param endtime: a datetime object
    :returns: a vector
    """
    tdelta = endtime - starttime;
    interval_years = tdelta.days / 365.24;  # number of years spanned by given velocity.
    Disps = [i * interval_years for i in LOS_rates];
    return Disps;


def convert_disps_to_rates(disps, starttime, endtime):
    """
    Compute displacement / time = rate

    :param disps: a vector
    :param starttime: a datetime object
    :param endtime: a datetime object
    :returns: a vector
    """
    tdelta = endtime - starttime;
    interval_years = tdelta.days / 365.24;  # number of years spanned by given velocity.
    LOS_rates = [i / interval_years for i in disps];
    return LOS_rates;


def get_file_dictionary(config_filename):
    """GET FILE NAMES"""
    this_dict = {};
    print("Reading file %s " % config_filename);
    ifile = open(config_filename);
    for line in ifile:
        data_type = line.split(':')[0];
        total_data_files = line.split()[1];  # assuming one file per list entry
        this_dict[data_type] = total_data_files;
    ifile.close();
    return this_dict;


def wrap_float(def_meas, wavelength):
    """
    Wrap a float or array (already referenced to refpixel) by a given wavelength. Tested.

    :param def_meas: float or array
    :param wavelength: float (in same units as def_meas)
    """
    wrapped_phase = np.mod(def_meas, wavelength/2) * (4*np.pi / wavelength) - np.pi;
    return wrapped_phase;