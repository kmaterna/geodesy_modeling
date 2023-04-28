"""
Mathematical functions independent of object specifics
"""

import numpy as np
from Tectonic_Utils.geodesy import haversine


def get_nearest_pixel_in_raster(raster_lon, raster_lat, target_lon, target_lat, min_dist_cutoff_km=0.25):
    """
    Find grid location closest to target in 2d arrays with lat and lon.
    This could be optimized with a distance formula or walking algorithm in the future.

    :param raster_lon: 2d array
    :param raster_lat: 2d array
    :param target_lon: float
    :param target_lat: float
    :param min_dist_cutoff_km: minimum distance for error codes, in km
    :returns: x_idx, y_idx, and minimum distance (km)
    """
    dist = np.zeros(np.shape(raster_lon));
    lon_shape = np.shape(raster_lon);
    for i in range(lon_shape[0]):
        for j in range(lon_shape[1]):
            mypt = [raster_lat[i][j], raster_lon[i][j]];
            dist[i][j] = haversine.distance((target_lat, target_lon), mypt);
    minimum_distance_km = np.nanmin(dist);
    if minimum_distance_km < min_dist_cutoff_km:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
        i_found = idx[0][0];
        j_found = idx[1][0];
        print(raster_lon[i_found][j_found], raster_lat[i_found][j_found]);
    else:
        i_found, j_found = np.nan, np.nan;  # error codes
    return i_found, j_found, minimum_distance_km;


def get_nearest_pixel_in_geocoded_array(array_lon, array_lat, target_lon, target_lat, min_dist_cutoff=0.005):
    """
    Find the nearest pixel to a latitude/longitude point in a 2D geocoded array (doesn't do distance formula).
    vector_lon and vector_lat don't have to be the same length, just representing a 2D grid of points.
    Throws an error if the target point is outside the domain.

    :param array_lon: 1d array, longitudes.
    :param array_lat: 1d array, latitudes (doesn't have to be the same length as vector_lon).
    :param target_lon: float, longitude of the point of interest.
    :param target_lat: float, latitude of the point of interest
    :param min_dist_cutoff: allowable border of the geocoded array, in degrees
    :returns: idx_lon, idx_lat, distance to nearest pixel, in km
    """
    idx_lon = np.abs(array_lon - target_lon).argmin();
    idx_lat = np.abs(array_lat - target_lat).argmin();
    # Throw an error if the recovered location is really far from the target (i.e., off the arrays)
    if np.abs(target_lon - array_lon[idx_lon]) > min_dist_cutoff:
        print("refidx[0]-InSAR_obj.lon[idx_lon]: %f degrees " % np.abs(target_lon - array_lon[idx_lon]));
        raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (longitude)");
    if np.abs(target_lat - array_lat[idx_lat]) > min_dist_cutoff:
        print("refidx[1]-InSAR_obj.lon[idx_lat]: %f degrees " % np.abs(target_lat - array_lat[idx_lat]));
        raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (latitude)");
    distance_km = haversine.distance((target_lat, target_lon), (array_lat[idx_lat], array_lon[idx_lon]));
    return idx_lon, idx_lat, distance_km;


def get_nearest_pixels_in_list(tuple_coord_list, target_lon, target_lat, close_threshold=0.0009, min_dist_cutoff=0.003):
    """
    Find the element closest to the target location in vector of lon and vector of lat.
    Fast function because of numpy math.
    This implementation seems unnecessarily complex and possibly bad error handling. Might want to replace it.

    :param tuple_coord_list: list of matching tuples of (lon, lat)
    :param target_lon: float
    :param target_lat: float
    :param close_threshold: distance in degrees where a pixel is "close".  Default 0.0009 degrees.
    :param min_dist_cutoff: distance wherein the closest pixel is still too far to be useful. Default 0.003 degrees.
    :returns: idx of closest pixel, distance of minimum pixel, list idx of all pixels where
    distance is < a certain value.
    """
    vector_lon = [x[0] for x in tuple_coord_list];
    vector_lat = [x[1] for x in tuple_coord_list];
    dist = np.sqrt(np.power(np.subtract(vector_lon, target_lon), 2) + np.power(np.subtract(vector_lat, target_lat), 2));
    minimum_distance = np.nanmin(dist);
    close_pixels = np.where(dist < close_threshold);  # the runner-up close pixels, about 100 of them
    if minimum_distance < min_dist_cutoff:  # if we're inside the domain.
        idx = np.where(dist == np.nanmin(dist));
        i_found = idx[0][0];
    else:  # if we're outside the domain, return an error code.
        i_found = np.nan;
    return i_found, minimum_distance, close_pixels;


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
