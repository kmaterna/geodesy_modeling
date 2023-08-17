"""
Mathematical functions independent of object specifics
"""

import numpy as np
import scipy
from Tectonic_Utils.geodesy import haversine, insar_vector_functions


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


def detrend_signal(x, y):
    """
    Remove an overall trend from a 1D timeseries signal, ignoring nans
    """
    x_nonans = x[np.where(~np.isnan(y))];
    y_nonans = y[np.where(~np.isnan(y))];
    A = np.vstack((x_nonans, np.ones(np.shape(x_nonans)))).T;
    coeffs = scipy.linalg.lstsq(np.array(A), np.array(y_nonans));  # the actual optimization step
    model_coeffs = coeffs[0];  # model: [z = ax + c]
    full_model = model_coeffs[0]*x + model_coeffs[1] * np.ones(np.shape(x));
    return y-full_model;


def wrap_float(def_meas, wavelength):
    """
    Wrap a float or array (already referenced to refpixel) by a given wavelength. Tested.

    :param def_meas: float or array
    :param wavelength: float (in same units as def_meas)
    """
    wrapped_phase = np.mod(def_meas, wavelength/2) * (4*np.pi / wavelength) - np.pi;
    return wrapped_phase;


def get_los_and_flight_vectors(flight_heading, look_dir='right'):
    """
    Return the unit vectors pointing in the flight direction and look direction

    :param flight_heading: degrees CW from north.
    :param look_dir: string, 'right' or 'left'.
    """
    [x_flight, y_flight] = insar_vector_functions.get_unit_vector_from_heading(flight_heading);
    if look_dir == 'right':
        [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading + 90);
    elif look_dir == 'left':
        [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading - 90);
    else:
        raise ValueError("Error! Invalid look_dir; look_dir must be right or left.");
    return x_flight, y_flight, x_los, y_los;
