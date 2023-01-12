"""
Utilities on 2D InSAR_Obj
InSAR_2D_Object = collections.namedtuple('InSAR_2D_Object',
['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime'])
LOS is in mm
"""

import numpy as np
from .class_model import InSAR_2D_Object
from .. import general_utils
from Tectonic_Utils.geodesy import insar_vector_functions


def impose_InSAR_bounding_box(_InSAR_obj, _bbox=(-180, 180, -90, 90)):
    """Impose a bounding box on some InSAR data. Not written yet. """
    return [];


def flip_los_sign(InSAR_obj):
    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=np.multiply(InSAR_obj.LOS, -1),
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;


def subtract_reference(InSAR_obj, refidx, tolerance=0.005):
    """
    Subtract the value of the reference pixel, essentially creating a referenced object.

    :param InSAR_obj: an InSAR object
    :param refidx: [int, int] for row/col or [float, float] for lon/lat
    :param tolerance: how close must the reference be to a viable pixel? In Degrees
    """
    if type(refidx[0]) is float:  # if refidx is lon, lat
        idx_lon, idx_lat = general_utils.get_nearest_pixel_in_vector(InSAR_obj.lon, InSAR_obj.lat, refidx[0],
                                                                     refidx[1], tolerance=tolerance);
        refvalue = InSAR_obj.LOS[idx_lon][idx_lat];
    else:  # if refidx is row, col
        refvalue = InSAR_2D_Object.LOS[refidx[0]][refidx[1]];
    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=np.subtract(InSAR_obj.LOS, refvalue),
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;


def get_look_vector_at_point(InSAR_obj, target_lon, target_lat):
    """"
    Return 3 component look vector, flight direction, and incidence angle at target location.
    Look vector is from ground to platform.

    :param InSAR_obj: an InSAR_2D_object
    :param target_lon: float
    :param target_lat: float
    :returns: E, N, U, flight, inc
    """
    # extract the look vector at a given spot:
    colnum = (np.abs(InSAR_obj.lon - target_lon)).argmin()
    rownum = (np.abs(InSAR_obj.lat - target_lat)).argmin()
    E = InSAR_obj.lkv_E[rownum][colnum];
    N = InSAR_obj.lkv_N[rownum][colnum];
    U = InSAR_obj.lkv_U[rownum][colnum];
    flight, inc = insar_vector_functions.look_vector2flight_incidence_angles(E, N, U);
    return E, N, U, flight, inc;


def get_incidence_grid(InSAR_obj):
    """
    Compute incidence angle across the grid, using the 3-component look vector.
    Currently not numpy-vectorized, so it takes a little while.
    """
    inc = np.zeros(np.shape(InSAR_obj.lkv_E));
    for y in range(len(InSAR_obj.lat)):
        for x in range(len(InSAR_obj.lon)):
            heading, inc_i = insar_vector_functions.look_vector2flight_incidence_angles(InSAR_obj.lkv_E[y][x],
                                                                                        InSAR_obj.lkv_N[y][x],
                                                                                        InSAR_obj.lkv_U[y][x]);
            inc[y][x] = inc_i;
    return inc;


def rewrap_InSAR(InSAR_obj, wavelength):
    """
    Take unwrapped LOS measurements (mm) and artificially wrap them around a certain radar wavelength (mm)

    :param InSAR_obj: 2D InSAR object
    :param wavelength: float, mm
    """
    rewrapped = wrap_float(InSAR_obj.LOS, wavelength);
    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=rewrapped,
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;


def wrap_float(def_meas, wavelength):
    """
    Wrap a float or array (already referenced to refpixel) by a given wavelength. Tested.

    :param def_meas: float or array
    :param wavelength: float (in same units as def_meas)
    """
    wrapped_phase = np.mod(def_meas, wavelength/2) * (4*np.pi / wavelength) - np.pi;
    return wrapped_phase;


def defensive_checks(InSAR_obj):
    """
    Check for array-size sanity in a newly created 2D InSAR object
    """
    (leny, lenx) = np.shape(InSAR_obj.LOS);
    if len(InSAR_obj.lon) != lenx:
        raise ValueError("length of InSAR_Obj lon array doesn't match shape of data.");
    if len(InSAR_obj.lat) != leny:
        raise ValueError("length of InSAR_Obj lat array doesn't match shape of data.");
    if np.shape(InSAR_obj.lkv_E) != np.shape(InSAR_obj.LOS):
        raise ValueError("Shape of InSAR data doesn't match shape of lkv arrays.");
    print("All arrays read in with size %s " % str(np.shape(InSAR_obj.LOS)));
    return;
