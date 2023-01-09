"""
Utilities on 2D InSAR_Obj
InSAR_2D_Object = collections.namedtuple('InSAR_2D_Object',
['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime'])
LOS is in mm
"""

import numpy as np
from .class_model import InSAR_2D_Object
from .. import multiSAR_utilities


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
        idx_lon, idx_lat = multiSAR_utilities.get_nearest_pixel_in_vector(InSAR_obj.lon, InSAR_obj.lat, refidx[0],
                                                                          refidx[1], tolerance=tolerance);
        refvalue = InSAR_obj.LOS[idx_lon][idx_lat];
    else:  # if refidx is row, col
        refvalue = InSAR_2D_Object.LOS[refidx[0]][refidx[1]];
    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=np.subtract(InSAR_obj.LOS, refvalue),
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;


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
