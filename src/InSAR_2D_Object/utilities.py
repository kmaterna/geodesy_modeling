"""
Utilities on 2D InSAR_Obj
InSAR_2D_Object = collections.namedtuple('InSAR_2D_Object',
['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime'])
LOS is in mm
"""

import numpy as np
from .class_model import InSAR_2D_Object


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
    if type(refidx[0]) is float:
        # refidx is lon, lat
        target_lon, target_lat = refidx[0], refidx[1];
        idx_lon = np.abs(InSAR_obj.lon - target_lon).argmin();
        idx_lat = np.abs(InSAR_obj.lat - target_lat).argmin();
        # Throw an error if the recovered location is really far from the target (i.e., off the arrays)
        if np.abs(refidx[0]-InSAR_obj.lon[idx_lon]) > tolerance:
            print("refidx[0]-InSAR_obj.lon[idx_lon]: %f degrees " % np.abs(refidx[0]-InSAR_obj.lon[idx_lon]));
            raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (longitude)");
        if np.abs(refidx[1]-InSAR_obj.lat[idx_lat]) > tolerance:
            print("refidx[1]-InSAR_obj.lon[idx_lat]: %f degrees " % np.abs(refidx[1] - InSAR_obj.lat[idx_lat]));
            raise ValueError("Error! Resolved Reference Pixel is not within tol of Target Reference Pixel (latitude)");
        refvalue = InSAR_obj.LOS[idx_lon][idx_lat];
    else:
        # refidx is row, col
        row, col = refidx[0], refidx[1];
        refvalue = InSAR_2D_Object.LOS[row][col];
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
