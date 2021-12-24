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


def rewrap_InSAR(InSAR_obj, wavelength, refidx):
    """
    Take unwrapped LOS measurements (mm) and artificially wrap them around a certain radar wavelength
    :param InSAR_obj: 2D InSAR object
    :param wavelength: float, mm
    :param refidx: [float, float] representing row, col
    """
    def_meas = InSAR_obj.LOS - InSAR_obj.LOS[refidx[0]][refidx[1]];
    rewrapped = wrap_float(def_meas, wavelength);
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
