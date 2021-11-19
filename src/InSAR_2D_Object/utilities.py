# July 2020
# Perform uniform downsampling on an InSAR_Obj
# Impose Bounding Box

import numpy as np
from .class_model import InSAR_2D_Object


# InSAR Object is my standard format:
# InSAR_2D_Object = collections.namedtuple('InSAR_2D_Object',
# ['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime']);
# where LOS is in mm

def impose_InSAR_bounding_box(InSAR_obj, bbox=(-180, 180, -90, 90)):
    """Impose a bounding box on some InSAR data. Not written yet. """
    return [];


def flip_los_sign(InSAR_obj):
    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=np.multiply(InSAR_obj.LOS, -1),
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;

def rewrap_InSAR(InSAR_obj, wavelength, refidx):
    """ Take unwrapped LOS measurements (mm) and artificially wrap them around a certain radar wavelength"""
    [numrows, numcols] = np.shape(InSAR_obj.LOS);
    rewrapped = np.zeros(np.shape(InSAR_obj.LOS));
    for i in range(numrows):
        for j in range(numcols):
            def_meas = InSAR_obj.LOS[i][j] - InSAR_obj.LOS[refidx[0]][refidx[1]];
            rewrapped[i][j] = wrap_float(def_meas, wavelength);

    new_InSAR_obj = InSAR_2D_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=rewrapped,
                                    LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                    lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return new_InSAR_obj;

def wrap_float(def_meas, wavelength):
    """Wrap a float by a given wavelength. Tested - seems to work."""
    wrapped_phase = np.mod(def_meas, wavelength/2) * (4*np.pi / wavelength) - np.pi;
    return wrapped_phase;
