
import numpy as np
from geodesy_modeling.datatypes.InSAR_1D_Object import class_model


def avg_uavsar_disp(TS, slicenum, row, col):
    """Average around a few pixels"""
    width_pixels = 10
    return np.nanmean(TS[slicenum, row - width_pixels:row + width_pixels, col - width_pixels:col + width_pixels])


def get_onetime_displacements(myGridTS, start_idx, end_idx):
    """ Turns a GridTS object into an InSAR_1D_Object.
    Turning a raster into a vector in the process."""
    los_raster = np.subtract(myGridTS.TS[end_idx][:, :], myGridTS.TS[start_idx][:, :])
    los_vector = np.reshape(los_raster, (np.size(los_raster),))
    lon_vector = np.reshape(myGridTS.lon, (np.size(los_raster),))
    lat_vector = np.reshape(myGridTS.lat, (np.size(los_raster),))
    InSAR_obj = class_model.Insar1dObject(lon=lon_vector, lat=lat_vector, LOS=los_vector,
                                          LOS_unc=np.nan*np.zeros(np.shape(los_vector)),
                                          lkv_E=np.nan*np.zeros(np.shape(los_vector)),
                                          lkv_N=np.nan*np.zeros(np.shape(los_vector)),
                                          lkv_U=np.nan*np.zeros(np.shape(los_vector)),
                                          starttime=myGridTS.dtarray[start_idx],
                                          endtime=myGridTS.dtarray[end_idx])
    return InSAR_obj
