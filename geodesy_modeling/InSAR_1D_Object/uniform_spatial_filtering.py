"""
InSAR_Obj spatial filtering algorithm
These functions do not appear to be implemented or used yet.
"""


from . import class_model


def uniform_filtering(InSAR_obj, _spatial_wavelength_x, _spatial_wavelength_y):
    print("Spatial Filtering: Starting with %d points " % (len(InSAR_obj.LOS)))
    print("STOP! Filtering algorithm is not written yet!!!! ")

    LOS_filt = []  # Here we would do filtering.
    filt_InSAR_obj = class_model.Insar1dObject(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=LOS_filt,
                                               LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E,
                                               lkv_N=InSAR_obj.lkv_N,
                                               lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime,
                                               endtime=InSAR_obj.endtime)
    print("Done with filtering: Ending with %d points " % (len(LOS_filt)))
    return filt_InSAR_obj
