# InSAR_Obj spatial filtering algorithm
# These functions do not appear to be implemented or used yet.


import numpy as np
from InSAR_Object import class_model


# InSAR Object is my standard format:
# InSAR_Object = collections.namedtuple('InSAR_Object',['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U',
#                                                       'starttime','endtime']);
# where LOS is in mm


def uniform_downsampling(InSAR_obj, spatial_wavelength_x, spatial_wavelength_y):
    print("Spatial Filtering: Starting with %d points " % (len(InSAR_obj.LOS)));
    print("STOP! Filtering algorithm is not written yet!!!! ");

    LOS_filt = [];  # Here we would do filtering.
    filt_InSAR_obj = class_model.InSAR_Object(lon=InSAR_obj.lon, lat=InSAR_obj.lat, LOS=LOS_filt,
                                              LOS_unc=InSAR_obj.LOS_unc, lkv_E=InSAR_obj.lkv_E, lkv_N=InSAR_obj.lkv_N,
                                              lkv_U=InSAR_obj.lkv_U, starttime=InSAR_obj.starttime,
                                              endtime=InSAR_obj.endtime);
    print("Done with filtering: Ending with %d points " % (len(LOS_filt)));
    return filt_InSAR_obj;


def get_average_within_box(lonlist, latlist, target_lon, target_lat, averaging_window, data):
    # averaging window in degrees
    # We search the averaging window in both directions from the target loc, and average the data
    new_data = [];
    for i in range(len(lonlist)):
        if target_lon - averaging_window <= lonlist[i] <= target_lon + averaging_window:
            if target_lat - averaging_window <= latlist[i] <= target_lat + averaging_window:
                new_data.append(data[i]);
    return np.nanmean(new_data);
