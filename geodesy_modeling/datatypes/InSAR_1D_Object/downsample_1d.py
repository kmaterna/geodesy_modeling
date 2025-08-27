
import numpy as np
from geodesy_modeling import general_utils
from . import class_model


def uniform_downsampling(InSAR_obj, sampling_interval, averaging_window=0):
    """
    Produce a regular grid of points, but don't keep pixels with no data / NaNs.

    :param InSAR_obj: an InSAR_1D_Object with 1D columns of data
    :param sampling_interval: degrees, float
    :param averaging_window: degrees, float
    """
    print("Uniform downsampling: Starting with %d points " % (len(InSAR_obj.lon)))

    # Step 1: Create uniform downsampled arrays
    x_array = np.arange(np.min(InSAR_obj.lon), np.max(InSAR_obj.lon), sampling_interval)
    y_array = np.arange(np.min(InSAR_obj.lat), np.max(InSAR_obj.lat), sampling_interval)
    [X, Y] = np.meshgrid(x_array, y_array)
    new_obs_array = np.zeros(np.shape(X))
    new_obs_unc = np.zeros(np.shape(X))
    new_e = np.zeros(np.shape(X))
    new_n = np.zeros(np.shape(X))
    new_u = np.zeros(np.shape(X))
    if len(x_array) * len(y_array) > len(InSAR_obj.lon):
        # Defensive programming
        print("ERROR!  Trying to uniformly downsample but the number of pixels actually increases.  Try again!")
        return InSAR_obj

    # Step 2: Populate uniform arrays
    for i in range(len(y_array)):
        for j in range(len(x_array)):
            if averaging_window == 0:  # If we just want to find THE nearest pixel
                idx, min_dist, _ = general_utils.get_nearest_pixels_in_list(InSAR_obj.get_coordinate_tuples(),
                                                                            x_array[j], y_array[i])
                if min_dist < sampling_interval * 110:  # rough degrees to km conversion
                    new_obs_array[i][j] = InSAR_obj.LOS[idx]
                    new_obs_unc[i][j] = InSAR_obj.LOS_unc[idx]
                    new_e[i][j] = InSAR_obj.lkv_E[idx]
                    new_n[i][j] = InSAR_obj.lkv_N[idx]
                    new_u[i][j] = InSAR_obj.lkv_U[idx]
                else:  # the nearest pixel was too far away
                    new_obs_array[i][j] = np.nan
                    new_obs_unc[i][j] = np.nan
                    new_e[i][j] = np.nan
                    new_n[i][j] = np.nan
                    new_u[i][j] = np.nan
            else:  # If we want to average over a spatial window
                avg_pixel = InSAR_obj.get_average_los_within_box(x_array[j], y_array[i], averaging_window)
                new_obs_array[i][j] = float(avg_pixel.LOS)
                new_obs_unc[i][j] = float(avg_pixel.LOS_unc)  # mathematically, not sure how to handle average unc
                new_e[i][j] = float(avg_pixel.lkv_E)
                new_n[i][j] = float(avg_pixel.lkv_N)
                new_u[i][j] = float(avg_pixel.lkv_U)

    ds_lon = np.reshape(X, (len(x_array) * len(y_array),))
    ds_lat = np.reshape(Y, (len(x_array) * len(y_array),))
    ds_LOS = np.reshape(new_obs_array, (len(x_array) * len(y_array),))
    ds_LOS_unc = np.reshape(new_obs_unc, (len(x_array) * len(y_array),))
    ds_lkv_e = np.reshape(new_e, (len(x_array) * len(y_array),))
    ds_lkv_n = np.reshape(new_n, (len(x_array) * len(y_array),))
    ds_lkv_u = np.reshape(new_u, (len(x_array) * len(y_array),))

    ds_InSAR_obj = class_model.Insar1dObject(lon=ds_lon, lat=ds_lat, LOS=ds_LOS, LOS_unc=ds_LOS_unc,
                                             lkv_E=ds_lkv_e, lkv_N=ds_lkv_n, lkv_U=ds_lkv_u,
                                             starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime)
    ds_InSAR_obj = ds_InSAR_obj.remove_nans(verbose=True)
    print("Done with downsampling: Ending with %d points " % (len(ds_InSAR_obj.lon)))
    return ds_InSAR_obj
