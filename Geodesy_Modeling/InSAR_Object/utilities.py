# July 2020
# Perform uniform downsampling on an InSAR_Obj
# Impose Bounding Box

import numpy as np
from .class_model import InSAR_Object


# InSAR Object is my standard format:
# InSAR_Object = collections.namedtuple('InSAR_Object',
# ['lon','lat','LOS','LOS_unc','lkv_E','lkv_N','lkv_U','starttime','endtime']);
# where LOS is in mm

def impose_InSAR_bounding_box(InSAR_obj, bbox=(-180, 180, -90, 90)):
    # Impose a bounding box on some InSAR data
    lon, lat, LOS = [], [], [];
    unit_E, unit_N, unit_U = [], [], [];
    LOS_unc = [];
    for i in range(len(InSAR_obj.lon)):
        if bbox[0] <= InSAR_obj.lon[i] <= bbox[1]:
            if bbox[2] <= InSAR_obj.lat[i] <= bbox[3]:
                if np.isnan(InSAR_obj.LOS[i]):
                    continue;
                else:
                    lon.append(InSAR_obj.lon[i]);
                    lat.append(InSAR_obj.lat[i]);
                    LOS.append(InSAR_obj.LOS[i]);
                    LOS_unc.append(InSAR_obj.LOS_unc[i]);
                    unit_E.append(InSAR_obj.lkv_E[i]);
                    unit_N.append(InSAR_obj.lkv_N[i]);
                    unit_U.append(InSAR_obj.lkv_U[i]);
    newInSAR_obj = InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N, lkv_U=unit_U,
                                starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return newInSAR_obj;


def remove_nans(InSAR_obj):
    # Remove Nans from some InSAR object
    lon, lat, LOS = [], [], [];
    unit_E, unit_N, unit_U = [], [], [];
    LOS_unc = [];
    for i in range(len(InSAR_obj.lon)):
        if np.isnan(InSAR_obj.LOS[i]):
            continue;
        else:
            lon.append(InSAR_obj.lon[i]);
            lat.append(InSAR_obj.lat[i]);
            LOS.append(InSAR_obj.LOS[i]);
            LOS_unc.append(InSAR_obj.LOS_unc[i]);
            unit_E.append(InSAR_obj.lkv_E[i]);
            unit_N.append(InSAR_obj.lkv_N[i]);
            unit_U.append(InSAR_obj.lkv_U[i]);
    newInSAR_obj = InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N, lkv_U=unit_U,
                                starttime=InSAR_obj.starttime, endtime=InSAR_obj.endtime);
    return newInSAR_obj;


def combine_objects(Obj1, Obj2):
    lon = np.hstack((Obj1.lon, Obj2.lon));
    lat = np.hstack((Obj1.lat, Obj2.lat));
    LOS = np.hstack((Obj1.LOS, Obj2.LOS));
    LOS_unc = np.hstack((Obj1.LOS_unc, Obj2.LOS_unc));
    lkv_E = np.hstack((Obj1.lkv_E, Obj2.lkv_E));
    lkv_N = np.hstack((Obj1.lkv_N, Obj2.lkv_N));
    lkv_U = np.hstack((Obj1.lkv_U, Obj2.lkv_U));
    newInSAR_obj = InSAR_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=lkv_E, lkv_N=lkv_N, lkv_U=lkv_U,
                                starttime=Obj1.starttime, endtime=Obj1.endtime);
    return newInSAR_obj;
