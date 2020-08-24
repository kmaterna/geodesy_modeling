# June 2020
# Remove a ramp, either naturally or by GPS

import numpy as np
from .class_model import InSAR_Object
from .inputs import inputs_txt
from .outputs import write_insar_invertible_format


def remove_ramp_filewise(insar_textfile, ramp_removed_file, ref_coord=None):
    # Here we will solve the least squares problem for the equation of a plane, and then remove it.
    # Then write out the data again.
    InSAR_Obj = inputs_txt(insar_textfile);
    noplane_Obj = remove_ramp_insarformat(InSAR_Obj, ref_coord);
    print("Writing ramp-removed data into file %s " % ramp_removed_file);
    write_insar_invertible_format(noplane_Obj, 0.0, ramp_removed_file);
    return;


def remove_ramp_insarformat(InSAR_Obj, ref_coord=None):
    # Plane equation: ax + by + c = z
    # Solving Ax = B
    Z = [];
    A = np.zeros((len(InSAR_Obj.lon), 3));
    for i in range(len(InSAR_Obj.lon)):
        A[i, :] = [InSAR_Obj.lon[i], InSAR_Obj.lat[i], 1];
        Z.append(InSAR_Obj.LOS[i]);
    model = np.linalg.lstsq(A, Z);
    model = model[0];

    # Removing the planar model
    new_disp = [];
    for i in range(len(InSAR_Obj.lon)):
        ramp_solution = model[0] * InSAR_Obj.lon[i] + model[1] * InSAR_Obj.lat[i] + model[2];
        new_disp.append(InSAR_Obj.LOS[i] - ramp_solution);

    # Re-reference if necessary
    if ref_coord is not None:
        ref_plane = model[0] * ref_coord[0] + model[1] * ref_coord[1] + model[2];
        new_disp = [x - ref_plane for x in new_disp];

    new_InSAR_Obj = InSAR_Object(lon=InSAR_Obj.lon, lat=InSAR_Obj.lat, LOS=new_disp, LOS_unc=InSAR_Obj.LOS_unc,
                                 lkv_E=InSAR_Obj.lkv_E, lkv_N=InSAR_Obj.lkv_N, lkv_U=InSAR_Obj.lkv_U,
                                 starttime=InSAR_Obj.starttime, endtime=InSAR_Obj.endtime);
    return new_InSAR_Obj;


def remove_ramp_with_GPS(insar_textfile, gps_textfile):
    return
