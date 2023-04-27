"""
June 2020
Remove a ramp or a constant from 1D InSAR object, either naturally or by GPS
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from .. import general_utils
from .class_model import InSAR_1D_Object
from .inputs import inputs_txt
from .outputs import write_insar_invertible_format


def remove_ramp_filewise(insar_textfile, ramp_removed_file, ref_coord=None):
    """
    Solve the least squares problem for the equation of a plane, and then remove it.
    Then write out the data again.
    """
    InSAR_Obj = inputs_txt(insar_textfile);
    noplane_Obj = remove_ramp(InSAR_Obj, ref_coord);
    print("Writing ramp-removed data into file %s " % ramp_removed_file);
    plotting_ramp_results(InSAR_Obj, noplane_Obj, insar_textfile+".png");
    write_insar_invertible_format(noplane_Obj, ramp_removed_file);
    return;


def remove_constant_filewise(insar_textfile, constant_removed_file, ref_coord=None):
    InSAR_Obj = inputs_txt(insar_textfile);
    noconst_Obj = remove_constant_insarformat(InSAR_Obj, ref_coord);
    print("Writing constant-removed data into file %s " % constant_removed_file);
    plotting_ramp_results(InSAR_Obj, noconst_Obj, insar_textfile+".png");
    write_insar_invertible_format(noconst_Obj, constant_removed_file);
    return;


def remove_constant_insarformat(InSAR_Obj, ref_coord=None):
    """
    Remove a constant from the InSAR_Obj.
    If ref_coord, then remove ref_coord.
    If not, then remove the median value.
    :param InSAR_Obj: 1D insar object
    :param ref_coord: [lon, lat] of point constrained to be zero.
    :returns: 1D insar object
    """
    if ref_coord:
        nearest_index, _, _ = general_utils.get_nearest_pixels_in_list(InSAR_Obj.get_coordinate_tuples(),
                                                                       ref_coord[0], ref_coord[1]);
        constant = InSAR_Obj.LOS[nearest_index];
    else:
        constant = np.nanmedian(InSAR_Obj.LOS);
    new_disp = [x - constant for x in InSAR_Obj.LOS];
    new_InSAR_Obj = InSAR_1D_Object(lon=InSAR_Obj.lon, lat=InSAR_Obj.lat, LOS=new_disp, LOS_unc=InSAR_Obj.LOS_unc,
                                    lkv_E=InSAR_Obj.lkv_E, lkv_N=InSAR_Obj.lkv_N, lkv_U=InSAR_Obj.lkv_U,
                                    starttime=InSAR_Obj.starttime, endtime=InSAR_Obj.endtime);
    return new_InSAR_Obj;


def remove_ramp(InSAR_Obj, ref_coord=None):
    """"
    Plane equation: ax + by + c = z
    Solving Ax = B
    We will re-reference if provided.
    Otherwise, we will remove the constant associated with the ramp.
    :param InSAR_Obj: 1D insar object
    :param ref_coord: [lon, lat] of point constrained to be zero.
    :returns: 1D insar object
    """
    nonan_obj = InSAR_Obj.remove_nans();
    Z = [];
    A = np.zeros((len(nonan_obj.lon), 3));
    for i in range(len(nonan_obj.lon)):
        A[i, :] = [nonan_obj.lon[i], nonan_obj.lat[i], 1];
        Z.append(nonan_obj.LOS[i]);
    model = np.linalg.lstsq(A, Z);
    model = model[0];

    # Removing the planar model
    new_disp = [];
    for i in range(len(InSAR_Obj.lon)):
        ramp_solution = model[0] * InSAR_Obj.lon[i] + model[1] * InSAR_Obj.lat[i] + model[2];
        new_disp.append(InSAR_Obj.LOS[i] - ramp_solution);

    # Re-reference if necessary
    if ref_coord:
        ref_plane = model[0] * ref_coord[0] + model[1] * ref_coord[1] + model[2];
        new_disp = [x - ref_plane for x in new_disp];

    new_InSAR_Obj = InSAR_1D_Object(lon=InSAR_Obj.lon, lat=InSAR_Obj.lat, LOS=new_disp, LOS_unc=InSAR_Obj.LOS_unc,
                                    lkv_E=InSAR_Obj.lkv_E, lkv_N=InSAR_Obj.lkv_N, lkv_U=InSAR_Obj.lkv_U,
                                    starttime=InSAR_Obj.starttime, endtime=InSAR_Obj.endtime);
    return new_InSAR_Obj;


# FUTURE WORK:
# def remove_ramp_with_GPS(insar_textfile, gps_textfile):
#     return


def plotting_ramp_results(Obj1, Obj2, filename):
    vmin = np.nanmin(Obj1.LOS);
    vmax = np.nanmax(Obj1.LOS);

    f, axarr = plt.subplots(1, 2, figsize=(12, 7), dpi=300);
    axarr[0].scatter(Obj1.lon, Obj1.lat, c=Obj1.LOS, cmap='rainbow', vmin=vmin, vmax=vmax);
    axarr[0].set_title('Object With Ramps');
    axarr[1].scatter(Obj2.lon, Obj2.lat, c=Obj2.LOS, cmap='rainbow', vmin=vmin, vmax=vmax);
    axarr[1].set_title('Object Without Ramps');

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(filename);
    return;
