"""
Input functions for 2D InSAR-format data
"""

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from .class_model import InSAR_2D_Object


def inputs_grd(los_grdfile):
    """
    Input function for netcdf file.
    :param los_grdfile: string, filename
    :returns InSAR_Obj: InSAR_2D_Object
    """
    [lon, lat, LOS] = netcdf_read_write.read_any_grd(los_grdfile);
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                                lkv_E=None, lkv_N=None, lkv_U=None,
                                starttime=None, endtime=None);
    return InSAR_Obj;


def inputs_from_synthetic_enu_grids(e_grdfile, n_grdfile, u_grdfile, flight_angle, constant_incidence_angle=None):
    """
    For synthetic models with three deformation components calculated.
    If constant_incidence_angle is provided, it uses one simple incidence angle and flight angle for the field.
    Future: Options for non-constant incidence angle have not yet been written. Range depends on track/orbit.

    :param e_grdfile: string, filename
    :param n_grdfile: string, filename
    :param u_grdfile: string, filename
    :param flight_angle: float, flight angle, degrees cw from n
    :param constant_incidence_angle: float, incidence angle, degrees from vertical
    """
    [lon, lat, e] = netcdf_read_write.read_any_grd(e_grdfile);
    [_, _, n] = netcdf_read_write.read_any_grd(n_grdfile);
    [_, _, u] = netcdf_read_write.read_any_grd(u_grdfile);
    look_vector = insar_vector_functions.flight_incidence_angles2look_vector(flight_angle, constant_incidence_angle);

    lkv_E = np.multiply(np.ones(np.shape(e)), look_vector[0]);  # constant incidence angle for now
    lkv_N = np.multiply(np.ones(np.shape(e)), look_vector[1]);  # constant incidence angle for now
    lkv_U = np.multiply(np.ones(np.shape(e)), look_vector[2]);  # constant incidence angle for now
    los = insar_vector_functions.def3D_into_LOS(e, n, u, flight_angle, constant_incidence_angle);
    los = np.multiply(los, 1000);  # convert from m to mm
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=los, LOS_unc=np.zeros(np.shape(los)),
                                lkv_E=lkv_E, lkv_N=lkv_N, lkv_U=lkv_U, starttime=None, endtime=None);
    print("Done with reading object");
    return InSAR_Obj;
