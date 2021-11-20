"""
March 2021
Definition of InSAR-format data
Input functions for InSAR-format data
"""

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from .class_model import InSAR_2D_Object


def inputs_grd(los_grdfile):
    """Input function for netcdf file"""
    [lon, lat, LOS] = netcdf_read_write.read_any_grd(los_grdfile);
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                                lkv_E=None, lkv_N=None, lkv_U=None,
                                starttime=None, endtime=None);
    return InSAR_Obj;


def inputs_enu_grids(e_grdfile, n_grdfile, u_grdfile, flight_angle, incidence_angle):
    """For synthetic models with three deformation components calculated.
    Uses a single incidence angle and flight angle right now."""
    [lon, lat, e] = netcdf_read_write.read_any_grd(e_grdfile);
    [_, _, n] = netcdf_read_write.read_any_grd(n_grdfile);
    [_, _, u] = netcdf_read_write.read_any_grd(u_grdfile);
    look_vector = insar_vector_functions.flight_incidence_angles2look_vector(flight_angle, incidence_angle);
    los = np.zeros(np.shape(e));
    [numrows, numcols] = np.shape(e);
    for i in range(numrows):
        for j in range(numcols):
            los[i][j] = 1000 * insar_vector_functions.def3D_into_LOS(e[i][j], n[i][j], u[i][j], flight_angle,
                                                                     incidence_angle);  # in mm
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=los, LOS_unc=np.zeros(np.shape(los)),
                                lkv_E=look_vector[0], lkv_N=look_vector[1], lkv_U=look_vector[2],
                                starttime=None, endtime=None);
    return InSAR_Obj;
