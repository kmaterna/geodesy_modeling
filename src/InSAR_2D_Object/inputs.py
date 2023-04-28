"""
Input functions for 2D InSAR-format data
"""

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from S1_batches.read_write_insar_utilities import isce_read_write
from .class_model import InSAR_2D_Object


def inputs_grd(los_grdfile, _rdrlosfile=None):
    """
    Read netcdf file.

    :param los_grdfile: string, filename
    :param _rdrlosfile: string, filename
    :returns InSAR_Obj: InSAR_2D_Object
    """
    [lon, lat, LOS] = netcdf_read_write.read_any_grd(los_grdfile);

    # Here, will write an input function for when there's a corresponding look vector file in GMTSAR format.
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                                lkv_E=None, lkv_N=None, lkv_U=None,
                                starttime=None, endtime=None);
    return InSAR_Obj;


def inputs_phase_isce(iscefile, los_rdr_file=None):
    """
    Create a 2D InSAR object from ISCE phase data.  Returns the scalar field in LOS, such as wrapped phase.
    """
    lon, lat = isce_read_write.get_xarray_yarray_from_xml(iscefile+'.xml');
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)));
    LOS = isce_read_write.read_phase_data(iscefile);
    if los_rdr_file:
        incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1);
        azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2);
    lkv_e, lkv_n, lkv_u = insar_vector_functions.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence);
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                                lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None);
    InSAR_Obj.defensive_checks();
    return InSAR_Obj;


def inputs_scalar_isce(iscefile, los_rdr_file=None):
    """
    Create a 2D InSAR object from ISCE data.  Returns the scalar field in LOS, such as unwrapped phase.
    """
    lon, lat = isce_read_write.get_xarray_yarray_from_xml(iscefile+'.xml');
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)));
    LOS = isce_read_write.read_scalar_data(iscefile);
    if los_rdr_file:
        incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1);
        azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2);
    lkv_e, lkv_n, lkv_u = insar_vector_functions.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence);
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                                lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None);
    InSAR_Obj.defensive_checks();
    return InSAR_Obj;


def inputs_from_synthetic_enu_grids(e_grdfile, n_grdfile, u_grdfile, flight_angle, constant_incidence_angle=None,
                                    convert_m_to_mm=True):
    """
    Read synthetic models with three deformation components.
    If constant_incidence_angle is provided, it uses one simple incidence angle and flight angle for the field.
    Future: Options for non-constant incidence angle have not yet been written. Range depends on track/orbit.

    :param e_grdfile: string, filename, with units of meters for default behavior
    :param n_grdfile: string, filename, with units of meters for default behavior
    :param u_grdfile: string, filename, with units of meters for default behavior
    :param flight_angle: float, flight angle, degrees cw from n
    :param constant_incidence_angle: float, incidence angle, degrees from vertical
    :param convert_m_to_mm: default True. Multiplies by 1000
    """
    [lon, lat, e] = netcdf_read_write.read_any_grd(e_grdfile);
    [_, _, n] = netcdf_read_write.read_any_grd(n_grdfile);
    [_, _, u] = netcdf_read_write.read_any_grd(u_grdfile);
    look_vector = insar_vector_functions.flight_incidence_angles2look_vector(flight_angle, constant_incidence_angle);

    lkv_E = np.multiply(np.ones(np.shape(e)), look_vector[0]);  # constant incidence angle for now, into 2D array
    lkv_N = np.multiply(np.ones(np.shape(e)), look_vector[1]);  # constant incidence angle for now, into 2D array
    lkv_U = np.multiply(np.ones(np.shape(e)), look_vector[2]);  # constant incidence angle for now, into 2D array
    los = insar_vector_functions.def3D_into_LOS(e, n, u, flight_angle, constant_incidence_angle);
    if convert_m_to_mm:
        los = np.multiply(los, 1000);  # convert from m to mm
    InSAR_Obj = InSAR_2D_Object(lon=lon, lat=lat, LOS=los, LOS_unc=np.zeros(np.shape(los)),
                                lkv_E=lkv_E, lkv_N=lkv_N, lkv_U=lkv_U, starttime=None, endtime=None);
    InSAR_Obj.defensive_checks();
    return InSAR_Obj;
