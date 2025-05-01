"""
Input functions for 2D InSAR-format data
"""

import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions as insar_vect
from cubbie.read_write_insar_utilities import isce_read_write, jpl_uav_read_write
from .class_model import Insar2dObject


def inputs_grd(los_grdfile, _rdrlosfile=None):
    """
    Read netcdf file.

    :param los_grdfile: string, filename
    :param _rdrlosfile: string, filename, will be used for incidence and azimuth in the future
    :returns InSAR_Obj: InSAR_2D_Object
    """
    [lon, lat, LOS] = netcdf_read_write.read_any_grd(los_grdfile)
    holder = np.zeros(np.shape(LOS))

    # Here, will write an input function for when there's a corresponding look vector file in GMTSAR format.
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=holder, lkv_E=holder, lkv_N=holder, lkv_U=holder)
    return InSAR_Obj


def inputs_phase_isce(iscefile, los_rdr_file=None):
    """
    Create a 2D InSAR object from ISCE phase data.  Returns the scalar field in LOS, such as wrapped phase.
    """
    lon, lat = isce_read_write.get_xarray_yarray_from_xml(iscefile+'.xml')
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)))
    LOS = isce_read_write.read_phase_data(iscefile)
    if los_rdr_file:
        incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1)
        azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2)
    lkv_e, lkv_n, lkv_u = insar_vect.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence)
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                              lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None)
    return InSAR_Obj


def inputs_scalar_isce(iscefile, los_rdr_file=None):
    """
    Create a 2D InSAR object from ISCE data.  Returns the scalar field in LOS, such as unwrapped phase.
    """
    lon, lat = isce_read_write.get_xarray_yarray_from_xml(iscefile+'.xml')
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)))
    LOS = isce_read_write.read_scalar_data(iscefile)
    if los_rdr_file:
        incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1)
        azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2)
    lkv_e, lkv_n, lkv_u = insar_vect.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence)
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                              lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None)
    return InSAR_Obj


def inputs_from_synthetic_enu_grids(e_grdfile, n_grdfile, u_grdfile, flight_angle, constant_incidence_angle=35,
                                    convert_m_to_mm=True, look_direction='right'):
    """
    Read synthetic models with three deformation components.
    If constant_incidence_angle is provided, it uses one simple incidence angle and flight angle for the field.
    Future: Options for non-constant incidence angle have not yet been written. Range depends on track/orbit.

    :param e_grdfile: string, filename, with units of meters for default behavior
    :param n_grdfile: string, filename, with units of meters for default behavior
    :param u_grdfile: string, filename, with units of meters for default behavior
    :param flight_angle: float, flight angle, degrees CW from n
    :param constant_incidence_angle: float, incidence angle, degrees from vertical, default 35 degrees
    :param convert_m_to_mm: default True. Multiplies by 1000
    :param look_direction: default 'right'. Can take 'left'.
    """
    [lon, lat, e] = netcdf_read_write.read_any_grd(e_grdfile)
    [_, _, n] = netcdf_read_write.read_any_grd(n_grdfile)
    [_, _, u] = netcdf_read_write.read_any_grd(u_grdfile)

    look_vector = insar_vect.flight_incidence_angles2look_vector(flight_angle, constant_incidence_angle, look_direction)

    lkv_E = np.multiply(np.ones(np.shape(e)), look_vector[0])  # constant incidence angle for now, into 2D array
    lkv_N = np.multiply(np.ones(np.shape(e)), look_vector[1])  # constant incidence angle for now, into 2D array
    lkv_U = np.multiply(np.ones(np.shape(e)), look_vector[2])  # constant incidence angle for now, into 2D array

    los = insar_vect.def3D_into_LOS(e, n, u, flight_angle, constant_incidence_angle, look_direction)

    if convert_m_to_mm:
        los = np.multiply(los, 1000)  # convert from m to mm
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=los, LOS_unc=np.zeros(np.shape(los)),
                              lkv_E=lkv_E, lkv_N=lkv_N, lkv_U=lkv_U, starttime=None, endtime=None,
                              look_direction=look_direction)
    return InSAR_Obj


def inputs_from_uavsar_igrams(data_file, ann_file):
    """
    Read an interferogram from UAVSAR from the NASA UAVSAR format, complete with look vector information

    :param data_file: string, name of binary data file
    :param ann_file: string, name of text file with metadata/annotation
    :return: InSAR2D object
    """
    x, y, phase, _ = jpl_uav_read_write.read_igram_data(data_file, ann_file, igram_type='ground')
    inc, az = jpl_uav_read_write.read_los_rdr_geo_from_ground_ann_file(ann_file, x, y)  # takes 20 seconds
    e, n, u = insar_vect.flight_incidence_angles2look_vector(az, inc, 'left')
    # Read JPL UAVSAR data into 2d insar object
    mydata = Insar2dObject(x, y, phase, LOS_unc=np.zeros(np.shape(phase)), lkv_E=e, lkv_N=n, lkv_U=u,
                           look_direction='left')
    return mydata
