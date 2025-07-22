"""
Input functions for 2D InSAR-format data
"""

import numpy as np
from tectonic_utils.read_write import netcdf_read_write
from tectonic_utils.geodesy import insar_vector_functions as insar_vect
from cubbie.read_write_insar_utilities import isce_read_write, jpl_uav_read_write
from .class_model import Insar2dObject


def inputs_grd(los_grdfile, lkvE_file=None, lkvN_file=None, lkvU_file=None, coh_file=None):
    """
    Read netcdf file.

    :param los_grdfile: string, filename of grd file
    :param lkvE_file: string, filename of grd file
    :param lkvN_file: string, filename of grd file
    :param lkvU_file: string, filename of grd file
    :param coh_file: string, filename of grd file
    :returns InSAR_Obj: InSAR_2D_Object
    """
    [lon, lat, LOS] = netcdf_read_write.read_any_grd(los_grdfile)
    lkvE, lkvN, lkvU = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)))
    coh = np.ones((len(lat), len(lon)))
    # Here, we read corresponding look vector information in GMTSAR format.
    if lkvE_file:
        _, _, lkvE = netcdf_read_write.read_netcdf4(lkvE_file)
    if lkvN_file:
        _, _, lkvN = netcdf_read_write.read_netcdf4(lkvN_file)
    if lkvU_file:
        _, _, lkvU = netcdf_read_write.read_netcdf4(lkvU_file)
    if coh_file:
        _, _, coh = netcdf_read_write.read_netcdf4(coh_file)
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                              lkv_E=lkvE, lkv_N=lkvN, lkv_U=lkvU, coherence=coh, starttime=None, endtime=None)
    return InSAR_Obj


def inputs_phase_isce(iscefile, los_rdr_file=None, coh_file=None):
    """
    Create a 2D InSAR object from ISCE phase data.  Returns the scalar field in LOS, such as wrapped phase.

    :param iscefile: string, filename
    :param los_rdr_file: string, filename, ISCE los.rdr file
    :param coh_file: string, filename, ISCE coherence file
    :returns: InSAR_2D_object
    """
    lon, lat, LOS = isce_read_write.read_phase_data(iscefile)
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)))
    coh = np.ones((len(lat), len(lon)))
    if los_rdr_file:
        _, _, incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1)
        _, _, azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2)
    if coh_file:
        _, _, coh = isce_read_write.read_scalar_data(coh_file, band=1)
    if lat[1] < lat[0]:
        lat = np.flipud(lat)
        LOS = np.flipud(LOS)
        incidence = np.flipud(incidence)
        azimuth = np.flipud(azimuth)
        coh = np.flipud(coh)
    lkv_e, lkv_n, lkv_u = insar_vect.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence)
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                              lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None,
                              coherence=coh)
    return InSAR_Obj


def inputs_scalar_isce(iscefile, los_rdr_file=None, coh_file=None):
    """
    Create a 2D InSAR object from ISCE data.  Returns the scalar field in LOS, such as unwrapped phase.

    :param iscefile: string, filename
    :param los_rdr_file: string, filename, ISCE los.rdr file
    :param coh_file: string, filename, ISCE coherence file
    :returns: InSAR_2D_object
    """
    lon, lat, LOS = isce_read_write.read_scalar_data(iscefile)
    incidence, azimuth = np.empty((len(lat), len(lon))), np.empty((len(lat), len(lon)))
    coh = np.ones((len(lat), len(lon)))
    if los_rdr_file:
        _, _, incidence = isce_read_write.read_scalar_data(los_rdr_file, band=1)
        _, _, azimuth = isce_read_write.read_scalar_data(los_rdr_file, band=2)
    if coh_file:
        _, _, coh = isce_read_write.read_scalar_data(coh_file, band=1)
    if lat[1] < lat[0]:
        lat = np.flipud(lat)
        LOS = np.flipud(LOS)
        incidence = np.flipud(incidence)
        azimuth = np.flipud(azimuth)
        coh = np.flipud(coh)
    lkv_e, lkv_n, lkv_u = insar_vect.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence)
    InSAR_Obj = Insar2dObject(lon=lon, lat=lat, LOS=LOS, LOS_unc=np.zeros(np.shape(LOS)),
                              lkv_E=lkv_e, lkv_N=lkv_n, lkv_U=lkv_u, starttime=None, endtime=None,
                              coherence=coh)
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


def inputs_from_uavsar_igrams(data_file, ann_file, corr_file=None):
    """
    Read an interferogram from UAVSAR from the NASA UAVSAR format, complete with look vector information

    :param data_file: string, name of binary data file
    :param ann_file: string, name of text file with metadata/annotation
    :param corr_file: string, optional, filename of binary coherence file
    :return: InSAR2D object
    """
    # Read JPL UAVSAR data into 2d insar object
    x, y, phase, _ = jpl_uav_read_write.read_igram_data(data_file, ann_file, igram_type='ground')
    inc, az = jpl_uav_read_write.read_los_rdr_geo_from_ground_ann_file(ann_file, x, y)  # takes 20 seconds
    if corr_file:
        coh = jpl_uav_read_write.read_corr_data(corr_file, ann_file)
    else:
        coh = None
    if y[1] < y[0]:
        y = np.flipud(y)
        phase = np.flipud(phase)
        inc = np.flipud(inc)
        az = np.flipud(az)
        coh = np.flipud(coh)
    e, n, u = insar_vect.flight_incidence_angles2look_vector(az, inc, 'left')
    mydata = Insar2dObject(x, y, phase, LOS_unc=np.zeros(np.shape(phase)), lkv_E=e, lkv_N=n, lkv_U=u,
                           look_direction='left', coherence=coh)
    return mydata
