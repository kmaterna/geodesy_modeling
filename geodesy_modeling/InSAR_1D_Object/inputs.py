"""
Input functions for 1D InSAR format
"""

import numpy as np
import datetime as dt
import sys
import pandas
import h5py
from cubbie.read_write_insar_utilities import isce_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from ..general_utils import convert_rates_to_disps
from .class_model import Insar1dObject
from elastic_stresses_py.PyCoulomb.fault_slip_triangle.file_io import io_other


def inputs_txt(insar_textfile, starttime=dt.datetime.strptime("19900101", "%Y%m%d"),
               endtime=dt.datetime.strptime("19900101", "%Y%m%d")):
    """
    Read slippy-invertible text file.
    Columns are: lon, lat, disp(m), sig(m), unit_e, unit_n, unit_u, with one header row.

    :param insar_textfile: string, file name
    :param starttime: optional, beginning of InSAR interval, dt.datetime object
    :param endtime: optional, end of InSAR interval, dt.datetime object
    """
    print("Reading file %s " % insar_textfile)
    [lon_meas, lat_meas, disp, sig, unit_e, unit_n, unit_u] = np.loadtxt(insar_textfile, unpack=True, skiprows=1)
    InSAR_Obj = Insar1dObject(lon=lon_meas, lat=lat_meas, LOS=disp * 1000, LOS_unc=sig * 1000, lkv_E=unit_e,
                              lkv_N=unit_n, lkv_U=unit_u, starttime=starttime, endtime=endtime)
    return InSAR_Obj


def inputs_simplest_txt(insar_textfile, starttime=dt.datetime.strptime("19900101", "%Y%m%d"),
                        endtime=dt.datetime.strptime("19900101", "%Y%m%d")):
    """
    Read text files with lon, lat, LOS (mm), with one header row.
    """
    [lon_meas, lat_meas, disp] = np.loadtxt(insar_textfile, unpack=True, skiprows=1)
    InSAR_Obj = Insar1dObject(lon=lon_meas, lat=lat_meas, LOS=disp,
                              LOS_unc=np.multiply(np.nan, np.zeros(np.shape(lon_meas))),
                              lkv_E=np.multiply(np.nan, np.zeros(np.shape(lon_meas))),
                              lkv_N=np.multiply(np.nan, np.zeros(np.shape(lon_meas))),
                              lkv_U=np.multiply(np.nan, np.zeros(np.shape(lon_meas))),
                              starttime=starttime, endtime=endtime)
    return InSAR_Obj


def inputs_TRE_vert_east(filename):
    """
    Read InSAR data from TRE data excel spreadsheets.
    """
    print("Reading in %s" % filename)

    # Vertical velocities
    vert_sheet = pandas.read_excel(filename, engine='openpyxl', sheet_name=2)
    zvel = vert_sheet['VEL'].values.tolist()
    zvel_std = vert_sheet['VEL_STD'].values.tolist()
    lat = vert_sheet['LAT'].values.tolist()
    lon = vert_sheet['LON'].values.tolist()

    # East velocities
    east_sheet = pandas.read_excel(filename, engine='openpyxl', sheet_name=3)
    evel = east_sheet['VEL'].values.tolist()
    evel_std = east_sheet['VEL_STD'].values.tolist()

    if "TSX" in filename:
        starttime = dt.datetime.strptime("2012-09-01", "%Y-%m-%d")
        endtime = dt.datetime.strptime("2013-09-01", "%Y-%m-%d")  # Hard coded for this experiment.
    elif "SNT1" in filename:
        starttime = dt.datetime.strptime("2015-04-01", "%Y-%m-%d")
        endtime = dt.datetime.strptime("2018-04-01", "%Y-%m-%d")  # Hard coded for this experiment.
    elif "SNT2" in filename:
        starttime = dt.datetime.strptime("2018-05-01", "%Y-%m-%d")
        endtime = dt.datetime.strptime("2019-08-01", "%Y-%m-%d")  # Hard coded for this experiment.
    elif "ENV_2005" in filename:
        starttime = dt.datetime.strptime("2005-08-01", "%Y-%m-%d")  # for Heber experiment, vertical/east
        endtime = dt.datetime.strptime("2010-08-01", "%Y-%m-%d")  # for Heber experiment, vertical/east
        # Heber dataset:
        # ascending: 12-2003 to 08-2010 (but denser data starts at 08-2005)
        # descending: 02-2003 to 09-2010
        # vertical/east : 08-2005 to 08-2010
    else:
        print("Unrecognized TRE filename option! Cannot find start and end times.  Exiting...")
        sys.exit(0)

    # Multiply by number of years of observation package it up.
    Vert_LOS = convert_rates_to_disps(zvel, starttime, endtime)
    East_LOS = convert_rates_to_disps(evel, starttime, endtime)
    zeros = np.zeros(np.shape(zvel))
    ones = np.ones(np.shape(zvel))
    Vert_obj = Insar1dObject(lon=lon, lat=lat, LOS=Vert_LOS, LOS_unc=zvel_std,
                             lkv_E=zeros, lkv_N=zeros, lkv_U=ones, starttime=starttime, endtime=endtime)
    East_obj = Insar1dObject(lon=lon, lat=lat, LOS=East_LOS, LOS_unc=evel_std,
                             lkv_E=ones, lkv_N=zeros, lkv_U=zeros, starttime=starttime, endtime=endtime)
    return Vert_obj, East_obj


def inputs_cornell_ou_velocities_hdf5(filename, lkv_filename, slicenum=0):
    """
    Read HDF5 file from Junle Jiang data format, one track at a time.
    HDF5 file contains 5 velocity measurements for each pixel, one at each time interval.
    Will return whichever slicenum (0-4) is supplied
    """
    print("Reading %s " % filename)
    f = h5py.File(filename, 'r')  # reads the hdf5 into a dictionary f

    # Unpacking
    dates = np.array(f.get('dates'))
    print(dates)
    lat = np.array(f.get('lat'))
    lon = np.array(f.get('lon'))
    rate = np.array(f.get('rate'))
    rstd = np.array(f.get('rstd'))

    f2 = h5py.File(lkv_filename, 'r')
    los = f2.get('los')
    lkv_E = los[0, :, :]
    lkv_N = los[1, :, :]
    lkv_U = los[2, :, :]

    num_data = np.shape(lon)[0] * np.shape(lon)[1]
    lon = np.reshape(lon, (num_data,))
    lat = np.reshape(lat, (num_data,))
    lkv_E = np.reshape(lkv_E, (num_data,))
    lkv_N = np.reshape(lkv_N, (num_data,))
    lkv_U = np.reshape(lkv_U, (num_data,))
    unc = np.reshape(rstd[:, :, slicenum], (num_data,))
    disps, dates = quick_convert_one_timeslice_to_disp(rate[:, :, slicenum], dates[slicenum])

    # Returning standard InSAR format of displacements in mm, etc.
    InSAR_data = Insar1dObject(lon=lon, lat=lat, LOS=disps, LOS_unc=unc, lkv_E=lkv_E, lkv_N=lkv_N, lkv_U=lkv_U,
                               starttime=dates[0], endtime=dates[1])
    return InSAR_data


def inputs_isce_unw_geo_losrdr(isce_unw_filename, los_filename, starttime=dt.datetime.strptime("19900101", "%Y%m%d"),
                               endtime=dt.datetime.strptime("19900101", "%Y%m%d")):
    """
    Read geocoded unw file and associated rdr.los.geo into InSAR_1D_Object.
    """
    xarray, yarray, data = isce_read_write.read_isce_unw_geo_alternative(isce_unw_filename)  # for brawley
    incidence = isce_read_write.read_scalar_data(los_filename, band=1)
    azimuth = isce_read_write.read_scalar_data(los_filename, band=2)
    lkv_e, lkv_n, lkv_u = insar_vector_functions.calc_lkv_from_rdr_azimuth_incidence(azimuth, incidence)

    [X, Y] = np.meshgrid(xarray, yarray)
    num_data = np.shape(data)[0] * np.shape(data)[1]
    lon = np.reshape(X, (num_data,))
    lat = np.reshape(Y, (num_data,))
    data = np.reshape(data, (num_data,))
    lkv_e = np.reshape(lkv_e, (num_data,))
    lkv_n = np.reshape(lkv_n, (num_data,))
    lkv_u = np.reshape(lkv_u, (num_data,))
    InSAR_data = Insar1dObject(lon=lon, lat=lat, LOS=data, LOS_unc=np.zeros(np.shape(data)), lkv_E=lkv_e,
                               lkv_N=lkv_n, lkv_U=lkv_u, starttime=starttime, endtime=endtime)
    return InSAR_data


def quick_convert_one_timeslice_to_disp(rateslice, date_intformat):
    """
    Compute displacement = rate * time

    :param rateslice: 2D array of velocities
    :param date_intformat: array of two strings, in [YYYYMMDD, YYYYMMDD] format.
    """
    numdata = np.shape(rateslice)[0] * np.shape(rateslice)[1]
    ratevector = np.reshape(rateslice, (numdata,))
    datevector = [dt.datetime.strptime(str(date_intformat[0]), "%Y%m%d"),
                  dt.datetime.strptime(str(date_intformat[1]), "%Y%m%d")]
    disps = convert_rates_to_disps(ratevector, datevector[0], datevector[1])
    return disps, datevector


def inputs_lohman_mcguire_2007(filename):
    """
    Read the Envisat data from the supplementary information of Lohman and McGuire 2007.
    Interferogram is from 21 August 2005 to 25 September 2005 (track 84, frame 2943).
    """
    xs, ys, los_list = np.loadtxt(filename, unpack=True)
    lon, lat, data = [], [], []
    for x, y, los in zip(xs, ys, los_list):
        latlontuple = io_other.convert_points_to_wgs84(x, y)
        lat.append(latlontuple[0])
        lon.append(latlontuple[1])
        data.append(10 * los)  # convert cm to mm.
    placeholder = np.zeros(np.shape(data))
    InSAR_data = Insar1dObject(lon=lon, lat=lat, LOS=data, LOS_unc=placeholder, lkv_E=placeholder, lkv_N=placeholder,
                               lkv_U=placeholder, starttime=dt.datetime.strptime('21-08-2005', '%d-%m-%Y'),
                               endtime=dt.datetime.strptime('25-09-2005', '%d-%m-%Y'))
    return InSAR_data
