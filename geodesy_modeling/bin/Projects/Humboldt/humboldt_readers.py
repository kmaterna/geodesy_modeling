# Reading functions for the humboldt bay vertical working group
# Tide gages and leveling surveys
# Some GNSS velocities too

import numpy as np
import pandas as pd
import datetime as dt
from gnss_timeseries_viewers.gps_tools import vel_functions
import elastic_stresses_py.PyCoulomb.coulomb_collections as cc
import elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from tectonic_utils.geodesy import euler_pole


def read_leveling(filename):
    """
    Read Humboldt Bay vertical working group data into GNSS Station_Vel object
    """
    print("Reading file %s " % filename)
    leveling_list = []
    df = pd.read_excel(filename)
    lat, lon = df['Latitude'], df['Longitude']
    vel, unc = df['uplift_rate'], df['rate_uncert']
    for i in range(len(lat)):
        newobj = vel_functions.Station_Vel(elon=lon[i], nlat=lat[i], e=np.nan, n=np.nan, u=vel[i], se=np.nan, sn=np.nan,
                                           su=unc[i], refframe='custom', meas_type='leveling', name='',
                                           first_epoch=None, last_epoch=None)
        leveling_list.append(newobj)
    return leveling_list


def read_tide_gage(filename):
    """
    Read Humboldt Bay vertical working group data
    """
    print("Reading file %s " % filename)
    valid_list = ['CC', 'TRD01', 'TRD02', 'MRS_01', 'SO', 'NS', 'HS_01']
    tide_gage_list = []
    df = pd.read_excel(filename)
    site = df['Site']
    lat = df['Latitude']
    lon = df['Longitude']
    vel = df['uplift_r']
    unc = df['uplift_u']
    start = df['epoch_1_beg']
    end = df['epoch_1_end']
    for i in range(len(site)):
        if site[i] in valid_list:
            newobj = vel_functions.Station_Vel(elon=lon[i], nlat=lat[i], e=np.nan, n=np.nan, u=vel[i], se=np.nan,
                                               sn=np.nan, su=unc[i], first_epoch=start[i].to_pydatetime(),
                                               last_epoch=end[i].to_pydatetime(), meas_type='tide_gage', name='')
            tide_gage_list.append(newobj)
    return tide_gage_list


def read_all_data_table(infile):
    """
    Reading a table of velocities from different types: continuous, campaign, leveling, and tide gage.
    The opposite of the function immediately below.
    Reads into cc.Displacement_Points objects (units of meters)
    """
    list_of_point_velocities = []
    ifile = open(infile, 'r')
    for line in ifile:
        temp = line.split()
        if temp[0] == '#':
            continue
        lon, lat, ve, vn = float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])
        vu, se, sn, su = float(temp[4]), float(temp[5]), float(temp[6]), float(temp[7])
        if temp[8] == 'nan':
            starttime = np.nan
        else:
            starttime = dt.datetime.strptime(temp[8], "%Y-%m-%d")
        if temp[9] == 'nan':
            endtime = np.nan
        else:
            endtime = dt.datetime.strptime(temp[9], "%Y-%m-%d")
        if temp[-1] == "tide_gage" and np.isnan(su):  # if unknown, set tide gage data to a large nominal uncertaint
            su = 0.5   # 0.5 mm/yr
        newobj = cc.Displacement_points(lon=lon, lat=lat, dE_obs=ve/1000, dN_obs=vn/1000, dU_obs=vu/1000,
                                        Se_obs=se/1000, Sn_obs=sn/1000, Su_obs=su/1000,
                                        starttime=starttime, endtime=endtime, meas_type=temp[-1], refframe='NA',
                                        name='')
        list_of_point_velocities.append(newobj)
    ifile.close()
    print("Reading file %s... %d points " % (infile, len(list_of_point_velocities)))
    return list_of_point_velocities


def read_correction_data_table(infile, _latlonfile):
    """
    Read the table of corrected velocities that Fred produced.
    Used for the ridge and transform faults and the southern SAF.
    Items listed in mm on the file but disp_points format takes meters.
    """
    list_of_point_velocities = []
    [lon, lat, dE, dN, dU] = np.loadtxt(infile, unpack=True)
    for i in range(len(lon)):
        newobj = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=dE[i] / 1000, dN_obs=dN[i] / 1000,
                                        dU_obs=dU[i] / 1000, Se_obs=0, Sn_obs=0, Su_obs=0, starttime=None,
                                        endtime=None, meas_type=None, name='', refframe=None)
        list_of_point_velocities.append(newobj)
    print("Reading file %s... %d points " % (infile, len(list_of_point_velocities)))
    return list_of_point_velocities


def read_corrected_data_table(infile):
    """Read a version of the data table produced by Fred (the version Fred corrected)"""
    list_of_point_velocities = []
    [lon, lat, dE, dN, dU] = np.loadtxt(infile, unpack=True, skiprows=1, usecols=(0, 1, 2, 3, 4))
    for i in range(len(lon)):
        newobj = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=dE[i] / 1000, dN_obs=dN[i] / 1000,
                                        dU_obs=dU[i] / 1000, Se_obs=0, Sn_obs=0, Su_obs=0, endtime=None,
                                        meas_type=None, name='', refframe=None, starttime=None)
        list_of_point_velocities.append(newobj)
    print("Reading file %s... %d points " % (infile, len(list_of_point_velocities)))
    return list_of_point_velocities


def read_addcp_correction_table(dispfile, coordfile):
    """Read a velocity table produced by strainCp and addC+Cp """
    [lat, lon] = np.loadtxt(coordfile, unpack=True, skiprows=1, usecols=(0, 1))
    list_of_point_velocities = []
    ifile = open(dispfile, 'r')
    counter = 0
    for line in ifile:
        de = float(line[20:33])
        dn = float(line[33:46])
        du = float(line[46:59])
        newobj = cc.Displacement_points(lon=lon[counter], lat=lat[counter], dE_obs=de/100, dN_obs=dn/100,
                                        dU_obs=du/100, Se_obs=0, Sn_obs=0, Su_obs=0, endtime=None, starttime=None,
                                        meas_type=None, name='', refframe=None)
        list_of_point_velocities.append(newobj)
        counter += 1
    ifile.close()
    print("Reading file %s... %d points " % (dispfile, len(list_of_point_velocities)))
    return list_of_point_velocities


def read_ghost_transient_table(infile, _latlonfile):
    velocity_disp_points = []
    [lon, lat, dE, dN] = np.loadtxt(infile, unpack=True)  # in mm/yr
    for i in range(len(lon)):
        newobj = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=dE[i] / 1000, dN_obs=dN[i] / 1000,
                                        dU_obs=0, Se_obs=0, Sn_obs=0, Su_obs=0, endtime=None, starttime=None,
                                        meas_type=None, name='', refframe=None)
        velocity_disp_points.append(newobj)
    print("Reading file %s... %d points " % (infile, len(velocity_disp_points)))
    return velocity_disp_points


def get_euler_pole_from_xyz(epx, epy, epz):
    """
    Some euler poles are expressed as a vector designed to be attached to the surface of the sphere.
    See McCaffrey 2007 for some examples. It's not my favorite, but it must be dealt with.

    :param epx: deg/ma
    :param epy: deg/ma
    :param epz: deg/ma
    """
    horizontal_projection = np.sqrt(np.square(epx) + np.square(epy))
    lat = np.rad2deg(np.arctan(-epz/horizontal_projection))
    lon = np.rad2deg(np.arctan(epy/epx))-180
    omega = np.sqrt(np.square(epz) + np.square(epy) + np.square(epz))
    return lon, lat, omega


def get_euler_pole_correction(ep, latlonfile):
    """Specific Euler Pole correction for Oregon block"""
    correction_dps = []
    [lat, lon] = np.loadtxt(latlonfile, unpack=True, skiprows=1, usecols=(0, 1))
    ep_lon, ep_lat, ep_omega = get_euler_pole_from_xyz(ep[0], ep[1], ep[2])
    for i in range(len(lon)):
        if lat[i] < 41.5:
            newobj = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=0, dN_obs=0,
                                            dU_obs=0, Se_obs=0, Sn_obs=0, Su_obs=0, endtime=None, starttime=None,
                                            meas_type=None, name='', refframe=None)   # in meters/yr
        else:
            Point = [lon[i], lat[i]]
            Euler_Pole = [ep_lon, ep_lat, ep_omega]
            [e_vel, n_vel, u_vel] = euler_pole.point_rotation_by_Euler_Pole(Point, Euler_Pole)
            newobj = cc.Displacement_points(lon=lon[i], lat=lat[i], dE_obs=-e_vel / 1000, dN_obs=-n_vel / 1000,
                                            dU_obs=-u_vel/1000, Se_obs=0, Sn_obs=0, Su_obs=0, endtime=None,
                                            starttime=None, meas_type=None, name='', refframe=None)   # meters/yr
        correction_dps.append(newobj)
    fso.plot_fault_slip.map_source_slip_distribution([], "ep_correction.png", disp_points=correction_dps,
                                                     region=[-125, -120, 38, 43],
                                                     scale_arrow=(1.0, 0.010, "10 mm"),)
    return correction_dps


def write_all_data_table(list_of_station_vels, outfile):
    """
    Write a complicated, combined list of stations, velocities, types, reference frames, etc.
    Input is list of station_vels
    """
    print("Writing file %s " % outfile)
    ofile = open(outfile, 'w')
    ofile.write("# lon, lat, Ve, Vn, Vu, Se, Sn, Su, starttime, endtime, reframe, type\n")
    for item in list_of_station_vels:  # both gps and tide gage / leveling
        ofile.write("%f %f %f %f %f %f %f %f " % (item.elon, item.nlat, item.e, item.n, item.u, item.se, item.sn,
                                                  item.su))
        if item.meas_type == 'gnss':   # gps velocity
            ofile.write("%s %s %s " % (dt.datetime.strftime(item.first_epoch, "%Y-%m-%d"),
                                       dt.datetime.strftime(item.last_epoch, "%Y-%m-%d"), item.refframe))
            if item.survey:
                ofile.write("survey")
            else:
                ofile.write("continuous")
        elif item.meas_type == 'leveling':  # tide gage or leveling
            ofile.write("nan nan lev %s" % item.meas_type)
        else:  # tide gage
            ofile.write("%s %s %s %s" % (dt.datetime.strftime(item.first_epoch, "%Y-%m-%d"),
                                         dt.datetime.strftime(item.last_epoch, "%Y-%m-%d"), "None", item.meas_type))
        ofile.write("\n")
    ofile.close()
    return
