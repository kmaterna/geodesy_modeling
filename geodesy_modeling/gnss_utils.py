"""
Higher-level conversion between GNSS Station_Vels and ESPy Disp_Points
Note that one is meant for velocities and one is meant for displacements
Be smart when you're using both.
"""

from gnss_timeseries_viewers.gps_tools import vel_functions
from Elastic_stresses_py.PyCoulomb import coulomb_collections


def disp_points_to_Station_Vels(list_of_disp_points):
    list_of_station_vels = []
    for dpo in list_of_disp_points:
        new_station_vel = vel_functions.Station_Vel(elon=dpo.lon, nlat=dpo.lat, e=dpo.dE_obs*1000, n=dpo.dN_obs*1000,
                                                    u=dpo.dU_obs*1000, se=dpo.Se_obs*1000, sn=dpo.Sn_obs*1000,
                                                    su=dpo.Su_obs*1000, name=dpo.name, first_epoch=dpo.starttime,
                                                    last_epoch=dpo.endtime, refframe=dpo.refframe,
                                                    meas_type=dpo.meas_type)  # in mm/yr
        list_of_station_vels.append(new_station_vel)
    return list_of_station_vels


def Station_Vels_to_disp_points(list_of_station_vels):
    list_of_disp_points = []
    for station_vel in list_of_station_vels:
        new_dpo = coulomb_collections.Displacement_points(lon=station_vel.elon, lat=station_vel.nlat,
                                                          dE_obs=station_vel.e/1000, dN_obs=station_vel.n/1000,
                                                          dU_obs=station_vel.u/1000, Se_obs=station_vel.se/1000,
                                                          Sn_obs=station_vel.sn/1000, Su_obs=station_vel.su/1000,
                                                          name=station_vel.name, starttime=station_vel.first_epoch,
                                                          endtime=station_vel.last_epoch, refframe=station_vel.refframe,
                                                          meas_type=station_vel.meas_type)   # in m
        list_of_disp_points.append(new_dpo)
    return list_of_disp_points
