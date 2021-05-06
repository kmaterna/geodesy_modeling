
import numpy as np
from scipy import stats
from . import leveling_inputs


def get_onetime_displacements(LevList, start_index, end_index):
    """
    Restrict a leveling object to only a certain timerange, setting the start to zero displacement.
    start_index and end_index refer to slices of dtarray.
    The data will have a bunch of arrays of [0, enddisp].
    Sign convention of (end - start) displacements.
    Should eventually re-write this for datetimes instead of indices.
    """
    dtarray, new_station_list = [], [];
    dtarray.append(LevList[0].dtarray[start_index]);
    dtarray.append(LevList[0].dtarray[end_index]);  # assuming they all have the same list of dates
    for station in LevList:
        referenced_data = [0, station.leveling[end_index] - station.leveling[start_index]];
        subsampled_station = leveling_inputs.LevStation(name=station.name, lon=station.lon, lat=station.lat,
                                                        dtarray=dtarray, leveling=referenced_data,
                                                        reflon=station.reflon, reflat=station.reflat);
        new_station_list.append(subsampled_station);
    return new_station_list;


def find_trend(LevList, start_time, end_time):
    """
    Returns the slope of the leveling time series observations between start_time and end_time
    Units: meters/yr
    LevList : leveling object
    start_time : datetime object
    end_time: datetime object
    """
    slope_solutions = [];
    for item in LevList:
        dtarray_limited = [];
        lev_obs_limited = [];
        for i in range(len(item.dtarray)):
            if start_time < item.dtarray[i] < end_time and ~np.isnan(item.leveling[i]):
                dtarray_limited.append(item.dtarray[i]);
                lev_obs_limited.append(item.leveling[i]);
        if len(dtarray_limited) > 1:
            date_int_array = [x.toordinal() for x in dtarray_limited];
            slope_0, intercept, _, _, _ = stats.linregress(date_int_array, lev_obs_limited);   # slope in mm per day
            slope_per_year = slope_0 * 365.24;
            slope_solutions.append(slope_per_year);
        else:
            slope_solutions.append(np.nan);
    return slope_solutions;


def detrend_leveling_object(LevList, slopes):
    """
    Remove a set of slopes (in m/yr) from leveling displacement objects
    """
    detrended_lev_list = [];
    for i, station in enumerate(LevList):
        date_int_array = [x.toordinal() for x in station.dtarray];
        detrended_data = np.subtract(station.leveling, [station.leveling[0] + slopes[i] * (1/365.24) *
                                                        (x-date_int_array[0]) for x in date_int_array]);
        detrended_station = leveling_inputs.LevStation(name=station.name, lon=station.lon, lat=station.lat,
                                                       dtarray=station.dtarray, leveling=detrended_data,
                                                       reflon=station.reflon, reflat=station.reflat);
        detrended_lev_list.append(detrended_station);
    return detrended_lev_list;
