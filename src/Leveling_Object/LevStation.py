import numpy as np
from scipy import stats


class LevStation:
    """
    LevStation: One object for each station. Units of meters.  Should be aggregated into lists-of-objects.
    """
    def __init__(self, name, lat, lon, dtarray, leveling, reflon, reflat):
        self.name = name;  # string
        self.lon = lon;  # float
        self.lat = lat;  # float
        self.dtarray = dtarray;  # list of dt.datetime objects
        self.leveling = leveling;  # list of floats, units of meters.
        self.reflon = reflon;  # float
        self.reflat = reflat;  # float

    def detrend_by_slope(self, slope):
        """
        Remove a slope (in m/yr) from leveling displacement object.

        :param slope: float, meters/year
        """
        date_int_array = [x.toordinal() for x in self.dtarray];
        detrended_data = np.subtract(self.leveling, [self.leveling[0] + slope * (1/365.24) *
                                                     (x-date_int_array[0]) for x in date_int_array]);
        detrended_station = LevStation(name=self.name, lon=self.lon, lat=self.lat, dtarray=self.dtarray,
                                       leveling=detrended_data, reflon=self.reflon, reflat=self.reflat);
        return detrended_station;

    def find_slope(self, start_time, end_time):
        """
        Returns the slope of leveling time series observations between start_time and end_time in meters/year.

        :param start_time: datetime object
        :param end_time: datetime object
        """
        dtarray_limited, lev_obs_limited = [], [];
        for i in range(len(self.dtarray)):
            if start_time < self.dtarray[i] < end_time and ~np.isnan(self.leveling[i]):
                dtarray_limited.append(self.dtarray[i]);
                lev_obs_limited.append(self.leveling[i]);
        if len(dtarray_limited) > 1:
            date_int_array = [x.toordinal() for x in dtarray_limited];
            slope_0, intercept, _, _, _ = stats.linregress(date_int_array, lev_obs_limited);   # slope in mm per day
            slope_per_year = slope_0 * 365.24;
            return slope_per_year;
        else:
            return np.nan;  # error code

    def get_onetime_displacements(self, start_index, end_index):
        """
        Restrict a leveling object to only a certain timerange, setting the start to zero displacement.
        start_index and end_index refer to slices of dtarray.
        Data will be [0, enddisp]. Sign convention of (end - start) displacements.
        Should eventually add the same function for datetimes instead of indices.
        """
        dtarray, new_station_list = [], [];
        dtarray.append(self.dtarray[start_index]);
        dtarray.append(self.dtarray[end_index]);
        referenced_data = [0, self.leveling[end_index] - self.leveling[start_index]];
        subsampled_station = LevStation(name=self.name, lon=self.lon, lat=self.lat, dtarray=dtarray,
                                        leveling=referenced_data, reflon=self.reflon, reflat=self.reflat);
        return subsampled_station;


# ---------------- ON LISTS OF OBJECTS --------------- #

def find_trend_list(LevList, start_time, end_time):
    """
    Returns the slope of the leveling time series observations between start_time and end_time, in meters/year

    :param LevList : leveling object
    :param start_time : datetime object
    :param end_time: datetime object
    """
    return [station.find_slope(start_time, end_time) for station in LevList];


def detrend_leveling_object_list(LevList, slopes):
    """
    Remove a set of slopes (in m/yr) from leveling displacement objects

    :param LevList: a list of Leveling objects
    :param slopes: list of floats, matching the list of stations in LevList
    :returns: a list of Leveling objects
    """
    detrended_lev_list = [];
    for i, station in enumerate(LevList):
        detrended_station = station.detrend_by_slope(slopes[i]);
        detrended_lev_list.append(detrended_station);
    return detrended_lev_list;


def get_onetime_displacements_list(LevList, start_index, end_index):
    """
    Restrict a leveling object to only a certain timerange, setting the start to zero displacement.
    start_index and end_index refer to slices of dtarray. Data will be [0, enddisp].
    Sign convention of (end - start) displacements.
    """
    return [station.get_onetime_displacements(start_index, end_index) for station in LevList];
