
from . import leveling_inputs


def get_onetime_displacements(LevList, start_index, end_index):
    """
    Restrict a leveling object to only a certain timerange, setting the start to zero displacement.
    start_index and end_index refer to slices of dtarray.
    The data will have a bunch of arrays of [0, enddisp].
    Sign convention of (end - start) displacements.
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
