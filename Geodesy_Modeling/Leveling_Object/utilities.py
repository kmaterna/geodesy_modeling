
from . import leveling_inputs


def get_onetime_displacements(data, start_index, end_index):
    """
    Restrict a leveling object to only a certain timerange, setting the start to zero displacement.
    start_index and end_index refer to slices of dtarray.
    The data will have a bunch of arrays of [0, enddisp].
    Sign convention of (end - start) displacements.
    """
    dtarray, leveling = [], [];
    dtarray.append(data.dtarray[start_index]);
    dtarray.append(data.dtarray[end_index]);
    for i in range(len(data.name)):
        referenced_data = [];
        referenced_data.append(0);
        referenced_data.append(data.leveling[i][end_index] - data.leveling[i][start_index])
        leveling.append(referenced_data);
    displacement_lev_data = leveling_inputs.LevData(name=data.name, lon=data.lon, lat=data.lat, dtarray=dtarray,
                                leveling=leveling);
    return displacement_lev_data;

