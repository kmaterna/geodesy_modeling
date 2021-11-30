

"""
Coulomb collections Displacement_points:
Displacement_points = collections.namedtuple('Disp_Points', [
    'lon', 'lat',
    'dE_obs', 'dN_obs', 'dU_obs',
    'Se_obs', 'Sn_obs', 'Su_obs',
    'name', 'starttime', 'endtime', 'refframe', 'meas_type'], defaults=(None,) * 13);
# Disp_points are now lists of individual disp_point elements
# Usually displacements here are in meters


GPS station_vel:
Station_Vel = collections.namedtuple("Station_Vel", ['name', 'nlat', 'elon', 'n', 'e', 'u', 'sn', 'se', 'su',
                                                     'first_epoch',
                                                     'last_epoch', 'refframe', 'proccenter', 'subnetwork', 'survey',
                                                     'meas_type'], defaults=(None,) * 15 + ('gnss',));
# Station_Vel are in mm/yr, with -180<lon<180, used for velfields
"""

import matplotlib.path

def filter_to_continuous_only(obs_points):
    """Filter a list of disp_points into a list of continuous GNSS disp_points"""
    continuous_obs_points = [];
    for item in obs_points:
        if item.meas_type == 'continuous':
            continuous_obs_points.append(item);
    return continuous_obs_points;


def filter_to_continuous_only_second_table(obs_points, second_table_obs_points):
    """Keep points on first table if they are continuous in second table"""
    continuous_obs_points = [];
    for i, item in enumerate(second_table_obs_points):
        if item.meas_type == 'continuous':
            continuous_obs_points.append(obs_points[i]);
    return continuous_obs_points;


def filter_to_remove_creep(obs_points, fault_points, radius_km=5):
    """
    Filter a list of disp_points to remove those that are very close to a creeping fault
    :param obs_points: list of disp_point objects
    :param fault_points: list of 2-tuples, lon-lat vertices of fault trace in creeping section
    :param radius_km: float, radius around fault trace, in km.  Ex: 5 km on each side of the fault.
    """
    far_field_points = [];
    newpath = matplotlib.path.Path(vertices=fault_points);
    radius_deg = 2 * radius_km * (1/111.000);  # 2x because we're taking 'radius' km on each side
    for item in obs_points:
        if newpath.contains_point((item.lon, item.lat), radius=radius_deg):
            continue;
        else:
            far_field_points.append(item);
    print("Filtering creep: Returning %d of %d points " % (len(far_field_points), len(obs_points)) );
    return far_field_points;
