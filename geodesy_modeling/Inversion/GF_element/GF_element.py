from elastic_stresses_py.PyCoulomb import disp_points_object as dpo
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
from Tectonic_Utils.geodesy import euler_pole


class GF_element:
    """
    GF_element is everything you would need to make a column of the Green's matrix and
    plot the impulse response function.
    'points' is coordinates of surface trace of the fault, if applicable.

    :param disp_points: modeled displacement_points due to unit activation of this GF_element
    :type disp_points: list of Displacement_points objects
    :param param_name: param_name
    :type param_name: string
    :param fault_dict_list: list of fault_slip_objects
    :type fault_dict_list: list
    :param upper_bound: highest allowed value of this GF_element
    :type upper_bound: float
    :param lower_bound: lowest allowed value of this GF_element
    :type lower_bound: float
    :param slip_penalty: a number that will be attached to minimum-norm smoothing in the G matrix
    :type slip_penalty: float
    :param units: what units is this 'unit activation' in?
    :type units: string
    :param points: coordinates of surface trace of fault, if provided
    :type points: np.array
    """

    def __init__(self, disp_points, param_name='', upper_bound=0, lower_bound=0, slip_penalty=0, units='',
                 fault_dict_list=(), points=()):
        self.disp_points = disp_points
        self.param_name = param_name
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound
        self.slip_penalty = slip_penalty
        self.units = units
        self.fault_dict_list = fault_dict_list
        self.points = points  # coordinates of surface trace of fault, if applicable

    def set_param_name(self, param_name):
        self.param_name = param_name

    def set_lower_bound(self, lower_bound):
        self.lower_bound = lower_bound

    def set_upper_bound(self, upper_bound):
        self.upper_bound = upper_bound

    def set_units(self, units_str):
        self.units = units_str

    def set_slip_penalty(self, slip_penalty):
        self.slip_penalty = slip_penalty


def get_GF_rotation_element(obs_disp_points, ep, target_region=(-180, 180, -90, 90), rot_name=''):
    """
    Build one GF_element for a horizontal rotation of GNSS velocities about some axis

    :param obs_disp_points: list of disp_points
    :param ep: [lon, lat, rate] of euler pole for desired rotation
    :param target_region: list of lon/lon/lat/lat for bounding box
    :param rot_name: string, optional metadata for naming the rotation (ex: ocb_)
    :returns: one GF_element with rotation displacements in x, y, and z directions
    """
    rot_disp_p = []
    for obs_item in obs_disp_points:
        coords = [obs_item.lon, obs_item.lat]
        if obs_item.is_within_bbox(target_region):
            mult = 1
        else:
            mult = 0
        response_to_rot = euler_pole.point_rotation_by_Euler_Pole(coords, ep)
        response = Displacement_points(lon=obs_item.lon, lat=obs_item.lat, dE_obs=mult*response_to_rot[0],
                                       dN_obs=mult*response_to_rot[1], dU_obs=mult*response_to_rot[2], Se_obs=0,
                                       Sn_obs=0, Su_obs=0, meas_type=obs_item.meas_type, refframe=obs_item.refframe,
                                       name=obs_item.name)
        rot_disp_p.append(response)
    rotation_gf = GF_element(disp_points=rot_disp_p, param_name=rot_name, upper_bound=1, lower_bound=-1,
                             slip_penalty=0, units='deg/Ma')
    return rotation_gf


def get_GF_rotation_elements(obs_disp_points, target_region=(-180, 180, -90, 90), rot_name=''):
    """
    Build 3 GF_elements for horizontal rotation of GNSS velocities due to reference frames
    X rotation: [0, 0, 1] Euler Pole
    Y rotation: [90, 0, 1] Euler Pole
    Z rotation: [0, 89.99, 1] Euler Pole

    :param obs_disp_points: list of disp_points
    :param target_region: list of lon/lon/lat/lat for bounding box
    :param rot_name: string, optional metadata for naming the rotation (ex: ocb_)
    :returns: list of 3 GF_elements with rotation displacements in x, y, and z directions
    """
    xresponse = get_GF_rotation_element(obs_disp_points, ep=[0, 0, 1], rot_name=rot_name + 'x_rot',
                                        target_region=target_region)  # X direction
    yresponse = get_GF_rotation_element(obs_disp_points, ep=[90, 0, 1], rot_name=rot_name + 'y_rot',
                                        target_region=target_region)  # Y direction
    zresponse = get_GF_rotation_element(obs_disp_points, ep=[0, 89.99, 1], rot_name=rot_name + 'z_rot',
                                        target_region=target_region)  # Z direction
    return [xresponse, yresponse, zresponse]


def get_GF_leveling_offset_element(obs_disp_points):
    """
    Build one GF_element for a reference frame leveling offset column of the GF matrix

    :param obs_disp_points: list of disp_point_objects
    :returns: a list of 1 GF_element, or an empty list if there is no leveling in this dataset
    """
    total_response_pts = []
    lev_count = 0
    for item in obs_disp_points:
        if item.meas_type == "leveling":
            response = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=1,
                                           Se_obs=0, Sn_obs=0, Su_obs=0, meas_type=item.meas_type)
            lev_count += 1
        else:
            response = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=0,
                                           Se_obs=0, Sn_obs=0, Su_obs=0, meas_type=item.meas_type)
        total_response_pts.append(response)
    lev_offset_gf = GF_element(disp_points=total_response_pts, param_name='lev_offset',
                               upper_bound=1, lower_bound=-1, slip_penalty=0, units='m/yr')
    if lev_count == 0:
        return []
    else:
        return [lev_offset_gf]


def add_gfs(GFs_list):
    """
    Take several gf_elements and add their green's functions together.

    :param GFs_list: list of GF_elements with the same modeled_disp_points lists
    :returns: list of disp_points
    """
    new_pts = GFs_list[0].disp_points
    for i in range(1, len(GFs_list)):
        new_pts = dpo.utilities.add_disp_points(new_pts, GFs_list[i].disp_points)
    return new_pts
