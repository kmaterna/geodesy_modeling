

import numpy as np
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso


def read_elliptical_params(filename, num_faults):
    """
    Read the parameter file produced by the output of the inversion.
    N slip values, N depth values, for each ellipse
    Slip is in cm
    Depth is in km, positive down
    """
    params = np.loadtxt(filename)
    slips = params[0:num_faults]
    depths = params[num_faults:2*num_faults]
    return slips, depths


def get_elliptical_slip_amount(surface_slip, crack_depth, target_depth):
    """ Implement the elliptical slip equation. crack_depth and target_depth in km, positive downward. """
    if target_depth > crack_depth:
        return 0
    else:
        return np.sqrt(surface_slip**2 * (1 - (np.square(target_depth) / (crack_depth ** 2))))  # Depths are km.


def create_column_of_faults(orig_fso, surface_slip, ellipse_depth, depth_increment):
    """
    :param orig_fso: fault_slip_object
    :param surface_slip: slip in cm, float
    :param ellipse_depth: depth in km, positive downward, float
    :param depth_increment: km, float
    """
    new_fault_list = []
    depth_divisions = np.arange(orig_fso.depth, ellipse_depth+depth_increment, depth_increment)
    for depth in depth_divisions:
        slip = get_elliptical_slip_amount(surface_slip, ellipse_depth, depth)
        new_fault = fso.FaultSlipObject(strike=orig_fso.strike,
                                        dip=orig_fso.dip,
                                        length=orig_fso.length,
                                        lon=orig_fso.lon,
                                        lat=orig_fso.lat,
                                        width=depth_increment,
                                        depth=depth,
                                        rake=orig_fso.rake,
                                        slip=slip/100,
                                        tensile=orig_fso.tensile,
                                        segment=orig_fso.segment)
        new_fault_list.append(new_fault)
    return new_fault_list
