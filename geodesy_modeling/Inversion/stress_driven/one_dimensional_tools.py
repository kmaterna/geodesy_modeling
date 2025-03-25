#!/usr/bin/env python

"""
Elements of a toolbox for 1-dimensional stress-driven inversion.
This toolbox can:
* Construct a subdivided one-dimensional fault
* Construct the shear stress Green's functions
* Plot the resulting stress, slip, or surface displacement across the fault
"""

import numpy as np
import Tectonic_Utils.geodesy.fault_vector_functions as fvf
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
import matplotlib.pyplot as plt


def build_depths(bottom_depth, top_depth, number_segments):
    """
    Subdivides a depth range into depth intervals and builds a list of top depths of the sub-fault patches.

    :param bottom_depth: float, km, positive downward
    :param top_depth: float, km, positive downward
    :param number_segments: integer
    :return: list of depths, depth_interval
    """
    depth_interval = (bottom_depth - top_depth) / number_segments
    depths = np.linspace(top_depth, bottom_depth-depth_interval, number_segments)
    return depths, depth_interval


def build_configs(config_params):
    """
    :param config_params: dictionary that contains "bottom_depth", "top_depth", and "number_segments", "dip"
    :return: list of depths, depth_interval
    """
    depths, depth_interval = build_depths(config_params["bottom_depth"], config_params["top_depth"],
                                          config_params["number_segments"])
    width = fvf.get_downdip_width(config_params["top_depth"], config_params["bottom_depth"], config_params["dip"])
    config_params["width"] = width
    config_params["depth_interval"] = depth_interval
    return config_params, depths


def construct_receiver_objects(configs):
    """
    Construct a list of one-dimensional fault patches by subdividing the fault into number_segments along dip.

    :param configs: a dictionary of config params
    :return: list of pyCoulomb fault objects
    """
    one_fault = fso.FaultSlipObject(strike=0, dip=configs["dip"], length=100, width=configs["width"],
                                    lon=-115.0, lat=32.8, depth=configs["top_depth"], rake=configs["rake"], slip=0)
    [receiver_object] = fso.fault_object_to_coulomb_fault([one_fault], configs["zerolon"], configs["zerolat"])
    receiver_object_list = receiver_object.split_single_fault(1, dip_num_split=configs["number_segments"])
    return receiver_object_list


def plot_stress_with_depth(depths, pred_stress, outname, actual_stress=None):
    """
    :param depths: 1-d array
    :param pred_stress: 1-d array in kPa
    :param actual_stress: 1-d array in kPa
    :param outname: string
    """
    # Plot the model predicted stress-drop distribution, and the input stress-drop distribution
    plt.figure(figsize=(6, 9), dpi=300)
    plt.plot(pred_stress, depths, marker='.', label='Inverted Stress')
    if actual_stress is not None:
        plt.plot(actual_stress, depths, label='Applied Stress')
    plt.legend(fontsize=14)
    plt.plot([0, 0], [0, 3], '--k')
    plt.xlabel("Stress (kPa)", fontsize=14)
    plt.ylabel("Depth (km)", fontsize=14)
    plt.gca().invert_yaxis()
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.savefig(outname)
    return


def plot_slip_with_depth(depths, slip, outname, pred_slip=None):
    """
    :param depths: 1d array
    :param slip: 1d array, in cm
    :param outname: string
    :param pred_slip: 1d array, in cm
    """
    # Plot the inverted slip distribution (in cm)
    plt.figure(figsize=(8, 9), dpi=300)
    plt.plot(slip, depths, marker='.', label='Inverted slip')
    plt.xlabel("Slip (cm)", fontsize=14)
    plt.ylabel("Depth (km)", fontsize=14)
    plt.gca().invert_yaxis()
    plt.plot([0, 0], [0, 3], '--k')
    if pred_slip is not None:
        plt.plot(pred_slip, depths, label='Original slip')
    plt.gca().legend(fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.savefig(outname)
    return
