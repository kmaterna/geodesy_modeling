#!/usr/bin/env python

"""
Elements of a toolbox for 1-dimensional stress-driven inversion.
This toolbox can:
* Construct a subdivided one-dimensional fault
* Construct the shear stress Green's functions
* Plot the resulting stress, slip, or surface displacement across the fault
"""

import numpy as np
import tectonic_utils.geodesy.fault_vector_functions as fvf
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
import matplotlib.pyplot as plt
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points


# Setting up the fixed 2D parameters that allow the experiment to operate in 1D.  These never change.
paramdict = {"zerolon": -115.0,
             "zerolat": 33.0,
             "bbox": [-116, -114, 32, 34],
             "length": 100}


class DislocArray_1D:
    def __init__(self, fault_list, elastic_params):
        """
        An object to facilitate one-dimensional modeling experiments.  Fault patches with slip in a 1D case.

        :param fault_list: list of PyCoulomb faults in the standard way
        :param elastic_params: a PyCoulomb Parameters class that contains mu, lame1, and B, options, and File-IO
        """
        self.fault_list = fault_list
        self.elastic_params = elastic_params
        self.shear_stress_drop, self.normal_stress_drop = self.forward_model_stress_drop()

    def get_depth_array(self):
        return [x.top for x in self.fault_list]

    def get_slip_array(self):
        return [x.get_fault_slip() for x in self.fault_list]

    def plot_slip_with_depth(self, outfile_name):
        """Extract the slip vs. depth profile for this list of fault objects. Plot it. """
        depth_values, slip_values = self.get_depth_array(), self.get_slip_array()
        plot_slip_with_depth(depth_values, slip_values, outfile_name)
        return

    def plot_stress_with_depth(self, outfile_name):
        """Extract the slip vs. depth profile for this list of fault objects. Plot it. """
        depth_values = [x.top for x in self.fault_list]
        plot_stress_with_depth(depth_values, self.shear_stress_drop, outfile_name)
        return

    def forward_model_stress_drop(self):
        """Calculate the shear and normal stress profiles from the slip currently prescribed onto each fault."""
        inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=self.fault_list,
                                                                                        zerolon=paramdict["zerolon"],
                                                                                        zerolat=paramdict["zerolat"],
                                                                                        bbox=paramdict["bbox"])
        inputs = inputs.modify_inputs_object(receiver_object=self.fault_list)
        out_object = PyCoulomb.run_dc3d.do_stress_computation(self.elastic_params, inputs)
        return out_object.receiver_shear, out_object.receiver_normal

    def build_shear_stress_G(self):
        """Compute the Green's matrix for self-stress from each element to each other element. """
        depths = self.get_depth_array()

        i = 0
        G = np.zeros((len(depths), len(depths)))

        for receiver in self.fault_list:
            # Source for Green's functions.
            rtlat, reverse = fvf.get_rtlat_dip_slip(slip=1, rake=receiver.rake)  # unit slip = 1 m
            source_object = [receiver.modify_fault_object(rtlat=rtlat, reverse=reverse)]  # modify slip on patch to 1 m
            inpts = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                           zerolon=paramdict["zerolon"],
                                                                                           zerolat=paramdict["zerolat"],
                                                                                           bbox=paramdict["bbox"])
            inpts = inpts.modify_inputs_object(receiver_object=self.fault_list)
            out_object = PyCoulomb.run_dc3d.do_stress_computation(self.elastic_params, inpts)
            G_column = out_object.receiver_shear  # only shear?  Combination of shear and normal?  Coulomb?
            G[i, :] = G_column
            i = i + 1
        return G

    def forward_model_surface_disps(self):
        """
        Does not produce any file-IO. Just computes displacements across the 1d profile

        :return: list of displacement_points objects with modeled displacements
        """
        lats = np.multiply(33.15, np.ones((1, 100)))[0]  # Create points along a line half-way up the fault
        lons = np.linspace(-115.5, -114.5, 100)  # Create points horizontally across the fault trace
        obs_disp_points = [Displacement_points(lon=x, lat=y) for x, y in zip(lons, lats)]
        inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=self.fault_list,
                                                                                        zerolon=paramdict["zerolon"],
                                                                                        zerolat=paramdict["zerolat"],
                                                                                        bbox=paramdict["bbox"],
                                                                                        domainsize=93.09,
                                                                                        num_points_x=60,
                                                                                        num_points_y=60)
        out_object = PyCoulomb.run_dc3d.do_stress_computation(self.elastic_params, inputs, disp_points=obs_disp_points)
        return out_object.model_disp_points


def build_depths_array(bottom_depth, top_depth, number_segments):
    """
    Subdivides a depth range into depth intervals and builds a list of top depths of the sub-fault patches.
    Helper function for creating DislocDepths object

    :param bottom_depth: float, km, positive downward
    :param top_depth: float, km, positive downward
    :param number_segments: integer
    :return: list of depths, depth_interval
    """
    depth_interval = (bottom_depth - top_depth) / number_segments
    depths = np.linspace(top_depth, bottom_depth-depth_interval, number_segments)
    return depths, depth_interval


def construct_receiver_objects(configs):
    """
    Construct a list of one-dimensional fault patches by subdividing the fault into number_segments along dip.

    :param configs: a dictionary of config params with certain fields
    :return: list of pyCoulomb fault objects
    """
    depths, depth_interval = build_depths_array(configs["bottom_depth"], configs["top_depth"],
                                                configs["number_segments"])
    width = fvf.get_downdip_width(configs["top_depth"], configs["bottom_depth"], configs["dip"])
    configs["width"] = width
    configs["depth_interval"] = depth_interval
    one_fault = fso.FaultSlipObject(strike=0, dip=configs["dip"], length=paramdict["length"], width=configs["width"],
                                    lon=paramdict["zerolon"], lat=32.7, depth=configs["top_depth"],
                                    rake=configs["rake"], slip=0)
    [receiver_object] = fso.fault_object_to_coulomb_fault([one_fault], zerolon_system=paramdict["zerolon"],
                                                          zerolat_system=paramdict["zerolat"])
    receiver_object_list = receiver_object.split_single_fault(1, dip_num_split=configs["number_segments"])
    return receiver_object_list


def construct_receivers_Disloc_Profile(configs, elastic_params):
    """Almost like a driver for the config stage of this experiment. """
    receivers = construct_receiver_objects(configs)  # list of pyc Faults, same zerolon/zerolat, no slip
    disloc_profile = DislocArray_1D(receivers, elastic_params)  # create the 1D object that will do the work
    return disloc_profile


def construct_source_Disloc_Profile(receiver_profile, slip_array, elastic_params):
    """
    Build a Disloc_profile object with given slip

    :param receiver_profile: Disloc_profile object
    :param slip_array: 1d list of slip, in m
    :param elastic_params: a PyCoulomb Parameters class that contains mu, lame1, and B, among other things
    :return: Disloc_profile object
    """
    slipping_faults = []
    for fault, slip_amount in zip(receiver_profile.fault_list, slip_array):
        rtlat, reverse = fvf.get_rtlat_dip_slip(slip_amount, rake=fault.rake)
        source_fault = fault.modify_fault_object(rtlat=rtlat, reverse=reverse)  # modify slip to specific value (m)
        slipping_faults.append(source_fault)
    return DislocArray_1D(slipping_faults, elastic_params)


def plot_stress_with_depth(depths, pred_stress, outname, actual_stress=None):
    """
    Simple plot of stress with depths, taking just the one-dimensional arrays of stress and depth.

    :param depths: 1-d array
    :param pred_stress: 1-d array in kPa
    :param actual_stress: 1-d array in kPa
    :param outname: string
    """
    # Plot the model predicted stress-drop distribution, and the input stress-drop distribution
    plt.figure(figsize=(6, 9), dpi=300)
    plt.plot(pred_stress, depths, marker='.', label='Stress')
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
    Simple plot of slip with depths, taking just the one-dimensional arrays of slip and depth.
    Are the units reasonable?

    :param depths: 1d array
    :param slip: 1d array, in cm
    :param outname: string
    :param pred_slip: 1d array, in cm
    """
    # Plot the inverted slip distribution (in cm)
    plt.figure(figsize=(8, 9), dpi=300)
    plt.plot(slip, depths, marker='.', label='Slip')
    plt.xlabel("Slip (m)", fontsize=14)
    plt.ylabel("Depth (km)", fontsize=14)
    plt.gca().invert_yaxis()
    plt.plot([0, 0], [0, 3], '--k')
    if pred_slip is not None:
        plt.plot(pred_slip, depths, label='Original slip')
    plt.gca().legend(fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    print("Saving slip with depth as %s" % outname)
    plt.savefig(outname)
    return
