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
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.inputs_object import io_intxt


class DislocArray_1D:
    def __init__(self, fault_list):
        """
        An object to facilitate one-dimensional modeling experiments.  Fault patches with slip in a 1D case.

        :param fault_list: list of PyCoulomb faults I guess
        """
        self.fault_list = fault_list
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
        """Calculate the stress drop from G."""
        # must configure calculation somehow
        default_params = PyCoulomb.configure_calc.configure_stress_calculation("my_config.txt")  # FIX THIS
        inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=self.fault_list,
                                                                                        zerolon=-115.0,
                                                                                        zerolat=33,
                                                                                        bbox=[-116, -114, 32, 34])
        inputs = inputs.modify_inputs_object(receiver_object=self.fault_list)
        out_object = PyCoulomb.run_dc3d.do_stress_computation(default_params, inputs)
        return out_object.receiver_shear, out_object.receiver_normal

    def build_shear_stress_G(self):
        """Compute the Green's matrix for self-stress from each element to each other element. """
        default_params = PyCoulomb.configure_calc.configure_stress_calculation("my_config.txt")  # FIX THIS
        depths = self.get_depth_array()

        i = 0
        G = np.zeros((len(depths), len(depths)))

        for receiver in self.fault_list:
            # Source for Green's functions.
            rtlat, reverse = fvf.get_rtlat_dip_slip(1, rake=receiver.rake)  # unit slip = 1 m
            source_object = [receiver.modify_fault_object(rtlat=rtlat, reverse=reverse)]  # modify slip on patch to 1 cm
            inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                            zerolon=-115.0,
                                                                                            zerolat=33,
                                                                                            bbox=[-116, -114, 32, 34])
            inputs = inputs.modify_inputs_object(receiver_object=self.fault_list)
            out_object = PyCoulomb.run_dc3d.do_stress_computation(default_params, inputs)
            G_column = out_object.receiver_shear  # only shear?  Combination of shear and normal?  Coulomb?
            G[i, :] = G_column
            i = i + 1
        return G

    def forward_model_surface_disps(self):
        default_params = PyCoulomb.configure_calc.configure_stress_calculation("my_config.txt")  # FIX THIS
        # Create points along a line
        lats = np.multiply(33.25, np.ones((1, 100)))
        lons = np.linspace(-115.5, -114.5, 100)
        save_data = np.vstack((lons, lats)).T
        np.savetxt("profile_coords.txt", save_data, "%f %f", header="Lon Lat")

        # Re-package the slip into the fault patches that they came from
        slipping_faults = []
        for fault in self.fault_list:
            rtlat, reverse = fvf.get_rtlat_dip_slip(fault.get_fault_slip(), rake=fault.rake)
            source_fault = fault.modify_fault_object(rtlat=rtlat, reverse=reverse)  # modify slip inverted value (m)
            slipping_faults.append(source_fault)
        inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=slipping_faults,
                                                                                        zerolon=-115.0,
                                                                                        zerolat=33,
                                                                                        bbox=[-116, -114, 32, 34])
        io_intxt.write_intxt(inputs, "resulting_patches.txt")

        # Do a displacement calculation
        default_params = default_params.modify_params_object(disp_points_file='profile_coords.txt', plot_grd_disp=True,
                                                             input_file='resulting_patches.txt')
        [inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(default_params)
        out_object = PyCoulomb.run_dc3d.do_stress_computation(default_params, inputs, disp_points=obs_disp_points)
        PyCoulomb.output_manager.produce_outputs(default_params, inputs, (), (), out_object)
        return


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
    one_fault = fso.FaultSlipObject(strike=0, dip=configs["dip"], length=100, width=configs["width"],
                                    lon=-115.0, lat=32.8, depth=configs["top_depth"], rake=configs["rake"], slip=0)
    [receiver_object] = fso.fault_object_to_coulomb_fault([one_fault], zerolon_system=-115.0,
                                                          zerolat_system=33.0)
    receiver_object_list = receiver_object.split_single_fault(1, dip_num_split=configs["number_segments"])
    return receiver_object_list


def construct_receivers_Disloc_Profile(configs):
    """Almost like a driver for the config stage of this experiment. """
    receivers = construct_receiver_objects(configs)  # list of pyc Faults, same zerolon/zerolat, no slip
    disloc_profile = DislocArray_1D(receivers)  # create the 1D object that will do the work
    return disloc_profile


def construct_source_Disloc_Profile(receiver_profile, slip_array):
    """
    Build a Disloc_profile object with given slip

    :param receiver_profile: Disloc_profile object
    :param slip_array: 1d list of slip, in m
    :return: Disloc_profile object
    """
    slipping_faults = []
    for fault, slip_amount in zip(receiver_profile.fault_list, slip_array):
        rtlat, reverse = fvf.get_rtlat_dip_slip(slip_amount, rake=fault.rake)
        source_fault = fault.modify_fault_object(rtlat=rtlat, reverse=reverse)  # modify slip inverted value (m)
        slipping_faults.append(source_fault)
    return DislocArray_1D(slipping_faults)


def plot_stress_with_depth(depths, pred_stress, outname, actual_stress=None):
    """
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


def mock_1d_slip_structure(disloc_array):
    """
    Please take a distribution of slip with depth (one dimensional), and create a surface slip profile using
    elastic dislocation modeling from Okada (1985).
    The fault can be dipping. The rake does not have to be strike-slip.
    Inputs: A list of fault objects. We will assume they strike north-south, and the profile runs east-west
    Outputs: A list of cross-fault positions,
    Goals: Make this run in code. No file-IO should be used here.

    :return:
    """

    return
