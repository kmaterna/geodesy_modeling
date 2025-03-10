#!/usr/bin/env python

"""
September 13, 2024.
Make a strike-slip fault. Prescribe its stress drop, and then solve for fault slip that releases the given stress drop.
Assumptions include long, strike-slip dislocations. The stress is resolved at the center of each element.
Only shear stress is considered.  Notes:
* Make sure the bounds on your inversion are wide enough (including negative numbers) to constrain the
slip that gives the right amount of stress drop.
* Make sure the code is able to compute self-stress, i.e., the stress drop on exactly the location of the
source patch.

Next: Convert the slip distribution with depth into a forward-model-able distribution of fault patches.
Could be through converting the whole tool into GF_elements.
"""

from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
import numpy as np
import scipy
import matplotlib.pyplot as plt
import Tectonic_Utils.geodesy.fault_vector_functions as fvf


start_configs = {"top_depth": 0,
                 "bottom_depth": 3,
                 "number_segments": 30,
                 "dip": 89.99,
                 "zerolon": -115.0,
                 "zerolat": 33}


def build_configs():
    depths = np.linspace(start_configs["top_depth"], start_configs["bottom_depth"], start_configs['number_segments'])
    config_params = start_configs
    depth_interval = (config_params["bottom_depth"] - config_params["top_depth"]) / config_params["number_segments"]
    width = fvf.get_downdip_width(config_params["top_depth"], config_params["bottom_depth"], config_params["dip"])
    config_params["width"] = width
    config_params["depth_interval"] = depth_interval
    return config_params, depths


def construct_receiver_objects(configs):
    """
    :param configs: a dictionary
    :return: list of pyCoulomb fault objects
    """
    one_fault = fso.FaultSlipObject(strike=0, dip=configs["dip"], length=100, width=configs["width"],
                                    lon=-115.0, lat=32.8, depth=configs["top_depth"], rake=180, slip=0)
    [receiver_object] = fso.fault_object_to_coulomb_fault([one_fault], configs["zerolon"], configs["zerolat"])
    receiver_object_list = receiver_object.split_single_fault(1, dip_num_split=configs["number_segments"])
    return receiver_object_list


def surface_rupture_crack_model(depths_1d, tau, crack_depth, mu):
    """
    Create the slip distribution of a surface rupturing mode III crack from a theoretical prediction

    :return:
    """
    depths_1d = np.multiply(depths_1d, 1000)  # convert km to m
    crack_depth = np.multiply(crack_depth, 1000)  # convert km to m
    theory_slip = 2 * tau / mu * np.sqrt(np.subtract(crack_depth**2, np.square(depths_1d)))  # Segall Page 96-ish
    return theory_slip


def build_G(configs, depths):
    default_params = PyCoulomb.configure_calc.configure_stress_calculation("my_config.txt")

    i = 0
    G = np.zeros((len(depths), len(depths)))

    receiver_object_list = construct_receiver_objects(configs)  # list of pyc Faults, same zerolon/zerolat

    for receiver in receiver_object_list:
        # Source for Green's functions.
        source_object = [receiver.modify_fault_object(rtlat=0.01)]  # modify slip on this patch to 1 cm
        inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                        zerolon=configs['zerolon'],
                                                                                        zerolat=configs['zerolat'],
                                                                                        bbox=[-116, -114, 32, 34])
        inputs = inputs.modify_inputs_object(receiver_object=receiver_object_list)
        out_object = PyCoulomb.run_dc3d.do_stress_computation(default_params, inputs)
        # PyCoulomb.output_manager.produce_outputs(default_params, inputs, (), (), out_object)
        G_column = out_object.receiver_shear
        G[i, :] = G_column
        i = i + 1

    return G, receiver_object_list


def plot_slip(depths, slip, G, actual_stress, config_params):

    # Forward model:
    pred_stress = np.dot(G, slip)
    print("Predicted stress:", pred_stress)  # The elliptical slip distribution should come out. And it does!

    # Plot the inverted slip distribution (in cm)
    plt.figure(figsize=(8, 9), dpi=300)
    plt.plot(slip, depths, marker='.', label='Inverted slip')
    plt.xlabel("Slip (cm)", fontsize=14)
    plt.ylabel("Depth (km)", fontsize=14)
    theory_slip = np.multiply(100, surface_rupture_crack_model(depths, mu=30e9, tau=-10e6,
                                                               crack_depth=config_params["bottom_depth"]))
    plt.plot(theory_slip, depths, color='green', label='Theory predicted')
    plt.gca().invert_yaxis()
    plt.gca().legend(fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.savefig("_inverted_slip.png")

    # Plot the model predicted stress-drop distribution, and the input stress-drop distribution
    plt.figure(figsize=(6, 9), dpi=300)
    pred_stress = np.divide(pred_stress, 1000)
    actual_stress = np.divide(actual_stress, 1000)
    plt.plot(pred_stress, depths, marker='.', label='Inverted Stress')
    plt.plot(actual_stress, depths, label='Applied Stress')
    plt.legend(fontsize=14)
    plt.xlabel("Stress (kPa)", fontsize=14)
    plt.ylabel("Depth (km)", fontsize=14)
    plt.gca().invert_yaxis()
    plt.gca().tick_params(axis='both', which='major', labelsize=14)
    plt.savefig("_stress_results.png")
    return


def main():
    config_params, depths = build_configs()
    G, receivers = build_G(config_params, depths)
    delta_tau = 10e3 * np.ones((len(depths),))  # stress drop constant within the crack (10 MPa)
    # delta_tau = np.add(delta_tau, np.multiply(depths, 1000))
    # delta_tau[-15:] = 0  # optionally: zero stress down at the bottom
    lb = -15000 * np.ones((len(depths),))  # make sure these amounts are big enough and the ranges are wide enough
    ub = 0 * np.ones((len(depths),))
    # lb[-15:] = -0.01  # optionally: pin the bottom.  This doesn't work very well.
    response = scipy.optimize.lsq_linear(G, delta_tau, bounds=(lb, ub), max_iter=1500, method='bvls')
    slip = response.x  # parameters of best-fitting model
    plot_slip(depths, slip, G, delta_tau, config_params)
    return


if __name__ == "__main__":
    main()
