#!/usr/bin/env python

"""
A model using the non-linear inversion that models slip as elliptical slip distributions.

# How to do this: create columns of fault elements.
# Create a function that takes Z and A and gives displacements
# A is non-linearly related to displacements.
Param vector is "SLIP -- DEPTH" for every fault patch.
"""

import numpy as np
from scipy.optimize import least_squares
import argparse
import json
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb.fault_slip_triangle import triangle_okada
from geodesy_modeling.datatypes.InSAR_1D_Object import covariance
from scipy.linalg import cholesky, solve_triangular
import elliptical_utilities  # local import
import inversion_utilities  # local import
import os


default_params = PyCoulomb.configure_calc.get_lightweight_config_params(mu=30e9, lame1=30e9, B=0)


def elastic_model(param_vector, cart_disp_points, faults, configs):
    """
    The forward model. The mesh geometry is known separately.

    :param cart_disp_points: used only for location of points, in cartesian coordinates
    :param param_vector: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :param configs: dictionary of configuration parameters for the experiment
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    fault_list = []
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        surface_slip, bottom_depth = param_vector[j], param_vector[j+configs["num_faults"]]  # for real faults. KEY LINE
        used_depths, widths, modeled_slip = elliptical_utilities.get_tuned_depth_arrays(configs['top_depth'],
                                                                                        configs['max_disloc_depth'],
                                                                                        configs["sampling_interval"],
                                                                                        surface_slip,
                                                                                        bottom_depth)
        for mslip, mwidth, depth in zip(modeled_slip, widths, used_depths):   # create all the source patches
            one_fault = fso.FaultSlipObject(strike=fault.strike, dip=89.99, length=fault.length,
                                            width=mwidth, lon=fault.lon,
                                            lat=fault.lat, depth=depth, rake=180, slip=mslip)
            # Change the slip to elliptical slip
            fault_list.append(one_fault)
    source_object = fso.fault_object_to_coulomb_fault(fault_list, zerolon_system=configs["zerolon"],
                                                      zerolat_system=configs["zerolat"])
    inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                    zerolon=configs['zerolon'],
                                                                                    zerolat=configs['zerolat'],
                                                                                    bbox=configs["bbox"])
    model_disp_points = triangle_okada.compute_cartesian_def_tris(inputs, default_params, cart_disp_points)  # run okada
    insar_1d_model = inversion_utilities.project_disp_points_into_insar1d(model_disp_points,
                                                                          configs["flight_angle"],
                                                                          configs["incidence_angle"])  # cart. space
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm


def set_up_initial_params_and_bounds(configs):
    # Set up constraints: lower bound on slip = from fieldwork and creepmeters
    _, fieldwork_lower_bounds = np.loadtxt(configs['fieldfile'], unpack=True)
    fieldwork_lower_bounds = np.multiply(fieldwork_lower_bounds, 0.001)  # convert mm to m
    lower_bound = np.zeros((configs["num_faults"]*2, ))   # lower bound on slip and depth is zero
    lower_bound[0:configs["num_faults"]] = fieldwork_lower_bounds
    upper_bound = np.multiply(5, np.ones((configs["num_faults"]*2, )))  # upper bound on slip is 50 mm, depth is 5 km
    upper_bound[0:configs["num_faults"]] = 0.050

    # Set up the original parameter vector and the upper bounds and lower bounds
    param0 = []  # param vector = [slip slip slip ..... slip depth depth depth.....]
    for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
        slip_guess = np.max([0.02, lower_bound[i]])  # the initial guess is either 20 mm or the minimum of range
        param0.append(slip_guess)
    for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
        param0.append(2.0)
    return param0, lower_bound, upper_bound


def invert_data(arguments):
    # We can vary the smoothing parameter and determine the most appropriate one through L-curve analysis
    # lam = 1  # Minimum norm Tikhonov smoothing regularization strength
    # gamma = 50  # Laplacian Regularization strength

    with open(arguments.config) as f:
        configs = json.load(f)
    if not isinstance(configs, dict):
        raise TypeError("Expected a JSON object at top level")

    data, cart_disp_points, faults = inversion_utilities.read_data_and_faults(configs, os.path.join(arguments.output,
                                                                                                    'used_faults.txt'))

    print("Experiment Setup: Arguments: ")
    print("Gamma Laplacian Smoothing Strength: ", arguments.laplacian)
    print("Lambda Tikhonov Smoothing Strength: ", arguments.tikhonov)
    print("Output Directory: ", arguments.output)
    os.makedirs(arguments.output, exist_ok=True)

    # Determine covariance matrix, and compute the inverse by triangular matrices in the Cholesky decomposition.
    L, sigma = covariance.read_covd_parameters(configs["cov_parameters"])
    cov = covariance.build_Cd(data, sigma, L)
    covariance.plot_full_covd(cov, os.path.join(arguments.output, 'covariance_matrix.png'))

    # Whitening the covariance matrix through Cholesky decomposition
    Lt_mat = cholesky(cov, lower=True, check_finite=False)
    covariance.plot_full_covd(Lt_mat, os.path.join(arguments.output, 'cholesky_Lt.png'))

    def make_data_whitener(Cd):
        # Cd must be symmetric-positive-definite. If near-singular, consider jitter or SVD.
        Ld = cholesky(Cd, lower=True, check_finite=False)  # Cd = Ld Ld^T

        # Wd @ x == solve(Ld, x) ⇒ Wd = Ld^{-1} without forming it
        def Wd_apply0(x):
            # handle 1D vector or 2D matrix (apply to columns)
            return solve_triangular(Ld, x, lower=True, check_finite=False)
        return Wd_apply0

    Wd_apply = make_data_whitener(cov)

    lam = arguments.tikhonov
    gamma = arguments.laplacian
    param0, lb, ub = set_up_initial_params_and_bounds(configs)  # create initial parameter vector

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults, configs)
        return insar_1d_model

    # Smoothing residuals (first differences) on first half of vector
    def smoothing_residuals(m2, gamma0):
        return gamma0 * np.diff(m2)

    def residuals(m, data0, gamma0, lam0):
        data_misfit = Wd_apply(forward_model(m).LOS - data0.LOS)  # normalize the misfit by the sqrt(cov_matrix)
        m2 = m[0:configs["num_faults"]]  # slip is the first half of the model vector
        mdepth = m[configs["num_faults"]:]
        tikhonov = lam0 * mdepth  # Tikhonov (minimum-norm) regularization
        smoothing = smoothing_residuals(m2, gamma0)
        return np.concatenate((data_misfit, smoothing, tikhonov))

    expname = 'laplacian_'+str(gamma)+'_tikhonov_'+str(lam)

    # NEXT: Put the additional three parameters for offset and planar fit.
    result = least_squares(residuals, x0=param0, verbose=True, bounds=[lb, ub], args=(data, gamma, lam))  # slip, z
    print(result.x)
    model_pred = forward_model(result.x)
    model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
    inversion_utilities.data_model_misfit_plot(data, model_pred, result.x, faults,
                                               os.path.join(arguments.output, expname+"_data_v_model.png"))
    inversion_utilities.write_outputs(data, model_pred, result.x, lam, gamma, arguments.output, expname, configs)

    # testcase_params = param0
    # simple_model = forward_model(testcase_params)  # InSAR1D object
    # simple_model = inversion_utilities.convert_xy_to_ll_insar1D(simple_model, configs)
    # inversion_utilities.data_model_misfit_plot(data, simple_model, testcase_params, faults,
    #                                            outname=os.path.join(arguments.output, "testcase_model.png"))
    # inversion_utilities.write_outputs(data, simple_model, testcase_params, lam, gamma, arguments.output,
    #                                   'testcase', configs)
    return


def do_main(parsed):
    invert_data(parsed)   # Build upon prior inversions, using InSAR1D object as data, invert for S0 and Z
    return


def parse_arguments():

    # 1. Create the parser
    parser = argparse.ArgumentParser(description="Invert data for elliptical slip distributions on a fault geometry.")

    # 2. Add arguments
    parser.add_argument(
        "-l", "--tikhonov",
        type=float,
        default=1,
        help="Lambda, Tikhonov smoothing (default: 1)"
    )
    parser.add_argument(
        "-g", "--laplacian",
        type=float,
        default=1,
        help="Gamma, Laplacian smoothing (default: 1)"
    )
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Path to the input config json file"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output directory"
    )

    # 3. Parse the command line
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    argms = parse_arguments()  # returns a dictionary? class?
    do_main(argms)
