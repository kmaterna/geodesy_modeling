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
from scipy.sparse import diags
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


def elastic_model(param_vector, data_1d, cart_disp_points, faults, configs):
    """
    The forward model. The mesh geometry is known separately.

    :param cart_disp_points: used only for location of points, in cartesian coordinates
    :param param_vector: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param data_1d: 1d_insar_object, which contains the look vector information for each pixel
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :param configs: dictionary of configuration parameters for the experiment
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    fault_list = []
    a, b, c = param_vector[-3], param_vector[-2], param_vector[-1]  # unpack parameters of the planar fit and refpix
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        surface_slip, bottom_depth = param_vector[j], param_vector[j+configs["num_faults"]]  # for real faults. KEY LINE
        surface_slip = surface_slip / 100  # converting from cm to m
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
    insar_1d_model = inversion_utilities.project_disp_points_into_insar1d(model_disp_points, data_1d)  # cart. space
    # Implement a planar fit to the data, which has three planar parameters for x and y and a constant offset
    insar_1d_model = insar_1d_model.subtract_ramp(a, b, c)
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm


def set_up_initial_params_and_bounds(configs, arguments):
    """
    Here we construct an initial parameter vector of length 2N + 3.
    The last three are plane and refpix parameters.
    The structure of the parameter vector is [N*slip_m, N*depth_km, a, b, refpix_offset]
    The structure of lower bounds is [2N+3]
    The structure of upper bounds is [2N+3]
    """
    # Set up constraints: lower bound on slip = from fieldwork and creepmeters
    numfaults = configs["num_faults"]
    lower_bound = np.zeros((numfaults*2 + 3, ))  # lower bound on slip and depth starts at zero
    upper_bound = np.zeros((numfaults*2 + 3,))

    _, fieldwork_data = np.loadtxt(configs['fieldfile'], unpack=True)
    fieldwork_data = np.multiply(fieldwork_data, 0.1)  # convert mm to cm

    # Setting proper bounds
    lower_bound[0:numfaults] = fieldwork_data  # lower bound on slip from field data
    upper_bound[0:numfaults] = 6  # upper bound on slip is 6 cm
    lower_bound[numfaults:2*numfaults] = 0  # lower bound on depth is zero
    upper_bound[numfaults:2*numfaults] = 5  # upper bound on depth is 5 km
    lower_bound[-3], lower_bound[-2], lower_bound[-1] = -50, -50, -5.0  # parameters for plane and reference pixel
    upper_bound[-3], upper_bound[-2], upper_bound[-1] = 50, 50, 5.0

    if configs['starting_point'] is not None:
        param0 = np.loadtxt(configs['starting_point'])
        print("Writing initial guess in %s" % os.path.join(arguments.output, 'initial_guess.txt'))
    else:
        param0 = []   # param vector = [slip slip slip ..... slip depth depth depth..... a, b, c]
        for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
            slip_guess = np.max([2.0, lower_bound[i]])  # the initial guess is either 20 mm or the minimum of range
            param0.append(slip_guess)
        for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
            param0.append(2.0)
        param0 = param0 + [0, 0, 0]  # plane, plane, and reference pixel

    # Set up the original parameter vector and the upper bounds and lower bounds
    x_scale = []  #
    for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
        x_scale.append(1.0)  # scale is mm
    for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
        x_scale.append(1.0)
    x_scale = x_scale + [5.0, 5.0, 1.0]
    value_and_bounds = np.vstack((param0, lower_bound, upper_bound)).T
    np.savetxt(os.path.join(arguments.output, 'initial_guess.txt'), value_and_bounds)
    return param0, lower_bound, upper_bound, x_scale


def invert_data(arguments):
    # We can vary the smoothing parameter and determine the most appropriate one through L-curve analysis

    with open(arguments.config) as f:
        configs = json.load(f)
    if not isinstance(configs, dict):
        raise TypeError("Expected a JSON object at top level")

    print("Experiment Setup: Arguments: ")
    print("Gamma Laplacian Smoothing Strength: ", arguments.laplacian)
    print("Lambda Tikhonov Smoothing Strength: ", arguments.tikhonov)
    print("Output Directory: ", arguments.output)
    os.makedirs(arguments.output, exist_ok=True)

    gamma = arguments.tikhonov  # tikhonov smoothing over depth range
    lam = arguments.laplacian  # laplacian smoothing minimizing differences in slip

    data, cart_disp_points, faults = inversion_utilities.read_data_and_faults(configs, os.path.join(arguments.output,
                                                                                                    'used_faults.txt'))
    # "data" is a insar1d object, so it contains all of its look vector information.

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size 81, 39 slips + 39 depths + plane + offset
        insar_1d_model = elastic_model(params, data, cart_disp_points, faults, configs)
        return insar_1d_model

    # Determine covariance matrix, and compute the inverse by triangular matrices in the Cholesky decomposition.
    L, sigma = covariance.read_covd_parameters(configs["cov_parameters"])
    cov = covariance.build_Cd(data, sigma, L)

    # Whitening the covariance matrix through Cholesky decomposition
    Lt_mat = cholesky(cov, lower=True, check_finite=False)

    def make_data_whitener(Cd):
        # Cd must be symmetric-positive-definite. If near-singular, consider jitter or SVD.
        Ld = cholesky(Cd, lower=True, check_finite=False)  # Cd = Ld Ld^T

        # Wd @ x == solve(Ld, x) ⇒ Wd = Ld^{-1} without forming it
        def Wd_apply0(x):
            # handle 1D vector or 2D matrix (apply to columns)
            return solve_triangular(Ld, x, lower=True, check_finite=False)
        return Wd_apply0

    Wd_apply = make_data_whitener(cov)

    def laplacian_v5(m):
        """
        Create two residual vectors, one for Laplacian smoothing of the stress drop and
        one for minimizing the magnitude of the stress drop.
        We apply the same Laplacian smoothing to the depth and slip parameters, and an additional Tikhonov to depth.

        :return: Two vectors, the laplacian smoothing penalty and the minimum-norm tikhonov penalty
        """
        m_slip = m[0:configs["num_faults"]]  # slip is first half of model vector
        m_depth = m[configs["num_faults"]:-3]  # depth is 2nd half of model vector, with last 3 reserved for plane

        # Doing the smoothing penalty on slip
        n = len(m_slip)
        main = -2 * np.ones(n)
        off = 1 * np.ones(n - 1)
        main[0] = main[-1] = 0  # Natural at boundaries
        Lapl = diags([off, main, off], offsets=[-1, 0, 1], format="csr")

        # The minimum norm penalty on depth
        A = np.eye(len(m_depth))

        return Lapl@m_slip, Lapl@m_depth, A@m_slip, A@m_depth

    param0, lb, ub, xscale = set_up_initial_params_and_bounds(configs, arguments)  # create initial parameter vector

    def residuals_double_L(m, data0, gamma0, lam0):   # if we're doing normal residuals
        data_misfit = Wd_apply(forward_model(m).LOS - data0.LOS)  # normalize the misfit by the sqrt(cov_matrix)
        l1, l2, l3, l4 = laplacian_v5(m)  # slip, depth, slip, depth
        beta = 2  # A factor for reasonable balance between the strength of slip smoothing and depth smoothing
        return np.concatenate((data_misfit,
                               np.multiply(lam0, l1),  # Laplacian on slip
                               np.multiply(lam0, l2),   # Laplacian on depth  # We can increase this.
                               np.multiply(gamma0*beta, l3),  # Tikhonov on slip
                               np.multiply(gamma0*beta, l4)))  # Tikhonov on depth

    expname = 'laplacian_'+str(lam)+'_tikhonov_'+str(gamma)

    residuals = residuals_double_L

    # If we're doing a forward model:
    if arguments.forward:
        print("Forward modeling from file %s." % arguments.forward)
        param_vector = np.loadtxt(arguments.forward)
        if np.size(param_vector) > 81:
            param_vector = param_vector[:, 3]  # use the column that represents the optimal solution
        model_pred = forward_model(param_vector)
        model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
        inversion_utilities.data_model_misfit_plot(data, model_pred, faults,
                                                   os.path.join(arguments.output, "test_data_v_model.png"),
                                                   region=(-115.88, -115.65, 32.90, 33.05), s=20)
        inversion_utilities.write_outputs(data, model_pred, param_vector, lam, gamma, arguments.output,
                                          "test_", configs)

        full_residuals = residuals(param_vector, data, gamma, lam)  # full residual vector, with applied coefficients

        total_misift, data_misfit = inversion_utilities.plot_complete_residual_vector_and_results(full_residuals,
                                                                                                  data,
                                                                                                  model_pred,
                                                                                                  param_vector,
                                                                                                  faults,
                                                                                                  arguments.output)

        return data_misfit

    else:
        # Visualize the covariance matrix
        covariance.plot_full_covd(cov, os.path.join(arguments.output, 'covariance_matrix.png'))
        covariance.plot_full_covd(Lt_mat, os.path.join(arguments.output, 'cholesky_Lt.png'))

        # The full least squares analysis
        result = least_squares(residuals, x0=param0, verbose=True, bounds=[lb, ub], args=(data, gamma, lam),
                               x_scale=xscale)  # slip, z, a, b, c

        print(result.x)
        model_pred = forward_model(result.x)
        model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
        inversion_utilities.data_model_misfit_plot(data, model_pred, faults,
                                                   os.path.join(arguments.output, expname+"_data_v_model.png"),
                                                   region=(-115.88, -115.65, 32.90, 33.05), s=20)
        inversion_utilities.write_outputs(data, model_pred, result.x, lam, gamma, arguments.output, expname, configs)
        residuals = data.LOS - model_pred.LOS
        data_misfit = np.sqrt(np.mean(residuals**2))

        full_residuals = residuals(result.x, data, gamma, lam)  # full residual vector, with applied coefficients

        _, _ = inversion_utilities.plot_complete_residual_vector_and_results(full_residuals,
                                                                             data,
                                                                             model_pred,
                                                                             result.x,
                                                                             faults,
                                                                             arguments.output)

        # testcase_params = param0
        # simple_model = forward_model(testcase_params)  # InSAR1D object
        # simple_model = inversion_utilities.convert_xy_to_ll_insar1D(simple_model, configs)
        # inversion_utilities.data_model_misfit_plot(data, simple_model, faults,
        #                                            outname=os.path.join(arguments.output, "testcase_model.png"))
        # inversion_utilities.write_outputs(data, simple_model, testcase_params, lam, gamma, arguments.output,
        #                                   'testcase', configs)
        return data_misfit


def do_main(parsed):
    misfit = invert_data(parsed)   # Build upon prior inversions, using InSAR1D object as data, invert for S0 and Z
    return misfit


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

    parser.add_argument(
        "-f", "--forward",
        default=None,
        help="Path to forward model file"
    )

    parser.add_argument(
        "-b", "--b",
        type=float,
        default=1,
        help="Scalar multiple between smoothing strengths. Set to -1 if you want to smooth the stress drop."
    )

    # 3. Parse the command line
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    argms = parse_arguments()  # returns a dictionary? class?
    print(argms)
    _misfit = do_main(argms)
