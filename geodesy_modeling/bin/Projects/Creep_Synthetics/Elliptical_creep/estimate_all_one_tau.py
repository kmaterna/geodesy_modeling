#!/usr/bin/env python

"""
A model using the non-linear inversion that models slip as elliptical slip distributions.
In this case, we are going to estimate slip at the surface and one constant stress drop, which gives us depth.
Param vector is "SLIP" for every fault patch, Tau, and then 3 planar fit parameters
"""

import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import diags
import argparse
import json
from geodesy_modeling.datatypes.InSAR_1D_Object import covariance
from scipy.linalg import cholesky, solve_triangular
import inversion_utilities  # local import
import experiment_specifics
import os


def const_tau_params_to_full_params(tau_params, n):
    """
    Convert n+4-parameter-vector into a 2n+3-parameter vector.
    Strain drop parameter is strain times 10^-5

    :param tau_params: vector of n slip, 1 tau, and 3 planar fit parameters
    :param n: number of faults, integer
    """
    slips = tau_params[0:n]
    strain_drop = tau_params[n]
    a, b, c = tau_params[n+1], tau_params[n+2], tau_params[n+3]
    depths = np.divide(slips, (2*strain_drop))
    full_params = np.concatenate((slips, depths))
    full_params = np.append(full_params, a)
    full_params = np.append(full_params, b)
    full_params = np.append(full_params, c)
    return full_params


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

    # Read InSAR data and fault geometry. "data" is an insar_1d object, so contains all of its look vector information.
    data, cart_disp_points, faults = inversion_utilities.read_data_and_faults(configs, os.path.join(arguments.output,
                                                                                                    'used_faults.txt'))

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size n + 4,
        # will convert to (2n+3): n slips + n depths + plane + offset
        full_params = const_tau_params_to_full_params(params, configs["num_faults"])  # convert 43 to 81
        print(full_params)
        insar_1d_model = experiment_specifics.elastic_model(full_params, data, cart_disp_points, faults, configs)
        return insar_1d_model

    # Determine covariance matrix, and compute the inverse by triangular matrices in the Cholesky decomposition.
    L, sigma = covariance.read_covd_parameters(configs["cov_parameters"])
    cov = covariance.build_Cd(data, sigma, L)
    cov = experiment_specifics.modify_cov_to_split_across_fault(cov, data, faults, jitter=30, lobotomy_value=1)

    # Whitening the data through its covariance matrix through Cholesky decomposition
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

        # Doing the smoothing penalty on slip -- second derivative smoothing
        n = len(m_slip)
        main = -2 * np.ones(n)
        off = 1 * np.ones(n - 1)
        Lapl = diags([off, main, off], offsets=[-1, 0, 1], format="csr")
        Lapl[0, 0:3] = 1, -2, 1  # Natural at boundaries
        Lapl[-1, -4:-1] = 1, -2, 1  # Natural at boundaries

        # The minimum norm penalty
        A = np.eye(len(m_slip))

        return Lapl@m_slip, A@m_slip

    # create initial parameter vector
    param0, lb, ub, xscale = experiment_specifics.set_up_initial_params_and_bounds_const_stress(configs, arguments)

    def residuals_double_L(m, data0, gamma0, lam0):   # if we're doing normal residuals
        data_misfit = Wd_apply(forward_model(m).LOS - data0.LOS)  # normalize the misfit by the sqrt(cov_matrix)
        l1, l2 = laplacian_v5(m)  # slip, depth, slip, depth
        beta = 1  # A factor for reasonable balance between the strength of slip smoothing and depth smoothing
        return np.concatenate((data_misfit,
                               np.multiply(lam0, l1),  # Laplacian on slip
                               np.multiply(gamma0*beta, l2)))  # Tikhonov on slip

    expname = 'laplacian_'+str(lam)+'_tikhonov_'+str(gamma)

    residuals = residuals_double_L

    # If we're doing a forward model:
    if arguments.forward:
        print("Forward modeling from file %s." % arguments.forward)
        # Visualize the covariance matrix
        covariance.plot_full_covd(cov, os.path.join(arguments.output, 'covariance_matrix.png'))
        covariance.plot_full_covd(Lt_mat, os.path.join(arguments.output, 'cholesky_Lt.png'))
        param_vector = np.loadtxt(arguments.forward)
        model_pred = forward_model(param_vector)
        model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
        inversion_utilities.data_model_misfit_plot(data, model_pred, faults,
                                                   os.path.join(arguments.output, "test_data_v_model.png"),
                                                   region=(-115.88, -115.65, 32.90, 33.05), s=20)

        np.savetxt(os.path.join(arguments.output, 'fitted_parameters.txt'), param_vector,
                   header="Params slip(cm), strain(1e-5), plane, plane, reference")  # original vector
        full_param_vector = const_tau_params_to_full_params(param_vector, configs["num_faults"])  # convert 43 to 81
        inversion_utilities.write_outputs(data, model_pred, full_param_vector, lam, gamma, arguments.output,
                                          "test_", configs)

        full_residuals = residuals(param_vector, data, gamma, lam)  # full residual vector, with applied coefficients

        total_misift, d_misfit = inversion_utilities.plot_complete_residual_vector_and_results(full_residuals,
                                                                                               data,
                                                                                               model_pred,
                                                                                               full_param_vector,
                                                                                               faults,
                                                                                               arguments.output,
                                                                                               arguments.laplacian)
        return d_misfit

    else:
        # Visualize the covariance matrix
        covariance.plot_full_covd(cov, os.path.join(arguments.output, 'covariance_matrix.png'))
        covariance.plot_full_covd(Lt_mat, os.path.join(arguments.output, 'cholesky_Lt.png'))

        # The full least squares analysis
        result = least_squares(residuals, x0=param0, verbose=True, bounds=[lb, ub], args=(data, gamma, lam),
                               x_scale=xscale)  # slip, strain, a, b, c

        print(result.x)
        model_pred = forward_model(result.x)
        model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
        inversion_utilities.data_model_misfit_plot(data, model_pred, faults,
                                                   os.path.join(arguments.output, expname+"_data_v_model.png"),
                                                   region=(-115.88, -115.65, 32.90, 33.05), s=20)
        np.savetxt(os.path.join(arguments.output, 'fitted_parameters.txt'), result.x,
                   header="Params slip(cm), strain(1e-5), plane, plane, reference")  # original vector

        full_param_vector = const_tau_params_to_full_params(result.x, configs["num_faults"])  # convert 43 to 81
        inversion_utilities.write_outputs(data, model_pred, full_param_vector, lam, gamma, arguments.output, expname,
                                          configs)
        d_misfit = data.LOS - model_pred.LOS
        rms_misfit = np.sqrt(np.mean(d_misfit**2))

        full_residuals = residuals(result.x, data, gamma, lam)  # full residual vector, with applied coefficients

        _, _ = inversion_utilities.plot_complete_residual_vector_and_results(full_residuals,
                                                                             data,
                                                                             model_pred,
                                                                             full_param_vector,
                                                                             faults,
                                                                             arguments.output,
                                                                             arguments.laplacian)
        return rms_misfit


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
    argms = parse_arguments()  # returns a Namespace
    print(argms)
    _misfit = do_main(argms)
