#!/usr/bin/env python

"""
A model using the non-linear inversion that models slip as elliptical slip distributions.

# How to do this: create columns of fault elements.
# Create a function that takes Z and A and gives displacements
# A is non-linearly related to displacements.
Param vector is "SLIP -- DEPTH" for every fault patch.
We have now included a reference pixel parameter.
"""
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb.fault_slip_triangle import triangle_okada
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.optimize import approx_fprime
from scipy.linalg import cholesky, solve_triangular
import elliptical_utilities  # local import
import inversion_utilities  # local import
from geodesy_modeling.datatypes.InSAR_1D_Object import covariance
import os


top_configs = {"datafile": "../../Prepare_Data/InSAR_Data/downsampled_data.txt",
               "faultfile": "../../../Get_Fault_Model/model_fault_patches_39.txt",
               "fieldfile": "../../Prepare_Data/Field_Data/slip_minimums.txt",
               "cov_parameters": "../../Prepare_Data/InSAR_Data/covariance_matrix_parameters.txt",
               "max_disloc_depth": 5.0,
               "num_faults": 39,
               "top_depth": 0,
               "zerolon": -115.0,
               "zerolat": 33,
               "fault_width": 5.0,
               "bbox": [-115.15, -114.9, 32.9, 33.15],
               "sampling_interval": 0.15,
               "flight_angle": 190,
               "incidence_angle": 37}
default_params = PyCoulomb.configure_calc.get_lightweight_config_params(mu=30e9, lame1=30e9, B=0)


def elastic_model(params, cart_disp_points, faults, configs):
    """
    The forward model. The mesh geometry is known separately.

    :param cart_disp_points: used only for location of points, in cartesian coordinates
    :param params: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :param configs: dictionary containing many geometry parameters
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    surface_slip, bottom_depth = params[0], params[1]  # for small inversion. KEY LINE.
    a, b, c = params[2], params[3], params[4]  # the three parameters for the plane fitting to the data
    fault_list = []
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        used_depths, widths, modeled_slip = elliptical_utilities.get_tuned_depth_arrays(configs['top_depth'],
                                                                                        configs['max_disloc_depth'],
                                                                                        configs["sampling_interval"],
                                                                                        surface_slip,
                                                                                        bottom_depth)
        for mslip, mwidth, depth in zip(modeled_slip, widths, used_depths):   # create all the source patches
            # Change the slip to elliptical slip on each fault patch
            one_fault = fso.FaultSlipObject(strike=fault.strike, dip=89.99, length=fault.length,
                                            width=mwidth, lon=fault.lon,
                                            lat=fault.lat, depth=depth, rake=180, slip=mslip)
            fault_list.append(one_fault)
    source_object = fso.fault_object_to_coulomb_fault(fault_list, zerolon_system=configs["zerolon"],
                                                      zerolat_system=configs["zerolat"])
    inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                    zerolon=configs['zerolon'],
                                                                                    zerolat=configs['zerolat'],
                                                                                    bbox=configs["bbox"])
    model_disp_points = triangle_okada.compute_cartesian_def_tris(inputs, default_params, cart_disp_points)  # run okada
    insar_1d_model = inversion_utilities.project_disp_points_into_insar1d(model_disp_points, configs["flight_angle"],
                                                                          configs["incidence_angle"])  # cartesian space
    # Implement a planar fit to the data, which has three planar parameters for x and y and a constant offset
    insar_1d_model = insar_1d_model.subtract_ramp(a, b, c)
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm, with lon/lat in cartesian space


def invert_data_2param(configs):
    data, cart_disp_points, faults = inversion_utilities.read_data_and_faults(configs, os.path.join('Two_Param_Tests',
                                                                                                    'used_faults.txt'))
    param0 = [0.020, 2.0, 0, 0, 0]
    lb, ub = [0.0, 0.0, -50, -50, -5.0], [0.050, 5.0, 50, 50, 5.0]  # last three parameters are a, b, and c for plane
    lam = 1.7  # Minimum-norm Tikhonov smoothing regularization strength

    # We are going to determine the real covariance matrix, and compute the inverse by
    # one of the two triangular matrices in the Cholesky decomposition.
    L, sigma = covariance.read_covd_parameters(configs["cov_parameters"])
    cov = covariance.build_Cd(data, sigma, L)
    covariance.plot_full_covd(cov, os.path.join('Two_Param_Tests', 'covariance_matrix.png'))

    # Whitening the covariance matrix through Cholesky decomposition
    Lt_mat = cholesky(cov, lower=True, check_finite=False)
    covariance.plot_full_covd(Lt_mat, os.path.join('Two_Param_Tests', 'cholesky_Lt.png'))

    def make_data_whitener(Cd):
        # Cd must be SPD. If near-singular, consider jitter or SVD.
        Ld = cholesky(Cd, lower=True, check_finite=False)  # Cd = Ld Ld^T

        # Wd @ x == solve(Ld, x) ⇒ Wd = Ld^{-1} without forming it
        def Wd_apply0(x):
            # handle 1D vector or 2D matrix (apply to columns)
            return solve_triangular(Ld, x, lower=True, check_finite=False)
        return Wd_apply0

    Wd_apply = make_data_whitener(cov)

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults, configs)
        return insar_1d_model

    def residuals(m, data0, lam0):
        data_misfit = Wd_apply(forward_model(m).LOS - data0.LOS)  # normalize the misfit by the sqrt(cov_matrix)
        tikhonov = lam0 * m  # Tikhonov (minimum-norm) regularization
        return np.concatenate((data_misfit, tikhonov))

    expname = "two_param_smooth7"
    result = least_squares(residuals, x0=param0, verbose=2, bounds=[lb, ub], args=(data, lam),
                           method='trf', x_scale=[0.01, 1.0, 5.0, 5.0, 1.0])  # slip, dep, ramp-ramp-refpix
    print(result.x)
    model_pred = forward_model(result.x)
    model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
    inversion_utilities.data_model_misfit_plot(data, model_pred, faults,
                                               os.path.join('Two_Param_Tests', expname+"_data_v_model.png"))
    inversion_utilities.write_outputs(data, model_pred, result.x, lam, 0, 'Two_Param_Tests', expname, configs)

    # testcase_params = param0
    # simple_model = forward_model(testcase_params)  # InSAR1D object
    # simple_model = inversion_utilities.convert_xy_to_ll_insar1D(simple_model, configs["zerolon"], configs["zerolat"])
    # inversion_utilities.data_model_misfit_plot(data, simple_model, faults, "two_param_model.png")
    # inversion_utilities.write_outputs(data, simple_model, testcase_params, 0, 0, 'two_param')
    return


def investigate_jacobian_2param(configs):
    data, cart_disp_points, faults = inversion_utilities.read_data_and_faults(configs)
    param0 = [0.020, 2.0]
    lam = 0  # Minimum norm Tikhonov smoothing regularization strength

    # Establish forward model and cost function
    def forward_model_rectangle(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults, configs)
        return insar_1d_model

    # Establish forward model and cost function
    def forward_model_ellipse(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults, configs)
        return insar_1d_model

    def residuals(m, data0, lam0):
        data_misfit = forward_model_ellipse(m).LOS - data0.LOS
        tikhonov = lam0 * m  # Tikhonov (minimum-norm) regularization
        return np.concatenate((data_misfit, tikhonov))

    model_pred = forward_model_ellipse([0.02, 2.75])
    model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
    inversion_utilities.data_model_misfit_plot(data, model_pred, faults, "ellipse_data_v_model.png")
    model_pred = forward_model_rectangle([0.02, 2.20])
    model_pred = inversion_utilities.convert_xy_to_ll_insar1D(model_pred, configs)
    inversion_utilities.data_model_misfit_plot(data, model_pred, faults, "rectangle_data_v_model.png")

    def wrapped_residual(x):
        return residuals(x, data, lam)

    r = residuals(param0, data, lam)
    print("residual norm:", np.linalg.norm(r))
    print("residual max:", np.max(np.abs(r)))

    J = approx_fprime(param0, wrapped_residual)
    print("Jacobian norm:", np.linalg.norm(J))
    print("Jacobian Shape:", np.shape(J))
    print("Jacobian:", J)

    f, axarr = plt.subplots(1, 2, figsize=(9, 5), dpi=300)
    slip_jacobian = J[0:1407, 0]
    depth_jacobian = J[0:1407, 1]
    x_vals = [x.lon for x in cart_disp_points]
    y_vals = [x.lat for x in cart_disp_points]
    print(np.shape(x_vals), np.shape(y_vals), np.shape(slip_jacobian))
    im1 = axarr[0].scatter(x_vals, y_vals, c=slip_jacobian, cmap='viridis')
    axarr[0].set_title("Jacobian with Slip")
    _cbar1 = f.colorbar(im1, ax=axarr[0])
    im2 = axarr[1].scatter(x_vals, y_vals, c=depth_jacobian, cmap='viridis')
    axarr[1].set_title("Jacobian with Depth")
    _cbar2 = f.colorbar(im2, ax=axarr[1])
    plt.savefig("Jacobian_spatial.png")

    f, axarr = plt.subplots(2, 1, figsize=(10, 8), dpi=300)
    axarr[0].plot(J[:, 0])
    axarr[1].plot(J[:, 1])
    plt.savefig("Jacobian.png")

    return


def do_main():
    invert_data_2param(top_configs)  # create a smaller inversion for testing speed-ups, like GPU and desktop work
    # elliptical_utilities.show_elliptical_distribution(top_configs["top_depth"], top_configs["max_disloc_depth"],
    #                                                   top_configs["sampling_interval"], surface_slip=0.02,
    #                                                   bottom_depth=2.3, outfile='Plot_rectangles_within_ellipse.png')
    # investigate_jacobian_2param(top_configs)
    return


if __name__ == "__main__":
    do_main()
