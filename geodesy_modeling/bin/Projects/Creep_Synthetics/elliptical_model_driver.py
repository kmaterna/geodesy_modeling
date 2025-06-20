#!/usr/bin/env python

"""
A model using the non-linear inversion that models slip as elliptical slip distributions.

# How to do this: create columns of fault elements.
# Create a function that takes Z and A and gives displacements
# A is non-linearly related to displacements.
Param vector is "SLIP -- DEPTH" for every fault patch.
"""
import elastic_stresses_py.PyCoulomb.fault_slip_object.file_io.io_slippy
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso
import numpy as np
import matplotlib.pyplot as plt
from Tectonic_Utils.geodesy import insar_vector_functions, fault_vector_functions
from geodesy_modeling.datatypes import InSAR_1D_Object
from scipy.optimize import least_squares
import sys


configs = {"datafile": "../Prepare_Data/downsampled_data.txt",
           "faultfile": "../../Get_Fault_Model/model_fault_patches_41.txt",
           "disloc_depth": 5.0,
           "num_faults": 41,
           "top_depth": 0,
           "zerolon": -115.0,
           "zerolat": 33,
           "fault_width": 5.0,
           "bbox": [-115.15, -114.9, 32.9, 33.15],
           "sampling_interval": 0.15,
           "flight_angle": 190,
           "incidence_angle": 37}
default_params = PyCoulomb.configure_calc.get_lightweight_config_params(mu=30e9, lame1=30e9, B=0)


def project_model_disp_points_into_insar1d(model_disp_points):
    """
    Project some modeled displacements points into the LOS and package them as InSAR1D object.

    :param model_disp_points: list of Disp Point objects
    :return: InSAR_1D object
    """
    u = [x.dE_obs for x in model_disp_points]
    v = [x.dN_obs for x in model_disp_points]
    w = [x.dU_obs for x in model_disp_points]
    lons = [x.lon for x in model_disp_points]
    lats = [x.lat for x in model_disp_points]
    los = insar_vector_functions.def3D_into_LOS(np.array(u), np.array(v), np.array(w),
                                                flight_angle=configs["flight_angle"],
                                                incidence_angle=configs["incidence_angle"],
                                                look_direction='right')
    los = np.multiply(los, -1000)  # insar1d object is in mm, sign convention is flipped
    lkvE, lkvN, lkvU = insar_vector_functions.flight_incidence_angles2look_vector(configs["flight_angle"],
                                                                                  configs["incidence_angle"],
                                                                                  look_direction='right')
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(lons, lats, los,
                                                           np.zeros(np.shape(los)),
                                                           np.multiply(np.ones(np.shape(lons)), lkvE),
                                                           np.multiply(np.ones(np.shape(lons)), lkvN),
                                                           np.multiply(np.ones(np.shape(lons)), lkvU))
    return insar_data


def applied_slip(depths, So, a):
    """
    Create an elliptical slip distribution based on surface slip So and bottom depth a.
    """
    # Depths are km.
    input_slip = np.sqrt(So*So * (1 - (np.square(depths)/(a**2))))
    input_slip[np.isnan(input_slip)] = 0
    return input_slip


def elastic_model(param_vector, disp_points, faults):
    """
    The mesh geometry is known separately.

    :param disp_points: used only for location of points
    :param param_vector: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    depths = np.arange(configs["top_depth"], configs["fault_width"], configs['sampling_interval'])
    fault_list = []
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        surface_slip, bottom_depth = param_vector[j], param_vector[j+configs["num_faults"]]
        modeled_slip = applied_slip(depths, surface_slip, bottom_depth)
        for i, depth in enumerate(depths):   # create all the source patches
            if modeled_slip[i] == 0:
                continue
            one_fault = fso.FaultSlipObject(strike=fault.strike, dip=89.99, length=fault.length,
                                            width=configs["sampling_interval"], lon=fault.lon,
                                            lat=fault.lat, depth=depth, rake=180, slip=modeled_slip[i])
            # Change the slip to elliptical slip
            fault_list.append(one_fault)
    source_object = fso.fault_object_to_coulomb_fault(fault_list, zerolon_system=configs["zerolon"],
                                                      zerolat_system=configs["zerolat"])
    inputs = PyCoulomb.inputs_object.input_obj.configure_default_displacement_input(source_object=source_object,
                                                                                    zerolon=configs['zerolon'],
                                                                                    zerolat=configs['zerolat'],
                                                                                    bbox=configs["bbox"])
    model_disp_points = PyCoulomb.run_dc3d.compute_ll_def(inputs, default_params, disp_points)  # run okada
    insar_1d_model = project_model_disp_points_into_insar1d(model_disp_points)
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm


def set_up_initial_params_and_bounds():
    # Set up the original parameter vector and the upper bounds and lower bounds
    param0 = []  # param vector = [slip slip slip ..... slip depth depth depth.....]
    for i in range(configs["num_faults"]):  # initial guess for the slip is 20 mm (m)
        param0.append(0.02)
    for i in range(configs["num_faults"]):  # initial guess for the lower-depth km is 1 km
        param0.append(2.0)
    lower_bound = np.zeros((configs["num_faults"]*2, ))   # lower bound on slip and depth is zero
    upper_bound = np.multiply(5, np.ones((configs["num_faults"]*2, )))  # upper bound on slip is 50 mm, depth is 5 km
    upper_bound[0:configs["num_faults"]] = 0.050
    return param0, lower_bound, upper_bound


def create_data_model_misfit_plot(datapts, modelpts, best_params, faults, outname):
    """ visualize the results of your inversion. """
    vmin, vmax = -15, 15
    residuals = datapts.LOS - modelpts.LOS
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
    _im1 = axes[0].scatter(datapts.lon, datapts.lat, c=datapts.LOS, vmin=vmin, vmax=vmax)
    axes[0].set_title("Input Data")
    im2 = axes[1].scatter(modelpts.lon, modelpts.lat, c=modelpts.LOS, vmin=vmin, vmax=vmax)
    axes[1].set_title("Model, Elliptical, depth0= %.2f km" % best_params[-1])
    for item in faults:
        lons, lats = item.get_updip_corners_lon_lat()
        axes[1].plot(lons, lats)
    axes[2].scatter(modelpts.lon, modelpts.lat, c=residuals, vmin=vmin, vmax=vmax)
    axes[2].set_title("Residuals, RMS=%.6f mm" % np.sqrt(np.mean(residuals**2)))
    cbar = fig.colorbar(im2, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('LOS Deformation (mm)')
    plt.savefig(outname)
    return


def write_outputs(data, model, fitted_params, lam, gamma, expname):
    fitted_params = np.array(fitted_params).reshape(2, configs["num_faults"]).T
    residuals = data.LOS - model.LOS
    rms = np.sqrt(np.mean(residuals**2))
    outdata = np.vstack((data.lon, data.lat, data.lkv_E, data.lkv_N, data.lkv_U, data.LOS, model.LOS)).T
    np.savetxt(expname+"_data_vs_model_predictions.txt", outdata, header="lon, lat, lkv_E, lkv_N, lkv_U, data, model")
    np.savetxt(expname+"_fitted_parameters.txt", fitted_params, header="slip(m) depth(km)")
    f, axarr = plt.subplots(2, 1, figsize=(14, 10), dpi=300)
    axarr[0].plot(fitted_params[:, 0])
    axarr[0].set_ylabel('Slip (m)')
    axarr[1].plot(fitted_params[:, 1])
    axarr[1].set_ylabel('Depth (km)')
    plt.savefig(expname+"_model_params_results.png")
    with open(expname+"_metrics.txt", 'w') as outfile:
        outfile.write("RMS: %f mm\n" % rms)
        outfile.write("lam (tikhonov): %f\n" % lam)
        outfile.write("gamma (smoothing): %f\n" % gamma)
    return


def invert_data():
    data = InSAR_1D_Object.inputs.inputs_txt(configs["datafile"])  # read as insar1d object
    disp_points = data.get_disp_points()
    faults = PyCoulomb.fault_slip_object.file_io.io_slippy.read_slippy_distribution(configs["faultfile"])
    PyCoulomb.fault_slip_object.file_io.io_slippy.write_slippy_distribution(faults, "used_faults.txt")
    lam = 10  # Minimum norm Tikhonov smoothing regularization strength
    gamma = 10  # Laplacian regularization strength
    param0, lb, ub = set_up_initial_params_and_bounds()

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size 102, 51 slips and 51 depths
        insar_1d_model = elastic_model(params, disp_points, faults)
        return insar_1d_model

    # Smoothing residuals (first differences) on first half of vector
    def smoothing_residuals(m2, gamma):
        return gamma * np.diff(m2)

    def residuals(m, data, gamma, lam):
        data_misfit = forward_model(m).LOS - data.LOS
        m2 = m[0:configs["num_faults"]]  # slip is the first half of the model vector
        tikhonov = lam * m  # Tikhonov (minimum-norm) regularization
        smoothing = smoothing_residuals(m2, gamma)
        return np.concatenate((data_misfit, smoothing, tikhonov))

    result = least_squares(residuals, x0=param0, verbose=True, bounds=[lb, ub], args=(data, gamma, lam))  # slip, depth
    print(result.x)
    model_pred = forward_model(result.x)
    create_data_model_misfit_plot(data, model_pred, result.x, faults, "inverse_data_v_model.png")
    write_outputs(data, model_pred, result.x, lam, gamma, 'inverse')

    testcase_params = param0
    simple_model = forward_model(testcase_params)
    create_data_model_misfit_plot(data, simple_model, testcase_params, faults, "testcase_data_v_model.png")
    write_outputs(data, simple_model, testcase_params, lam, gamma, 'testcase')
    return


def do_main():
    invert_data()   # Build upon the Inversion stress_driven, using InSAR1D object as data, invert for S0 and Z
    return


if __name__ == "__main__":
    do_main()
