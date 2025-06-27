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
from elastic_stresses_py.PyCoulomb.fault_slip_triangle import triangle_okada
import numpy as np
import matplotlib.pyplot as plt
from Tectonic_Utils.geodesy import insar_vector_functions, fault_vector_functions
from geodesy_modeling.datatypes import InSAR_1D_Object
from scipy.optimize import least_squares
from scipy.optimize import approx_fprime


configs = {"datafile": "../Prepare_Data/InSAR_Data/downsampled_data.txt",
           "faultfile": "../../Get_Fault_Model/model_fault_patches_39.txt",
           "fieldfile": "../Prepare_Data/Field_Data/slip_minimums.txt",
           "disloc_depth": 5.0,
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


def project_model_disp_points_into_insar1d(model_disp_points):
    """
    Project some modeled displacements points into the LOS and package them as InSAR1D object.

    :param model_disp_points: list of Disp Point objects
    :return: InSAR_1D object
    """
    u = np.array([x.dE_obs for x in model_disp_points])
    v = np.array([x.dN_obs for x in model_disp_points])
    w = np.array([x.dU_obs for x in model_disp_points])
    lons = [x.lon for x in model_disp_points]
    lats = [x.lat for x in model_disp_points]
    los = insar_vector_functions.def3D_into_LOS(u, v, w, flight_angle=configs["flight_angle"],
                                                incidence_angle=configs["incidence_angle"],
                                                look_direction='right')
    los = np.multiply(los, -1000)  # insar1d object is in mm; sign convention is flipped
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(lons, lats, los,
                                                           np.zeros(np.shape(los)), np.zeros(np.shape(los)),
                                                           np.zeros(np.shape(los)), np.zeros(np.shape(los)))
    return insar_data


def convert_xy_to_ll_insar1D(insar1D, zerolon, zerolat):
    """
    Do the last step, turning the model results into a useful InSAR1D object in lat-lon space
    """
    [real_lons, real_lats] = fault_vector_functions.xy2lonlat(insar1D.lon, insar1D.lat, zerolon, zerolat)
    lkvE, lkvN, lkvU = insar_vector_functions.flight_incidence_angles2look_vector(configs["flight_angle"],
                                                                                  configs["incidence_angle"],
                                                                                  look_direction='right')
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(real_lons, real_lats, insar1D.LOS,
                                                           insar1D.LOS_unc,
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvE),
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvN),
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvU))
    return insar_data


def applied_slip(depths, So, a):
    """
    Create an elliptical slip distribution based on surface slip So and bottom depth a.
    """
    # Depths are km.
    input_slip = np.sqrt(So*So * (1 - (np.square(depths)/(a**2))))
    used_depths = depths[~np.isnan(input_slip)]
    input_slip = input_slip[~np.isnan(input_slip)]
    return used_depths, input_slip


def get_width_of_depth(starting_width, depth):
    """ The function that increases size of rectangles going deeper. """
    width = starting_width + 0.27 * depth
    return width


def get_tuned_depth_array(top_depth, bottom_depth, starting_width, surface_slip, bottom_depth_slip):
    """
    Produce an array of depth values that are associated with the
    """
    counter = top_depth
    top_depths, widths = [], []
    while counter <= bottom_depth:
        top_depths.append(counter)
        width = get_width_of_depth(starting_width, counter)
        widths.append(width)
        counter += width
    top_depths, widths = np.array(top_depths), np.array(widths)
    middle_depths = np.array([x+y/2 for x, y in zip(top_depths, widths)])
    exp_depths, exp_slip = applied_slip(middle_depths, surface_slip, bottom_depth_slip)  # middle depths
    exp_widths = widths[0:len(exp_depths)]
    exp_depths = top_depths[0:len(exp_depths)]
    return exp_depths, exp_widths, exp_slip


def plot_rectangles(ax, top_depths, widths, slips, color):
    for i in range(len(top_depths)):
        x = [0, slips[i], slips[i], 0]
        y = [top_depths[i], top_depths[i], top_depths[i]+widths[i], top_depths[i]+widths[i]]
        ax.plot(x, y, color=color)
    return ax


def show_sampling():
    surface_slip = 0.02
    bottom_depth = 2.3
    sample_depths = np.arange(configs["top_depth"], configs["disloc_depth"], configs["sampling_interval"])

    # # Create a model of the fault using non-equally-spaced rectangles
    exp_depths, exp_widths, exp_slip = get_tuned_depth_array(configs['top_depth'], configs['disloc_depth'],
                                                             configs["sampling_interval"], surface_slip, bottom_depth)

    # Create a model of the fault using equally-spaced rectangles
    rect_depths, rect_slip = applied_slip(sample_depths, surface_slip, bottom_depth)  # top depths
    rect_widths = configs["sampling_interval"]*np.ones(np.shape(rect_depths))

    model_depths = np.arange(configs["top_depth"], configs["disloc_depth"], 0.05)
    model_depths, model_slip = applied_slip(model_depths, surface_slip, bottom_depth)
    f, axarr = plt.subplots(1, 2, dpi=300, figsize=(12, 8))
    axarr[0].plot(model_slip, model_depths)
    plot_rectangles(axarr[0], rect_depths, rect_widths, rect_slip, color='red')
    axarr[0].invert_yaxis()
    axarr[0].set_title("N=" + str(len(rect_depths)))
    axarr[0].set_xlabel("Slip (m)")
    axarr[0].set_ylabel("Depth (km)")

    axarr[1].plot(model_slip, model_depths)
    plot_rectangles(axarr[1], exp_depths, exp_widths, exp_slip, color='black')
    axarr[1].invert_yaxis()
    axarr[1].set_title("N=" + str(len(exp_depths)))
    axarr[1].set_xlabel("Slip (m)")
    axarr[1].set_ylabel("Depth (km)")

    plt.savefig("Modeling_ellipse.png")
    return


def elastic_model(param_vector, cart_disp_points, faults):
    """
    The mesh geometry is known separately.

    :param cart_disp_points: used only for location of points, in cartesian coordinates
    :param param_vector: vector of surface_slip values and depth-values for the various elliptical slip distributions
    :param faults: used for sources; the fault slip will be re-set to an elliptical slip distribution
    :return: InSAR_1D_object, a matching object to the data structure. The LOS values contain the model predictions.
    """
    # depths = np.arange(configs["top_depth"], configs["fault_width"], configs['sampling_interval'])
    fault_list = []
    for j, fault in enumerate(faults):
        # Here we could use one value for the entire fault length, or use a spatial distribution of values.
        surface_slip, bottom_depth = param_vector[j], param_vector[j+configs["num_faults"]]  # for real faults
        # surface_slip, bottom_depth = param_vector[0], param_vector[1]  # for 2-param inversion
        used_depths, widths, modeled_slip = get_tuned_depth_array(configs['top_depth'], configs['disloc_depth'],
                                                                  configs["sampling_interval"], surface_slip,
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
    insar_1d_model = project_model_disp_points_into_insar1d(model_disp_points)  # this will be in cartesian space
    print("Faults and points: %d and %d" % (len(fault_list), len(insar_1d_model.LOS)))
    return insar_1d_model  # in mm


def set_up_initial_params_and_bounds():
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
    if len(fitted_params) > 2:
        fitted_params = np.array(fitted_params).reshape(2, configs["num_faults"]).T
        f, axarr = plt.subplots(2, 1, figsize=(14, 10), dpi=300)
        axarr[0].plot(fitted_params[:, 0])
        axarr[0].set_ylabel('Slip (m)')
        axarr[1].plot(fitted_params[:, 1])
        axarr[1].set_ylabel('Depth (km)')
        plt.savefig(expname+"_model_params_results.png")
    residuals = data.LOS - model.LOS
    rms = np.sqrt(np.mean(residuals**2))
    outdata = np.vstack((data.lon, data.lat, data.lkv_E, data.lkv_N, data.lkv_U, data.LOS, model.LOS)).T
    np.savetxt(expname+"_data_vs_model_predictions.txt", outdata, header="lon, lat, lkv_E, lkv_N, lkv_U, data, model")
    np.savetxt(expname+"_fitted_parameters.txt", fitted_params, header="slip(m) depth(km)")
    with open(expname+"_metrics.txt", 'w') as outfile:
        outfile.write("RMS: %f mm\n" % rms)
        outfile.write("lam (tikhonov): %f\n" % lam)
        outfile.write("gamma (smoothing): %f\n" % gamma)
    return


def read_data_and_faults(config):
    data = InSAR_1D_Object.inputs.inputs_txt(config["datafile"])  # read as insar1d object
    disp_points = data.get_disp_points()
    cart_disp_points = PyCoulomb.utilities.convert_ll2xy_disp_points(disp_points, config["zerolon"],
                                                                     config["zerolat"])  # ll to cart.
    faults = PyCoulomb.fault_slip_object.file_io.io_slippy.read_slippy_distribution(config["faultfile"])
    PyCoulomb.fault_slip_object.file_io.io_slippy.write_slippy_distribution(faults, "used_faults.txt")
    return data, cart_disp_points, faults


def invert_data():
    # We can vary the smoothing parameter and determine the most appropriate one through L-curve analysis
    data, cart_disp_points, faults = read_data_and_faults(configs)

    # lam_array = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10, 50, 100]
    lam_array = [1]
    # lam = 1  # Minimum norm Tikhonov smoothing regularization strength
    param0, lb, ub = set_up_initial_params_and_bounds()  # create initial parameter vector

    for lam in lam_array:
        gamma = 50  # Laplacian Regularization strength

        # Establish forward model and cost function
        def forward_model(params):  # params = vector of size 82, 41 slips and 41 depths
            insar_1d_model = elastic_model(params, cart_disp_points, faults)
            return insar_1d_model

        # Smoothing residuals (first differences) on first half of vector
        def smoothing_residuals(m2, gamma):
            return gamma * np.diff(m2)

        def residuals(m, data, gamma, lam):
            data_misfit = forward_model(m).LOS - data.LOS
            m2 = m[0:configs["num_faults"]]  # slip is the first half of the model vector
            mdepth = m[configs["num_faults"]:]
            tikhonov = lam * mdepth  # Tikhonov (minimum-norm) regularization
            smoothing = smoothing_residuals(m2, gamma)
            return np.concatenate((data_misfit, smoothing, tikhonov))

        expname = 'Multi_Param/smoothing_'+str(lam)
        result = least_squares(residuals, x0=param0, verbose=True, bounds=[lb, ub], args=(data, gamma, lam))  # slip, z
        print(result.x)
        model_pred = forward_model(result.x)
        model_pred = convert_xy_to_ll_insar1D(model_pred, configs["zerolon"], configs["zerolat"])
        create_data_model_misfit_plot(data, model_pred, result.x, faults, expname+"_data_v_model.png")
        write_outputs(data, model_pred, result.x, lam, gamma, expname)

    # testcase_params = param0
    # simple_model = forward_model(testcase_params)  # InSAR1D object
    # simple_model = convert_xy_to_ll_insar1D(simple_model, configs["zerolon"], configs["zerolat"])
    # create_data_model_misfit_plot(data, simple_model, testcase_params, faults, "testcase_data_v_model.png")
    # write_outputs(data, simple_model, testcase_params, lam, gamma, 'testcase')
    return


def invert_data_2param():
    data, cart_disp_points, faults = read_data_and_faults(configs)
    param0 = [0.020, 2.0]
    lb, ub = [0.0, 0.0], [0.050, 5.0]
    lam = 7  # Minimum norm Tikhonov smoothing regularization strength

    # Establish forward model and cost function
    def forward_model(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults)
        return insar_1d_model

    def residuals(m, data, lam):
        data_misfit = forward_model(m).LOS - data.LOS
        tikhonov = lam * m  # Tikhonov (minimum-norm) regularization
        return np.concatenate((data_misfit, tikhonov))

    expname = "two_param_smooth8"
    result = least_squares(residuals, x0=param0, verbose=2, bounds=[lb, ub], args=(data, lam),
                           method='trf', x_scale=[0.01, 1.0])  # slip, dep
    print(result.x)
    model_pred = forward_model(result.x)
    model_pred = convert_xy_to_ll_insar1D(model_pred, configs["zerolon"], configs["zerolat"])
    create_data_model_misfit_plot(data, model_pred, result.x, faults, expname+"_data_v_model.png")
    write_outputs(data, model_pred, result.x, lam, 0, expname)

    # testcase_params = param0
    # simple_model = forward_model(testcase_params)  # InSAR1D object
    # simple_model = convert_xy_to_ll_insar1D(simple_model, configs["zerolon"], configs["zerolat"])
    # create_data_model_misfit_plot(data, simple_model, testcase_params, faults, "two_param_data_v_model.png")
    # write_outputs(data, simple_model, testcase_params, 0, 0, 'two_param')
    return


def investigate_jacobian_2param():
    data, cart_disp_points, faults = read_data_and_faults(configs)
    param0 = [0.020, 2.0]
    lb, ub = [0.0, 0.0], [0.050, 5.0]
    lam = 0  # Minimum norm Tikhonov smoothing regularization strength

    # Establish forward model and cost function
    def forward_model_rectangle(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults)
        return insar_1d_model

    # Establish forward model and cost function
    def forward_model_ellipse(params):  # params = vector of size 82, 41 slips and 41 depths
        insar_1d_model = elastic_model(params, cart_disp_points, faults)
        return insar_1d_model

    def residuals(m, data, lam):
        data_misfit = forward_model_ellipse(m).LOS - data.LOS
        tikhonov = lam * m  # Tikhonov (minimum-norm) regularization
        return np.concatenate((data_misfit, tikhonov))

    model_pred = forward_model_ellipse([0.02, 2.75])
    model_pred = convert_xy_to_ll_insar1D(model_pred, configs["zerolon"], configs["zerolat"])
    create_data_model_misfit_plot(data, model_pred, [0.02, 2.75], faults, "ellipse_data_v_model.png")
    model_pred = forward_model_rectangle([0.02, 2.20])
    model_pred = convert_xy_to_ll_insar1D(model_pred, configs["zerolon"], configs["zerolat"])
    create_data_model_misfit_plot(data, model_pred, [0.02, 2.20], faults, "rectangle_data_v_model.png")

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
    invert_data()   # Build upon the Inversion stress_driven, using InSAR1D object as data, invert for S0 and Z
    # invert_data_2param()  # create a smaller inversion for testing speed-ups, like GPU and desktop work
    # investigate_jacobian_2param()
    # show_sampling()
    return


if __name__ == "__main__":
    do_main()
