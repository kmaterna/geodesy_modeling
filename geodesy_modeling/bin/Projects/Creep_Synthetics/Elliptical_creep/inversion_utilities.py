
import numpy as np
import matplotlib.pyplot as plt
import os
from tectonic_utils.geodesy import insar_vector_functions, fault_vector_functions
from geodesy_modeling.datatypes import InSAR_1D_Object
from elastic_stresses_py import PyCoulomb
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso


def convert_xy_to_ll_insar1D(insar1D, configs):
    """
    Do the last step of the inversion, turning the model results into a useful InSAR1D object in lat-lon space
    """
    [real_lons, real_lats] = fault_vector_functions.xy2lonlat(insar1D.lon, insar1D.lat,
                                                              configs["zerolon"],
                                                              configs["zerolat"])
    lkvE, lkvN, lkvU = insar_vector_functions.flight_incidence_angles2look_vector(configs["flight_angle"],
                                                                                  configs["incidence_angle"],
                                                                                  look_direction='right')
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(real_lons, real_lats, insar1D.LOS,
                                                           insar1D.LOS_unc,
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvE),
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvN),
                                                           np.multiply(np.ones(np.shape(real_lons)), lkvU))
    return insar_data


def project_disp_points_into_insar1d(model_disp_points, flight_angle, incidence_angle):
    """
    Project some modeled displacements points into the LOS and package them as InSAR1D object, matching the data.

    :param model_disp_points: list of Disp Point objects
    :param flight_angle: float, in degrees
    :param incidence_angle: float, in degrees, angle of the look vector away from the vertical
    :return: InSAR_1D object
    """
    u = np.array([x.dE_obs for x in model_disp_points])
    v = np.array([x.dN_obs for x in model_disp_points])
    w = np.array([x.dU_obs for x in model_disp_points])
    lons = [x.lon for x in model_disp_points]
    lats = [x.lat for x in model_disp_points]
    los = insar_vector_functions.def3D_into_LOS(u, v, w, flight_angle=flight_angle, incidence_angle=incidence_angle,
                                                look_direction='right')
    los = np.multiply(los, -1000)  # insar1d object is in mm; sign convention is flipped
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(lons, lats, los,
                                                           np.zeros(np.shape(los)), np.zeros(np.shape(los)),
                                                           np.zeros(np.shape(los)), np.zeros(np.shape(los)))
    return insar_data


def read_data_and_faults(config, write_out_faults_file='used_faults.txt'):
    """ Read all the required data files and fault geometry files. """
    data = InSAR_1D_Object.inputs.inputs_txt(config["datafile"])  # read as insar1d object
    disp_points = data.get_disp_points()
    cart_disp_points = PyCoulomb.utilities.convert_ll2xy_disp_points(disp_points, config["zerolon"],
                                                                     config["zerolat"])  # ll to cart, for performance.
    faults = PyCoulomb.fault_slip_object.file_io.io_slippy.read_slippy_distribution(config["faultfile"])
    PyCoulomb.fault_slip_object.file_io.io_slippy.write_slippy_distribution(faults, write_out_faults_file)
    return data, cart_disp_points, faults


def data_model_misfit_plot(datapts, modelpts, best_params, faults, outname):
    """ Visualize the results of your inversion. """
    vmin, vmax = -15, 15
    residuals = datapts.LOS - modelpts.LOS
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
    # Show the data
    _im1 = axes[0].scatter(datapts.lon, datapts.lat, c=datapts.LOS, vmin=vmin, vmax=vmax)
    axes[0].set_title("Input Data")
    # Show the model
    im2 = axes[1].scatter(modelpts.lon, modelpts.lat, c=modelpts.LOS, vmin=vmin, vmax=vmax)
    axes[1].set_title("Model, Elliptical, depth0= %.2f km" % best_params[-1])  # this could be printing something diff.
    for item in faults:
        lons, lats = item.get_updip_corners_lon_lat()
        axes[1].plot(lons, lats)
    # Show the residuals
    axes[2].scatter(modelpts.lon, modelpts.lat, c=residuals, vmin=vmin, vmax=vmax)
    axes[2].set_title("Residuals, RMS=%.6f mm" % np.sqrt(np.mean(residuals**2)))
    cbar = fig.colorbar(im2, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('LOS Deformation (mm)')
    plt.savefig(outname)
    return


def write_outputs(data, model, fitted_params, lam, gamma, outdir, expname: str, configs):
    if len(fitted_params) > 6:
        # If there are a lot of parameters, let's make a plot of them. Should formalize this later.
        fitted_params = np.array(fitted_params).reshape(2, configs["num_faults"]).T
        f, axarr = plt.subplots(2, 1, figsize=(14, 10), dpi=300)
        axarr[0].plot(fitted_params[:, 0])
        axarr[0].set_ylabel('Slip (m)')
        axarr[1].plot(fitted_params[:, 1])
        axarr[1].set_ylabel('Depth (km)')
        plt.savefig(os.path.join(outdir, expname+"_model_params_results.png"))
    residuals = data.LOS - model.LOS
    rms = np.sqrt(np.mean(residuals**2))
    outdata = np.vstack((data.lon, data.lat, data.lkv_E, data.lkv_N, data.lkv_U, data.LOS, model.LOS)).T
    outfilename = os.path.join(outdir, expname+'_data_vs_model_predictions.txt')
    np.savetxt(outfilename, outdata, header="lon, lat, lkv_E, lkv_N, lkv_U, data, model")
    outfilename = os.path.join(outdir, expname + '_fitted_parameters.txt')
    np.savetxt(outfilename, fitted_params, header="Params slip(m) depth(km) reference")
    with open(os.path.join(outdir, expname+"_metrics.txt"), 'w') as outfile:
        outfile.write("RMS: %f mm\n" % rms)
        outfile.write("lam (tikhonov): %f\n" % lam)
        outfile.write("gamma (smoothing): %f\n" % gamma)
    return
