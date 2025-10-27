
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
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(real_lons,
                                                           real_lats,
                                                           insar1D.LOS,
                                                           insar1D.LOS_unc,
                                                           insar1D.lkv_E,
                                                           insar1D.lkv_N,
                                                           insar1D.lkv_U)
    return insar_data


def project_disp_points_into_insar1d(model_disp_points, data_1d):
    """
    Project some modeled displacements points into the LOS and package them as InSAR1D object, matching the data.

    :param model_disp_points: list of Disp Point objects
    :param data_1d: insar1d object that contains look vectors, with the same shape as model_disp_points
    :return: InSAR_1D object
    """
    u = np.array([x.dE_obs for x in model_disp_points])
    v = np.array([x.dN_obs for x in model_disp_points])
    w = np.array([x.dU_obs for x in model_disp_points])
    lons = [x.lon for x in model_disp_points]
    lats = [x.lat for x in model_disp_points]
    lkvE = data_1d.lkv_E
    lkvN = data_1d.lkv_N
    lkvU = data_1d.lkv_U
    flight_angles, incidence_angles = insar_vector_functions.look_vector2flight_incidence_angles(lkvE, lkvN, lkvU,
                                                                                                 look_direction='right')
    los = insar_vector_functions.def3D_into_LOS(u, v, w, flight_angle=flight_angles, incidence_angle=incidence_angles,
                                                look_direction='right')
    los = np.multiply(los, -1000)  # insar1d object is in mm; sign convention is flipped
    insar_data = InSAR_1D_Object.class_model.Insar1dObject(lons, lats, los,
                                                           np.zeros(np.shape(los)), data_1d.lkv_E,
                                                           data_1d.lkv_N, data_1d.lkv_U)
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


def data_model_misfit_plot(datapts, modelpts, faults, outname, region, s=20):
    """ Visualize the results of your inversion. """
    vmin, vmax = -15, 15
    residuals = datapts.LOS - modelpts.LOS
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
    # Show the data
    _im1 = axes[0].scatter(datapts.lon, datapts.lat, c=datapts.LOS, vmin=vmin, vmax=vmax, s=s)
    axes[0].set_title("Input Data")
    axes[0].set_xlim([region[0], region[1]])
    axes[0].set_ylim([region[2], region[3]])
    # Show the model
    im2 = axes[1].scatter(modelpts.lon, modelpts.lat, c=modelpts.LOS, vmin=vmin, vmax=vmax, s=s)
    axes[1].set_title("Model, Elliptical")
    # for item in faults:
    #     lons, lats = item.get_updip_corners_lon_lat()
    #     axes[1].plot(lons, lats, linewidth=0.5)
    axes[1].set_xlim([region[0], region[1]])
    axes[1].set_ylim([region[2], region[3]])
    # Show the residuals
    axes[2].scatter(modelpts.lon, modelpts.lat, c=residuals, vmin=vmin, vmax=vmax, s=s)
    axes[2].set_title("Residuals, RMS=%.6f mm" % np.sqrt(np.mean(residuals**2)))
    # for item in faults:
    #     lons, lats = item.get_updip_corners_lon_lat()
    #     axes[2].plot(lons, lats, linewidth=0.5)
    axes[2].set_xlim([region[0], region[1]])
    axes[2].set_ylim([region[2], region[3]])
    cbar = fig.colorbar(im2, ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
    cbar.set_label('LOS Deformation (mm)')
    plt.savefig(outname)
    return


def write_outputs(data, model, fitted_params, lam, gamma, outdir, expname: str, configs):
    num_faults = configs["num_faults"]
    if len(fitted_params) > 6:
        # If there are a lot of parameters, let's make a plot of them. Should formalize this later.
        fault_fitted_params = np.array(fitted_params)[0:2*num_faults].reshape(2, num_faults).T
        f, axarr = plt.subplots(2, 1, figsize=(14, 10), dpi=300)
        axarr[0].plot(fault_fitted_params[:, 0])
        axarr[0].set_ylabel('Slip (cm)')
        axarr[1].plot(fault_fitted_params[:, 1])
        axarr[1].set_ylabel('Depth (km)')
        plt.savefig(os.path.join(outdir, expname+"_model_params_results.png"))
    fitted_params = np.array(fitted_params)
    residuals = data.LOS - model.LOS
    rms = np.sqrt(np.mean(residuals**2))
    outdata = np.vstack((data.lon, data.lat, data.lkv_E, data.lkv_N, data.lkv_U, data.LOS, model.LOS)).T
    outfilename = os.path.join(outdir, expname+'_data_vs_model_predictions.txt')
    np.savetxt(outfilename, outdata, header="lon, lat, lkv_E, lkv_N, lkv_U, data, model")
    outfilename = os.path.join(outdir, expname + '_fitted_parameters.txt')
    np.savetxt(outfilename, fitted_params, header="Params slip(cm), depth(km), plane, plane, reference")  # whole vector
    with open(os.path.join(outdir, expname+"_metrics.txt"), 'w') as outfile:
        outfile.write("RMS: %f mm\n" % rms)
        outfile.write("lam (tikhonov): %f\n" % lam)
        outfile.write("gamma (smoothing): %f\n" % gamma)
    return


def plot_complete_residual_vector_and_results(full_residuals, data, model_pred, param_vector, faults, outdir):
    slip_values = param_vector[0:39]
    depth_values = param_vector[39:39 * 2]
    normalized_data_resid = full_residuals[0:len(data.LOS)]  # the data part divided by L(Cd)
    data_residuals = data.LOS - model_pred.LOS  # the pure data residuals, in mm
    laplacian_slip = full_residuals[len(data.LOS):len(data.LOS) + 39]
    laplacian_depth = full_residuals[len(data.LOS) + 39:len(data.LOS) + 2 * 39]
    tikhonov_slip = full_residuals[len(data.LOS) + 2 * 39:len(data.LOS) + 3 * 39]
    tikhonov_depth = full_residuals[len(data.LOS) + 3 * 39:-1]
    data_misfit = np.sqrt(np.mean(data_residuals ** 2))
    chi2 = np.sqrt(np.mean(normalized_data_resid ** 2))
    lslip = np.sqrt(np.mean(laplacian_slip ** 2))
    ldepth = np.sqrt(np.mean(laplacian_depth ** 2))
    tslip = np.sqrt(np.mean(tikhonov_slip ** 2))
    tdepth = np.sqrt(np.mean(tikhonov_depth ** 2))
    total_misfit = np.sqrt(np.mean(full_residuals ** 2))

    print(f"  Data misfit: {data_misfit} mm")
    print(f"  Normalized data misfit: {chi2}")
    print(f"  Laplacian part: {np.sqrt(np.mean(np.concatenate((laplacian_slip, laplacian_depth)) ** 2))}")
    print(f"  Tikhonov part: {np.sqrt(np.mean(np.concatenate((tikhonov_slip, tikhonov_depth)) ** 2))}")
    print(f"  Total Loss Function RMS: {total_misfit}")

    plt.figure(figsize=(5, 5), dpi=300)
    plt.scatter(data.lon, data.lat, c=abs(normalized_data_resid), s=20, vmin=0, vmax=3.1)
    plt.colorbar()
    for item in faults:
        lons, lats = item.get_updip_corners_lon_lat()
        plt.plot(lons, lats, linewidth=0.5, color='firebrick')
    plt.title(f"RMS Data Misfit {data_misfit:.3f} mm and Chi2 {chi2:.3f}")
    plt.savefig(os.path.join(outdir, 'data_misfit_contributions.png'))

    # --- Figure & layout ---
    fig = plt.figure(figsize=(10, 7), constrained_layout=False, dpi=300)
    gs = fig.add_gridspec(
        nrows=2, ncols=2,
        width_ratios=[1, 2.2],  # left (square) : right (wide)
        height_ratios=[1, 1],
        wspace=0.3, hspace=0.35
    )

    ax00 = fig.add_subplot(gs[0, 0])  # top-left (square)
    ax01 = fig.add_subplot(gs[0, 1])  # top-right (rectangular)
    ax10 = fig.add_subplot(gs[1, 0])  # bottom-left (square)
    ax11 = fig.add_subplot(gs[1, 1])  # bottom-right (rectangular)

    # --- Left: square plots with colorbars ---
    im0 = ax00.scatter(data.lon, data.lat, c=abs(data_residuals), s=17, vmin=-3.1, vmax=3.1)
    ax00.set_title(f'g(m) - d: {data_misfit:.3f} mm rms')
    ax00.set_aspect('equal', adjustable='box')
    cbar0 = fig.colorbar(im0, ax=ax00, fraction=0.046, pad=0.04)
    cbar0.set_label('Value')

    im1 = ax10.scatter(data.lon, data.lat, c=abs(normalized_data_resid), s=17, vmin=0, vmax=3.1)
    ax10.set_title(f'L^-1 (g(m)-d): {chi2:.3f} rms')
    ax10.set_aspect('equal', adjustable='box')
    cbar1 = fig.colorbar(im1, ax=ax10, fraction=0.046, pad=0.04)
    cbar1.set_label('Value')

    # --- Right: long rectangular plots ---
    ax01.plot(slip_values, lw=2)
    ax01.set_title(f'Slip Values: Laplacian = {lslip:.3f}, Tikhonov = {tslip:.3f}')
    ax01.set_xlabel('Fault segment')
    ax01.set_ylim([0, 3.0])
    ax01.set_ylabel('Slip (cm)')

    ax11.plot(depth_values, lw=2)
    ax11.set_title(f'Depth Values: Laplacian = {ldepth:.3f}, Tikhonov = {tdepth:.3f}')
    ax11.set_xlabel('Fault segment')
    ax11.set_ylim([0, 5.0])
    ax11.set_ylabel('Depth (km)')

    # --- Overall title ---
    fig.suptitle(f'Total Loss Function Phi: {total_misfit:.3f}', fontsize=14, y=0.98)

    # Leave space for the suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(os.path.join(outdir, 'all_loss_contributions.png'))

    return total_misfit, data_misfit
