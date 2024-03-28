
import numpy as np
from Elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library
import Tectonic_Utils.seismo.moment_calculations as mo
from .GF_element.GF_element import GF_element
import os


def pair_obs_model(obs_disp_pts, model_disp_pts, tol=0.001):
    """
    Filters two lists of disp_points objects, just pairing the objects together where their locations match

    :param obs_disp_pts: list of disp_points objects, length n
    :param model_disp_pts: list of disp_points objects, length m
    :param tol: tolerance, default 0.001 degrees
    :returns: list of disp_point_obj, list of disp_point_obj, with matching lengths p
    """
    paired_model, paired_obs = [], []
    for obs_item in obs_disp_pts:
        for model_item in model_disp_pts:
            if abs(obs_item.lon - model_item.lon) < tol and abs(obs_item.lat - model_item.lat) < tol:
                paired_model.append(model_item)
                paired_obs.append(obs_item)
                break
    return paired_obs, paired_model


def pair_gf_elements_with_obs(obs_disp_points, gf_elements, tol=0.001):
    """
    Take list of GF_elements, and list of obs_disp_points. Pare them down to a matching set of points in same order.
    The assumption is that all gf_elements have same points inside them (because we take first one as representative)

    :param obs_disp_points: list of disp_points
    :param gf_elements: a list of gf_elements with all the same points inside them
    :param tol: tolerance for pairing station with station, in degrees
    :returns: paired_obs (list of disp_points), paired_gf_elements (list of gf_elements)
    """
    paired_gf_elements = []  # a list of GF_element objects
    paired_obs, _ = pair_obs_model(obs_disp_points, gf_elements[0].disp_points, tol=tol)  # get paired obs disp_points
    target_len = len(paired_obs)
    for gf_model in gf_elements:
        _, paired_gf = pair_obs_model(obs_disp_points, gf_model.disp_points, tol=tol)  # one fault or CSZ patch
        paired_gf_elements.append(GF_element(disp_points=paired_gf, param_name=gf_model.param_name,
                                             fault_dict_list=gf_model.fault_dict_list, lower_bound=gf_model.lower_bound,
                                             upper_bound=gf_model.upper_bound,
                                             slip_penalty=gf_model.slip_penalty, units=gf_model.units,
                                             points=gf_model.points))
        if len(paired_gf) != target_len:
            raise ValueError("ERROR! Not all points have green's functions.")
    return paired_obs, paired_gf_elements


def get_displacement_directions(obs_disp_point, model_point):
    """
    Code up the logic for which components we model for each GNSS/leveling/tidegage/insar point, etc
    """
    if obs_disp_point.meas_type == "continuous":
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs])
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Su_obs])
    elif obs_disp_point.meas_type == "survey":
        disps = np.array([model_point.dE_obs, model_point.dN_obs])
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs])
    elif obs_disp_point.meas_type == "leveling":
        disps = np.array([model_point.dU_obs])
        sigmas = np.array([model_point.Su_obs])
    elif obs_disp_point.meas_type == "tide_gage":
        disps = np.array([model_point.dU_obs])
        sigmas = np.array([model_point.Su_obs])
    elif obs_disp_point.meas_type == "insar":
        disps = np.array([model_point.dE_obs])
        sigmas = np.array([model_point.Se_obs])
    else:
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs])
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Su_obs])
    return disps, sigmas


def unpack_model_pred_vector(model_pred, paired_obs):
    """
    Unpack a model vector into a bunch of disp_point objects. Same logic implemented here as in the functions above.
    In the future, this can be done with a dictionary or look-up table passed into the function.

    :param model_pred: long vector of model parameters, corresponding to each component being used
    :param paired_obs: list of disp_point_objects (shorter than model_pred vector)
    """
    disp_points_list = []
    counter = 0
    for i in range(len(paired_obs)):
        if paired_obs[i].meas_type == "survey":
            E = model_pred[counter]
            N = model_pred[counter+1]
            U = np.nan
            counter = counter+2
        elif paired_obs[i].meas_type == "leveling":
            E, N = np.nan, np.nan
            U = model_pred[counter]
            counter = counter+1
        elif paired_obs[i].meas_type == "tide_gage":
            E, N = np.nan, np.nan
            U = model_pred[counter]
            counter = counter+1
        elif paired_obs[i].meas_type == "insar":
            E = model_pred[counter]
            N, U = np.nan, np.nan
            counter = counter+1
        else:
            E = model_pred[counter]
            N = model_pred[counter+1]
            U = model_pred[counter+2]
            counter = counter+3

        new_disp_point = Displacement_points(lon=paired_obs[i].lon, lat=paired_obs[i].lat, dE_obs=E, dN_obs=N, dU_obs=U,
                                             Se_obs=0, Sn_obs=0, Su_obs=0, name=paired_obs[i].name,
                                             meas_type=paired_obs[i].meas_type, refframe=paired_obs[i].refframe)
        disp_points_list.append(new_disp_point)
    return disp_points_list


def buildG_column(GF_disp_points, obs_disp_points):
    """
    Green's functions for a single model element at each observation point (i.e., slip on one fault).
    Returns a single column of nx1.
    """
    GF_col = []
    if len(GF_disp_points) != len(obs_disp_points):
        raise ValueError("Error! Length of modeled and observe vectors do not agree.")
    for g_item, obs in zip(GF_disp_points, obs_disp_points):
        new_values, _ = get_displacement_directions(obs, g_item)
        GF_col = np.concatenate((GF_col, new_values))
    GF_col = np.reshape(GF_col, (len(GF_col), 1))
    return GF_col


def build_obs_vector(obs_disp_points):
    """
    Build observation 1D-vector.
    """
    obs, sigmas = [], []
    for item in obs_disp_points:
        new_values, new_sigmas = get_displacement_directions(item, item)
        obs = np.concatenate((obs, new_values), axis=0)
        sigmas = np.concatenate((sigmas, new_sigmas), axis=0)
    return obs, sigmas


def forward_disp_points_predictions(G, m, sigmas, paired_obs):
    """Create a convenient list of disp_points from a forward prediction based on G and m and sigma matrices/vectors."""
    model_pred = G.dot(m) * sigmas
    model_disp_points = unpack_model_pred_vector(model_pred, paired_obs)
    return model_disp_points


def unpack_model_of_target_param(M_vector, parameter_names, target_names=()):
    """
    Simplify a model vector into only those components from particular target_name (string), such as 'CSZ_dist'
    Useful for investigating a forward model.

    :param M_vector: 1D array of numbers, model parameter values
    :param parameter_names: list of parameters names for each model parameter
    :param target_names: list, names of desired model parameter
    :returns: vector of zeros except for the model values corresponding to target name (e.g., 5 mm/yr on one fault)
    """
    target_params_vector = np.zeros(np.shape(M_vector))
    for i in range(len(parameter_names)):
        if parameter_names[i] in target_names:
            target_params_vector[i] = 1
    M_target = np.multiply(M_vector, target_params_vector)
    return M_target


def unpack_model_without_target_param(M_vector, parameter_names, exclude_names=()):
    """
    Simplify a model vector into only those components excluding target_names (string)
    Useful for investigating a forward model.

    :param M_vector: 1D array of numbers, model parameter values
    :param parameter_names: list of parameters names for each model parameter
    :param exclude_names: list, names of undesired model parameter
    :returns: vector of parameter values except for the model values corresponding to exclude name
    """
    target_params_vector = np.ones(np.shape(M_vector))
    for i in range(len(parameter_names)):
        if parameter_names[i] in exclude_names:
            target_params_vector[i] = 0
    M_target = np.multiply(M_vector, target_params_vector)
    return M_target


def build_smoothing(gf_elements, param_name_list, strength, lengthscale, G, obs, sigmas, distance_3d=True,
                    laplacian_operator=-1/4):
    """
    Make a weighted connectivity matrix that has the same number of columns as G, that can be appended to the bottom.
    Any gf_element that has param_name will have its immediate neighbors subtracted for smoothing.
    Append the matching number of zeros to the obs_vector.
    Assumes similar-sized patches throughout the slip distribution

    :param gf_elements: list of gf_element objects
    :type gf_elements: list
    :param param_name_list: which fault elements are we smoothing, tuple of strings
    :param strength: lambda parameter in smoothing equation
    :param lengthscale: distance over which fault elements are smooth (i.e., correlated)
    :param G: already existing G matrix
    :param obs: already existing obs vector
    :param sigmas: already existing sigma vector
    :param distance_3d: bool, do you compute distance between fault patches in 3d way, YES or NO?
    :param laplacian_operator: how strong do you smooth the neighbor? -1/4 (compared to 1 for base element) is default.
    """
    print("G and obs before smoothing:", np.shape(G), np.shape(obs))
    if strength == 0:
        print("No change, smoothing set to 0")
        return G, obs, sigmas   # returning unaltered G if there is no smoothing

    G_smoothing = np.zeros((len(gf_elements), len(gf_elements)))

    # Get critical distance, a typical small distance between neighboring fault patches.
    # Operates on the first patch that it finds.
    distances = []
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name in param_name_list:
            for j in range(i+1, len(gf_elements)):
                if gf_elements[j].param_name in param_name_list:
                    distances.append(gf_elements[i].fault_dict_list[0].get_fault_element_distance(
                        gf_elements[j].fault_dict_list[0]))
            break
    critical_distance = sorted(distances)[2] + lengthscale  # smooth adjacent patches with some wiggle room

    # Build the parts of the matrix for smoothing
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name in param_name_list:
            G_smoothing[i][i] = 1
            for j in range(len(gf_elements)):
                if gf_elements[j].param_name in param_name_list:
                    if i != j and gf_elements[i].fault_dict_list[0].get_fault_element_distance(
                            gf_elements[j].fault_dict_list[0], threedimensional=distance_3d) < critical_distance:
                        G_smoothing[i][j] = laplacian_operator

    G_smoothing = G_smoothing * strength  # multiplying by lambda factor

    # observation vector of zeros
    zero_vector = np.zeros((len(gf_elements),))

    G_smoothing = np.vstack((G, G_smoothing))    # appending smoothing matrix
    smoothed_obs = np.concatenate((obs, zero_vector))   # appending smoothing components to data
    smoothed_sigmas = np.concatenate((sigmas, zero_vector))  # appending smoothing components to sigmas
    print("G and obs after smoothing:", np.shape(G_smoothing), np.shape(smoothed_obs))
    return G_smoothing, smoothed_obs, smoothed_sigmas


def build_slip_penalty(gf_elements, penalty, G, obs, sigmas):
    """
    Minimum-norm smoothing constraint.
    Build a square diagonal matrix to go at the bottom of G, with 1's along the elements that will be slip-penalized.
    Zeros along obs and sigmas.
    Overwrites old G, obs, and sigma variables.
    """
    print("G and obs before slip penalty:", np.shape(G), np.shape(obs))
    if penalty == 0:
        print("No change, slip penalty set to 0")
        return G, obs, sigmas   # returning unaltered G if there is no smoothing

    G_penalty = np.zeros((len(gf_elements), len(gf_elements)))
    for i in range(len(gf_elements)):
        if gf_elements[i].slip_penalty > 0:
            G_penalty[i][i] = gf_elements[i].slip_penalty

    G_penalty = G_penalty * penalty  # multiplying by lambda factor

    # observation vector of zeros
    zero_vector = np.zeros((len(gf_elements),))

    G_penalty = np.vstack((G, G_penalty))    # appending smoothing matrix
    smoothed_obs = np.concatenate((obs, zero_vector))   # appending smoothing components to data
    smoothed_sigmas = np.concatenate((sigmas, zero_vector))  # appending smoothing components to sigmas
    print("G and obs after penalty:", np.shape(G_penalty), np.shape(smoothed_obs))
    return G_penalty, smoothed_obs, smoothed_sigmas


def filter_out_smoothing_lines(pred_vector, obs_vector, sigma_vector, tol=1e-9):
    """
    Remove lines of zeros automatically added to obs/sigma vectors for smoothing and slip penalty parameters.
    All inputs are expected to be vectors with the same length.

    :param pred_vector: array, in m
    :param obs_vector: array, in m
    :param sigma_vector: array, in m
    :param tol: optional, 1e-9
    """
    new_pred_vector, new_obs_vector, new_sigma_vector = [], [], []
    for i in range(len(pred_vector)):
        if abs(pred_vector[i]) < tol and abs(obs_vector[i]) < tol and abs(sigma_vector[i]) < tol:
            continue
        else:
            new_pred_vector.append(pred_vector[i])
            new_obs_vector.append(obs_vector[i])
            new_sigma_vector.append(sigma_vector[i])
    print("During RMS calc., filter vectors from %d to %d for "
          "smoothing and penalty" % (len(pred_vector), len(new_pred_vector)))
    return new_pred_vector, new_obs_vector, new_sigma_vector


def write_fault_traces(M_vector, paired_gf_elements, outfile, ignore_faults=()):
    """Write a gmt multi-segment file with fault traces (GF_element.points) and colors from values.
    Will only write when the GF_element has the field "points".

    :param M_vector: vector containing value of model parameters, floats (e.g., in mm/yr)
    :param paired_gf_elements: matching list of GF_elements
    :param outfile: filename, string
    :param ignore_faults: parameter names to ignore when printing gmt file
    """
    print("Writing %s" % outfile)
    with open(outfile, 'w') as ofile:
        for i in range(len(M_vector)):
            if paired_gf_elements[i].param_name in ignore_faults:
                continue
            if len(paired_gf_elements[i].points) == 0:
                continue
            ofile.write("> -Z%f \n" % M_vector[i])
            for coord in paired_gf_elements[i].points:
                ofile.write("%f %f\n" % (coord[0], coord[1]))
    return


def write_raw_model_params(outfile, v, rms_residual=0, GF_elements=None):
    """
    The most general way of reporting the model vector, not especially human-readable.

    :param outfile: string, filename
    :param v: vector of model parameters, floats
    :param rms_residual: optional, float, mm/yr
    :param GF_elements: optional, list of paired GF_element objects to write the parameter names
    """
    print("Writing raw model outputs in %s" % outfile)
    ofile = open(outfile, 'w')
    for item in v:
        ofile.write(str(item)+"\n")
    report_string = "\nRMS: %f mm/yr\n" % rms_residual
    ofile.write(report_string)
    if GF_elements:
        for item in GF_elements:
            ofile.write(item.param_name+" ")
    ofile.close()
    return


def write_summary_params(v, outfile, GF_elements, ignore_faults=(), message=''):
    """
    Write a human-readable results file, with the potential to ignore some faults for clarity.

    :param v: vector of model parameters, floats
    :param outfile: string
    :param GF_elements: list of GF_element objects
    :param ignore_faults: list of strings
    :param message: optional message from inverse routine about solution quality
    """
    print("Writing %s" % outfile)
    ofile = open(outfile, 'w')
    for i in range(len(v)):
        if GF_elements[i].param_name in ignore_faults:
            continue
        ofile.write(GF_elements[i].param_name + ": ")
        if GF_elements[i].units == "cm/yr":   # converting fault slip rates into mm/yr for convenience
            units = 'mm/yr'
            multiplier = 10
        else:
            units = GF_elements[i].units
            multiplier = 1
        ofile.write("%.5f %s" % (v[i]*multiplier, units))
        ofile.write('  [within %.3f to %.3f]' % (GF_elements[i].lower_bound*multiplier,
                                                 GF_elements[i].upper_bound*multiplier))
        ofile.write("\n")
    report_string = "\nWith %d observations\n" % (len(GF_elements[0].disp_points))
    ofile.write(report_string)
    ofile.write("Message: "+message+"\n")
    ofile.close()
    return


def write_standard_misfit_report(outdir, obs_disp_pts, model_disp_pts, total_moment):
    """
    Write a standard type of misfit report into a text file.

    :param outdir: string
    :param obs_disp_pts: list of disp_point_objects
    :param model_disp_pts: list of disp_point_objects
    :param total_moment: float, total moment of the slip
    """
    [_all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_disp_pts, model_disp_pts)
    with open(outdir+'/metrics.txt', 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm)
        print("total moment (N-m): ", total_moment)
        print("Equivalent to:", mo.mw_from_moment(total_moment))
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm)
        ofile.write("total moment (N-m): %f\n" % total_moment)
        ofile.write("Equivalent to: %f\n" % mo.mw_from_moment(total_moment))
    return


def view_full_results(exp_dict, paired_obs, modeled_disp_points, residual_pts, rotation_pts, norot_pts, title, region):
    # Plot the data, model, and residual in separate plots
    # Not plotting the fault patches because it takes a long time.
    scale_arrow = (1.0, 0.010, "1 cm/yr")
    library.plot_fault_slip.map_source_slip_distribution([], os.path.join(exp_dict["outdir"], "data_only.png"),
                                                         disp_points=paired_obs, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], os.path.join(exp_dict["outdir"], "residuals.png"),
                                                         disp_points=residual_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001,
                                                         title=title)
    library.plot_fault_slip.map_source_slip_distribution([], os.path.join(exp_dict["outdir"], "model_pred.png"),
                                                         disp_points=modeled_disp_points, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], os.path.join(exp_dict["outdir"], "rotation_pred.png"),
                                                         disp_points=rotation_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], os.path.join(exp_dict["outdir"], "faults_only_pred.png"),
                                                         disp_points=norot_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    return


def remove_nearfault_pts(obs_points, fault_trace_file):
    """
    Wraps dpo.utilities.filter_to_remove_near_fault so that we can pass a filename.

    :param obs_points: a list of disp_point objects
    :param fault_trace_file: filename that contains two columns of numbers for fault trace: lon, lat
    """
    trace_pts = np.loadtxt(fault_trace_file)
    print("Removing near-fault points from file %s " % fault_trace_file)
    obs_disp_points = dpo.utilities.filter_to_remove_near_fault(obs_points, trace_pts, radius_km=10)
    return obs_disp_points
