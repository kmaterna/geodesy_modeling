
import numpy as np
import collections
from Tectonic_Utilities.Tectonic_Utils.geodesy import euler_pole, haversine
import Elastic_stresses_py.PyCoulomb.coulomb_collections as cc
import Elastic_stresses_py.PyCoulomb.fault_slip_object.disp_points_object as disp_points
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library

"""
GF_element is everything you would need to make a column of the Green's matrix and plot the impulse response function. 
"""
GF_element = collections.namedtuple('GF_element', ['disp_points', 'fault_name',
                                                   'fault_dict_list', 'upper_bound', 'lower_bound',
                                                   'slip_penalty_flag']);


def pair_obs_gf(obs_disp_points, model_disp_points):
    """
    Operates on two disp_points objects, just pairing the 3-component objects together
    """
    paired_gf, paired_obs = [], [];
    for obs_item in obs_disp_points:
        for gf_item in model_disp_points:
            if abs(obs_item.lon - gf_item.lon) < 0.001 and abs(obs_item.lat - gf_item.lat) < 0.001:
                paired_gf.append(gf_item);
                paired_obs.append(obs_item);
                break;
    return paired_obs, paired_gf;


def pair_gf_elements_with_obs(obs_disp_points, gf_elements):
    """
    Take list of GF_elements, and list of obs_disp_points. Pare them down to a matching set of points in same order.
    The assumption is that all gf_elements have same points inside them (because we take first one as representative)
    Returns:
        paired_obs (list of disp_points)
        paired_gf_elements (list of gf_elements)
    """
    paired_gf_elements = [];  # a list of GF_element objects
    paired_obs, _ = pair_obs_gf(obs_disp_points, gf_elements[0].disp_points);  # get paired obs disp_points
    target_len = len(paired_obs);
    for gf_model in gf_elements:
        _, paired_gf = pair_obs_gf(obs_disp_points, gf_model.disp_points);  # one fault or CSZ patch
        paired_gf_elements.append(GF_element(disp_points=paired_gf, fault_name=gf_model.fault_name,
                                             fault_dict_list=gf_model.fault_dict_list, lower_bound=gf_model.lower_bound,
                                             upper_bound=gf_model.upper_bound,
                                             slip_penalty_flag=gf_model.slip_penalty_flag));
        if len(paired_gf) != target_len:
            raise ValueError("ERROR! Not all points have green's functions.");
    return paired_obs, paired_gf_elements;


def add_all_csz_patches(one_patch_gfs):
    """take n patches of the subduction interface and add their green's functions together """
    new_pts = one_patch_gfs[0];
    for i in range(1, len(one_patch_gfs)):
        new_pts = disp_points.add_disp_points(new_pts, one_patch_gfs[i]);
    return new_pts;


def get_GF_rotation_elements(obs_disp_points):
    """
    Build 3 GF_elements for horizontal rotation of GNSS velocities due to reference frames
    X rotation: [0, 0, 1] Euler Pole
    Y rotation: [90, 0, 1] Euler Pole
    Z rotation: [0, 89.99, 1] Euler Pole
    Returns list of GF_elements with theoretical displacements in all 3 directions.
    """
    x_disp_p, y_disp_p, z_disp_p = [], [], [];
    for obs_item in obs_disp_points:
        coords = [obs_item.lon, obs_item.lat];
        response_to_rot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 0, 1]);   # X direction
        response = cc.Displacement_points(lon=obs_item.lon, lat=obs_item.lat, dE_obs=response_to_rot[0],
                                          dN_obs=response_to_rot[1], dU_obs=response_to_rot[2],
                                          Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, meas_type=obs_item.meas_type);
        x_disp_p.append(response);
        response_to_rot = euler_pole.point_rotation_by_Euler_Pole(coords, [90, 0, 1]);  # Y direction
        response = cc.Displacement_points(lon=obs_item.lon, lat=obs_item.lat, dE_obs=response_to_rot[0],
                                          dN_obs=response_to_rot[1], dU_obs=response_to_rot[2],
                                          Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, meas_type=obs_item.meas_type);
        y_disp_p.append(response);
        response_to_rot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 89.99, 1]);  # Z direction
        response = cc.Displacement_points(lon=obs_item.lon, lat=obs_item.lat, dE_obs=response_to_rot[0],
                                          dN_obs=response_to_rot[1], dU_obs=response_to_rot[2],
                                          Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, meas_type=obs_item.meas_type);
        z_disp_p.append(response);

    xresponse = GF_element(disp_points=x_disp_p, fault_name='x_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1,
                           slip_penalty_flag=0);
    yresponse = GF_element(disp_points=y_disp_p, fault_name='y_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1,
                           slip_penalty_flag=0);
    zresponse = GF_element(disp_points=z_disp_p, fault_name='z_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1,
                           slip_penalty_flag=0);
    return [xresponse, yresponse, zresponse];


def get_GF_leveling_offset_element(obs_disp_points):
    """
    Build a GF_element for a reference frame leveling offset column of the GF matrix
    Input: a list of disp_points
    Output: a list of 1 GF_element, or an empty list if there is no leveling in this dataset
    """
    total_response_pts = [];
    lev_count = 0;
    for item in obs_disp_points:
        if item.meas_type == "leveling":
            response = cc.Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=1,
                                              Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, meas_type=item.meas_type);
            lev_count += 1;
        else:
            response = cc.Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=0,
                                              Se_obs=np.nan, Sn_obs=np.nan, Su_obs=np.nan, meas_type=item.meas_type);
        total_response_pts.append(response);
    lev_offset_gf = GF_element(disp_points=total_response_pts, fault_name='lev_offset', fault_dict_list=[],
                               upper_bound=1, lower_bound=-1, slip_penalty_flag=0);
    if lev_count == 0:
        return [];
    else:
        return [lev_offset_gf];


def get_displacement_directions(obs_disp_point, model_point):
    """
    Code up the logic for which components we model for each GNSS/leveling/tidegage point, etc
    """
    if obs_disp_point.meas_type == "continuous":
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs]);
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Su_obs]);
    elif obs_disp_point.meas_type == "survey":
        disps = np.array([model_point.dE_obs, model_point.dN_obs]);
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs]);
    elif obs_disp_point.meas_type == "leveling":
        disps = np.array([model_point.dU_obs]);
        sigmas = np.array([model_point.Su_obs]);
    elif obs_disp_point.meas_type == "tide_gage":
        disps = np.array([model_point.dU_obs]);
        sigmas = np.array([model_point.Su_obs]);
    else:
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs]);
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Su_obs]);
    return disps, sigmas;


def buildG_column(GF_disp_points, obs_disp_points):
    """
    Green's functions for a single model element at each observation point (i.e., slip on one fault).
    Returns a single column of nx1.
    """
    GF_col = [];
    for g_item, obs in zip(GF_disp_points, obs_disp_points):
        new_values, _ = get_displacement_directions(obs, g_item);
        GF_col = np.concatenate((GF_col, new_values));
    GF_col = np.reshape(GF_col, (len(GF_col), 1));
    return GF_col;


def build_obs_vector(obs_disp_points):
    """
    Build observation vector (will eventually have leveling here too)
    """
    obs, sigmas = [], [];
    for item in obs_disp_points:
        new_values, new_sigmas = get_displacement_directions(item, item);
        obs = np.concatenate((obs, new_values), axis=0);
        sigmas = np.concatenate((sigmas, new_sigmas), axis=0);
    return obs, sigmas;


def unpack_model_pred_vector(model_pred, paired_obs):
    """
    Unpack a model vector into a bunch of disp_point objects. Same logic implemented here as in the functions above.
    """
    disp_points_list = [];
    counter = 0;
    for i in range(len(paired_obs)):
        if paired_obs[i].meas_type == "survey":
            E = model_pred[counter];
            N = model_pred[counter+1];
            U = np.nan;
            counter = counter+2;
        elif paired_obs[i].meas_type == "leveling":
            E, N = np.nan, np.nan;
            U = model_pred[counter];
            counter = counter+1;
        elif paired_obs[i].meas_type == "tide_gage":
            E, N = np.nan, np.nan;
            U = model_pred[counter];
            counter = counter+1;
        else:
            E = model_pred[counter];
            N = model_pred[counter+1];
            U = model_pred[counter+2];
            counter = counter+3;

        new_disp_point = cc.Displacement_points(lon=paired_obs[i].lon, lat=paired_obs[i].lat,
                                                dE_obs=E, dN_obs=N, dU_obs=U,
                                                Se_obs=0, Sn_obs=0, Su_obs=0,
                                                name=paired_obs[i].name, meas_type=paired_obs[i].meas_type);
        disp_points_list.append(new_disp_point);
    return disp_points_list;


def unpack_model_of_rotation_only(M_vector, parameter_names):
    """
    1) Zeros out all model parameters except the rotation (leaving the last three)
    2) Zeros out all the rotation (canceling the last three)
    :parameter M_vector: vector of model parameters
    :parameter parameter_names: list of strings, derived from fault_name of a GF_element
    """
    rot_multiplier = np.zeros(np.shape(M_vector));
    fault_multiplier = np.zeros(np.shape(M_vector));
    for i in range(len(parameter_names)):
        if parameter_names[i] == 'lev_offset':
            continue;
        elif parameter_names[i] in ["x_rot", "y_rot", "z_rot"]:
            rot_multiplier[i] = 1;
        else:
            fault_multiplier[i] = 1;
    M_rot_only = np.multiply(M_vector, rot_multiplier);
    M_no_rot = np.multiply(M_vector, fault_multiplier);
    return M_rot_only, M_no_rot;


def unpack_model_of_particular_fault(M_vector, parameter_names, target_fault):
    """
    Simplify a model vector into only those components from particular target_fault (string), such as 'CSZ_dist'
    Useful for investigating a forward model.
    """
    fault_params_vector = np.zeros(np.shape(M_vector));
    for i in range(len(parameter_names)):
        if parameter_names[i] == target_fault:
            fault_params_vector[i] = 1;
    M_fault = np.multiply(M_vector, fault_params_vector);
    return M_fault;

def get_fault_element_distance(fault_dict1, fault_dict2):
    distance = haversine.distance([fault_dict1["lat"], fault_dict1["lon"]], [fault_dict2["lat"], fault_dict2["lon"]]);
    return distance;


def build_smoothing(gf_elements, fault_name, strength, G, obs, sigmas):
    """
    Make a weighted connectivity matrix that has the same number of columns as G, that can be appended to the bottom.
    Any element within gf_element that has fault_name will have its immediate neighbors subtracted for smoothing.
    Append the matching number of zeros to the obs_vector.
    Assumes similar-sized patches throughout the slip distribution

    :param gf_elements: list of gf_element objects
    :type gf_elements: list
    :param fault_name: which fault elements are we smoothing, string
    :param strength: lambda parameter in smoothing equation
    :param G: already existing G matrix
    :param obs: already existing obs vector
    :param sigmas: already existing sigma vector
    """
    print("G and obs before smoothing:", np.shape(G), np.shape(obs));
    if strength == 0:
        print("No change, smoothing set to 0");
        return G, obs, sigmas;   # returning unaltered G if there is no smoothing

    G_smoothing = np.zeros((len(gf_elements), len(gf_elements)));

    # Get critical distance, a typical small distance between neighboring fault patches.
    # Operates on the first patch that it finds.
    distances = [];
    for i in range(len(gf_elements)):
        if gf_elements[i].fault_name == fault_name:
            for j in range(i+1, len(gf_elements)):
                if gf_elements[j].fault_name == fault_name:
                    distances.append(get_fault_element_distance(gf_elements[i].fault_dict_list[0],
                                                                gf_elements[j].fault_dict_list[0]));
            break;
    critical_distance = sorted(distances)[2] + 5;  # take adjacent patches with some wiggle room

    # Build the parts of the matrix for smoothing
    for i in range(len(gf_elements)):
        if gf_elements[i].fault_name == fault_name:
            G_smoothing[i][i] = 1;
            for j in range(len(gf_elements)):
                if gf_elements[j].fault_name == fault_name:
                    if i != j and get_fault_element_distance(gf_elements[i].fault_dict_list[0],
                                                             gf_elements[j].fault_dict_list[0]) < critical_distance:
                        G_smoothing[i][j] = -1/4;

    G_smoothing = G_smoothing * strength;  # multiplying by lambda factor

    # observation vector of zeros
    zero_vector = np.zeros((len(gf_elements),));

    G_smoothing = np.vstack((G, G_smoothing));    # appending smoothing matrix
    smoothed_obs = np.concatenate((obs, zero_vector));   # appending smoothing components to data
    smoothed_sigmas = np.concatenate((sigmas, zero_vector));  # appending smoothing components to sigmas
    print("G and obs after smoothing:", np.shape(G_smoothing), np.shape(smoothed_obs));
    return G_smoothing, smoothed_obs, smoothed_sigmas;

def build_slip_penalty(gf_elements, penalty, G, obs, sigmas):
    """
    Minimum-norm smoothing constraint.
    Build a square diagonal matrix to go at the bottom of G, with 1's along the elements that will be slip-penalized.
    Zeros along obs and sigmas.
    Overwrites old G, obs, and sigma variables.
    """
    print("G and obs before slip penalty:", np.shape(G), np.shape(obs));
    if penalty == 0:
        print("No change, slip penalty set to 0");
        return G, obs, sigmas;   # returning unaltered G if there is no smoothing

    G_penalty = np.zeros((len(gf_elements), len(gf_elements)));
    for i in range(len(gf_elements)):
        if gf_elements[i].slip_penalty_flag == 1:
            G_penalty[i][i] = 1;

    G_penalty = G_penalty * penalty;  # multiplying by lambda factor

    # observation vector of zeros
    zero_vector = np.zeros((len(gf_elements),));

    G_penalty = np.vstack((G, G_penalty));    # appending smoothing matrix
    smoothed_obs = np.concatenate((obs, zero_vector));   # appending smoothing components to data
    smoothed_sigmas = np.concatenate((sigmas, zero_vector));  # appending smoothing components to sigmas
    print("G and obs after penalty:", np.shape(G_penalty), np.shape(smoothed_obs));
    return G_penalty, smoothed_obs, smoothed_sigmas;


def filter_out_smoothing_lines(pred_vector, obs_vector, sigma_vector):
    """
    Remove lines of zeros automatically added to obs/sigma vectors for smoothing and slip penalty parameters.
    All inputs are expected to be vectors with the same length.
    """
    tol = 1e-9;
    new_pred_vector, new_obs_vector, new_sigma_vector = [], [], [];
    for i in range(len(pred_vector)):
        if abs(pred_vector[i]) < tol and abs(obs_vector[i]) < tol and abs(sigma_vector[i]) < tol:
            continue;
        else:
            new_pred_vector.append(pred_vector[i]);
            new_obs_vector.append(obs_vector[i]);
            new_sigma_vector.append(sigma_vector[i]);
    print("During RMS calc., filter vectors from %d to %d for "
          "smoothing and penalty" % (len(pred_vector), len(new_pred_vector)) );
    return new_pred_vector, new_obs_vector, new_sigma_vector;


def rms_from_model_pred_vector(pred_vector, obs_vector, sigma_vector):
    """
    Various ways to compute the RMS value
    All inputs are expected to be vectors with the same length.
    """
    # Filter out lines added for smoothing
    pred_vector, obs_vector, sigma_vector = filter_out_smoothing_lines(pred_vector, obs_vector, sigma_vector);
    residuals = np.subtract(obs_vector, pred_vector);
    rms_mm = np.multiply(np.sqrt(np.mean(np.square(residuals))), 1000);
    added_sum = 0;
    for i in range(len(pred_vector)):
        chi2 = np.square(obs_vector[i]-pred_vector[i]) / np.square(sigma_vector[i]);
        added_sum = added_sum + chi2
    reported_chi2 = added_sum / len(sigma_vector);
    return rms_mm, reported_chi2;


def write_model_params(v, residual, outfile, GF_elements=None):
    """
    :param v: vector of model parameters, floats
    :param residual: float, mm/yr
    :param outfile: string
    :param GF_elements: optional, list of GF_element objects
    """
    print("Writing %s" % outfile);
    ofile = open(outfile, 'w');
    for item in v:
        ofile.write(str(item)+"\n");
    report_string = "\nRMS: %f mm/yr\n" % residual;
    ofile.write(report_string);
    if GF_elements:
        for item in GF_elements:
            ofile.write(item.fault_name+" ");
    ofile.close();
    return;

def write_summary_params(v, residual, outfile, GF_elements, ignore_faults=(), message=''):
    """
    Write a human-readable results file, with the potential to ignore faults with distributed models for clarity
    :param v: vector of model parameters, floats
    :param residual: float, mm/yr
    :param outfile: string
    :param GF_elements: list of GF_element objects
    :param ignore_faults: list of strings
    :param message: optional message from inverse routine about solution quality
    """
    print("Writing %s" % outfile);
    ofile = open(outfile, 'w');
    for i in range(len(v)):
        if GF_elements[i].fault_name in ignore_faults:
            continue;
        if GF_elements[i].fault_name in ["x_rot", "y_rot", "z_rot"]:
            unit = "deg/Ma"
        else:
            unit = "cm/yr"
        ofile.write(GF_elements[i].fault_name + ": ");
        ofile.write(str(v[i])+" "+unit+"\n");
    report_string = "\nRMS: %f mm/yr, on %d observations\n" % (residual, len(GF_elements[0].disp_points));
    ofile.write(report_string);
    ofile.write("Message: "+message+"\n");
    ofile.close();
    return;


def view_full_results(exp_dict, paired_obs, modeled_disp_points, residual_pts, rotation_pts, norot_pts, title, region,
                      special_pts=()):
    # Plot the data, model, and residual in separate plots
    # Not plotting the fault patches because it takes a long time.
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"] + "/data_only.png",
                                                         disp_points=paired_obs, region=region,
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/residuals.png",
                                                         disp_points=residual_pts, region=region,
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001,
                                                         title=title);
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/model_pred.png",
                                                         disp_points=modeled_disp_points, region=region,
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/rotation_pred.png",
                                                         disp_points=rotation_pts, region=region,
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/faults_only_pred.png",
                                                         disp_points=norot_pts, region=region,
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001)
    if special_pts:
        library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"] + "/csz_only_pred.png",
                                                             disp_points=special_pts, region=region,
                                                             scale_arrow=(1.0, 0.010, "1 cm/yr"),
                                                             v_labeling_interval=0.001)
    return;


def visualize_GF_elements(GF_elements_list, outdir, exclude_list=()):
    """
    Aside to the main calculation, just view each GF.
    Inputs: list of GF_elements objects
    string for outdir
    """
    for GF_element in GF_elements_list:
        if GF_element.fault_name in exclude_list:   # plot these elements separately, like individual CSZ patches
            continue;
        print(GF_element.fault_name);
        if GF_element.fault_name == "CSZ":
            scale_arrow = (1.0, 0.010, "1 cm");
        else:
            scale_arrow = (1.0, 0.001, "1 mm");
        library.plot_fault_slip.map_source_slip_distribution(GF_element.fault_dict_list, outdir + "/gf_" +
                                                             GF_element.fault_name + "_only.png",
                                                             disp_points=GF_element.disp_points,
                                                             region=[-127, -119.7, 37.7, 43.3],
                                                             scale_arrow=scale_arrow,
                                                             v_labeling_interval=0.001);
    return;


def view_csz(csz_patches, one_patch_gfs):
    """
    View different parts of the distributed CSZ Green's functions.
    Currently not called. Probably.
    """
    collective_csz_gf = add_all_csz_patches(one_patch_gfs);
    library.plot_fault_slip.map_source_slip_distribution(csz_patches, "CSZ_patches.png", disp_points=collective_csz_gf,
                                                         region=[-127, -119.7, 37.7, 43.3],
                                                         scale_arrow=(1.0, 0.010, "1 cm"), v_labeling_interval=0.001);
    idx = 0;
    patch = csz_patches[idx];
    one_patch_gf = one_patch_gfs[idx];
    library.plot_fault_slip.map_source_slip_distribution([patch], "CSZ_patch.png", disp_points=one_patch_gf,
                                                         region=[-127, -119.7, 37.7, 43.3],
                                                         scale_arrow=(1.0, 0.001, "1 mm"), v_labeling_interval=0.001);
    return;
