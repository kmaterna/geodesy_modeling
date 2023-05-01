
import numpy as np
from Tectonic_Utilities.Tectonic_Utils.geodesy import euler_pole
from Elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library


class GF_element:
    """
    GF_element is everything you would need to make a column of the Green's matrix and
    plot the impulse response function.
    'points' is coordinates of surface trace of the fault, if applicable.

    :param disp_points: modeled displacement_points due to unit activation of this GF_element
    :type disp_points: list of Displacement_points objects
    :param param_name: param_name
    :type param_name: string
    :param fault_dict_list: list of fault_slip_objects
    :type fault_dict_list: list
    :param upper_bound: highest allowed value of this GF_element
    :type upper_bound: float
    :param lower_bound: lowest allowed value of this GF_element
    :type lower_bound: float
    :param slip_penalty: a number that will be attached to smoothing in the G matrix
    :type slip_penalty: float
    :param units: what units is this 'unit activation' in?
    :type units: string
    :param points: coordinates of surface trace of fault, if provided
    :type points: np.array
    """

    def __init__(self, disp_points, param_name, upper_bound, lower_bound, slip_penalty, units,
                 fault_dict_list=(), points=()):
        self.disp_points = disp_points;
        self.param_name = param_name;
        self.fault_dict_list = fault_dict_list;
        self.upper_bound = upper_bound;
        self.lower_bound = lower_bound;
        self.slip_penalty = slip_penalty;
        self.units = units;
        self.points = points;  # coordinates of surface trace of fault, if applicable


def pair_obs_model(obs_disp_pts, model_disp_pts):
    """
    Filters two lists of disp_points objects, just pairing the objects together where their locations match

    :param obs_disp_pts: list of disp_points objects
    :param model_disp_pts: list of disp_points objects
    :returns: list of disp_point_obj, list of disp_point_obj
    """
    tol = 0.001;
    paired_model, paired_obs = [], [];
    for obs_item in obs_disp_pts:
        for model_item in model_disp_pts:
            if abs(obs_item.lon - model_item.lon) < tol and abs(obs_item.lat - model_item.lat) < tol:
                paired_model.append(model_item);
                paired_obs.append(obs_item);
                break;
    return paired_obs, paired_model;


def pair_gf_elements_with_obs(obs_disp_points, gf_elements):
    """
    Take list of GF_elements, and list of obs_disp_points. Pare them down to a matching set of points in same order.
    The assumption is that all gf_elements have same points inside them (because we take first one as representative)

    :param obs_disp_points: list of disp_points
    :param gf_elements: a list of gf_elements with all the same points inside them
    :returns: paired_obs (list of disp_points), paired_gf_elements (list of gf_elements)
    """
    paired_gf_elements = [];  # a list of GF_element objects
    paired_obs, _ = pair_obs_model(obs_disp_points, gf_elements[0].disp_points);  # get paired obs disp_points
    target_len = len(paired_obs);
    for gf_model in gf_elements:
        _, paired_gf = pair_obs_model(obs_disp_points, gf_model.disp_points);  # one fault or CSZ patch
        paired_gf_elements.append(GF_element(disp_points=paired_gf, param_name=gf_model.param_name,
                                             fault_dict_list=gf_model.fault_dict_list, lower_bound=gf_model.lower_bound,
                                             upper_bound=gf_model.upper_bound,
                                             slip_penalty=gf_model.slip_penalty, units=gf_model.units,
                                             points=gf_model.points));
        if len(paired_gf) != target_len:
            raise ValueError("ERROR! Not all points have green's functions.");
    return paired_obs, paired_gf_elements;


def add_gfs(single_gfs):
    """
    Take some gf_elements and add their green's functions together.

    :param single_gfs: list of gf_objects with the same modeled_disp_points lists
    :returns: list of disp_points
    """
    new_pts = single_gfs[0];
    for i in range(1, len(single_gfs)):
        new_pts = dpo.utilities.add_disp_points(new_pts, single_gfs[i]);
    return new_pts;


def get_GF_rotation_element(obs_disp_points, ep, target_region=(-180, 180, -90, 90), rot_name=''):
    """
    Build one GF_elements for a horizontal rotation of GNSS velocities about some axis

    :param obs_disp_points: list of disp_points
    :param ep: [lon, lat, rate] of euler pole for desired rotation
    :param target_region: list of lon/lon/lat/lat for bounding box
    :param rot_name: string, optional metadata for naming the rotation (ex: ocb_)
    :returns: one GF_element with rotation displacements in x, y, and z directions
    """
    rot_disp_p = [];
    for obs_item in obs_disp_points:
        coords = [obs_item.lon, obs_item.lat];
        if obs_item.is_within_bbox(target_region):
            mult = 1;
        else:
            mult = 0;
        response_to_rot = euler_pole.point_rotation_by_Euler_Pole(coords, ep);
        response = Displacement_points(lon=obs_item.lon, lat=obs_item.lat, dE_obs=mult*response_to_rot[0],
                                       dN_obs=mult*response_to_rot[1], dU_obs=mult*response_to_rot[2], Se_obs=0,
                                       Sn_obs=0, Su_obs=0, meas_type=obs_item.meas_type, refframe=obs_item.refframe,
                                       name=obs_item.name);
        rot_disp_p.append(response);
    rot_response = GF_element(disp_points=rot_disp_p, param_name=rot_name, upper_bound=1, lower_bound=-1,
                              slip_penalty=0, units='deg/Ma');
    return rot_response;


def get_GF_rotation_elements(obs_disp_points, target_region=(-180, 180, -90, 90), rot_name=''):
    """
    Build 3 GF_elements for horizontal rotation of GNSS velocities due to reference frames
    X rotation: [0, 0, 1] Euler Pole
    Y rotation: [90, 0, 1] Euler Pole
    Z rotation: [0, 89.99, 1] Euler Pole

    :param obs_disp_points: list of disp_points
    :param target_region: list of lon/lon/lat/lat for bounding box
    :param rot_name: string, optional metadata for naming the rotation (ex: ocb_)
    :returns: list of 3 GF_elements with rotation displacements in x, y, and z directions
    """
    xresponse = get_GF_rotation_element(obs_disp_points, ep=[0, 0, 1],  rot_name=rot_name+'x_rot',
                                        target_region=target_region);  # X direction
    yresponse = get_GF_rotation_element(obs_disp_points, ep=[90, 0, 1],  rot_name=rot_name+'y_rot',
                                        target_region=target_region);  # Y direction
    zresponse = get_GF_rotation_element(obs_disp_points, ep=[0, 89.99, 1],  rot_name=rot_name+'z_rot',
                                        target_region=target_region);  # Z direction
    return [xresponse, yresponse, zresponse];


def get_GF_leveling_offset_element(obs_disp_points):
    """
    Build a GF_element for a reference frame leveling offset column of the GF matrix

    :param obs_disp_points: list of disp_point_objects
    :returns: a list of 1 GF_element, or an empty list if there is no leveling in this dataset
    """
    total_response_pts = [];
    lev_count = 0;
    for item in obs_disp_points:
        if item.meas_type == "leveling":
            response = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=1,
                                           Se_obs=0, Sn_obs=0, Su_obs=0, meas_type=item.meas_type);
            lev_count += 1;
        else:
            response = Displacement_points(lon=item.lon, lat=item.lat, dE_obs=0, dN_obs=0, dU_obs=0,
                                           Se_obs=0, Sn_obs=0, Su_obs=0, meas_type=item.meas_type);
        total_response_pts.append(response);
    lev_offset_gf = GF_element(disp_points=total_response_pts, param_name='lev_offset',
                               upper_bound=1, lower_bound=-1, slip_penalty=0, units='m/yr');
    if lev_count == 0:
        return [];
    else:
        return [lev_offset_gf];


def get_displacement_directions(obs_disp_point, model_point):
    """
    Code up the logic for which components we model for each GNSS/leveling/tidegage/insar point, etc
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
    elif obs_disp_point.meas_type == "insar":
        disps = np.array([model_point.dE_obs]);
        sigmas = np.array([model_point.Se_obs]);
    else:
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs]);
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Su_obs]);
    return disps, sigmas;


def unpack_model_pred_vector(model_pred, paired_obs):
    """
    Unpack a model vector into a bunch of disp_point objects. Same logic implemented here as in the functions above.

    :param model_pred: long vector of model parameters, corresponding to each component being used
    :param paired_obs: list of disp_point_objects (shorter than model_pred vector)
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
        elif paired_obs[i].meas_type == "insar":
            E = model_pred[counter];
            N, U = np.nan, np.nan;
            counter = counter+1;
        else:
            E = model_pred[counter];
            N = model_pred[counter+1];
            U = model_pred[counter+2];
            counter = counter+3;

        new_disp_point = Displacement_points(lon=paired_obs[i].lon, lat=paired_obs[i].lat, dE_obs=E, dN_obs=N, dU_obs=U,
                                             Se_obs=0, Sn_obs=0, Su_obs=0, name=paired_obs[i].name,
                                             meas_type=paired_obs[i].meas_type, refframe=paired_obs[i].refframe);
        disp_points_list.append(new_disp_point);
    return disp_points_list;


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
    Build observation 1D-vector
    """
    obs, sigmas = [], [];
    for item in obs_disp_points:
        new_values, new_sigmas = get_displacement_directions(item, item);
        obs = np.concatenate((obs, new_values), axis=0);
        sigmas = np.concatenate((sigmas, new_sigmas), axis=0);
    return obs, sigmas;


def forward_disp_points_predictions(G, m, sigmas, paired_obs):
    """Create a convenient list of disp_points from a forward prediction based on G and m and sigma matrices/vectors."""
    model_pred = G.dot(m) * sigmas;
    model_disp_points = unpack_model_pred_vector(model_pred, paired_obs);
    return model_disp_points;


def unpack_model_of_target_param(M_vector, parameter_names, target_names=()):
    """
    Simplify a model vector into only those components from particular target_name (string), such as 'CSZ_dist'
    Useful for investigating a forward model.

    :param M_vector: 1D array of numbers, model parameter values
    :param parameter_names: list of parameters names for each model parameter
    :param target_names: list, names of desired model parameter
    :returns: vector of zeros except for the model values corresponding to target name (e.g., 5 mm/yr on one fault)
    """
    target_params_vector = np.zeros(np.shape(M_vector));
    for i in range(len(parameter_names)):
        if parameter_names[i] in target_names:
            target_params_vector[i] = 1;
    M_target = np.multiply(M_vector, target_params_vector);
    return M_target;


def unpack_model_without_target_param(M_vector, parameter_names, exclude_names=()):
    """
    Simplify a model vector into only those components excluding target_names (string)
    Useful for investigating a forward model.

    :param M_vector: 1D array of numbers, model parameter values
    :param parameter_names: list of parameters names for each model parameter
    :param exclude_names: list, names of undesired model parameter
    :returns: vector of parameter values except for the model values corresponding to exclude name
    """
    target_params_vector = np.ones(np.shape(M_vector));
    for i in range(len(parameter_names)):
        if parameter_names[i] in exclude_names:
            target_params_vector[i] = 0;
    M_target = np.multiply(M_vector, target_params_vector);
    return M_target;


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
    print("G and obs before smoothing:", np.shape(G), np.shape(obs));
    if strength == 0:
        print("No change, smoothing set to 0");
        return G, obs, sigmas;   # returning unaltered G if there is no smoothing

    G_smoothing = np.zeros((len(gf_elements), len(gf_elements)));

    # Get critical distance, a typical small distance between neighboring fault patches.
    # Operates on the first patch that it finds.
    distances = [];
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name in param_name_list:
            for j in range(i+1, len(gf_elements)):
                if gf_elements[j].param_name in param_name_list:
                    distances.append(gf_elements[i].fault_dict_list[0].get_fault_element_distance(
                        gf_elements[j].fault_dict_list[0]));
            break;
    critical_distance = sorted(distances)[2] + lengthscale;  # smooth adjacent patches with some wiggle room

    # Build the parts of the matrix for smoothing
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name in param_name_list:
            G_smoothing[i][i] = 1;
            for j in range(len(gf_elements)):
                if gf_elements[j].param_name in param_name_list:
                    if i != j and gf_elements[i].fault_dict_list[0].get_fault_element_distance(
                            gf_elements[j].fault_dict_list[0], threedimensional=distance_3d) < critical_distance:
                        G_smoothing[i][j] = laplacian_operator;

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
        if gf_elements[i].slip_penalty > 0:
            G_penalty[i][i] = gf_elements[i].slip_penalty;

    G_penalty = G_penalty * penalty;  # multiplying by lambda factor

    # observation vector of zeros
    zero_vector = np.zeros((len(gf_elements),));

    G_penalty = np.vstack((G, G_penalty));    # appending smoothing matrix
    smoothed_obs = np.concatenate((obs, zero_vector));   # appending smoothing components to data
    smoothed_sigmas = np.concatenate((sigmas, zero_vector));  # appending smoothing components to sigmas
    print("G and obs after penalty:", np.shape(G_penalty), np.shape(smoothed_obs));
    return G_penalty, smoothed_obs, smoothed_sigmas;


def filter_out_smoothing_lines(pred_vector, obs_vector, sigma_vector, tol=1e-9):
    """
    Remove lines of zeros automatically added to obs/sigma vectors for smoothing and slip penalty parameters.
    All inputs are expected to be vectors with the same length.

    :param pred_vector: array, in m
    :param obs_vector: array, in m
    :param sigma_vector: array, in m
    :param tol: optional, 1e-9
    """
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


def write_fault_traces(M_vector, paired_gf_elements, outfile, ignore_faults=()):
    """Write a gmt multi-segment file with fault traces (GF_element.points) and colors from values.
    Will only write when the GF_element has the field "points".

    :param M_vector: vector containing value of model parameters, floats (e.g., in mm/yr)
    :param paired_gf_elements: matching list of GF_elements
    :param outfile: filename, string
    :param ignore_faults: parameter names to ignore when printing gmt file
    """
    print("Writing %s" % outfile);
    with open(outfile, 'w') as ofile:
        for i in range(len(M_vector)):
            if paired_gf_elements[i].param_name in ignore_faults:
                continue;
            if len(paired_gf_elements[i].points) == 0:
                continue;
            ofile.write("> -Z%f \n" % M_vector[i]);
            for coord in paired_gf_elements[i].points:
                ofile.write("%f %f\n" % (coord[0], coord[1]) );
    return;


def write_raw_model_params(outfile, v, rms_residual=0, GF_elements=None):
    """
    The most general way of reporting the model vector, not especially human-readable.

    :param outfile: string, filename
    :param v: vector of model parameters, floats
    :param rms_residual: optional, float, mm/yr
    :param GF_elements: optional, list of paired GF_element objects to write the parameter names
    """
    print("Writing raw model outputs in %s" % outfile);
    ofile = open(outfile, 'w');
    for item in v:
        ofile.write(str(item)+"\n");
    report_string = "\nRMS: %f mm/yr\n" % rms_residual;
    ofile.write(report_string);
    if GF_elements:
        for item in GF_elements:
            ofile.write(item.param_name+" ");
    ofile.close();
    return;


def write_summary_params(v, outfile, GF_elements, ignore_faults=(), message=''):
    """
    Write a human-readable results file, with the potential to ignore some faults for clarity.

    :param v: vector of model parameters, floats
    :param outfile: string
    :param GF_elements: list of GF_element objects
    :param ignore_faults: list of strings
    :param message: optional message from inverse routine about solution quality
    """
    print("Writing %s" % outfile);
    ofile = open(outfile, 'w');
    for i in range(len(v)):
        if GF_elements[i].param_name in ignore_faults:
            continue;
        ofile.write(GF_elements[i].param_name + ": ");
        if GF_elements[i].units == "cm/yr":   # converting fault slip rates into mm/yr for convenience
            units = 'mm/yr';
            multiplier = 10;
        else:
            units = GF_elements[i].units;
            multiplier = 1;
        ofile.write("%.5f %s" % (v[i]*multiplier, units));
        ofile.write('  [within %.3f to %.3f]' % (GF_elements[i].lower_bound*multiplier,
                                                 GF_elements[i].upper_bound*multiplier) );
        ofile.write("\n");
    report_string = "\nWith %d observations\n" % (len(GF_elements[0].disp_points));
    ofile.write(report_string);
    ofile.write("Message: "+message+"\n");
    ofile.close();
    return;


def view_full_results(exp_dict, paired_obs, modeled_disp_points, residual_pts, rotation_pts, norot_pts, title, region):
    # Plot the data, model, and residual in separate plots
    # Not plotting the fault patches because it takes a long time.
    scale_arrow = (1.0, 0.010, "1 cm/yr");
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"] + "/data_only.png",
                                                         disp_points=paired_obs, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/residuals.png",
                                                         disp_points=residual_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001,
                                                         title=title);
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/model_pred.png",
                                                         disp_points=modeled_disp_points, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/rotation_pred.png",
                                                         disp_points=rotation_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"]+"/faults_only_pred.png",
                                                         disp_points=norot_pts, region=region,
                                                         scale_arrow=scale_arrow, v_labeling_interval=0.001)
    return;


def visualize_GF_elements(GF_elements_list, outdir, exclude_list=()):
    """
    Aside to the main calculation, just view each GF.

    :param GF_elements_list: list of GF_elements objects
    :param outdir: string for outdir
    :param exclude_list: optional list of GF_element.param_name to exclude from visualizing
    """
    if exclude_list == 'all':
        return;
    for GF_el in GF_elements_list:
        if GF_el.param_name in exclude_list:   # plot these elements separately, like individual CSZ patches
            continue;
        print(GF_el.param_name);
        if GF_el.param_name == "CSZ":
            scale_arrow = (1.0, 0.010, "1 cm");
        else:
            scale_arrow = (1.0, 0.001, "1 mm");
        library.plot_fault_slip.map_source_slip_distribution(GF_el.fault_dict_list, outdir + "/gf_" +
                                                             GF_el.param_name + "_only.png",
                                                             disp_points=GF_el.disp_points,
                                                             region=[-127, -119.7, 37.7, 43.3],
                                                             scale_arrow=scale_arrow,
                                                             v_labeling_interval=0.001);
    return;


def remove_nearfault_pts(obs_points, fault_trace_file):
    """
    Wraps dpo.utilities.filter_to_remove_near_fault so that we can pass a filename.

    :param obs_points: a list of disp_point objects
    :param fault_trace_file: filename that contains two columns of numbers for fault trace: lon, lat
    """
    trace_pts = np.loadtxt(fault_trace_file);
    print("Removing near-fault points from file %s " % fault_trace_file);
    obs_disp_points = dpo.utilities.filter_to_remove_near_fault(obs_points, trace_pts, radius_km=10);
    return obs_disp_points;


def write_insar_greens_functions(GF_elements, outfile):
    """
    Serialize a bunch of InSAR Green's Functions into written text file, in meters, with rows for each InSAR point:
    For each observation point: lon, lat, dLOS1, dLOS2, dLOS3.... [n fault patches].

    :param GF_elements: list of GF_elements with InSAR displacements in disp_points.
    :param outfile: string
    """
    print("Writing file %s " % outfile);
    for item in GF_elements:
        if item.meas_type != 'insar':
            raise ValueError("Error! In this function, we can only serialize GF_elements with meas_type == insar.");
    ofile = open(outfile, 'w');
    for i, pt in enumerate(GF_elements[0].disp_points):
        ofile.write('%f %f ' % (pt.lon, pt.lat));
        point_displacements = [GF_el.disp_points[i].dE_obs for GF_el in GF_elements];
        for x in point_displacements:
            ofile.write(str(x)+" ");
        ofile.write("\n");
    ofile.close();
    return;


def read_insar_greens_functions(gf_file, fault_patches, param_name='', lower_bound=0, upper_bound=0):
    """
    Read pre-computed green's functions in the format matching the write-function.
    Currently, only works for triangles.

    :param gf_file: string, filename
    :param fault_patches: list of fault_slip_objects or fault_slip_triangles
    :param param_name: string
    :param lower_bound: float
    :param upper_bound: float
    """
    GF_elements = [];
    gf_data_array = np.loadtxt(gf_file);
    lons, lats = gf_data_array[:, 0], gf_data_array[:, 1];
    model_disp_pts = [];
    for tlon, tlat in zip(lons, lats):
        mdp = Displacement_points(lon=tlon, lat=tlat, dE_obs=0, dN_obs=0, dU_obs=0, Se_obs=0, Sn_obs=0, Su_obs=0,
                                  meas_type='insar');
        model_disp_pts.append(mdp);

    for i, tri in enumerate(fault_patches):  # right now, only the triangle version is written.
        changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0);  # triangle-specific
        changed_slip = changed_slip.change_reference_loc();  # triangle-specific interface
        index = i+2;  # moving to the correct column in the GF file, skipping lon and lat.
        los_defo = gf_data_array[:, index];
        model_disp_pts = dpo.utilities.set_east(model_disp_pts, los_defo);
        GF_elements.append(GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                      param_name=param_name, lower_bound=lower_bound, upper_bound=upper_bound,
                                      slip_penalty=0));
    return GF_elements;
