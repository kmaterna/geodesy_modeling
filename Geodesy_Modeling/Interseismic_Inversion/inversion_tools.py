
import numpy as np
import collections
from Tectonic_Utilities.Tectonic_Utils.geodesy import euler_pole
import Elastic_stresses_py.PyCoulomb.coulomb_collections as cc
import Elastic_stresses_py.PyCoulomb.fault_slip_object.disp_points_object as disp_points
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library

"""
GF_element is everything you would need to make a column of the Green's matrix and plot the impulse response function. 
"""
GF_element = collections.namedtuple('GF_element', ['disp_points', 'fault_name',
                                                   'fault_dict_list', 'upper_bound', 'lower_bound']);


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
    Take a list of GF_elements, and a list of obs_disp_points. Pare them down to a matching set of points.
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
                                             upper_bound=gf_model.upper_bound));
        if len(paired_gf) != target_len:
            raise ValueError("ERROR! Not all points have green's functions.");
    return paired_obs, paired_gf_elements;


def add_all_csz_patches(one_patch_gfs):
    """take n patches of the subduction interface and add their green's functions together """
    new_pts = one_patch_gfs[0];
    for i in range(1, len(one_patch_gfs)):
        new_pts = disp_points.add_disp_points(new_pts, one_patch_gfs[i]);
    return new_pts;


def get_G_rotation_elements(obs_disp_points):
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

    xresponse = GF_element(disp_points=x_disp_p, fault_name='x_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1);
    yresponse = GF_element(disp_points=y_disp_p, fault_name='y_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1);
    zresponse = GF_element(disp_points=z_disp_p, fault_name='z_rot', fault_dict_list=[], upper_bound=1, lower_bound=-1);
    return [xresponse, yresponse, zresponse];


def get_displacement_directions(obs_disp_point, model_point):
    """
    Code up the logic for which components we model for each GNSS/leveling/tidegage point, etc
    """
    if obs_disp_point.meas_type == "continuous":
        disps = np.array([model_point.dE_obs, model_point.dN_obs, model_point.dU_obs]);
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Sn_obs]);
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
        sigmas = np.array([model_point.Se_obs, model_point.Sn_obs, model_point.Sn_obs]);
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
        sigmas = np.concatenate((sigmas, new_values), axis=0);
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


def unpack_model_of_rotation_only(M_vector):
    """
    1) Zeros out all model parameters except the rotation (leaving the last three)
    2) Zeros out all the rotation (canceling the last three)
    """
    multiplier = np.zeros(np.shape(M_vector));
    multiplier[-3:] = [1, 1, 1];
    M_rot_only = np.multiply(M_vector, multiplier);
    M_no_rot = np.zeros(np.shape(M_vector));
    M_no_rot[0:-3] = M_vector[0:-3];
    return M_rot_only, M_no_rot;

def rms_from_model_pred_vector(pred_vector, obs_vector):
    """Various ways to compute the RMS value"""
    residuals = np.subtract(obs_vector, pred_vector);
    rms_mm = np.multiply(np.sqrt(np.mean(np.square(residuals))), 1000);
    return rms_mm;


def write_model_params(v, residual, outfile, GF_elements=None):
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


def view_full_results(exp_dict, paired_obs, modeled_disp_points, residual_pts, rotation_pts, norot_pts, title, region):
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
