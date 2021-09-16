
import numpy as np
from Tectonic_Utilities.Tectonic_Utils.geodesy import euler_pole
import Elastic_stresses_py.PyCoulomb.coulomb_collections as cc
import Elastic_stresses_py.PyCoulomb.fault_slip_object.disp_points_object as disp_points


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


def add_all_csz_patches(one_patch_gfs):
    """take n patches of the subduction interface and add their green's functions together """
    new_pts = one_patch_gfs[0];
    for i in range(1, len(one_patch_gfs)):
        new_pts = disp_points.add_disp_points(new_pts, one_patch_gfs[i]);
    return new_pts;


def add_G_rotation_columns(G, obs_disp_points):
    """
    Build 3 columns of the G matrix for horizontal rotation of GNSS velocities due to reference frames
    """
    xrot_col, yrot_col, zrot_col = [], [], [];
    for obs_item in obs_disp_points:
        coords = [obs_item.lon, obs_item.lat];
        response_to_xrot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 0, 1]);   # X-dir EP
        response_to_yrot = euler_pole.point_rotation_by_Euler_Pole(coords, [90, 0, 1]);  # Y-dir EP
        response_to_zrot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 89.99, 1]);   # Z-dir EP

        response_to_xrot = clip_response_to_rot(response_to_xrot, obs_item);
        response_to_yrot = clip_response_to_rot(response_to_yrot, obs_item);
        response_to_zrot = clip_response_to_rot(response_to_zrot, obs_item);

        xrot_col = np.concatenate((xrot_col, np.array(response_to_xrot)), axis=0);
        yrot_col = np.concatenate((yrot_col, np.array(response_to_yrot)), axis=0);
        zrot_col = np.concatenate((zrot_col, np.array(response_to_zrot)), axis=0);

    # Once you're done, add the columns for Euler Pole rotations
    xrot_col = np.reshape(xrot_col, (len(xrot_col), 1));
    yrot_col = np.reshape(yrot_col, (len(xrot_col), 1));
    zrot_col = np.reshape(zrot_col, (len(xrot_col), 1));
    G_with_rot = np.concatenate((G, xrot_col, yrot_col, zrot_col), axis=1)
    return G_with_rot;


def get_displacement_directions(obs_disp_point, model_point):
    """Code up the logic for which components we model for each GNSS/leveling/tidegage point, etc """
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


def clip_response_to_rot(partials, obs_disp_point):
    """
    Extract the partial derivatives that we will use for a certain disp_point.
    Should match the logic from the function above.
    """
    if obs_disp_point.meas_type == "continuous":
        return np.array([partials[0], partials[1], partials[2]]);
    elif obs_disp_point.meas_type == "survey":
        return np.array([partials[0], partials[1]]);
    elif obs_disp_point.meas_type == "leveling":
        return np.array([partials[2]]);
    elif obs_disp_point.meas_type == "tide_gage":
        return np.array([partials[2]]);
    else:
        return np.array([partials[0], partials[1], partials[2]]);


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


def write_model_params(v, residual, outfile, fault_names=None):
    print("Writing %s" % outfile);
    ofile = open(outfile, 'w');
    for item in v:
        ofile.write(str(item)+"\n");
    report_string = "\nRMS: %f mm/yr\n" % residual;
    ofile.write(report_string);
    if fault_names:
        for item in fault_names:
            ofile.write(item+" ");
    ofile.close();
    return;
