
import numpy as np
from Tectonic_Utilities.Tectonic_Utils.geodesy import euler_pole
import Elastic_stresses_py.PyCoulomb.coulomb_collections as cc


def pair_obs_gf(obs_disp_points, model_disp_points):
    """
    Operates on two disp_points objects
    """
    paired_gf, paired_obs = [], [];
    for obs_item in obs_disp_points:
        for gf_item in model_disp_points:
            if abs(obs_item.lon - gf_item.lon) < 0.001 and abs(obs_item.lat - gf_item.lat) < 0.001:
                paired_gf.append(gf_item);
                paired_obs.append(obs_item);
                break;
    return paired_obs, paired_gf;


def add_G_rotation_columns(G, obs_disp_points):
    """
    Build 3 columns of the G matrix for horizontal rotation of GNSS velocities due to reference frames
    """
    xrot_col, yrot_col, zrot_col = [], [], [];
    for item in obs_disp_points:
        coords = [item.lon, item.lat];
        response_to_xrot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 0, 1]);   # X-dir EP
        response_to_yrot = euler_pole.point_rotation_by_Euler_Pole(coords, [90, 0, 1]);  # Y-dir EP
        response_to_zrot = euler_pole.point_rotation_by_Euler_Pole(coords, [0, 89.99, 1]);   # Z-dir EP
        xrot_col = np.concatenate((xrot_col, np.array(response_to_xrot)), axis=0);
        yrot_col = np.concatenate((yrot_col, np.array(response_to_yrot)), axis=0);
        zrot_col = np.concatenate((zrot_col, np.array(response_to_zrot)), axis=0);

    # Once you're done, add the columns for Euler Pole rotations
    xrot_col = np.reshape(xrot_col, (len(xrot_col), 1));
    yrot_col = np.reshape(yrot_col, (len(xrot_col), 1));
    zrot_col = np.reshape(zrot_col, (len(xrot_col), 1));
    G_with_rot = np.concatenate((G, xrot_col, yrot_col, zrot_col), axis=1)
    return G_with_rot;


def buildG_column(GF_disp_points):
    """
    Green's functions for a single model element (i.e., slip on one fault).
    Returns a single column of nx1.
    """
    GF_col = [];
    for item in GF_disp_points:
        GF_col = np.concatenate((GF_col, np.array([item.dE_obs])));
        GF_col = np.concatenate((GF_col, np.array([item.dN_obs])));
        GF_col = np.concatenate((GF_col, np.array([item.dU_obs])));
    GF_col = np.reshape(GF_col, (len(GF_col), 1));
    return GF_col;


def build_obs_vector(obs_disp_points):
    """
    Build observation vector (will eventually have leveling here too)
    """
    obs, sigmas = [], [];
    for item in obs_disp_points:
        obs = np.concatenate((obs, np.array([item.dE_obs])), axis=0);
        obs = np.concatenate((obs, np.array([item.dN_obs])), axis=0);
        obs = np.concatenate((obs, np.array([item.dU_obs])), axis=0);
        sigmas = np.concatenate((sigmas, np.array([item.Se_obs])), axis=0);
        sigmas = np.concatenate((sigmas, np.array([item.Sn_obs])), axis=0);
        sigmas = np.concatenate((sigmas, np.array([item.Su_obs])), axis=0);
    return obs, sigmas;


def unpack_model_pred_vector(model_pred, paired_gf):
    """
    When the model vector changes and has lots of leveling, this will also change.
    """
    disp_points_list = [];
    for i in range(len(paired_gf)):
        new_disp_point = cc.Displacement_points(lon=paired_gf[i].lon, lat=paired_gf[i].lat,
                                                dE_obs=model_pred[i*3],
                                                dN_obs=model_pred[i*3 + 1],
                                                dU_obs=model_pred[i*3 + 2],
                                                Se_obs=0,
                                                Sn_obs=0,
                                                Su_obs=0, name=paired_gf[i].name);
        disp_points_list.append(new_disp_point);
    return disp_points_list;


def write_model_params(v, residual, outfile):
    ofile = open(outfile, 'w');
    for item in v:
        ofile.write(str(item)+"\n");
    report_string = "\nRMS: %f mm/yr" % residual;
    ofile.write(report_string);
    ofile.close();
    return;
