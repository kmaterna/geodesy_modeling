import numpy as np
from Elastic_stresses_py.PyCoulomb.disp_points_object import utilities, compute_rms

# The usual metrics for the whole project are done with operations on displacement point objects

def print_typical_uncs(obs_points):
    """
    Print a median uncertainty in mm for each data type
    """
    unc_tide_gage = [1000*x.Su_obs for x in obs_points if x.meas_type == 'tide_gage']
    unc_leveling = [1000*x.Su_obs for x in obs_points if x.meas_type == 'leveling']
    unc_east_campaign = [1000*x.Se_obs for x in obs_points if x.meas_type == 'survey']
    unc_north_campaign = [1000 * x.Sn_obs for x in obs_points if x.meas_type == 'survey']
    unc_survey = unc_east_campaign + unc_north_campaign
    unc_east_continuous = [1000*x.Se_obs for x in obs_points if x.meas_type == 'continuous']
    unc_north_continuous = [1000 * x.Sn_obs for x in obs_points if x.meas_type == 'continuous']
    unc_up_continuous = [1000 * x.Su_obs for x in obs_points if x.meas_type == 'continuous']
    unc_continuous = unc_east_continuous + unc_north_continuous + unc_up_continuous
    print("continuous:", np.median(unc_continuous), np.mean(unc_continuous))
    print("survey:", np.median(unc_survey), np.mean(unc_survey))
    print("leveling:", np.median(unc_leveling), np.mean(unc_leveling))
    print("tide gage:", np.median(unc_tide_gage), np.mean(unc_tide_gage))
    return


# ----------- IMPLEMENTATIONS OF MISFIT CALCULATIONS ------------------ #

def obs_vs_model_L2_aggregate(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points)   # make residual points
    res_cgps = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"),
                                                     ("E", "N", "U"))
    res_sgps = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "survey"), ("E", "N"))
    res_lev = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "leveling"), ("U",))
    res_tg = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "tide_gage"), ("U",))
    resid_data = res_cgps[0] + res_sgps[0] + res_lev[0] + res_tg[0]
    resid_sigma = res_cgps[1] + res_sgps[1] + res_lev[1] + res_tg[1]
    rms_m, reported_chi2 = compute_rms.L2_on_vector(resid_data, resid_sigma)
    rms_mm = np.multiply(rms_m, 1000)
    return rms_mm, reported_chi2


def obs_vs_model_L2_horiz(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points)   # make residual points
    res_cgps = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"), ("E", "N"))
    res_sgps = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "survey"), ("E", "N"))
    resid_data = res_cgps[0] + res_sgps[0]
    resid_sigma = res_cgps[1] + res_sgps[1]
    rms_m, reported_chi2 = compute_rms.L2_on_vector(resid_data, resid_sigma)
    rms_mm = np.multiply(rms_m, 1000)
    return rms_mm, reported_chi2


def obs_vs_model_L2_vertical(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points)   # make residual points
    res_cgps = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"), ("U",))
    res_lev = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "leveling"), ("U",))
    res_tg = compute_rms.filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "tide_gage"), ("U",))
    resid_data = res_cgps[0] + res_lev[0] + res_tg[0]
    resid_sigma = res_cgps[1] + res_lev[1] + res_tg[1]
    rms_m, reported_chi2 = compute_rms.L2_on_vector(resid_data, resid_sigma)
    rms_mm = np.multiply(rms_m, 1000)
    return rms_mm, reported_chi2
