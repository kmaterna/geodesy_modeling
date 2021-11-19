#!/usr/bin/env python
"""
An example driver for interseismic velocity inversion for fault slip rates.
Humboldt Bay application.
Mostly research code.
Reaches into Humboldt Bay code for some functions.
Second half of 2021.
"""


import numpy as np
import scipy.optimize
import subprocess, json, sys
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library
import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Geodesy_Modeling.src.Interseismic_Inversion as Interseismic_Inversion
sys.path.append("/Users/kmaterna/Documents/B_Research/Mendocino_Geodesy/Humboldt/_Project_Code");   # add local code to path
import humboldt_readers as HR


def configure(config_file):
    """Set up an experiment directory."""
    with open(config_file) as f:
        exp_dict = json.load(f);
    with open(exp_dict['fault_data']) as f:
        metadata = json.load(f);
    exp_dict["faults"] = metadata["faults"];
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);
    subprocess.call(['cp', config_file, exp_dict["outdir"]+"/config.txt"]);
    print("Success resolving imports")
    return exp_dict;

def correct_for_far_field_terms(exp_dict, obs_disp_points):
    """
    Velocity corrections, Pollitz & Evans, 2017
    """
    csz_correction_disp_points = HR.read_correction_data_table(exp_dict["csz_correction"]);
    ridge_correction_disp_points = HR.read_correction_data_table(exp_dict["ridge_correction"]);  # Fred correction
    ssaf_correction_disp_points = HR.read_correction_data_table(exp_dict["ssaf_correction"]);  # Fred correction
    obs_disp_points = library.disp_points_object.subtract_disp_points(obs_disp_points, csz_correction_disp_points);
    obs_disp_points = library.disp_points_object.subtract_disp_points(obs_disp_points, ridge_correction_disp_points);
    obs_disp_points = library.disp_points_object.subtract_disp_points(obs_disp_points, ssaf_correction_disp_points);
    return obs_disp_points;

def read_fault_gf_elements(exp_dict):
    """
    Input: a config dictionary
    Return: a list of inversion_tools.GF_elements,
    which are the building blocks for the columns of the Green's function inversion.
    """
    # Read a list of green's functions for the modeled faults in this experiment.  One for each model parameter.
    gf_elements = [];  # building blocks for columns in the Green's matrix
    for i in range(len(exp_dict["exp_faults"])):  # for each fault
        fault_name = exp_dict["exp_faults"][i];
        if fault_name == "CSZ":  # Reading for single-parameter CSZ case
            one_patch_gfs, csz_patches = Interseismic_Inversion.readers.read_distributed_GF(exp_dict["faults"]["CSZ"]["GF"],
                                                                         exp_dict["faults"]["CSZ"]["geometry"],
                                                                         exp_dict["lonlatfile"]);
            collective_csz_gf = Interseismic_Inversion.inversion_tools.add_all_csz_patches(one_patch_gfs);
            one_gf_element = Interseismic_Inversion.inversion_tools.GF_element(disp_points=collective_csz_gf, fault_name=fault_name,
                                                        fault_dict_list=csz_patches,
                                                        lower_bound=exp_dict["faults"][fault_name]["slip_min"],
                                                        upper_bound=exp_dict["faults"][fault_name]["slip_max"]);
            gf_elements.append(one_gf_element);
        elif fault_name == "CSZ_dist":  # Reading for distributed CSZ patches as unit slip.
            one_patch_gfs, csz_patches = Interseismic_Inversion.readers.read_distributed_GF(exp_dict["faults"]["CSZ"]["GF"],
                                                                         exp_dict["faults"]["CSZ"]["geometry"],
                                                                         exp_dict["lonlatfile"], unit_slip=True);
            for gf_disp_points, patch in zip(one_patch_gfs, csz_patches):
                one_gf_element = Interseismic_Inversion.inversion_tools.GF_element(disp_points=gf_disp_points, fault_name=fault_name,
                                                            fault_dict_list=[patch],
                                                            lower_bound=exp_dict["faults"]["CSZ"]["slip_min"],
                                                            upper_bound=exp_dict["faults"]["CSZ"]["slip_max"]);
                gf_elements.append(one_gf_element);
        else:  # Reading for LSF, MRF, other fault cases
            fault_gf = exp_dict["faults"][fault_name]["GF"];
            fault_geom = exp_dict["faults"][fault_name]["geometry"];
            temp, _ = library.io_static1d.read_static1D_source_file(fault_geom, headerlines=1);
            mod_disp_points = library.io_static1d.read_static1D_output_file(fault_gf, exp_dict["lonlatfile"]);
            one_gf_element = Interseismic_Inversion.inversion_tools.GF_element(disp_points=mod_disp_points, fault_name=fault_name,
                                                        fault_dict_list=temp,
                                                        lower_bound=exp_dict["faults"][fault_name]["slip_min"],
                                                        upper_bound=exp_dict["faults"][fault_name]["slip_max"]);
            gf_elements.append(one_gf_element);
    return gf_elements;

def run_humboldt_inversion(config_file):
    # Starting program.  Configure stage
    exp_dict = configure(config_file);

    # # INPUT stage: Read obs velocities as cc.Displacement_Points
    obs_disp_points = HR.read_all_data_table(exp_dict["data_file"]);
    obs_disp_points = correct_for_far_field_terms(exp_dict, obs_disp_points);  # needed from Fred's work
    # Experimental options:
    obs_disp_points = HR.filter_to_continuous_only(obs_disp_points);  # an experimental design step.

    # INPUT stage: Read GF models based on the configuration parameters
    gf_elements = read_fault_gf_elements(exp_dict);  # list of GF_elements, one for each fault-related column of G.

    # COMPUTE STAGE: PREPARE ROTATION GREENS FUNCTIONS AND LEVELING OFFSET
    gf_elements_rotation = Interseismic_Inversion.inversion_tools.get_GF_rotation_elements(obs_disp_points);  # 3 elements: rot_x, rot_y, rot_z
    gf_elements = gf_elements + gf_elements_rotation;  # add rotation elements to matrix
    gf_element_lev = Interseismic_Inversion.inversion_tools.get_GF_leveling_offset_element(obs_disp_points);  # 1 element: lev reference frame
    gf_elements = gf_elements + gf_element_lev;

    # COMPUTE STAGE: PAIRING DATA AND GREENS FUNCTION POINTS.
    paired_obs, paired_gf_elements = Interseismic_Inversion.inversion_tools.pair_gf_elements_with_obs(obs_disp_points, gf_elements);

    exclude_list = ["CSZ_dist", "NSAF", "x_rot", "y_rot", "z_rot", "CSZ", "lev_offset", "MadRiver", "LSF", "BSF",
                    "Maacama"]
    Interseismic_Inversion.inversion_tools.visualize_GF_elements(paired_gf_elements, exp_dict["outdir"], exclude_list=exclude_list);

    # COMPUTE STAGE: INVERSE.  Reduces certain points to only-horizontal, only-vertical, etc.
    list_of_gf_columns = [];
    for paired_gf in paired_gf_elements:
        G_one_col = Interseismic_Inversion.inversion_tools.buildG_column(paired_gf.disp_points, paired_obs);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);

    # Build observation vector
    obs, _sigmas = Interseismic_Inversion.inversion_tools.build_obs_vector(paired_obs);
    sigmas = np.ones(np.shape(obs));  # placeholder.
    G /= sigmas[:, None];
    weighted_obs = obs / sigmas;

    # Add optional smoothing penalty, overwriting old variables
    if 'smoothing' in exp_dict.keys():
        G, weighted_obs, sigmas = Interseismic_Inversion.inversion_tools.build_smoothing(paired_gf_elements, 'CSZ_dist',
                                                                  exp_dict["smoothing"], G, weighted_obs, sigmas);

    # Money line: Constrained inversion
    lower_bound = [x.lower_bound for x in paired_gf_elements];
    upper_bound = [x.upper_bound for x in paired_gf_elements];
    response = scipy.optimize.lsq_linear(G, weighted_obs, bounds=(lower_bound, upper_bound), max_iter=1500,
                                         method='bvls');
    M_opt = response.x;  # parameters of best-fitting model
    print(response.message);
    print("First fault slip rate:", M_opt[0])
    if response.message == "The maximum number of iterations is exceeded.":
        print("Maximum number of iterations exceeded. Cannot trust this inversion. Exiting");
        sys.exit(0);

    def forward_disp_points_predictions(G, m, sigmas, paired_obs):
        model_pred = G.dot(m) * sigmas;
        model_disp_points = Interseismic_Inversion.inversion_tools.unpack_model_pred_vector(model_pred, paired_obs);
        return model_disp_points;

    # Make forward predictions
    M_rot_only, M_no_rot = Interseismic_Inversion.inversion_tools.unpack_model_of_rotation_only(M_opt,
                                                                         [x.fault_name for x in paired_gf_elements]);
    modeled_disp_points = forward_disp_points_predictions(G, M_opt, sigmas, paired_obs);
    rot_modeled_pts = forward_disp_points_predictions(G, M_rot_only, sigmas, paired_obs);
    norot_modeled_pts = forward_disp_points_predictions(G, M_no_rot, sigmas, paired_obs);

    # Output stage
    fault_dict_lists = [item.fault_dict_list for item in paired_gf_elements];
    rms_mm, rms_chi2 = Interseismic_Inversion.inversion_tools.rms_from_model_pred_vector(G.dot(M_opt) * sigmas, weighted_obs * sigmas, sigmas);
    rms_title = "RMS: %f mm/yr" % rms_mm;
    print(" ", rms_title);
    residual_pts = library.disp_points_object.subtract_disp_points(paired_obs, modeled_disp_points);
    PyCoulomb.io_additionals.write_disp_points_results(modeled_disp_points,
                                                       exp_dict["outdir"] + '/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(paired_obs, exp_dict["outdir"] + '/simple_obs_file.txt');
    Interseismic_Inversion.inversion_tools.write_model_params(M_opt, rms_mm, exp_dict["outdir"] + '/' + exp_dict["model_file"],
                                       paired_gf_elements);
    Interseismic_Inversion.inversion_tools.write_summary_params(M_opt, rms_mm, exp_dict["outdir"] + '/model_results_human.txt',
                                         paired_gf_elements,
                                         ignore_faults=['CSZ_dist']);
    Interseismic_Inversion.readers.write_csz_dist_fault_patches(fault_dict_lists, M_opt, exp_dict["outdir"] + '/csz_model.gmt');
    Interseismic_Inversion.inversion_tools.view_full_results(exp_dict, paired_obs, modeled_disp_points, residual_pts, rot_modeled_pts,
                                      norot_modeled_pts, rms_title, region=[-127, -119.7, 37.7, 43.5]);
    library.plot_fault_slip.plot_data_model_residual(exp_dict["outdir"] + "/results.png", paired_obs,
                                                     modeled_disp_points, residual_pts, [-126, -119.7, 37.7, 43.3],
                                                     scale_arrow=(0.5, 0.020, "2 cm"), v_labeling_interval=0.001,
                                                     fault_dict_list=[], rms=rms_mm);
    return;


if __name__ == "__main__":
    config_file = sys.argv[1];  # take config file from runstring.
    run_humboldt_inversion(config_file);
