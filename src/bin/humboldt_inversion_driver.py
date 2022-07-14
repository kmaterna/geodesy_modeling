#!/usr/bin/env python

"""
A driver for interseismic velocity inversion for fault slip rates.
Humboldt Bay application.
Mostly research code.  It is included in the repo as an example of what to do with this toolbox.
Depends on Humboldt Bay project code for some functions. It will not work on a different general system.
"""

import numpy as np
import scipy.optimize
import subprocess, json, sys, argparse, os
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library
import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import Geodesy_Modeling.src.Inversion.readers as readers
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.disp_points_object.outputs as dpo_out
sys.path.append("/Users/kmaterna/Documents/B_Research/Mendocino_Geodesy/Humboldt/_Project_Code");  # add local code
import humboldt_readers as HR  # Also had to add into pycharm project settings.

reader_dictionary = {
    "csz": HR.read_addcp_correction_table,
    "ridge": HR.read_correction_data_table,
    "fred": HR.read_correction_data_table,
    "ghost": HR.read_ghost_transient_table,
    "ep": HR.get_euler_pole_correction,
    "bc": library.io_static1d.read_static1D_output_file
}

def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('configfile', default='config.json');  # required first argument
    p.add_argument('faultconfigfile', default='fault_config.json');  # required second argument
    p.add_argument('--data_file', type=str, help='''filename with velocities''');
    p.add_argument('--lonlatfile', type=str, help='''filename with lonlats''');
    p.add_argument('--corrections', type=list, help='''List of correction dictionaries''');
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''');
    p.add_argument('--slip_penalty', type=float, help='''strength of minimum-norm slip constraint''');
    p.add_argument('--max_depth_csz_slip', type=float, help='''Maximum depth of slip distribution on CSZ''');
    p.add_argument('--bbox', type=list, help='''Area being included in the data [W, E, S, N]''');
    p.add_argument('--exclude_regions', type=list, help='''Areas optionally excluded from the data [[W, E, S, N],]''');
    p.add_argument('--continuous_only', type=float, help='''Flag to use only continuous GNSS''');
    p.add_argument('--outdir', type=str, help='''Output directory''');
    p.add_argument('--model_file', type=str, help='''Name of model file''');
    p.add_argument('--exp_faults', type=list, help='''Faults used in this inversion''');
    p.add_argument('--inverse_dir', type=str, help='''Way to get to the home directory of inverses''');
    exp_dict = vars(p.parse_args())

    if os.path.exists(exp_dict["configfile"]):
        config_file = open(exp_dict["configfile"], 'r');
        exp_dict = json.load(config_file)
    else:
        exp_dict = {}
    p.set_defaults(**exp_dict);
    exp_dict = vars(p.parse_args())
    exp_dict["lonlatfile"] = exp_dict["inverse_dir"] + exp_dict["lonlatfile"];
    exp_dict["data_file"] = exp_dict["inverse_dir"] + exp_dict["data_file"];
    with open(exp_dict["faultconfigfile"]) as f:   # """ Load the faults"""
        exp_dict["faults"] = json.load(f)["faults"];
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4);
    print("Success resolving imports");
    return exp_dict;


def correct_for_far_field_terms(exp_dict, obs_disp_points):
    """
    Velocity corrections, to set boundary conditions, some from Pollitz & Evans, 2017.
    """
    for correction in exp_dict["corrections"]:
        correction_dps = reader_dictionary[correction["type"]](exp_dict["inverse_dir"]+correction["file"],
                                                               exp_dict["lonlatfile"]);
        correction_dps = dpo.utilities.mult_disp_points_by(correction_dps, correction["scale"]);
        if correction["type"] == "csz":
            obs_disp_points = dpo.utilities.subtract_disp_points(obs_disp_points, correction_dps, tol=0.001,
                                                                 target='horizontal');
        else:
            obs_disp_points = dpo.utilities.subtract_disp_points(obs_disp_points, correction_dps, tol=0.001);
    return obs_disp_points;


def read_hb_fault_gf_elements(exp_dict):
    """
    Input: a config dictionary
    Return: a list of inversion_tools.GF_elements,
    which are the building blocks for the columns of the Green's function inversion.
    Read a list of green's functions for the modeled faults in this experiment.  One for each model parameter.
    """
    gf_elements = [];  # building blocks for columns in the Green's matrix
    for i in range(len(exp_dict["exp_faults"])):  # for each fault
        fault_name = exp_dict["exp_faults"][i];
        if fault_name == "CSZ_dist":  # Reading for distributed CSZ patches as unit slip.
            one_patch_dps, csz_patches, maxslip = readers.read_distributed_GF(exp_dict["inverse_dir"]+exp_dict["faults"]["CSZ"]["GF"],
                                                                              exp_dict["inverse_dir"]+exp_dict["faults"]["CSZ"]["geometry"],
                                                                              exp_dict["lonlatfile"], unit_slip=True,
                                                                              latlonbox=(-127, -120, 38, 44.5));
            for gf_disp_points, patch, max0 in zip(one_patch_dps, csz_patches, maxslip):
                amount_of_slip_penalty = 1;
                if patch["depth"] > exp_dict["max_depth_csz_slip"]:
                    amount_of_slip_penalty = 100;   # optionally: force CSZ slip to be above a certain depth
                one_gf_element = inv_tools.GF_element(disp_points=gf_disp_points, fault_name=fault_name,
                                                      fault_dict_list=[patch],
                                                      lower_bound=exp_dict["faults"]["CSZ"]["slip_min"],
                                                      upper_bound=max0*130,  # max slip from geometry, in units of cm
                                                      slip_penalty=amount_of_slip_penalty, units='cm/yr',
                                                      points=[]);
                gf_elements.append(one_gf_element);
        else:  # Reading for LSF, MRF, other fault cases
            fault_gf = exp_dict["inverse_dir"]+exp_dict["faults"][fault_name]["GF"];
            fault_geom = exp_dict["inverse_dir"]+exp_dict["faults"][fault_name]["geometry"];
            temp, _ = library.io_static1d.read_static1D_source_file(fault_geom, headerlines=1);
            mod_disp_points = library.io_static1d.read_static1D_output_file(fault_gf, exp_dict["lonlatfile"]);
            fault_points = np.loadtxt(exp_dict["inverse_dir"]+exp_dict["faults"][fault_name]["points"]);
            if "creep_multiplier" in exp_dict["faults"][fault_name].keys():
                correction_strength = exp_dict["faults"][fault_name]["creep_multiplier"];
                if correction_strength > 0:
                    correction_file1 = exp_dict["inverse_dir"]+exp_dict["faults"][fault_name]["GF_15km_visco"]
                    correction_file2 = exp_dict["inverse_dir"]+exp_dict["faults"][fault_name]["GF_15km_stat"]
                    mod_dpo1 = library.io_static1d.read_static1D_output_file(correction_file1, exp_dict["lonlatfile"]);
                    mod_dpo1 = dpo.utilities.mult_disp_points_by(mod_dpo1, 1/500);
                    mod_dpo2 = library.io_static1d.read_static1D_output_file(correction_file2, exp_dict["lonlatfile"]);
                    mod_disp_points = dpo.utilities.add_disp_points(mod_disp_points, mod_dpo1);
                    mod_disp_points = dpo.utilities.add_disp_points(mod_disp_points, mod_dpo2);
            one_gf_element = inv_tools.GF_element(disp_points=mod_disp_points, fault_name=fault_name,
                                                  fault_dict_list=temp,
                                                  lower_bound=exp_dict["faults"][fault_name]["slip_min"],
                                                  upper_bound=exp_dict["faults"][fault_name]["slip_max"],
                                                  slip_penalty=0, units='cm/yr', points=fault_points);
            gf_elements.append(one_gf_element);
    return gf_elements;


def run_humboldt_inversion():
    # Starting program.  Configure stage
    exp_dict = configure();

    # # INPUT stage: Read obs velocities as cc.Displacement_Points
    obs_disp_pts = HR.read_all_data_table(exp_dict["data_file"]);   # all 783 points
    # inv_tools.print_typical_uncs(obs_disp_pts);
    obs_disp_pts = correct_for_far_field_terms(exp_dict, obs_disp_pts);  # needed from Fred's work
    # Experimental options:
    if exp_dict["continuous_only"] == 1:
        obs_disp_pts = dpo.utilities.filter_to_meas_type(obs_disp_pts, 'continuous');  # experimental design step
    obs_disp_pts = inv_tools.remove_nearfault_pts(obs_disp_pts, exp_dict["inverse_dir"]+exp_dict["faults"]["Maa"]["points"])
    obs_disp_pts = inv_tools.remove_nearfault_pts(obs_disp_pts, exp_dict["inverse_dir"]+exp_dict["faults"]["BSF"]["points"])
    for excluded_region in exp_dict["exclude_regions"]:
        obs_disp_pts = dpo.utilities.filter_to_exclude_bounding_box(obs_disp_pts, excluded_region);  # Lassen etc.
    obs_disp_pts = dpo.utilities.filter_by_bounding_box(obs_disp_pts, exp_dict["bbox"]);  # north of 38.5

    # INPUT stage: Read GF models based on the configuration parameters
    gf_elements = read_hb_fault_gf_elements(exp_dict);  # list of GF_elements, one for each fault-related column of G.

    # COMPUTE STAGE: PREPARE ROTATION GREENS FUNCTIONS AND LEVELING OFFSET
    gf_elements_rotation = inv_tools.get_GF_rotation_elements(obs_disp_pts);  # 3 elements: rot_x, rot_y, rot_z
    gf_elements = gf_elements + gf_elements_rotation;  # add rotation elements to matrix
    gf_elements_rotation2 = inv_tools.get_GF_rotation_elements(obs_disp_pts, target_region=[-126, -119, 40.4, 46]);
    gf_elements = gf_elements + gf_elements_rotation2;  # add second rotation elements (Oregon Coast Block)
    gf_element_lev = inv_tools.get_GF_leveling_offset_element(obs_disp_pts);  # 1 element: lev reference frame
    gf_elements = gf_elements + gf_element_lev;

    # COMPUTE STAGE: Pairing is necessary in case you've filtered out any observations along the way.
    paired_obs, paired_gf_elements = inv_tools.pair_gf_elements_with_obs(obs_disp_pts, gf_elements);

    inv_tools.visualize_GF_elements(paired_gf_elements, exp_dict["outdir"], exclude_list='all');

    # COMPUTE STAGE: INVERSE.  Reduces certain points to only-horizontal, only-vertical, etc.
    list_of_gf_columns = [];
    for paired_gf in paired_gf_elements:
        G_one_col = inv_tools.buildG_column(paired_gf.disp_points, paired_obs);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);

    # Build observation vector
    obs, sigmas = inv_tools.build_obs_vector(paired_obs);
    sigmas = np.divide(sigmas, np.nanmean(sigmas));  # normalizing so smoothing has same order-of-magnitude
    if exp_dict["unc_weighted"] == 0:
        sigmas = np.ones(np.shape(obs));
    G /= sigmas[:, None];
    weighted_obs = obs / sigmas;

    # Add optional smoothing penalty, overwriting old variables
    if 'smoothing' in exp_dict.keys():
        G, weighted_obs, sigmas = inv_tools.build_smoothing(paired_gf_elements, 'CSZ_dist',
                                                            exp_dict["smoothing"], G, weighted_obs, sigmas);
    # Add optional slip weighting penalty, overwriting old variables
    if 'slip_penalty' in exp_dict.keys():
        G, weighted_obs, sigmas = inv_tools.build_slip_penalty(paired_gf_elements,
                                                               exp_dict["slip_penalty"], G, weighted_obs, sigmas);

    # Money line: Constrained inversion
    lb = [x.lower_bound for x in paired_gf_elements];
    ub = [x.upper_bound for x in paired_gf_elements];
    response = scipy.optimize.lsq_linear(G, weighted_obs, bounds=(lb, ub), max_iter=1500, method='bvls');
    M_opt = response.x;  # parameters of best-fitting model
    print(response.message);
    if response.message == "The maximum number of iterations is exceeded.":
        print("Maximum number of iterations exceeded. Cannot trust this inversion. Exiting");
        sys.exit(0);

    # Make forward predictions.  Work in disp_pts as soon as possible, not matrices.
    M_rot_only, M_no_rot = inv_tools.unpack_model_of_rotation_only(M_opt, [x.fault_name for x in paired_gf_elements]);
    M_csz = inv_tools.unpack_model_of_target_param(M_opt, [x.fault_name for x in paired_gf_elements], 'CSZ_dist');
    M_LSF = inv_tools.unpack_model_of_target_param(M_opt, [x.fault_name for x in paired_gf_elements], 'LSFRev');
    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, paired_obs);
    rot_modeled_pts = inv_tools.forward_disp_points_predictions(G, M_rot_only, sigmas, paired_obs);
    norot_modeled_pts = inv_tools.forward_disp_points_predictions(G, M_no_rot, sigmas, paired_obs);
    csz_modeled_pts = inv_tools.forward_disp_points_predictions(G, M_csz, sigmas, paired_obs);
    lsf_modeled_pts = inv_tools.forward_disp_points_predictions(G, M_LSF, sigmas, paired_obs);

    # Output stage
    fault_dict_lists = [item.fault_dict_list for item in paired_gf_elements];
    rms_mm_h, rms_chi2_h = dpo.compute_rms.obs_vs_model_L2_horiz(paired_obs, model_disp_pts);
    rms_mm_v, rms_chi2_v = dpo.compute_rms.obs_vs_model_L2_vertical(paired_obs, model_disp_pts);
    rms_mm_t, rms_chi2_t = dpo.compute_rms.obs_vs_model_L2_aggregate(paired_obs, model_disp_pts);
    rms_obj = [rms_mm_h, rms_mm_v, rms_mm_t, rms_chi2_h, rms_chi2_v, rms_chi2_t];
    rms_title = "RMS: %f mm/yr" % rms_mm_t;
    print(" ", rms_title);

    residual_pts = dpo.utilities.subtract_disp_points(paired_obs, model_disp_pts);
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, exp_dict["outdir"] + '/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(residual_pts, exp_dict["outdir"] + '/resid_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(paired_obs, exp_dict["outdir"] + '/simple_obs_file.txt');

    dpo_out.write_disp_points_gmt(model_disp_pts, exp_dict["outdir"] + '/model_pred_gmt.txt', write_meas_type=True)
    dpo_out.write_disp_points_gmt(residual_pts, exp_dict["outdir"] + '/resid_file_gmt.txt', write_meas_type=True)
    dpo_out.write_disp_points_gmt(paired_obs, exp_dict["outdir"] + '/simple_obs_file_gmt.txt', write_meas_type=True)
    dpo_out.write_disp_points_gmt(csz_modeled_pts, exp_dict["outdir"] + '/csz_model_pred.txt', write_meas_type=True)
    dpo_out.write_disp_points_gmt(lsf_modeled_pts, exp_dict["outdir"] + '/lsf_model_pred.txt', write_meas_type=True)

    inv_tools.write_model_params(M_opt, rms_mm_t, exp_dict["outdir"] + '/' + exp_dict["model_file"], paired_gf_elements)
    inv_tools.write_summary_params(M_opt, rms_obj, exp_dict["outdir"] + '/model_results_human.txt',
                                   paired_gf_elements, ignore_faults=['CSZ_dist'], message=response.message);
    inv_tools.write_fault_traces(M_opt, paired_gf_elements, exp_dict["outdir"] + '/fault_output.txt',
                                 ignore_faults=['CSZ_dist', 'x_rot', 'y_rot', 'z_rot', 'lev_offset']);
    readers.write_csz_dist_fault_patches(fault_dict_lists, M_opt, exp_dict["outdir"] + '/csz_model.gmt');
    inv_tools.view_full_results(exp_dict, paired_obs, model_disp_pts, residual_pts, rot_modeled_pts,
                                norot_modeled_pts, rms_title, region=[-127, -119.7, 37.7, 43.5]);
    library.plot_fault_slip.map_source_slip_distribution([], exp_dict["outdir"] + "/csz_only_pred.png",
                                                         disp_points=csz_modeled_pts, region=[-127, -119.7, 37.7, 43.5],
                                                         scale_arrow=(1.0, 0.010, "1 cm/yr"), v_labeling_interval=0.001)
    library.plot_fault_slip.plot_data_model_residual(exp_dict["outdir"] + "/results.png", paired_obs,
                                                     model_disp_pts, residual_pts, [-126, -119.7, 37.7, 43.3],
                                                     scale_arrow=(0.5, 0.020, "2 cm"), v_labeling_interval=0.003,
                                                     fault_dict_list=[], rms=rms_mm_t);
    return;


if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Error! Please provide the name of config files. "
              " Ex: humboldt_inversion_driver.py config.json fault_config.json");
        sys.exit(0);
    else:
        run_humboldt_inversion();
