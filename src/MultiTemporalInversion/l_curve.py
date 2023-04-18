# Tools for L-curve analysis.

import json, glob, os
from Geodesy_Modeling.src.Inversion.l_curve_plots import plot_l_curve_coordinator
from Geodesy_Modeling.src.Inversion.post_inversion_tools import read_misfits_from_list_of_files


def collect_curve_points(config):
    """ Harvest the parameter values and misfit values from a bunch of l-curve directories """
    exp_dir = config["output_dir_lcurve"];
    other_dirs = glob.glob(exp_dir+"/*");
    params_array, misfit_array = [], [];
    for i in other_dirs:
        if os.path.isdir(i):
            print("Reading ", i);
            params = read_params_from_dir(i);
            resultsfile = i + "/summary_stats_simple.txt";
            misfit = read_misfits_from_list_of_files([resultsfile], 'Average normalized misfit')[0]
            params_array.append(params); misfit_array.append(misfit);
    return [params_array, misfit_array];


def read_params_from_dir(dirname):
    """ Utility Function: Read json parameters like smoothing and penalties from a results directory """
    config_file = dirname + "/config.json"
    config_fh = open(config_file, 'r');
    config_object = json.load(config_fh);
    alpha = config_object['alpha'];
    faults = list(config_object["faults"].keys())
    penalty = config_object["faults"][faults[0]]["penalty"];
    return [alpha, penalty];


def main_driver(config):
    [params, misfits] = collect_curve_points(config);
    plot_l_curve_coordinator(params, misfits, config["output_dir_lcurve"] + "/l_curve.png");
    return;
