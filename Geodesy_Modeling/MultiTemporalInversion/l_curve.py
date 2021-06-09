# Tools for L-curve analysis.

import matplotlib.pyplot as plt
import json
import glob, os


def collect_curve_points(config):
    """ Harvest the parameter values and misfit values from a bunch of l-curve directories """
    exp_dir = config["output_dir_lcurve"];
    other_dirs = glob.glob(exp_dir+"/*");
    params_array, misfit_array = [], [];
    for i in other_dirs:
        if os.path.isdir(i):
            print("Reading ", i);
            params = read_params_from_dir(i);
            misfit = read_misfits_from_dir(i);
            params_array.append(params);
            misfit_array.append(misfit);
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


def read_misfits_from_dir(dirname):
    """
    Utility Function: Read misfits from a results directory.
    Depends on the structure of the misfit file.
    """
    ifile = dirname + "/summary_stats.txt";
    with open(ifile, 'r') as fp:
        for line in fp:
            if 'Average misfit' in line:
                misfit = float(line.split()[2]);  # misfit metric somehow
    return misfit;


def plot_l_curve(params, misfits, outfile):
    """ Coordiantor function for driving l-curve plots """
    all_alphas = [x[0] for x in params];
    all_penalties = [x[1] for x in params];
    if len(set(all_alphas)) == 1:
        plot_1d_curve(all_penalties, misfits, 'Smoothing Penalty', outfile);
    elif len(set(all_penalties)) == 1:
        plot_1d_curve(all_alphas, misfits, 'Slip Penalty, alpha', outfile);
    else:
        plot_2d_curve(all_alphas, all_penalties, misfits, 'alpha', 'smoothing', outfile);
    return;


def plot_1d_curve(param, misfit, axis_name, outfile):
    """ Making 1D plot for L-curve """
    plt.figure();
    plt.plot(param, misfit, '.', markersize=14);
    plt.xlabel(axis_name, fontsize=14);
    plt.ylabel('Misfit (mm)', fontsize=14);
    plt.gca().invert_xaxis();
    plt.savefig(outfile);
    return;


def plot_2d_curve(alphas, penalties, misfits, param1_name, param2_name, outfile):
    """ Making 2D surface plot for L-curve """
    plt.figure();
    plt.savefig(outfile);
    return;


def main_driver(config):
    [params, misfits] = collect_curve_points(config);
    plot_l_curve(params, misfits, config["output_dir_lcurve"]+"/l_curve.png");
    return;
