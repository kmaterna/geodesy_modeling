# Tools for L-curve analysis.

import matplotlib.pyplot as plt
import json
import glob, os
import numpy as np


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
        plot_1d_curve(all_penalties, misfits, 'Smoothing Penalty', outfile.split('.')[0]+"_smoothing.png");
    elif len(set(all_penalties)) == 1:
        plot_1d_curve(all_alphas, misfits, 'Slip Penalty, alpha', outfile.split('.')[0]+"_slip.png");
    else:
        plot_2d_curve(all_alphas, all_penalties, misfits, '1/alpha (slip)', '1/smoothing (smoothing)',
                      outfile.split('.')[0]+"_2d.png");
    return;


def plot_1d_curve(param, misfit, axis_name, outfile):
    """ Making 1D plot for L-curve """
    # Linear scale
    plt.figure(figsize=(8, 8), dpi=300);
    plt.plot(param, misfit, '.', markersize=14);
    plt.xlabel(axis_name, fontsize=14);
    plt.ylabel('Misfit (mm)', fontsize=14);
    plt.gca().invert_xaxis();
    plt.savefig(outfile);
    # Log scale
    plt.figure(figsize=(8, 8), dpi=300);
    plt.plot(param, misfit, '.', markersize=14);
    plt.gca().set_yscale('log');
    plt.gca().set_xscale('log');
    plt.xlabel(axis_name, fontsize=14);
    plt.ylabel('Log Misfit (mm)', fontsize=14);
    plt.gca().invert_xaxis();
    plt.savefig(outfile.split('.')[0]+"_log.png");
    return;


def chosen_axis_annotations(ax):
    # Manual annotations for certain projects
    # Example: Heber base inversion: alpha 2, penalty 0.5.
    # WARNING: This is mixing code and data.
    chosen_alpha = 2;        # HEBER
    chosen_penalty = 0.5;    # HEBER
    alpha_range = ax.get_xlim();
    penalty_range = ax.get_ylim();
    [bottom, _] = ax.get_zlim();
    ax.plot([1/chosen_alpha, 1/chosen_alpha], penalty_range, [bottom, bottom], color='black');
    ax.plot(alpha_range, [1/chosen_penalty, 1/chosen_penalty], [bottom, bottom], color='black');
    ax.plot([1/chosen_alpha], [1/chosen_penalty], [bottom], marker='.', color='red', markersize=15);
    return ax;


def plot_2d_curve(alphas, penalties, misfits, param1_name, param2_name, outfile):
    """ Making 2D surface plot for L-curve """
    fig = plt.figure(figsize=(4, 4), dpi=200);
    ax = fig.add_subplot(111, projection='3d')
    alphas = [1/x for x in alphas];
    penalties = [1/x for x in penalties];
    misfits = [np.log(x) for x in misfits];
    ax.scatter(alphas, penalties, misfits, '.', c=misfits, cmap='viridis');
    ax.set_xlabel(param1_name);
    ax.set_ylabel(param2_name);
    ax.set_zlabel('Misfit (mm)')
    ax = chosen_axis_annotations(ax);
    ax.set_title("Misfit vs Regularization Parameters");
    plt.show();  # good for playing.
    # plt.savefig(outfile);
    return;


def main_driver(config):
    [params, misfits] = collect_curve_points(config);
    plot_l_curve(params, misfits, config["output_dir_lcurve"]+"/l_curve.png");
    return;
