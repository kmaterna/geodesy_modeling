#!/usr/bin/env python

"""
A generic driver for multiple projects
Run Slippy across multiple choices of parameters, for l-curve analysis
"""

import sys, json, os, shutil
from Geodesy_Modeling.src import MultiTemporalInversion
from Geodesy_Modeling.src.Inversion import l_curve, l_curve_plots


def welcome_and_parse(argv):
    print("Welcome to the INVERSION L-CURVE.")
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ")
        sys.exit(0)
    else:
        config = argv[1]
    config_file = open(config, 'r')
    config1 = json.load(config_file)
    os.makedirs(config1['output_dir_lcurve'], exist_ok=True)
    for i, key in enumerate(config1["faults"].keys()):
        fault_file_name = config1["faults"][key]["filename"]
        fault_name = os.path.split(fault_file_name)[1]
        shutil.copy2(fault_file_name, os.path.join(config1['output_dir_lcurve'], fault_name))  # save fault files
    if 'G' not in config1.keys():
        config1['G'] = 30e9  # default shear modulus is 30 GPa
    if 'resolution_test' not in config1.keys():
        config1['resolution_test'] = ''   # default resolution test is none
    output_json = os.path.join(config1['output_dir_lcurve'], 'config.json')
    with open(output_json, 'w') as fp:
        json.dump(config1, fp, indent="  ")   # save master config file, record-keeping
    return config1

def iterate_many_inversions(config):
    """
    A driver for looping multiple inversions depending on the experiment, testing the impact of alpha or smoothing.
    """
    if not config["switch_alpha"] and not config["switch_penalty"]:   # no search at all.
        print("Check your configuration. Not searching through alpha or lambda.")
        sys.exit(0)
    elif config["switch_alpha"] and not config["switch_penalty"]:   # 1d search in slip penalty
        for alpha in config['range_alpha']:
            config["alpha"] = alpha    # set the alpha
            config["output_dir"] = os.path.join(config["output_dir_lcurve"], "alpha_"+str(alpha), "")
            os.makedirs(config['output_dir'], exist_ok=True)  # set output dir
            MultiTemporalInversion.buildG.beginning_calc(config)  # the guts of drive_slippy_multitemporal.py
            MultiTemporalInversion.metrics.main_function(config)
    elif config["switch_penalty"] and not config["switch_alpha"]:   # 1d search in smoothing penalty
        for penalty in config['range_penalty']:
            for key in config["faults"].keys():
                config["faults"][key]["penalty"] = penalty    # set the smoothing penalty
            config["output_dir"] = os.path.join(config['output_dir_lcurve'], "penalty_"+str(penalty), "")
            os.makedirs(config['output_dir'], exist_ok=True)  # set output dir
            MultiTemporalInversion.buildG.beginning_calc(config)
            MultiTemporalInversion.metrics.main_function(config)
    else:  # 2d search for smoothing and slip penalty
        for alpha in config['range_alpha']:
            for penalty in config['range_penalty']:
                for key in config["faults"].keys():
                    config["faults"][key]["penalty"] = penalty  # set the smoothing penalty
                config["alpha"] = alpha  # set the alpha
                config["output_dir"] = os.path.join(config["output_dir_lcurve"],
                                                    "alpha_"+str(alpha)+"_"+str(penalty), "")
                os.makedirs(config['output_dir'], exist_ok=True)
                MultiTemporalInversion.buildG.beginning_calc(config)
                MultiTemporalInversion.metrics.main_function(config)
    return


def make_lcurve_driver(config):
    """
    A driver for creating an L-curve from multiple runs of an inversion, with config syntax matching slippy.metrics.
    """
    [params, misfits] = l_curve.collect_curve_points_slippy(top_level_dir=config["output_dir_lcurve"],
                                                            config_file_name='config.json',
                                                            results_file_name="summary_stats_simple.txt",
                                                            misfitname='Average normalized misfit')
    l_curve_plots.plot_l_curve_coordinator(params, misfits, config["output_dir_lcurve"])
    return


if __name__ == "__main__":
    config = welcome_and_parse(sys.argv)
    iterate_many_inversions(config)
    make_lcurve_driver(config)
