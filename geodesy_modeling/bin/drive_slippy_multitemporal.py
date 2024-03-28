#!/usr/bin/env python

"""
Generic Driver for multiple projects
Run Slippy with multiple time intervals in the input (a really big G matrix)
"""

import sys
import json
import os
import shutil
from geodesy_modeling import MultiTemporalInversion


def welcome_and_parse(argv):
    print("Welcome to the MULTITEMPORAL INVERSION.")
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ")
        sys.exit(0)
    else:
        config = argv[1]
    config_file = open(config, 'r')
    config1 = json.load(config_file)
    os.makedirs(config1['output_dir'], exist_ok=True)
    for i, key in enumerate(config1["faults"].keys()):
        fault_file_name = config1["faults"][key]["filename"]
        fault_name = os.path.split(fault_file_name)[1]
        shutil.copy2(fault_file_name, os.path.join(config1['output_dir'], fault_name))  # save faults, record-keeping
    if 'G' not in config1.keys():
        config1['G'] = 30e9  # default shear modulus is 30 GPa
    return config1


if __name__ == "__main__":
    config = welcome_and_parse(sys.argv)
    MultiTemporalInversion.buildG.beginning_calc(config)
    MultiTemporalInversion.metrics.main_function(config)
