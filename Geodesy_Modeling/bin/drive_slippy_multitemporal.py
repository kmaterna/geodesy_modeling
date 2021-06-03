#!/usr/bin/env python
# Run Slippy with multiple time intervals in the input
# (A really big G matrix)

import sys, json, subprocess
from Geodesy_Modeling import MultiTemporalInversion


def welcome_and_parse(argv):
    print("Welcome to the MULTITEMPORAL INVERSION.");
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ");
        sys.exit(0);
    else:
        config = argv[1];
    config_file = open(config, 'r');
    config1 = json.load(config_file);
    subprocess.call(['mkdir', '-p', config1["output_dir"]], shell=False);
    for i, key in enumerate(config1["faults"].keys()):
        fault_name = config1["faults"][key]["filename"]
        subprocess.call(['cp', fault_name, config1['output_dir']]);  # save fault files, record-keeping
    with open(config1['output_dir']+'/config.json', 'w') as fp:
        json.dump(config1, fp, indent="  ");   # save config file, record-keeping
    return config1;


if __name__ == "__main__":
    config = welcome_and_parse(sys.argv);
    MultiTemporalInversion.buildG.beginning_calc(config);
