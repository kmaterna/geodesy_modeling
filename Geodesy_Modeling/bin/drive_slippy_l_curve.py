#!/usr/bin/env python
# Run Slippy across multiple choices of parameters, for l-curve analysis

import sys, json, subprocess
from Geodesy_Modeling import MultiTemporalInversion


def welcome_and_parse(argv):
    print("Welcome to the INVERSION L-CURVE.");
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ");
        sys.exit(0);
    else:
        config = argv[1];
    config_file = open(config, 'r');
    config1 = json.load(config_file);
    subprocess.call(['mkdir', '-p', config1["output_dir_lcurve"]], shell=False);
    for i, key in enumerate(config1["faults"].keys()):
        fault_name = config1["faults"][key]["filename"]
        subprocess.call(['cp', fault_name, config1['output_dir']]);  # save fault files, record-keeping
    with open(config1['output_dir_lcurve']+'/config.json', 'w') as fp:
        json.dump(config1, fp, indent="  ");   # save master config file, record-keeping
    return config1;


if __name__ == "__main__":
    config = welcome_and_parse(sys.argv);
    for alpha in config['range_alpha']:
        config["alpha"] = alpha;    # set the alpha
        config["output_dir"] = config["output_dir_lcurve"]+"/alpha_"+str(alpha)+"/";   # set the output dir
        subprocess.call(['mkdir', '-p', config["output_dir"]], shell=False);
        MultiTemporalInversion.buildG.beginning_calc(config);
        MultiTemporalInversion.metrics.main_function(config);
