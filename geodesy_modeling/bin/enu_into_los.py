#!/usr/bin/env python

"""
Project 3 ENU grids into a certain LOS. Command line tool.
"""

import json
import sys
from geodesy_modeling.datatypes.InSAR_2D_Object import model_enu_grids_into_los


def parse_config(argv):
    if len(argv) < 2:
        print("Error! Please provide the name of a config json. Exiting. ")
        sys.exit(0)
    else:
        config = argv[1]
    config_file = open(config, 'r')
    config1 = json.load(config_file)
    return config1


if __name__ == "__main__":
    my_params = parse_config(sys.argv)  # a command line application
    model_enu_grids_into_los.do_synthetic_grid_LOS(my_params)
