
import glob, json, os
import numpy as np
from . import l_curve_plots


def glob_and_drive_1d_lcurve(target_dir='.', name_of_printed_config="configs_used.txt", paramname='smoothing',
                             name_of_results_file="model_results_human.txt", misfitname="RMS",
                             outname="smoothing_curve.png", xlabel="Smoothing", corner_point=None):
    """
    Get every smoothing parameter in a directory where smoothing experiment has been run multiple times.

    :param target_dir: directory name
    :param name_of_printed_config: file name
    :param paramname: string, found within config_file
    :param name_of_results_file: file name
    :param misfitname: string, found within results_file
    :param outname: string
    :param xlabel: string
    :param corner_point: float, optional x-location where an annotation will be drawn
    """
    config_files = glob.glob(target_dir+"/**/"+name_of_printed_config);
    results_files = glob.glob(target_dir + "/**/" + name_of_results_file);
    smoothings = read_param_from_list_of_config_files(config_files, paramname);
    misfits = read_misfits_from_list_of_files(results_files, misfitname);
    l_curve_plots.plot_1d_curve(smoothings, misfits, xlabel, outname, corner_point);
    l_curve_plots.write_1d_curve(smoothings, misfits, target_dir+"/lcurve_points.txt");
    return;


def collect_curve_points_slippy(top_level_dir, config_file_name, results_file_name, misfitname):
    """ Harvest parameter values and misfit values from a bunch of l-curve directories. """
    other_dirs = glob.glob(top_level_dir+"/*");
    params_array, misfit_array = [], [];
    for i in other_dirs:
        if os.path.isdir(i):
            print("Reading ", i);
            params = read_alpha_and_penalty_from_file(i+'/'+config_file_name);  # get config params
            resultsfile = i + "/"+results_file_name;
            misfit = read_misfit_value_from_text_file(resultsfile, misfitname)  # get results
            params_array.append(params); misfit_array.append(misfit);
    return [params_array, misfit_array];


def read_alpha_and_penalty_from_file(config_filename):
    """ Read json parameters smoothing (lamda) and slip penalty (alpha) from a config file, Slippy format. """
    alpha = read_param_from_config_file(config_filename, paramname='alpha');
    penalty = read_fault_param_from_config_file(config_filename, fault_name='fault1', paramname='penalty');
    return [alpha, penalty];


def read_param_from_config_file(config_filename, paramname):
    """Simplest dictionary reader for one parameter, paramname."""
    config_file = open(config_filename, 'r');
    tempdict = json.load(config_file)
    retval = tempdict[paramname];
    config_file.close();
    return retval;


def read_value_from_results_file(filename, paramname):
    """Assumes the line has a certain format, like Maacama: 10.3 """
    with open(filename, 'r') as fp:
        for line in fp:
            if paramname in line:
                value = float(line.split()[1]);  # slip rate value
    return value;


def read_fault_param_from_config_file(config_filename, fault_name='fault1', paramname='penalty'):
    """Simplest dictionary reader for one parameter, paramname, associated with a fault, fault_name."""
    config_file = open(config_filename, 'r');
    tempdict = json.load(config_file)
    faults = tempdict["faults"];
    retval = faults[fault_name][paramname];
    config_file.close();
    return retval;


def read_misfit_value_from_text_file(filename, fieldname):
    """Assumes the line has a certain format. Like: Normalized misfit: 11.851312 sigma """
    misfit = np.nan;
    with open(filename, 'r') as fp:
        for line in fp:  # if 'Average misfit' in line:
            if fieldname in line:  # Have found that misfit should consider sigmas
                misfit = float(line.split()[-2]);  # misfit metric somehow
    return misfit;


def read_param_from_list_of_config_files(filelist, paramname):
    """
    Read list of parameter value from a list of config files. Ex: a series of smoothing values from experiments.
    """
    return [read_param_from_config_file(file, paramname) for file in filelist];


def read_fault_param_from_list_of_config_files(filelist, fault_name, param_name):
    """
    Read list of fault-specific parameter value from a list of config files. Ex: a series of fault minimum slip rates.
    """
    return [read_fault_param_from_config_file(file, fault_name, param_name) for file in filelist];


def read_misfits_from_list_of_files(filelist, paramname):
    """
    Read list of misfit values from a list of human-readable results files.
    """
    return [read_misfit_value_from_text_file(ifile, paramname) for ifile in filelist];


def read_values_from_list_of_files(filelist, paramname):
    """
    Read parameter values from a list of human-readable results files. Ex: a series of experiments' slip rate values.
    """
    return [read_value_from_results_file(ifile, paramname) for ifile in filelist];
