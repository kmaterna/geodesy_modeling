# After you've done an inversion, what are the results?
# How much moment is created, and how big is the misfit?
# MultiT and SingleT
# Writes into a summary text file, and prints to screen. 
# Useful for L-curve analysis.

import numpy as np
from Tectonic_Utils.seismo import moment_calculations
import slippy.io


# -------- READ FUNCTIONS ----------- #
def read_obs_vs_predicted_object(config):
    """Read data and model prediction files. There may be many files. """
    obs_pos_column = np.zeros((0, 3));  # no predicted pos because it's the same as obs
    obs_disp_column = np.zeros((0,));
    pred_disp_column = np.zeros((0,));
    obs_sigma_column = np.zeros((0,));  # no predicted sigma since it's 0 for a model
    obs_type_column = [];
    for data_category in config["data_files"].keys():
        obs_file = config["data_files"][data_category]["data_file"];  # infile expected of all data files
        pred_file = config["output_dir"]+config["data_files"][data_category]["outfile"];  # predicted outfile
        data_type = config["data_files"][data_category]["type"];  # type expected of all data files
        print("Reading data and model prediction from %s: %s, %s" % (data_category, obs_file, pred_file))
        if data_type == 'gps':
            pos_geodetic, disp, sigma = slippy.io.read_gps_data(obs_file);  # format: [llh], [enu], [sigs3]
            _, pred_disp, _ = slippy.io.read_gps_data(pred_file);
            obs_disp_fi = disp.reshape((len(pos_geodetic) * 3,))   # reshaping 3-component data into 1d vector
            pred_disp_fi = pred_disp.reshape((len(pos_geodetic) * 3,))
            sigma_fi = sigma.reshape((len(pos_geodetic) * 3,))
            obs_pos_column = np.concatenate((obs_pos_column, pos_geodetic));
            obs_disp_column = np.concatenate((obs_disp_column, obs_disp_fi));
            pred_disp_column = np.concatenate((pred_disp_column, pred_disp_fi));
            obs_sigma_column = np.concatenate((obs_sigma_column, sigma_fi));
            obs_type_column += ["gps"] * len(pos_geodetic);
        else:
            pos_geodetic, disp, sigma, _ = slippy.io.read_insar_data(obs_file);   # ignored elements are basis vectors
            _, pred_disp, _, _ = slippy.io.read_insar_data(pred_file);
            obs_pos_column = np.concatenate((obs_pos_column, pos_geodetic));
            obs_disp_column = np.concatenate((obs_disp_column, disp));
            pred_disp_column = np.concatenate((pred_disp_column, pred_disp));
            obs_sigma_column = np.concatenate((obs_sigma_column, sigma));
            if data_type == 'insar':
                obs_type_column += ["insar"] * len(pos_geodetic);
            else:
                obs_type_column += ["leveling"] * len(pos_geodetic);
    return [obs_pos_column, obs_disp_column, pred_disp_column, obs_sigma_column, obs_type_column];


# -------- SIMPLE MISFIT FUNCTIONS ----------- #
def simple_misfit_driver(config):
    # Compute simple misfit - one column for each data.
    [obs_pos, obs_disp, pred_disp, obs_sigma, obs_type] = read_obs_vs_predicted_object(config);
    metrics = compute_simple_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type);
    write_simple_misfit(metrics, config["summary_file"]);
    return;


def compute_simple_misfit(_obs_pos, obs_disp, pred_disp, obs_sigma, obs_type, data_type='all'):
    """
    The simplest misfit calculation, from any or all types of data
    Options for data_type: ['all', 'gps', 'insar', 'leveling'];
    Want in both absolute numbers and relative to the respective uncertainties.
    """
    if data_type == 'all':
        idx = [index for index, element in enumerate(obs_type) if element != ''];  # take everything in obs vector
    else:
        idx = [index for index, element in enumerate(obs_type) if element == data_type];  # filter on data type
    npts = len(obs_disp[idx]);

    # The L1 norm
    # abs_misfit = np.abs(obs_disp[idx]-pred_disp[idx]);
    # norm_misfit = np.divide(abs_misfit, obs_sigma[idx]);  # divide by sigma

    # The L2 norm: RMS and Chi Squared Values
    sq_misfit = np.square(abs(obs_disp[idx] - pred_disp[idx]));
    norm_misfit = np.divide(sq_misfit, np.square(obs_sigma[idx]));
    rms_misfit = np.sqrt(np.nanmean(sq_misfit));
    chisquared = np.sqrt(np.nanmean(norm_misfit));

    return [rms_misfit, chisquared, npts];


def write_simple_misfit(metrics, outfile):
    ofile = open(outfile, 'w');  # cleaning the file from last times
    print("Writing %s " % outfile);
    print("Average total misfit: %f mm" % (1000 * metrics[0]));
    print("Average normalized misfit: %f sigma \n" % (metrics[1]));
    ofile.write("Average misfit: %f mm\n" % (1000 * metrics[0]));
    ofile.write("Average normalized misfit: %f sigma \n" % (metrics[1]));
    ofile.write("Total npts: %d \n" % metrics[2]);
    ofile.close();
    return;


# -------- COMPOUND MISFIT FUNCTIONS ----------- #
def compound_misfit_driver(config):
    # Compute three metrics for gps, insar, and leveling respectively
    print("Calculating metrics for inversion results.");
    [obs_pos, obs_disp, pred_disp, obs_sigma, obs_type] = read_obs_vs_predicted_object(config);
    metrics = compute_compound_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type);
    write_compound_misfit(metrics, config["summary_file"]);  # matching write function
    return;


def compute_compound_misfit(_obs_pos, obs_disp, pred_disp, obs_sigma, obs_type):
    [gps_misfit, gps_x2, gps_npts] = compute_simple_misfit(_obs_pos, obs_disp, pred_disp, obs_sigma, obs_type,
                                                           data_type='gps');
    [ins_misfit, ins_x2, ins_npts] = compute_simple_misfit(_obs_pos, obs_disp, pred_disp, obs_sigma, obs_type,
                                                           data_type='insar');
    [lev_misfit, lev_x2, lev_npts] = compute_simple_misfit(_obs_pos, obs_disp, pred_disp, obs_sigma, obs_type,
                                                           data_type='leveling');
    return [gps_misfit, gps_x2, gps_npts, ins_misfit, ins_x2, ins_npts, lev_misfit, lev_x2, lev_npts];


def write_compound_misfit(metrics, outfile):
    ofile = open(outfile, 'w');  # cleaning the file from last times
    print("Writing %s " % outfile);

    if metrics[0]:
        print("Average GPS misfit: %f mm" % (1000*metrics[0]) );
        print("Average normalized GPS misfit: %f sigma \n" % (metrics[1]) );
        ofile.write("Average GPS misfit: %f mm\n" % (1000*metrics[0]) );
        ofile.write("Average normalized GPS misfit: %f sigma \n" % (metrics[1]) );
        ofile.write("GPS npts: %d \n" % metrics[2]);

    if metrics[3]:
        print("Average InSAR misfit: %f mm" % (1000*metrics[3]) );
        print("Average normalized InSAR misfit: %f sigma \n" % (metrics[4]) );
        ofile.write("Average InSAR misfit: %f mm\n" % (1000*metrics[3]) );
        ofile.write("Average normalized InSAR misfit: %f sigma \n" % (metrics[4]) );
        ofile.write("InSAR npts: %d \n" % metrics[5]);

    if metrics[6]:
        print("Average Leveling misfit: %f mm" % (1000*metrics[6]) );
        print("Average normalized Leveling misfit: %f sigma \n" % (metrics[7]) );
        ofile.write("Average Leveling misfit: %f mm\n" % (1000*metrics[6]) );
        ofile.write("Average normalized Leveling misfit: %f sigma \n" % (metrics[7]) );
        ofile.write("Leveling npts: %d \n" % metrics[8]);

    ofile.close();
    return;


# -------- COMPOUND MISFIT FUNCTIONS ----------- #
def brawley_misfit_driver(config):
    # Compute three metrics for gps, insar, and leveling respectively
    print("Calculating metrics for inversion results.");
    [obs_pos, obs_disp, pred_disp, obs_sigma, obs_type] = read_obs_vs_predicted_object(config);
    metrics = compute_brawley_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type);
    write_simple_misfit(metrics, config["summary_file"]);  # matching write function
    return;

def compute_brawley_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type):
    """ Ignore the western part of the domain, specific to Brawley. """
    idx = np.where(obs_pos[:, 0] > -115.7);
    idx = idx[0]
    obs_type = np.array(obs_type);
    metrics = compute_simple_misfit(obs_pos[idx], obs_disp[idx], pred_disp[idx], obs_sigma[idx], obs_type[idx]);
    return metrics;


# -------- SLIP COMPUTE FUNCTIONS ----------- #
def slip_metrics_driver(config):
    moments = get_slip_moments(config["output_dir"]+config['epochs']['EpochA']["slip_output_file"], config["G"]);  # A
    write_slip_moments(moments, config["G"], config["summary_file"]);
    return;


def get_slip_moments(slip_filename, mu=30e9):
    # From the inversion results, what is the moment of the slip distribution?
    # Later, it might be important to describe how much of this slip is strike-slip vs reverse.
    # Question for a later time.
    # mu is shear modulus, in Pa.
    moment_total = 0;
    length, width, leftlat, thrust, _ = np.loadtxt(slip_filename, skiprows=1, unpack=True, usecols=(5, 6, 7, 8, 9));
    for i in range(len(length)):
        slip = np.sqrt(leftlat[i]*leftlat[i] + thrust[i]*thrust[i]);
        area = length[i]*width[i];  # m^2
        momenti = moment_calculations.moment_from_muad(mu, area, slip);
        moment_total = moment_total+momenti;
    mw = moment_calculations.mw_from_moment(moment_total);
    print("Calculating moment from %s" % slip_filename);
    return [moment_total, mw];


def write_slip_moments(moments, G, filename):
    # moments
    ofile = open(filename, 'a');  # appending to the file
    scinot = "{:e}".format(moments[0]);
    ofile.write("G = %f GPa\n" % (G/1e9));
    print("Total Slip Moment is %s N-m, equivalent to mw=%f \n" % (scinot, moments[1]) );
    ofile.write("Total Slip Moment is %s N-m, equivalent to mw=%f \n" % (scinot, moments[1]) );
    ofile.close();
    return;


# -------- ACCESS FUNCTIONS ----------- # 
def main_function(config):
    """Misfit can be defined in many ways. Here we point to them. """
    config["summary_file"] = config["output_dir"]+"summary_stats.txt";  # Creates an output file
    brawley_misfit_driver(config);
    slip_metrics_driver(config);   # for the amount of slip.
    return;
