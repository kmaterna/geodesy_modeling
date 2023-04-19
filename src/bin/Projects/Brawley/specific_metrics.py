import numpy as np
from ....MultiTemporalInversion import metrics


# -------- BRAWLEY EXPT MISFIT FUNCTIONS ----------- #
def compute_brawley_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type):
    """ Ignore the western part of the domain, specific to Brawley. """
    idx = np.where(obs_pos[:, 0] > -115.7);  # -116 is pretty global for Brawley
    idx = idx[0]
    obs_type = np.array(obs_type);
    metric = metrics.compute_simple_misfit(obs_pos[idx], obs_disp[idx], pred_disp[idx], obs_sigma[idx], obs_type[idx]);
    return metric;


def brawley_misfit_driver(config, outfile):
    """Compute three metrics for gps, insar, and leveling respectively"""
    print("Calculating metrics for Brawley inversion results.");
    [obs_pos, obs_disp, pred_disp, obs_sigma, obs_type] = metrics.read_obs_vs_predicted_object(config);
    metric = compute_brawley_misfit(obs_pos, obs_disp, pred_disp, obs_sigma, obs_type);
    metrics.write_simple_misfit(metric, outfile);  # matching write function
    return;
