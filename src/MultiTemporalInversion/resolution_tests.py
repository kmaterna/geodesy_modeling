# Tools for resolution tests on inversions
# Option 1: View model resolution matrix, R
# Option 2: how much displacement is caused by a unit displacement at each model cell?
# Option 3: invert 100 iterations of the data, and take the standard deviation of that distribution
#     Can only do this once we have defined G, obs, etc.
#     This would involve a somewhat refactor of this script (configure, input, compute, output)
#     Or would do for only a simple inversion, slippy-style, not a compound inversion

import numpy as np
import matplotlib.pyplot as plt
import scipy
import slippy.basis

def analyze_model_resolution_matrix(G, num_obs, outdir):
    """
    Simply analyze the resolution matrix R (Menke, 1989)
    Before leveling offsets have been added.
    """
    # Show big-G matrix for all times, all data
    plt.figure(figsize=(12, 8), dpi=300);
    plt.imshow(G, vmin=-0.02, vmax=0.02, aspect=1);
    plt.savefig(outdir + "/G_resolution.png");

    # Model Resolution Matrix
    Ggi = scipy.linalg.pinv(G);  # Replace all smoothing with zeros for second half of this expression
    G_nosmoothing = G.copy();
    G_nosmoothing[num_obs:, :] = 0;
    Rmatrix = np.dot(Ggi, G_nosmoothing);  # R matrix

    # Model Covariance values
    covb_slip = np.dot(Ggi, Ggi.T);
    sig_slipb = np.sqrt(np.diag(covb_slip));

    # Viewing the diagonal elements of R (might be helpful?)
    plt.figure();
    plt.plot(np.diag(Rmatrix));
    plt.savefig(outdir + '/model_resolution_diagonal.png');

    # Viewing the total picture of R: shows the model resolution along the diagonal.
    plt.figure(figsize=(12, 8), dpi=300);
    plt.imshow(Rmatrix, aspect=1);
    plt.colorbar();
    plt.savefig(outdir + "/Rmatrix.png");
    return np.diag(Rmatrix), sig_slipb;


def empirical_slip_resolution(G, total_fault_slip_basis):
    """
    A relative measure of model resolution for strike slip and dip slip motion on faults.
    Analyze model resolution by putting unit slip on each fault model parameter, calculating average geodetic response.
    This function does not take into account observation uncertainties or smoothing parameters
    """
    num_model_params = np.shape(G)[1];  # number of model parameters
    resolution_vector = np.zeros((num_model_params, 1));
    if len(total_fault_slip_basis[0]) != 2:
        print("Error, I don't know how to do empirical uncertainties for 1 or 3 basis vectors yet.");
        return resolution_vector;
    # Reshaping the slip basis into a long list of vectors, one for each slip component
    for i in range(len(total_fault_slip_basis)):   # for each fault patch in the model
        mhat = np.zeros((num_model_params, 1));   # create a model vector of zeros
        angle_with_reverse1 = np.arctan(total_fault_slip_basis[0][0][0] / total_fault_slip_basis[0][0][1]);
        angle_with_reverse2 = np.arctan(total_fault_slip_basis[0][1][0] / total_fault_slip_basis[0][1][1]);
        transform_matrix = np.array([[np.sin(angle_with_reverse1), -np.sin(angle_with_reverse1)],
                                     [-np.cos(angle_with_reverse2), -np.cos(angle_with_reverse2)]]);  # transform
        cardinal_strike_slip = np.dot(np.linalg.inv(transform_matrix), np.array([1, 0]));
        cardinal_dip_slip = np.dot(np.linalg.inv(transform_matrix), np.array([0, -1]));  # coordinate transform
        mhat[2*i] = cardinal_strike_slip[0];    # set particular model parameters to 1m strike slip
        mhat[2*i + 1] = cardinal_strike_slip[1];
        dhat = np.dot(G, mhat);    # predicted displacement from unit slip
        resolution_vector[2*i] = np.nanmean(np.abs(dhat));   # measure a metric from unit strike slip
        mhat[2*i] = cardinal_dip_slip[0];    # set particular model parameters to 1m dip slip
        mhat[2*i + 1] = cardinal_dip_slip[1];
        dhat = np.dot(G, mhat);    # predicted displacement from unit slip
        resolution_vector[2*i + 1] = np.nanmean(np.abs(dhat));   # measure a metric from unit dip slip
    return resolution_vector;


def parse_empirical_res_outputs(res_f, Ns_total, Ds,  num_lev_offsets):
    """
    Rotate the resolution into dip slip and strike slip components.
    res_f: resolution vector
    Ns_total: number of segments
    Ds: number of dimensions in the basis
    total_fault_slip_basis:
    num_lev_offsets: the last n model parameters are leveling offsets that don't get plotted on faults
    """
    n_params = int(len(res_f) - num_lev_offsets);  # num fault params
    res = res_f[0:n_params].reshape((Ns_total, Ds))  # ASSUMES ALL FAULTS HAVE SAME NUMBER OF BASIS VECTORS
    cardinal_res = np.hstack((res, np.zeros((len(res), 1))));  # adding the "tensile"
    return cardinal_res;


def parse_checkerboard_res_outputs(res_f, Ns_total, Ds, total_fault_slip_basis, num_lev_offsets):
    """
    Rotate the slip into dip slip and strike slip components.
    res_f: resolution vector
    Ns_total: number of segments
    Ds: number of dimensions in the basis
    total_fault_slip_basis:
    num_lev_offsets: the last n model parameters are leveling offsets that don't get plotted on faults
    """
    n_params = int(len(res_f) - num_lev_offsets);  # num fault params
    res_f = res_f[0:n_params].reshape((Ns_total, Ds))  # ASSUMES ALL FAULTS HAVE SAME NUMBER OF BASIS VECTORS
    cardinal_res = slippy.basis.cardinal_components(res_f, total_fault_slip_basis)
    return cardinal_res;


def bootstrapped_model_resolution(_G_total, _G_nosmooth, _d, _sig, _weights):
    """
    An absolute measure (in mm) of model resolution on faults
    Run the model a hundred times with random noise realizations. Get the noise floor.
    We probably need G_total to have smoothing here, while G_nosmooth does not.
    Useful lines:
    slip_f = reg_nnls(G_total, d_total)   # slip_f is the model
    pred_disp_f = G_nosmooth_total.dot(slip_f) * sig_total * weight_total;   # the forward prediction
    I'm now thinking this is an interesting problem - not as useful as I thought
    This depends on the input data
    """
    return;


def checkerboard_test(patches_f, Ds, num_leveling_params, num_width, fault_num_array, checker_width=3, fault_num=0):
    """
    Basic checkerboard utility to build checkers on a single fault. Making a checkerboard input pattern.
    """

    checkerboard_vector = np.zeros((len(patches_f)+num_leveling_params,));
    is_target_patch = np.where(np.array(fault_num_array) == fault_num);
    num_patches = len(is_target_patch[0]);   # ex: 169
    patch_vector = np.zeros((num_patches,));
    start_patch_idx = is_target_patch[0][0];  # ex: 0

    # Define the rows for checker-0 and checker-1
    type0_slice = np.zeros((num_width,));
    type1_slice = np.zeros((num_width,));
    for i in range(len(type0_slice)):
        if np.mod(i, checker_width*2) < checker_width:
            type0_slice[i] = 0.4;
    for i in range(len(type1_slice)):
        if np.mod(i, checker_width*2) >= checker_width:
            type1_slice[i] = 0.4;

    # Fill in the rows for checker-0 and checker-1
    for i in np.arange(len(patch_vector), step=num_width):
        if np.mod(i, num_width*checker_width*2) < num_width*checker_width:
            patch_vector[i:i+num_width] = type0_slice;
        else:
            patch_vector[i:i + num_width] = type1_slice;

    for i in range(len(patch_vector)):
        checkerboard_vector[start_patch_idx + (i*Ds)] = patch_vector[i];
        checkerboard_vector[start_patch_idx + (i*Ds) + 1] = patch_vector[i];

    return checkerboard_vector;
