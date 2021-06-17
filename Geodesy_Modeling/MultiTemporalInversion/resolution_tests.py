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

def analyze_model_resolution_matrix(G, outdir):
    """ Analyze the resolution matrix R (Menke, 1989) """
    Ggi = scipy.linalg.pinv(G);
    Rmatrix = np.dot(Ggi, G);
    for i in range(np.shape(Rmatrix)[0]):
        for j in range(np.shape(Rmatrix)[1]):
            if abs(Rmatrix[i][j]) < 1e-10:  # cleaning up the really small entries
                Rmatrix[i][j] = np.nan;

    # Viewing the diagonal elements of R (might be helpful?)
    plt.figure();
    plt.plot(np.diag(Rmatrix));
    plt.savefig(outdir + '/model_resolution_diagonal.png');

    # Viewing the total picture of R: shows the model resolution along the diagonal.
    # Which ones are rows and columns? Assuming we go from shallow to deep.
    plt.figure(figsize=(12, 8), dpi=300);
    plt.imshow(np.log10(Rmatrix), aspect=1);
    plt.colorbar();
    plt.savefig(outdir + "/model_resolution.png");
    return;


def empirical_slip_resolution(G, m):
    """
    Analyze model resolution by putting unit slip on each model parameter and calculating average response
    This function does not take into account observation uncertainties.
    """


    return;
