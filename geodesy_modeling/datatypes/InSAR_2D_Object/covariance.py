"""Create covariance matrix routine for 2d insar object."""


from ..InSAR_1D_Object import covariance
import matplotlib.pyplot as plt
import numpy as np


def compute_insar_varigram(insarobj, rmin, rmax, dr, Nmax, rng=None):
    """
    Estimate the semi-varigram structure within an InSAR interferogram, assuming the noise structure is isotropic,
    and has an exponential function with distance.
    Note that if your lat/lon are in degrees, this could cause some distortion in the Cartesian distance metric,
    especially if you're working in the high latitudes where lon and lat do not have similar length scales.

    C(r) = sigma * exp(-r/L)
    From Lohman and Simons, 2005.

    :param insarobj: InSAR_2D object
    :param rmin: inner radius, in the same units as insarobj.lon and insarobj.lat
    :param rmax: outer radius, in the same units as insarobj.lon and insarobj.lat
    :param dr: step-size in radius, in the same units as insarobj.lon and insarobj.lat
    :param Nmax: number of samples, a large integer
    :param rng: specific type of random number generate you wish to use
    :return: radius, covariance at each radius, float sigma, and float L
    """
    insar1d = insarobj.convert_to_insar1D()
    rnew, znew, sigma, L = covariance.compute_insar_varigram(insar1d, rmin, rmax, dr, Nmax, rng=rng)
    return rnew, znew, sigma, L


def plot_cov_results(rnew, znew, sigma, L, outfile):
    """ Plot the covariance vs. distance, both empirical/observed and exponential model fit. """
    print("Plotting file %s " % outfile)
    predicted = np.multiply(sigma, np.exp(-rnew / L))  # Estimating the covariance structure as an exponential
    plt.figure(figsize=(8, 7), dpi=300)
    plt.plot(rnew, znew, '.k', label='Observed')
    plt.plot(rnew, predicted, '-r', label='Predicted, L=%.3f, sigma=%.3f' % (L, sigma))
    plt.xlabel('Distance')
    plt.ylabel('Covariance (mm^2)')
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


def write_covd_parameters(L, sigma, outfile):
    """
    :param L: lengthscale, float
    :param sigma: amplitude of the variance near zero distance, float
    :param outfile: string, filename
    """
    print("Writing covariance exponential parameters in file %s" % outfile)
    with open(outfile, 'w') as ofile:
        ofile.write("# Lengthscale, Sigma\n")
        ofile.write("%f %f" % (L, sigma))
    return


def read_covd_parameters(infile):
    print("Reading file %s " % infile)
    with open(infile) as ifile:
        ifile.readline()
        data = ifile.readline()
        L = float(data.split()[0])
        sigma = float(data.split()[1])
    return L, sigma
