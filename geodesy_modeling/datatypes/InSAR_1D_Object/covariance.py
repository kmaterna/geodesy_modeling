"""
function [rnew,znew,sigma,L]=get_insar_varigram(xx,yy,zz,rmin,rmax,dr,Nmax);
%  get the insar semi-varigram assuming the noise structure is isotropic,
%
%  and fit it as exponential function
%
%  Usage: [rout,zout,sigma,L]=get_insar_varigram(xin,yin,zin,rmin,rmax,dr,Nmax);
%
%     input:   xin ---- x coordinates of the noise [m]
%              yin ---- y coordinates of the noise [m]
%              zin ---- noise [cm]
%              rmin ---- minimum distance
%              rmax ---- maximum distance
%              dr   ---- distance step
%              Nmax ---- number pairs to use
%
%    C(r) = sigma * exp(-r/L);
%
%    Matlab version by Kang Wang in  Feb. 2018
%    Last updated by Kang Wang in Sep. 2018
%
%%%%% relationship between semivarigram and covariance
% http://pro.arcgis.com/en/pro-app/help/analysis/geostatistical-analyst/semivariogram-and-covariance-functions.htm
%
%  gamma(si,sj)=sill-Cov(si,sj);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def compute_insar_varigram(insarobj, rmin, rmax, dr, Nmax, rng=None):
    """
    Estimate the semi-varigram structure within an InSAR interferogram, assuming the noise structure is isotropic,
    and has an exponential function with distance.

    C(r) = sigma * exp(-r/L)
    From Lohman and Simons, 2005.

    :param insarobj: InSAR_1d object
    :param rmin: inner radius, in the same units as insarobj.lon and insarobj.lat
    :param rmax: outer radius, in the same units as insarobj.lon and insarobj.lat
    :param dr: step-size in radius, in the same units as insarobj.lon and insarobj.lat
    :param Nmax: number of samples, a large integer
    :param rng: specific type of random number generate you wish to use
    :return: radius, covariance at each radius, float sigma, and float L
    """
    print("Computing data covariance")
    insarobj = insarobj.remove_nans()
    zin = insarobj.LOS
    # unbiased variance like MATLAB var(...), i.e., divide by n-1
    var_all = np.var(zin, ddof=1)
    print("Total unbiased variance: %f" % var_all)
    npt = len(zin)
    if npt < 3:
        raise ValueError("Not enough valid points after removing NaNs.")
    # random unique index pairs (rows unique, order not canonicalized)
    rng = np.random.default_rng() if rng is None else rng
    idx = rng.integers(0, npt, size=(Nmax, 2))
    idx = np.unique(idx, axis=0)
    i1, i2 = idx[:, 0], idx[:, 1]

    # differences and pairwise distances
    dz = zin[i1] - zin[i2]
    dx = insarobj.lon[i2] - insarobj.lon[i1]
    dy = insarobj.lat[i2] - insarobj.lat[i1]
    r = np.sqrt(dx * dx + dy * dy)

    # distance bins
    rin = np.arange(rmin, rmax + dr, dr)
    if rin.size < 2:
        raise ValueError("Need at least two bin edges.")

    rout_list, zout_list = [], []

    for k in range(rin.size - 1):
        r1, r2 = rin[k], rin[k + 1]
        in_bin = (r >= r1) & (r < r2)
        N = np.count_nonzero(in_bin)
        if N > 2:  # bins with at least 3 pairs
            # structure function Sr ≈ var(dz) with unbiased variance
            Sr = np.var(dz[in_bin], ddof=1)  # % structure function in eq. (2) of Lohman and Simons (2005)
            # covariance from semivariogram relation: C = var_all - Sr/2
            rout_list.append(0.5 * (r1 + r2))
            zout_list.append(var_all - 0.5 * Sr)   # % eq. (3) of Lohman and Simons (2005)

    if not rout_list:
        raise RuntimeError("No bins had at least 3 pairs. Try lowering rmin/rmax, increasing Nmax, or enlarging dr.")

    rout = np.asarray(rout_list)
    zout = np.asarray(zout_list)

    # sort by r
    order = np.argsort(rout)
    rnew = rout[order]
    znew = zout[order]

    # fit C(r) = sigma * exp(-r/L)  #
    def model(r_, sigma0, L0):
        return sigma0 * np.exp(-r_ / L0)

    # initial guess similar to MATLAB x0=[1, 1e4]
    p0 = (max(znew.max(), 1e-12), max((rmax - rmin) / 3.0, 1.0))
    # bounds: sigma >= 0, L > 0
    popt, _ = curve_fit(model, rnew, znew, p0=p0, bounds=([0.0, 1e-12], [np.inf, np.inf]))
    sigma, L = map(float, popt)

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


def build_Cd(insarobj, sigma, L):
    """
    Construct a formal covariance matrix given the empirically-derived exponential approximation of the cov matrix.

    :param insarobj: an object of type InSAR1Dobj
    :param sigma: the variance at zero distance
    :param L: the decay constant associated with distance between points
    :return: matrix of size
    """
    sizeof = len(insarobj.LOS)
    covd = np.zeros((sizeof, sizeof))
    for i in range(sizeof):
        for j in range(i, sizeof):
            x1, y1 = insarobj.lon[i], insarobj.lat[i]
            x2, y2 = insarobj.lon[j], insarobj.lat[j]
            distance = np.sqrt((x1-x2)**2 + (y1-y2)**2)
            cov = sigma * np.exp(-distance / L)  # this is sigma^2, but it comes from varigram function as sigma
            covd[i][j] = cov
            covd[j][i] = cov
    return covd


def plot_full_covd(covd, outfilename):
    """ Plot the covariance matrix to see its structure. """
    print("Plotting file %s " % outfilename)
    plt.figure(figsize=(7, 7), dpi=300)
    plt.imshow(covd)
    plt.colorbar(label='Covariance (mm^2)')
    plt.savefig(outfilename)
    plt.close()
    return
