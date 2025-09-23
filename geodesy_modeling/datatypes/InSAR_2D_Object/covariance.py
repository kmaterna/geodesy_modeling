"""Create covariance matrix routine for 2d insar object."""


from ..InSAR_1D_Object import covariance


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
