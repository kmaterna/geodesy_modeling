
import numpy as np


def read_stress_tensors_from_espy(filename):
    """
    Read a file with a list of stress tensors from ll_strains in elastic_stresses_py
    # Format: lon lat depth_km sigma_xx sigma_xy sigma_xz sigma_yy sigma_yz sigma_zz (kPa)

    :param filename: string
    :return: list of 3x3 matrices, in pascals
    """
    print("Reading filename %s " % filename)
    tau_list = []
    [lon, _, _, s11, s12, s13, s22, s23, s33] = np.loadtxt(filename, unpack=True)
    s11 = np.multiply(s11, 1000)  # convert from kPa to pascals
    s12 = np.multiply(s12, 1000)
    s13 = np.multiply(s13, 1000)
    s22 = np.multiply(s22, 1000)
    s23 = np.multiply(s23, 1000)
    s33 = np.multiply(s33, 1000)
    for i in range(len(lon)):
        tau = np.array([[s11[i], s12[i], s13[i]], [s12[i], s22[i], s23[i]], [s13[i], s23[i], s33[i]]])
        tau_list.append(tau)
    return tau_list
