"""
Remove a ramp or a constant from 1D InSAR object
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from .class_model import Insar1dObject


def fit_ramp(InSAR_Obj: Insar1dObject):
    """
    Fit the best-fitting plane parameters to an object of InSAR data (1d).

    :param InSAR_Obj: 1D InSAR object
    :return: a, b, c
    """
    nonan_obj = InSAR_Obj.remove_nans()  # remove nans before fitting
    # Solve for the best-fitting ramp of equation ax + by + c = z
    Z = []
    A = np.zeros((len(nonan_obj.lon), 3))
    for i in range(len(nonan_obj.lon)):
        A[i, :] = [nonan_obj.lon[i], nonan_obj.lat[i], 1]
        Z.append(nonan_obj.LOS[i])
    model = np.linalg.lstsq(A, Z)
    model = model[0]
    a, b, c = model[0], model[1], model[2]
    return a, b, c


def remove_best_fit_ramp(InSAR_Obj: Insar1dObject):
    """
    Find the best-fitting ramp from a set of InSAR observations and remove it from the data.
    Plane equation: ax + by + c = z
    Solving Ax = B

    :param InSAR_Obj: 1D insar object
    :returns: 1D insar object
    """
    a, b, c = fit_ramp(InSAR_Obj)  # Solve for the best-fitting ramp of equation ax + by + c = z
    new_insar = InSAR_Obj.subtract_ramp(a, b, c)  # Removing the planar model
    return new_insar


def plot_ramp_results(Obj1, Obj2, filename):
    """
    Plot difference between the object with and without the ramp removal.

    :param Obj1: Original InSAR 1d object
    :param Obj2: De-ramped InSAR 1d object
    :param filename: string
    :return:
    """
    vmin = np.nanmin(Obj1.LOS)
    vmax = np.nanmax(Obj1.LOS)

    f, axarr = plt.subplots(1, 2, figsize=(12, 7), dpi=300)
    axarr[0].scatter(Obj1.lon, Obj1.lat, c=Obj1.LOS, cmap='rainbow', vmin=vmin, vmax=vmax)
    axarr[0].set_title('Object With Ramps')
    axarr[1].scatter(Obj2.lon, Obj2.lat, c=Obj2.LOS, cmap='rainbow', vmin=vmin, vmax=vmax)
    axarr[1].set_title('Object Without Ramps')

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False)
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='rainbow')
    custom_cmap.set_array(np.arange(vmin, vmax))
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical')
    cb.set_label('Displacement (mm)', fontsize=18)
    cb.ax.tick_params(labelsize=12)

    plt.savefig(filename)
    return
