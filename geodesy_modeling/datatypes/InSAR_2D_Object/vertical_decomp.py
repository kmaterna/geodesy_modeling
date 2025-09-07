"""
From two InSAR LOS objects with co-registered pixels, perform a decomposition.
"""

import numpy as np
import sys


def east_vertical_decomp(object1, object2):
    """
    Take scenes from two different look vectors and convert them into vertical and east,
    assuming north motion is none.
    This is a very common decomposition in InSAR processing.

    :param object1: an object of type InSAR_2D_object
    :param object2: an object of type InSAR_2D_object, with equal array size as Object1
    :return: east array, vertical array
    """
    # Step 1: sanity check.
    if np.shape(object1.LOS) != np.shape(object2.LOS):
        raise ValueError("Error! Object1 and Object2 do not have compatible shapes: " + str(np.shape(object1.LOS)) +
                         " " + str(np.shape(object2.LOS)))
    y, x = np.shape(object1.LOS)
    east, vert = np.zeros(np.shape(object1.LOS)), np.zeros(np.shape(object1.LOS))
    for i in range(y):
        for j in range(x):
            lkv1e, lkv1u = object1.lkv_E[i][j], object1.lkv_U[i][j]
            lkv2e, lkv2u = object2.lkv_E[i][j], object2.lkv_U[i][j]
            disps = np.array([object1.LOS[i][j], object2.LOS[i][j]])
            matrix = np.array([[lkv1e, lkv1u], [lkv2e, lkv2u]])
            solution = np.linalg.solve(matrix, disps)
            east[i][j], vert[i][j] = solution[0], solution[1]
    return east, vert
