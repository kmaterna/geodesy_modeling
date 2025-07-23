"""
Create a quadtree downsampling of an InSAR_2D object.
Start with the notion that your data fits into a square of side-length 2^n, where you solve for n.
Then split up the leaves of the quad-tree over and over again.
"""

import numpy as np
import matplotlib.pyplot as plt


def get_square_2n_size(LOS):
    """
    Find the appropriate power of 2 that is needed to cover your dataset. The result will be a square.

    :param LOS: 2D array
    :return: integer, power of two that is needed to fit this data
    """
    ny, nx = np.shape(LOS)
    i = 1
    while np.power(2, i) < np.max((nx, ny)):
        i += 1
    return i


def create_padded_boxes(insar2d_obj, nx_padding, ny_padding):
    """
    Pad the box of data into a square of 2^n, surrounded by nans.

    :param insar2d_obj: data of type InSAR2d
    :param nx_padding: integer, number of places away from the origin to start the x-padding
    :param ny_padding: integer, number of places away from the origin to start the y-padding
    :return:
    """
    ny_data, nx_data = np.shape(insar2d_obj.LOS)
    max_square = get_square_2n_size(insar2d_obj.LOS)
    padded_box = np.multiply(np.ones((np.power(2, max_square), np.power(2, max_square))), np.nan)
    LOS = padded_box.copy()
    coh = padded_box.copy()
    lkvE = padded_box.copy()
    lkvN = padded_box.copy()
    lkvU = padded_box.copy()
    LOS[ny_padding:ny_padding+ny_data, nx_padding:nx_padding+nx_data] = insar2d_obj.LOS
    coh[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.coherence
    lkvE[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_E
    lkvN[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_N
    lkvU[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_U
    return LOS, coh, lkvE, lkvN, lkvU


def do_quadtree(insar2d_obj, var_max=0.5, nx_padding=0, ny_padding=0):
    """
    Data structure: IMIN IMAX JMIN JMAX VAR NGOOD NTOTAL

    :return:
    """
    LOS, coh, lkvE, lkvN, lkvU = create_padded_boxes(insar2d_obj, nx_padding, ny_padding)
    var_max = 7.5
    iter_max = 100

    # Goal: Split each leaf if the variance is above a certain value and there are >50% good pixels
    def split_condition(A):
        return (A[:, 4] > var_max) & ((A[:, 5]/A[:, 6]) > 0.5) & (A[:, 6] > 4)

    # Initialize the routine with a single leaf
    IMIN, JMIN = 0, 0
    IMAX, JMAX = np.shape(LOS)
    VAR = np.nanvar(LOS)
    NGOOD = np.sum(~np.isnan(LOS))
    NTOTAL = np.size(LOS)
    A = np.array([[IMIN, IMAX, JMIN, JMAX, VAR, NGOOD, NTOTAL]])  # the first leaf

    print("\n\n\n")
    IOsplit = np.array([True])  # force the code to split the first leaf

    iterations = 0
    while np.any(IOsplit) and iterations < iter_max:
        iterations = iterations + 1
        print("Iteration ", iterations)
        print("Number of leaves: ", len(A))
        Asave = A[~IOsplit, :]  # Save any leaves which we are NOT splitting.
        A = A[IOsplit, :]  # Extract leaves which we are going to split.

        # During each iteration, take every leaf that needs splitting and split it into 4 leaves.
        Anew = None
        for k in range(len(A)):
            Imin, Imax = int(A[k, 0]), int(A[k, 1])
            Jmin, Jmax = int(A[k, 2]), int(A[k, 3])
            Imean = np.mean((Imin, Imax))
            Jmean = np.mean((Jmin, Jmax))
            Asplit = np.zeros((4, 7))  # Create four new rows
            Asplit[:, 0:4] = np.array([[Imin, int(np.floor(Imean)), Jmin, int(np.floor(Jmean))],
                                       [int(np.ceil(Imean)), Imax, Jmin, int(np.floor(Jmean))],
                                       [int(np.ceil(Imean)), Imax, int(np.ceil(Jmean)), Jmax],
                                       [Imin, int(np.floor(Imean)), int(np.ceil(Jmean)), Jmax]])
            for n in range(0, 4):
                chip = LOS[int(Asplit[n, 0]):int(Asplit[n, 1]), int(Asplit[n, 2]):int(Asplit[n, 3])]
                Asplit[n, 4] = np.nanvar(chip[:])
                Asplit[n, 5] = np.sum(~np.isnan(chip[:]))  # number of good values in the chip
                Asplit[n, 6] = np.size(chip)   # number of total values in the chip
            if Anew is None:
                Anew = Asplit
            else:
                Anew = np.vstack([Anew, Asplit])
        if len(Asave) == 0:
            A = Anew
        else:
            A = np.vstack([Asave, Anew])

        IOsplit = split_condition(A)
        # VAR=A(:,5); NGOOD=A(:,6); NTOTAL=A(:,7);

    print('Ending after iter ', iterations, 'with ', len(A), 'leaves')

    def draw_box(Aline):
        x = [Aline[2], Aline[3], Aline[3], Aline[2], Aline[2]]
        y = [Aline[0], Aline[0], Aline[1], Aline[1], Aline[0]]
        return x, y

    plt.figure(dpi=300, figsize=(10, 10))
    plt.imshow(LOS)
    for i in range(len(A)):
        x, y = draw_box(A[i, :])
        plt.plot(x, y, color='black')
    plt.savefig("Quadtree_downsampling.png")
