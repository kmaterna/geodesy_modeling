"""
Create a quadtree downsampling of an InSAR_2D object.
Start with the notion that your data fits into a square of side-length 2^n, where you solve for n.
Then split up the leaves of the quad-tree over and over again.
NOTE: This will basically turn your InSAR1d object into InSAR1d object
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from ..InSAR_1D_Object.class_model import Insar1dObject
from tectonic_utils.geodesy import haversine


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


def nearest_distance_to_fault(target_lon, target_lat, fault_pts_lon, fault_pts_lat, distance=100):
    """
    Find the nearest point on the fault trace to a target point of interest.
    This is useful in case you want a customized downsampling scheme involving distance to a major fault.

    :param target_lon: float, longitude
    :param target_lat: float, latitude
    :param fault_pts_lon: list of longitudes
    :param fault_pts_lat: list of latitudes
    :param distance: a starting value larger than the maximum distance we are likely to ever encounter, in km.
    """
    for xf, yf in zip(fault_pts_lon, fault_pts_lat):
        d = haversine.distance((yf, xf), (target_lat, target_lon))
        if d < distance:
            distance = d
    return distance


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
    LOS[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.LOS
    coh[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.coherence
    lkvE[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_E
    lkvN[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_N
    lkvU[ny_padding:ny_padding + ny_data, nx_padding:nx_padding + nx_data] = insar2d_obj.lkv_U
    lats = np.multiply(np.ones((np.power(2, max_square),)), np.nan)
    lons = np.multiply(np.ones((np.power(2, max_square),)), np.nan)
    lats[ny_padding:ny_padding + ny_data] = insar2d_obj.lat
    lons[nx_padding:nx_padding + nx_data] = insar2d_obj.lon
    return lons, lats, LOS, coh, lkvE, lkvN, lkvU


def get_center_of_box(Aline):
    """ Get the indices of the center point of a quadtree leaf. """
    y_coord = int(np.mean((Aline[0], Aline[1])))
    x_coord = int(np.mean((Aline[2], Aline[3])))
    return x_coord, y_coord


def get_mean_value_in_box(Aline, grid_data):
    chip = grid_data[int(Aline[0]):int(Aline[1]), int(Aline[2]):int(Aline[3])]
    mean_value = np.nanmean(chip)
    return mean_value


def draw_box_coordinates(Aline):
    """ Get the vertices for plotting a box that surrounds a quadtree leaf. """
    x = [Aline[2], Aline[3], Aline[3], Aline[2], Aline[2]]
    y = [Aline[0], Aline[0], Aline[1], Aline[1], Aline[0]]
    return x, y


def draw_box_coordinates_geographic(Aline, lons, lats):
    """ Get the vertices for plotting a box that surrounds a quadtree leaf, in Lon/Lat coordinates. """
    x = [lons[int(Aline[2])], lons[int(Aline[3])], lons[int(Aline[3])], lons[int(Aline[2])], lons[int(Aline[2])]]
    y = [lats[int(Aline[0])], lats[int(Aline[0])], lats[int(Aline[1])], lats[int(Aline[1])], lats[int(Aline[0])]]
    return x, y


def summarize_quadtree_pixels(Atable, lons, lats, LOS, coh, lkvE, lkvN, lkvU):
    """
    Get the relevant information for each leaf that becomes one pixel. Turns them into InSAR1D pixels.
    Also removes any leaves that have mostly nans, i.e., outside of the main data box.

    :param Atable: list of quadtree leaves in format [IMIN IMAX JMIN JMAX VAR NGOOD NTOTAL]
    :param lons: 1d list
    :param lats: 1d list
    :param LOS: 2d array of LOS values, nan-padded
    :param coh: 2d array of coherence values, nan-padded
    :param lkvE: 2d array of LKV east component, nan-padded
    :param lkvN: 2d array of LKV north component, nan-padded
    :param lkvU: 2d array of LKV up component, nan-padded
    :return: list of pixels in InSAR1D format, array of valid quadtree leaves
    """
    valid_A_list = None
    lon_list, lat_list, LOS_list, coh_list, lkvE_list, lkvN_list, lkvU_list = [], [], [], [], [], [], []
    for leaf in Atable:
        xcoord, ycoord = get_center_of_box(leaf)
        if ~np.isnan(lons[xcoord]) and ~np.isnan(lats[ycoord]):
            lon_list.append(lons[xcoord])
            lat_list.append(lats[ycoord])
            LOS_list.append(get_mean_value_in_box(leaf, LOS))
            coh_list.append(get_mean_value_in_box(leaf, coh))
            lkvE_list.append(get_mean_value_in_box(leaf, lkvE))
            lkvN_list.append(get_mean_value_in_box(leaf, lkvN))
            lkvU_list.append(get_mean_value_in_box(leaf, lkvU))
            if valid_A_list is None:
                valid_A_list = leaf
            else:
                valid_A_list = np.vstack([valid_A_list, leaf])
    insar1d_pixels = Insar1dObject(lon=np.array(lon_list), lat=np.array(lat_list),
                                   LOS=np.array(LOS_list), LOS_unc=np.zeros(np.shape(LOS_list)),
                                   lkv_E=np.array(lkvE_list), lkv_N=np.array(lkvN_list), lkv_U=np.array(lkvU_list),
                                   coherence=np.array(coh_list))
    return insar1d_pixels, valid_A_list


def do_quadtree_downsampling(insar2d_obj, var_max=0.5, nx_padding=0, ny_padding=0, min_pixels=4,
                             custom_function=None):
    """
    Data structure: [IMIN IMAX JMIN JMAX VAR NGOOD NTOTAL]

    :param insar2d_obj:
    :param var_max: maximum variance allowed in a single leaf
    :param nx_padding: option to move the dataset horizontally within the square of power of two
    :param ny_padding: option to move the dataset vertically within the square of power of two
    :param min_pixels: the size of the smallest leaf, default 4.
    :param custom_function: optional conditional function that can be used to split or not-split the leaves
    :return:
    """
    print("Entering quadtree downsampling routine")

    # Create square box with side length equal to a power of two
    lons, lats, LOS, coh, lkvE, lkvN, lkvU = create_padded_boxes(insar2d_obj, nx_padding, ny_padding)
    iter_max = 100

    # Goal: Split each leaf if the variance is above a certain value and there are >50% good pixels
    def split_condition(Amatrix):
        return (Amatrix[:, 4] > var_max) & ((Amatrix[:, 5]/Amatrix[:, 6]) > 0.5) & (Amatrix[:, 6] > min_pixels)

    # If the user wants, you can redefine the conditional function used to split the leaves.
    # Otherwise, we will just use the default split condition based on variance, nans, and minimum size.
    if custom_function is None:
        split_condition_function = split_condition
    else:
        split_condition_function = custom_function

    # Initialize the routine by a single leaf over the entire box
    IMIN, JMIN = 0, 0
    IMAX, JMAX = np.shape(LOS)
    IMAX = IMAX - 1  # converting matlab to python, zero-indexed
    JMAX = JMAX - 1  # converting matlab to python, zero-indexed
    VAR = np.nanvar(LOS)
    NGOOD = np.sum(~np.isnan(LOS))
    NTOTAL = np.size(LOS)
    A = np.array([[IMIN, IMAX, JMIN, JMAX, VAR, NGOOD, NTOTAL]])  # the first leaf

    IOsplit = np.array([True])  # force the code to split the first leaf

    iterations = 0
    while np.any(IOsplit) and iterations < iter_max:
        iterations = iterations + 1
        print("Iteration ", iterations)
        print("Number of leaves: ", len(A))
        Asave = A[~IOsplit, :]  # Extract any leaves which we are NOT splitting.
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

        IOsplit = split_condition_function(A)

    print('Ending after iter ', iterations, 'with ', len(A), 'leaves')

    # Summarize results
    insar1d_object, valid_A = summarize_quadtree_pixels(A, lons, lats, LOS, coh, lkvE, lkvN, lkvU)
    print("Resulting in %d valid pixels" % len(valid_A))

    # Visualize the results really quickly
    plt.figure(dpi=300, figsize=(10, 10))
    plt.imshow(LOS)
    for i in range(len(A)):
        x, y = draw_box_coordinates(A[i, :])
        plt.plot(x, y, color='black')
    plt.xlabel("Array Index")
    plt.ylabel("Array Index")
    plt.title("Preliminary Quadtree Downsampling Scheme")
    print("Plotting initial results in quadtree_downsampling_preliminary.png")
    plt.savefig("quadtree_downsampling_preliminary.png")

    return insar1d_object, valid_A, lons, lats


def quadtree_outputs(insar2d_obj, insar1d_object, valid_A, lons, lats, plotting_file="quadtree_ds_results.png",
                     outfile="qt_ds_pixels.txt"):
    """
    Plot the outputs of a quadtree calcultion in a 2-panel figure showing all the leaves.
    Also write the output leaves into a table showing each pixel.

    :param insar2d_obj: the original insar_2d object
    :param insar1d_object: the resulting insar_1d object
    :param valid_A: the quadtree leaves matrix that tells how big each leaf is
    :param lons: the lon array that corresponds with the A matrix
    :param lats: the lat array that corresponds with the A matrix
    :param plotting_file: string, filename for png
    :param outfile: string, filename for txt pixels
    :return:
    """
    print("Plotting the results of quadtree downsampling in %s" % plotting_file)
    fig, axarr = plt.subplots(1, 2, dpi=300, figsize=(14, 10), constrained_layout=True)
    vmin = np.nanmin(insar2d_obj.LOS)
    vmax = np.nanmax(insar2d_obj.LOS)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis

    # First panel: Just the original data + boxes for each leaf
    im = axarr[0].imshow(np.flipud(insar2d_obj.LOS), cmap=cmap, norm=norm, extent=(np.nanmin(lons), np.nanmax(lons),
                                                                                   np.nanmin(lats), np.nanmax(lats)))
    for i in range(len(valid_A)):
        x, y = draw_box_coordinates_geographic(valid_A[i, :], lons, lats)
        axarr[0].plot(x, y, color='black')

    # Second panel: Colored values of LOS for each leaf.
    patches, rect_vals = [], []
    for i in range(len(valid_A)):
        w = lons[int(valid_A[i, 3])] - lons[int(valid_A[i, 2])]
        h = lats[int(valid_A[i, 1])] - lats[int(valid_A[i, 0])]
        if ~np.isnan(w) and ~np.isnan(h):
            patches.append(Rectangle((lons[int(valid_A[i, 2])], lats[int(valid_A[i, 0])]), w, h))
            rect_vals.append(insar1d_object.LOS[i])

    for i in range(len(valid_A)):
        x, y = draw_box_coordinates_geographic(valid_A[i, :], lons, lats)
        axarr[1].plot(x, y, color='black', linewidth=0.01)
    coll = PatchCollection(patches, cmap=cmap, norm=norm, edgecolor='k', linewidth=0.01)
    coll.set_array(np.array(rect_vals))
    axarr[1].add_collection(coll)
    axarr[1].set_aspect('equal')
    axarr[1].set_title('%d Quadtree rectangles' % (len(valid_A)))
    axarr[1].set_xlabel('Longitude')
    axarr[1].set_ylabel('Latitude')

    _cbar = fig.colorbar(im, ax=axarr.ravel().tolist(), label='Value')
    plt.savefig(plotting_file)

    write_insar_ds_invertible_format(insar1d_object, valid_A, lons, lats, outfile)
    return


def write_insar_ds_invertible_format(InSAR_obj, valid_A, lons, lats, filename):
    """
    Write InSAR displacements from quadtree downsampling into file that can be inverted.
    Write one header line and multiple data lines.
    InSAR_obj is in mm

    :param InSAR_obj: insar1d object
    :param valid_A: matrix of leaves in quadtree downsampling
    :param lons: the lon array that corresponds with the A matrix
    :param lats: the lat array that corresponds with the A matrix
    :param filename: string, filename for output text file
    """
    print("Writing InSAR displacements into file %s " % filename)
    ofile = open(filename, 'w')
    ofile.write("# Displacements: Lon, Lat, disp(mm), sigma, coherence, unitE, unitN, unitU, TLx, BRx, TLy, BRy\n")
    for i in range(len(InSAR_obj.lon)):
        x, y = draw_box_coordinates_geographic(valid_A[i, :], lons, lats)
        tlx, brx = x[0], x[2]
        tly, bry = y[0], y[2]
        if np.isnan(InSAR_obj.LOS[i]):
            continue
        if np.isnan(tlx) or np.isnan(brx) or np.isnan(tly) or np.isnan(bry):
            continue
        else:
            ofile.write('%f %f ' % (InSAR_obj.lon[i], InSAR_obj.lat[i]))
            ofile.write('%f %f %f ' % (InSAR_obj.LOS[i], InSAR_obj.LOS_unc[i], InSAR_obj.coherence[i]))  # mm
            ofile.write('%f %f %f ' % (InSAR_obj.lkv_E[i], InSAR_obj.lkv_N[i], InSAR_obj.lkv_U[i]))
            ofile.write('%f %f %f %f\n' % (tlx, brx, tly, bry))
    ofile.close()
    return


def qtree_compute(insar2d_obj, var_max=0.5, nx_padding=0, ny_padding=0, min_pixels=4,
                  plotting_file="quadtree_ds_results.png", outfile="qt_ds_pixels.txt",
                  custom_function=None):
    insar1d_object, valid_A, lons, lats = do_quadtree_downsampling(insar2d_obj, var_max, nx_padding,
                                                                   ny_padding, min_pixels, custom_function)
    quadtree_outputs(insar2d_obj, insar1d_object, valid_A, lons, lats, plotting_file, outfile)
    return
