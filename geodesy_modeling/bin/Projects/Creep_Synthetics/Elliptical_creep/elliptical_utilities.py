"""
A set of utilities used in creating elliptical slip distributions out of Okada rectangular fault patches.
This library allows for increasing the size of the rectangular patches as the fault patches go deeper.
"""

import numpy as np
import matplotlib.pyplot as plt


def create_elliptical_slip_array(depths_array, So, a):
    """
    Create an elliptical slip distribution based on surface slip (So) and bottom depth (a).
    The parameter depths_array is usually linearly spaced. It can be larger than the intended crack depth.

    :param depths_array: list, in km.
    :param So: slip at surface, in whatever units
    :param a: bottom depth of the elliptical crack, in km
    :returns: list of depths that are used in the crack, list of slip values that are used.
    """
    input_slip = np.sqrt(So*So * (1 - (np.square(depths_array)/(a**2))))  # Depths are km.
    used_depths = depths_array[~np.isnan(input_slip)]
    input_slip = input_slip[~np.isnan(input_slip)]
    return used_depths, input_slip


def get_width_from_depth(starting_width, depth, thickening_slope=0.27):
    """
    Implements a linear function to determine the desired width of a fault patch as a function of depth.
    We increase the width of rectangle fault patches as we go deeper.
    This is used for approximating elliptical slip distributions with bigger rectangles in the deeper depths.
    It is important for reducing the computational cost with identical results, and for producing smooth Jacobians.

    :param starting_width: float, in km, user-defined parameter to determine the start of the thickening rectangles
    :param depth: float, depth of the target rectangular fault patch
    :param thickening_slope: float, user-defined parameter to determine how quickly width of faults increases with depth
    """
    width = starting_width + thickening_slope * depth
    return width


def get_tuned_depth_arrays(top_depth, max_depth, starting_width, ell_surface_slip, ell_bottom_depth):
    """
    Produce an array of calibrated depth, width, and slip values associated with a column of rectangular fault patches.
    The size of the fault patches increases with depth.
    The slip values conform to an elliptical slip distribution.

    :param top_depth: float, km, top of the column of faults.
    :param max_depth: float, km, bottom of the entire column of faults that could be slipping.
    :param starting_width: float, initial guess for the width of the top rectangular fault element.
    :param ell_surface_slip: in meters, the slip at the top of the crack.
    :param ell_bottom_depth: float, bottom depth of the elliptical crack
    :returns: lists of top_depths, widths, and slip associated with each calibrated rectangle within elliptical crack
    """
    counter = top_depth  # start at the top of the column of fault patches
    top_depths, widths = [], []  # initialize empty arrays for top depth and width of each patch.
    while counter <= max_depth:
        top_depths.append(counter)
        width = get_width_from_depth(starting_width, counter)
        widths.append(width)
        counter += width
    top_depths, widths = np.array(top_depths), np.array(widths)
    middle_depths = np.array([x+y/2 for x, y in zip(top_depths, widths)])  # extract the middle depth
    # for getting elliptical slip values in the next step
    exp_depths, exp_slip = create_elliptical_slip_array(middle_depths, ell_surface_slip, ell_bottom_depth)  # middles
    exp_widths = widths[0:len(exp_depths)]
    exp_top_depths = top_depths[0:len(exp_depths)]
    return exp_top_depths, exp_widths, exp_slip


def plot_rectangles(ax, top_depths, widths, slips, color):
    """ Small function to plot the edges of many rectangular fault patches in profile. """
    for i in range(len(top_depths)):
        x = [0, slips[i], slips[i], 0]
        y = [top_depths[i], top_depths[i], top_depths[i]+widths[i], top_depths[i]+widths[i]]
        ax.plot(x, y, color=color)
    return ax


def show_elliptical_distribution(top_depth, max_depth, sampling_interval, surface_slip=0.02, bottom_depth=2.3,
                                 outfile="Modeling_ellipse.png"):
    """
    Create a plot to show how a given elliptical slip distribution will be produced from a number of
    rectangular fault elements of increasing size with increasing depth.

    :param top_depth: float, in km, depth of the top of the column of fault patches.
    :param max_depth: float, in km, depth of the fault patch at bottom of the column. Not all patches must have slip.
    :param sampling_interval: in km, default width of rectangular fault patches in the constant-width case
    :param surface_slip: float, in meters
    :param bottom_depth: float, in km
    :param outfile: string, name of file where PNG will be written.
    """
    # for initial rectangular distributions:
    sample_depths = np.arange(top_depth, max_depth, sampling_interval)
    model_depths = np.arange(top_depth, max_depth, 0.05)  # very fine sampling, just for 1d line on the plot

    # Create a model of the fault using equally-spaced rectangles
    rect_depths, rect_slip = create_elliptical_slip_array(sample_depths, surface_slip, bottom_depth)  # gets top depths
    rect_widths = sampling_interval*np.ones(np.shape(rect_depths))

    # Create a model of the fault using non-equally-spaced rectangles
    exp_depths, exp_widths, exp_slip = get_tuned_depth_arrays(top_depth, max_depth, sampling_interval, surface_slip,
                                                              bottom_depth)  # gets top depths

    # Create a densely sampled elliptical curve for plotting on top of the Okada rectangles
    model_depths, model_slip = create_elliptical_slip_array(model_depths, surface_slip, bottom_depth)

    print("Plotting file %s " % outfile)
    f, axarr = plt.subplots(1, 2, dpi=300, figsize=(12, 8))
    axarr[0].plot(model_slip, model_depths)
    plot_rectangles(axarr[0], rect_depths, rect_widths, rect_slip, color='red')
    axarr[0].invert_yaxis()
    axarr[0].set_title("Constant-width: N=" + str(len(rect_depths)))
    axarr[0].set_xlabel("Slip (m)")
    axarr[0].set_ylabel("Depth (km)")
    axarr[1].plot(model_slip, model_depths)
    plot_rectangles(axarr[1], exp_depths, exp_widths, exp_slip, color='black')
    axarr[1].invert_yaxis()
    axarr[1].set_title("Variable-width: N=" + str(len(exp_depths)))
    axarr[1].set_xlabel("Slip (m)")
    axarr[1].set_ylabel("Depth (km)")
    plt.savefig(outfile)
    return
