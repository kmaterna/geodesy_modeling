"""
Write and output functions for 1D InSAR data format
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from tectonic_utils.geodesy import insar_vector_functions
from geodesy_modeling import general_utils


def write_insar_invertible_format(InSAR_obj, filename, unc_min=0):
    """
    Write InSAR displacements into insar file that can be inverted.
    Write one header line and multiple data lines.
    InSAR_obj is in mm, and written out is in meters
    """
    print("Writing InSAR displacements into file %s " % filename)
    ofile = open(filename, 'w')
    ofile.write("# InSAR Displacements: Lon, Lat, disp(m), sigma, unitE, unitN, unitU \n")
    for i in range(len(InSAR_obj.lon)):
        if np.isnan(InSAR_obj.LOS[i]):
            continue
        else:
            if InSAR_obj.LOS_unc is not None:
                std = InSAR_obj.LOS_unc[i] * 0.001  # in m
                if std < unc_min:
                    std = unc_min
            else:   # sometimes there's an error code in LOS_unc field
                std = unc_min
            ofile.write('%f %f ' % (InSAR_obj.lon[i], InSAR_obj.lat[i]))
            ofile.write('%f %f ' % (0.001 * InSAR_obj.LOS[i], std))  # writing in m
            ofile.write('%f %f %f\n' % (InSAR_obj.lkv_E[i], InSAR_obj.lkv_N[i], InSAR_obj.lkv_U[i]))
    ofile.close()
    return


def plot_insar(InSAR_Obj, plotname, vmin=None, vmax=None, lons_annot=None, lats_annot=None, refpix=None,
               title=None, colormap='viridis'):
    """lons_annot and lat_annot are for lines to annotate the plot, such as faults or field boundaries"""
    print("Plotting insar in file %s " % plotname)
    fig = plt.figure(dpi=300, figsize=(5, 5))
    if vmin is not None:
        im = plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=28, cmap=colormap, vmin=vmin, vmax=vmax)
    else:
        im = plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=28, cmap=colormap)
    if lons_annot:
        plt.plot(lons_annot, lats_annot, color='darkred')
    if refpix:  # expect a list or tuple, (lon, lat)
        plt.plot(refpix[0], refpix[1], '.', color='red')
    if title:
        plt.title(title)
    else:
        if InSAR_Obj.starttime is None:
            starttime, endtime = ' ', ' '
        else:
            starttime = dt.datetime.strftime(InSAR_Obj.starttime, "%Y-%m")
            endtime = dt.datetime.strftime(InSAR_Obj.endtime, "%Y-%m")
        plt.title("InSAR with %d Points from %s to %s" % (len(InSAR_Obj.LOS), starttime, endtime))
    cb = fig.colorbar(im, ax=plt.gca())
    cb.set_label('Displacement (mm)', fontsize=16)
    plt.savefig(plotname)
    return


def plot_look_vectors(Data, plotname):
    """
    Visualize the look vectors for gut-checking.
    """
    print("Plotting look vectors in %s " % plotname)
    f, axarr = plt.subplots(1, 3, figsize=(11, 4), dpi=300)
    title_fontsize = 14
    flight_heading, _ = insar_vector_functions.look_vector2flight_incidence_angles(Data.lkv_E[0], Data.lkv_N[0],
                                                                                   Data.lkv_U[0])
    x_flight, y_flight, x_los, y_los = general_utils.get_los_and_flight_vectors(flight_heading, 'right')

    cm = plt.cm.get_cmap('viridis')

    im1 = axarr[0].scatter(Data.lon, Data.lat, c=Data.lkv_E, cmap=cm)
    cbaxes = inset_axes(axarr[0], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[0].bbox)
    plt.colorbar(im1, cax=cbaxes, orientation='vertical')
    axarr[0].set_title("Look Vector East (ground to sat)", fontsize=title_fontsize)
    axarr[0].quiver(0.1, 0.85, 2*x_flight, 2*y_flight, transform=axarr[0].transAxes, scale=8, scale_units='inches')
    axarr[0].quiver(0.1, 0.85, x_los, y_los, transform=axarr[0].transAxes, scale=8, scale_units='inches')

    im2 = axarr[1].scatter(Data.lon, Data.lat, c=Data.lkv_N, cmap=cm)
    axarr[1].set_yticks([])
    cbaxes = inset_axes(axarr[1], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[1].bbox)
    plt.colorbar(im2, cax=cbaxes, orientation='vertical')
    axarr[1].set_title("Look Vector North", fontsize=title_fontsize)

    im3 = axarr[2].scatter(Data.lon, Data.lat, c=Data.lkv_U, cmap=cm)
    axarr[2].set_yticks([])
    cbaxes = inset_axes(axarr[2], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[2].bbox)
    plt.colorbar(im3, cax=cbaxes, orientation='vertical')
    axarr[2].set_title("Look Vector Up", fontsize=title_fontsize)

    plt.savefig(plotname)
    return
