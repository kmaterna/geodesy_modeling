"""
Write/output/plotting functions for InSAR 2D data format
"""

import pygmt
import numpy as np
import matplotlib.pyplot as plt
from tectonic_utils.read_write import netcdf_read_write
from . import inputs
from geodesy_modeling import general_utils
from elastic_stresses_py.PyCoulomb import utilities as overall_utils


def write_InSAR2D(InSAR_Obj, filename):
    """Write grdfile"""
    netcdf_read_write.write_netcdf4(InSAR_Obj.lon, InSAR_Obj.lat, InSAR_Obj.LOS, filename)
    return


def write_insar2D_invertible_format(InSAR_obj, filename):
    """
    Write InSAR 2D displacements into insar text file that can be inverted.
    Write one header line and multiple data lines, with different look vectors for each pixel.
    InSAR_2D_obj is in mm, and written out is in meters
    """
    print("Writing InSAR displacements into file %s" % filename)
    with open(filename, 'w') as ofile:
        ofile.write("# Lon, Lat, Disp, LkvE, LkvN, LkvU, Coh\n")
        ofile.write("# Look direction = %s \n" % InSAR_obj.look_direction)
        X, Y = np.meshgrid(InSAR_obj.lon, InSAR_obj.lat)
        X = X.reshape(-1)
        Y = Y.reshape(-1)
        LOS = InSAR_obj.LOS.reshape(-1)
        lkvE = InSAR_obj.lkv_E.reshape(-1)
        lkvN = InSAR_obj.lkv_N.reshape(-1)
        lkvU = InSAR_obj.lkv_U.reshape(-1)
        coherence = InSAR_obj.coherence.reshape(-1)
        for i in range(len(X)):
            ofile.write("%f %f %f %f %f %f %f\n" % (X[i], Y[i], LOS[i], lkvE[i], lkvN[i], lkvU[i], coherence[i]))
    return


def plot_incidence_azimuth_angle(InSAR_2D_obj, plotname):
    """
    Plot the incidence angle (degrees from the vertical) and the flight azimuth (CW from north)
    for the 2D InSAR object.

    :param InSAR_2D_obj: an InSAR_2D_object
    :param plotname: string
    """
    print("Plotting %s " % plotname)
    az, inc = InSAR_2D_obj.get_azimuth_incidence_grid(InSAR_2D_obj.look_direction)
    inc[inc > 90] = np.nan  # ignore incidence angles greater than 90
    inc[inc < 0] = np.nan

    fig, axs = plt.subplots(1, 2, dpi=300, figsize=(8, 5))
    fig.subplots_adjust(wspace=0.5)  # Increase horizontal space
    d1 = axs[0].imshow(inc, extent=(InSAR_2D_obj.lon.min(), InSAR_2D_obj.lon.max(),
                                    InSAR_2D_obj.lat.min(), InSAR_2D_obj.lat.max()))
    axs[0].set_xlabel('Longitude', fontsize=14)
    axs[0].set_ylabel('Latitude', fontsize=14)
    axs[0].set_title('Incidence')
    _cb = fig.colorbar(d1, label="Incidence (degrees from vertical)", ax=axs[0])
    d2 = axs[1].imshow(az, extent=(InSAR_2D_obj.lon.min(), InSAR_2D_obj.lon.max(),
                                   InSAR_2D_obj.lat.min(), InSAR_2D_obj.lat.max()))
    axs[1].set_xlabel('Longitude', fontsize=14)
    axs[1].set_title('Azimuth, '+InSAR_2D_obj.look_direction+'-looking')
    _cb = fig.colorbar(d2, label="Azimuth (degrees CW from north)", ax=axs[1])

    plt.savefig(plotname)
    plt.close()
    return


def plot_disp_point_annotations(fig, disp_points, disp_points_color=None, cmap=None, region=None):
    """
    :param fig: pygmt/python figure, for plotting
    :param disp_points: list of disp_point objects
    :param disp_points_color: optional 1d array of floats, used for coloring the fill of the symbols using cmap
    :param cmap: optional string, a file for the color map plotting the symbol fill
    :param region: optional array, [W, E, S, N], for plotting the scale vector
    """

    # Unpack modeled displacements
    model_lon = np.array([x.lon for x in disp_points])
    model_lat = np.array([x.lat for x in disp_points])
    model_dE = np.array([x.dE_obs for x in disp_points])
    model_dN = np.array([x.dN_obs for x in disp_points])

    if disp_points_color:  # plot vectors and color-fill
        [scale_arrow, vectext] = overall_utils.define_vector_scale_size(model_dE, model_dN)
        scale = 1.5 / scale_arrow  # empirical scaling, convenient display
        fig.plot(x=model_lon, y=model_lat, style='t0.10i', fill=disp_points_color, pen="thin,black", cmap=cmap)
        fig.plot(x=model_lon, y=model_lat, style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                 direction=[model_dE, model_dN], pen="thin,black")
        fig.plot(x=[region[0]+0.05], y=[region[2]+0.55],  style='v0.2c+e+gblack+h0+p1p,black+z'+str(scale),
                 direction=[[scale_arrow], [0]],  pen="thin,black")  # scale vector
        fig.text(x=[region[0]+0.20], y=[region[2]+0.65], text=vectext+" obs")  # scale label
    else:  # just plot triangles for location
        fig.plot(x=model_lon, y=model_lat, style='t0.07i', fill='black', pen="thin,black")
    return fig


def plot_reference_annotations(fig, refpoint):
    # Plot the reference point
    fig.plot(x=refpoint[0], y=refpoint[1], style='s0.13i', pen="thin,red")
    return fig


def add_text_annotation(fig, text_annot, region, proj):
    fig.text(position='TL', text=text_annot, region=region, projection=proj, font="14p,Helvetica,black",
             offset="0.2/-0.2", pen="0.5p,black")
    return fig


def add_flight_vector(fig, flight_heading, look_dir, region):
    x_flight, y_flight, x_los, y_los = general_utils.get_los_and_flight_vectors(flight_heading, look_dir)
    fig.text(x=region[0], y=region[2], text='LOS', offset=str(0.6 + 0.4 * x_los) + "i/0.45i",
             fill='white', font="10p,Helvetica,black")  # LOS text annotation
    fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
             direction=[[2*x_flight], [2*y_flight]], pen="thin,black", offset="0.6i/0.45i")  # flight vector
    fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
             direction=[[x_los], [y_los]], pen="thin,black", offset="0.6i/0.45i")  # los vector
    return fig


def map_wrapped_insar(grd_filename, plotname, text_annot=None, flight_heading=None, look_dir="right",
                      disp_points=None, region=None, refloc=None, wrapped=False,
                      disp_points_color=None, vmin=-10, vmax=10, label_inc=2):
    """
    Pygmt map of wrapped InSAR deformation with optional annotations

    :param grd_filename: string, filename to be mapped
    :param plotname: string, output file to be created
    :param text_annot: optional string, additional annotation
    :param flight_heading: optional float, degrees CW from north
    :param look_dir: optional string, either 'right' or 'left'
    :param disp_points: optional list of disp_point objects
    :param disp_points_color: optional 1d array of floats, used for coloring the fill of the symbols using cmap
    :param refloc: optional lon/lat for plotting reference pixel
    :param wrapped: default False. If True, the plot is assumed to be wrapped (-pi to pi) on a cycle color scale
    :param region: optional array, [W, E, S, N] for the plot
    :param vmin: optional, minimum for unwrapped colorbar
    :param vmax: optional, maximum for unwrapped colorbar
    :param label_inc: optional, labeling increment for unwrapped colorbar
    """
    InSAR_Obj = inputs.inputs_grd(grd_filename)
    print("Plotting file %s " % plotname)
    proj = 'M4i'
    if not region:
        region = [np.min(InSAR_Obj.lon), np.max(InSAR_Obj.lon), np.min(InSAR_Obj.lat), np.max(InSAR_Obj.lat)]

    if wrapped:  # Options for wrapped phase plots (wrapped is always on the same scale no need to change this)
        pygmt.makecpt(cmap="cyclic", series="-3.14/3.14/0.01", background="o", output="mycpt.cpt")
        title = "Wrapped LOS Displacement"
        label_inc = 1.57  # for wrapped phase
        label = "Phase"
    else:  # currently hard-coded for different applications
        pygmt.makecpt(cmap="polar", series=str(vmin)+"/"+str(vmax)+"/0.01", background="o", output="mycpt.cpt")
        title = "Unwrapped LOS Displacement"
        label = "LOS Deformation (mm)"

    # Build a PyGMT plot
    fig = pygmt.Figure()
    fig.basemap(region=region, projection=proj, frame="+t\""+title+"\"")
    fig.grdimage(grd_filename, region=region, cmap="mycpt.cpt")
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.43/0.06+c" + str(region[2]) + "+w20", frame="1.0")
    if text_annot:
        fig = add_text_annotation(fig, text_annot, region, proj)
    if flight_heading:  # draw vectors for flight direction and look direction
        fig = add_flight_vector(fig, flight_heading, look_dir, region)
    if disp_points:
        fig = plot_disp_point_annotations(fig, disp_points, disp_points_color, "mycpt.cpt", region=region)
    if refloc:  # reference location
        fig = plot_reference_annotations(fig, refloc)

    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", frame=["x"+str(label_inc), "y+L\""+label+"\""])
    fig.savefig(plotname)
    return
