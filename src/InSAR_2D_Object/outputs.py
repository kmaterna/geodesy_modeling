"""
Write and output functions for InSAR 2D data format
"""

import pygmt
import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from Tectonic_Utils.geodesy import insar_vector_functions
from . import inputs


def write_InSAR2D(InSAR_Obj, filename):
    """Write grdfile"""
    netcdf_read_write.produce_output_netcdf(InSAR_Obj.lon, InSAR_Obj.lat, InSAR_Obj.LOS, '', filename);
    return;


def write_insar2D_invertible_format(_InSAR_obj, _unc_min, filename):
    """
    Write InSAR 2D displacements into insar text file that can be inverted.
    Write one header line and multiple data lines, with different look vectors for each pixel.
    InSAR_2D_obj is in mm, and written out is in meters
    """
    print("Writing InSAR displacements into file %s - function not yet written!" % filename);
    return;


def plot_disp_point_annotations(fig, disp_points):
    # Unpack modeled displacements
    model_lon = np.array([x.lon for x in disp_points]);
    model_lat = np.array([x.lat for x in disp_points]);
    fig.plot(x=model_lon, y=model_lat, style='t0.07i', color='black', pen="thin,black");
    return fig;


def plot_reference_annotations(fig, refpoint):
    # Plot the reference point
    fig.plot(x=refpoint[0], y=refpoint[1], style='t0.04i', color='red', pen="thin,red");
    return fig;


def add_text_annotation(fig, text_annot, region, proj):
    fig.text(position='TL', text=text_annot, region=region, projection=proj, font="14p,Helvetica,black",
             offset="0.2/-0.2", pen="0.5p,black");
    return fig;


def add_flight_vector(fig, flight_heading, look_dir, region):
    [x_flight, y_flight] = insar_vector_functions.get_unit_vector_from_heading(flight_heading);
    if look_dir == 'right':
        [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading + 90);
    else:
        [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading - 90);
    fig.text(x=region[0], y=region[2], text='LOS', offset=str(0.6 + 0.4 * x_los) + "i/0.45i",
             fill='white', font="10p,Helvetica,black");  # LOS text annotation
    fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
             direction=[[x_flight], [y_flight]], pen="thin,black", offset="0.6i/0.45i");  # flight vector
    fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
             direction=[[x_los / 2], [y_los / 2]], pen="thin,black", offset="0.6i/0.45i");  # los vector
    return fig;


def map_wrapped_insar(grd_filename, plotname, text_annot=None, flight_heading=None, look_dir="right",
                      disp_points=None, region=None, refloc=None, wrapped=False):
    """
    Pygmt map of wrapped InSAR deformation with optional annotations
    """
    InSAR_Obj = inputs.inputs_grd(grd_filename);
    print("Plotting file %s " % grd_filename);
    proj = 'M4i'
    if not region:
        region = [np.min(InSAR_Obj.lon), np.max(InSAR_Obj.lon), np.min(InSAR_Obj.lat), np.max(InSAR_Obj.lat)];

    if wrapped:  # Hard-coded options for wrapped phase plots
        pygmt.makecpt(cmap="cyclic", series="-3.14/3.14/0.01", background="o", output="mycpt.cpt");
        title = "Wrapped LOS Displacement"
        label_inc = 1.57;  # for wrapped phase
        label = "Phase";
    else:  # currently hard-coded for different applications
        pygmt.makecpt(cmap="roma", series="-30/30/0.01", background="o", output="mycpt.cpt");
        title = "Unwrapped LOS Displacement"
        label_inc = 5;  # For unwrapped displacement
        label = "LOS Deformation (mm)";

    # Build a PyGMT plot
    fig = pygmt.Figure();
    fig.basemap(region=region, projection=proj, frame="+t\""+title+"\"");
    fig.grdimage(grd_filename, region=region, cmap="mycpt.cpt");
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.43/0.06+c" + str(region[2]) + "+w20", frame="1.0");
    if text_annot:
        fig = add_text_annotation(fig, text_annot, region, proj);
    if flight_heading:  # draw vectors for flight direction and look direction
        fig = add_flight_vector(fig, flight_heading, look_dir, region)
    if disp_points:
        fig = plot_disp_point_annotations(fig, disp_points);
    if refloc:  # reference location
        fig = plot_reference_annotations(fig, refloc);

    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", frame=["x"+str(label_inc), "y+L\""+label+"\""]);
    fig.savefig(plotname);
    return;
