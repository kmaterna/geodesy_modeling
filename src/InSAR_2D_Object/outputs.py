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


def map_wrapped_insar(grd_filename, plotname, text_annot=None, flight_heading=None, look_dir="right"):
    """
    Pygmt map of wrapped InSAR deformation with optional annotations
    """
    print("Mapping file %s " % grd_filename);
    InSAR_2D_Obj = inputs.inputs_grd(grd_filename);
    proj = 'M4i'
    region = [np.min(InSAR_2D_Obj.lon), np.max(InSAR_2D_Obj.lon), np.min(InSAR_2D_Obj.lat), np.max(InSAR_2D_Obj.lat)];
    # Build a PyGMT plot
    fig = pygmt.Figure();
    pygmt.makecpt(cmap="rainbow", series="-3.14/3.14/0.01", background="o", output="mycpt.cpt");
    fig.basemap(region=region, projection=proj, frame="+t\"Wrapped LOS Displacement\"");
    fig.grdimage(grd_filename, region=region, cmap="mycpt.cpt");
    fig.coast(region=region, projection=proj, borders='1', shorelines='1.0p,black', water='lightblue',
              map_scale="n0.43/0.06+c" + str(region[2]) + "+w20", frame="1.0");
    if text_annot is not None:
        fig.text(position='TL', text=text_annot, region=region, projection=proj, font="14p,Helvetica,black",
                 offset="0.2/-0.2", pen="0.5p,black");
    if flight_heading is not None:  # draw vectors for flight direction and look direction
        [x_flight, y_flight] = insar_vector_functions.get_unit_vector_from_heading(flight_heading);
        if look_dir == 'right':
            [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading+90);
        else:
            [x_los, y_los] = insar_vector_functions.get_unit_vector_from_heading(flight_heading-90);
        fig.text(x=region[0], y=region[2], text='LOS', offset=str(0.6 + 0.4*x_los)+"i/0.45i",
                 fill='white', font="10p,Helvetica,black");  # LOS text annotation
        fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
                 direction=[[x_flight], [y_flight]], pen="thin,black", offset="0.6i/0.45i");  # flight vector
        fig.plot(x=[region[0]], y=[region[2]], style='v0.2c+e+gblack+h0+p1p,black+z' + str(1.0),
                 direction=[[x_los/2], [y_los/2]], pen="thin,black", offset="0.6i/0.45i");  # los vector
    fig.colorbar(position="JCR+w4.0i+v+o0.7i/0i", cmap="mycpt.cpt", truncate="-3.14/3.14",
                 frame=["x1.57", "y+L\"Phase\""]);
    fig.savefig(plotname);
    return;
