"""
# Write and output functions for InSAR 2D data format
"""

import pygmt
import numpy as np
from Tectonic_Utils.read_write import netcdf_read_write
from . import inputs

def write_InSAR2D(InSAR_Obj, filename):
    netcdf_read_write.produce_output_netcdf(InSAR_Obj.lon, InSAR_Obj.lat, InSAR_Obj.LOS, '', filename);
    return;

def plot_wrapped_insar(grd_filename, plotname):
    print("Mapping file %s " % grd_filename);
    InSAR_2D_Obj = inputs.inputs_grd(grd_filename);
    proj = 'M4i'
    region = [np.min(InSAR_2D_Obj.lon), np.max(InSAR_2D_Obj.lon), np.min(InSAR_2D_Obj.lat), np.max(InSAR_2D_Obj.lat)];
    # Build a PyGMT plot
    fig = pygmt.Figure();
    pygmt.makecpt(C="rainbow", T="-3.14/3.14/0.01", D="o", H="mycpt.cpt");
    fig.basemap(region=region, projection=proj, B="+t\"Wrapped LOS Displacement\"");
    fig.grdimage(grd_filename, region=region, C="mycpt.cpt");
    fig.coast(region=region, projection=proj, N='1', W='1.0p,black', S='lightblue',
              L="n0.23/0.06+c" + str(region[2]) + "+w20", B="1.0");
    fig.colorbar(D="JCR+w4.0i+v+o0.7i/0i", C="mycpt.cpt", G="-3.14/3.14", B=["x1.57", "y+L\"Phase\""]);
    fig.savefig(plotname);
    return;
