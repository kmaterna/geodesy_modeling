# Downsample a geocoded interferogram with quadtree;
# Write the outputs to a geojson.
# The script below is the equivalent of fiddling with the file using the gui:
# spool --load unwrap_ll.grd
# Notes: The data files and los files should be in the directory where this is called. 
# Notes for GMTSAR LOS BINARY:
# If there's a .los.enu file in the directory, it will be assumed as the valid one for this experiment
# Be careful what's in your directory. 
# example: gmt grd2xyz unwrap_ll.grd | gmt grdtrack -Gdsdem.grd | awk {'print $1, $2, $4'} |
#      SAT_look S1A20170330_165508_F3.PRM -bos > bot2017.los.enu
# Not used: fig = qt.plot()  # this will open a plot.show(), won't save it. 

from kite import Scene
from Tectonic_Utils.geodesy import geojson2txt
from . import plot_geojson


def kite_downsample_isce_unw(datafile, outname,
                             epislon=1, nan_allowed=0.99, tile_size_min=0.002, tile_size_max=0.010):
    # -------- quadtree downsample an interferogram --------- #
    # epsilon - variance cutoff before the quadtree splits
    # nan_allowed - fraction of pixels that can be nan and still get used
    # tile_size_min - degrees
    # tile_size_max - degrees
    # datafile: .unw.geo file with a matching .xml in the same directory
    # outname: the geojson produced
    # los_rdr_file: los.rdr.geo as produced by isce must be in the same directory
    print("Quadtree Downsampling the file %s into geojson %s " % (datafile, outname));
    sc = Scene.import_data(datafile);
    qt = sc.quadtree
    qt.epsilon = epislon
    qt.nan_allowed = nan_allowed
    qt.tile_size_min = tile_size_min
    qt.tile_size_max = tile_size_max
    qt.export_geojson(outname);
    return;


def geojson_to_outputs(geojsonfile, plotfile, textfile, bbox=(-180, 180, -90, 90), std_min=0.001):
    # ----------------------------
    # This plots downsampled data and standard deviation.
    # It also writes a text file for inversion.
    # The functions here live in the Utility_functions repository
    pixel_list = geojson2txt.read_geojson(geojsonfile);
    plot_geojson.plot_downsampled_InSAR(pixel_list, plotfile, vmin=-120, vmax=20);
    geojson2txt.pixels_to_txt(pixel_list, textfile, bbox, std_min);  # can take a bbox optionally
    return;
