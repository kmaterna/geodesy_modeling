# Downsample a geocoded interferogram with quadtree;
# Write the outputs to a geojson.
# The script below is the equivalent of fiddling with the file using the gui:
# spool --load unwrap_ll.grd
# Notes: The data files and los files should be in the directory where this is called. 
# Notes for GMTSAR LOS BINARY:
# If there's a .los.enu file in the directory, it will be assumed as the valid one for this experiment
# Be careful what's in your directory. 
# example: gmt grd2xyz unwrap_ll.grd | gmt grdtrack -Gdsdem.grd | awk {'print $1, $2, $4'} | SAT_look S1A20170330_165508_F3.PRM -bos > bot2017.los.enu

import numpy as np 
import netcdf_read_write
import subprocess
import matplotlib.pyplot as plt 
from kite import Scene
import geojson2txt

# ----------------------------
# This downsamples an interferogram.
# sc = Scene.import_data(datafile);
# qt = sc.quadtree; set parameters...;
# fig = qt.plot()  # this will open a plot.show(), won't save it. 
# qt.export_geojson(outdir+"/programmatic.geojson");

def kite_downsample_isce_unw(datafile, los_rdr_file, outname):
	# datafile: some .unw.geo file with a matching .xml in the same directory
	# los_rdr_file: los.rdr.geo as produced by isce
	# outname: the geojson produced
	# If there's not matching .xml, it should fail. 
	# The kite parameters should eventually be fed in. 
	print("Quadtree Downsampling the file %s into geojson %s " % (datafile, outname) );
	sc = Scene.import_data(datafile);
	qt = sc.quadtree
	qt.epsilon = 8
	qt.nan_allowed = 0.99
	qt.tile_size_min = 0.002
	qt.tile_size_max = 0.009
	qt.export_geojson(outname);
	return;

def geojson_to_outputs(geojsonfile, plotfile, textfile):
	# ----------------------------
	# This plots downsampled data and standard deviation.
	# It also writes a text file for inversion.
	bbox = [-115.75, -115.48, 32.98, 33.09];
	pixel_list = geojson2txt.read_geojson(geojsonfile);
	geojson2txt.plot_downsampled_InSAR(pixel_list,plotfile,vmin=-120, vmax=20);
	geojson2txt.pixels_to_txt(pixel_list, textfile, bbox);  # can take a bbox optionally
	return;
