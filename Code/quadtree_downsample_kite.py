# A test to read a geocoded interferogram
# Use quadtree downsampling
# Write the outputs to a geojson
# The script below is the equivalent of fiddling with the file using the gui:
# spool --load unwrap_ll.grd
# Notes: The data files and los files should be in the directory 
# where this is called. 
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
# def setup():
# 	i=6;
# 	working_dir = "isce_geocode/scene_"+str(i)
# 	working_file = "ts_slice_"+str(i)+".unw.geo";
# 	subprocess.call(["cp","isce_geocode/los.rdr.geo","."],shell=False);
# 	datafile=working_dir+"/"+working_file
# 	subprocess.call(["cp",datafile,"."],shell=False);
# 	subprocess.call(["cp",datafile+".xml","."],shell=False);
# 	return working_file, working_dir;

# datafile, outdir = setup();

# sc = Scene.import_data(datafile);
# qt = sc.quadtree
# qt.epsilon = 1
# qt.nan_allowed = 0.99
# qt.tile_size_min = 0.00025
# qt.tile_size_max = 0.001

# fig = qt.plot()  # this will open a plot.show(), won't save it. 
# qt.export_geojson(outdir+"/programmatic.geojson");

# subprocess.call(["rm",datafile],shell=False);
# subprocess.call(["rm",datafile+".xml"],shell=False);
# subprocess.call(["rm","los.rdr.geo"],shell=False);

# ----------------------------
# This plots and writes a text file
# Once you've downsampled, you can plot and write the text file. 
folder = "isce_geocode/scene_1"
datafile = folder+"/programmatic.geojson"
plotfile = folder+"/downsampled.png"
textfile = folder+"/insar_disps.txt"
bbox = [-115.59, -115.47, 32.98, 33.06];
pixel_list = geojson2txt.read_geojson(datafile);
geojson2txt.plot_downsampled_InSAR(pixel_list,plotfile,vmax=120);
geojson2txt.pixels_to_txt(pixel_list, textfile, bbox);


