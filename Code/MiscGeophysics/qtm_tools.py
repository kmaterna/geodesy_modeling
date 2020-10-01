

import numpy as np
import matplotlib.pyplot as plt  
import collections
import datetime as dt
import subprocess
import general_python_io_functions as readers


Catalog = collections.namedtuple("Catalog", ["dtarray", "lon", "lat", "depth", "Mag"]);


def input_qtm(filename):
	# Designed to use the txt file format of Ross et al. (2019)'s QTM catalog
	# downloaded from https://scedc.caltech.edu/research-tools/altcatalogs.html
	print("Reading file %s " % filename);
	print(dt.datetime.now());
	
	# This takes 15 seconds (surprisingly, np.loadtxt takes 30 seconds to do the same thing)
	year = []; month = []; day = []; hour = []; minute = []; second = []; lat = []; lon = []; depth = []; mag = [];
	ifile = open(filename, 'r');
	ifile.readline();
	for line in ifile:
		temp = line.split();
		year.append(temp[0]);
		month.append(temp[1]);
		day.append(temp[2]);
		hour.append(temp[3]);
		minute.append(temp[4]);
		second.append(temp[5]);
		lat.append(float(temp[7]));
		lon.append(float(temp[8]));
		depth.append(float(temp[9]));
		mag.append(float(temp[10]));
	ifile.close();

	dtarray = [];
	for i in range(len(year)):
		try:
			newdate = dt.datetime.strptime(year[i]+month[i]+day[i]+hour[i]+minute[i], "%Y%m%d%H%M");
		except ValueError:  # WE ACTUALLY GOT AN EARTHQUAKE DURING A LEAP SECOND!!!! 
			print("You may have a problem at: ");
			print(year[i]+month[i]+day[i]+hour[i]+minute[i]);
			newdate = dt.datetime.strptime(year[i]+month[i]+day[i]+hour[i], "%Y%m%d%H");
		dtarray.append(newdate);
	MyCat = Catalog(dtarray=dtarray, lon=lon, lat=lat, depth=depth, Mag=mag);
	print("done at : ", dt.datetime.now());
	return MyCat;


def read_simple_catalog_txt(filename):
	# Reading a very simple .txt format for earthquake catalogs
	print("Reading Catalog in %s " % filename);
	[datestrs, lon, lat, depth, Mag] = np.loadtxt(filename, dtype={'names': ('datestr', 'lon', 'lat', 'depth', 'mag'),
																   'formats': ('U19', np.float, np.float, np.float, np.float)
																   }, unpack=True);
	dtarray = [dt.datetime.strptime(i, "%Y-%m-%d-%H-%M-%S") for i in datestrs];
	MyCat = Catalog(dtarray=dtarray, lon=lon, lat=lat, depth=depth, Mag=Mag);
	return MyCat;


def write_simple_catalog_txt(MyCat, outfile):
	# Writing a very simple .txt format for earthquake catalogs
	print("Writing Catalog in %s " % outfile );
	ofile = open(outfile, 'w');
	for i in range(len(MyCat.dtarray)):
		datestr = dt.datetime.strftime(MyCat.dtarray[i], "%Y-%m-%d-%H-%M-%S");
		ofile.write("%s %f %f %.3f %.2f\n" % (datestr, MyCat.lon[i], MyCat.lat[i], MyCat.depth[i], MyCat.Mag[i]) );
	ofile.close();
	return;


def restrict_cat_box(catalog, bbox):
	# A function on earthquake catalogs
	# Limit a catalog to the provided bounding box in lon, lat, depth, and optionally time
	new_dtarray = []; new_lon = []; new_lat = []; new_depth = []; new_Mag = [];
	print("Restricting catalog to box ", bbox)
	for i in range(len(catalog.dtarray)):
		if bbox[0] <= catalog.lon[i] <= bbox[1]:
			if bbox[2] <= catalog.lat[i] <= bbox[3]:
				if bbox[4] <= catalog.depth[i] <= bbox[5]:
					if bbox[6] <= catalog.dtarray[i] <= bbox[7]:
						new_dtarray.append(catalog.dtarray[i]);
						new_lon.append(catalog.lon[i]);
						new_lat.append(catalog.lat[i]);
						new_depth.append(catalog.depth[i]);
						new_Mag.append(catalog.Mag[i]);
	MyCat = Catalog(dtarray=new_dtarray, lon=new_lon, lat=new_lat, depth=new_depth, Mag=new_Mag);
	print("Returning %d out of %d events" % (len(MyCat.depth), len(catalog.depth)) );
	return MyCat;


def make_cumulative_stack(MyCat):
	# Take a catalog
	# Returns two arrays: time and seismicity
	# They can be plotted for a cumulative seismicity plot
	dt_total = []; eq_total = []; adding_sum = 0;
	dt_total.append(MyCat.dtarray[0]);
	eq_total.append(0);
	for i in range(len(MyCat.lon)):	
		dt_total.append(MyCat.dtarray[i]);
		eq_total.append(adding_sum);
		adding_sum = adding_sum+1;
		eq_total.append(adding_sum);
		dt_total.append(MyCat.dtarray[i]);
	return dt_total, eq_total;


def mapping_plot(MyCat, boundary_file, outdir):
	boundary_lons, boundary_lats = readers.read_gmt_multisegment_latlon(boundary_file, ',');
	fault_lon, fault_lat = readers.read_gmt_multisegment_latlon('../Injection_Data/Data/SwarmFault1.txt_gmt');
	n_fault_lon, n_fault_lat = readers.read_gmt_multisegment_latlon('../Injection_Data/Data/Wei_normalfault.txt_gmt');
	fig = plt.figure(figsize=(14, 12), dpi=300);
	plt.scatter(MyCat.lon, MyCat.lat, s=MyCat.Mag, c=MyCat.depth, cmap='viridis');
	for i in range(len(boundary_lons)):
		plt.plot(boundary_lons[i], boundary_lats[i], marker=None, linewidth=1, color='indianred');
	for i in range(len(fault_lon)):
		plt.plot(fault_lon[i], fault_lat[i], marker=None, linewidth=0.5, color='gray');
	for i in range(len(n_fault_lon)):
		plt.plot(n_fault_lon[i], n_fault_lat[i], marker=None, linewidth=0.5, color='gray');
	plt.colorbar();
	plt.title("QTM Catalog: %d events " % len(MyCat.lon), fontsize=20 );
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	plt.xlabel("Longitude", fontsize=18);
	plt.ylabel("Latitude", fontsize=18);
	plt.xlim([-115.66, -115.43]);
	plt.ylim([32.9, 33.1]);
	plt.savefig(outdir+"/QTM_events.png");
	return;

def make_lollipop_plot(MyCat):
	plt.figure(dpi=300, figsize=(10,7));
	for i in range(len(MyCat.dtarray)):
		plt.plot(MyCat.dtarray[i], MyCat.Mag[i], marker='o', markersize=10, linewidth=0,color='black');
		plt.plot([MyCat.dtarray[i], MyCat.dtarray[i]],[0,MyCat.Mag[i]], color='black',linewidth=1);
	plt.ylabel('Magnitude',fontsize=20);
	plt.xlabel('Time',fontsize=20);
	plt.gca().tick_params(axis='both', which='major', labelsize=16)
	plt.ylim([2.5, 5.0])
	plt.savefig('lollipop.png');
	return;

def add_Brawley_timeseries_annotations(ax):
	# For Brawley
	[top, bottom] = ax.get_ylim();
	emc = dt.datetime.strptime("2010-04-03", "%Y-%m-%d");
	swarm = dt.datetime.strptime("2012-08-26", "%Y-%m-%d");
	s_14 = dt.datetime.strptime("2014-06-01", "%Y-%m-%d");
	e_14 = dt.datetime.strptime("2014-10-15", "%Y-%m-%d");
	ax.plot([emc, emc], [bottom, top], '--k');
	ax.plot([swarm, swarm], [bottom, top], '--k');
	ax.plot([s_14, s_14], [bottom, top], '--g');
	ax.plot([e_14, e_14], [bottom, top], '--g');
	return ax;


def cumulative_seismicity_plot(MyCat, outdir):
	dt_total, eq_total = make_cumulative_stack(MyCat);
	fig = plt.figure(figsize=(18, 9), dpi=300);
	plt.plot_date(dt_total, eq_total, linewidth=2, linestyle='solid', marker=None, color='blue');
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	ax = add_Brawley_timeseries_annotations(plt.gca());
	plt.xlabel("Time", fontsize=18);
	plt.ylabel("Cumulative Earthquakes", fontsize=18);
	plt.title("Cumulative Seismicity In Brawley", fontsize=20);
	plt.savefig(outdir+"/StackedSeismicity.png");
	return;


def main_function(catfilename, field_filename, bbox, label):
	outdir = 'Steps/'+label+'_depth_'+str(bbox[4])+'_'+str(bbox[5])+'_'+str(bbox[2])+'_'+str(bbox[3]);
	outdir = outdir+'_'+dt.datetime.strftime(bbox[6], "%Y%m%d")+'_'+dt.datetime.strftime(bbox[7], '%Y%m%d');  # making a nice outdir
	subprocess.call(['mkdir', '-p', outdir], shell=False);

	# Inputs
	MyCat = input_qtm(catfilename);
	
	# Compute and Output
	BrawleyCat = restrict_cat_box(MyCat, bbox);
	mapping_plot(BrawleyCat, field_filename, outdir);
	cumulative_seismicity_plot(BrawleyCat, outdir);
	write_simple_catalog_txt(BrawleyCat, outdir+"/Brawley_QTM.txt");
	return;


if __name__ == "__main__":
	catfilename = filename = "../../../Misc/General_Data/QTM/qtm_final_12dev.hypo"
	field_filename = "../_Project_Data/TRE_Data/CEC_Data/Data/DOGGR_GIS/Fields_Boundaries.txt"
	starttime = dt.datetime.strptime('20080101', "%Y%m%d");
	endtime = dt.datetime.strptime('20220101', "%Y%m%d");
	bbox = [-115.66, -115.43, 32.9, 33.1, 0, 8, starttime, endtime];  # lon, lat, depth, time ranges  # general Brawley Area
	label = 'test'
	main_function(catfilename, field_filename, bbox, label);

