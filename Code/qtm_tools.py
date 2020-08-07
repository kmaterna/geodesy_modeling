

import numpy as np
import matplotlib.pyplot as plt  
import collections
import datetime as dt
import subprocess
import multiSAR_input_functions


Catalog = collections.namedtuple("Catalog",["dtarray","lon","lat","depth","Mag"]);


def input_qtm(filename):
	print("Reading file %s " % filename);
	print(dt.datetime.now());
	
	# This takes 15 seconds
	year=[]; month=[]; day=[]; hour=[]; minute=[]; second=[]; lat=[]; lon=[]; depth=[]; mag=[];
	ifile=open(filename,'r');
	ifile.readline();
	for line in ifile:
		temp=line.split();
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


	# # This takes 30 seconds 
	# year, month, day, hour, minute, second, lat, lon, depth, mag = np.loadtxt(filename, skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 7, 8, 9, 10), unpack=True,
	# 	dtype={'names':('year','month','day','hour','minute','second','lat','lon','depth','mag'),'formats':('U4','U2','U2','U2','U2',np.float, np.float, np.float, np.float, np.float)} );


	dtarray=[];
	for i in range(len(year)):
		try:
			newdate = dt.datetime.strptime(year[i]+month[i]+day[i]+hour[i]+minute[i],"%Y%m%d%H%M");
		except ValueError:  # WE ACTUALLY GOT AN EARTHQUAKE DURING A LEAP SECOND!!!! 
			print("You may have a problem at: ");
			print(year[i]+month[i]+day[i]+hour[i]+minute[i]);
			newdate = dt.datetime.strptime(year[i]+month[i]+day[i]+hour[i],"%Y%m%d%H");
		dtarray.append(newdate);
	MyCat = Catalog(dtarray=dtarray, lon=lon, lat=lat, depth=depth, Mag=mag);
	print("done at : ", dt.datetime.now());
	return MyCat;

def restrict_cat_box(catalog, bbox):
	# Limit a catalog to the provided bounding box in lon, lat, depth, and optionally time
	new_dtarray=[]; new_lon=[]; new_lat=[]; new_depth=[]; new_Mag=[]; 
	if len(bbox)==6:
		starttime=min(catalog.dtarray);
		endtime=max(catalog.dtarray);
	else:
		starttime=bbox[6];
		endtime=bbox[7];
	print("Restricting catalog to box ", bbox, starttime, endtime)		
	for i in range(len(catalog.dtarray)):
		if catalog.lon[i]>=bbox[0] and catalog.lon[i]<=bbox[1]:
			if catalog.lat[i]>=bbox[2] and catalog.lat[i]<=bbox[3]:
				if catalog.depth[i]>=bbox[4] and catalog.depth[i]<=bbox[5]:
					if catalog.dtarray[i]>=starttime and catalog.dtarray[i]<=endtime:
						new_dtarray.append(catalog.dtarray[i]);
						new_lon.append(catalog.lon[i]);
						new_lat.append(catalog.lat[i]);
						new_depth.append(catalog.depth[i]);
						new_Mag.append(catalog.Mag[i]);
	MyCat = Catalog(dtarray=new_dtarray, lon=new_lon, lat=new_lat, depth=new_depth, Mag=new_Mag);
	print("Returning %d out of %d events" % (len(MyCat.depth), len(catalog.depth)) );
	return MyCat;

def make_cumulative_stack(MyCat):
	dt_total = []; eq_total=[]; adding_sum=0;
	dt_total.append(MyCat.dtarray[0]);
	eq_total.append(0);
	for i in range(len(MyCat.lon)):	
		dt_total.append(MyCat.dtarray[i]);
		eq_total.append(adding_sum);
		adding_sum=adding_sum+1;
		eq_total.append(adding_sum);
		dt_total.append(MyCat.dtarray[i]);
	return dt_total, eq_total;

def mapping_plot(MyCat, boundary_file, outdir):
	boundary_lons, boundary_lats = multiSAR_input_functions.read_gmt_multisegment_latlon(boundary_file,',');
	fault_lon, fault_lat = multiSAR_input_functions.read_gmt_multisegment_latlon('../Injection_Data/Data/SwarmFault1.txt_gmt');
	n_fault_lon, n_fault_lat = multiSAR_input_functions.read_gmt_multisegment_latlon('../Injection_Data/Data/Wei_normalfault.txt_gmt');	
	fig = plt.figure(figsize=(14,12), dpi=300);
	plt.scatter(MyCat.lon, MyCat.lat, s=MyCat.Mag, c=MyCat.depth, cmap='viridis');
	for i in range(len(boundary_lons)):
		plt.plot(boundary_lons[i], boundary_lats[i], marker=None, linewidth=1, color='indianred');
	for i in range(len(fault_lon)):
		plt.plot(fault_lon[i], fault_lat[i], marker=None, linewidth=0.5, color='gray');
	for i in range(len(n_fault_lon)):
		plt.plot(n_fault_lon[i], n_fault_lat[i], marker=None, linewidth=0.5, color='gray');
	plt.colorbar();
	plt.title("QTM Catalog: %d events " % len(MyCat.lon),fontsize=20 );
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	plt.xlabel("Longitude",fontsize=18);
	plt.ylabel("Latitude",fontsize=18);
	plt.xlim([-115.66, -115.43]);
	plt.ylim([32.9, 33.1]);
	plt.savefig(outdir+"/QTM_events.png");
	return;

def add_Brawley_annotations(ax):
	# For Brawley
	[top,bottom]=ax.get_ylim();
	emc = dt.datetime.strptime("2010-04-03","%Y-%m-%d");
	swarm = dt.datetime.strptime("2012-08-26","%Y-%m-%d");
	s_14 = dt.datetime.strptime("2014-06-01","%Y-%m-%d");
	e_14 = dt.datetime.strptime("2014-10-15","%Y-%m-%d");
	ax.plot([emc, emc],[bottom, top],'--k');
	ax.plot([swarm, swarm],[bottom, top],'--k');	
	ax.plot([s_14, s_14],[bottom, top],'--g');
	ax.plot([e_14, e_14],[bottom, top],'--g');
	return ax;

def cumulative_seismicity_plot(MyCat, outdir):
	dt_total, eq_total = make_cumulative_stack(MyCat);
	fig=plt.figure(figsize=(18,9),dpi=300);
	plt.plot_date(dt_total, eq_total, linewidth=2, linestyle='solid',marker=None, color='blue');
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	ax = add_Brawley_annotations(plt.gca());
	plt.xlabel("Time",fontsize=18);
	plt.ylabel("Cumulative Earthquakes",fontsize=18);
	plt.title("Cumulative Seismicity In Brawley",fontsize=20);
	plt.savefig(outdir+"/StackedSeismicity.png");
	return;

def write_catalog(MyCat, outfile):
	print("Writing Catalog in %s " % (outfile) );
	ofile=open(outfile,'w');
	for i in range(len(MyCat.dtarray)):
		datestr=dt.datetime.strftime(MyCat.dtarray[i],"%Y-%m-%d-%H-%M-%S");
		ofile.write("%s %f %f %.3f %.2f\n" % (datestr, MyCat.lon[i], MyCat.lat[i], MyCat.depth[i], MyCat.Mag[i]) );
	ofile.close();
	return;


def main_function(catfilename, field_filename, bbox):
	outdir='Steps/depth_'+str(bbox[4])+'_'+str(bbox[5])+'_'+str(bbox[2])+'_'+str(bbox[3]);
	if len(bbox)>6:
		outdir=outdir+'_'+dt.datetime.strftime(bbox[6],"%Y%m%d")+'_'+dt.datetime.strftime(bbox[7],'%Y%m%d');
	subprocess.call(['mkdir','-p',outdir],shell=False);

	# Inputs
	MyCat = input_qtm(catfilename);
	
	# Compute and Output
	BrawleyCat = restrict_cat_box(MyCat, bbox);
	mapping_plot(BrawleyCat, field_filename, outdir);
	cumulative_seismicity_plot(BrawleyCat, outdir);
	write_catalog(BrawleyCat,outdir+"/Brawley_QTM.txt");
	return;


if __name__=="__main__":
	catfilename = filename = "../../../Misc/General_Data/QTM/qtm_final_12dev.hypo"
	field_filename = "../_Project_Data/TRE_Data/CEC_Data/Data/DOGGR_GIS/Fields_Boundaries.txt"
	# bbox = [-115.66, -115.43, 33.0, 33.04, 0, 6]; # lon, lat, depth ranges, NBGF-specific
	bbox = [-115.66, -115.43, 32.9, 33.1, 0, 8]; # lon, lat, depth ranges  # general Brawley Area	
	main_function(catfilename, field_filename, bbox);


