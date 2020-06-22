

import numpy as np
import matplotlib.pyplot as plt  
import collections
import datetime as dt

Catalog = collections.namedtuple("Catalog",["dtarray","lon","lat","depth","Mag"]);

def configure():
	filename = "../../../Misc/General_Data/QTM/qtm_final_12dev.hypo"
	field_filename = "../_Project_Data/TRE_Data/CEC_Data/Data/DOGGR_GIS/Fields_Boundaries.txt"
	# bbox = [-115.66, -115.43, 33.0, 33.04, 0, 6]; # lon, lat, depth ranges, NBGF-specific
	bbox = [-115.66, -115.43, 32.9, 33.1, 0, 8]; # lon, lat, depth ranges  # general Brawley Area
	return filename, field_filename, bbox;

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
	# year, month, day, hour, minute, second, lat, lon, depth, mag = np.loadtxt(filename, 
	# 	skiprows=1, usecols=(0, 1, 2, 3, 4, 5, 7, 8, 9, 10), unpack=True,
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

def read_fields(filename):
	ifile=open(filename,'r');
	boundary_lons=[]; boundary_lats=[];
	for line in ifile:
		if line.split()[0]==">>":
			boundary_lons.append(np.nan);
			boundary_lats.append(np.nan);
			continue;
		boundary_lons.append(float(line.split(',')[0]));
		boundary_lats.append(float(line.split(',')[1]));
	ifile.close();
	return boundary_lons, boundary_lats;

def restrict_cat_box(catalog, bbox):
	# Limit a catalog to the provided bounding box in lon, lat, depth
	new_dtarray=[]; new_lon=[]; new_lat=[]; new_depth=[]; new_Mag=[]; 
	for i in range(len(catalog.dtarray)):
		if catalog.lon[i]>=bbox[0] and catalog.lon[i]<=bbox[1]:
			if catalog.lat[i]>=bbox[2] and catalog.lat[i]<=bbox[3]:
				if catalog.depth[i]>=bbox[4] and catalog.depth[i]<=bbox[5]:
					new_dtarray.append(catalog.dtarray[i]);
					new_lon.append(catalog.lon[i]);
					new_lat.append(catalog.lat[i]);
					new_depth.append(catalog.depth[i]);
					new_Mag.append(catalog.Mag[i]);
	MyCat = Catalog(dtarray=new_dtarray, lon=new_lon, lat=new_lat, depth=new_depth, Mag=new_Mag);
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

def mapping_plot(MyCat, boundary_lons, boundary_lats):
	fig = plt.figure(figsize=(14,12), dpi=300);
	plt.scatter(MyCat.lon, MyCat.lat, s=MyCat.Mag, c=MyCat.depth, cmap='viridis');
	plt.plot(boundary_lons, boundary_lats, marker=None, linewidth=1, color='indianred');
	plt.colorbar();
	plt.title("QTM Catalog: %d events " % len(MyCat.lon),fontsize=20 );
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	plt.xlabel("Longitude",fontsize=18);
	plt.ylabel("Latitude",fontsize=18);
	plt.xlim([-115.66, -115.43]);
	plt.ylim([32.9, 33.1]);
	plt.savefig("QTM_12dev.png");
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

def cumulative_seismicity_plot(MyCat):
	dt_total, eq_total = make_cumulative_stack(MyCat);
	print(len(dt_total));
	print(len(eq_total));
	fig=plt.figure(figsize=(18,9),dpi=300);
	plt.plot_date(dt_total, eq_total, linewidth=2, linestyle='solid',marker=None, color='blue');
	plt.gca().tick_params(axis='both', which='major', labelsize=16);
	ax = add_Brawley_annotations(plt.gca());
	plt.xlabel("Time",fontsize=18);
	plt.ylabel("Cumulative Earthquakes",fontsize=18);
	plt.title("Cumulative Seismicity In Brawley",fontsize=20);
	plt.savefig("StackedSeismicity.png");
	return;

def write_catalog(MyCat, outfile):
	print("Writing Catalog in %s " % (outfile) );
	ofile=open(outfile,'w');
	for i in range(len(MyCat.dtarray)):
		datestr=dt.datetime.strftime(MyCat.dtarray[i],"%Y-%m-%d-%H-%M-%S");
		ofile.write("%s %f %f %.3f %.2f\n" % (datestr, MyCat.lon[i], MyCat.lat[i], MyCat.depth[i], MyCat.Mag[i]) );
	ofile.close();
	return;

if __name__=="__main__":
	filename, field_filename, bbox =configure();
	MyCat = input_qtm(filename);
	boundary_lons, boundary_lats = read_fields(field_filename);
	BrawleyCat = restrict_cat_box(MyCat, bbox);
	mapping_plot(BrawleyCat, boundary_lons, boundary_lats);
	cumulative_seismicity_plot(BrawleyCat);
	write_catalog(BrawleyCat,"Brawley_QTM.txt");



