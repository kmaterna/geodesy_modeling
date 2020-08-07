# August 2020
# Python scripts to analyze injection at specific wells
# The sum of injected fluid is 2x higher than extracted fluid in 2009-2019.  Might want to check that out. 

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import datetime as dt
import xlrd
import pygmt

Wells = collections.namedtuple('Wells',['api','lon','lat','masspermonth','dtarray','welltype']);


def driver(file_injection, file_production, file_mass, fields_file):
	InjectionWells = read_excel_wells(file_injection);
	ProductionWells = read_excel_wells(file_production);
	field_boundaries_lon, field_boundaries_lat = read_gmt_multisegment(fields_file,',');
	MassWells = read_well_masses(file_mass);
	InjectionWells, ProductionWells = combine_masses_wells(InjectionWells, ProductionWells, MassWells);
	map_well_locations(InjectionWells, ProductionWells, field_boundaries_lon, field_boundaries_lat);
	return;

def read_gmt_multisegment(fields_file, split_delimiter=' '):
	print("reading gmt multisegment file %s" % fields_file);
	ifile=open(fields_file);
	lon_collection=[];
	lat_collection=[];
	lon_temp=[];
	lat_temp=[];
	for line in ifile:
		if line.split()[0]=='>>' or line.split()[0]=='>':
			if lon_temp !=[]:
				lon_collection.append(lon_temp);
				lat_collection.append(lat_temp);
			lon_temp=[];
			lat_temp=[];
			continue;
		else:
			temp=line.split(split_delimiter);
			lon_temp.append(float(temp[0]));
			lat_temp.append(float(temp[1]));
	lon_collection.append(lon_temp);
	lat_collection.append(lat_temp);
	return lon_collection, lat_collection;

def read_excel_wells(filename):
	print("Reading in %s" % filename);
	wb = xlrd.open_workbook(filename);
	sheet = wb.sheet_by_index(0);
	names = sheet.col_values(0)[1:];
	latitudes = sheet.col_values(1)[1:];
	longitudes = sheet.col_values(2)[1:];
	Well_collection = [];
	for i in range(len(longitudes)):  # make a list of wells
		New_Well = Wells(api=names[i], lon=longitudes[i], lat=latitudes[i], masspermonth=None, dtarray=None, welltype=None);
		Well_collection.append(New_Well);
	return Well_collection;

def read_well_masses(filename):
	# Make a bunch of well objects with name, months, and volumes
	print("Reading in %s" % filename);
	wb = xlrd.open_workbook(filename);
	sheet = wb.sheet_by_index(0);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];

	# Get the dates on the left column 
	dates = [data[i][0] for i in range(3,np.shape(data)[0])];
	dtarray=[];
	for i in range(len(dates)):
		dtarray.append(dt.datetime(*xlrd.xldate_as_tuple(dates[i], wb.datemode)))

	# Get the well information, one for each column
	wellArray = [];
	for i in range(2,numcols):
		api = data[1][i];
		injection_type = data[2][i];
		mass = [data[j][i] for j in range(3,numrows)];
		OneWell = Wells(api=api, lon=None, lat=None, masspermonth=mass, dtarray=dtarray, welltype=injection_type);
		wellArray.append(OneWell);
	print("Read masspermonth from %d wells " % (len(wellArray)/2) );
	return wellArray;

def combine_masses_wells(InjectionWells, ProductionWells, MassWells):
	# Associate each production and injection well with its right volume. 
	print("Combining well location and monthly mass information.");
	new_injection_wells = []; 
	new_production_wells = [];
	for well in InjectionWells:
		found_well_mass = 0;
		for sample_well in MassWells:
			if sample_well.api == well.api:
				if sample_well.welltype=="inject":
					found_well_mass=1
					temp = Wells(api=well.api, lon=well.lon, lat=well.lat, dtarray=sample_well.dtarray, 
						masspermonth=sample_well.masspermonth, welltype=sample_well.welltype);
					new_injection_wells.append(temp);
		if found_well_mass==0 and well.lat>32.8 and well.lat<33.1 and well.lon<-115.4 and well.lon>-115.7:
			print("Error! Could not find mass information for Brawley injection well %d at %f, %f " % (well.api,well.lon, well.lat) );
	print("Found lon/lat/mass information for %d injection wells" % (len(new_injection_wells)) );
	
	for well in ProductionWells:
		found_well_mass=0;
		for sample_well in MassWells:
			if sample_well.api == well.api:
				if sample_well.welltype=="product":
					found_well_mass=1
					temp = Wells(api=well.api, lon=well.lon, lat=well.lat, dtarray=sample_well.dtarray, 
						masspermonth=sample_well.masspermonth, welltype=sample_well.welltype);
					new_production_wells.append(temp);
		if found_well_mass==0 and well.lat>32.8 and well.lat<33.1 and well.lon<-115.4 and well.lon>-115.7:
			print("Error! Could not find mass information for Brawley production well %d at %f, %f " % (well.api, well.lon, well.lat) );
	print("Found lon/lat/mass information for %d production wells" % (len(new_production_wells)) );	

	# Identify missing wells (no lat/lon information)
	missing_production=0; 
	missing_injection =0;
	api_lonlat = [x.api for x in InjectionWells];
	for sample_well in MassWells:
		if sample_well.welltype=='inject' and sample_well.api not in api_lonlat:
			temp = get_total_production(sample_well.dtarray, sample_well.masspermonth);
			if temp > 0.0001:
				print("Error! Brawley Injection Well %s not found in lat/lon database" % (sample_well.api) );
				print("  Accounting for %f injection " % get_total_production(sample_well.dtarray, sample_well.masspermonth));
				missing_injection = missing_injection+temp;
	api_lonlat = [x.api for x in ProductionWells];
	for sample_well in MassWells:
		if sample_well.welltype=='product' and sample_well.api not in api_lonlat:
			temp = get_total_production(sample_well.dtarray, sample_well.masspermonth);
			if temp > 0.0001: 
				print("Error! Brawley Production Well %s not found in lat/lon database" % (sample_well.api) );
				print("  Accounting for %f production " % get_total_production(sample_well.dtarray, sample_well.masspermonth));
				missing_production = missing_production+temp;

	print("Missing Injection is %f 1e6" % (missing_injection/1e6));
	print("Missing Production is %f 1e6" % (missing_production/1e6));

	return new_injection_wells, new_production_wells;


def get_total_production(dtarray, masspermonth, 
	starttime=dt.datetime.strptime("20090101","%Y%m%d"), 
	endtime=dt.datetime.strptime("20200101","%Y%m%d")):
	sum_production = 0;
	for i in range(len(dtarray)):
		if dtarray[i]>=starttime and dtarray[i]<=endtime:
			temp = masspermonth[i];
			if masspermonth[i]=='':
				temp=0;
			sum_production=sum_production+temp;
	return sum_production;


def map_well_locations(InjWells, ProWells, field_boundaries_lon, field_boundaries_lat):

	rupture_lon, rupture_lat = np.loadtxt('Data/M4p9_surface_rupture.txt',unpack=True);
	eq_lon, eq_lat = np.loadtxt('../QTM_exploring/Brawley_QTM.txt',unpack=True,usecols=(1,2));
	fault_lon, fault_lat = read_gmt_multisegment('Data/SwarmFault1.txt_gmt');
	n_fault_lon, n_fault_lat = read_gmt_multisegment('Data/Wei_normalfault.txt_gmt');

	region = [-115.64, -115.42, 32.95, 33.08];  # close-ish
	# region = [-116.0, -114.3, 32.3, 33.5];  # big range
	proj = "M6i"

	fig = pygmt.Figure()
	fig.basemap(region=region,projection=proj,B="+t\"Wells at North Brawley Geothermal Field\"");
	fig.coast(shorelines="1.0p,black",region=region,projection=proj,N=[1,2],L="g-115.45/33.07+c1.5+w5",B="0.1")
	for i in range(len(field_boundaries_lat)):
		fig.plot(x=field_boundaries_lon[i], y=field_boundaries_lat[i], W='thick,red');
	for i in range(len(fault_lat)):     # strike slip fault
		fig.plot(x=fault_lon[i], y=fault_lat[i], W='thin,gainsboro');
	for i in range(len(n_fault_lat)):   # normal fault
		fig.plot(x=n_fault_lon[i], y=n_fault_lat[i], W='thin,gainsboro');
	fig.plot(x=rupture_lon, y=rupture_lat, W='thick,black')  # M4.7 surface rupture 
	fig.plot(x=eq_lon, y=eq_lat, S='c0.9p',G='darkgray');  # QTM seismicity

	pygmt.makecpt(C="hot",T="0/20.0/0.1",I=True,H="mycpt.cpt");
	
	manywell_production=0;
	manywell_injection=0;
	for well in InjWells:
		total_production = get_total_production(well.dtarray, well.masspermonth);
		total_production = total_production/1e6;
		if total_production > 0.000:
			# fig.plot(x=well.lon, y=well.lat, S='i0.1i',G='royalblue2',W='thin,black');
			fig.plot(x=well.lon, y=well.lat, G=[total_production], S='i0.13i',W='thin,black',C="mycpt.cpt");
		manywell_injection=manywell_injection+total_production;
	for well in ProWells: # production wells
		total_production = get_total_production(well.dtarray, well.masspermonth);
		total_production = total_production/1e6;
		if total_production > 0.000:
			fig.plot(x=well.lon, y=well.lat, G=[total_production], S='t0.13i',W='thin,black',C="mycpt.cpt");
		manywell_production=manywell_production+total_production;
	
	fig.legend(region=region, projection=proj,spec="legendfile.txt",D="jBL+o0.2c",F="+gantiquewhite+pthick,black");
	fig.colorbar(D="jBr+w2.2i/0.1i+o1.2c/1.0c+h",C="mycpt.cpt",I="0.8",G="0/20.0",B=["x"+str(4),"y+L\"Vol\""]); 
	fig.savefig('well_locations.png')

	print("manywell_injection: %f x 1e6" % manywell_injection);
	print("manywell_production: %f x 1e6" % manywell_production);
	print("From excel, would expect total injection to be 117");
	print("From excel, would expect total production to be 130");
	return;

