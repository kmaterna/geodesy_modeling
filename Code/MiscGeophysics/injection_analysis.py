# August 2020
# Python scripts to analyze injection at specific wells
# This script creates a pygmt map of injection wells color coded by amount of injection and production
# It also prints out the mismatch between reporting systems: some wells are missing lat/lon information, 
# and some wells are missing volume information 

import numpy as np 
import matplotlib.pyplot as plt 
import collections
import datetime as dt
import xlrd
import csv
import general_python_io_functions as readers

Wells = collections.namedtuple('Wells', ['api', 'lon', 'lat', 'masspermonth', 'dtarray', 'welltype']);


def driver(mass_file, ll_file):
	WellsLatLon = read_doggr_csv_wells(ll_file);
	MassWells = read_well_masses(mass_file);
	TotalWells = combine_masses_wells(WellsLatLon, MassWells);
	map_well_locations(TotalWells);  # reads extra files for annotations as well
	return;

def read_nbgf_injection_totals(filename):
	# Reading a text file with injection and production totals for the NBGF
	ifile = open(filename,'r');
	dtarray = [];
	injection = [];
	production = [];
	for line in ifile:
		temp = line.split();
		if temp[0] == "#":
			continue;
		dtarray.append(dt.datetime.strptime(temp[0]+' '+temp[1], "%Y %m"));
		injection.append(float(temp[3].replace(",", ""))/1000);
		production.append(float(temp[2].replace(",", ""))/1000);
	return dtarray, injection, production;


def read_doggr_csv_wells(filename):
	# Making a list of wells lat/lon from a second format of DOGGR data
	# This will be ALL OF CALIFORNIA
	WellsLatLon = [];
	print("Reading in %s" % filename);
	with open(filename) as csvfile:
		next(csvfile);  # skip the header row
		reader = csv.reader(csvfile);
		for row in reader:
			onewell = Wells(api=int(row[1]), lon=float(row[19]), lat=float(row[18]), masspermonth=None, dtarray=None, welltype=None);
			WellsLatLon.append(onewell);
	return WellsLatLon;


def read_well_masses(filename):
	# Make a bunch of well objects with name, months, and volumes
	# This excel file is from Doggr/CEC-supplied data
	print("Reading in %s" % filename);
	wb = xlrd.open_workbook(filename);
	sheet = wb.sheet_by_index(0);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r, c) for c in range(numcols)] for r in range(numrows)];

	# Get the dates on the left column 
	dates = [data[i][0] for i in range(3, np.shape(data)[0])];
	dtarray = [];
	for i in range(len(dates)):
		dtarray.append(dt.datetime(*xlrd.xldate_as_tuple(dates[i], wb.datemode)))

	# Get the well information, one for each column
	wellArray = [];
	for i in range(2, numcols):
		api = data[1][i];
		injection_type = data[2][i];
		mass = [data[j][i] for j in range(3, numrows)];
		OneWell = Wells(api=api, lon=None, lat=None, masspermonth=mass, dtarray=dtarray, welltype=injection_type);
		wellArray.append(OneWell);
	print("Read masspermonth from %d wells " % (len(wellArray)/2) );
	return wellArray;


def combine_masses_wells(WellsLatLon, MassWells):
	# Associate each production and injection well with its right volume. 
	print("Combining well location and monthly mass information.");
	full_located_wells = [];
	missing_production = 0;
	missing_injection = 0;

	for mass_sample_well in MassWells: 
		found_well = 0;
		for ll_sample_well in WellsLatLon:
	
			if mass_sample_well.api == ll_sample_well.api:  # if there's full information for this well (mass, lat/lon)
				found_well = 1
				temp = Wells(api=mass_sample_well.api, lon=ll_sample_well.lon, lat=ll_sample_well.lat, dtarray=mass_sample_well.dtarray, 
					masspermonth=mass_sample_well.masspermonth, welltype=mass_sample_well.welltype);

				full_located_wells.append(temp);

		# Identify if injection or production wells are missing location
		if found_well == 0 and get_total_production(mass_sample_well.dtarray, mass_sample_well.masspermonth) > 0.001:
			print("Error! Brawley %s Well %s not found in lat/lon database" % (mass_sample_well.welltype, mass_sample_well.api) );
			print("  Accounting for %f volume " % get_total_production(mass_sample_well.dtarray, mass_sample_well.masspermonth));
			temp = get_total_production(mass_sample_well.dtarray, mass_sample_well.masspermonth);
			if mass_sample_well.welltype == 'inject':
				missing_injection = missing_injection+temp;
			else:
				missing_production = missing_production+temp;
	print("Found lon/lat/mass information for %d wells" % (len(full_located_wells)) );
	print("Missing Lat/Lon Injection of %f 1e6" % (missing_injection/1e6));
	print("Missing Lat/Lon Production of %f 1e6" % (missing_production/1e6));

	return full_located_wells; 


def get_total_production(dtarray, masspermonth, 
	starttime=dt.datetime.strptime("20090101", "%Y%m%d"),
	endtime=dt.datetime.strptime("20200101", "%Y%m%d")):
	# Integrate the well total from startdate to enddate
	sum_production = 0;
	for i in range(len(dtarray)):
		if starttime <= dtarray[i] <= endtime:
			temp = masspermonth[i];
			if masspermonth[i] == '':
				temp = 0;
			sum_production = sum_production+temp;
	return sum_production;


def map_well_locations(TotalWells):
	# Makes a pygmt overview map of the wells at NBGF
	import pygmt

	rupture_lon, rupture_lat = np.loadtxt('Data/M4p9_surface_rupture.txt', unpack=True);
	eq_lon, eq_lat = np.loadtxt('../QTM_exploring/Brawley_QTM_shallower8.txt', unpack=True, usecols=(1, 2));
	fault_lon, fault_lat = readers.read_gmt_multisegment_latlon('Data/SwarmFault1.txt_gmt');
	n_fault_lon, n_fault_lat = readers.read_gmt_multisegment_latlon('Data/Wei_normalfault.txt_gmt');
	field_boundaries_lon, field_boundaries_lat = readers.read_gmt_multisegment_latlon('Data/Fields_Boundaries.txt', ',');

	region = [-115.64, -115.42, 32.95, 33.08];  # close-ish
	# region = [-116.0, -114.3, 32.3, 33.5];  # big range
	proj = "M6i"

	fig = pygmt.Figure()
	fig.basemap(region=region, projection=proj, B="+t\"Wells at North Brawley Geothermal Field\"");
	fig.coast(shorelines="1.0p,black", region=region, projection=proj, N=[1, 2], L="g-115.45/33.07+c1.5+w5", B="0.1")
	for i in range(len(field_boundaries_lat)):
		fig.plot(x=field_boundaries_lon[i], y=field_boundaries_lat[i], W='thick,red');
	for i in range(len(fault_lat)):     # strike slip fault
		fig.plot(x=fault_lon[i], y=fault_lat[i], W='thin,gainsboro');
	for i in range(len(n_fault_lat)):   # normal fault
		fig.plot(x=n_fault_lon[i], y=n_fault_lat[i], W='thin,gainsboro');
	fig.plot(x=rupture_lon, y=rupture_lat, W='thick,black')  # M4.7 surface rupture 
	fig.plot(x=eq_lon, y=eq_lat, S='c0.9p', G='darkgray');  # QTM seismicity

	pygmt.makecpt(C="hot", T="0/20.0/0.1", I=True, H="mycpt.cpt");
	
	manywell_production = 0;
	manywell_injection = 0;
	
	for well in TotalWells:
		if well.welltype == 'inject':  # injection wells
			total_production = get_total_production(well.dtarray, well.masspermonth);
			total_production = total_production/1e6;
			if total_production > 0.000:
				fig.plot(x=well.lon, y=well.lat, G=[total_production], S='i0.13i', W='thin,black', C="mycpt.cpt");
			manywell_injection = manywell_injection+total_production;
		if well.welltype == 'product':
			total_production = get_total_production(well.dtarray, well.masspermonth);
			total_production = total_production/1e6;
			if total_production > 0.000:
				fig.plot(x=well.lon, y=well.lat, G=[total_production], S='t0.13i', W='thin,black', C="mycpt.cpt");
			manywell_production = manywell_production+total_production;
	
	fig.legend(region=region, projection=proj, spec="legendfile.txt", D="jBL+o0.2c", F="+gantiquewhite+pthick,black");
	fig.colorbar(D="jBr+w2.2i/0.1i+o1.2c/1.0c+h", C="mycpt.cpt", I="0.8", G="0/20.0", B=["x"+str(4), "y+L\"Vol\""]);
	fig.savefig('well_locations.png')

	print("manywell_injection: %f x 1e6" % manywell_injection);
	print("manywell_production: %f x 1e6" % manywell_production);
	print("From excel well-by-well tally, would expect total injection to be 117");
	print("From excel well-by-well tally, would expect total production to be 130");
	return;

