# A series of functions to  input various data types
# Leveling
# UAVSAR
# TSX
# S1
# Envisat

import numpy as np 
import datetime as dt 
import collections
import xlrd
import netcdf_read_write

# Collections
UavsarData = collections.namedtuple("UavsarData",["dtarray","lon","lat","TS"]);
LevData = collections.namedtuple("LevData",["name","lat","lon","dtarray", "leveling"]);
TREData = collections.namedtuple("TREData",["lon","lat","vvel","vvel_std",
	"evel","evel_std","starttime","endtime"]);

# GET FILE NAMES
def get_file_dictionary():
	this_dict = {};
	ifile=open("directory_points.txt");
	for line in ifile:
		data_type = line.split(":")[0];
		total_data_files = line.split(":")[1][1:-1];
		this_dict[data_type]=total_data_files;
	ifile.close();
	return this_dict;


# UAVSAR INPUT FUNCTIONS
def inputs_uavsar(filename):
	# Reads a TS file from my interferogram processing
	print("Reading UAVSAR files in directory %s" % filename);
	[tdata, xdata, ydata, zdata] = netcdf_read_write.read_3D_netcdf(filename+"TS.nc");
	print("tdata:",tdata);
	print("zdata:",np.shape(zdata));
	lon = netcdf_read_write.read_grd(filename+"cut_lon.nc");
	print("lon:",np.shape(lon));
	lat = netcdf_read_write.read_grd(filename+"cut_lat.nc");
	print("lat:",np.shape(lat));
	dtarray = [];
	day0 = dt.datetime.strptime("2009-04-24","%Y-%m-%d"); # Hard-coded
	for i in range(len(tdata)):
		dtarray.append(day0+dt.timedelta(days=int(tdata[i])));
	myUAVSAR = UavsarData(dtarray=dtarray, lon=lon, lat=lat, TS=zdata);
	return myUAVSAR;


# LEVELING INPUT FUNCTIONS
def inputs_leveling(data_filename, errors_filename):
	# New data structure: 
	# names, lat, lon, datetimes, lists of deformation values

	print("Reading in %s" % data_filename);
	wb = xlrd.open_workbook(data_filename);
	sheet = wb.sheet_by_index(0);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];
	[rownum, colnum, old_values, new_values] = read_errors(errors_filename);
	data = implement_changes(data, rownum, colnum, old_values, new_values);
	
	dtarray = get_datetimes(data[0][1:-1]);	
	names = [data[i][0] for i in range(1,numrows)];
	
	# # Get the latitude and longitude information
	latsheet = wb.sheet_by_index(2);
	latitudes = latsheet.col_values(1)[1:];
	longitudes = latsheet.col_values(2)[1:];
	ll_names = latsheet.col_values(0)[1:];	
	names, lons, lats = match_lon_lat(names, latitudes, longitudes, ll_names);

	leveling_array = [];
	for i in range(1,numrows):
		single_leveling_array=data[i][1:13];
		single_leveling_array=clean_single_ts(single_leveling_array);
		leveling_array.append(single_leveling_array);

	myLev = LevData(name=names, lon=lons, lat=lats, dtarray=dtarray, leveling=leveling_array);
	return myLev;

def read_errors(filename):
	print("Reading documented errors in %s " % (filename) )
	rownum=[]; colnum=[]; old_values=[]; new_values=[];
	ifile=open(filename, 'r');
	for  line in ifile:
		temp=line.split("::");
		if temp[0]=="Row":
			continue;
		rownum.append(int(temp[0]));
		colnum.append(int(temp[1]));
		old_values.append(temp[2]);
		new_values.append(temp[3]);
	ifile.close();
	return [rownum, colnum, old_values, new_values];

def implement_changes(data, rownum, colnum, old_values, new_values):
	print("Implementing changes to data:");
	for i in range(len(rownum)):
		# Catching any bugs and errors
		if str(data[rownum[i]][colnum[i]]) != str(old_values[i]):
			print("PROBLEM at row %d! Attempting to modify data. Value found %s does not match expected value %s " 
				% (i, data[rownum[i]][colnum[i]], old_values[i]) );
			print("Skipping.");
			continue;

		else:
			# Over-writing the data
			print("MODIFYING data at %d, %d: %s goes to %s" % (rownum[i], colnum[i], old_values[i], new_values[i]) );
			if type(data[rownum[i]][colnum[i]]) == str:
				data[rownum[i]][colnum[i]] = new_values[i];
			elif type(data[rownum[i]][colnum[i]]) == float:
				data[rownum[i]][colnum[i]] = float(new_values[i]);
			elif type(data[rownum[i]][colnum[i]]) == int:
				data[rownum[i]][colnum[i]] = int(new_values[i]);
	return data; 

def get_datetimes(timestrings):
	dtarray=[];
	for i in range(len(timestrings)):
		# Normal dates
		if " 88 " in timestrings[i]:
			temp = timestrings[i].split(" 88 ");
			temp2 = temp[1].split();
			mmm = temp2[0];
			year = temp2[1];
			dtarray.append(dt.datetime.strptime(year+"-"+mmm+"-01","%Y-%b-%d"));  # issue here, but not too bad. Fix at work. 
		else: # For the special cases
			if "NOLTE 2008" in timestrings[i]:
				dtarray.append(dt.datetime.strptime("2008-Nov-01","%Y-%b-%d"));
	return dtarray;

def match_lon_lat(names, lats, lons, ll_names):
	# Pair up the latlon info with the timeseries info
	matched_lons=[]; matched_lats=[];
	for i in range(len(names)):
		find_name=names[i];
		if names[i]=="Y-1225 Datum":
			find_name = "Y 1225";
		idx = ll_names.index(find_name);
		matched_lats.append(lats[idx]);
		matched_lons.append(lons[idx]);
	return [names, matched_lons, matched_lats];

def clean_single_ts(array):
	newarray=[];
	for i in range(len(array)):
		if str(array[i])=="-" or str(array[i])=="DESTROYED" or str(array[i])=="DAMAGED" or str(array[i])=="NOT" or str(array[i])=="FOUND":
			newarray.append(np.nan);
		else:
			newarray.append(array[i]);
	return newarray;

# LEVELING COMPUTE FUNCITON (REFERENCE TO DATUM)
def compute_rel_to_datum_nov_2009(data):
	# Skips the 2008 measurement. Returns an object that is 83x10
	arrays_of_ref_leveling=[];
	for i in range(len(data.name)):
		
		# Automatically find the first day that matters.  Either after 2008 or has data. 
		for j in range(len(data.dtarray)):
			if ~np.isnan(data.leveling[i][j]) and data.dtarray[j]>dt.datetime.strptime("2009-01-01","%Y-%m-%d"):
				idx = j;  # this is the first date after 2009 that has data
				break;

		# Accounting for a change in Datum height in 2014
		idx_early = 6; # the placement of 2014 before adjustment on the spreadsheet
		idx_late = 7;  # the placement of 2014 after adjustment on the spreadsheet
		step = data.leveling[i][idx_early] - data.leveling[i][idx_late];

		referenced_dates = []; referenced_data = [];

		for j in range(1,len(data.dtarray)):  #skipping 2008 anyway. 
			if j==6:
				continue;  # passing over the 2014 measurement before re-referencing. 
			if data.dtarray[j]>dt.datetime.strptime("2014-01-01","%Y-%m-%d"):
				referenced_dates.append(data.dtarray[j]);
				referenced_data.append(data.leveling[i][j] - data.leveling[i][idx] + step);
			else:
				referenced_dates.append(data.dtarray[j]);
				referenced_data.append(data.leveling[i][j] - data.leveling[i][idx]);

		arrays_of_ref_leveling.append(referenced_data);
	referenced_object = LevData(name=data.name, lon=data.lon, lat = data.lat, dtarray=referenced_dates, leveling=arrays_of_ref_leveling);
	return referenced_object; 


# TSX INPUT FUNCTIONS
def inputs_tsx(filename):
	# Reading data from TSX TRE data. 
	print("Reading in %s" % filename);
	wb = xlrd.open_workbook(filename);

	# Vertical velocities
	sheet = wb.sheet_by_index(2);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];
	length, width = np.shape(data);
	zvel = [data[i][0] for i in range(1,length)];
	zvel_std = [data[i][1] for i in range(1,length)];
	lat = [data[i][2] for i in range(1,length)];
	lon = [data[i][3] for i in range(1,length)];

	# East velocities
	sheet = wb.sheet_by_index(3);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];
	length, width = np.shape(data);
	evel = [data[i][0] for i in range(1,length)];
	evel_std = [data[i][1] for i in range(1,length)];

	starttime=dt.datetime.strptime("2012-09-01","%Y-%m-%d");
	endtime=dt.datetime.strptime("2013-09-01","%Y-%m-%d");  # Hard coded for this experiment. 

	# Packaging it up
	myTSX = TREData(lon=lon, lat=lat, vvel=zvel, vvel_std=zvel_std, evel=evel, evel_std=evel_std, 
		starttime=starttime, endtime=endtime);

	return myTSX;


# S1 INPUT FUNCTIONS
def inputs_S1(filename):
	# Reading data from SNT1 TRE data. 
	print("Reading in %s" % filename);
	wb = xlrd.open_workbook(filename);

	# Vertical velocities
	sheet = wb.sheet_by_index(2);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];
	length, width = np.shape(data);
	zvel = [data[i][0] for i in range(1,length)];
	zvel_std = [data[i][1] for i in range(1,length)];
	lat = [data[i][2] for i in range(1,length)];
	lon = [data[i][3] for i in range(1,length)];

	# East velocities
	sheet = wb.sheet_by_index(3);
	numcols = sheet.ncols;
	numrows = sheet.nrows; 
	data = [[sheet.cell_value(r,c) for c in range(numcols)] for r in range(numrows)];
	length, width = np.shape(data);
	evel = [data[i][0] for i in range(1,length)];
	evel_std = [data[i][1] for i in range(1,length)];

	if "SNT1" in filename:
		starttime=dt.datetime.strptime("2014-04-01","%Y-%m-%d");
		endtime=dt.datetime.strptime("2018-04-01","%Y-%m-%d");  # Hard coded for this experiment. 
	if "SNT2" in filename:
		starttime=dt.datetime.strptime("2018-05-01","%Y-%m-%d");
		endtime=dt.datetime.strptime("2019-08-01","%Y-%m-%d");  # Hard coded for this experiment. 

	# Packaging it up
	myS1 = TREData(lon=lon, lat=lat, vvel=zvel, vvel_std=zvel_std, evel=evel, evel_std=evel_std, 
		starttime=starttime, endtime=endtime);

	return myS1;





