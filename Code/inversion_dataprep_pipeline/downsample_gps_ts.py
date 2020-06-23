# June 2020
# A series of functions to help get GPS TS into 

import numpy as np 
import datetime as dt 
import matplotlib.pyplot as plt 
import collections
import stations_within_radius
import gps_input_pipeline
import gps_ts_functions

Timeseries = collections.namedtuple("Timeseries",['name','coords','dtarray','dN', 'dE','dU','Sn','Se','Su','EQtimes']);  # in mm

def read_station_ts(gps_bbox, gps_reference):
	blacklist=[];
	station_names,_,_ = stations_within_radius.get_stations_within_box(gps_bbox);
	[dataobj_list, offsetobj_list, eqobj_list, _] = gps_input_pipeline.multi_station_inputs(station_names, blacklist, "pbo","NA");
	
	# Now we are doing a bit of adjustments, like for coseismic offsets and base stations, etc.
	ref_dataobjlist = [];
	for one_object in dataobj_list:
		if one_object.name==gps_reference:
			reference_station = one_object;

	for one_object in dataobj_list:
		refobj = gps_ts_functions.get_referenced_data(one_object, reference_station);
		ref_dataobjlist.append(refobj);

	return ref_dataobjlist;
	# return dataobj_list;

def subsample_in_time(station, starttime, endtime):
	# Take a station and give us the data points corresponding to the starttime and endtime
	# return E0, N0, U0, E1, N1, U1;
	dE_array_start=[]; dN_array_start=[]; dU_array_start=[];
	dE_array_end=[]; dN_array_end=[]; dU_array_end=[];
	for i in range(len(station.dtarray)):
		if abs((station.dtarray[i] - starttime).days) < 30:
			dE_array_start.append(station.dE[i]);
			dN_array_start.append(station.dN[i]);
			dU_array_start.append(station.dU[i]);
		if abs((station.dtarray[i] - endtime).days) < 30:
			dE_array_end.append(station.dE[i]);
			dN_array_end.append(station.dN[i]);
			dU_array_end.append(station.dU[i]);
	if len(dE_array_start)>2:
		E0 = np.nanmean(dE_array_start);
		N0 = np.nanmean(dN_array_start);
		U0 = np.nanmean(dU_array_start);
	else:
		E0 = np.nan; N0 = np.nan; U0 = np.nan;
	if len(dE_array_end)>2:
		E1 = np.nanmean(dE_array_end);
		N1 = np.nanmean(dN_array_end);
		U1 = np.nanmean(dU_array_end);
	else:
		E1 = np.nan; N1 = np.nan; U1 = np.nan;

	return E0, N0, U0, E1, N1, U1;

def get_displacements_show_ts(stations, starttime, endtime, gps_sigma, prep_dir):
	# Get the values of TS at starttime and endtime
	startlim = starttime - dt.timedelta(days=365);
	endlim = endtime + dt.timedelta(days=365);
	gps_displacements_object=[];

	for station in stations:
		E0, N0, U0, E1, N1, U1 = subsample_in_time(station, starttime, endtime);
		one_object = Timeseries(name=station.name, coords=station.coords, dtarray=[starttime, endtime], 
			dN=[0, N1-N0], dE=[0, E1-E0], dU=[0,U1-U0],
			Sn=[gps_sigma,gps_sigma], Se=[gps_sigma,gps_sigma], Su=[3*gps_sigma,3*gps_sigma],
			EQtimes=station.EQtimes);
		gps_displacements_object.append(one_object);

		f,axarr = plt.subplots(3,1,figsize=(12,8),dpi=300);
		axarr[0].plot(station.dtarray, station.dE,'.');
		axarr[0].set_xlim([startlim, endlim]);
		axarr[0].plot(starttime, E0, '.', color='red',markersize=15);
		axarr[0].plot(endtime, E1, '.', color='red',markersize=15);
		axarr[0].set_ylabel('East (mm)')
		axarr[1].plot(station.dtarray, station.dN,'.');
		axarr[1].set_xlim([startlim, endlim]);
		axarr[1].plot(starttime, N0, '.', color='red',markersize=15);
		axarr[1].plot(endtime, N1, '.', color='red',markersize=15);
		axarr[1].set_ylabel('North (mm)');
		axarr[2].plot(station.dtarray, station.dU,'.');
		axarr[2].set_xlim([startlim, endlim]);
		axarr[2].plot(starttime, U0, '.', color='red',markersize=15);
		axarr[2].plot(endtime, U1, '.', color='red',markersize=15);		
		axarr[2].set_ylabel('Up (mm)');
		plt.savefig(prep_dir+"gps_"+station.name+"_ts.png");
		plt.close();

	return gps_displacements_object;
