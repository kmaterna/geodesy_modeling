# June 2020
# A series of functions to help get GPS TS into 

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import collections
import stations_within_radius
import gps_input_pipeline
import gps_ts_functions
import gps_seasonal_removals
import gps_postseismic_remove
import offsets

# Global variables
Timeseries = collections.namedtuple("Timeseries", ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su',
                                                   'EQtimes']);  # in mm
gps_data_dir = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/"
gps_data_config_file = gps_data_dir + "config.txt"


def read_station_ts(gps_bbox, gps_reference, remove_coseismic=0):
    blacklist = [];
    network = 'pbo'
    station_names, _, _ = stations_within_radius.get_stations_within_box(gps_data_config_file, coord_box=gps_bbox,
                                                                         network=network);
    print(station_names);
    [dataobj_list, offsetobj_list, eqobj_list, _] = gps_input_pipeline.multi_station_inputs(station_names, blacklist,
                                                                                            network, "NA",
                                                                                            gps_data_config_file);
    # Importing BRAW
    [myData, offset_obj, eq_obj] = gps_input_pipeline.get_station_data("BRAW", "unr", gps_data_config_file,
                                                                       refframe="NA")
    myData = gps_ts_functions.impose_time_limits(myData, dt.datetime.strptime("20080505", "%Y%m%d"),
                                                 dt.datetime.strptime("20200101", "%Y%m%d"));
    dataobj_list.append(myData);
    offsetobj_list.append(offset_obj);
    eqobj_list.append(eq_obj);

    # Now we are doing a bit of adjustments, for seasonal corrections and base station.
    cleaned_objects = [];
    for i in range(len(dataobj_list)):
        one_object = dataobj_list[i];
        newobj = offsets.remove_offsets(one_object, offsetobj_list[i]);  # will remove antenna offsets from everything

        if newobj.name == 'BRAW':
            newobj = gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1, seasonals_type="lssq",
                                                             data_config_file=gps_data_config_file, remove_trend=0);
        else:
            newobj = gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1, seasonals_type="nldas",
                                                             data_config_file=gps_data_config_file, remove_trend=0);

        if remove_coseismic:
            print("Removing coseismic offsets");
            newobj = offsets.remove_offsets(newobj, eqobj_list[i]);

        newobj = gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

        # Here we detrend using pre-2010 velocities,
        # assuming tectonic strain accumulation won't contribute to geothermal field deformation.
        # Remove the postseismic by the Hines model
        endtime = dt.datetime.strptime("2010-04-01", "%Y-%m-%d");
        [east_slope, north_slope, vert_slope, _, _, _] = gps_ts_functions.get_slope(newobj, endtime=endtime,
                                                                                    missing_fraction=0.2);
        east_params = [east_slope, 0, 0, 0, 0];
        north_params = [north_slope, 0, 0, 0, 0];
        vert_params = [vert_slope, 0, 0, 0, 0];
        newobj = gps_postseismic_remove.remove_by_model(newobj, gps_data_config_file);
        # This will actually remove the coseismic offset if within the window.
        newobj = gps_ts_functions.detrend_data_by_value(newobj, east_params, north_params, vert_params);
        cleaned_objects.append(newobj)

    # Subtracting the reference GPS station.
    ref_dataobjlist = [];
    reference_station = [x for x in cleaned_objects if x.name == gps_reference][0];
    for one_object in cleaned_objects:
        refobj = gps_ts_functions.get_referenced_data(one_object, reference_station);
        ref_dataobjlist.append(refobj);

    return ref_dataobjlist;


def get_displacements_show_ts(stations, starttime, endtime, gps_sigma, prep_dir):
    """Get the values of TS at starttime and endtime"""
    startlim = starttime - dt.timedelta(days=1065);
    endlim = endtime + dt.timedelta(days=1065);
    gps_displacements_object = [];

    for station in stations:
        E0, N0, U0, E1, N1, U1 = subsample_in_time(station, starttime, endtime);
        one_object = Timeseries(name=station.name, coords=station.coords, dtarray=[starttime, endtime],
                                dN=[0, N1 - N0], dE=[0, E1 - E0], dU=[0, U1 - U0],
                                Sn=[gps_sigma, gps_sigma], Se=[gps_sigma, gps_sigma], Su=[3 * gps_sigma, 3 * gps_sigma],
                                EQtimes=station.EQtimes);
        gps_displacements_object.append(one_object);

        f, axarr = plt.subplots(3, 1, figsize=(12, 8), dpi=300);
        axarr[0].plot(station.dtarray, station.dE, '.');
        axarr[0].set_xlim([startlim, endlim]);
        axarr[0].plot(starttime, E0, '.', color='red', markersize=15);
        axarr[0].plot(endtime, E1, '.', color='red', markersize=15);
        axarr[0].plot([starttime, endtime], [E0, E1], color='red');
        axarr[0].set_ylabel('East (mm)')
        axarr[1].plot(station.dtarray, station.dN, '.');
        axarr[1].set_xlim([startlim, endlim]);
        axarr[1].plot(starttime, N0, '.', color='red', markersize=15);
        axarr[1].plot(endtime, N1, '.', color='red', markersize=15);
        axarr[1].set_ylabel('North (mm)');
        axarr[2].plot(station.dtarray, station.dU, '.');
        axarr[2].set_xlim([startlim, endlim]);
        axarr[2].plot(starttime, U0, '.', color='red', markersize=15);
        axarr[2].plot(endtime, U1, '.', color='red', markersize=15);
        axarr[2].set_ylabel('Up (mm)');
        plt.savefig(prep_dir + "gps_" + station.name + "_ts.png");
        plt.close();

    return gps_displacements_object;


def add_gps_constant_offset(displacement_object, enu_constant_offset):
    """
    In case your reference station undergoes some offset, you might want to put that back
    to all data because it's been referenced-out already.
    Allows for more realistic GNSS inversions
    Offset and displacement obj in mm
    """
    new_gps_displacements_object = [];

    for one_object in displacement_object:
        dE = [i + enu_constant_offset[0] for i in one_object.dE];
        dN = [i + enu_constant_offset[1] for i in one_object.dN];
        dU = [i + enu_constant_offset[2] for i in one_object.dU];
        object_after_offset = Timeseries(name=one_object.name, coords=one_object.coords, dtarray=one_object.dtarray,
                                         dN=dN, dE=dE, dU=dU, Sn=one_object.Sn, Se=one_object.Se, Su=one_object.Su,
                                         EQtimes=one_object.EQtimes);
        new_gps_displacements_object.append(object_after_offset);
    return new_gps_displacements_object;


def subsample_in_time(station, starttime, endtime):
    """
    Take a station and give us the data points corresponding to the starttime and endtime
    return E0, N0, U0, E1, N1, U1
    """
    dE_array_start = [];
    dN_array_start = [];
    dU_array_start = [];
    dE_array_end = [];
    dN_array_end = [];
    dU_array_end = [];
    for i in range(len(station.dtarray)):
        if abs((station.dtarray[i] - starttime).days) < 30:
            dE_array_start.append(station.dE[i]);
            dN_array_start.append(station.dN[i]);
            dU_array_start.append(station.dU[i]);
        if abs((station.dtarray[i] - endtime).days) < 30:
            dE_array_end.append(station.dE[i]);
            dN_array_end.append(station.dN[i]);
            dU_array_end.append(station.dU[i]);
    if len(dE_array_start) > 2:
        E0 = np.nanmean(dE_array_start);
        N0 = np.nanmean(dN_array_start);
        U0 = np.nanmean(dU_array_start);
    else:
        E0 = np.nan;
        N0 = np.nan;
        U0 = np.nan;
    if len(dE_array_end) > 2:
        E1 = np.nanmean(dE_array_end);
        N1 = np.nanmean(dN_array_end);
        U1 = np.nanmean(dU_array_end);
    else:
        E1 = np.nan;
        N1 = np.nan;
        U1 = np.nan;

    return E0, N0, U0, E1, N1, U1;
