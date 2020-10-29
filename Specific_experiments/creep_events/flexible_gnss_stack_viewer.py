#!/usr/bin/env python
# September 2020
# Examine a few apparent offsets that could be creep events related to the Imperial Fault

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import gps_input_pipeline
import gps_seasonal_removals
import offsets
import gps_ts_functions
import gps_postseismic_remove


data_config_file = "../../../../Mendocino_Geodesy/GPS_POS_DATA/config.txt"

def driver(myExp):
    Data_objlist, offset_objlist, eq_objlist, seasonals_list = inputs(myExp);
    cleaned_objects, start_end_values = compute(Data_objlist, offset_objlist, eq_objlist, seasonals_list, myExp);
    # outputs(cleaned_objects, offset_objlist, eq_objlist, myExp);
    outputs_singlepanel(cleaned_objects, offset_objlist, eq_objlist, start_end_values, myExp);
    return;

def unpack_stations_and_options(myExp):
    stations, networks, refframes, seasonal_types = [], [], [], [];
    for station in myExp.station_codes.keys():
        for control_string in myExp.station_codes[station]:
            [network, refframe, seasonal_type] = control_string.split('_');
            stations.append(station);
            networks.append(network);
            refframes.append(refframe);
            seasonal_types.append(seasonal_type);
    return [stations, networks, refframes, seasonal_types];

def inputs(datasources):
    Data_objlist, offset_objlist, eq_objlist, seasonals_list = [], [], [], [];
    [stations, networks, refframes, seasonal_types] = unpack_stations_and_options(datasources);
    for i in range(len(stations)):
        [Data, offset_obj, eq_obj] = gps_input_pipeline.get_station_data(stations[i], networks[i], 
            data_config_file, refframe=refframes[i]);
        Data_objlist.append(Data);
        offset_objlist.append(offset_obj);
        eq_objlist.append(eq_obj);
    return Data_objlist, offset_objlist, eq_objlist, seasonal_types;

def compute(Data_objlist, offset_objlist, eq_objlist, seasonals_list, ExpParams):
    # Preferred method: remove offsets/earthquakes, remove PS, remove_seasonals, get slope, detrend.

    cleaned_objects = []; 
    start_end_values = [];
    for i in range(len(Data_objlist)):
        newobj = Data_objlist[i];
        newobj = offsets.remove_offsets(newobj, offset_objlist[i]);
        newobj = offsets.remove_offsets(newobj, eq_objlist[i]);	
        newobj = gps_ts_functions.remove_outliers(newobj, 20);  
        newobj = gps_postseismic_remove.remove_by_model(newobj, data_config_file);
        newobj = gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1, seasonals_type=seasonals_list[i],
                                                         data_config_file=data_config_file, remove_trend=0); 

        endtime = dt.datetime.strptime("2010-04-01", "%Y-%m-%d");
        [east_slope, north_slope, vert_slope, _, _, _] = gps_ts_functions.get_slope(newobj, endtime=endtime,
                                                                                    missing_fraction=0.2);
        east_params = [east_slope, 0, 0, 0, 0];
        north_params = [north_slope, 0, 0, 0, 0];
        vert_params = [vert_slope, 0, 0, 0, 0];
        newobj = gps_ts_functions.detrend_data_by_value(newobj, east_params, north_params, vert_params);
        cleaned_objects.append(newobj);

        # Get the average values
        starttime = dt.datetime.strptime(ExpParams.close_starttime[0],"%Y%m%d");
        endtime = dt.datetime.strptime(ExpParams.close_starttime[1],"%Y%m%d");
        east_val_start, north_val_start, up_val_start = gps_ts_functions.get_means(newobj, starttime, endtime);
        starttime = dt.datetime.strptime(ExpParams.close_endtime[0],"%Y%m%d");
        endtime = dt.datetime.strptime(ExpParams.close_endtime[1],"%Y%m%d");
        east_val_end, north_val_end, up_val_end = gps_ts_functions.get_means(newobj, starttime, endtime);
        means_package = [east_val_start, north_val_start, up_val_start, east_val_end, north_val_end, up_val_end];
        start_end_values.append(means_package);

    return cleaned_objects, start_end_values;

def get_colors(station):
    if station == 'P498':
        color = 'dodgerblue'
    elif station == 'P744':
        color = 'darkcyan'
    elif station == 'IVCO':
        color = 'purple'  
    elif station == 'P503':
        color = 'indianred'
    else:
        color = 'black';  
    return color;


def outputs(Data_objlist2, offset_objlist, eq_objlist, myExp):

    [stations, networks, refframes, seasonal_types] = unpack_stations_and_options(myExp);
    plt.figure(figsize=(12,12),dpi=300);
    for i in range(len(Data_objlist2)):
        offset=i*12; 
        color = get_colors(Data_objlist2[i].name); 
        plt.plot(Data_objlist2[i].dtarray, offset+Data_objlist2[i].dE, linewidth=0, marker='.', color=color);
        text_string = stations[i]+'/'+networks[i]+'/'+seasonal_types[i];
        plt.text(dt.datetime.strptime(myExp.labeltime,'%Y%m%d'), offset+6, text_string, color='black',fontsize=12);
    plt.gca().set_xlim([dt.datetime.strptime(myExp.starttime,'%Y%m%d'), dt.datetime.strptime(myExp.endtime,'%Y%m%d')]);
    plt.gca().set_title('East Residuals');
    plt.gca().set_ylabel('Position (mm)'); 
    plt.gca().grid(True);
    plt.savefig("East.png");
    plt.close();

    plt.figure(figsize=(12,12),dpi=300);
    for i in range(len(Data_objlist2)):
        offset=i*12; 
        color = get_colors(Data_objlist2[i].name);
        plt.plot(Data_objlist2[i].dtarray, offset+Data_objlist2[i].dN, linewidth=0, marker='.', color=color);
        text_string = stations[i]+'/'+networks[i]+'/'+seasonal_types[i];
        plt.text(dt.datetime.strptime(myExp.labeltime,'%Y%m%d'), offset+6, text_string, color='black',fontsize=12);
    plt.gca().set_xlim([dt.datetime.strptime(myExp.starttime,'%Y%m%d'), dt.datetime.strptime(myExp.endtime,'%Y%m%d')]);
    plt.gca().set_title('North Residuals');
    plt.gca().set_ylabel('Position (mm)'); 
    plt.gca().grid(True);
    plt.savefig("North.png");
    plt.close();

    # plt.figure(figsize=(12,12),dpi=300);
    # for i in range(len(cleaned_objects)):
    #     offset=i*35; 
    #     color = get_colors(datasources[i][0]);
    #     plt.plot(cleaned_objects[i].dtarray, offset+cleaned_objects[i].dU, linewidth=0, marker='.', color=color);
    #     text_string = datasources[i][0]+'/'+datasources[i][1]+'/'+datasources[i][3]
    #     plt.text(labeltime, offset+17, text_string, color='black',fontsize=12);
    # plt.gca().set_xlim([starttime, endtime])
    # plt.gca().set_title('Vertical Residuals');
    # plt.gca().set_ylabel('Position (mm)'); 
    # plt.gca().grid(True);
    # plt.savefig("Up.png");
    # plt.close();
    return;


def outputs_singlepanel(Data_objlist, offset_objlist, eq_objlist, start_end_values, myExp):
    fig, axarr = plt.subplots(2,3,figsize=(15, 10), dpi=300);
    [stations, networks, refframes, seasonal_types] = unpack_stations_and_options(myExp);
    for i in range(len(Data_objlist)):
        offset=i*12; 
        color = get_colors(Data_objlist[i].name); 
        axarr[0][0].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dE, linewidth=0, marker='.', color=color);
        axarr[1][0].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dE, linewidth=0, marker='.', color=color);
        text_string = stations[i]+'/'+networks[i]+'/'+seasonal_types[i];
        axarr[0][0].text(dt.datetime.strptime(myExp.labeltime,'%Y%m%d'), offset+3, text_string, color='black',fontsize=12);
        if i==0:
            axarr[1][0].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_starttime[-1],'%Y%m%d')],[offset+start_end_values[0][0], offset+start_end_values[0][0]],'-k');
            axarr[1][0].plot([dt.datetime.strptime(myExp.close_endtime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[offset+start_end_values[0][3], offset+start_end_values[0][3]],'-k');

    avg_east_offset = [x[3]-x[0] for x in start_end_values];
    avg_east_offset = np.nanmean(avg_east_offset);

    [top, bottom] = axarr[0][0].get_ylim();
    axarr[0][0].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][0].plot([dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d'), dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][0].set_xlim([dt.datetime.strptime(myExp.starttime,'%Y%m%d'), dt.datetime.strptime(myExp.endtime,'%Y%m%d')]);
    axarr[0][0].set_title('East Residuals: Avg %.3fmm' % (avg_east_offset) );
    axarr[0][0].set_ylabel('Position (mm)');
    axarr[0][0].tick_params(axis='x', rotation=40);
    axarr[0][0].grid(True);

    axarr[1][0].set_xlim([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
        dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')]);
    axarr[1][0].set_ylabel('Position (mm)');
    axarr[1][0].grid(True);
    axarr[1][0].tick_params(axis='x', rotation=40);



    for i in range(len(Data_objlist)):
        offset=i*12; 
        color = get_colors(Data_objlist[i].name); 
        axarr[0][1].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dN, linewidth=0, marker='.', color=color);
        axarr[1][1].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dN, linewidth=0, marker='.', color=color);
        if i==0:
            axarr[1][1].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_starttime[-1],'%Y%m%d')],[offset+start_end_values[0][1], offset+start_end_values[0][1]],'-k');
            axarr[1][1].plot([dt.datetime.strptime(myExp.close_endtime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[offset+start_end_values[0][4], offset+start_end_values[0][4]],'-k');        

    avg_north_offset = [x[4]-x[1] for x in start_end_values];
    avg_north_offset = np.nanmean(avg_north_offset);

    [top, bottom] = axarr[0][1].get_ylim();
    axarr[0][1].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][1].plot([dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d'), dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][1].set_xlim([dt.datetime.strptime(myExp.starttime,'%Y%m%d'), dt.datetime.strptime(myExp.endtime,'%Y%m%d')]);
    axarr[0][1].set_title('North Residuals: Avg %.3fmm' % (avg_north_offset) );
    axarr[0][1].set_ylabel('Position (mm)');
    axarr[0][1].tick_params(axis='x', rotation=40);
    axarr[0][1].grid(True);

    axarr[1][1].set_xlim([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
        dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')]);
    axarr[1][1].set_ylabel('Position (mm)');
    axarr[1][1].grid(True);
    axarr[1][1].tick_params(axis='x', rotation=40);


    for i in range(len(Data_objlist)):
        offset=i*30; 
        color = get_colors(Data_objlist[i].name); 
        axarr[0][2].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dU, linewidth=0, marker='.', color=color);
        axarr[1][2].plot(Data_objlist[i].dtarray, offset+Data_objlist[i].dU, linewidth=0, marker='.', color=color);
        if i==0:
            axarr[1][2].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_starttime[-1],'%Y%m%d')],[offset+start_end_values[0][2], offset+start_end_values[0][2]],'-k');
            axarr[1][2].plot([dt.datetime.strptime(myExp.close_endtime[0],'%Y%m%d'), 
                dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[offset+start_end_values[0][5], offset+start_end_values[0][5]],'-k');        

    avg_vert_offset = [x[5]-x[2] for x in start_end_values];
    avg_vert_offset = np.nanmean(avg_vert_offset);

    [top, bottom] = axarr[0][2].get_ylim();
    axarr[0][2].plot([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][2].plot([dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d'), dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')],[top, bottom],'--k');
    axarr[0][2].set_xlim([dt.datetime.strptime(myExp.starttime,'%Y%m%d'), dt.datetime.strptime(myExp.endtime,'%Y%m%d')]);
    axarr[0][2].set_title('Vertical Residuals: Avg %.3fmm' % (avg_vert_offset) );
    axarr[0][2].set_ylabel('Position (mm)');
    axarr[0][2].tick_params(axis='x', rotation=40);
    axarr[0][2].grid(True);

    axarr[1][2].set_xlim([dt.datetime.strptime(myExp.close_starttime[0],'%Y%m%d'), 
        dt.datetime.strptime(myExp.close_endtime[-1],'%Y%m%d')]);
    axarr[1][2].set_ylabel('Position (mm)');
    axarr[1][2].grid(True);
    axarr[1][2].tick_params(axis='x', rotation=40);

    fig.savefig(myExp.outdir+'singlepanel.png');

    return;

if __name__=="__main__":
    driver();

