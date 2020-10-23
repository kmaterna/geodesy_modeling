#!/usr/bin/env python
# September 2020
# Examine a few apparent offsets that could be creep events related to the Imperial Fault

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt
import collections
import gps_input_pipeline
import gps_seasonal_removals
import offsets
import gps_ts_functions
import gps_postseismic_remove

Experiment = collections.namedtuple('Experiment',['station_codes','starttime','endtime','labeltime','outdir']);
data_config_file = "../../../../Mendocino_Geodesy/GPS_POS_DATA/config.txt"

# Experiment parameters for controlling the plots
imperial_2015_station_codes = {'P744':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_nldas','unr_NA_lssq','nmt_NA_lssq'), 
                               'P498':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_nldas','unr_NA_lssq','nmt_NA_lssq')};
shf_2017_station_codes      = {'P503':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_lssq','unr_NA_lssq','cwu_NA_nldas')};
imperial_2006_station_codes = {'IVCO':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_lssq','unr_NA_lssq','cwu_NA_nldas')};
imperial_2004_station_codes = {'IVCO':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_lssq','unr_NA_lssq','cwu_NA_nldas')};
imperial_2010_station_codes = {'P498':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_nldas','unr_NA_lssq','nmt_NA_lssq'),
                               'P744':('cwu_NA_lssq','nmt_NA_lssq','pbo_NA_lssq','unr_NA_lssq','cwu_NA_nldas')};

IF_2015_Expt = Experiment(station_codes=imperial_2015_station_codes, starttime='20170303', endtime='20180401', 
               labeltime='20160606', outdir='IF_2015/');
SHF_2017_Expt = Experiment(station_codes=shf_2017_station_codes, starttime='20090508', endtime='20121101', 
                labeltime='20060606', outdir='SHF_2017/');
IF_2006_Expt = Experiment(station_codes=imperial_2006_station_codes, starttime='20060601', endtime='20060606', 
                labeltime='20060606', outdir='IF_2006/');
IF_2004_Expt = Experiment(station_codes=imperial_2004_station_codes, starttime='20040101', endtime='20050101', 
                labeltime='20050606', outdir='IF_2004/');
IF_2010_Expt = Experiment(station_codes=imperial_2010_station_codes, starttime='20090101', endtime='20201001',
                labeltime='20100606', outdir='IF_2010/');

def driver():
    myExp = IF_2010_Expt;
    Data_objlist, offset_objlist, eq_objlist, seasonals_list = inputs(myExp);
    cleaned_objects = compute(Data_objlist, offset_objlist, eq_objlist, seasonals_list);
    outputs(cleaned_objects, offset_objlist, eq_objlist, myExp);
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


def compute(Data_objlist, offset_objlist, eq_objlist, seasonals_list):
    # Preferred method: remove offsets/earthquakes, remove PS, remove_seasonals, get slope, detrend.

    cleaned_objects = []; 
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

    return cleaned_objects;


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

if __name__=="__main__":
    driver();
