"""
Research-specific read functions to preprocess GNSS time series.
NOTE: We use the Timeseries object from the GNSS repo for all uses of this object.
import GNSS_TimeSeries_Viewers.gps_tools.gps_io_functions as gps_io_functions
Timeseries = GNSS_TimeSeries_Viewers.gps_tools.gps_io_functions.Timeseries;
Reference: Timeseries = collections.namedtuple("Timeseries",
['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su', 'EQtimes']);  # in mm
"""


import GNSS_TimeSeries_Viewers.gps_tools as gpstools
import datetime as dt


def read_station_ts_NBGF(gps_bbox, gps_reference, remove_coseismic=0, network='pbo', blacklist=()):
    """
    Read a set of GNSS stations. Specific to North Brawley Geothermal Field.
    Good to have around as an example.
    The Methods section of the paper can basically be read straight from this function.
    """
    gps_data_config_file = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/config.txt"  # sys path
    starttime1 = dt.datetime.strptime("20100403", "%Y%m%d");
    endtime1 = dt.datetime.strptime("20100405", "%Y%m%d");
    starttime2 = dt.datetime.strptime("20200328", "%Y%m%d");
    endtime2 = dt.datetime.strptime("20200330", "%Y%m%d");  # parameters for removal of EMC postseismic transient
    station_names, _, _ = gpstools.stations_within_radius.get_stations_within_box(gps_data_config_file,
                                                                                  coord_box=gps_bbox, network=network);
    print(station_names);
    [dataobj_list, offsetobj_list, eqobj_list, _] = \
        gpstools.gps_input_pipeline.multi_station_inputs(station_names, blacklist, network, "NA",
                                                         gps_data_config_file);
    # Next 5 lines: Importing BRAW from UNR
    [myData, offset_obj, eq_obj] = gpstools.gps_input_pipeline.get_station_data("BRAW", "unr", gps_data_config_file,
                                                                                refframe="NA")
    myData = gpstools.gps_ts_functions.impose_time_limits(myData, dt.datetime.strptime("20080505", "%Y%m%d"),
                                                          dt.datetime.strptime("20200101", "%Y%m%d"));
    dataobj_list.append(myData);
    offsetobj_list.append(offset_obj);
    eqobj_list.append(eq_obj);

    # Now we are doing a bit of adjustments, for seasonal corrections and base station.
    cleaned_objects = [];
    for i in range(len(dataobj_list)):
        one_object = dataobj_list[i];
        newobj = gpstools.offsets.remove_offsets(one_object, offsetobj_list[i]);  # remove antenna offsets

        if newobj.name == 'BRAW':
            newobj = gpstools.gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=0, seasonals_type="lssq",
                                                                      data_config_file=gps_data_config_file,
                                                                      remove_trend=0);
        else:
            newobj = gpstools.gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1,
                                                                      seasonals_type="nldas",
                                                                      data_config_file=gps_data_config_file,
                                                                      remove_trend=0);

        if remove_coseismic:
            print("Removing coseismic offsets");
            newobj = gpstools.offsets.remove_offsets(newobj, eqobj_list[i]);

        newobj = gpstools.gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

        # Here we detrend using pre-2010 velocities,
        # assuming tectonic strain accumulation won't contribute to geothermal field deformation.
        endtime = dt.datetime.strptime("2010-04-01", "%Y-%m-%d");
        [east_slope, north_slope, vert_slope, _, _, _] = gpstools.gps_ts_functions.get_slope(newobj, endtime=endtime,
                                                                                             missing_fraction=0.2);
        east_params = [east_slope, 0, 0, 0, 0];
        north_params = [north_slope, 0, 0, 0, 0];
        vert_params = [vert_slope, 0, 0, 0, 0];
        # Remove postseismic transient by Hines model
        model_obj = gpstools.gps_postseismic_remove.get_station_hines(newobj.name, gps_data_config_file);
        newobj = gpstools.gps_postseismic_remove.remove_by_model(newobj, model_obj, starttime1, endtime1, starttime2,
                                                                 endtime2);
        # Remove coseismic offset if within window
        newobj = gpstools.gps_ts_functions.detrend_data_by_value(newobj, east_params, north_params, vert_params);
        cleaned_objects.append(newobj)

    # Subtracting the reference GPS station, if desired.
    if gps_reference == "NO_REF":
        ref_dataobjlist = cleaned_objects;
    else:
        ref_dataobjlist = [];
        reference_station = [x for x in cleaned_objects if x.name == gps_reference][0];
        for one_object in cleaned_objects:
            refobj = gpstools.gps_ts_functions.get_referenced_data(one_object, reference_station);
            ref_dataobjlist.append(refobj);

    return ref_dataobjlist;
