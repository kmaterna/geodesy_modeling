import GNSS_TimeSeries_Viewers.gps_tools as gps_tools
import datetime as dt
from . import gnss_object


def read_station_ts(gps_bbox, gps_reference, remove_coseismic=0):
    """A reading function specific to Brawley right now. """
    blacklist = [];
    network = 'pbo'
    station_names, _, _ = gps_tools.stations_within_radius.get_stations_within_box(gnss_object.gps_data_config_file,
                                                                         coord_box=gps_bbox, network=network);
    print(station_names);
    [dataobj_list, offsetobj_list, eqobj_list, _] = \
        gps_tools.gps_input_pipeline.multi_station_inputs(station_names, blacklist, network, "NA",
                                                gnss_object.gps_data_config_file);
    # Importing BRAW
    [myData, offset_obj, eq_obj] = gps_tools.gps_input_pipeline.get_station_data("BRAW", "unr", gnss_object.gps_data_config_file,
                                                                       refframe="NA")
    myData = gps_tools.gps_ts_functions.impose_time_limits(myData, dt.datetime.strptime("20080505", "%Y%m%d"),
                                                 dt.datetime.strptime("20200101", "%Y%m%d"));
    dataobj_list.append(myData);
    offsetobj_list.append(offset_obj);
    eqobj_list.append(eq_obj);

    # Now we are doing a bit of adjustments, for seasonal corrections and base station.
    cleaned_objects = [];
    for i in range(len(dataobj_list)):
        one_object = dataobj_list[i];
        newobj = gps_tools.offsets.remove_offsets(one_object, offsetobj_list[i]);  # will remove antenna offsets from everything

        if newobj.name == 'BRAW':
            newobj = gps_tools.gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1, seasonals_type="lssq",
                                                             data_config_file=gnss_object.gps_data_config_file,
                                                             remove_trend=0);
        else:
            newobj = gps_tools.gps_seasonal_removals.make_detrended_ts(newobj, seasonals_remove=1, seasonals_type="nldas",
                                                             data_config_file=gnss_object.gps_data_config_file,
                                                             remove_trend=0);

        if remove_coseismic:
            print("Removing coseismic offsets");
            newobj = gps_tools.offsets.remove_offsets(newobj, eqobj_list[i]);

        newobj = gps_tools.gps_ts_functions.remove_outliers(newobj, 20);  # 20mm outlier definition

        # Here we detrend using pre-2010 velocities,
        # assuming tectonic strain accumulation won't contribute to geothermal field deformation.
        # Remove the postseismic by the Hines model
        endtime = dt.datetime.strptime("2010-04-01", "%Y-%m-%d");
        [east_slope, north_slope, vert_slope, _, _, _] = gps_tools.gps_ts_functions.get_slope(newobj, endtime=endtime,
                                                                                    missing_fraction=0.2);
        east_params = [east_slope, 0, 0, 0, 0];
        north_params = [north_slope, 0, 0, 0, 0];
        vert_params = [vert_slope, 0, 0, 0, 0];
        newobj = gps_tools.gps_postseismic_remove.remove_by_model(newobj, gnss_object.gps_data_config_file);
        # This will actually remove the coseismic offset if within the window.
        newobj = gps_tools.gps_ts_functions.detrend_data_by_value(newobj, east_params, north_params, vert_params);
        cleaned_objects.append(newobj)

    # Subtracting the reference GPS station.
    ref_dataobjlist = [];
    reference_station = [x for x in cleaned_objects if x.name == gps_reference][0];
    for one_object in cleaned_objects:
        refobj = gps_tools.gps_ts_functions.get_referenced_data(one_object, reference_station);
        ref_dataobjlist.append(refobj);

    return ref_dataobjlist;
