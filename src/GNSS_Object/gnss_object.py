import GNSS_TimeSeries_Viewers.gps_tools.gps_io_functions as gps_io_functions

# Global variables
Timeseries = gps_io_functions.Timeseries;
# For reference: Timeseries = collections.namedtuple("Timeseries",
# ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su', 'EQtimes']);  # in mm
gps_data_dir = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/"
gps_data_config_file = gps_data_dir + "config.txt"
