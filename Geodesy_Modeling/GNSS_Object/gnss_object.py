import collections

# Global variables
Timeseries = collections.namedtuple("Timeseries", ['name', 'coords', 'dtarray', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su',
                                                   'EQtimes']);  # in mm
gps_data_dir = "/Users/kmaterna/Documents/B_Research/GEOPHYS_DATA/GPS_POS_DATA/"
gps_data_config_file = gps_data_dir + "config.txt"
