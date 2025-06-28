import datetime as dt
import numpy as np
from tectonic_utils.read_write import netcdf_read_write
from geodesy_modeling.datatypes.InSAR_1D_Object import class_model


# DATA STRUCTURE
class GrdTSData:
    def __init__(self, dtarray, lon, lat, TS):
        self.dtarray = dtarray  # one-dimensional array
        self.lon = lon  # one-dimensional array
        self.lat = lat  # one-dimensional array
        self.TS = TS   # three-dimensional array

    def get_onetime_displacements(self, start_idx, end_idx):
        """
        Convert a particular epoch of a GridTS object into an InSAR_1D_Object of dots.
        Converts a 2d raster into a 1d vector in the process.
        """
        los_raster = np.subtract(self.TS[end_idx][:, :], self.TS[start_idx][:, :])
        los_vector = np.reshape(los_raster, (np.size(los_raster),))
        lon_vector = np.reshape(self.lon, (np.size(los_raster),))
        lat_vector = np.reshape(self.lat, (np.size(los_raster),))
        nodata_object = np.nan * np.zeros(np.shape(los_vector))
        InSAR_obj = class_model.Insar1dObject(lon=lon_vector, lat=lat_vector, LOS=los_vector, LOS_unc=nodata_object,
                                              lkv_E=nodata_object, lkv_N=nodata_object, lkv_U=nodata_object,
                                              starttime=self.dtarray[start_idx], endtime=self.dtarray[end_idx])
        return InSAR_obj

    def avg_disp(self, num, row, col, width_pixels=10):
        """
        Average around a few pixels.

        :param num: slice index number of the time series.
        :param row: row number at center of the sampling region
        :param col: column number at center of the sampling region
        :param width_pixels: integer, default 10
        """
        return np.nanmean(self.TS[num, row - width_pixels:row + width_pixels, col - width_pixels:col + width_pixels])


# INPUT FUNCTIONS FOR NETCDF FORMAT
def inputs_TS_grd(filename, lonfile, latfile, day0=dt.datetime.strptime("2009-04-24", "%Y-%m-%d")):
    """
    Reads a TS file with associated lat/lon files
    The files generally are not orthorectified grids
    GRDnetcdf has tdata (days since day0), x, y, and zdata (3D cube)
    lon and lat files are 2D arrays with corresponding lon and lat for each point
    day0 is the day of the first acquisition in the time series (hard coded for a InSAR_Timeseries track default)
    """
    print("Reading TS Grid file  %s" % filename)
    [tdata, _, _, zdata] = netcdf_read_write.read_3D_netcdf(filename)
    print("tdata:", tdata)
    print("   where Day0 of this time series is %s " % dt.datetime.strftime(day0, "%Y-%m-%d"))
    [_, _, lon] = netcdf_read_write.read_any_grd(lonfile)
    [_, _, lat] = netcdf_read_write.read_any_grd(latfile)
    print("lon and lat:", np.shape(lon))
    print("zdata:", np.shape(zdata))
    zdata_correct_size = []
    if np.shape(zdata[0])[0] == np.shape(lon)[0] + 1 and np.shape(zdata[0])[1] == np.shape(lat)[1] + 1:
        for i in range(len(zdata)):
            zdata_correct_size.append(zdata[i][0:-1, 0:-1])  # cutting off one pixel on each end for pixel node problem
    else:
        zdata_correct_size = zdata
    dtarray = []
    for i in range(len(tdata)):
        dtarray.append(day0 + dt.timedelta(days=int(tdata[i])))

    assert (np.shape(lon) == np.shape(zdata_correct_size[0])), ValueError("Lon and Data size don't match")
    assert (np.shape(lat) == np.shape(zdata_correct_size[0])), ValueError("Lat and Data size don't match")
    assert (np.shape(zdata_correct_size)[0] == len(dtarray)), ValueError("dtarray and zdata size don't match")
    myGridTS = GrdTSData(dtarray=dtarray, lon=lon, lat=lat, TS=zdata_correct_size)
    return myGridTS
