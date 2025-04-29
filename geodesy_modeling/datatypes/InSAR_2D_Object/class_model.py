import numpy as np
from Tectonic_Utils.geodesy import insar_vector_functions as ivf
from geodesy_modeling import general_utils
from cubbie.math_tools import grid_tools


class Insar2dObject:
    """
    A generalized 2D Grid InSAR format where all data fields are 2D grids.
    Displacements in mm (if LOS is a displacement measurement instead of phase or other)
    """
    def __init__(self, lon, lat, LOS, LOS_unc, lkv_E, lkv_N, lkv_U, starttime=None, endtime=None):
        self.lon = lon  # 1d array
        self.lat = lat  # 1d array
        self.LOS = LOS  # displacements, in mm, 2d Grid
        self.LOS_unc = LOS_unc   # grid, in mm, 2d Grid
        self.lkv_E = lkv_E  # Ground to Satellite, 2d grid
        self.lkv_N = lkv_N  # Ground to Satellite, 2d grid
        self.lkv_U = lkv_U  # Ground to Satellite, 2d grid
        self.starttime = starttime  # just metadata, datetime object
        self.endtime = endtime  # just metadata, datetime object
        (leny, lenx) = np.shape(self.LOS)
        if len(self.lon) != lenx:
            raise ValueError("length of InSAR_Obj lon array doesn't match shape of data.")
        if len(self.lat) != leny:
            raise ValueError("length of InSAR_Obj lat array doesn't match shape of data.")
        if np.shape(self.lkv_E) != np.shape(self.LOS):
            raise ValueError("Shape of InSAR data doesn't match shape of lkv arrays:", np.shape(self.LOS),
                             " versus ", np.shape(self.lkv_E))

    def impose_InSAR_bounding_box(self, bbox=(-180, 180, -90, 90)):
        """Impose a bounding box on some InSAR data. """
        xc, yc, los_c = grid_tools.clip_array_by_bbox(self.lon, self.lat, self.LOS, bbox)
        _, _, unc_c = grid_tools.clip_array_by_bbox(self.lon, self.lat, self.LOS_unc, bbox)
        _, _, lkve_c = grid_tools.clip_array_by_bbox(self.lon, self.lat, self.lkv_E, bbox)
        _, _, lkvn_c = grid_tools.clip_array_by_bbox(self.lon, self.lat, self.lkv_N, bbox)
        _, _, lkvu_c = grid_tools.clip_array_by_bbox(self.lon, self.lat, self.lkv_U, bbox)
        smaller_object = Insar2dObject(lon=xc, lat=yc, LOS=los_c, LOS_unc=unc_c, lkv_E=lkve_c, lkv_N=lkvn_c,
                                       lkv_U=lkvu_c, starttime=self.starttime, endtime=self.endtime)
        return smaller_object

    def flip_los_sign(self):
        new_InSAR_obj = Insar2dObject(lon=self.lon, lat=self.lat, LOS=np.multiply(self.LOS, -1),
                                      LOS_unc=self.LOS_unc, lkv_E=self.lkv_E, lkv_N=self.lkv_N,
                                      lkv_U=self.lkv_U, starttime=self.starttime, endtime=self.endtime)
        return new_InSAR_obj

    def subtract_value(self, value):
        new_InSAR_obj = Insar2dObject(lon=self.lon, lat=self.lat, LOS=np.subtract(self.LOS, value),
                                      LOS_unc=self.LOS_unc, lkv_E=self.lkv_E, lkv_N=self.lkv_N,
                                      lkv_U=self.lkv_U, starttime=self.starttime, endtime=self.endtime)
        return new_InSAR_obj

    def subtract_reference(self, refidx, tolerance=0.005):
        """
        Subtract the value of the reference pixel, essentially creating a referenced object.

        :param refidx: [int, int] for row/col of reference pixel or [float, float] for lon/lat
        :param tolerance: how close must the reference be to a viable pixel? In Degrees
        """
        if type(refidx[0]) is float:  # if refidx is lon, lat
            idx_lon, idx_lat = general_utils.get_nearest_pixel_in_geocoded_array(self.lon, self.lat, refidx[0],
                                                                                 refidx[1], min_dist_cutoff=tolerance)
            refvalue = self.LOS[idx_lon][idx_lat]
        else:  # if refidx is row, col
            refvalue = self.LOS[refidx[0]][refidx[1]]
        new_InSAR_obj = self.subtract_value(refvalue)
        return new_InSAR_obj

    def rewrap_InSAR(self, wavelength):
        """
        Take unwrapped LOS measurements (mm) and artificially wrap them around a certain radar wavelength (mm)

        :param wavelength: float, mm
        """
        rewrapped = general_utils.wrap_float(self.LOS, wavelength)
        new_InSAR_obj = Insar2dObject(lon=self.lon, lat=self.lat, LOS=rewrapped,
                                      LOS_unc=self.LOS_unc, lkv_E=self.lkv_E, lkv_N=self.lkv_N,
                                      lkv_U=self.lkv_U, starttime=self.starttime, endtime=self.endtime)
        return new_InSAR_obj

    def get_look_vector_at_point(self, target_lon, target_lat):
        """
        Return 3 component look vector, flight direction, and incidence angle at target location.
        Look vector is from ground to platform.

        :param target_lon: float, value of longitude
        :param target_lat: float, value of latitude
        :returns: E, N, U, flight, inc
        """
        # extract the look vector at a given spot:
        colnum = (np.abs(self.lon - target_lon)).argmin()
        rownum = (np.abs(self.lat - target_lat)).argmin()
        E, N, U = self.lkv_E[rownum][colnum], self.lkv_N[rownum][colnum], self.lkv_U[rownum][colnum]
        flight, inc = ivf.look_vector2flight_incidence_angles(E, N, U)
        return E, N, U, flight, inc

    def get_azimuth_incidence_grid(self, lookdir='right'):
        """
        Compute incidence angle across the grid, using the 3-component look vector.
        We assume the look direction is 'right', but if it's left, we flip the flight heading.
        """
        az, inc = ivf.look_vector2flight_incidence_angles(self.lkv_E, self.lkv_N, self.lkv_U)
        if lookdir == 'left':
            az = np.subtract(az, 180)
        return az, inc
