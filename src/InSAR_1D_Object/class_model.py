import numpy as np


class InSAR_1D_Object:
    """
    A generalized 1D InSAR format where all data fields are 1D-vectors.
    The vector entries represent one pixel each.
    Displacements are Displacements, not velocities (for now)
    """
    def __init__(self, lon, lat, LOS, LOS_unc, lkv_E, lkv_N, lkv_U, starttime, endtime):
        self.lon = lon;  # 1d vector
        self.lat = lat;  # 1d vector
        self.LOS = LOS;  # 1d vector of displacements, in mm
        self.LOS_unc = LOS_unc;   # 1d vector, in mm
        self.lkv_E = lkv_E;  # Ground to Satellite, 1d vector
        self.lkv_N = lkv_N;  # Ground to Satellite, 1d vector
        self.lkv_U = lkv_U;  # Ground to Satellite, 1d vector
        self.starttime = starttime;  # just metadata, datetime object
        self.endtime = endtime;  # just metadata, datetime object

    def impose_bounding_box(self, bbox=(-180, 180, -90, 90)):
        """Impose a bounding box on InSAR data"""
        lon, lat, LOS, LOS_unc, unit_E, unit_N, unit_U = [], [], [], [], [], [], [];
        for i in range(len(self.lon)):
            if bbox[0] <= self.lon[i] <= bbox[1] and bbox[2] <= self.lat[i] <= bbox[3]:
                lon.append(self.lon[i]);
                lat.append(self.lat[i]);
                LOS.append(self.LOS[i]);
                LOS_unc.append(self.LOS_unc[i]);
                unit_E.append(self.lkv_E[i]);
                unit_N.append(self.lkv_N[i]);
                unit_U.append(self.lkv_U[i]);
        newInSAR_obj = InSAR_1D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N,
                                       lkv_U=unit_U, starttime=self.starttime, endtime=self.endtime);
        return newInSAR_obj;

    def remove_nans(self, verbose=False):
        """Remove Nans from 1D InSAR object"""
        lon, lat, LOS, LOS_unc, unit_E, unit_N, unit_U = [], [], [], [], [], [], [];
        for i in range(len(self.lon)):
            if np.isnan(self.LOS[i]):
                continue;
            else:
                lon.append(self.lon[i]);
                lat.append(self.lat[i]);
                LOS.append(self.LOS[i]);
                LOS_unc.append(self.LOS_unc[i]);
                unit_E.append(self.lkv_E[i]);
                unit_N.append(self.lkv_N[i]);
                unit_U.append(self.lkv_U[i]);
        newInSAR_obj = InSAR_1D_Object(lon=lon, lat=lat, LOS=LOS, LOS_unc=LOS_unc, lkv_E=unit_E, lkv_N=unit_N,
                                       lkv_U=unit_U, starttime=self.starttime, endtime=self.endtime);
        if verbose:
            print("Removing nans from 1D InSAR obj. Starting with %d, ending with %d pixels" % (len(self.lon),
                                                                                                len(newInSAR_obj.lon)));
        return newInSAR_obj;

    def flip_los_sign(self):
        new_InSAR_obj = InSAR_1D_Object(lon=self.lon, lat=self.lat, LOS=np.multiply(self.LOS, -1), LOS_unc=self.LOS_unc,
                                        lkv_E=self.lkv_E, lkv_N=self.lkv_N, lkv_U=self.lkv_U, starttime=self.starttime,
                                        endtime=self.endtime);
        return new_InSAR_obj

    def get_average_los_within_box(self, target_lon, target_lat, averaging_window):
        """
        Averaging window in degrees.
        Search the averaging window in both directions from target loc, and average the data.
        """
        new_los, new_unc, new_lkvE, new_lkvN, new_lkvU = [], [], [], [], [];
        for i in range(len(self.lon)):
            if target_lon - averaging_window <= self.lon[i] <= target_lon + averaging_window:
                if target_lat - averaging_window <= self.lat[i] <= target_lat + averaging_window:
                    new_los.append(self.LOS[i]);
                    new_unc.append(self.LOS_unc[i]);
                    new_lkvE.append(self.lkv_E[i]);
                    new_lkvN.append(self.lkv_N[i]);
                    new_lkvU.append(self.lkv_U[i]);
        new_InSAR_obj = InSAR_1D_Object(lon=target_lon, lat=target_lat, LOS=np.nanmean(new_los),
                                        LOS_unc=np.nanmean(new_unc), lkv_E=np.nanmean(new_lkvE),
                                        lkv_N=np.nanmean(new_lkvN), lkv_U=np.nanmean(new_lkvU),
                                        starttime=self.starttime, endtime=self.endtime);
        return new_InSAR_obj;

    def get_coordinate_tuples(self):
        """
        Return a list of tuples containing (lon, lat) for each pixel in the 1d list of pixels.
        """
        tuple_list = [];
        for x, y in zip(self.lon, self.lat):
            tuple_list.append((x, y));
        return tuple_list;
