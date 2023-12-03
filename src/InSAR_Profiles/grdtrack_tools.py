"""
Tools to manipulate profiles from GMT GRDTRACK
"""

import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from Tectonic_Utils.geodesy import haversine
from .. import general_utils


class GMTProfile:
    """The internal object for a single GMT GRDTRACK profile."""
    def __init__(self, center_lon, center_lat, azimuth, pt_lons, pt_lats, pt_dist, pt_azimuth, pt_value, offset=0):
        self.center_lon = center_lon  # float
        self.center_lat = center_lat  # float
        self.azimuth = azimuth  # float
        self.pt_lons = pt_lons  # 1d array
        self.pt_lats = pt_lats  # 1d array
        self.pt_dist = pt_dist  # 1d array
        self.pt_azimuth = pt_azimuth  # 1d array
        self.pt_value = pt_value  # 1d array
        self.offset = offset  # float.  essentially a type of metadata associated with each profile.

    def set_offset(self, offset):
        self.offset = offset

    def is_surface_breaking(self, critical_distance=0.25, critical_slope=0.30):
        """Is the slope at the center region of the profile sharp enough to be considered surface-breaking?"""
        idx = np.where(np.abs(self.pt_dist) < critical_distance)
        central_region = self.pt_value[idx]
        derivative = np.diff(central_region)
        if sum(np.isnan(derivative)) == len(derivative):
            return 0
        if np.nanmax(np.abs(derivative)) > critical_slope:
            return 1
        else:
            return 0

    def get_offsets_surface_breaking_simple(self, medfilt_window=101, neg_sample_idx=-25, pos_sample_idx=25):
        """The methods for determining slip amount when surface does break."""
        negative_side = self.pt_value[np.where(self.pt_dist < 0)]
        positive_side = self.pt_value[np.where(self.pt_dist > 0)]
        negative_filt = scipy.signal.medfilt(negative_side, medfilt_window)
        positive_filt = scipy.signal.medfilt(positive_side, medfilt_window)
        offset = negative_filt[neg_sample_idx] - positive_filt[pos_sample_idx]
        return offset

    def get_offsets_subsurface_profile_simple(self, medfilt_window=101, critical_distance=0.22):
        """The methods for determining slip amount when surface doesn't break."""
        filtered = scipy.signal.medfilt(self.pt_value, medfilt_window)  # FILTER THE PROFILE
        filtered = general_utils.detrend_signal(self.pt_dist, filtered)  # THEN WE DE-TREND THE SIGNAL
        filtered_zoomed = filtered[np.where(np.abs(self.pt_dist) <= critical_distance)]
        if sum(np.isnan(filtered_zoomed)) == len(filtered_zoomed):
            return 0
        offset = np.nanmax(filtered_zoomed) - np.nanmin(filtered_zoomed)  # gives a reasonable fit.
        return offset


def read_gmt_profiles(infile):
    """
    Read a set of profiles automatically created from GMT GRDTRACK.  File-IO function.
    Returns a list of GMTProfile objects.
    """
    profiles = []
    pt_lons, pt_lats, pt_dist, pt_azimuth, pt_value = [], [], [], [], []
    with open(infile, 'r') as ifile:
        print("Reading file %s " % infile)
        for line in ifile:
            if '>' in line:
                temp = line.split()
                center_lon = float(temp[6].split('/')[0])
                center_lat = float(temp[6].split('/')[1])
                azimuth = float(temp[7].split('=')[1])
                if len(pt_lons) > 0:
                    new_profile = GMTProfile(center_lon=center_lon, center_lat=center_lat, azimuth=azimuth,
                                             pt_lons=np.array(pt_lons), pt_lats=np.array(pt_lats),
                                             pt_dist=np.array(pt_dist), pt_azimuth=np.array(pt_azimuth),
                                             pt_value=np.array(pt_value))
                    profiles.append(new_profile)
                    pt_lons, pt_lats, pt_dist, pt_azimuth, pt_value = [], [], [], [], []
            else:
                temp = line.split()
                pt_lons.append(float(temp[0]))
                pt_lats.append(float(temp[1]))
                pt_dist.append(float(temp[2]))
                pt_azimuth.append(float(temp[3]))
                pt_value.append(float(temp[4]))
    return profiles


def get_nearest_profile(profile_list, target_pt):
    """Return the profile whose center point is closest to a target point, such as a creep-meter location."""
    distance_list = []
    for item in profile_list:
        distance_list.append(haversine.distance((item.center_lat, item.center_lon), (target_pt[1], target_pt[0])))
    return profile_list[np.argmin(distance_list)], np.argmin(distance_list)


def visualize_offsets(profile_list, outfile, vmin=0, vmax=30):
    """Draw a plot with the fault trace color-coded by offset."""
    plt.figure(dpi=300)
    cm = plt.cm.get_cmap('hot')
    colors = [x.offset for x in profile_list]
    sc = plt.scatter([x.center_lon for x in profile_list], [x.center_lat for x in profile_list], c=colors,
                     vmin=vmin, vmax=vmax, cmap=cm)
    plt.colorbar(sc, label="LOS (mm)")
    plt.savefig(outfile)
    return


def visualize_profile_displacement(profile_list, outfile):
    if len(profile_list) > 0:
        plt.figure(dpi=300)
        for i, item in enumerate(profile_list):
            plt.plot(item.pt_dist, item.pt_value, '.')
        plt.savefig(outfile)
    return


def draw_surface_trace(profile_list, outfile, surface_breaking=()):
    plt.figure()
    plt.plot([x.center_lon for x in profile_list], [x.center_lat for x in profile_list], '.b', label='Profiles')
    plt.plot([x.center_lon for x in surface_breaking], [x.center_lat for x in surface_breaking], '.r',
             label='Surface breaking')
    plt.legend()
    plt.savefig(outfile)
    return


def write_profile_offsets(profile_list, outfile):
    """Write text files with user-calculated profile offsets."""
    print("Writing file %s " % outfile)
    ofile = open(outfile, 'w')
    ofile.write("# center_lon center_lat azimuth offset\n")
    for item in profile_list:
        ofile.write("%f %f %f %f\n" % (item.center_lon, item.center_lat, item.azimuth, item.offset) )
    ofile.close()
    return


def read_profile_offsets(infile):
    """Read text files with user-calculated profile offsets."""
    print("Reading file %s " % infile)
    profiles = []
    [lons, lats, azimuths, offsets] = np.loadtxt(infile, skiprows=1, unpack=True)
    for i in range(len(lons)):
        new_profile = GMTProfile(center_lon=lons[i], center_lat=lats[i], azimuth=azimuths[i], offset=offsets[i],
                                 pt_azimuth=(), pt_value=(), pt_lons=(), pt_lats=(), pt_dist=())
        profiles.append(new_profile)
    return profiles
