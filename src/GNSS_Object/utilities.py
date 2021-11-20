"""
Operations on GNSS objects
"""

import numpy as np
from . import gnss_object

def add_gps_constant_offset(displacement_object, enu_constant_offset):
    """
    In case your reference station undergoes some offset, you might want to put that back
    to all data because it's been referenced-out already.
    Allows for more realistic GNSS inversions
    Offset and displacement obj in mm

    :param displacement_object: a list of gnss time series objects (in mm)
    :param enu_constant_offset: list of 3 numbers (in mm)
    :returns: a list of gnss time series objects
    """
    new_gps_displacements_object = [];
    for one_object in displacement_object:
        dE = [i + enu_constant_offset[0] for i in one_object.dE];
        dN = [i + enu_constant_offset[1] for i in one_object.dN];
        dU = [i + enu_constant_offset[2] for i in one_object.dU];
        object_after_offset = gnss_object.Timeseries(name=one_object.name, coords=one_object.coords,
                                                     dtarray=one_object.dtarray, dN=dN, dE=dE, dU=dU, Sn=one_object.Sn,
                                                     Se=one_object.Se, Su=one_object.Su, EQtimes=one_object.EQtimes);
        new_gps_displacements_object.append(object_after_offset);
    return new_gps_displacements_object;


def write_gps_invertible_format(gps_object_list, filename):
    """One header line
    GPS Data in meters"""
    print("Writing GPS displacements into file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# Header: lon, lat, dE, dN, dU, Se, Sn, Su (m)\n");
    for station in gps_object_list:
        if np.isnan(station.dE[1]):
            continue;
        else:
            ofile.write('%f %f ' % (station.coords[0], station.coords[1]));
            ofile.write("%f %f %f " % (0.001 * station.dE[1], 0.001 * station.dN[1], 0.001 * station.dU[1]));
            ofile.write("%f %f %f\n" % (station.Se[1], station.Sn[1], station.Su[1]));
    ofile.close();
    return;
