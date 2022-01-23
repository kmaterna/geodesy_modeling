"""
A few rare operations on GNSS objects. Utilities are mostly pushed into GNSS repo.
"""

import numpy as np


def write_gps_invertible_format(gps_object_list, filename):
    """Write displacements from a special format associated with an interval
    (i.e., disp object has only 2 elements, start=0 and finish=finish).
    Format has one header line. GPS displacements are in meters."""
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
