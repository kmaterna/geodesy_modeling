# Some useful functions for comparing pixels across object types
import numpy as np


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
