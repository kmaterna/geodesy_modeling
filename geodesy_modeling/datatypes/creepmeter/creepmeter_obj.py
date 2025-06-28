# Get the locations of points adjacent to the fault, i.e., a synthetic creep-meter

import tectonic_utils.geodesy.haversine as haversine
import tectonic_utils.geodesy.fault_vector_functions as fault_vector_functions
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points


class Synthetic_creepmeter:
    def __init__(self, center, west_point, east_point):
        self.center = center  # disp_points object
        self.west_point = west_point  # disp_points object
        self.east_point = east_point  # disp_points object
        self.displacement = get_synthetic_creepmeter_displacement([west_point, east_point])  # disp_points objects


def get_synthetic_creepmeter(lon0, lat0, strike, creepmeter_length):
    """Calculate the two points that would make a synthetic creepmeter orthogonal to a particular fault.
    The coordinates lon0 and lat0 are the desired coordinates of the creepmeter.
    The total length is creepmeter_length (given in m)

    :param lon0: lon of the creepmeter center
    :param lat0: lat of the creepmeter center
    :param strike: strike of the fault, in degrees CW from north
    :param creepmeter_length: meters
    :returns: a Synthetic_creepmeter object"""

    west_vector = fault_vector_functions.get_strike_vector(strike - 90)
    east_vector = fault_vector_functions.get_strike_vector(strike + 90)

    west_coord = haversine.add_vector_to_coords(lon0, lat0, dx=(creepmeter_length/2000) * west_vector[0],
                                                dy=(creepmeter_length/2000) * west_vector[1])
    east_coord = haversine.add_vector_to_coords(lon0, lat0, dx=(creepmeter_length/2000) * east_vector[0],
                                                dy=(creepmeter_length/2000) * east_vector[1])
    east_pt = Displacement_points(lon=east_coord[0], lat=east_coord[1])
    west_pt = Displacement_points(lon=west_coord[0], lat=west_coord[1])
    center_pt = Displacement_points(lon=lon0, lat=lat0)
    return Synthetic_creepmeter(center_pt, west_pt, east_pt)


def get_synthetic_creepmeter_displacement(disp_points):
    """
    :param disp_points: list of two disp_points objects with modeled displacements
    :return: displacement, in mm
    """
    # Returns synthetic creepmeter displacements in mm
    mag_v1 = disp_points[0].get_magnitude()
    mag_v2 = disp_points[1].get_magnitude()
    return 1000 * (mag_v1 + mag_v2)


def write_creepmeter_results(list_of_cm_objects, outfile):
    print("Writing synthetic creepmeter results in file %s " % outfile)
    with open(outfile, 'w') as ofile:
        ofile.write("# center_lon center_lat disp(mm)\n")
        for cm in list_of_cm_objects:
            ofile.write("%.9f %.9f %.7f\n" % (cm.center.lon, cm.center.lat, cm.displacement))
    return
