# Get the locations of points adjacent to the fault, i.e., a synthetic creep-meter

import tectonic_utils.geodesy.haversine as haversine
import tectonic_utils.geodesy.fault_vector_functions as fault_vector_functions


def get_synthetic_creepmeter_points(lon0, lat0, strike, length_km, creepmeter_length):
    """Calculate the two points that would make a synthetic creepmeter orthogonal to a particular fault.
    The coordinates lon0 and lat0 are the back edge corner of the fault at the surface, not the center of the fault.
    The total length is creepmeter_length (given in m)"""

    strike_vector = fault_vector_functions.get_strike_vector(strike)
    west_vector = fault_vector_functions.get_strike_vector(strike - 90)
    east_vector = fault_vector_functions.get_strike_vector(strike + 90)

    center_lon, center_lat = haversine.add_vector_to_coords(lon0, lat0, (length_km/2) * strike_vector[0],
                                                            (length_km/2) * strike_vector[1])
    west_coord = haversine.add_vector_to_coords(center_lon, center_lat, dx=(creepmeter_length/2000) * west_vector[0],
                                                dy=(creepmeter_length/2000) * west_vector[1])
    east_coord = haversine.add_vector_to_coords(center_lon, center_lat, dx=(creepmeter_length/2000) * east_vector[0],
                                                dy=(creepmeter_length/2000) * east_vector[1])
    return west_coord, east_coord


def write_cm_points(west_coord, east_coord, filename):
    """Write the simple two points into a text file for espy calculations"""
    with open(filename, 'w') as ofile:
        ofile.write("%.9f %.9f \n" % (west_coord[0], west_coord[1]))
        ofile.write("%.9f %.9f \n" % (east_coord[0], east_coord[1]))
    return
