"""
A small project to take a fault trace (lon/lat paris) and mesh it into chunks
The result will be a series of rectangular mesh elements.
The format will be the internal format used in PyCoulomb, a class with geometry and slip for each patch.
"""

import numpy as np
import matplotlib.pyplot as plt
from Tectonic_Utils.geodesy import haversine, fault_vector_functions
from elastic_stresses_py.PyCoulomb.fault_slip_object import fault_slip_object as fso


def read_surface_fault_trace(infile):
    """
    Requires a single continuous trace.
    Returns fault_trace = [lon, lat]
    Ordered from south to north.
    """
    lon, lat = np.loadtxt(infile, unpack=True)
    if lat[0] < lat[-1]:
        return [lon, lat]
    elif lat[0] > lat[-1]:
        return [np.flipud(lon), np.flipud(lat)]
    else:
        raise Exception("Error! Lat[0] = Lat[-1]. Don't know how to order this fault trace south-to-north. ")


def split_fault_trace(fault_trace, typical_spacing_km):
    """
    Algorithm: Walk up the fault in chunks.  If the chunk is smaller than typical_spacing_km, it becomes one segment.
    Otherwise we split the segment into something similar to typical_spacing_km.
    Returns: (starting lon, starting_lat, strike, length, ending_lon, ending_lat) for each fault trace segment
    """
    all_fault_segments = []
    for i in range(1, len(fault_trace[0])):
        start_point = (fault_trace[1][i-1], fault_trace[0][i-1])
        end_point = (fault_trace[1][i], fault_trace[0][i])
        segment_distance = haversine.distance(start_point, end_point)
        strike = haversine.calculate_initial_compass_bearing(start_point, end_point)

        if segment_distance < typical_spacing_km:    # for the really small segments
            one_fault_segment = (fault_trace[0][i-1], fault_trace[1][i-1], strike, segment_distance,
                                 fault_trace[0][i], fault_trace[1][i])
            all_fault_segments.append(one_fault_segment)

        else:   # the segments greater than 2x the characteristic spacing
            num_subsegments = int(np.ceil(segment_distance / typical_spacing_km))  # making segments < typical_spacing
            counter = [x for x in range(0, num_subsegments+1)]

            lon_difference = fault_trace[0][i] - fault_trace[0][i-1]
            lon_step = lon_difference / num_subsegments
            lon_subarray = [fault_trace[0][i-1] + lon_step*x for x in counter]
            # for x subsegments, we have x+1 elements in lon_subarray

            lat_difference = fault_trace[1][i] - fault_trace[1][i-1]
            lat_step = lat_difference / num_subsegments
            lat_subarray = [fault_trace[1][i-1] + lat_step*x for x in counter]

            for j in range(1, len(lon_subarray)):
                one_fault_segment = (lon_subarray[j-1], lat_subarray[j-1], strike, segment_distance/num_subsegments,
                                     lon_subarray[j], lat_subarray[j])
                all_fault_segments.append(one_fault_segment)
    print("Meshed fault trace of %d segments into %d segments " % (len(fault_trace[0]), len(all_fault_segments)))
    return all_fault_segments


def convert_2d_segments_to_fault_dictionary(fault_segments, dip, dip_direction, top_depth, bottom_depth, slip_cm, rake):
    """
    fault_segment includes: (starting lon, starting_lat, strike, length, ending_lon, ending_lat)
    where starting_lon is probably south of ending_lon.
    We are converting to the internal fault format from Elastic_stresses_py
    """
    fault_list = []
    for item in fault_segments:
        if dip_direction == 'south':  # Flip strike/dip if necessary
            strike = item[2] - 180
            # Do I flip the lon/lat to the other corner here too?  Not totally sure.  Should test this.
        else:
            strike = item[2]
        downdip_width = fault_vector_functions.get_downdip_width(top_depth, bottom_depth, dip)
        new_fault = fso.FaultSlipObject(strike=strike, dip=dip, length=item[3], rake=rake, slip=slip_cm/100,
                                        tensile=0, depth=top_depth, width=downdip_width, lon=item[0], lat=item[1])
        fault_list.append(new_fault)
    return fault_list


def plot_surface_fault_and_segments(surface_fault, all_fault_segments, outfile="surface_segments.png"):
    """ A little function to help see how you've done """
    plt.figure()
    plt.plot(surface_fault[0], surface_fault[1], '-k', linewidth=4)
    for segment in all_fault_segments:
        plt.plot([segment[0], segment[4]], [segment[1], segment[5]], linewidth=2, color='red')
    plt.title("Fault segments: %d " % len(all_fault_segments))
    plt.savefig(outfile)
    return
