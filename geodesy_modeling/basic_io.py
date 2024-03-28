
import shapefile


def extract_polyline_from_shapefile(s):
    """
    :param s: shapefile.shape
    returns: list of lon, list of lat
    """
    return [coord[0] for coord in s.points], [coord[1] for coord in s.points]


def read_one_fault_from_shapefile(shape_file, fault_index):
    """
    A wrapper for extracting lon/lat from a fault in a shapefile.
    One could also extract shapes = sf.shapes()
    """
    sf = shapefile.Reader(shape_file)
    lonarray1, latarray1 = extract_polyline_from_shapefile(sf.shape(fault_index))
    return lonarray1, latarray1
