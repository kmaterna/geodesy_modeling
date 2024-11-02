#!/usr/bin/env python

import random
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
from elastic_stresses_py.PyCoulomb import io_additionals


def do_main():
    pts = []
    lon0, lat0, width_deg = -121.28, 36.7, 0.16
    minlon, maxlon, minlat, maxlat = lon0 - width_deg, lon0 + width_deg, lat0 - width_deg, lat0 + width_deg
    for i in range(50):
        newlon = random.uniform(minlon, maxlon)
        newlat = random.uniform(minlat, maxlat)
        newpt = Displacement_points(lon=newlon, lat=newlat)
        pts.append(newpt)
    io_additionals.write_disp_points_locations(pts, "espy_files/gnss_points.txt")
    return


if __name__ == "__main__":
    do_main()
