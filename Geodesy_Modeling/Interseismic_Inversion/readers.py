import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc


def inside_lonlat_box(bbox, lonlat):
    """Return boolean condition testing whether [lon, lat] point is within a bounding box [w, e, s, n]"""
    if bbox[0] <= lonlat[0] <= bbox[1] and bbox[2] <= lonlat[1] <= bbox[3]:
        return 1;
    else:
        return 0;


def read_distributed_GF(gf_file, geom_file, latlonfile, latlonbox=(-127, -120, 38, 43)):
    """
    Read the results of Fred's Static1D file
    For example: stat2C.outCascadia
    We also restrict the range of fault elements using a bounding box
    NEXT: Here we should divide by the imposed slip rate to get a 1cm/yr Green's Function
    Returns a list of lists of disp_point objects, and a matching list of fault patches
    """
    fault_patches = fso.io_static1d.read_stat2C_geometry(geom_file);
    gps_disp_locs = fso.io_static1d.read_disp_points_from_static1d(latlonfile);

    print("Reading file %s " % gf_file);
    xdisps, ydisps, zdisps = [], [], [];
    ifile = open(gf_file, 'r');
    for line in ifile:
        x = line[20:33];
        y = line[33:46];
        z = line[46:59];
        xdisps.append(float(x));
        ydisps.append(float(y));
        zdisps.append(float(z));
    ifile.close();

    counter = 0;
    disp_points_all_patches = [];
    all_patches = [];
    for i in range(len(fault_patches)):
        if not inside_lonlat_box(latlonbox, [fault_patches[i]["lon"], fault_patches[i]["lat"]]):  # southern patches
            counter = counter + len(gps_disp_locs);
            continue;
        disp_points_one_patch = [];
        for j in range(len(gps_disp_locs)):
            # Build a list of GF disp_points for each patch in inversion
            disp_point = cc.Displacement_points(lon=gps_disp_locs[j].lon, lat=gps_disp_locs[j].lat,
                                                dE_obs=xdisps[counter], dN_obs=ydisps[counter], dU_obs=zdisps[counter],
                                                Se_obs=0, Sn_obs=0, Su_obs=0, name="");
            counter = counter + 1;
            disp_points_one_patch.append(disp_point);
            if counter == len(gps_disp_locs):
                break;
        disp_points_all_patches.append(disp_points_one_patch);
        all_patches.append(fault_patches[i]);

    return disp_points_all_patches, all_patches;
