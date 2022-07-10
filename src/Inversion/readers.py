import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from Elastic_stresses_py.PyCoulomb import coulomb_collections as cc


def inside_lonlat_box(bbox, lonlat):
    """Return boolean condition testing whether [lon, lat] point is within a bounding box [w, e, s, n]"""
    if bbox[0] <= lonlat[0] <= bbox[1] and bbox[2] <= lonlat[1] <= bbox[3]:
        return 1;
    else:
        return 0;


def read_distributed_GF(gf_file, geom_file, latlonfile, latlonbox=(-127, -120, 38, 52), unit_slip=False):
    """
    Read the results of Fred's Static1D file
    For example: stat2C.outCascadia
    We also restrict the range of fault elements using a bounding box
    If unit_slip, we divide by the imposed slip rate to get a 1 cm/yr Green's Function.
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
        xdisps.append(float(x)/100);
        ydisps.append(float(y)/100);
        zdisps.append(float(z)/100);  # convert from cm to m
    ifile.close();

    counter = 0;
    norm_factor = 1;
    disp_points_all_patches, all_patches, given_slip = [], [], [];
    for i in range(len(fault_patches)):
        if not inside_lonlat_box(latlonbox, [fault_patches[i]["lon"], fault_patches[i]["lat"]]):  # southern patches
            counter = counter + len(gps_disp_locs);
            continue;
        disp_points_one_patch = [];
        for j in range(len(gps_disp_locs)):
            # Build a list of GF disp_points for each patch in inversion
            if unit_slip:
                norm_factor = 0.010 / fault_patches[i]["slip"];  # normalizing to 1 cm/yr Green's Function
            else:
                norm_factor = 1;
            disp_point = cc.Displacement_points(lon=gps_disp_locs[j].lon, lat=gps_disp_locs[j].lat,
                                                dE_obs=-xdisps[counter] * norm_factor,  # negative means backslip
                                                dN_obs=-ydisps[counter] * norm_factor,
                                                dU_obs=-zdisps[counter] * norm_factor,
                                                Se_obs=0, Sn_obs=0, Su_obs=0, name="", meas_type='model',
                                                starttime=None, endtime=None, refframe=None);
            counter = counter + 1;
            disp_points_one_patch.append(disp_point);
            if counter == len(gps_disp_locs):
                break;
        disp_points_all_patches.append(disp_points_one_patch);
        [fault_slip_patch] = fso.fault_slip_object.change_fault_slip([fault_patches[i]],
                                                                     fault_patches[i]["slip"] * norm_factor);
        all_patches.append(fault_slip_patch);
        given_slip.append(fault_patches[i]["slip"]);  # in mm

    return disp_points_all_patches, all_patches, given_slip;


def write_csz_dist_fault_patches(all_fault_patches, model_vector, outfile):
    """Write out the slip results for a distributed CSZ model into GMT format"""
    modeled_slip_patches = [];
    for i in range(len(all_fault_patches)):
        if len(all_fault_patches[i]) == 1:   # CSZ patches are 1 patch per model element.
            [new_patch] = fso.fault_slip_object.change_fault_slip(all_fault_patches[i], model_vector[i]*10);  # mm/yr
            modeled_slip_patches.append(new_patch);
    if len(modeled_slip_patches) > 0:
        fso.fault_slip_object.write_gmt_fault_file(modeled_slip_patches, outfile);
    return;
