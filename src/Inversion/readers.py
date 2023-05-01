import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
from Elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points


def read_distributed_GF(gf_file, geom_file, latlonfile, latlonbox=(-127, -120, 38, 52), unit_slip=False):
    """
    Read results of Fred's Static1D file
    For example: stat2C.outCascadia
    We also restrict the range of fault elements using a bounding box
    If unit_slip, we divide by the imposed slip rate to get a 1 cm/yr Green's Function.
    Returns a list of lists of disp_point objects, and a matching list of fault patches
    """
    fault_patches = fso.file_io.io_static1d.read_stat2C_geometry(geom_file);
    gps_disp_locs = fso.file_io.io_static1d.read_disp_points_from_static1d(latlonfile);

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
        if not fault_patches[i].is_within_bbox(latlonbox):  # can restrict to the southern patches if desired
            counter = counter + len(gps_disp_locs);
            continue;
        disp_points_one_patch = [];
        for j in range(len(gps_disp_locs)):
            # Build a list of GF disp_points for each patch in inversion
            if unit_slip:
                norm_factor = 0.010 / fault_patches[i].slip;  # normalizing to 1 cm/yr Green's Function
            else:
                norm_factor = 1;
            disp_point = Displacement_points(lon=gps_disp_locs[j].lon, lat=gps_disp_locs[j].lat,
                                             dE_obs=-xdisps[counter] * norm_factor,  # negative means backslip
                                             dN_obs=-ydisps[counter] * norm_factor,
                                             dU_obs=-zdisps[counter] * norm_factor,
                                             Se_obs=0, Sn_obs=0, Su_obs=0, meas_type='model');
            counter = counter + 1;
            disp_points_one_patch.append(disp_point);
            if counter == len(gps_disp_locs):
                break;
        disp_points_all_patches.append(disp_points_one_patch);
        fault_slip_patch = fault_patches[i].change_fault_slip(fault_patches[i].slip * norm_factor);
        all_patches.append(fault_slip_patch);
        given_slip.append(fault_patches[i].slip);  # in mm

    return disp_points_all_patches, all_patches, given_slip;


def write_csz_dist_fault_patches(gf_elements, model_vector, outfile_gmt, outfile_txt):
    """Write out the slip results for a distributed CSZ model into GMT format"""
    modeled_slip_patches = [];
    fault_dict_lists = [item.fault_dict_list for item in gf_elements];
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name == 'CSZ_dist':   # CSZ patches are 1 patch per model element.
            new_patch = fault_dict_lists[i][0].change_fault_slip(model_vector[i]*10);  # mm/yr
            modeled_slip_patches.append(new_patch);
    if len(modeled_slip_patches) > 0:
        fso.fault_slip_object.write_gmt_fault_file(modeled_slip_patches, outfile_gmt,
                                                   color_mappable=lambda x: x.get_total_slip());

    fso.file_io.io_slippy.write_slippy_distribution(modeled_slip_patches, outfile_txt, slip_units='mm/yr');
    return;
