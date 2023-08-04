import numpy as np
from Elastic_stresses_py.PyCoulomb import fault_slip_object as fso, disp_points_object as dpo
from Elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
from src.Inversion.GF_element import GF_element


def write_csz_dist_fault_patches(gf_elements, model_results_vector, outfile_gmt, outfile_txt):
    """Write out slip results for a distributed CSZ model into GMT format"""
    modeled_slip_patches = [];
    fault_dict_lists = [item.fault_dict_list for item in gf_elements];
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name == 'CSZ_dist':   # CSZ patches are 1 patch per model element.
            new_patch = fault_dict_lists[i][0].change_fault_slip(model_results_vector[i]*10);  # mm/yr
            modeled_slip_patches.append(new_patch);
    if len(modeled_slip_patches) > 0:
        fso.fault_slip_object.write_gmt_fault_file(modeled_slip_patches, outfile_gmt,
                                                   color_mappable=lambda x: x.get_total_slip());

    fso.file_io.io_slippy.write_slippy_distribution(modeled_slip_patches, outfile_txt, slip_units='mm/yr');
    return;


def read_distributed_GF_static1d(gf_file, geom_file, latlonfile, latlonbox=(-127, -120, 38, 52), unit_slip=False):
    """
    Read results of Fred's Static1D file (e.g., stat2C.outCascadia), getting partially into the GF_element object.
    We also restrict the range of fault elements using a bounding box
    If unit_slip, we divide by the imposed slip rate to get a 1 cm/yr Green's Function.
    Returns a list of lists of disp_point objects, and a matching list of fault patches.
    We get partially into the GF_element object; the rest of the way research-specific code in the Humboldt driver.
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


def read_insar_greens_functions(gf_file, fault_patches, param_name='', lower_bound=0, upper_bound=0):
    """
    Read pre-computed green's functions in the format matching the write-function.
    Currently, only works for triangles.

    :param gf_file: string, filename
    :param fault_patches: list of fault_slip_objects or fault_slip_triangles
    :param param_name: string
    :param lower_bound: float
    :param upper_bound: float
    """
    GF_elements = [];
    gf_data_array = np.loadtxt(gf_file);
    lons, lats = gf_data_array[:, 0], gf_data_array[:, 1];
    model_disp_pts = [];
    for tlon, tlat in zip(lons, lats):
        mdp = Displacement_points(lon=tlon, lat=tlat, dE_obs=0, dN_obs=0, dU_obs=0, Se_obs=0, Sn_obs=0, Su_obs=0,
                                  meas_type='insar');
        model_disp_pts.append(mdp);

    for i, tri in enumerate(fault_patches):  # right now, only the triangle version is written.
        changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0);  # triangle-specific
        changed_slip = changed_slip.change_reference_loc();  # triangle-specific interface
        index = i+2;  # moving to the correct column in the GF file, skipping lon and lat.
        los_defo = gf_data_array[:, index];
        model_disp_pts = dpo.utilities.set_east(model_disp_pts, los_defo);
        GF_elements.append(GF_element.GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                                 param_name=param_name, lower_bound=lower_bound,
                                                 upper_bound=upper_bound, slip_penalty=0));
    return GF_elements;


def write_insar_greens_functions(GF_elements, outfile):
    """
    Serialize a bunch of InSAR Green's Functions into written text file, in meters, with rows for each InSAR point:
    For each observation point: lon, lat, dLOS1, dLOS2, dLOS3.... [n fault patches].

    :param GF_elements: list of GF_elements with InSAR displacements in disp_points.
    :param outfile: string
    """
    print("Writing file %s " % outfile);
    ofile = open(outfile, 'w');
    for i, pt in enumerate(GF_elements[0].disp_points):
        ofile.write('%f %f ' % (pt.lon, pt.lat));
        point_displacements = [GF_el.disp_points[i].dE_obs for GF_el in GF_elements];
        for x in point_displacements:
            ofile.write(str(x)+" ");
        ofile.write("\n");
    ofile.close();
    return;
