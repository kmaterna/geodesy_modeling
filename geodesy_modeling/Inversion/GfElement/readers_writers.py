from elastic_stresses_py.PyCoulomb import fault_slip_object as fso
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
import elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
from .GfElement import GfElement
import scipy.io


def write_csz_dist_fault_patches(gf_elements, model_results_vector, outfile_gmt, outfile_txt):
    """Write out slip results for a distributed CSZ model into GMT format"""
    modeled_slip_patches = []
    fault_dict_lists = [item.fault_dict_list for item in gf_elements]
    for i in range(len(gf_elements)):
        if gf_elements[i].param_name == 'CSZ_dist':   # CSZ patches are 1 patch per model element.
            new_patch = fault_dict_lists[i][0].change_fault_slip(model_results_vector[i]*10)  # mm/yr
            modeled_slip_patches.append(new_patch)
    if len(modeled_slip_patches) > 0:
        fso.file_io.outputs.write_gmt_fault_file(modeled_slip_patches, outfile_gmt,
                                                 color_mappable=lambda x: x.get_total_slip())

    fso.file_io.io_slippy.write_slippy_distribution(modeled_slip_patches, outfile_txt, slip_units='mm/yr')
    return


def read_distributed_GF_static1d(gf_file, geom_file, latlonfile, latlonbox=(-127, -120, 38, 52), unit_slip=False):
    """
    Read results of Fred's Static1D file (e.g., stat2C.outCascadia), getting into a minimal GfElement object.
    We also restrict the range of fault elements using a bounding box
    If unit_slip, we divide by the imposed slip rate to get a 1 cm/yr Green's Function.
    Returns a list of lists of disp_point objects, and a matching list of fault patches.
    We get into a minimal GfElement object the rest of the way research-specific code in the Humboldt driver.
    """
    fault_patches = fso.file_io.io_static1d.read_stat2C_geometry(geom_file)
    gps_disp_locs = fso.file_io.io_static1d.read_disp_points_from_static1d(latlonfile)

    print("Reading file %s " % gf_file)
    xdisps, ydisps, zdisps = [], [], []
    ifile = open(gf_file, 'r')
    for line in ifile:
        x = line[20:33]
        y = line[33:46]
        z = line[46:59]
        xdisps.append(float(x)/100)
        ydisps.append(float(y)/100)
        zdisps.append(float(z)/100)  # convert from cm to m
    ifile.close()

    counter = 0
    norm_factor = 1
    disp_points_all_patches, all_patches, given_slip = [], [], []
    GF_elements = []
    for i in range(len(fault_patches)):
        if not fault_patches[i].is_within_bbox(latlonbox):  # can restrict to the southern patches if desired
            counter = counter + len(gps_disp_locs)
            continue
        disp_points_one_patch = []
        for j in range(len(gps_disp_locs)):
            # Build a list of GF disp_points for each patch in inversion
            if unit_slip:
                norm_factor = 0.010 / fault_patches[i].slip  # normalizing to 1 cm/yr Green's Function
            else:
                norm_factor = 1
            disp_point = Displacement_points(lon=gps_disp_locs[j].lon, lat=gps_disp_locs[j].lat,
                                             dE_obs=-xdisps[counter] * norm_factor,  # negative means backslip
                                             dN_obs=-ydisps[counter] * norm_factor,
                                             dU_obs=-zdisps[counter] * norm_factor,
                                             Se_obs=0, Sn_obs=0, Su_obs=0, meas_type='model')
            counter = counter + 1
            disp_points_one_patch.append(disp_point)
            if counter == len(gps_disp_locs):
                break
        fault_slip_patch = fault_patches[i].change_fault_slip(fault_patches[i].slip * norm_factor)

        new_gf = GfElement(fault_dict_list=[fault_slip_patch], disp_points=disp_points_one_patch)
        GF_elements.append(new_gf)
        given_slip.append(fault_patches[i].slip)  # in mm

    return GF_elements, given_slip


def read_GFs_matlab_CSZ(gf_file):
    """
    Read the Green's functions for the CSZ calculated in Materna et al., 2019 by Noel Bartlow.
    Returns a list of Green's Functions elements.
    """

    print("Reading file %s " % gf_file)
    data_structure = scipy.io.loadmat(gf_file)  # a large dictionary object
    kern = data_structure['Kern']  # 165 x 303 (E, N, U for each grid element)
    fault_patches, nodes = fst.file_io.io_other.read_csz_bartlow_2019(gf_file)

    num_gf_elements = len(fault_patches)
    gf_elements = []
    for patch_number in range(num_gf_elements):
        model_disp_points = []
        for i in range(len(data_structure['Lons'])):
            new_item = Displacement_points(lon=data_structure['Lons'][i][0], lat=data_structure['Lats'][i][0],
                                           dE_obs=kern[(3*i)][patch_number], dN_obs=kern[(3*i)+1][patch_number],
                                           dU_obs=kern[(3*i)+2][patch_number])
            model_disp_points.append(new_item)  # Read model disp_points associated with one fault patch.
        fault_patches[patch_number] = fault_patches[patch_number].change_fault_slip(1.0, 0, 0)

        one_GF = GfElement(disp_points=model_disp_points, param_name=str(patch_number), units='meters',
                           fault_dict_list=[fault_patches[patch_number]])
        gf_elements.append(one_GF)
    return gf_elements
