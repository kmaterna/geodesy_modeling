import numpy as np
from elastic_stresses_py.PyCoulomb import disp_points_object as dpo
from elastic_stresses_py.PyCoulomb.fault_slip_triangle import fault_slip_triangle
from elastic_stresses_py.PyCoulomb.disp_points_object.disp_points_object import Displacement_points
from .GfElement import GfElement


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
    GF_elements = []
    gf_data_array = np.loadtxt(gf_file)
    lons, lats = gf_data_array[:, 0], gf_data_array[:, 1]
    model_disp_pts = []
    for tlon, tlat in zip(lons, lats):
        model_disp_pts.append(Displacement_points(lon=tlon, lat=tlat, dE_obs=0, dN_obs=0, dU_obs=0, meas_type='insar'))

    if isinstance(fault_patches[0], fault_slip_triangle.TriangleFault):  # triangle version
        for i, tri in enumerate(fault_patches):
            changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0)  # triangle-specific
            changed_slip = changed_slip.change_reference_loc()  # triangle-specific interface
            index = i+2  # moving to the correct column in the GF file, skipping lon and lat.
            los_defo = gf_data_array[:, index]
            model_disp_pts = dpo.utilities.with_easts_as(model_disp_pts, los_defo)
            GF_elements.append(GfElement(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                         param_name=param_name, lower_bound=lower_bound, upper_bound=upper_bound))

    else:
        for i, patch in enumerate(fault_patches):  # Rectangular version
            changed_slip = patch.change_fault_slip(new_slip=1, new_rake=180, new_tensile=0)
            index = i+2  # moving to the correct column in the GF file, skipping lon and lat.
            los_defo = gf_data_array[:, index]
            model_disp_pts = dpo.utilities.with_easts_as(model_disp_pts, los_defo)
            GF_elements.append(GfElement(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                         param_name=param_name, lower_bound=lower_bound, upper_bound=upper_bound))
    return GF_elements


def write_insar_greens_functions(GF_elements, outfile):
    """
    Serialize a bunch of InSAR Green's Functions into written text file, in meters, with rows for each InSAR point:
    For each observation point: lon, lat, dLOS1, dLOS2, dLOS3.... [n fault patches].

    :param GF_elements: list of GF_elements with InSAR displacements in disp_points.
    :param outfile: string
    """
    print("Writing file %s " % outfile)
    ofile = open(outfile, 'w')
    for i, pt in enumerate(GF_elements[0].disp_points):
        ofile.write('%f %f ' % (pt.lon, pt.lat))
        point_displacements = [GF_el.disp_points[i].dE_obs for GF_el in GF_elements]
        for x in point_displacements:
            ofile.write(str(x)+" ")
        ofile.write("\n")
    ofile.close()
    return
