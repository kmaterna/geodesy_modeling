#!/usr/bin/env python

"""
Inverting gnss and InSAR data in the Salton Sea for slip on the Kalin Fault.
"""

import elastic_stresses_py.PyCoulomb as PyCoulomb
import elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
import elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import geodesy_modeling.Inversion.inversion_tools as inv_tools
import geodesy_modeling.Inversion.GfElement.GfElement as GF_element
import geodesy_modeling.Inversion.GfElement.rw_insar_gfs as rw_insar_gfs
from geodesy_modeling.datatypes.InSAR_1D_Object import Insar1dObject
import tectonic_utils.seismo.moment_calculations as mo
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import argparse
import json
import subprocess
import os


def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''')
    p.add_argument('--outdir', type=str, help='''Output directory''')
    one_exp_dict = vars(p.parse_args())
    one_exp_dict["obs_desc"] = "../../../InSAR_Exps/sample_InSAR_on_roads/output/avg_desc.txt"
    one_exp_dict["obs_asc"] = "../../../InSAR_Exps/sample_InSAR_on_roads/output/avg_asc.txt"
    one_exp_dict["fault_file"] = "../../../../_Data/Lohman_Fault_Geom/forK.mat"
    one_exp_dict["smoothing_length"] = 2  # smooth adjacent patches with some wiggle room (2 for salton sea)
    os.makedirs(one_exp_dict['outdir'], exist_ok=True)  # """Set up an experiment directory."""
    with open(one_exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(one_exp_dict, fp, indent=4)
    return one_exp_dict


def compute_insar_gf_elements_kalin(fault_file: str, insar_object: Insar1dObject):
    """Create a list of fault triangle elements and their associated InSAR GF's for use in inversion. """
    tri_gf_elements = []
    fault_tris = fst.file_io.io_other.read_brawley_lohman_2005(fault_file)
    all_disp_points = insar_object.get_disp_points()   # get the locations of InSAR points
    for tri in fault_tris:
        changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0)
        changed_slip = changed_slip.change_reference_loc()
        model_pts = fst.triangle_okada.compute_disp_points_from_triangles([changed_slip], all_disp_points,
                                                                          poisson_ratio=0.25)
        # PROJECT 3D DISPLACEMENTS INTO LOS
        los_defo = []
        for j in range(len(all_disp_points)):
            x = model_pts[j].project_into_los(insar_object.lkv_E[j], insar_object.lkv_N[j], insar_object.lkv_U[j])
            los_defo.append(x)

        model_pts = dpo.utilities.with_easts_as(model_pts, los_defo)
        model_pts = dpo.utilities.mult_disp_points_by(model_pts, -1)  # sign convention
        tri_gf_elements.append(
            GF_element.GfElement(disp_points=model_pts, fault_dict_list=[changed_slip], units='m',
                                 param_name='kalin', lower_bound=-1, upper_bound=0, slip_penalty=0))
    print("Computed Green's functions for %d triangles" % len(tri_gf_elements))
    return tri_gf_elements


def read_gf_elements_kalin(fault_file: str, gf_file: str):
    """Read a list of fault triangle elements and their associated GF's for use in inversion. """
    print("Reading pre-computed InSAR Green's functions.")
    fault_tris = fst.file_io.io_other.read_brawley_lohman_2005(fault_file)
    kalin_gf_elements = rw_insar_gfs.read_insar_greens_functions(gf_file, fault_tris, param_name='kalin',
                                                                 lower_bound=-1, upper_bound=0)
    return kalin_gf_elements


def write_misfit_report(outfile, obs_pts, model_pts):
    [_all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_pts, model_pts)
    with open(outfile, 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm)
        print("total moment (N-m): ", total_moment)
        print("Equivalent to:", mo.mw_from_moment(total_moment))
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm)
        ofile.write("total moment (N-m): %f\n" % total_moment)
        ofile.write("Equivalent to: %f\n" % mo.mw_from_moment(total_moment))
    return


def compute_all_the_GFS(desc_pts: Insar1dObject, asc_pts: Insar1dObject):
    """ One-time function to compute and write ascending and descending green's functions (takes a few minutes)."""
    GF_elements_descend = compute_insar_gf_elements_kalin(exp_dict['fault_file'], desc_pts)
    rw_insar_gfs.write_insar_greens_functions(GF_elements_descend, "desc_insar_gfs.txt")
    GF_elements_ascend = compute_insar_gf_elements_kalin(exp_dict['fault_file'], asc_pts)
    rw_insar_gfs.write_insar_greens_functions(GF_elements_ascend, "asc_insar_gfs.txt")
    return


def combine_two_matching_lists_of_GF_elements(gf1, gf2):
    new_gf_elements = []
    for i in range(len(gf1)):
        x = GF_element.GfElement(disp_points=gf1[i].disp_points + gf2[i].disp_points, param_name=gf1[i].param_name,
                                 fault_dict_list=gf1[i].fault_dict_list, upper_bound=gf1[i].upper_bound,
                                 lower_bound=gf1[i].lower_bound, slip_penalty=gf1[i].slip_penalty, units=gf1[i].units)
        new_gf_elements.append(x)
    return new_gf_elements


if __name__ == "__main__":
    exp_dict = configure()
    outdir = exp_dict['outdir']

    # INPUT STAGE:
    desc_insar = geodesy_modeling.datatypes.InSAR_1D_Object.inputs.inputs_txt(exp_dict['obs_desc'])  # get list of insar data
    asc_insar = geodesy_modeling.datatypes.InSAR_1D_Object.inputs.inputs_txt(exp_dict['obs_asc'])
    obs_disp_pts_desc = desc_insar.get_disp_points()  # get locations and data of InSAR points
    obs_disp_pts_asc = asc_insar.get_disp_points()
    GF_elements_desc = read_gf_elements_kalin(exp_dict['fault_file'], "desc_insar_gfs.txt")  # get GFs
    GF_elements_asc = read_gf_elements_kalin(exp_dict['fault_file'], "asc_insar_gfs.txt")

    # SWITCH: Determine which data goes inside the inversion
    GF_elements = combine_two_matching_lists_of_GF_elements(GF_elements_desc, GF_elements_asc)  # if multiple datasets
    obs_disp_pts = obs_disp_pts_desc + obs_disp_pts_asc  # if multiple datasets

    # # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = []
    for gf in GF_elements:
        G_one_col = inv_tools.buildG_column(gf.disp_points, gf.disp_points)  # for one fault model parameter
        list_of_gf_columns.append(G_one_col)
    G = np.concatenate(tuple(list_of_gf_columns), axis=1)
    obs, sigmas = inv_tools.build_obs_vector(obs_disp_pts)
    G /= sigmas[:, None]
    w_obs = obs / sigmas
    G, w_obs, sigmas = inv_tools.build_smoothing(GF_elements, ('kalin',), exp_dict["smoothing"],
                                                 exp_dict["smoothing_length"], G, w_obs, sigmas)
    plt.imshow(G, vmin=-3, vmax=3)
    plt.savefig(outdir+"/G_matrix.png")

    # Money line: Constrained inversion
    lb, ub = [x.lower_bound for x in GF_elements], [x.upper_bound for x in GF_elements]
    response = scipy.optimize.lsq_linear(G, w_obs, bounds=(lb, ub), max_iter=1500, method='bvls')
    M_opt = response.x  # parameters of best-fitting model
    #
    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_disp_pts)
    resid = dpo.utilities.subtract_disp_points(obs_disp_pts, model_disp_pts)   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, outdir + '/model_pred_file.txt')
    PyCoulomb.io_additionals.write_disp_points_results(obs_disp_pts, outdir + '/obs_file.txt')
    #
    # Unpack into a collection of fault triangles with optimal slip values
    modeled_faults = []
    for _i, gf_element in enumerate(GF_elements):
        [new_fault] = gf_element.fault_dict_list
        new_fault = new_fault.change_fault_slip(rtlat=M_opt[_i])
        modeled_faults.append(new_fault)
    total_moment = fst.fault_slip_triangle.get_total_moment(modeled_faults)

    fst.file_io.tri_outputs.write_gmt_plots_geographic(modeled_faults, outdir+"/temp-outfile.txt",
                                                       color_mappable=lambda x: x.get_rtlat_slip())
    # fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/model_disps.png',
    #                                                  fault_traces_from_file=outdir+"/temp-outfile.txt")
    fst.file_io.tri_outputs.write_gmt_vertical_fault_file(modeled_faults, outdir+'/vertical_fault.txt',
                                                          color_mappable=lambda x: x.get_rtlat_slip(), strike=225)

    write_misfit_report(outdir+'/metrics.txt', obs_disp_pts, model_disp_pts)
    subprocess.call(['./gmt_plot.sh', outdir])   # plotting the results in pretty GMT plots
