#!/usr/bin/env python

"""
A cute little script starting off the process of inverting gnss data in the Salton Sea for slip on the Kalin Fault
"""

import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import Tectonic_Utils.seismo.moment_calculations as mo
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import argparse, json, subprocess, os
from Geodesy_Modeling.src.Inversion.GF_element import GF_element

def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''')
    p.add_argument('--outdir', type=str, help='''Output directory''')
    exp_dict = vars(p.parse_args())
    exp_dict["obs_disp_points"] = "Input_Data/ssgf_vectors_manual_m.txt"
    exp_dict["fault_file"] = "../../../_Data/Lohman_Fault_Geom/forK.mat"
    exp_dict["smoothing_length"] = 2  # smooth adjacent patches with some wiggle room (2 for salton sea)
    os.makedirs(exp_dict['outdir'], exist_ok=True)  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4)
    return exp_dict


def read_gf_elements_kalin(exp_dict, obs_disp_pts):
    """Create a list of fault triangle elements and their associated GF's for use in inversion. """
    GF_elements = []
    fault_tris = fst.file_io.io_other.read_brawley_lohman_2005(exp_dict['fault_file'])
    for tri in fault_tris:
        changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0)
        changed_slip = changed_slip.change_reference_loc()
        model_disp_pts = fst.triangle_okada.compute_disp_points_from_triangles([changed_slip], obs_disp_pts,
                                                                               poisson_ratio=0.25)
        GF_elements.append(
            GF_element.GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                  param_name='kalin', lower_bound=-1, upper_bound=0, slip_penalty=0))
    return GF_elements


def write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts):
    [_all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_disp_pts, model_disp_pts)
    with open(exp_dict["outdir"]+'/metrics.txt', 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm)
        print("total moment (N-m): ", total_moment)
        print("Equivalent to:", mo.mw_from_moment(total_moment))
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm)
        ofile.write("total moment (N-m): %f\n" % total_moment)
        ofile.write("Equivalent to: %f\n" % mo.mw_from_moment(total_moment))
    return


if __name__ == "__main__":
    exp_dict = configure()
    outdir = exp_dict['outdir']
    obs_disp_pts = PyCoulomb.io_additionals.read_disp_points(exp_dict["obs_disp_points"])
    GF_elements = read_gf_elements_kalin(exp_dict, obs_disp_pts)

    # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = []
    for gf in GF_elements:
        G_one_col = inv_tools.buildG_column(gf.disp_points, obs_disp_pts)  # for one fault model parameter
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

    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_disp_pts)
    resid = dpo.utilities.subtract_disp_points(obs_disp_pts, model_disp_pts)   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, outdir + '/model_pred_file.txt')
    PyCoulomb.io_additionals.write_disp_points_results(obs_disp_pts, outdir + '/obs_file.txt')

    # Unpack into a collection of fault triangles with optimal slip values
    modeled_faults = []
    for i in range(len(GF_elements)):
        [new_fault] = GF_elements[i].fault_dict_list
        new_fault = new_fault.change_fault_slip(rtlat=M_opt[i])
        modeled_faults.append(new_fault)
    total_moment = fst.fault_slip_triangle.get_total_moment(modeled_faults)

    fst.file_io.tri_outputs.write_gmt_plots_geographic(modeled_faults, outdir+"/temp-outfile.txt",
                                                       color_mappable=lambda x: x.get_rtlat_slip())
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/model_disps.png', disp_points=model_disp_pts,
                                                     fault_traces_from_file=outdir+"/temp-outfile.txt")
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/obs_disps.png', disp_points=obs_disp_pts)
    fst.file_io.tri_outputs.write_gmt_vertical_fault_file(modeled_faults, outdir+'/vertical_fault.txt',
                                                          color_mappable=lambda x: x.get_rtlat_slip(), strike=225)

    write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts)
    subprocess.call(['./gmt_plot.sh', outdir])   # plotting the results in pretty GMT plots

    # Do a special experiment: while the fault elements are in computer memory, let's forward predict their
    # displacement everywhere.
    newfaults = []
    for fault in modeled_faults:  # have to change the reference location.
        newfaults.append(fault.change_reference_loc((-115.7, 33.1)))

    grid_pts = dpo.utilities.generate_grid_of_disp_points(-115.90, -115.2, 32.9, 33.3, 0.005, 0.005)
    model_grid = fst.triangle_okada.compute_disp_points_from_triangles(newfaults, grid_pts, poisson_ratio=0.25)
    PyCoulomb.io_additionals.write_disp_points_results(model_grid, outdir+'/modeled_grid_pts.txt')
