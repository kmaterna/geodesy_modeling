#!/usr/bin/env python

"""
A cute little script starting off the process of inverting gnss data in the Salton Sea for slip on the Kalin Fault
"""

import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import Tectonic_Utilities.Tectonic_Utils.seismo.moment_calculations as mo
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import argparse, json, subprocess


def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''');
    p.add_argument('--outdir', type=str, help='''Output directory''');
    exp_dict = vars(p.parse_args())
    exp_dict["obs_disp_points"] = "ssgf_vectors_manual_m.txt";
    exp_dict["fault_file"] = "../../_Data/Lohman_Fault_Geom/forK.mat";
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4);
    return exp_dict;


def read_gf_elements_kalin(exp_dict, obs_disp_pts):
    GF_elements = [];
    fault_tris = fst.io_other.read_brawley_lohman_2005(exp_dict['fault_file']);
    for tri in fault_tris:
        changed_slip = fst.fault_slip_triangle.change_fault_slip(tri, rtlat=1, dipslip=0, tensile=0);
        changed_slip = fst.fault_slip_triangle.change_reference_loc(changed_slip);
        model_disp_pts = fst.triangle_okada.compute_disp_points_from_triangles([changed_slip], obs_disp_pts,
                                                                               poisson_ratio=0.25);
        GF_elements.append(inv_tools.GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                                fault_name='kalin', lower_bound=-1, upper_bound=0, points=(),
                                                slip_penalty=0))
    return GF_elements;


if __name__ == "__main__":
    exp_dict = configure();
    obs_disp_pts = PyCoulomb.io_additionals.read_disp_points(exp_dict["obs_disp_points"]);
    GF_elements = read_gf_elements_kalin(exp_dict, obs_disp_pts);

    # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = [];
    for gf in GF_elements:
        G_one_col = inv_tools.buildG_column(gf.disp_points, obs_disp_pts);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);
    obs, sigmas = inv_tools.build_obs_vector(obs_disp_pts);
    G /= sigmas[:, None];
    weighted_obs = obs / sigmas;
    G, weighted_obs, sigmas = inv_tools.build_smoothing(GF_elements, ('kalin',), exp_dict["smoothing"], G, weighted_obs,
                                                        sigmas);
    plt.imshow(G, vmin=-3, vmax=3); plt.savefig(exp_dict['outdir']+"/G_matrix.png");

    # Money line: Constrained inversion
    lb = [x.lower_bound for x in GF_elements];
    ub = [x.upper_bound for x in GF_elements];
    response = scipy.optimize.lsq_linear(G, weighted_obs, bounds=(lb, ub), max_iter=1500, method='bvls');
    M_opt = response.x;  # parameters of best-fitting model

    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_disp_pts);
    resid = dpo.utilities.subtract_disp_points(obs_disp_pts, model_disp_pts);   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, exp_dict['outdir']+'/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(obs_disp_pts, exp_dict['outdir'] + '/obs_file.txt');

    # Unpack into a collection of fault triangles with optimal slip values
    modeled_faults = [];
    for i in range(len(GF_elements)):
        [new_fault] = GF_elements[i].fault_dict_list;
        new_fault = fst.fault_slip_triangle.change_fault_slip(new_fault, rtlat=M_opt[i]);
        modeled_faults.append(new_fault);
    total_moment = fst.fault_slip_triangle.get_total_moment(modeled_faults)

    fst.fault_slip_triangle.write_gmt_plots_geographic(modeled_faults, exp_dict['outdir']+"/temp-outfile.txt",
                                                       plotting_function=fst.fault_slip_triangle.get_rtlat_slip);
    fso.plot_fault_slip.map_source_slip_distribution([], exp_dict['outdir']+'/model_disps.png',
                                                     disp_points=model_disp_pts,
                                                     fault_traces_from_file=exp_dict['outdir']+"/temp-outfile.txt");
    fso.plot_fault_slip.map_source_slip_distribution([], exp_dict['outdir']+'/obs_disps.png', disp_points=obs_disp_pts);

    [all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_disp_pts, model_disp_pts);
    with open(exp_dict["outdir"]+'/metrics.txt', 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm);
        print("total moment (N-m): ", total_moment);
        print("Equivalent to:", mo.mw_from_moment(total_moment));
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm);
        ofile.write("total moment (N-m): %f\n" % total_moment);
        ofile.write("Equivalent to: %f\n" % mo.mw_from_moment(total_moment));
