#!/usr/bin/env python

"""
Inverting gnss and InSAR data in the Salton Sea for slip on the Kalin Fault.
"""

import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import Geodesy_Modeling.src.InSAR_1D_Object as InSAR_1D
from Geodesy_Modeling.src.InSAR_1D_Object.class_model import InSAR_1D_Object
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
    exp_dict["obs_desc"] = "../../../InSAR_Exps/sample_InSAR_on_roads/output/avg_desc.txt";
    exp_dict["obs_asc"] = "../../../InSAR_Exps/sample_InSAR_on_roads/output/avg_asc.txt";
    exp_dict["fault_file"] = "../../../../_Data/Lohman_Fault_Geom/forK.mat";
    exp_dict["smoothing_length"] = 2;  # smooth adjacent patches with some wiggle room (2 for salton sea)
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4);
    return exp_dict;


def compute_insar_gf_elements_kalin(fault_file: str, insar_object: InSAR_1D_Object):
    """Create a list of fault triangle elements and their associated InSAR GF's for use in inversion. """
    GF_elements = [];
    fault_tris = fst.file_io.io_other.read_brawley_lohman_2005(fault_file);
    all_disp_points = insar_object.get_disp_points();   # get the locations of InSAR points
    for tri in fault_tris:
        changed_slip = tri.change_fault_slip(rtlat=1, dipslip=0, tensile=0);
        changed_slip = changed_slip.change_reference_loc();
        model_disp_pts = fst.triangle_okada.compute_disp_points_from_triangles([changed_slip], all_disp_points,
                                                                               poisson_ratio=0.25);
        # PROJECT 3D DISPLACEMENTS INTO LOS
        los_defo = [];
        for i in range(len(all_disp_points)):
            x = dpo.utilities.project_into_los(model_disp_pts[i], insar_object.lkv_E[i], insar_object.lkv_N[i],
                                               insar_object.lkv_U[i]);
            los_defo.append(x);

        model_disp_pts = dpo.utilities.set_east(model_disp_pts, los_defo);
        model_disp_pts = dpo.utilities.mult_disp_points_by(model_disp_pts, -1);  # sign convention
        GF_elements.append(inv_tools.GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                                param_name='kalin', lower_bound=-1, upper_bound=0, slip_penalty=0));
    print("Computed Green's functions for %d triangles" % len(GF_elements));
    return GF_elements;


def read_gf_elements_kalin(fault_file: str, gf_file: str):
    """Read a list of fault triangle elements and their associated GF's for use in inversion. """
    print("Reading pre-computed InSAR Green's functions.");
    fault_tris = fst.file_io.io_other.read_brawley_lohman_2005(fault_file);
    GF_elements = inv_tools.read_insar_greens_functions(gf_file, fault_tris, param_name='kalin', lower_bound=-1,
                                                        upper_bound=0);
    return GF_elements;


def write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts):
    [_all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_disp_pts, model_disp_pts);
    with open(exp_dict["outdir"]+'/metrics.txt', 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm);
        print("total moment (N-m): ", total_moment);
        print("Equivalent to:", mo.mw_from_moment(total_moment));
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm);
        ofile.write("total moment (N-m): %f\n" % total_moment);
        ofile.write("Equivalent to: %f\n" % mo.mw_from_moment(total_moment));
    return;


if __name__ == "__main__":
    exp_dict = configure();
    outdir = exp_dict['outdir'];

    desc_insar = InSAR_1D.inputs.inputs_txt(exp_dict['obs_desc']);
    obs_disp_pts = desc_insar.get_disp_points();  # get the locations and data of InSAR points
    GF_elements_desc = compute_insar_gf_elements_kalin(exp_dict['fault_file'], desc_insar);
    inv_tools.write_insar_greens_functions(GF_elements_desc, "desc_insar_gfs.txt");
    GF_elements_desc = read_gf_elements_kalin(exp_dict['fault_file'], "desc_insar_gfs.txt");

    # asc_insar = InSAR_1D.inputs.inputs_txt(exp_dict['obs_asc']);
    # obs_disp_pts = asc_insar.get_disp_points();  # get the locations and data of InSAR points
    # GF_elements_asc = compute_insar_gf_elements_kalin(exp_dict['fault_file'], asc_insar);
    # inv_tools.write_insar_greens_functions(GF_elements_desc, "asc_insar_gfs.txt");
    # GF_elements_asc = read_gf_elements_kalin(exp_dict['fault_file'], "asc_insar_gfs.txt");

    # # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = [];
    for gf in GF_elements_desc:
        G_one_col = inv_tools.buildG_column(gf.disp_points, gf.disp_points);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);
    obs, sigmas = inv_tools.build_obs_vector(obs_disp_pts);
    G /= sigmas[:, None];
    w_obs = obs / sigmas;
    G, w_obs, sigmas = inv_tools.build_smoothing(GF_elements_desc, ('kalin',), exp_dict["smoothing"],
                                                 exp_dict["smoothing_length"], G, w_obs, sigmas);
    plt.imshow(G, vmin=-3, vmax=3); plt.savefig(outdir+"/G_matrix.png");

    # Money line: Constrained inversion
    lb, ub = [x.lower_bound for x in GF_elements_desc], [x.upper_bound for x in GF_elements_desc];
    response = scipy.optimize.lsq_linear(G, w_obs, bounds=(lb, ub), max_iter=1500, method='bvls');
    M_opt = response.x;  # parameters of best-fitting model
    #
    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_disp_pts);
    resid = dpo.utilities.subtract_disp_points(obs_disp_pts, model_disp_pts);   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, outdir + '/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(obs_disp_pts, outdir + '/obs_file.txt');
    #
    # Unpack into a collection of fault triangles with optimal slip values
    modeled_faults = [];
    for i in range(len(GF_elements_desc)):
        [new_fault] = GF_elements_desc[i].fault_dict_list;
        new_fault = new_fault.change_fault_slip(rtlat=M_opt[i]);
        modeled_faults.append(new_fault);
    total_moment = fst.fault_slip_triangle.get_total_moment(modeled_faults)

    fst.fault_slip_triangle.write_gmt_plots_geographic(modeled_faults, outdir+"/temp-outfile.txt",
                                                       color_mappable=lambda x: x.get_rtlat_slip());
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/model_disps.png',
                                                     fault_traces_from_file=outdir+"/temp-outfile.txt");
    fst.fault_slip_triangle.write_gmt_vertical_fault_file(modeled_faults, outdir+'/vertical_fault.txt',
                                                          color_mappable=lambda x: x.get_rtlat_slip(), strike=225);

    write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts);
    subprocess.call(['./gmt_plot.sh', outdir]);   # plotting the results in pretty GMT plots
