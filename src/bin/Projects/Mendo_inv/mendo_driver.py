#!/usr/bin/env python

from Geodesy_Modeling.src.Inversion.GF_element import readers_writers as gf_rw
from Elastic_stresses_py.PyCoulomb import disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Elastic_stresses_py.PyCoulomb.fault_slip_triangle as fst
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import json, argparse, subprocess
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt

filedict = {
    "gf_file": "../../../../GEOPHYS_DATA/Meshes/CSZ_Bartlow/2014_inversion_het.mat",
    "data_file": "../../stack/calculated_displacements.txt",
    "sample_gf_file": "../../Visualize_GFs/gf_pred_file.txt"
};
"""
Knobs that can be changed: 
    FIXED: Vertical uncertainty, Number of fault patches, Range of GNSS stations
    TUNABLE: Laplacian Smoothing, Laplacian Smoothing length, Tikhonov regularization
"""

def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''', default=20);
    p.add_argument('--outdir', type=str, help='''Output directory''', default='test_output/');
    p.add_argument('--tikhonov0', type=float, help='''strength of minimum-norm penalty''', default=30)
    exp_dict = vars(p.parse_args())
    exp_dict["smoothing_length"] = 10;  # smooth adjacent patches with some wiggle room
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4);
    return exp_dict;


def run_main():
    exp_dict = configure();
    outdir = exp_dict['outdir'];

    # Reading step
    GF_elements = gf_rw.read_GFs_matlab_CSZ(filedict['gf_file']);
    GF_elements = [x for x in GF_elements if x.fault_dict_list[0].depth > -50];
    for item in GF_elements:
        item.set_lower_bound(0); item.set_upper_bound(1);
        item.set_slip_penalty(exp_dict['tikhonov0'] - 0.1*item.fault_dict_list[0].depth);  # - item.fault_dict_list[0].depth a tunable parameter
    obs_data_points = dpo.io_gmt.read_disp_points_gmt(filedict['data_file']);
    obs_data_points = dpo.utilities.filter_to_remove_nans(obs_data_points);
    obs_data_points = dpo.utilities.filter_by_bounding_box(obs_data_points, (-126, -122, 39.65, 41.5));
    obs_data_points = dpo.utilities.filter_to_remove_outliers(obs_data_points, 0.02, verbose=True);  # remove P335/P794

    # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = [];
    obs_data_points, GF_elements = inv_tools.pair_gf_elements_with_obs(obs_data_points, GF_elements, tol=0.014);
    for gf in GF_elements:
        G_one_col = inv_tools.buildG_column(gf.disp_points, obs_data_points);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);
    obs, sigmas = inv_tools.build_obs_vector(obs_data_points);
    G /= sigmas[:, None];
    w_obs = obs / sigmas;
    smoothing_list = [x.param_name for x in GF_elements]
    G, w_obs, sigmas = inv_tools.build_smoothing(GF_elements, smoothing_list, exp_dict["smoothing"],
                                                 exp_dict["smoothing_length"], G, w_obs, sigmas);
    G, w_obs, sigmas = inv_tools.build_slip_penalty(GF_elements, 1, G, w_obs, sigmas);
    plt.imshow(G, vmin=-3, vmax=3); plt.savefig(outdir+"/G_matrix.png");

    # Money line: Constrained inversion
    lb, ub = [x.lower_bound for x in GF_elements], [x.upper_bound for x in GF_elements];
    response = scipy.optimize.lsq_linear(G, w_obs, bounds=(lb, ub), max_iter=1500, method='bvls');
    M_opt = response.x;  # parameters of best-fitting model

    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_data_points);
    _resid = dpo.utilities.subtract_disp_points(obs_data_points, model_disp_pts);   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, outdir + '/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(obs_data_points, outdir + '/obs_file.txt');

    # Unpack into a collection of fault triangles with optimal slip values
    modeled_faults = [];
    for i in range(len(GF_elements)):
        [new_fault] = GF_elements[i].fault_dict_list;
        new_fault = new_fault.change_fault_slip(rtlat=M_opt[i]);
        modeled_faults.append(new_fault);
    total_moment = fst.fault_slip_triangle.get_total_moment(modeled_faults)

    fst.fault_slip_triangle.write_gmt_plots_geographic(modeled_faults, outdir+"/temp-outfile.txt",
                                                       color_mappable=lambda x: x.get_rtlat_slip());
    scale_arrow = (1.0, 0.001, "1 mm");
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/model_disps.png', disp_points=model_disp_pts,
                                                     fault_traces_from_file=outdir+"/temp-outfile.txt",
                                                     region=(-125.0, -122, 39.05, 41.8), vert_mult=1000,
                                                     vert_disp_units='mm', vmin=-2, vmax=6, scale_arrow=scale_arrow,
                                                     slip_cbar_opts=(-0.06, 0.06, 0.001));
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/obs_disps.png', disp_points=obs_data_points,
                                                     region=(-125.0, -122, 39.05, 41.8), vert_mult=1000,
                                                     vert_disp_units='mm', vmin=-2, vmax=6, scale_arrow=scale_arrow);

    inv_tools.write_standard_misfit_report(exp_dict["outdir"], obs_data_points, model_disp_pts, total_moment);
    return;


if __name__ == "__main__":
    run_main();
