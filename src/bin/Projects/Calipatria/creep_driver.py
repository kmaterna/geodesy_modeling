#!/usr/bin/env python

"""
A cute little script starting off the process of inverting gnss data in the Salton Sea for slip
"""

import Geodesy_Modeling.src.InSAR_1D_Object as InSAR_1D
from Geodesy_Modeling.src.InSAR_1D_Object.class_model import InSAR_1D_Object
import Elastic_stresses_py.PyCoulomb as PyCoulomb
import Elastic_stresses_py.PyCoulomb.disp_points_object as dpo
import Elastic_stresses_py.PyCoulomb.fault_slip_object as fso
import Geodesy_Modeling.src.Inversion.inversion_tools as inv_tools
import Geodesy_Modeling.src.Inversion.GF_element.rw_insar_gfs as rw_gf
import Tectonic_Utilities.Tectonic_Utils.seismo.moment_calculations as mo
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import argparse, json, subprocess
from Geodesy_Modeling.src.Inversion.GF_element import GF_element

class inputs_obj:
    def __init__(self, zerolon, zerolat, source):
        self.zerolon = zerolon
        self.zerolat = zerolat
        self.source_object = source
        return;

class params_obj:
    def __init__(self, alpha, nu):
        self.alpha = alpha
        self.nu = nu

def configure():
    p = argparse.ArgumentParser(description='''Inversion of geodetic data''')
    p.add_argument('--smoothing', type=float, help='''strength of Laplacian smoothing constraint''');
    p.add_argument('--outdir', type=str, help='''Output directory''');
    exp_dict = vars(p.parse_args())
    exp_dict["obs_desc"] = "../_4_Downsample_unw/igram_sum/pixels_filtered_m.txt";
    exp_dict["fault_file"] = "../Get_Fault_Model/model_fault_patches.txt";
    exp_dict["smoothing_length"] = 6;  # smooth adjacent patches with some wiggle room
    subprocess.call(['mkdir', '-p', exp_dict["outdir"]]);  # """Set up an experiment directory."""
    with open(exp_dict["outdir"] + "/configs_used.txt", 'w') as fp:
        json.dump(exp_dict, fp, indent=4);
    return exp_dict;

def compute_insar_gf_elements(fault_file: str, insar_object: InSAR_1D_Object):
    """Create a list of fault elements and their associated InSAR GF's for use in inversion. This is cool! """
    GF_elements = [];
    fault_patches = fso.file_io.io_slippy.read_slippy_distribution(fault_file);
    all_disp_points = insar_object.get_disp_points();   # get the locations of InSAR points

    for patch in fault_patches:
        changed_slip = patch.change_fault_slip(new_slip=1, new_rake=180, new_tensile=0);
        pycoulomb_fault = changed_slip.fault_object_to_coulomb_fault(-115.8, 33.1)
        inputs = inputs_obj(zerolon=-115.8, zerolat=33.1, source=[pycoulomb_fault]);
        params = params_obj(2/3, 0.25);
        model_disp_pts = PyCoulomb.run_dc3d.compute_ll_def(inputs, params, all_disp_points);

        # PROJECT 3D DISPLACEMENTS INTO LOS
        los_defo = [];
        for i in range(len(all_disp_points)):
            x = model_disp_pts[i].project_into_los(insar_object.lkv_E[i], insar_object.lkv_N[i], insar_object.lkv_U[i]);
            los_defo.append(x);

        model_disp_pts = dpo.utilities.set_east(model_disp_pts, los_defo);
        model_disp_pts = dpo.utilities.mult_disp_points_by(model_disp_pts, -1);  # sign convention
        GF_elements.append(GF_element.GF_element(disp_points=model_disp_pts, fault_dict_list=[changed_slip], units='m',
                                                 param_name='shf'));
    print("Computed Green's functions for %d patches" % len(GF_elements));
    return GF_elements;

def read_gf_elements(fault_file: str, gf_file: str):
    """Read a list of fault triangle elements and their associated GF's for use in inversion. """
    print("Reading pre-computed InSAR Green's functions.");
    fault_patches = fso.file_io.io_slippy.read_slippy_distribution(fault_file);
    GF_elements = rw_gf.read_insar_greens_functions(gf_file, fault_patches, param_name='shf', lower_bound=0,
                                                    upper_bound=0.3);
    return GF_elements;

def write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts):
    [_all_L2_norm, avg_misfit_norm, _, _] = dpo.compute_rms.obs_vs_model_L2_misfit(obs_disp_pts, model_disp_pts);
    with open(exp_dict["outdir"]+'/metrics.txt', 'w') as ofile:
        print('Avg misfit (mm):', avg_misfit_norm);
        print("total moment (N-m): ", total_moment);
        print("Equivalent to Mw:", mo.mw_from_moment(total_moment));
        ofile.write('Avg misfit: %f mm\n' % avg_misfit_norm);
        ofile.write("total moment (N-m): %f\n" % total_moment);
        ofile.write("Equivalent to Mw: %f\n" % mo.mw_from_moment(total_moment));
    return;


if __name__ == "__main__":
    exp_dict = configure();
    outdir = exp_dict['outdir'];
    desc_insar = InSAR_1D.inputs.inputs_txt(exp_dict['obs_desc']);  # get list of insar data
    obs_disp_pts = desc_insar.get_disp_points();  # get locations and data of InSAR points

    # ONE TIME: Compute the Green's Functions
    # GF_elements = compute_insar_gf_elements(exp_dict['fault_file'], desc_insar);
    # rw_gf.write_insar_greens_functions(GF_elements, "desc_insar_gfs.txt");

    GF_elements = read_gf_elements(exp_dict['fault_file'], "desc_insar_gfs.txt");  # get GFs

    # COMPUTE STAGE: INVERSE.
    list_of_gf_columns = [];
    for gf in GF_elements:
        G_one_col = inv_tools.buildG_column(gf.disp_points, obs_disp_pts);  # for one fault model parameter
        list_of_gf_columns.append(G_one_col);
    G = np.concatenate(tuple(list_of_gf_columns), axis=1);
    obs, sigmas = inv_tools.build_obs_vector(obs_disp_pts);
    G /= sigmas[:, None];
    w_obs = obs / sigmas;
    G, w_obs, sigmas = inv_tools.build_smoothing(GF_elements, ('shf',), exp_dict["smoothing"],
                                                 exp_dict["smoothing_length"], G, w_obs, sigmas);
    plt.imshow(G, vmin=-3, vmax=3); plt.savefig(outdir+"/G_matrix.png");

    # Money line: Constrained inversion
    lb, ub = [x.lower_bound for x in GF_elements], [x.upper_bound for x in GF_elements];
    response = scipy.optimize.lsq_linear(G, w_obs, bounds=(lb, ub), max_iter=1500, method='bvls');
    M_opt = response.x;  # parameters of best-fitting model

    model_disp_pts = inv_tools.forward_disp_points_predictions(G, M_opt, sigmas, obs_disp_pts);
    resid = dpo.utilities.subtract_disp_points(obs_disp_pts, model_disp_pts);   # make residual points
    PyCoulomb.io_additionals.write_disp_points_results(model_disp_pts, outdir + '/model_pred_file.txt');
    PyCoulomb.io_additionals.write_disp_points_results(obs_disp_pts, outdir + '/obs_file.txt');

    # Unpack into a collection of fault patches with optimal slip values
    modeled_faults = [];
    for i in range(len(GF_elements)):
        [new_fault] = GF_elements[i].fault_dict_list;
        new_fault = new_fault.change_fault_slip(new_slip=M_opt[i]);
        modeled_faults.append(new_fault);
    total_moment = fso.fault_slip_object.get_total_moment(modeled_faults)

    fso.fault_slip_object.write_gmt_vertical_fault_file(modeled_faults, outdir+"/vertical-outfile.txt",
                                                        color_mappable=lambda x: x.get_rtlat_slip());
    fso.file_io.io_slippy.write_slippy_distribution(modeled_faults, outdir+'/model_fault_patches.txt');
    for x in model_disp_pts:
        x.set_vert_value(x.dE_obs);
    for x in obs_disp_pts:
        x.set_vert_value(x.dE_obs);
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/model_disps.png', disp_points=model_disp_pts,
                                                     fault_traces_from_dict=modeled_faults, vmin=-0.015, vmax=0.015);
    fso.plot_fault_slip.map_source_slip_distribution([], outdir+'/obs_disps.png', disp_points=obs_disp_pts,
                                                     fault_traces_from_dict=modeled_faults, vmin=-0.015, vmax=0.015);
    write_misfit_report(exp_dict, obs_disp_pts, model_disp_pts);
    subprocess.call(['./gmt_plot.sh', outdir]);   # plotting the results in pretty GMT plots
