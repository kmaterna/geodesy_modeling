#!/usr/bin/env python

import os
import sys
from elastic_stresses_py.PyCoulomb import configure_calc, input_values, run_dc3d, output_manager
from geodesy_modeling.InSAR_2D_Object import model_enu_grids_into_los
from geodesy_modeling.Creepmeter import slip_1d_profiles, creepmeter_geometry


# 29 May 2024
# Steps:
# Read 1D depth profile; plot it
# Build dislocation file, assuming a certain length and the creepmeter is in the middle of the dislocation
# Run Okada
# Predict displacements at creepmeter site
# Predict InSAR displacements, project into LOS.
# Plot InSAR displacements.


def configure_experiment(txt_filename):
    _, filename = os.path.split(txt_filename)
    exp_name = filename.split('.')[0]
    os.makedirs("Outputs", exist_ok=True)
    os.makedirs(os.path.join("Outputs", exp_name), exist_ok=True)
    os.makedirs("espy_files", exist_ok=True)
    return exp_name, os.path.join("Outputs", exp_name)


def build_okada_files(txt_filename, expname):
    lon0, lat0, width_deg = -121.3, 36.7, 0.16
    strike, rake, dip, length_km = 321, 180, 89.9, 8
    tops, bottoms, slips = slip_1d_profiles.read_1d_profile(txt_filename)
    minlon, maxlon, minlat, maxlat = lon0-width_deg, lon0+width_deg, lat0-width_deg, lat0+width_deg
    espy_filename = os.path.join('espy_files', expname + '.intxt')
    with open(espy_filename, 'w') as ofile:
        print("Writing %s " % espy_filename)
        ofile.write("# Top line\n")
        ofile.write("General: 0.250 0.40 ")
        ofile.write(str(minlon)+" "+str(maxlon)+" "+str(lon0)+" "+str(minlat)+" "+str(maxlat)+" "+str(lat0)+"\n")
        for x1, x2, slip in zip(tops, bottoms, slips):
            # Format: strike rake dip length_km width_km lon lat depth_km slip_m (opt: tensile_m)
            width_km = x2 - x1
            ofile.write("Source_Patch: "+str(strike)+" "+str(rake)+" "+str(dip)+" "+str(length_km)+" ")
            ofile.write(str(width_km)+' '+str(lon0)+' '+str(lat0)+' '+str(x1)+' '+str(slip/1000)+' 0\n')

    # Make the synthetic creepmeter locations
    pt1, pt2 = creepmeter_geometry.get_synthetic_creepmeter_points(lon0, lat0, strike, length_km, creepmeter_length=10)
    creepmeter_geometry.write_cm_points(pt1, pt2, filename=os.path.join('espy_files', 'target_pts.txt'))
    return


def get_synthetic_creepmeter_magnitude(disp_points):
    # Returns synthetic creepmeter displacements in mm
    mag_v1 = disp_points[0].get_magnitude()
    mag_v2 = disp_points[1].get_magnitude()
    return 1000 * (mag_v1 + mag_v2)


def run_okada_calc(exp_name, outdir):
    espy_filename = os.path.join('espy_files', exp_name + '.intxt')
    params = configure_calc.Params(input_file=espy_filename,
                                   disp_points_file=os.path.join("espy_files", "target_pts.txt"),
                                   plot_stress=0, outdir=outdir)
    [inputs, obs_disp_points, obs_strain_points] = input_values.read_inputs(params)
    out_object = run_dc3d.do_stress_computation(params, inputs, obs_disp_points, obs_strain_points)
    output_manager.produce_outputs(params, inputs, obs_disp_points, obs_strain_points, out_object)
    synthetic_disp = get_synthetic_creepmeter_magnitude(out_object.model_disp_points)
    print("SYNTHETIC DISP = ", synthetic_disp, " mm")
    with open(os.path.join(outdir, 'synthetic_creepmeter.txt'), 'w') as ofile:
        ofile.write("SYNTHETIC DISP = " + str(synthetic_disp) + " mm\n")
    return out_object.model_disp_points


def project_to_los(outdir, disp_points=()):
    los_params = {  # REQUIRED PARAMETERS
              'east_grdfile': os.path.join(outdir, 'east.grd'),  # with units of meters by default
              'north_grdfile': os.path.join(outdir, 'north.grd'),  # with units of meters by default
              'up_grdfile': os.path.join(outdir, 'vert.grd'),  # with units of meters by default
              'flight_direction': 190.0,
              'incidence_angle': 37.5,
              'wavelength_mm': 56,  # with units of mm
              'look_direction': 'right',
              'plot_wrapped': True,                 # OPTIONAL PARAMETERS:
              'plot_unwrapped': True,
              'wrapped_plot_name': 'pred.jpg',
              'unwrapped_plot_name': 'unw_pred.jpg',
              'refidx': [0, 0],   # either [float, float] for lon, lat, or [int, int] for row, col
              'plot_range': None,
              'plot_wrapped_annot': None,
              'plot_unwrapped_annot': None,
              'outdir': outdir,
              'label': 's1_descending_'}
    model_enu_grids_into_los.do_synthetic_grid_LOS(los_params, disp_points=disp_points)

    los_params['flight_direction'] = 85
    los_params['incidence_angle'] = 52
    los_params['wavelength_mm'] = 240
    los_params['look_direction'] = 'left'
    los_params['label'] = 'uavsar_eastward_'
    model_enu_grids_into_los.do_synthetic_grid_LOS(los_params, disp_points=disp_points)

    los_params['flight_direction'] = 265
    los_params['incidence_angle'] = 52
    los_params['wavelength_mm'] = 240
    los_params['look_direction'] = 'left'
    los_params['label'] = 'uavsar_westward_'
    model_enu_grids_into_los.do_synthetic_grid_LOS(los_params, disp_points=disp_points)
    return


if __name__ == "__main__":
    given_filename = sys.argv[1]   # the path to a 1d slip-distribution file

    exper_name, out_dir = configure_experiment(given_filename)
    slip_1d_profiles.plot_1d_profile(given_filename, os.path.join(out_dir, 'slip_dist_' + exper_name + '.png'))
    build_okada_files(given_filename, exper_name)
    model_disp_points = run_okada_calc(exper_name, out_dir)
    project_to_los(out_dir, model_disp_points)
