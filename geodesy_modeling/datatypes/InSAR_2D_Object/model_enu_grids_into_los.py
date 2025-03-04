"""
Project 3 grd files of modeled east/north/up deformation into the LOS of an assumed satellite flight.
Simple calculation - single incidence angle for now.
"""

import os


def parse_config(params_provided):
    defaults = {  # REQUIRED PARAMETERS
              'east_grdfile': 'east.grd',  # with units of meters by default
              'north_grdfile': 'north.grd',  # with units of meters by default
              'up_grdfile': 'vert.grd',  # with units of meters by default
              'flight_direction': 190.0,
              'incidence_angle': 37.5,
              'look_direction': 'right',
              'wavelength_mm': 56,  # with units of mm
              'plot_wrapped': True,                 # OPTIONAL PARAMETERS:
              'plot_unwrapped': True,
              'wrapped_plot_name': 'pred.jpg',
              'unwrapped_plot_name': 'unw_pred.jpg',
              'refidx': [0, 0],   # either [float, float] for lon, lat, or [int, int] for row, col
              'plot_range': None,
              'plot_wrapped_annot': None,
              'plot_unwrapped_annot': None,
              'outdir': '',
              'label': ''}

    # Combine dictionary with priority order using ** operator
    prio_dict = {1: params_provided, 2: defaults}
    params = {**prio_dict[2], **prio_dict[1]}
    return params


def read_grd_inputs(params):
    """
    Reading an InSAR_2D_Object
    """
    myobj = geodesy_modeling.datatypes.InSAR_2D_Object.inputs.inputs_from_synthetic_enu_grids(params['east_grdfile'], params['north_grdfile'],
                                                                                              params['up_grdfile'], params['flight_direction'],
                                                                                              params['incidence_angle'],
                                                                                              look_direction=params['look_direction'])
    return myobj


def plot_synthetic_grid_los(params, insarobj, disp_points=None, disp_points_color=None):
    """
    :param params: dictionary
    :param insarobj: a 2D InSAR object
    :param disp_points: a list of disp_points for optional plotting annotations
    :param disp_points_color: a 1d array of floats to be plotted as colors in the disp_points fill
    """
    if params['plot_unwrapped']:  # the outputs are in los deformation (mm)
        myobj_ref = insarobj.subtract_reference(params['refidx'])  # Subtract reference pix
        geodesy_modeling.datatypes.InSAR_2D_Object.outputs.write_InSAR2D(myobj_ref, os.path.join(params['outdir'],
                                                                                                 params['label'] +"unw_los.grd"))
        geodesy_modeling.datatypes.InSAR_2D_Object.outputs.map_wrapped_insar(os.path.join(params['outdir'], params['label'] + "unw_los.grd"),
                                                                             os.path.join(params['outdir'],
                                                               params['label']+params['unwrapped_plot_name']),
                                                                             text_annot=params['plot_unwrapped_annot'],
                                                                             flight_heading=params['flight_direction'],
                                                                             look_dir=params['look_direction'],
                                                                             disp_points=disp_points, region=params['plot_range'],
                                                                             refloc=params['refidx'], disp_points_color=disp_points_color)

    if params['plot_wrapped']:  # the outputs are in phase (radians)
        myobj_wrapped = insarobj.rewrap_InSAR(params['wavelength_mm'])
        geodesy_modeling.datatypes.InSAR_2D_Object.outputs.write_InSAR2D(myobj_wrapped, os.path.join(params['outdir'],
                                                                                                     params['label'] +"phase.grd"))
        geodesy_modeling.datatypes.InSAR_2D_Object.outputs.map_wrapped_insar(os.path.join(params['outdir'], params['label'] + "phase.grd"),
                                                                             os.path.join(params['outdir'],
                                                               params['label']+params["wrapped_plot_name"]),
                                                                             text_annot=params['plot_wrapped_annot'],
                                                                             look_dir=params['look_direction'],
                                                                             flight_heading=params['flight_direction'],
                                                                             disp_points=disp_points, region=params['plot_range'], wrapped=True)
    return


def do_synthetic_grid_LOS(params_provided, disp_points=None, disp_points_color=None):
    # SAMPLE DRIVER: CONFIG, INPUT, OUTPUT
    params = parse_config(params_provided)
    insarobj = read_grd_inputs(params)
    plot_synthetic_grid_los(params, insarobj, disp_points, disp_points_color)
    return


if __name__ == "__main__":
    my_params = parse_config({})  # could consider making this a command line application
    do_synthetic_grid_LOS(my_params)
