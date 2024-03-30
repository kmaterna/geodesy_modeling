from elastic_stresses_py.PyCoulomb import fault_slip_object as library
import os


def visualize_GF_elements(GF_elements_list, outdir, exclude_list=()):
    """
    Aside to the main calculation, just view each GF.

    :param GF_elements_list: list of GF_elements objects
    :param outdir: string for outdir
    :param exclude_list: optional list of GfElement.param_name to exclude from visualizing
    """
    if exclude_list == 'all':
        return
    for GF_el in GF_elements_list:
        if GF_el.param_name in exclude_list:   # plot these elements separately, like individual CSZ patches
            continue
        print(GF_el.param_name)
        if GF_el.param_name == "CSZ":
            scale_arrow = (1.0, 0.010, "1 cm")
        else:
            scale_arrow = (1.0, 0.001, "1 mm")
        library.plot_fault_slip.map_source_slip_distribution(GF_el.fault_dict_list, os.path.join(outdir, "gf_" +
                                                             GF_el.param_name + "_only.png"),
                                                             disp_points=GF_el.disp_points,
                                                             region=[-127, -119.7, 37.7, 43.3],
                                                             scale_arrow=scale_arrow,
                                                             v_labeling_interval=0.001)
    return
