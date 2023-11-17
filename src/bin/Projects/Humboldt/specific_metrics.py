# Metrics and print statements specific to the Humboldt modeling project


import Tectonic_Utils.seismo.moment_calculations as moment_calcs
import Elastic_stresses_py.PyCoulomb.fault_slip_object as library


def write_custom_humboldt_metrics(outfile, values, GF_elements):
    """
    Accounting of moment rate on the CSZ and other faults, specific to Humboldt project

    :param outfile: existing file name
    :param values: model vector
    :param GF_elements: list of green's functions elements associated with vector of model values
    """

    with open(outfile, 'a') as ofile:
        csz_fault_dicts = [];
        for value, gf_element in zip(values, GF_elements):
            if gf_element.param_name == 'CSZ_dist':
                if gf_element.fault_dict_list[0].lat > 43:  # considering southern section only
                    continue;
                new_patches = library.fault_slip_object.change_fault_slip_list(gf_element.fault_dict_list, value/100);
                csz_fault_dicts.append(new_patches[0]);
        csz_mo = library.fault_slip_object.get_total_moment(csz_fault_dicts);
        ofile.write("\nCSZ Over 300 years: equivalent to Mw");
        ofile.write("%f \n" % moment_calcs.mw_from_moment(csz_mo*300));
        ofile.write("%f N-m\n" % csz_mo);

        new_patches = [];
        for value, gf_element in zip(values, GF_elements):
            if gf_element.param_name == 'LSFRev':
                new_patches = library.fault_slip_object.change_fault_slip_list(gf_element.fault_dict_list, value / 100);
        lsf_mo = library.fault_slip_object.get_total_moment(new_patches);
        ofile.write("\nLSFRev Over 300 years: equivalent to Mw");
        ofile.write("%f \n" % moment_calcs.mw_from_moment(lsf_mo*300));
        ofile.write("%f N-m\n" % lsf_mo);

        for value, gf_element in zip(values, GF_elements):
            if gf_element.param_name == 'MadRiverRev':
                new_patches = library.fault_slip_object.change_fault_slip_list(gf_element.fault_dict_list, value / 100);
        lsf_mo = library.fault_slip_object.get_total_moment(new_patches);
        ofile.write("\nMadRiverRev Over 300 years: equivalent to Mw");
        ofile.write("%f \n" % moment_calcs.mw_from_moment(lsf_mo*300));
        ofile.write("%f N-m\n" % lsf_mo);
    return;


def write_custom_misfit_metrics(outfile, rms_list):
    """
    The "rms_list" field is a bit specific to the experiment.

    :param outfile: filename of existing file, string.
    :param rms_list: array of floats, mm/yr and normalized.  This might be a little project-specific
    """
    with open(outfile, 'a') as ofile:
        report_string = "RMS misfit [h, v, t]: %f %f %f mm/yr\n" % (rms_list[0], rms_list[1], rms_list[2]);
        ofile.write(report_string);
        report_string = "RMS normalized [h, v, t]: %f %f %f \n" % (rms_list[3], rms_list[4], rms_list[5]);
        ofile.write(report_string);
    return;
