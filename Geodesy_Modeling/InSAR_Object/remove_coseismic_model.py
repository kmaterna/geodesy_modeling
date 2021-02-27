# Remove coseismic model relative to a reference pixel. 
# Most of the subtraction math is here. 
# June 2020

import numpy as np


def remove_model_los(los_file, model_disps_file, adjusted_file):
    # Pseudocode:
    # What is the modeled delta-e, delta-n, and delta-u of each point relative to the reference point?
    # Then, how does that project into the local LOS?
    # Finally, subtract that offset from each point
    # Then we write everything except the reference line
    # Into a file with "_updated" in its name
    # Assumes the reference pixel is the last row of the data file.
    [_, _, u_pred, v_pred, w_pred] = np.loadtxt(model_disps_file, unpack=True, skiprows=1);
    [lon_meas, lat_meas, disp, sig, unit_e, unit_n, unit_u] = np.loadtxt(los_file, unpack=True, skiprows=1);
    corrected_los = [];
    los_reference_pixel = [u_pred[-1], v_pred[-1], w_pred[-1]];  # last model is reference
    for i in range(len(lon_meas)):
        model_prediction = [u_pred[i], v_pred[i], w_pred[i]];
        model_deltas = np.subtract(model_prediction, los_reference_pixel);

        los_unitv = [unit_e[i], unit_n[i], unit_u[i]];
        los_view = np.dot(model_deltas, los_unitv);

        corrected_los.append(disp[i] - los_view);

    ofile = open(adjusted_file, 'w');
    ofile.write("# Header: lon, lat, disp(m), sig(m), unitE, unitN, unitU from ground to satellite\n");
    for i in range(len(lon_meas)):
        ofile.write("%f %f %f %f %f %f %f \n" % (
            lon_meas[i], lat_meas[i], corrected_los[i], sig[i], unit_e[i], unit_n[i], unit_u[i]));
    ofile.close();
    return;


def remove_model_gps(gps_file, model_disps_file):
    # Take the prediction at each station,
    # subtract the reference prediction
    # and then adjust by that amount.
    # Assumes the reference pixel is the first row of the data file.
    [_, _, u_pred, v_pred, w_pred] = np.loadtxt(model_disps_file, unpack=True, skiprows=1);
    [lon_meas, lat_meas, u_meas, v_meas, w_meas, sige_meas, sign_meas, sigu_meas] = np.loadtxt(gps_file, unpack=True,
                                                                                               skiprows=1);
    disp_reference_pixel = [u_pred[0], v_pred[0], w_pred[0]];
    namestem = gps_file.split(".txt")[0];
    ofile = open(namestem + "_cos_corrected.txt", 'w');
    ofile.write("# Header: lon, lat, dE, dN, dU, Se, Sn, Su (m)\n");
    for i in range(len(lon_meas)):
        disp_at_station = [u_meas[i], v_meas[i], w_meas[i]];
        model_for_station = [u_pred[i], v_pred[i], w_pred[i]];
        model_delta = np.subtract(model_for_station, disp_reference_pixel);
        new_data = np.subtract(disp_at_station, model_delta);
        ofile.write("%f %f %f %f %f %f %f %f\n" % (
            lon_meas[i], lat_meas[i], new_data[0], new_data[1], new_data[2], sige_meas[i], sign_meas[i], sigu_meas[i]));
    ofile.close();
    return;
