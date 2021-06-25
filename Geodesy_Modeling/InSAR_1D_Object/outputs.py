# Definition of InSAR-format data
# Write and output functions for the InSAR data format

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


def write_insar_invertible_format(InSAR_obj, unc_min, filename):
    """
    Writes InSAR displacements into insar file that can be inverted.
    Writes one header line and multiple data lines.
    InSAR_obj is in mm, and written out is in meters
    """
    print("Writing InSAR displacements into file %s " % filename);
    ofile = open(filename, 'w');
    ofile.write("# InSAR Displacements: Lon, Lat, disp(m), sigma, unitE, unitN, unitN \n");
    for i in range(len(InSAR_obj.lon)):
        if np.isnan(InSAR_obj.LOS[i]):
            continue;
        else:
            if InSAR_obj.LOS_unc:
                std = InSAR_obj.LOS_unc[i] * 0.001;  # in m
                if std < unc_min:
                    std = unc_min;
            else:   # sometimes there's an error code in LOS_unc field
                std = unc_min;
            ofile.write('%f %f ' % (InSAR_obj.lon[i], InSAR_obj.lat[i]));
            ofile.write('%f %f ' % (0.001 * InSAR_obj.LOS[i], std));  # writing in m
            ofile.write('%f %f %f\n' % (InSAR_obj.lkv_E[i], InSAR_obj.lkv_N[i], InSAR_obj.lkv_U[i]));
    ofile.close();
    return;


def plot_insar(InSAR_Obj, plotname, vmin=None, vmax=None):
    print("Plotting insar in file %s " % plotname);
    plt.figure(dpi=300, figsize=(8, 8));
    if vmin is not None:
        plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=18, cmap='rainbow', vmin=vmin, vmax=vmax);
    else:
        plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=18, cmap='rainbow');
    cb = plt.colorbar();
    cb.set_label('Displacement (mm)', fontsize=16);
    starttime = dt.datetime.strftime(InSAR_Obj.starttime, "%Y-%m");
    endtime = dt.datetime.strftime(InSAR_Obj.endtime, "%Y-%m");
    plt.title("InSAR with %d Points from %s to %s" % (len(InSAR_Obj.LOS), starttime, endtime));
    plt.savefig(plotname);
    return;
