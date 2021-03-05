import datetime as dt
import matplotlib
import numpy as np
import collections
import stacking_utilities
from .. import multiSAR_utilities
from Tectonic_Utils.read_write import netcdf_read_write
from matplotlib import pyplot as plt, cm as cm

# Collections
GrdTSData = collections.namedtuple("GrdTSData", ["dtarray", "lon", "lat", "TS"]);

# INPUT FUNCTIONS FOR NETCDF FORMAT
def inputs_TS_grd(filename, lonfile, latfile, day0=dt.datetime.strptime("2009-04-24", "%Y-%m-%d")):
    """
    Reads a TS file with associated lat/lon files
    The files generally are not orthorectified grids
    GRDnetcdf has tdata (days since day0), x, y, and zdata (3D cube)
    lon and lat files are 2D arrays with corresponding lon and lat for each point
    day0 is the day of the first acquisition in the time series (hard coded for a UAVSAR track default)
    """
    print("Reading TS Grid file  %s" % filename);
    [tdata, _, _, zdata] = netcdf_read_write.read_3D_netcdf(filename);
    print("tdata:", tdata);
    print("   where Day0 of this time series is %s " % dt.datetime.strftime(day0, "%Y-%m-%d"));
    [_, _, lon] = netcdf_read_write.read_any_grd(lonfile);
    [_, _, lat] = netcdf_read_write.read_any_grd(latfile);
    print("lon and lat:", np.shape(lon));
    print("zdata:", np.shape(zdata));
    zdata_correct_size = [];
    if np.shape(zdata[0])[0] == np.shape(lon)[0]+1 and np.shape(zdata[0])[1] == np.shape(lat)[1]+1:
        for i in range(len(zdata)):
            zdata_correct_size.append(zdata[i][0:-1,0:-1]);  # cutting off one pixel on each end for pixel node problem
    else:
        zdata_correct_size = zdata;
    dtarray = [];
    for i in range(len(tdata)):
        dtarray.append(day0 + dt.timedelta(days=int(tdata[i])));

    assert(np.shape(lon) == np.shape(zdata_correct_size[0])), ValueError("Lon and Data size don't match");
    assert(np.shape(lat) == np.shape(zdata_correct_size[0])), ValueError("Lat and Data size don't match");
    assert(np.shape(zdata_correct_size)[0] == len(dtarray)), ValueError("dtarray and zdata size don't match");
    myGridTS = GrdTSData(dtarray=dtarray, lon=lon, lat=lat, TS=zdata_correct_size);
    return myGridTS;


def plot_grid_TS_redblue(TS_GRD_OBJ, outfile, vmin=-50, vmax=200, aspect=1, incremental=False, gps_i=None,
                         gps_j=None, selected=None):
    """Make a nice time series plot.
    With incremental or cumulative displacement data."""
    num_rows_plots = 3;
    num_cols_plots = 4;

    # With selected, you can select certain intervals and combine others
    xdates = TS_GRD_OBJ.dtarray;
    if selected is None:
        selected = np.arange(len(xdates));
    TS_array = TS_GRD_OBJ.TS[selected, :, :];
    xdates = [xdates[i] for i in range(max(selected) + 1) if i in selected];

    f, axarr = plt.subplots(num_rows_plots, num_cols_plots, figsize=(16, 10), dpi=300);

    for i in range(0, len(xdates)):
        # Collect the data for each plot.
        if incremental:
            data = np.subtract(TS_array[i, :, :], TS_array[i - 1, :, :]);
        else:
            data = TS_array[i, :, :];
        data = data.T;  # for making a nice map with west approximately to the left.
        data = np.fliplr(data);  # for making a nice map with west approximately to the left.

        gps_j_flipped = gps_j;
        gps_i_flipped = [np.shape(data)[1] - x for x in gps_i];

        # Plotting now
        rownum, colnum = stacking_utilities.get_axarr_numbers(num_cols_plots, i);
        if i == 0 and incremental is True:
            axarr[rownum][colnum].set_visible(False);
            continue;
        else:
            axarr[rownum][colnum].imshow(data, aspect=aspect, cmap='RdYlBu_r', vmin=vmin, vmax=vmax);
            titlestr = dt.datetime.strftime(xdates[i], "%Y-%m-%d");
            axarr[rownum][colnum].plot(gps_i_flipped, gps_j_flipped, 'v', markersize=6, color='black');
            axarr[rownum][colnum].get_xaxis().set_visible(False);
            axarr[rownum][colnum].set_title(titlestr, fontsize=20);

    _ = f.add_axes([0.75, 0.35, 0.2, 0.3], visible=False);
    color_boundary_object = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    custom_cmap = cm.ScalarMappable(norm=color_boundary_object, cmap='RdYlBu_r');
    custom_cmap.set_array(np.arange(vmin, vmax));
    cb = plt.colorbar(custom_cmap, aspect=12, fraction=0.2, orientation='vertical');
    cb.set_label('Displacement (mm)', fontsize=18);
    cb.ax.tick_params(labelsize=12);

    plt.savefig(outfile);
    return;


def plot_pixel_ts(TS, dtarray, i, j, name, outdir):
    pixel_value, pixel_value2 = [], [];
    width_pixels = 80;
    for date in range(len(dtarray)):
        pixel_value.append(np.nanmean(TS[date, i - width_pixels:i + width_pixels, j - width_pixels:j + width_pixels]));
        pixel_value2.append(TS[date, i, j]);
    plt.figure(figsize=(8, 8));
    plt.plot(dtarray, pixel_value, '.--', markersize=12);
    plt.plot(dtarray, pixel_value2, '.--', color='red', markersize=12);
    plt.xlabel("Time");
    plt.ylabel("Displacement (mm)");
    plt.savefig(outdir + "/" + name + "_onepixel.png");
    return;


def total_ts_visualizing(myUAVSAR, gps_lon, gps_lat, gps_names, selected, outdir):
    """Vizualizations for grid time series format
    GPS points can be added by index. For future work. """
    plot_grid_TS_redblue(myUAVSAR, outdir + "/increments.png", vmin=-100, vmax=100, aspect=4,
                         incremental=True, gps_i=None, gps_j=None, selected=selected);
    plot_grid_TS_redblue(myUAVSAR, outdir + "/full_TS.png", vmin=-160, vmax=160, aspect=4,
                         incremental=False, gps_i=None, gps_j=None, selected=selected);
    # Comparing InSAR TS with GPS
    for i in range(len(gps_lon)):
        ipt, jpt = multiSAR_utilities.get_nearest_pixel_in_raster(myUAVSAR.lon, myUAVSAR.lat, gps_lon[i], gps_lat[i]);
        plot_pixel_ts(myUAVSAR.TS, myUAVSAR.dtarray, ipt, jpt, gps_names[i], outdir);
    return;
