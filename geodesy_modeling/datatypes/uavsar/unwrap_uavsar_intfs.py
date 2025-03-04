#!/usr/bin/env python
"""
July 2020
Manipulate a uavsar file from the jpl website into a properly unwrapped and geocoded data file.
It takes .ann, .int.grd, and .cor.grd files.
This is a complicated multi-step process that needs both ISCE and GMTSAR functions.
It writes the filtered, unwrapped, masked, geocoded interferograms in isce format,
and writes los.rdr.geo in isce format too.  Useful for Kite downsampling next.
"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess
from cubbie.read_write_insar_utilities import jpl_uav_read_write, isce_read_write, netcdf_plots
from cubbie.intf_generating import isce_geocode_tools
from cubbie.math_tools import mask_and_interpolate
from Tectonic_Utils.read_write import netcdf_read_write


def read_jpl_ground_range_data(data_file_slant, corr_file_slant, ann_file):
    real, imag = jpl_uav_read_write.read_igram_data(data_file_slant, ann_file, igram_type='ground',
                                                    return_type='real_imag')
    corr = jpl_uav_read_write.read_corr_data(corr_file_slant, ann_file, igram_type='ground')
    return real, imag, corr


def cut_and_write_out_igram(real, imag, corr, cut_rowcol):
    """Cut imag, real, and corr arrays write them in lots of formats."""
    real = real[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]]
    imag = imag[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]]
    corr = corr[cut_rowcol[0]:cut_rowcol[1], cut_rowcol[2]:cut_rowcol[3]]
    phase = np.arctan2(imag, real)
    complex_numbers = np.float32(real) + 1j * np.float32(imag)
    cor32 = np.float32(corr)
    xdata = range(0, np.shape(phase)[1])
    ydata = range(0, np.shape(phase)[0])
    (ny, nx) = np.shape(phase)
    netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, 'radians', 'phase.grd', dtype=float)
    netcdf_read_write.produce_output_netcdf(xdata, ydata, corr, 'corr', 'corr.grd', dtype=float)
    netcdf_plots.produce_output_plot('phase.grd', 'Phase', 'phase.png', 'phase', aspect=1.0, invert_yaxis=False)
    netcdf_plots.produce_output_plot('corr.grd', 'Coherence', 'corr.png', 'corr', aspect=1.0, cmap='binary_r',
                                     invert_yaxis=False)
    isce_read_write.write_isce_data(complex_numbers, nx, ny, 'CFLOAT', 'cut_slc.int')
    isce_read_write.write_isce_data(cor32, nx, ny, 'FLOAT', 'cut_cor.cor')
    return


def filter_with_looks_py(after_filtering, after_filtering_corr, looks_x, looks_y):
    # # FILTER STEP: Filter using looks.py.
    subprocess.call(['looks.py', '-i', 'cut_slc.int', '-o', after_filtering, '-r', str(looks_x), '-a',
                     str(looks_y)], shell=False)
    subprocess.call(
        ['looks.py', '-i', 'cut_cor.cor', '-o', after_filtering_corr, '-r', str(looks_x), '-a', str(looks_y)],
        shell=False)
    plotting_filtering('cut_slc.int', after_filtering)
    return


def plotting_filtering(before_file, after_file):
    before_phase = isce_read_write.read_phase_data(before_file)
    after_phase = isce_read_write.read_phase_data(after_file)
    f, axarr = plt.subplots(1, 2, figsize=(10, 8), dpi=300)
    axarr[0].imshow(before_phase, cmap='rainbow')
    axarr[0].set_title('Before')
    axarr[1].imshow(after_phase, cmap='rainbow')
    axarr[1].set_title('After')
    plt.savefig("before_after_filtering.png")
    return


def multiply_igram_by_coherence_mask(after_filtering, after_filtering_corr, cutoff):
    """Multiply by coherence mask"""
    phase = isce_read_write.read_phase_data(after_filtering)
    corr = isce_read_write.read_scalar_data(after_filtering_corr)
    coherence_mask = mask_and_interpolate.make_coherence_mask(corr, cutoff)
    masked_phase = mask_and_interpolate.apply_coherence_mask(phase, coherence_mask)
    xdata = range(0, np.shape(phase)[1])
    ydata = range(0, np.shape(phase)[0])
    netcdf_read_write.produce_output_netcdf(xdata, ydata, phase, 'radians', 'phase_filtered.grd')
    netcdf_read_write.produce_output_netcdf(xdata, ydata, masked_phase, 'radians', 'phase_masked.grd', dtype=float)
    netcdf_read_write.produce_output_netcdf(xdata, ydata, corr, 'corr', 'corr.grd', dtype=float)
    netcdf_plots.produce_output_plot('phase_filtered.grd', 'Phase', 'phase_filtered.png', 'phase', aspect=1.0,
                                     invert_yaxis=False)
    netcdf_plots.produce_output_plot('phase_masked.grd', 'Phase', 'phase_masked.png', 'phase', aspect=1.0,
                                     invert_yaxis=False)
    netcdf_plots.produce_output_plot('corr.grd', 'Coherence', 'corr.png', 'corr', aspect=1.0, cmap='binary_r',
                                     invert_yaxis=False)
    return masked_phase, coherence_mask


def phase_interpolation(phase):
    """Perform phase interpolation"""
    interp_array = mask_and_interpolate.interpolate_2d(phase)
    xdata = range(0, np.shape(phase)[1])
    ydata = range(0, np.shape(phase)[0])
    netcdf_read_write.produce_output_netcdf(xdata, ydata, interp_array, 'radians', 'phase_interp.grd', dtype=float)
    netcdf_plots.produce_output_plot('phase_interp.grd', 'Phase', 'phase_interp.png', 'phase', aspect=1.0,
                                     invert_yaxis=False)
    return interp_array


def read_and_reapply_mask(mask):  # re-apply the mask
    unw_grd = netcdf_read_write.read_netcdf4("unwrap.grd")
    unw_grd = np.multiply(unw_grd, mask)
    xdata = range(0, np.shape(unw_grd)[1])
    ydata = range(0, np.shape(unw_grd)[0])
    netcdf_read_write.produce_output_netcdf(xdata, ydata, unw_grd, 'radians', 'unwrap_masked.grd')
    netcdf_plots.produce_output_plot('unwrap_masked.grd', 'Unwrapped Phase', 'unw_masked.png', 'phase', aspect=1.0,
                                     invert_yaxis=False)
    return


def main(ann_file, data_file, corr_file, after_filtering, after_filtering_corr, cut_rowcol, cor_cutoff, looks_x,
         looks_y, wavelength):
    # Example: ann_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.ann"
    # data_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.int.grd"  # 1 GB
    # corr_file = "Downloads/SanAnd_08508_11073-010_12083-007_0321d_s01_L090HH_01.cor.grd"  # 500 Mb
    # after_filtering = "cut_filtered_slc.int"
    # after_filtering_corr = "cut_filtered_cor.cor"
    # cut_rowcol = [2500, 5100, 7800, 13000]  # how to cut the frame for unwrapping and filtering.
    # cor_cutoff = 0.21
    # wavelength=237.9  # mm
    # looks_y = 5
    # looks_x = 5

    # # WE BEGIN WITH CUTTING, WRITING, AND FILTERING THE INTERFEROGRAM.
    real, imag, corr = read_jpl_ground_range_data(data_file, corr_file, ann_file)
    cut_and_write_out_igram(real, imag, corr, cut_rowcol)  # cut and write in many formats
    filter_with_looks_py(after_filtering, after_filtering_corr, looks_x,
                         looks_y)  # Filter the igrams using ISCE's toolkit

    # # MASK, INTERPOLATE, AND UNWRAP STEP
    masked_phase, mask = multiply_igram_by_coherence_mask(after_filtering, after_filtering_corr,
                                                          cor_cutoff)  # multiply by coherence mask
    phase_interpolation(masked_phase)  # perform phase interpolation
    subprocess.call(['/Users/kmaterna/Documents/B_Research/Salton/Brawley_multiSAR_project/Code/custom_unwrap.sh'],
                    shell=True)  # THEN UNWRAP
    read_and_reapply_mask(mask)  # re-apply the mask and write

    # # THEN, GEOCODE BASED ON ANNOTATION FILE, RANGES, AND MULTILOOKING PARAMETERS
    x_axis, y_axis = isce_geocode_tools.get_geocoded_axes_from_ann(ann_file, cut_rowcol, looks_x, looks_y)
    isce_geocode_tools.write_unwrapped_ground_range_displacements("unwrap_masked.grd", "uavsar.unw.geo", x_axis, y_axis,
                                                                  wavelength)
    isce_geocode_tools.create_los_rdr_geo_from_ground_ann_file(ann_file, x_axis, y_axis)  # produce los.rdr.geo
