"""
Write and output functions for 1D InSAR data format
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import zipfile
import tempfile
import datetime as dt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from tectonic_utils.geodesy import insar_vector_functions
from geodesy_modeling import general_utils


def write_insar_invertible_format(InSAR_obj, filename, unc_min=0):
    """
    Write InSAR displacements into insar file that can be inverted.
    Write one header line and multiple data lines.
    InSAR_obj is in mm, and written out is in meters
    """
    print("Writing InSAR displacements into file %s " % filename)
    ofile = open(filename, 'w')
    ofile.write("# InSAR Displacements: Lon, Lat, disp(m), sigma, unitE, unitN, unitU \n")
    for i in range(len(InSAR_obj.lon)):
        if np.isnan(InSAR_obj.LOS[i]):
            continue
        else:
            if InSAR_obj.LOS_unc is not None:
                std = InSAR_obj.LOS_unc[i] * 0.001  # in m
                if std < unc_min:
                    std = unc_min
            else:   # sometimes there's an error code in LOS_unc field
                std = unc_min
            ofile.write('%f %f ' % (InSAR_obj.lon[i], InSAR_obj.lat[i]))
            ofile.write('%f %f ' % (0.001 * InSAR_obj.LOS[i], std))  # writing in m
            ofile.write('%f %f %f\n' % (InSAR_obj.lkv_E[i], InSAR_obj.lkv_N[i], InSAR_obj.lkv_U[i]))
    ofile.close()
    return


def plot_insar(InSAR_Obj, plotname, vmin=None, vmax=None, lons_annot=(), lats_annot=(), refpix=None,
               title=None, colormap='viridis', symbolsize=28):
    """lons_annot and lat_annot are for lines to annotate the plot, such as faults or field boundaries"""
    print("Plotting insar in file %s " % plotname)
    fig = plt.figure(dpi=300, figsize=(5, 5))
    if vmin is not None:
        im = plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=symbolsize,
                         cmap=colormap, vmin=vmin, vmax=vmax)
    else:
        im = plt.scatter(InSAR_Obj.lon, InSAR_Obj.lat, c=InSAR_Obj.LOS, s=symbolsize, cmap=colormap)
    if len(lons_annot) > 0:
        plt.plot(lons_annot, lats_annot, color='darkred')
    if refpix:  # expect a list or tuple, (lon, lat)
        plt.plot(refpix[0], refpix[1], '.', color='red')
    if title:
        plt.title(title)
    else:
        if InSAR_Obj.starttime is None:
            starttime, endtime = ' ', ' '
        else:
            starttime = dt.datetime.strftime(InSAR_Obj.starttime, "%Y-%m")
            endtime = dt.datetime.strftime(InSAR_Obj.endtime, "%Y-%m")
        plt.title("InSAR with %d Points from %s to %s" % (len(InSAR_Obj.LOS), starttime, endtime))
    cb = fig.colorbar(im, ax=plt.gca())
    cb.set_label('Displacement (mm)', fontsize=16)
    plt.savefig(plotname)
    return


def plot_look_vectors(Data, plotname):
    """
    Visualize the look vectors for gut-checking.
    """
    print("Plotting look vectors in %s " % plotname)
    f, axarr = plt.subplots(1, 3, figsize=(11, 4), dpi=300)
    title_fontsize = 14
    flight_heading, _ = insar_vector_functions.look_vector2flight_incidence_angles(Data.lkv_E[0], Data.lkv_N[0],
                                                                                   Data.lkv_U[0])
    x_flight, y_flight, x_los, y_los = general_utils.get_los_and_flight_vectors(flight_heading, 'right')

    cm = plt.cm.get_cmap('viridis')

    im1 = axarr[0].scatter(Data.lon, Data.lat, c=Data.lkv_E, cmap=cm)
    cbaxes = inset_axes(axarr[0], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[0].bbox)
    plt.colorbar(im1, cax=cbaxes, orientation='vertical')
    axarr[0].set_title("Look Vector East (ground to sat)", fontsize=title_fontsize)
    axarr[0].quiver(0.1, 0.85, 2*x_flight, 2*y_flight, transform=axarr[0].transAxes, scale=8, scale_units='inches')
    axarr[0].quiver(0.1, 0.85, x_los, y_los, transform=axarr[0].transAxes, scale=8, scale_units='inches')

    im2 = axarr[1].scatter(Data.lon, Data.lat, c=Data.lkv_N, cmap=cm)
    axarr[1].set_yticks([])
    cbaxes = inset_axes(axarr[1], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[1].bbox)
    plt.colorbar(im2, cax=cbaxes, orientation='vertical')
    axarr[1].set_title("Look Vector North", fontsize=title_fontsize)

    im3 = axarr[2].scatter(Data.lon, Data.lat, c=Data.lkv_U, cmap=cm)
    axarr[2].set_yticks([])
    cbaxes = inset_axes(axarr[2], width="5%", height="90%", loc=1, bbox_to_anchor=axarr[2].bbox)
    plt.colorbar(im3, cax=cbaxes, orientation='vertical')
    axarr[2].set_title("Look Vector Up", fontsize=title_fontsize)

    plt.savefig(plotname)
    return


def write_kml_points(InSAR_obj, outfile, cmap='RdBu_r', vmin=None, vmax=None, n_bins=11,
                     alpha=1.0, icon_href='http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png',
                     icon_scale=0.9, include_name=False, include_extended=True, kmz=False, unit='mm'):
    """
    Write a KML (or KMZ) with point Placemarks colored by LOS displacement.

    Parameters
    ----------
    InSAR_obj : Insar1dObject
        1D InSAR object containing lon, lat, LOS and optional coherence.
    outfile : str
        Path to output .kml or .kmz file.
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : float or None
        Color range limits; if None, uses robust percentiles (2, 98).
    n_bins : int
        Number of discrete color bins/styles.
    alpha : float
        Overall alpha in [0,1] applied to styles.
    icon_href : str
        Icon URL for point styling.
    icon_scale : float
        Icon scale for point styling.
    include_name : bool
        If True, include the value of the point as a text attribute in the KML.
    include_extended : bool
        If True, include ExtendedData with LOS (and coherence if available).
    kmz : bool
        If True, write a KMZ; otherwise write plain KML.
    unit : str
        Unit string to annotate LOS values (e.g., 'mm').
    """
    print("Writing InSAR displacements into file %s " % outfile)

    los = np.asarray(InSAR_obj.LOS, dtype=float)
    lon = np.asarray(InSAR_obj.lon, dtype=float)
    lat = np.asarray(InSAR_obj.lat, dtype=float)

    valid = ~np.isnan(los)
    los_valid = los[valid]
    lon_valid = lon[valid]
    lat_valid = lat[valid]

    if vmin is None or vmax is None:
        p2, p98 = np.nanpercentile(los_valid, [2, 98])
        if vmin is None:
            vmin = p2
        if vmax is None:
            vmax = p98
    if vmin == vmax:
        vmin -= 1.0
        vmax += 1.0

    cmap_obj = cm.get_cmap(cmap, n_bins)
    edges = np.linspace(vmin, vmax, n_bins + 1)
    idx = np.digitize(los_valid, edges, right=False) - 1
    idx[idx < 0] = 0
    idx[idx >= n_bins] = n_bins - 1

    def rgba_to_kml_color(r, g, b, a):
        aa = int(round(a * alpha * 255))
        rr = int(round(r * 255))
        gg = int(round(g * 255))
        bb = int(round(b * 255))
        return f"{aa:02x}{bb:02x}{gg:02x}{rr:02x}"

    styles = []
    for i in range(n_bins):
        rgba = cmap_obj(i / max(1, n_bins - 1))
        styles.append(rgba_to_kml_color(rgba[0], rgba[1], rgba[2], rgba[3]))

    name = os.path.basename(outfile)

    # Build KML content
    parts = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append('<kml xmlns="http://www.opengis.net/kml/2.2">')
    parts.append('<Document>')
    parts.append(f'<name>{name}</name>')

    # Styles
    for i, col in enumerate(styles):
        parts.append(f'<Style id="s{i}">')
        parts.append('<IconStyle>')
        parts.append(f'<color>{col}</color>')
        parts.append(f'<scale>{icon_scale}</scale>')
        parts.append('<Icon>')
        parts.append(f'<href>{icon_href}</href>')
        parts.append('</Icon>')
        parts.append('</IconStyle>')
        parts.append('</Style>')

    # Placemarks
    coh = getattr(InSAR_obj, 'coherence', None)
    for i in range(len(lon_valid)):
        style_id = idx[i]
        val = los_valid[i]
        x = lon_valid[i]
        y = lat_valid[i]
        parts.append('<Placemark>')
        if include_name:
            parts.append(f'<name>{val:.3f} {unit}</name>')
        if include_extended:
            parts.append('<ExtendedData>')
            parts.append(f'<Data name="LOS"><value>{val:.6f}</value></Data>')
            if coh is not None and not np.isnan(coh[valid][i]):
                parts.append(f'<Data name="coherence"><value>{float(coh[valid][i]):.4f}</value></Data>')
            parts.append('</ExtendedData>')
        parts.append(f'<styleUrl>#s{style_id}</styleUrl>')
        parts.append('<Point>')
        parts.append(f'<coordinates>{x:.6f},{y:.6f},0</coordinates>')
        parts.append('</Point>')
        parts.append('</Placemark>')

    parts.append('</Document>')
    parts.append('</kml>')
    kml_text = "\n".join(parts)

    if kmz:
        with tempfile.TemporaryDirectory() as td:
            doc_path = os.path.join(td, 'doc.kml')
            with open(doc_path, 'w') as f:
                f.write(kml_text)
            with zipfile.ZipFile(outfile, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
                zf.write(doc_path, arcname='doc.kml')
    else:
        with open(outfile, 'w') as f:
            f.write(kml_text)
    return
