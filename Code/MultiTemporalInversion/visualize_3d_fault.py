import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
import matplotlib
import slippy.io
import slippy.xyz2geo
import slippy.patch
import qtm_tools


def input_slippy_patches(slip_file):
    inputs = slippy.io.read_slip_data(slip_file);
    patch_pos_geo = inputs[0]  # 3-column pos
    bm = slippy.xyz2geo.create_default_basemap(patch_pos_geo[:, 0], patch_pos_geo[:, 1], resolution='i')
    pos_cart = slippy.xyz2geo.geodetic_to_cartesian(patch_pos_geo, bm)
    strike = inputs[1]
    dip = inputs[2]
    length = inputs[3]
    width = inputs[4]
    slip = inputs[5]
    patches = [slippy.patch.Patch(p, l, w, s, d) for p, l, w, s, d in zip(pos_cart, length, width, strike, dip)]
    return patches, slip, bm;


def plot_3d_patches(patches, slip, bm, qtm_filename):
    # The function that makes a 3D plot for exploring.
    polys = []
    for p in patches:
        vert = p.get_3d_polygon();
        polys += [vert];  # get the 3D vertex positions in meters

    # Create a 3D patch collection
    pc = Poly3DCollection(polys, edgecolor='black', zorder=0, alpha=0.9, cmap='viridis')
    total_slip = np.array(np.sqrt(np.add(np.square(slip[:, 0]), np.square(slip[:, 1]))));
    pc.set_array(total_slip)

    # Get some earthquakes into the right coordinate system
    catalog = qtm_tools.read_simple_catalog_txt(qtm_filename);
    eq_pos_geo = [];
    for i in range(len(catalog.lon)):
        eq_pos_geo.append([catalog.lon[i], catalog.lat[i], -1000*catalog.depth[i]]);
    eq_pos_cart = slippy.xyz2geo.geodetic_to_cartesian(eq_pos_geo, bm);

    fig = plt.figure(dpi=200);
    ax = fig.add_subplot(111, projection='3d')
    ax.add_collection(pc)
    # for i in range(len(catalog.lon)):
    #     if catalog.Mag[i]>1.0:
    #         ax.plot([eq_pos_cart[i,0]], [eq_pos_cart[i,1]], [eq_pos_cart[i,2]], marker='.', markersize=catalog.Mag[i]*catalog.Mag[i]/2, color='black');
    ax.plot(eq_pos_cart[:,0], eq_pos_cart[:,1], eq_pos_cart[:,2], marker='.', markersize=2, color='black', linewidth=0);  # every event at once.
    ax.set_xlim([-5000, 5000])
    ax.set_ylim([-5000, 5000])
    ax.set_zlim([-8000, 0])
    ax.set_xlabel('Easting (m)');
    ax.set_ylabel('Northing (m)');
    ax.set_zlabel('Depth (m)');
    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.3)
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap='viridis'), ax=ax)
    cbar.set_label('Total Slip (m)');

    plt.show()
    return;


if __name__ == "__main__":
    slip_output_file = sys.argv[1];
    qtm_filename = sys.argv[2];
    patches, slip, bm = input_slippy_patches(slip_output_file);
    plot_3d_patches(patches, slip, bm, qtm_filename);
