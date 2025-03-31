
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def heatmap_plot(df, plotname, plottitle=''):
    """
    Draw a heat map of a stress component's positive and negative values over all geometries.

    :param df: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb']
    :param plotname: string to name the figure
    :param plottitle: string to put as title of the figure
    :return:
    """
    # Step 2: Make plot of shear stress from this stress tensor
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.05, hspace=0.1)  # Adjust spacing
    fig.suptitle(plottitle, fontsize=18, y=0.92)
    dip = -10
    levels = np.arange(-8.6, 8.6, 0.1)
    for i in range(5):
        for j in range(2):
            dip += 10
            ax = fig.add_subplot(gs[i, j])

            # Create tricontour plots
            cf = ax.tricontourf(df[df['dip'] == dip]['strike'], df[df['dip'] == dip]['rake'],
                                df[df['dip'] == dip]['shear'], levels=levels, cmap='magma')

            ax.text(5, 120, f"Dip={dip}")
            if j == 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Rake')
            if i != 4:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Strike')

    # Add a single colorbar on the right
    cbar_ax = fig.add_axes([0.92, 0.2, 0.02, 0.6])  # [left, bottom, width, height]
    fig.colorbar(cf, cax=cbar_ax, label="$\mu = 0.1$ Shear Stress (kPa)")
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
    return
