
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def heatmap_plot(df, plotname, plottitle='', stresstype='shear'):
    """
    Draw a heat map of a stress component's positive and negative values over all geometries.

    :param df: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb']
    :param plotname: string to name the figure
    :param plottitle: string to put as title of the figure
    :param stresstype: string, 'shear', 'normal', or 'coulomb'
    :return:
    """
    # Step 2: Make plot of stress component from this stress tensor
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.05, hspace=0.1)  # Adjust spacing
    fig.suptitle(plottitle, fontsize=18, y=0.92)
    dip = -10
    levels = np.arange(-18.6, 18.6, 0.1)
    cf = []
    for i in range(5):
        for j in range(2):
            dip += 10
            ax = fig.add_subplot(gs[i, j])

            # Create tricontour plots
            cf = ax.tricontourf(df[df['dip'] == dip]['strike'], df[df['dip'] == dip]['rake'],
                                df[df['dip'] == dip][stresstype], levels=levels, cmap='magma')

            ax.text(5, 120, f"Dip={dip}", backgroundcolor='0.75')
            if j == 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Rake')
            if i != 4:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Strike')

    # Add a single colorbar on the right
    cbar_ax = fig.add_axes((0.92, 0.2, 0.02, 0.6))  # [left, bottom, width, height]
    fig.colorbar(cf, cax=cbar_ax, label="$\mu = 0.1$ "+stresstype+" stress (kPa)")
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
    plt.close()
    return


def two_heatmaps(df1, df2, plotname, plottitle='', stresstype='shear', nextype='shear_nex'):
    """
    Draw a heat map of a stress component's positive and negative values over all geometries.

    :param df1: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear_nex', 'coulomb_nex']
    :param df2: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb']
    :param plotname: string to name the figure
    :param plottitle: string to put as title of the figure
    :param stresstype: string, 'shear', 'normal', or 'coulomb'
    :param nextype: string, 'shear_nex' or 'coulomb_nex'
    :return:
    """
    # Step 2: Make plot of stress component from this stress tensor
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.05, hspace=0.1)  # Adjust spacing
    fig.suptitle(plottitle, fontsize=18, y=0.92)
    dip = -10
    levels = np.arange(0, 1.0, 0.01)
    cf = []
    for i in range(5):
        for j in range(2):
            dip += 10
            ax = fig.add_subplot(gs[i, j])

            x_values1 = df1[df1['dip'] == dip]['strike']
            y_values1 = df1[df1['dip'] == dip]['rake']
            z_values1 = df1[df1['dip'] == dip][nextype]
            z_values2 = df2[df2['dip'] == dip][stresstype]

            z_values1[z_values1 < 0] = 0
            z_values2[z_values2 > 0] = 0
            z_values_total = np.multiply(z_values1, z_values2)  # Product positive when both conditions are met

            # Create tricontour plots
            cf = ax.tricontourf(x_values1, y_values1, z_values_total, levels=levels, cmap='magma')

            ax.text(5, 120, f"Dip={dip}", backgroundcolor='0.75')
            if j == 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Rake')
            if i != 4:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Strike')

    # Add a single colorbar on the right
    cbar_ax = fig.add_axes((0.92, 0.2, 0.02, 0.6))  # [left, bottom, width, height]
    fig.colorbar(cf, cax=cbar_ax, label="$\mu = 0.1$ " + stresstype + " stress (kPa)")
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
    plt.close()
    return


def three_heatmaps(df1, df2, df3, plotname, plottitle='', stresstype='shear', nextype='shear_nex'):
    """
    Draw a heat map of a stress component's positive and negative values over all geometries.

    :param df1: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear_nex', 'coulomb_nex']
    :param df2: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb']
    :param df3: a pandas dataframe that contains columns=['strike', 'dip', 'rake', 'shear', 'normal', 'coulomb']
    :param plotname: string to name the figure
    :param plottitle: string to put as title of the figure
    :param stresstype: string, 'shear', 'normal', or 'coulomb'
    :param nextype: string, 'shear_nex' or 'coulomb_nex'
    :return:
    """
    # Step 2: Make plot of stress component from this stress tensor
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(5, 2, figure=fig, wspace=0.05, hspace=0.1)  # Adjust spacing
    fig.suptitle(plottitle, fontsize=18, y=0.92)
    dip = -10
    levels = np.arange(0, 1.0, 0.01)
    cf = []
    for i in range(5):
        for j in range(2):
            dip += 10
            ax = fig.add_subplot(gs[i, j])

            x_values1 = df1[df1['dip'] == dip]['strike']
            y_values1 = df1[df1['dip'] == dip]['rake']
            z_values1 = df1[df1['dip'] == dip][nextype]
            z_values2 = df2[df2['dip'] == dip][stresstype]
            z_values3 = df3[df3['dip'] == dip][stresstype]

            z_values1[z_values1 < 0] = 0
            z_values2[z_values2 > -1] = 0  # coseismic shear stress change must be less than -1 kPa
            z_values3[z_values3 < 0] = 0   # interseismic stress rate must be positive
            z_values_total = np.multiply(z_values1, z_values2)  # Product positive when both conditions are met
            z_values_total = np.multiply(z_values_total, z_values3)  # Product positive when both conditions are met

            # Create tricontour plots
            cf = ax.tricontourf(x_values1, y_values1, z_values_total, levels=levels, cmap='magma')

            ax.text(5, 120, f"Dip={dip}", backgroundcolor='0.75')
            if j == 1:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel('Rake')
            if i != 4:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel('Strike')

    # Add a single colorbar on the right
    cbar_ax = fig.add_axes((0.92, 0.2, 0.02, 0.6))  # [left, bottom, width, height]
    fig.colorbar(cf, cax=cbar_ax, label="$\mu = 0.1$ " + stresstype + " stress (kPa)")
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
    plt.close()
    return
