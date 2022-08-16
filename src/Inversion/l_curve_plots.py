
import numpy as np
from matplotlib import pyplot as plt


def plot_1d_curve(param_values, misfit, axis_name, outfile, corner_point=None):
    """
    Make 1D plot for L-curve

    :param param_values: 1D array
    :param misfit: 1D array
    :param axis_name: string
    :param outfile: string
    :param corner_point: float, optional x-location where an annotation will be drawn
    """
    param = [1/x for x in param_values];

    def tick_function(X):  # for the upper axis labels
        V = 1 / X
        return ["%.4f" % z for z in V]

    # Linear scale
    plt.figure(figsize=(6, 6), dpi=300);
    plt.plot(param, misfit, '.', markersize=14);
    plt.xlabel("1/"+axis_name, fontsize=16);
    plt.ylabel('Misfit (mm/yr)', fontsize=16);
    ax1 = plt.gca();
    top, bottom = plt.gca().get_ylim();
    if corner_point is not None:  # draw optional annotation
        ax1.plot([1/corner_point, 1/corner_point], [bottom, top], '--r');
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax2 = plt.twiny(ax1);
    ax2.set_xlim(ax1.get_xlim())
    new_tick_locations = np.array(np.divide(1, param_values));
    # cutting off some labels
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations), rotation=70)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.set_xlabel(axis_name, fontsize=16)
    ax2.grid(True)
    plt.tight_layout()
    plt.savefig(outfile);

    # Log scale
    plt.figure(figsize=(6, 6), dpi=300);
    plt.plot(param, misfit, '.', markersize=14);
    ax1 = plt.gca()
    # ax1.set_yscale('log');
    ax1.set_xscale('log');
    top, bottom = plt.gca().get_ylim();
    if corner_point is not None:   # draw optional annotation
        ax1.plot([1/corner_point, 1/corner_point], [bottom, top], '--r');

    ax1.set_xlabel("1/"+axis_name, fontsize=16);
    ax1.set_ylabel('Misfit (mm/yr)', fontsize=16);
    # ax1.grid(True)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax2 = plt.twiny(ax1);
    ax2.set_xscale('log');
    ax2.set_xlim(ax1.get_xlim())
    new_tick_locations = np.array(np.divide(1, param_values));
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations), rotation=70)
    ax2.set_xlabel(axis_name, fontsize=16)
    ax2.tick_params(axis='x', labelsize=16)
    plt.tight_layout()
    plt.savefig(outfile.split('.')[0]+"_log.png");
    return;


def plot_2d_curve(alphas, penalties, misfits, param1_name, param2_name, outfile):
    """
    Make 2D surface plot for L-curve

    :param alphas: 1D array
    :param penalties: 1D array
    :param misfits: 2D array
    :param param1_name: string
    :param param2_name: string
    :param outfile: string
    """
    fig = plt.figure(figsize=(4, 4), dpi=200);
    ax = fig.add_subplot(111, projection='3d')
    alphas = [1/x for x in alphas];
    penalties = [1/x for x in penalties];
    misfits = [np.log(x) for x in misfits];
    ax.scatter(alphas, penalties, misfits, '.', c=misfits, cmap='viridis');
    ax.set_xlabel(param1_name);
    ax.set_ylabel(param2_name);
    ax.set_zlabel('Misfit (mm)')
    ax = chosen_axis_annotations(ax);
    ax.set_title("Misfit vs Regularization Parameters");
    plt.show();  # good for playing.
    plt.savefig(outfile);
    return;


def chosen_axis_annotations(ax):
    # Manual annotations for certain projects
    # Example: Heber base inversion: alpha 2, penalty 0.5.
    # WARNING: This is mixing code and data.
    chosen_alpha = 2;        # HEBER
    chosen_penalty = 0.5;    # HEBER
    alpha_range = ax.get_xlim();
    penalty_range = ax.get_ylim();
    [bottom, _] = ax.get_zlim();
    ax.plot([1/chosen_alpha, 1/chosen_alpha], penalty_range, [bottom, bottom], color='black');
    ax.plot(alpha_range, [1/chosen_penalty, 1/chosen_penalty], [bottom, bottom], color='black');
    ax.plot([1/chosen_alpha], [1/chosen_penalty], [bottom], marker='.', color='red', markersize=15);
    return ax;


def plot_l_curve_coordinator(params, misfits, outfile):
    """ Coordiantor function for driving l-curve plots """
    all_alphas = [x[0] for x in params];
    all_penalties = [x[1] for x in params];
    if len(set(all_alphas)) == 1:
        plot_1d_curve(all_penalties, misfits, 'Smoothing Penalty', outfile.split('.')[0] + "_smoothing.png");
    elif len(set(all_penalties)) == 1:
        plot_1d_curve(all_alphas, misfits, 'Slip Penalty, alpha', outfile.split('.')[0] + "_slip.png");
    else:
        plot_2d_curve(all_alphas, all_penalties, misfits, '1/alpha (slip)', '1/smoothing (smoothing)',
                      outfile.split('.')[0] + "_2d.png");
    return;
