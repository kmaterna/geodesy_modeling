#!/usr/bin/env python

from geodesy_modeling.Inversion import l_curve

l_curve.glob_and_drive_1d_lcurve(target_dir='.', name_of_printed_config="configs_used.txt", paramname='tikhonov0',
                                 name_of_results_file="metrics.txt", misfitname="Avg misfit",
                                 outname="smoothing_curve.png", xlabel="Smoothing", corner_point=30)
