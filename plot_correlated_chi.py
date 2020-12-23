######################################################
# plot_correlated_chi.py
#
# Plots chi1 correlations between amino acids in GNSRV
# Currently uses combined cluster data
######################################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math

# Loop over all clusters of interest
for cluster_loop in range(0, 3):
    folder = 's1s2/'

    # Load data
    if (cluster_loop == 0):
        folder = folder + 'cluster1/'
        data = np.loadtxt(folder + 'GNSRV_both_cluster1.xvg')
    if (cluster_loop == 1):
        folder = folder + 'cluster2/'
        data = np.loadtxt(folder + 'GNSRV_both_cluster2.xvg')
    if (cluster_loop == 2):
        folder = folder + 'cluster3/'
        data = np.loadtxt(folder + 'GNSRV_both_cluster3.xvg')
    if (cluster_loop == 3):
        folder = folder + 'cluster4/'
        data = np.loadtxt(folder + 'sGNSRV_both_cluster4.xvg')
    if (cluster_loop == 4):
        folder = folder + 'cluster5/'
        data = np.loadtxt(folder + 'GNSRV_both_cluster5.xvg')

    # Loop over all possible amino acid pairs (excluding Gly)
    for aa_pairs in range(0, 3):

        if (aa_pairs == 0):
            # Plot Asn, Ser
            res1_chi1 = data[:, 12].copy()
            res2_chi1 = data[:, 14].copy()
            xlabel = 'Asn $\chi_1$'
            ylabel = 'Ser $\chi_1$'
            save_str = 'Asn_Ser'
        if (aa_pairs == 1):
            # Plot Ser, Arg
            res1_chi1 = data[:, 14].copy()
            res2_chi1 = data[:, 15].copy()
            xlabel = 'Ser $\chi_1$'
            ylabel = 'Arg $\chi_1$'
            save_str = 'Ser_Arg'

        if (aa_pairs == 2):
            # Plot Arg, Val
            res1_chi1 = data[:, 15].copy()
            res2_chi1 = data[:, 19].copy()
            xlabel = 'Arg $\chi_1$'
            ylabel = 'Val $\chi_1$'
            save_str = 'Arg_Val'

        # Make all chi > 0
        ind0 = res1_chi1 < 0
        res1_chi1[ind0] = res1_chi1[ind0] + 360
        ind0 = res2_chi1 < 0
        res2_chi1[ind0] = res2_chi1[ind0] + 360

        # Make 2d data
        square_chi1 = np.zeros([36, 36])
        r1_P = np.zeros([36, 1])
        r2_P = np.zeros([36, 1])
        for i in range(0, res1_chi1.shape[0]):
            chi_X = res1_chi1[i]
            chi_Z = res2_chi1[i]
            chi_X_index = (math.floor(chi_X / 10))
            chi_Z_index = (math.floor(chi_Z / 10))
            square_chi1[chi_Z_index, chi_X_index] = square_chi1[chi_Z_index, chi_X_index] + 1

            r1_P[chi_X_index] = r1_P[chi_X_index] + 1
            r2_P[chi_Z_index] = r2_P[chi_Z_index] + 1

        square_chi1 = square_chi1 / res1_chi1.shape[0] / 100

        # Calculate amount in each bin and save to file
        percent_bins = np.zeros([9, 1])
        percent_bins[0] = sum(sum(square_chi1[0:12, 0:12])) * 100 * 100
        percent_bins[1] = sum(sum(square_chi1[0:12, 12:24])) * 100 * 100
        percent_bins[2] = sum(sum(square_chi1[0:12, 24:36])) * 100 * 100
        percent_bins[3] = sum(sum(square_chi1[12:24, 0:12])) * 100 * 100
        percent_bins[4] = sum(sum(square_chi1[12:24, 12:24])) * 100 * 100
        percent_bins[5] = sum(sum(square_chi1[12:24, 24:36])) * 100 * 100
        percent_bins[6] = sum(sum(square_chi1[24:36, 0:12])) * 100 * 100
        percent_bins[7] = sum(sum(square_chi1[24:36, 12:24])) * 100 * 100
        percent_bins[8] = sum(sum(square_chi1[24:36, 24:36])) * 100 * 100
        np.savetxt(folder + save_str + '_chi1_vs_chi1.txt', percent_bins)

        # Set tiny values to show up as white
        ind0 = square_chi1 <= 1E-8
        square_chi1[ind0] = 'nan'

        # Make plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        plt.imshow(square_chi1, origin='lower', cmap='jet')
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel(ylabel, fontsize=20)
        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        plt.clim(0, 0.0005)
        plt.axis([-.5, 35.5, -.5, 35.5])
        plt.colorbar()
        plt.savefig(folder + save_str + '_chi1_vs_chi1.png', bbox_inches='tight')
        plt.close()

        # Make uncorrelated Plot
        r1_P = r1_P / res1_chi1.shape[0] / 10
        r2_P = r2_P / res1_chi1.shape[0] / 10
        square_chi1 = np.zeros([36, 36])
        for i in range(0, 36):
            for j in range(0, 36):
                square_chi1[i, j] = r2_P[i] * r1_P[j]

        # Calculate amount in each bin and save to file
        percent_bins = np.zeros([9, 1])
        percent_bins[0] = sum(sum(square_chi1[0:12, 0:12])) * 100 * 100
        percent_bins[1] = sum(sum(square_chi1[0:12, 12:24])) * 100 * 100
        percent_bins[2] = sum(sum(square_chi1[0:12, 24:36])) * 100 * 100
        percent_bins[3] = sum(sum(square_chi1[12:24, 0:12])) * 100 * 100
        percent_bins[4] = sum(sum(square_chi1[12:24, 12:24])) * 100 * 100
        percent_bins[5] = sum(sum(square_chi1[12:24, 24:36])) * 100 * 100
        percent_bins[6] = sum(sum(square_chi1[24:36, 0:12])) * 100 * 100
        percent_bins[7] = sum(sum(square_chi1[24:36, 12:24])) * 100 * 100
        percent_bins[8] = sum(sum(square_chi1[24:36, 24:36])) * 100 * 100
        np.savetxt(folder + save_str + '_chi1_uncoorelated.txt', percent_bins)

        # Set tiny values to show up as white
        ind0 = square_chi1 <= 1E-8
        square_chi1[ind0] = 'nan'

        # Make plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        plt.imshow(square_chi1, origin='lower', cmap='jet')
        plt.xlabel(xlabel, fontsize=20)
        plt.ylabel(ylabel, fontsize=20)
        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
        plt.clim(0, 0.0005)
        plt.axis([-.5, 35.5, -.5, 35.5])
        plt.colorbar()
        plt.savefig(folder + save_str + '_chi1_uncoorelated.png', bbox_inches='tight')
        plt.close()