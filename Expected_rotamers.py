########################################################################
# Create Dunbrack predicted values based on phi/psi of GNSRV
#
# Input:
# xvg file
# column ids for phi, psi and chi values
# list of rotamer values corresponding to dunbrack rotamer centers
# rotamer bin widths
# number of unique rotamers
#
# Output:
# Columns
#   1: phi bin
#   2: psi bin
#   3: # of frames in rotamer 1
#   4: # of frames in rotamer 2, etc
########################################################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math


folder = 's1s2/cluster2/'
cluster_number = '2'

# Load data
data_CP = np.loadtxt(folder + 'GNSRV_both_cluster' + cluster_number + '.xvg')

# Loop over all AA locations
for res_loop in range(0, 4):

    if (res_loop == 0):
        CP_chi1 = data_CP[:, 12].copy()
        CP_phi = data_CP[:, 4].copy()
        CP_psi = data_CP[:, 5].copy()
        data = np.loadtxt('../Sidechain_bbDependent/Dunbrack_smoothed_5/asn.bbdep.rotamers.lib', skiprows=28,
                          usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9))

    elif (res_loop == 1):
        CP_chi1 = data_CP[:, 14].copy()
        CP_phi = data_CP[:, 6].copy()
        CP_psi = data_CP[:, 7].copy()
        data = np.loadtxt('../Sidechain_bbDependent/Dunbrack_smoothed_5/ser.bbdep.rotamers.lib', skiprows=28,
                          usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9))

    elif (res_loop == 2):
        CP_chi1 = data_CP[:, 15].copy()
        CP_phi = data_CP[:, 8].copy()
        CP_psi = data_CP[:, 9].copy()
        data = np.loadtxt('../Sidechain_bbDependent/Dunbrack_smoothed_5/arg.bbdep.rotamers.lib', skiprows=28,
                          usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9))

    elif (res_loop == 3):
        CP_chi1 = data_CP[:, 19].copy()
        CP_phi = data_CP[:, 10].copy()
        CP_psi = data_CP[:, 11].copy()
        data = np.loadtxt('../Sidechain_bbDependent/Dunbrack_smoothed_5/val.bbdep.rotamers.lib', skiprows=28,
                          usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9))

    sidechain_bin_width = 30
    Nrot = 3
    # Columns of Dunbrack input
    # 0: phi
    # 1: psi
    # 2: count
    # 3: rotamer number
    # 4: probability
    # 5: chi1 Value

    # Creat phi / psi probability distribution for CP data

    CP_phi_psi = np.zeros([36, 36])
    counts = np.zeros([36, 36])
    for i in range(0, CP_phi.shape[0]):
        phi = CP_phi[i]
        phi_index = (math.floor(phi / 10) + 17)
        psi = CP_psi[i]
        psi_index = (math.floor(psi / 10) + 17)
        CP_phi_psi[psi_index, phi_index] = CP_phi_psi[psi_index, phi_index] + 1

    # Use CP data to get predicted Dunbrack distrbutions
    # Loop over all rotamers
    p_each_rot = np.zeros([Nrot, 1])
    for rotamer_loop in range(0, Nrot):

        # Get Dunbrack data for this rotamer
        ind0 = data[:, 3] == (rotamer_loop + 1)
        this_rot_Dun = data[ind0, :].copy()

        expected_values = np.zeros([37 * 37, 3])    # Expected phi/psi data
        to_plot = np.zeros([36, 36])
        for bb_loop in range(0, (this_rot_Dun.shape[0])):

            # Find phi/psi indexes based on Dunbrack data
            phi = this_rot_Dun[bb_loop, 0]
            phi_index = math.floor(phi / 10) + 17
            psi = this_rot_Dun[bb_loop, 1]
            psi_index = math.floor(psi / 10) + 17
            CP_count = CP_phi_psi[psi_index, phi_index]  # Frequency of this phi/psi in CP MD

            prob = this_rot_Dun[bb_loop, 7]  # Probability from Dunbrack data

            if (phi < 180 and psi < 180):
                #  Multiply phi/psi frequency in MD by Dunbrack Probablity
                to_plot[psi_index, phi_index] = to_plot[psi_index, phi_index] + prob * CP_count
                if (prob * CP_count > 0):
                    p_each_rot[rotamer_loop] = p_each_rot[rotamer_loop] + prob * CP_count

        to_plot = to_plot / CP_psi.shape[0]  # Normalize by number of time steps in MD

        ind0 = to_plot <= 4.0E-7
        to_plot[ind0] = 'nan'

        # Plot data
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        plt.imshow(to_plot / 100, origin='lower', cmap='jet')
        plt.colorbar()
        plt.clim(0, 0.0001)
        plt.xlabel('$\phi$', fontsize=20)
        plt.ylabel('$\psi$', fontsize=20)

        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)

        if (res_loop == 0):
            save_start = folder + 'Asn_2'
            plt.title('Asn 2', fontsize=20)
        elif (res_loop == 1):
            plt.title('Ser 3', fontsize=20)
            save_start = folder + 'Ser_3'
        elif (res_loop == 2):
            plt.title('Arg 4', fontsize=20)
            save_start = folder + 'Arg_4'
        elif (res_loop == 3):
            plt.title('Val 5', fontsize=20)
            save_start = folder + 'Val 5'

        if (rotamer_loop == 0):
            plt.title('$\chi_1$ = 60', fontsize=20)
            plt.savefig(save_start + '_60_PChi_Dun.png', bbox_inches='tight')
        if (rotamer_loop == 1):
            plt.title('$\chi_1$ = 180', fontsize=20)
            plt.savefig(save_start + '_180_PChi_Dun.png', bbox_inches='tight')
        if (rotamer_loop == 2):
            plt.title('$\chi_1$ = 300', fontsize=20)
            plt.savefig(save_start + '_300_PChi_Dun.png', bbox_inches='tight')
        plt.close()
    print(p_each_rot / CP_psi.shape[0])

    # Additional analysis of Asn
    # 12 chi 2 bins
    # if (res_loop == 0):
    #     p_each_rot = np.zeros([12 * 3, 3])
    #     for chi1_loop in range(0, 3):
    #         ind0 = data[:, 3] == (chi1_loop + 1)
    #         this_rot_Dun_c1 = data[ind0, :].copy()

    #         for chi2_loop in range(0, 12):
    #             ind0 = this_rot_Dun_c1[:, 4] == (chi2_loop + 1)
    #             this_rot_Dun = this_rot_Dun_c1[ind0, :]

    #             expected_values = np.zeros([37 * 37, 3])
    #             to_plot = np.zeros([36, 36])
    #             p_each_rot[(chi1_loop * 12) + chi2_loop, 0] = chi1_loop + 1
    #             p_each_rot[(chi1_loop * 12) + chi2_loop, 1] = chi2_loop + 1
    #             for bb_loop in range(0, (this_rot_Dun.shape[0])):

    #                 phi = this_rot_Dun[bb_loop, 0]
    #                 phi_index = math.floor(phi / 10) + 17
    #                 psi = this_rot_Dun[bb_loop, 1]
    #                 psi_index = math.floor(psi / 10) + 17
    #                 CP_count = CP_phi_psi[psi_index, phi_index]

    #                 prob = this_rot_Dun[bb_loop, 7]

    #                 if (phi < 180 and psi < 180):

    #                     to_plot[psi_index, phi_index] = to_plot[psi_index, phi_index] + prob * CP_count
    #                     if (prob * CP_count > 0):
    #                         p_each_rot[(chi1_loop * 12) + chi2_loop, 2] = p_each_rot[(chi1_loop * 12) + chi2_loop, 2] + prob * CP_count

    #             to_plot = to_plot / CP_psi.shape[0]
    #     p_each_rot[:, 2] = p_each_rot[:, 2] / CP_psi.shape[0] * 100
