########################################################################
# CP_MD_cluster_phi_psi_chi.py
#
# Creates all phi/psi/chi plots for each cluster and amino acid
# from a given GNSRV MD simulation
#
########################################################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import math

folder = 's1s2/'

# outlines of phi/psi regions
beta_x = np.array([-180, -180, -50, -50, -115, -115, -180]) / 10
beta_y = np.array([80, 180, 180, 80, 80, 100, 100]) / 10
beta_x1 = np.array([-180, -180, -50, -50, -180]) / 10
beta_y1 = np.array([-180, -170, -170, -180, -180]) / 10

pII_x = np.array([-110, -110, -50, -50, -110]) / 10
pII_y = np.array([120, 180, 180, 120, 120]) / 10

pII_x1 = np.array([-180, -180, -115, -115, -180]) / 10
pII_y1 = np.array([50, 100, 100, 50, 50]) / 10

alphaL_x = np.array([5, 5, 75, 75, 5]) / 10
alphaL_y = np.array([25, 120, 120, 25, 25]) / 10

# grouped alphaR and near alphaR
alphaR_x = np.array([-180, -180, -30, -30, -180]) / 10
alphaR_y = np.array([-80, -5, -5, -80, -80]) / 10


# data = np.loadtxt('s1/GNSRV_s1_angle.xvg', skiprows=17)
    # data2 = np.loadtxt('s2/GNSRV_s2_angle.xvg', skiprows=17)
    # data = np.concatenate((data, data2), axis=0)
    # np.savetxt('s1s2/GNSRV_both_angles.xvg', data)


# Loop over all clusteres
# cluster 5 and 6 are full data for s1 and s2
# cluster 7 is full data for combined s1 and s2

for cluster_loop in range(0, 1):
    folder = 's1s2/'

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
        data = np.loadtxt(folder + 'GNSRV_both_cluster4.xvg')
    if (cluster_loop == 4):
        folder = folder + 'cluster5/'
        data = np.loadtxt(folder + 'GNSRV_both_cluster5.xvg')
    if (cluster_loop == 5):
        data = np.loadtxt(folder + 'GNSRV_s1_angle.xvg', skiprows=17)
    if (cluster_loop == 6):
        data = np.loadtxt(folder + 'GNSRV_s2_angle.xvg', skiprows=17)
    if (cluster_loop == 7):
        data = np.loadtxt(folder + 'GNSRV_both_angles.xvg')

    print(data.shape[0])

    # Loop over all residues in the cyclic peptide
    for val_loop in range(0, 4):
        if (val_loop == 0):  # Asn
            Val_chi1 = data[:, 12].copy()
            Val_chi2 = data[:, 13].copy()
            Val_phi = data[:, 4].copy()
            Val_psi = data[:, 5].copy()
            nchi = 2
        elif (val_loop == 1):  # Ser
            Val_chi1 = data[:, 14].copy()
            Val_phi = data[:, 6].copy()
            Val_psi = data[:, 7].copy()
            nchi = 1
        elif (val_loop == 2):  # Arg
            Val_chi1 = data[:, 15].copy()
            Val_chi2 = data[:, 16].copy()
            Val_chi3 = data[:, 17].copy()
            Val_chi4 = data[:, 18].copy()
            Val_phi = data[:, 8].copy()
            Val_psi = data[:, 9].copy()
            nchi = 4
        elif (val_loop == 3):  # Val
            Val_chi1 = data[:, 19].copy()
            Val_phi = data[:, 10].copy()
            Val_psi = data[:, 11].copy()
            nchi = 1

        # Make chi values between 0 and 360
        ind0 = Val_chi1 < 0
        Val_chi1[ind0] = Val_chi1[ind0] + 360

        if (nchi >= 2):
            ind0 = Val_chi2 < 0
            Val_chi2[ind0] = Val_chi2[ind0] + 360
        if (nchi >= 3):
            ind0 = Val_chi3 < 0
            Val_chi3[ind0] = Val_chi3[ind0] + 360
        if (nchi >= 4):
            ind0 = Val_chi4 < 0
            Val_chi4[ind0] = Val_chi4[ind0] + 360

        num_steps = Val_chi1.shape[0]

        # Plot phi/psi heatmap of full trajectory
        square_chi1 = np.zeros([36, 36])
        counts = np.zeros([36, 36])
        max_phi = 0
        max_psi = 0
        for i in range(0, Val_phi.shape[0]):
            phi = Val_phi[i]
            psi = Val_psi[i]
            phi_index = (math.floor(phi / 10) + 17)
            psi_index = (math.floor(psi / 10) + 17)
            square_chi1[psi_index, phi_index] = square_chi1[psi_index, phi_index] + 1
        square_chi1 = square_chi1 / Val_phi.shape[0] / 100
        ind0 = square_chi1 == 0
        square_chi1[ind0] = 'nan'

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.set_aspect('equal', adjustable='box')
        plt.imshow(square_chi1, origin='lower', cmap='jet')
        plt.plot(beta_x + 17.5, (beta_y + 17.5), 'r')
        plt.plot(beta_x1 + 17.5, (beta_y1 + 17.5), 'r')
        plt.plot(pII_x + 17.5, (pII_y + 17.5), 'r')
        plt.plot(pII_x1 + 17.5, (pII_y1 + 17.5), 'r')
        plt.plot(alphaL_x + 17.5, (alphaL_y + 17.5), 'r')
        plt.plot(alphaR_x + 17.5, (alphaR_y + 17.5), 'r')
        plt.colorbar()
        plt.xlabel('$\phi$', fontsize=20)
        plt.ylabel('$\psi$', fontsize=20)

        locs, labels = plt.xticks()            # Get locations and labels
        plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
        locs, labels = plt.yticks()            # Get locations and labels
        plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
        plt.clim(0, 0.0005)
        plt.axis([-.5, 35.5, -.5, 35.5])

        if (val_loop == 0):
            plt.title('Asn 2', fontsize=20)
            plt.savefig(folder + 'Asn_2_phipsiAll.png', bbox_inches='tight')
        elif (val_loop == 1):
            plt.title('Ser 3', fontsize=20)
            plt.savefig(folder + 'Ser_3_phipsiAll.png', bbox_inches='tight')
        elif (val_loop == 2):
            plt.title('Arg 4', fontsize=20)
            plt.savefig(folder + 'Arg_4_phipsiAll.png', bbox_inches='tight')
        elif (val_loop == 3):
            plt.title('Val 5', fontsize=20)
            plt.savefig(folder + 'Val_5_phipsiAll.png', bbox_inches='tight')
        plt.close()
        plt.close()

        # Bin phi/psi by chi1 rotamer bins
        ind0 = Val_chi1 < 120
        ind1 = Val_chi1 >= 0
        Val_60_phi = Val_phi[ind0 & ind1]
        Val_60_psi = Val_psi[ind0 & ind1]
        print(Val_60_phi.shape[0] / Val_phi.shape[0])

        ind0 = Val_chi1 < 240
        ind1 = Val_chi1 >= 120
        Val_180_phi = Val_phi[ind0 & ind1]
        Val_180_psi = Val_psi[ind0 & ind1]
        print(Val_180_phi.shape[0] / Val_phi.shape[0])

        ind0 = Val_chi1 < 360
        ind1 = Val_chi1 >= 240
        Val_300_phi = Val_phi[ind0 & ind1]
        Val_300_psi = Val_psi[ind0 & ind1]
        print(Val_300_phi.shape[0] / Val_phi.shape[0])

        # Plot each 1D chi distribution
        bins = np.arange(0, 360, 10) + 5

        for loop_chi in range(0, nchi):
            if (loop_chi == 0):
                this_chi = Val_chi1
                save_str = 'chi1'
            elif (loop_chi == 1):
                this_chi = Val_chi2
                save_str = 'chi2'
            elif (loop_chi == 2):
                this_chi = Val_chi3
                save_str = 'chi3'
            elif (loop_chi == 3):
                this_chi = Val_chi4
                save_str = 'chi4'

            n = np.zeros([36, 1])
            for i in range(0, this_chi.shape[0]):
                chi = this_chi[i]
                chi_index = math.floor(chi / 10)
                n[chi_index] = n[chi_index] + 1
            n = n / this_chi.shape[0]

            plt.plot(bins[0:36], n, 'k')
            plt.xticks(np.arange(0, 360, step=60))

            plt.axis([0, 360, 0, 0.4])
            if (loop_chi == 0):
                plt.xlabel('$\chi_1$')
                plt.ylabel('P($\chi_1$)')
            elif (loop_chi == 1):
                plt.xlabel('$\chi_2$')
                plt.ylabel('P($\chi_2$)')
            elif (loop_chi == 2):
                plt.xlabel('$\chi_3$')
                plt.ylabel('P($\chi_3$)')
            elif (loop_chi == 3):
                plt.xlabel('$\chi_4$')
                plt.ylabel('P($\chi_4$)')
            if (val_loop == 0):
                plt.title('Asn 2', fontsize=20)
                plt.savefig(folder + 'Asn_2_P' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 1):
                plt.title('Ser 3', fontsize=20)
                plt.savefig(folder + 'Ser_3_P' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 2):
                plt.title('Arg 4', fontsize=20)
                plt.savefig(folder + 'Arg_4_P' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 3):
                plt.title('Val 5', fontsize=20)
                plt.savefig(folder + 'Val_5_P' + save_str + '.png', bbox_inches='tight')
            plt.close()
            plt.close()

        # Plot chi1/chi2 if relevant
        if (nchi >= 2):
            n = np.zeros([36, 36])
            for i in range(0, this_chi.shape[0]):
                chi1 = Val_chi1[i]
                chi2 = Val_chi2[i]
                chi1_index = math.floor(chi1 / 10)
                chi2_index = math.floor(chi2 / 10)
                n[chi2_index, chi1_index] = n[chi2_index, chi1_index] + 1
            n = n / this_chi.shape[0]
            square_chi1 = n.copy()
            ind0 = square_chi1 == 0
            square_chi1[ind0] = 'nan'
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111)
            ax.set_aspect('equal', adjustable='box')
            plt.imshow(square_chi1, origin='lower', cmap='jet')
            plt.colorbar()
            plt.xlabel('$\chi_1$', fontsize=20)
            plt.ylabel('$\chi_2$', fontsize=20)

            locs, labels = plt.xticks()            # Get locations and labels
            plt.xticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
            locs, labels = plt.yticks()            # Get locations and labels
            plt.yticks(np.arange(0, 36, 6), np.arange(0, 360, 60), fontsize=20)
            plt.axis([-.5, 35.5, -.5, 35.5])
            if (val_loop == 0):
                plt.title('Asn 2', fontsize=20)
                plt.plot([-0.5, 35.5], [2.3, 2.3], 'r')
                plt.plot([-0.5, 35.5], [5.3, 5.3], 'r')
                plt.plot([-0.5, 35.5], [8.3, 8.3], 'r')
                plt.plot([-0.5, 35.5], [11.3, 11.3], 'r')
                plt.plot([-0.5, 35.5], [14.3, 14.3], 'r')
                plt.plot([-0.5, 35.5], [17.3, 17.3], 'r')
                plt.plot([-0.5, 35.5], [20.3, 20.3], 'r')
                plt.plot([-0.5, 35.5], [23.3, 23.3], 'r')
                plt.plot([-0.5, 35.5], [26.3, 26.3], 'r')
                plt.plot([-0.5, 35.5], [29.3, 29.3], 'r')
                plt.plot([-0.5, 35.5], [32.3, 32.3], 'r')
                plt.plot([-0.5, 35.5], [35.3, 35.3], 'r')
                plt.plot([12, 12], [-0.5, 35.5], 'r')
                plt.plot([24, 24], [-0.5, 35.5], 'r')

                plt.savefig(folder + 'Asn_2_chi1_chi2.png', bbox_inches='tight')

                # get rotamer values
                p_each_rot = np.zeros([12 * 3, 3])
                for chi1_loop in range(0, 3):
                    min_c1 = chi1_loop * 120
                    max_c1 = min_c1 + 120

                    for chi2_loop in range(0, 12):
                        min_c2 = (chi2_loop - 1) * 30 + 23
                        max_c2 = min_c2 + 30

                        ind0 = Val_chi1 > min_c1
                        ind1 = Val_chi1 < max_c1
                        ind2 = Val_chi2 > min_c2
                        ind3 = Val_chi2 < max_c2
                        p_each_rot[(chi1_loop * 12) + chi2_loop, 2] = sum(ind0 & ind1 & ind2 & ind3)

                p_each_rot[:, 2] = p_each_rot[:, 2] / Val_chi1.shape[0] * 100
            elif (val_loop == 1):
                plt.title('Ser 3', fontsize=20)
                plt.savefig(folder + 'Ser_3_chi1_chi2.png', bbox_inches='tight')
            elif (val_loop == 2):
                plt.title('Arg 4', fontsize=20)
                plt.savefig(folder + 'Arg_4_chi1_chi2.png', bbox_inches='tight')
            elif (val_loop == 3):
                plt.title('Val 5', fontsize=20)
                plt.savefig(folder + 'Val_5_chi1_chi2.png', bbox_inches='tight')
            plt.close()

        # Loop over chi1 rotamers and plot phi/psi
        for chi_rotamers in range(0, 3):
            if (chi_rotamers == 0):    # 60
                this_phi = Val_60_phi
                this_psi = Val_60_psi
                save_str = '60'
            elif (chi_rotamers == 1):  # 180
                this_phi = Val_180_phi
                this_psi = Val_180_psi
                save_str = '180'
            elif (chi_rotamers == 2):  # 300
                this_phi = Val_300_phi
                this_psi = Val_300_psi
                save_str = '300'

            square_chi1 = np.zeros([36, 36])
            counts = np.zeros([36, 36])
            for i in range(0, this_phi.shape[0]):
                phi = this_phi[i]
                phi_index = (math.floor(phi / 10) + 17)
                psi = this_psi[i]
                psi_index = (math.floor(psi / 10) + 17)
                square_chi1[psi_index, phi_index] = square_chi1[psi_index, phi_index] + 1

            if (chi_rotamers == 0):
                square_chi1_60 = square_chi1.copy()
            elif (chi_rotamers == 1):
                square_chi1_180 = square_chi1.copy()
            elif (chi_rotamers == 2):
                square_chi1_300 = square_chi1.copy()
            square_chi1 = square_chi1
            ind0 = square_chi1 == 0
            square_chi1[ind0] = 'nan'

            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111)
            ax.set_aspect('equal', adjustable='box')
            plt.imshow(square_chi1 / num_steps / 100, origin='lower', cmap='jet')
            plt.plot(beta_x + 17.5, (beta_y + 17.5), 'r')
            plt.plot(beta_x1 + 17.5, (beta_y1 + 17.5), 'r')
            plt.plot(pII_x + 17.5, (pII_y + 17.5), 'r')
            plt.plot(pII_x1 + 17.5, (pII_y1 + 17.5), 'r')
            plt.plot(alphaL_x + 17.5, (alphaL_y + 17.5), 'r')
            plt.plot(alphaR_x + 17.5, (alphaR_y + 17.5), 'r')
            plt.colorbar()
            plt.xlabel('$\phi$', fontsize=20)
            plt.ylabel('$\psi$', fontsize=20)

            locs, labels = plt.xticks()            # Get locations and labels
            plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
            locs, labels = plt.yticks()            # Get locations and labels
            plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
            plt.clim(0, 0.0001)
            plt.axis([-.5, 35.5, -.5, 35.5])

            plt.title('$\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)

            if (val_loop == 0):
                plt.title('Asn 2', fontsize=20)
                plt.savefig(folder + 'Asn_2_' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 1):
                plt.title('Ser 3', fontsize=20)
                plt.savefig(folder + 'Ser_3_' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 2):
                plt.title('Arg 4', fontsize=20)
                plt.savefig(folder + 'Arg_4_' + save_str + '.png', bbox_inches='tight')
            elif (val_loop == 3):
                plt.title('Val 5', fontsize=20)
                plt.savefig(folder + 'Val_5_' + save_str + '.png', bbox_inches='tight')

            plt.close()

            # Normalize in plot
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111)
            ax.set_aspect('equal', adjustable='box')
            plt.imshow(square_chi1 / this_phi.shape[0] / 100, origin='lower', cmap='jet')
            plt.colorbar()
            plt.xlabel('$\phi$', fontsize=20)
            plt.ylabel('$\psi$', fontsize=20)
            plt.title('$\chi_1$ = 60$^{\circ}$', fontsize=20)
            locs, labels = plt.xticks()            # Get locations and labels
            plt.xticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
            locs, labels = plt.yticks()            # Get locations and labels
            plt.yticks(np.arange(0, 36, 6), np.arange(-180, 180, 60), fontsize=20)
            plt.clim(0, 0.0007)

            plt.axis([-.5, 35.5, -.5, 35.5])

            if (val_loop == 0):
                plt.title('Asn 2 $\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)
                plt.savefig(folder + 'Asn_2_' + save_str + '_normIn.png', bbox_inches='tight')
            elif (val_loop == 1):
                plt.title('Ser 3 $\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)
                plt.savefig(folder + 'Ser_3_' + save_str + '_normIn.png', bbox_inches='tight')
            elif (val_loop == 2):
                plt.title('Arg 4 $\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)
                plt.savefig(folder + 'Arg_4_' + save_str + '_normIn.png', bbox_inches='tight')
            elif (val_loop == 3):
                plt.title('Val 5 $\chi_1$ = ' + save_str + '$^{\circ}$', fontsize=20)
                plt.savefig(folder + 'Val_5_' + save_str + '_normIn.png', bbox_inches='tight')

            plt.close()
