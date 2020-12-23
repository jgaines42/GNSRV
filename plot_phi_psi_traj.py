###################################################################
# plot_phi_psi_traj.py
#
# Plots the phi/psi trajectory for all 5 amino acids in GNSRV
# from MD simulation specified in "file_name"
###################################################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

file_name = 's2/GNSRV_s2_angle.xvg'
folder = 's2/'

# Load data
data = np.loadtxt(file_name, skiprows=17)
time = data[:, 0] / 100

# Loop over all amino acids
for aa_loop in range(0, 5):
    if (aa_loop == 0):
        Val_chi1 = data[:, 12].copy()
        Val_phi = data[:, 4].copy()
        Val_psi = data[:, 5].copy()
    elif (aa_loop == 1):
        Val_chi1 = data[:, 14].copy()
        Val_phi = data[:, 6].copy()
        Val_psi = data[:, 7].copy()
    elif (aa_loop == 2):
        Val_chi1 = data[:, 15].copy()
        Val_phi = data[:, 8].copy()
        Val_psi = data[:, 9].copy()
    elif (aa_loop == 3):
        Val_chi1 = data[:, 19].copy()
        Val_phi = data[:, 10].copy()
        Val_psi = data[:, 11].copy()
    elif (aa_loop == 4):
        Val_phi = data[:, 2].copy()
        Val_psi = data[:, 3].copy()

    # Set chi1 to be from 0 to 360
    ind0 = Val_chi1 < 0
    Val_chi1[ind0] = Val_chi1[ind0] + 360

    # Make plot
    plt.rcParams.update({'font.size': 15})
    fig, axs = plt.subplots(3, 1, figsize=(15, 10))

    if (aa_loop < 4):
        axs[0].plot(time, Val_chi1, '.', markersize=2)
        axs[0].axis([500, 1000, 0, 360])
    axs[1].plot(time, Val_phi, '.', markersize=2)
    axs[2].plot(time, Val_psi, '.', markersize=2)
    axs[1].axis([500, 1000, -180, 180])
    axs[2].axis([500, 1000, -180, 180])
    axs[0].set_ylabel('$\chi_1$')
    axs[1].set_ylabel('$\phi$')
    axs[2].set_ylabel('$\psi$')
    axs[2].set_xlabel('Time (ns)')
    axs[0].set_yticks(np.arange(0, 360, step=60))
    axs[1].set_yticks(np.arange(-180, 180, step=60))
    axs[2].set_yticks(np.arange(-180, 180, step=60))

    if (aa_loop == 0):
        axs[0].set_title('Asn 2', fontsize=20)
        plt.savefig(folder + 'Asn_2_traj.png', bbox_inches='tight')
    elif (aa_loop == 1):
        axs[0].set_title('Ser 3', fontsize=20)
        plt.savefig(folder + 'Ser_3_trajl.png', bbox_inches='tight')
    elif (aa_loop == 2):
        axs[0].set_title('Arg 4', fontsize=20)
        plt.savefig(folder + 'Arg_4_traj.png', bbox_inches='tight')
    elif (aa_loop == 3):
        axs[0].set_title('Val 5', fontsize=20)
        plt.savefig(folder + 'Val_5_traj.png', bbox_inches='tight')
    elif (aa_loop == 4):
        axs[0].set_title('Gly 1', fontsize=20)
        plt.savefig(folder + 'Gly_1_traj.png', bbox_inches='tight')
    plt.show()
