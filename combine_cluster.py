# Combine data from both s1 and s2 clusters

import numpy as np  # type: ignore

for cluster_loop in range(1, 8):

    data = np.loadtxt('s1/cluster' + str(cluster_loop) + '/s1_cluster' + str(cluster_loop) + '.xvg', skiprows=17)
    data2 = np.loadtxt('s2/s2_cluster' + str(cluster_loop) + '.xvg', skiprows=17)
    data = np.concatenate((data, data2), axis=0)
    np.savetxt('s1s2/cluster' + str(cluster_loop) + '/GNSRV_both_cluster' + str(cluster_loop) + '.xvg', data, fmt='%10.3f')
