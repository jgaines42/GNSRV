# analyze h-bonds from results of hbond command in gromacs

import numpy as np  # type: ignore

folder = 's1/'
cluster = 'c3'
aa_list = ['G1', 'N2', 'S3', 'R4', 'V5']
aa_hbond_data = np.zeros([5, 3])
for aa_loop in range(0, 5):
    aa = aa_list[aa_loop]

    data = np.loadtxt(folder + 'hbond_' + cluster + '_' + aa + '.xvg', skiprows=25)
    print(data[0, :])
    n_hbonds = sum(data[:, 1])
    n_contacts = sum(data[:, 2])
    x = data[:, 1] + data[:, 2]
    aa_hbond_data[aa_loop, 0] = n_hbonds / data.shape[0]
    aa_hbond_data[aa_loop, 1] = n_contacts / data.shape[0]
    aa_hbond_data[aa_loop, 2] = (n_hbonds + n_contacts) / data.shape[0]
    print((n_hbonds + n_contacts))

print(aa_hbond_data)
print(data.shape[0])
