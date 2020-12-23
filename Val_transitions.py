import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

data = np.loadtxt('s1/GNSRV_s1_angle.xvg',skiprows=17)

Val_chi1 = data[:,19]
ind0 = Val_chi1 < 0
Val_chi1[ind0] = Val_chi1[ind0]+360

Val_phi = data[:,10]
Val_psi = data[:,11]


Val_chi_binned = np.zeros(Val_chi1.shape,dtype=int)

ind0 = Val_chi1 <= 120
ind1 = Val_chi1 >= 0
Val_chi_binned[ind0&ind1]= 0
Val_60 = sum(ind0&ind1)

ind0 = Val_chi1 <= 240
ind1 = Val_chi1 >= 120
Val_chi_binned[ind0&ind1] = 1
Val_180 = sum(ind0&ind1)

ind0 = Val_chi1 <= 360
ind1 = Val_chi1 >= 240
Val_chi_binned[ind0&ind1] = 2
Val_300 = sum(ind0&ind1)

print(Val_60 + Val_180+Val_300)
print(Val_chi1.shape)

Val_transitions = np.zeros([3,3])
for i in range(0,Val_chi1.shape[0]-1):
    if ((i)%2500 > 10):   # i+1 because 0 indexed
        x = Val_chi_binned[i]
        y = Val_chi_binned[i+1]
        Val_transitions[x,y] = Val_transitions[x,y]+1

print(Val_transitions/250000)
#np.savetxt('Val_chi_binned.txt', Val_chi_binned)


# All:
#  15642.    721.    216.
#    719. 179333.   2305.
#    218.   2303.  48542.

# no (i+1)%2500
# 15638.    719.    216.
#   716. 179268.   2299.
#   218.   2296.  48530.

# no (i+1)%2500 < 10
# 15571.    717.    215.
#   714. 178526.   2288.
#   217.   2285.  48367.