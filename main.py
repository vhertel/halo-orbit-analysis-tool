"""
File    : main.py
Author  : Victor Hertel
Date    : 24.04.2018

Computation of halo orbits
"""

# Imports
import numpy as np
from Orbit import Orbit
from OrbitFamily import OrbitFamily

from Utility import NumericalMethods, Plot
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Initial Conditions
l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
test = np.array([1.10156010, 0, -0.19840472, 0, -0.21676273, 0])

mFirstPrimary = 5.97237e+24  # Earth
mSecondPrimary = 7.342e+22  # Moon
mu = mSecondPrimary / (mSecondPrimary + mFirstPrimary)



L1 = np.array([0.8627102723, 0, -0.4702996003, 0, 0.1351811306, 0])
LTEST = np.array([0.9214943426, 0, -0.3238011029, 0, 0.0848204819, 0])

test = np.array([0.9325284505, 0, -0.2384555359, 0, 0.0961225566, 0])


fam = OrbitFamily(test, 0.0045, 5, mu)
fam.getHaloFamily()
fam.writeData()
fam.plot(background="on", haloFamily="southern")




# period 2.7159986012
L1 = np.array([0.8627102723, 0, -0.4702996003, 0, 0.1351811306, 0])
# period 2.3534714874
bla = np.array([0.8631244803, 0, -0.1859405036, 0, 0.2500788516, 0])

halo = NumericalMethods.rk4System(test, 2.5, mu, 2000)
x = halo[:, 0]
y = halo[:, 1]
z = halo[:, 2]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, color='blue', linewidth=1)
ax.scatter(1-mu, 0, 0, color='black', s = 8, label='Moon')
Plot.setAxesEqual(ax)
plt.show()
