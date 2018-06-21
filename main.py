"""
File    : main.py
Author  : Victor Hertel
Date    : 28.05.2018

Computation of halo orbits
"""

# Imports
from GUI import HaloTool
import numpy as np
from Orbit import Orbit
from OrbitFamily import OrbitFamily
from Utility import Plot, System



import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from Utility import Utility, Plot
from scipy.integrate import odeint



app = HaloTool()
app.mainloop()

























# # Initial Conditions
# l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
# l1_middle = np.array([0.8235990912, 0, -0.0399866715, 0, 0.1492106867, 0])
# l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
# l2_middle = np.array([1.1542349115, 0, -0.1379744940, 0, -0.2147411949, 0])

# earth_moon = System(nameFP="Earth", massFP=5.97237e+24, nameSP="Moon", massSP=7.342e+22, distance=384402 * 1.0e3)
# o = Orbit(l1_middle, "x", earth_moon)
# o.stableManifold(numberOfPoints=30, durationFactor=1.35)
# o.unstableManifold(numberOfPoints=30, durationFactor=1.35)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# t = np.linspace(0, 1.25*o.period, num=10000)
# for i in range(len(o.stableManifolds)):
#     manifold = odeint(Utility.backwards, o.stableManifolds[i], t, args=(o.system.mu,), rtol=2.5e-13, atol=1e-22)
#     x = manifold[:, 0]
#     y = manifold[:, 1]
#     z = manifold[:, 2]
#     ax.plot(x, y, z, color='green', linewidth=0.5)
# for i in range(len(o.unstableManifolds)):
#     manifold = odeint(Utility.sysEquations, o.unstableManifolds[i], t, args=(o.system.mu,), rtol=2.5e-13, atol=1e-22)
#     x = manifold[:, 0]
#     y = manifold[:, 1]
#     z = manifold[:, 2]
#     ax.plot(x, y, z, color='blue', linewidth=0.5)
# orbitStates = odeint(Utility.sysEquations, o.x0, t, args=(o.system.mu,), rtol=2.5e-13, atol=1e-22)
# xOrbit = orbitStates[:, 0]
# yOrbit = orbitStates[:, 1]
# zOrbit = orbitStates[:, 2]
# ax.plot(xOrbit, yOrbit, zOrbit, color='black', linewidth=1)
# ax.scatter((1 - o.system.mu), 0, 0, color='grey', s=1, label="Moon")
# #ax.scatter((-o.system.mu), 0, 0, color='blue', s=3, label="Moon")
# #Plot.setAxesEqual(ax)
# ax.set_xlim3d(0.6, 0.9)
# ax.set_ylim3d(-0.15,0.15)
# ax.set_zlim3d(-0.15,0.15)
# plt.show()


