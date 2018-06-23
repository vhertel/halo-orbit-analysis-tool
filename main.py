"""
File    : main.py
Author  : Victor Hertel
Date    : 28.05.2018

Computation of halo orbits
"""

# Imports
from GUI import HaloTool
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from Utility import System, Utility, Plot
from Orbit import Orbit
from scipy.integrate import odeint


app = HaloTool()
app.mainloop()


# # Initial Conditions
# test = np.array([8.23390186e-01, 0, 0, 0, 1.26668772e-01, 0])
# l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
# l1_middle = np.array([0.8235990912, 0, -0.0399866715, 0, 0.1492106867, 0])
# l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
# l2_middle = np.array([1.1542349115, 0, -0.1379744940, 0, -0.2147411949, 0])
#
# earth_moon = System(nameFP="Earth", massFP=5.97237e+24, nameSP="Moon", massSP=7.342e+22, distance=384402 * 1.0e3)
#
# o = Orbit(l2, "x", earth_moon)
# o.getClosestNRHO()
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# t = np.linspace(0, 1.25*o.period, num=300)
# orbitStates = odeint(Utility.sysEquations, o.x0, t, args=(o.system.mu,), rtol=2.5e-13, atol=1e-22)
# xOrbit = orbitStates[:, 0]
# yOrbit = orbitStates[:, 1]
# zOrbit = orbitStates[:, 2]
# ax.plot(xOrbit, yOrbit, zOrbit, color='black', linewidth=1)
# ax.scatter((1 - o.system.mu), 0, 0, color='grey', s=1, label="Moon")
# #ax.scatter((-o.system.mu), 0, 0, color='blue', s=3, label="Moon")
# Plot.setAxesEqual(ax)
# plt.show()
