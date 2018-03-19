
"""
File    : main.py
Author  : Victor Hertel
Date    : 18.03.2018

Computation of halo orbits
"""

# Imports
import numpy as np
import matplotlib as mpl
# necessary to use matplotlib for Mac
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from library import *



# Initial Condition
l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
l2 = np.array([1.1808881373,  0, -0.0032736457,  0, -0.1559184478,  0])

t0 = 0
mFirstPrimary = 5.97237e+24           # Earth
mSecondPrimary = 7.342e+22            # Moon
mu = mSecondPrimary / (mSecondPrimary + mFirstPrimary)



# Prepares figure for plot
fig = plt.figure()
ax = fig.gca(projection='3d')

#--------------------------------------------------------------------------
# SINGLE HALO ORBIT
calculation.singleHalo(l1, t0, mu, 1.0e-8, "x", "southern", ax)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# NATURAL PARAMETER CONTINUATION
#calculation.natParaConti(l1, t0, mu, 1.0e-6, "all", 0.002, "L1", "northern", ax)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# PSEUDO-ARCLENGTH CONTINUATION
#calculation.pseudoArcLenConti(l2, t0, mu, 1.0e-6, 20, 0.05, "vertical", "southern", ax)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# PLOT OF PRIMARIES AND LAGRANGIAN POINTS L1 AND L2
calculation.primaries(mu, "second", ax)
calculation.lagrangianPoints(mu, ax)
#--------------------------------------------------------------------------

# Equals axis lengths and plots figure
utilities.setAxesEqual(ax)
utilities.plotTraj(ax)
