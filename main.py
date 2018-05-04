"""
File    : main.py
Author  : Victor Hertel
Date    : 24.04.2018

Computation of halo orbits
"""

# Imports
import numpy as np
from Orbit import Orbit
from OrbitFamily import OrbitFamily, L1Family, L2Family

# Initial Conditions
l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
test = np.array([1.10156010, 0, -0.19840472, 0, -0.21676273, 0])

mFirstPrimary = 5.97237e+24  # Earth
mSecondPrimary = 7.342e+22  # Moon
mu = mSecondPrimary / (mSecondPrimary + mFirstPrimary)




Orbit = Orbit(test, "x", mu)
Orbit.plot(background="on", haloFamily="southern")

