"""
File    : main.py
Author  : Victor Hertel
Date    : 28.05.2018

Computation of halo orbits
"""

# Imports
import numpy as np
from Orbit import Orbit
from OrbitFamily import OrbitFamily
from Utility import Plot, System



# Initial Conditions
l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
l2_middle = np.array([1.1542349115, 0, -0.1379744940, 0, -0.2147411949, 0])


earth_moon = System(nameFP="Earth", massFP=5.97237e+24, nameSP="Moon", massSP=7.342e+22, distance=384402*1.0e3)

orbit = Orbit(l1, "x", earth_moon)
orbit.plot()

#fam = OrbitFamily(l1, 0.007, earth_moon, orbitNumber=5)
#fam.getHaloFamily()
#fam.plot()

#Plot.plotFromTable("L2")
