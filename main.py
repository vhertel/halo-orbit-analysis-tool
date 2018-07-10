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
from Orbit import Orbit, InitialGuess
from scipy.integrate import odeint


app = HaloTool()
app.mainloop()


#Plot.plotFromTable("L1&L2")















