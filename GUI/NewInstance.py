"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import OrbitInputPage, FamilyInputPage

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
from tkinter import filedialog, ttk
from Utility import Utility, System, Plot, NumericalMethods
from Orbit import Orbit, InitialGuess
from OrbitFamily import OrbitFamily
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from PIL import ImageTk, Image
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import os

# global variables
LARGE_FONT = ("Verdana", 12)



# ----------------------------------------------------------------------------------------------------------------------
# This class specifies whether a single orbit or a family is to be calculated.
# ----------------------------------------------------------------------------------------------------------------------
class NewInstance(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # The GUI structure is loaded during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        # single orbit
        button = Button(self, text="Single Orbit", command=lambda: NewInstance.singleOrbit(self, parent, controller, dynamicalSystem))
        button.pack()
        # orbit family
        button2 = Button(self, text="Orbit Family", command=lambda: NewInstance.orbitFamily(self, parent, controller, dynamicalSystem))
        button2.pack()
        # back to home
        homeButton = Button(self, text="Home", command=lambda: [controller.show_frame("StartPage"), self.destroy()])
        homeButton.pack()
    # ------------------------------------------------------------------------------------------------------------------
    # This method calls the page to enter the orbit parameters.
    # ------------------------------------------------------------------------------------------------------------------
    def singleOrbit(self, parent, controller, dynamicalSystem):
        controller.frames["OrbitInputPage"] = OrbitInputPage.OrbitInputPage(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem)
        controller.frames["OrbitInputPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitInputPage")
        self.update()
        time.sleep(0.5)
        self.destroy()
    # ------------------------------------------------------------------------------------------------------------------
    # This method calls the page to enter the orbit family parameters.
    # ------------------------------------------------------------------------------------------------------------------
    def orbitFamily(self, parent, controller, dynamicalSystem):
        controller.frames["FamilyInputPage"] = FamilyInputPage.FamilyInputPage(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem)
        controller.frames["FamilyInputPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("FamilyInputPage")
        self.update()
        time.sleep(0.5)
        self.destroy()
