"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import OrbitFamilyPage


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
# This class contains the configuration of the Orbit family.
# ----------------------------------------------------------------------------------------------------------------------
class FamilyInputPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # The GUI structure is loaded during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        # initial state
        stateFrame = LabelFrame(self, text=" LAGRANGIAN ")
        stateFrame.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        self.lagrangian = StringVar(value="L1")
        Radiobutton(stateFrame, text="L1", variable=self.lagrangian, value="L1").grid(row=0, sticky='W', padx=5, pady=2)
        Radiobutton(stateFrame, text="L2", variable=self.lagrangian, value="L2").grid(row=0, column=1, sticky='W', padx=5, pady=2)
        self.family = StringVar(value="Northern")
        Radiobutton(stateFrame, text="Northern", variable=self.family, value="Northern").grid(row=1, sticky='W', padx=5, pady=2)
        Radiobutton(stateFrame, text="Southern", variable=self.family, value="Southern").grid(row=1, column=1, sticky='W', padx=5, pady=2)

        Label(stateFrame, text="Orbit Distance:").grid(row=2, sticky='W', padx=5)
        self.distanceEntry = Entry(stateFrame, width=7)
        self.distanceEntry.grid(row=2, column=1, sticky='E', padx=5, pady=2)
        self.distance = 0.001
        self.distanceEntry.insert(END, self.distance)

        Label(self, font=("TkDefaultFont", 11), text="Status:",).grid(row=2, sticky='NW', padx=5, pady=(15,5))
        self.statusBar = Text(self, height=15, relief=RIDGE, bd=3, spacing1=3)
        self.statusBar.grid(row=2, column=1, columnspan=3, sticky='W', padx=5, pady=(15,5))
        self.statusBar.insert(INSERT, time.strftime("%Y-%m-%dT%H.%M.%S:\n\n>>>"))

        # button for calculation of entry data
        Button(self, text="Calculate", command=lambda: FamilyInputPage.calculate(self, dynamicalSystem)).grid(row=3, column=1, sticky='E', padx=5, pady=2)
        Button(self, text="Next", command=lambda: FamilyInputPage.nextPage(self, parent, controller, dynamicalSystem)).grid(row=3, column=2, sticky='E', padx=5, pady=2)
        # back to home
        Button(self, text="Home", command=lambda: [controller.show_frame("StartPage"), self.destroy()]).grid(row=3, column=0, sticky='W', padx=5, pady=2)
    # ------------------------------------------------------------------------------------------------------------------
    # This method calculates the Orbit family.
    # ------------------------------------------------------------------------------------------------------------------
    def calculate(self, dynamicalSystem):
        guess = InitialGuess(dynamicalSystem, self.lagrangian.get(), self.family.get(), "z", 0)
        self.orbitFamily = OrbitFamily(guess.x0, self.family.get(), self.lagrangian.get(), float(self.distanceEntry.get()), dynamicalSystem, statusBar=self.statusBar)
    # ------------------------------------------------------------------------------------------------------------------
    # This method calls the detail page for the orbit family found.
    # ------------------------------------------------------------------------------------------------------------------
    def nextPage(self, parent, controller, dynamicalSystem):
            controller.frames["OrbitFamilyPage"] = OrbitFamilyPage.OrbitFamilyPage(parent=parent, controller=controller, orbitFamily=self.orbitFamily, dynamicalSystem=dynamicalSystem)
            controller.frames["OrbitFamilyPage"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("OrbitFamilyPage")
            self.update()
            time.sleep(0.5)
            self.destroy()
