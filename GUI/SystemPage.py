"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import NewInstance

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
# This class configures the dynamic system.
# ----------------------------------------------------------------------------------------------------------------------
class SystemPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # The GUI structure is loaded during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        # dynamical system
        stepOne = LabelFrame(self, text=" DYNAMICAL SYSTEM ")
        stepOne.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        describtion = Label(stepOne, text="Describtion")
        describtion.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        nameFPLabel = Label(stepOne, text="Name of first Primary:")
        nameFPLabel.grid(row=1, column=0, sticky='W', padx=5, pady=2)
        self.nameFPEntry = Entry(stepOne)
        self.nameFPEntry.grid(row=1, column=1, sticky='W', padx=5, pady=2)
        massFPLabel = Label(stepOne, text="Mass of first Primary:")
        massFPLabel.grid(row=2, column=0, sticky='W', padx=5, pady=2)
        self.massFPEntry = Entry(stepOne)
        self.massFPEntry.grid(row=2, column=1, sticky='W', padx=5, pady=2)
        unit = Label(stepOne, text="kg")
        unit.grid(row=2, column=2, sticky='W')
        nameSPLabel = Label(stepOne, text="Name of second Primary:")
        nameSPLabel.grid(row=3, column=0, sticky='W', padx=5, pady=2)
        self.nameSPEntry = Entry(stepOne)
        self.nameSPEntry.grid(row=3, column=1, sticky='W', padx=5, pady=2)
        massSPLabel = Label(stepOne, text="Mass of second Primary:")
        massSPLabel.grid(row=4, column=0, sticky='W', padx=5, pady=2)
        self.massSPEntry = Entry(stepOne)
        self.massSPEntry.grid(row=4, column=1, sticky='W', padx=5, pady=2)
        unit = Label(stepOne, text="kg")
        unit.grid(row=4, column=2, sticky='W')
        distanceLabel = Label(stepOne, text="Distance of Primaries:")
        distanceLabel.grid(row=5, column=0, sticky='W', padx=5, pady=2)
        self.distanceEntry = Entry(stepOne)
        self.distanceEntry.grid(row=5, column=1, sticky='W', padx=5, pady=2)
        unit = Label(stepOne, text="m")
        unit.grid(row=5, column=2, sticky='W')
        # defaul systems
        defaultSystems = LabelFrame(self, text=" DEFAULT SYSTEMS ")
        defaultSystems.grid(row=0, column=9, columnspan=2, rowspan=7, sticky='NS', padx=5, pady=5)
        var = IntVar()
        earthMoon = Radiobutton(defaultSystems, text="Earth - Moon System", variable=var, value=0, command=lambda: SystemPage.earthMoon(self))
        earthMoon.grid(row=0, sticky='W', padx=5, pady=2)
        sunEarth = Radiobutton(defaultSystems, text="Sun - Earth System", variable=var, value=1, command=lambda: SystemPage.sunEarth(self))
        sunEarth.grid(row=1, sticky='W', padx=5, pady=2)
        sunMars = Radiobutton(defaultSystems, text="Sun - Mars System", variable=var, value=2, command=lambda: SystemPage.sunMars(self))
        sunMars.grid(row=2, sticky='W', padx=5, pady=2)
        sunJupiter = Radiobutton(defaultSystems, text="Sun - Jupiter System", variable=var, value=3, command=lambda: SystemPage.sunJupiter(self))
        sunJupiter.grid(row=3, sticky='W', padx=5, pady=2)
        # control
        confirmButton = Button(self, text="Confirm", command=lambda: SystemPage.confirm(self, parent, controller))
        confirmButton.grid(row=7, column=1, sticky='E', padx=5, pady=2)
        homeButton = Button(self, text="Home", command=lambda: [controller.show_frame("StartPage"), self.destroy()])
        homeButton.grid(row=7, column=0, sticky='W', padx=5, pady=2)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets corresponding default values in the input.
    # ------------------------------------------------------------------------------------------------------------------
    def earthMoon(self):
        self.nameFPEntry.delete(0, "end")
        self.massFPEntry.delete(0, "end")
        self.nameSPEntry.delete(0, "end")
        self.massSPEntry.delete(0, "end")
        self.distanceEntry.delete(0, "end")
        self.nameFPEntry.insert(0, "Earth")
        self.massFPEntry.insert(0, 5.97237e+24)
        self.nameSPEntry.insert(0, "Moon")
        self.massSPEntry.insert(0, 7.342e+22)
        self.distanceEntry.insert(0, 384402 * 1.0e3)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets corresponding default values in the input.
    # ------------------------------------------------------------------------------------------------------------------
    def sunEarth(self):
        self.nameFPEntry.delete(0, "end")
        self.massFPEntry.delete(0, "end")
        self.nameSPEntry.delete(0, "end")
        self.massSPEntry.delete(0, "end")
        self.distanceEntry.delete(0, "end")
        self.nameFPEntry.insert(0, "Sun")
        self.massFPEntry.insert(0, 1.98892e+30)
        self.nameSPEntry.insert(0, "Earth")
        self.massSPEntry.insert(0, 5.97237e+24)
        self.distanceEntry.insert(0, 149597870.7 * 1.0e3)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets corresponding default values in the input.
    # ------------------------------------------------------------------------------------------------------------------
    def sunMars(self):
        self.nameFPEntry.delete(0, "end")
        self.massFPEntry.delete(0, "end")
        self.nameSPEntry.delete(0, "end")
        self.massSPEntry.delete(0, "end")
        self.distanceEntry.delete(0, "end")
        self.nameFPEntry.insert(0, "Sun")
        self.massFPEntry.insert(0, 1.98892e+30)
        self.nameSPEntry.insert(0, "Mars")
        self.massSPEntry.insert(0, 6.419e+23)
        self.distanceEntry.insert(0, 227900000 * 1.0e3)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets corresponding default values in the input.
    # ------------------------------------------------------------------------------------------------------------------
    def sunJupiter(self):
        self.nameFPEntry.delete(0, "end")
        self.massFPEntry.delete(0, "end")
        self.nameSPEntry.delete(0, "end")
        self.massSPEntry.delete(0, "end")
        self.distanceEntry.delete(0, "end")
        self.nameFPEntry.insert(0, "Sun")
        self.massFPEntry.insert(0, 1.98892e+30)
        self.nameSPEntry.insert(0, "Jupiter")
        self.massSPEntry.insert(0, 1.89813e+27)
        self.distanceEntry.insert(0, 778547200 * 1.0e3)
    # ------------------------------------------------------------------------------------------------------------------
    # This method creates an object of the dynamical system from the input.
    # ------------------------------------------------------------------------------------------------------------------
    def confirm(self, parent, controller):
        try:
            nameFP = self.nameFPEntry.get()
            massFP = float(self.massFPEntry.get())
            nameSP = self.nameSPEntry.get()
            massSP = float(self.massSPEntry.get())
            distance = float(self.distanceEntry.get())
        except ValueError:
            self.nameFPEntry.config(highlightbackground='red')
            self.massFPEntry.config(highlightbackground='red')
            self.nameSPEntry.config(highlightbackground='red')
            self.massSPEntry.config(highlightbackground='red')
            self.distanceEntry.config(highlightbackground='red')
        else:
            dynamicalSystem = System(nameFP=nameFP, massFP=massFP, nameSP=nameSP, massSP=massSP, distance=distance)
            controller.frames["NewInstance"] = NewInstance.NewInstance(parent=parent, controller=controller,
                                                           dynamicalSystem=dynamicalSystem)
            controller.frames["NewInstance"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("NewInstance")
            self.update()
            time.sleep(0.5)
            self.destroy()
