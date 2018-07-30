"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import SystemPage, OrbitPage, OrbitFamilyPage
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
# This method contains the start page.
# ----------------------------------------------------------------------------------------------------------------------
class StartPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # The GUI structure is loaded during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        # main frame
        frameStartPage = Frame(self)
        frameStartPage.pack(expand=True, fill=BOTH, padx=10, pady=10)
        # button for creating new instance of orbit or orbit family
        instanceButton = Button(frameStartPage, text="Create New Instance", command=lambda: StartPage.createNewInstance(self, parent, controller))
        instanceButton.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        # button for creating new instance by loading data
        loadButton = Button(frameStartPage, text="Load From Data", command=lambda: StartPage.loadFromData(self, parent, controller))
        loadButton.grid(row=0, column=1, sticky='E', padx=5, pady=2)
        # button to quit program
        quitButton = Button(frameStartPage, text="Exit", command=self.client_exit)
        quitButton.grid(row=1, column=0, sticky='W', padx=5, pady=2)
    # ------------------------------------------------------------------------------------------------------------------
    # This method is called when the user wants to instantiate a new object.
    # ------------------------------------------------------------------------------------------------------------------
    def createNewInstance(self, parent, controller):
        controller.frames["SystemPage"] = SystemPage.SystemPage(parent=parent, controller=controller)
        controller.frames["SystemPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("SystemPage")
    # ------------------------------------------------------------------------------------------------------------------
    # This method is called when the user wants to load input data from a file.
    # ------------------------------------------------------------------------------------------------------------------
    def loadFromData(self, parent, controller):
        defaultPath = os.path.dirname(os.path.abspath(__file__)) + "/Output/"
        filePath = filedialog.askopenfilename(initialdir=defaultPath, title="Select file", filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
        file = open(filePath, "r")
        words = []
        type = None
        nameFP = None
        massFP = None
        nameSP = None
        massSP = None
        distance = None
        for lines in file:
            line = lines.split()
            if "TYPE" in lines:
                type = line[2]
            if "NAME FIRST PRIMARY" in lines:
                nameFP = line[4]
            if "MASS FIRST PRIMARY" in lines:
                massFP = float(line[4])
            if "NAME SECOND PRIMARY" in lines:
                nameSP = line[4]
            if "MASS SECOND PRIMARY" in lines:
                massSP = float(line[4])
            if "PRIMARY DISTANCE" in lines:
                distance = float(line[3])
            if "LAGRANGIAN" in lines:
                lagrangian = line[2]
            if "ORBIT NUMBER" in lines:
                orbitNumber = int(line[3])
            if "ORBIT DISTANCE" in lines:
                orbitDistance = float(line[3])
            # if "DATA_START" in lines:
            #     for i in range(2):
            #         line = file.readline()
            #     words = line.split()

            if "DATA_START" in lines:
                for i in range(2):
                    line = file.readline()
                words = line.split()
                x0 = np.array([float(words[2]), 0, float(words[3]), 0, float(words[4]), 0])
                for i in range(3):
                    words.insert(2 * i + 3, 0)
                data = words
                while "DATA_STOP" not in line:
                    words = line.split()
                    for i in range(3):
                        words.insert(2 * i + 3, 0)
                    data = np.vstack([data, words])
                    line = file.readline()
                data = data.astype(np.float)
        file.close()
        dynamicalSystem = System(nameFP, massFP, nameSP, massSP, distance)

        if type == "ORBIT":
            orbit = Orbit(x0, "x", dynamicalSystem)
            controller.frames["OrbitPage"] = OrbitPage.OrbitPage(parent=parent, controller=controller, orbit=orbit, dynamicalSystem=dynamicalSystem)
            controller.frames["OrbitPage"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("OrbitPage")
        elif type == "FAMILY":
            if data[0,4] > data[-1,4]:
                direction = "Northern"
            else:
                direction = "Southern"
            orbitFamily = OrbitFamily(x0, direction, lagrangian, orbitDistance, dynamicalSystem, familyData=data)
            controller.frames["OrbitFamilyPage"] = OrbitFamilyPage.OrbitFamilyPage(parent=parent, controller=controller, orbitFamily=orbitFamily, dynamicalSystem=dynamicalSystem)
            controller.frames["OrbitFamilyPage"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("OrbitFamilyPage")
    # ------------------------------------------------------------------------------------------------------------------
    # This method terminates the program.
    # ------------------------------------------------------------------------------------------------------------------
    def client_exit(self):
        exit()
