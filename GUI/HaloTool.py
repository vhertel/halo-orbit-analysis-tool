"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import Window, StartPage

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
# The class HaloTool is the main class of the GUI.
# ----------------------------------------------------------------------------------------------------------------------
class HaloTool(Tk):
    # ------------------------------------------------------------------------------------------------------------------
    # The start page is called during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        # set GUI title
        Tk.wm_title(self, "HALO TOOL")
        # center window and set GUI size
        Window.Window.centerWindow(self, 700, 600)
        container = Frame(self)
        container.pack(side="left", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}
        self.frames["StartPage"] = StartPage.StartPage(parent=container, controller=self)
        self.frames["StartPage"].grid(row=0, column=0, sticky="nsew")
        self.setMenu()
        self.show_frame("StartPage")
    # ------------------------------------------------------------------------------------------------------------------
    # This method contains the structure of the menu.
    # ------------------------------------------------------------------------------------------------------------------
    def setMenu(self):
        menubar = Menu(self)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="New", command=lambda: HaloTool.donothing(self))
        filemenu.add_command(label="Open", command=lambda: HaloTool.donothing(self))
        filemenu.add_command(label="Save", command=lambda: HaloTool.donothing(self))
        filemenu.add_command(label="Save as...", command=lambda: HaloTool.donothing(self))
        filemenu.add_command(label="Close", command=lambda: HaloTool.donothing(self))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        editmenu = Menu(menubar, tearoff=0)
        editmenu.add_command(label="Undo", command=lambda: HaloTool.donothing(self))
        editmenu.add_separator()
        editmenu.add_command(label="Cut", command=lambda: HaloTool.donothing(self))
        editmenu.add_command(label="Copy", command=lambda: HaloTool.donothing(self))
        editmenu.add_command(label="Paste", command=lambda: HaloTool.donothing(self))
        editmenu.add_command(label="Delete", command=lambda: HaloTool.donothing(self))
        editmenu.add_command(label="Select All", command=lambda: HaloTool.donothing(self))
        menubar.add_cascade(label="Edit", menu=editmenu)
        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="Help Index", command=lambda: HaloTool.donothing(self))
        helpmenu.add_command(label="About...", command=lambda: HaloTool.donothing(self))
        menubar.add_cascade(label="Help", menu=helpmenu)
        self.config(menu=menubar)
    # ------------------------------------------------------------------------------------------------------------------
    # Dummy Method
    # ------------------------------------------------------------------------------------------------------------------
    def donothing(self):
        print("DO NOTHING")
    # ------------------------------------------------------------------------------------------------------------------
    # This method shows an other frame.
    # ------------------------------------------------------------------------------------------------------------------
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
