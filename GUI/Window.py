"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
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
# The class Window contains methods for changing the window size.
# ----------------------------------------------------------------------------------------------------------------------
class Window():
    # ------------------------------------------------------------------------------------------------------------------
    # This method centers the window and sets its size.
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def centerWindow(window, width, height):
        # get screen width and height
        screen_width = window.winfo_screenwidth()
        screen_height = window.winfo_screenheight()
        # calculate position x and y coordinates
        x = (screen_width/2) - (width/2)
        y = (screen_height/2) - (height/2)
        window.geometry('%dx%d+%d+%d' % (width, height, x, y))
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets the window to full screen.
    # ------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def fullScreen(window):
        # get screen width and height
        screen_width = window.winfo_screenwidth()
        screen_height = window.winfo_screenheight()
        Window.centerWindow(window, screen_width, screen_height)
