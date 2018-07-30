"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import Window

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
# This OrbitFamily main page contains all functionalities for orbit families.
# ----------------------------------------------------------------------------------------------------------------------
class OrbitFamilyPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller, orbitFamily, dynamicalSystem):
        Frame.__init__(self, parent)
        frameLeft = Frame(self)
        frameLeft.pack(fill=BOTH, expand=True, side=LEFT, padx=10, pady=10)
        frameRight = Frame(self)
        frameRight.pack(fill=BOTH, expand=False, side=RIGHT, padx=10, pady=10)
        # orbit attributes
        orbitAttributes = LabelFrame(frameRight, text=" ORBIT ATTRIBUTES ")
        orbitAttributes.pack(fill=BOTH, expand=True, pady=(0, 10))
        orbitAttributes.columnconfigure(0, weight=1)

        # initial state
        Label(orbitAttributes, text="Calculated Orbits:").grid(row=0, column=0, columnspan=2, padx=5, sticky='W')
        self.calculatedOrbitNumber = Label(orbitAttributes, text=len(orbitFamily.familyData))
        self.calculatedOrbitNumber.grid(row=0, column=2, padx=5, sticky='E')
        Label(orbitAttributes, text="Displayed Orbits:").grid(row=1, column=0, columnspan=2, padx=5, sticky='W')
        self.displayedOrbitNumber = Label(orbitAttributes, text=len(orbitFamily.familyData))
        self.displayedOrbitNumber.grid(row=1, column=2, padx=5, sticky='E')
        Label(orbitAttributes, text="Lagrangian:").grid(row=2, column=0, columnspan=2, padx=5, sticky='W')
        self.lagrangianLabel = Label(orbitAttributes, text=orbitFamily.lagrangian)
        self.lagrangianLabel.grid(row=2, column=2, padx=5, sticky='E')
        Label(orbitAttributes, text="Distance calculated Orbits:").grid(row=3, column=0, columnspan=2, padx=5, sticky='W')
        self.calDistanceLabel = Label(orbitAttributes, text=orbitFamily.orbitDistance)
        self.calDistanceLabel.grid(row=3, column=2, padx=5, sticky='E')
        Label(orbitAttributes, text="Distance displayed Orbits:").grid(row=4, column=0, columnspan=2, padx=5, sticky='W')
        self.disDistanceLabel = Label(orbitAttributes, text=orbitFamily.orbitDistance)
        self.disDistanceLabel.grid(row=4, column=2, padx=5, sticky='E')
        Label(orbitAttributes, text="Unit:", anchor=W, font=("TkDefaultFont", 11)).grid(row=5, column=0, padx=5, sticky='W', ipadx=5)
        self.unit = IntVar(value=0)
        Radiobutton(orbitAttributes, text="Normalized", variable=self.unit, value=0, command=lambda: OrbitFamilyPage.updateAttributes(self, orbitFamily, dynamicalSystem, 0), font=("TkDefaultFont", 11)).grid(row=5, column=1, sticky='W', padx=5)
        Radiobutton(orbitAttributes, text="Physical", variable=self.unit, value=1, command=lambda: OrbitFamilyPage.updateAttributes(self, orbitFamily, dynamicalSystem, 1), font=("TkDefaultFont", 11)).grid(row=5, column=2, sticky='W', padx=5)


        # system attributes
        systemAttributes = LabelFrame(frameRight, text=" SYSTEM ATTRIBUTES ")
        systemAttributes.pack(fill=BOTH, expand=True, pady=(0, 10))
        systemAttributes.columnconfigure(0, weight=1)
        # first primary
        Label(systemAttributes, text="First Primary:", width=17, anchor=W).grid(row=0, column=0, padx=5, sticky='W')
        Label(systemAttributes, text=dynamicalSystem.nameFP).grid(row=0, column=1, padx=5, sticky='E')
        Label(systemAttributes, text="Mass of " + dynamicalSystem.nameFP + ":").grid(row=1, column=0, padx=5, sticky='W')
        Label(systemAttributes, text="%e kg" % (dynamicalSystem.massFP)).grid(row=1, column=1, padx=5, sticky='E')
        # second primary
        Label(systemAttributes, text="Second Primary:").grid(row=2, column=0, padx=5, sticky='W')
        Label(systemAttributes, text=dynamicalSystem.nameSP).grid(row=2, column=1, padx=5, sticky='E')
        Label(systemAttributes, text="Mass of " + dynamicalSystem.nameSP + ":").grid(row=3, column=0, padx=5, sticky='W')
        Label(systemAttributes, text="%e kg" % (dynamicalSystem.massSP)).grid(row=3, column=1, padx=5, sticky='E')
        # distance between primaries
        Label(systemAttributes, text="Distance:").grid(row=4, column=0, padx=5, sticky='W')
        Label(systemAttributes, text="%.0f km" % (dynamicalSystem.distance/1.0e3)).grid(row=4, column=1, padx=5, sticky='E')

        # plot options
        plotOptions = LabelFrame(frameRight, text=" PLOT CONFIGURATION ")
        plotOptions.pack(fill=BOTH, expand=True, pady=(0, 10))

        self.actualPlot = IntVar(value=0)
        Label(plotOptions, text="Display:").grid(row=0, column=0, padx=5, sticky='W')
        Radiobutton(plotOptions, text="Orbits", variable=self.actualPlot, value=0, command=lambda: OrbitFamilyPage.activateOrbitProperties(self)).grid(row=0, column=1, columnspan=2, sticky='W', padx=5, pady=2)
        Radiobutton(plotOptions, text="Period", variable=self.actualPlot, value=1, command=lambda: OrbitFamilyPage.activateOrbitProperties(self)).grid(row=0, column=3, columnspan=2, sticky='W', padx=5, pady=2)
        Radiobutton(plotOptions, text="Jacobi", variable=self.actualPlot, value=2, command=lambda: OrbitFamilyPage.activateOrbitProperties(self)).grid(row=1, column=1, columnspan=2, sticky='W', padx=5, pady=2)
        Radiobutton(plotOptions, text="Stability", variable=self.actualPlot, value=3, command=lambda: OrbitFamilyPage.activateOrbitProperties(self)).grid(row=1, column=3, columnspan=2, sticky='W', padx=5, pady=2)

        Label(plotOptions, text="Number of Orbits:").grid(row=2, column=0, columnspan=3, padx=5, sticky='W')
        self.numberofOrbitsEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.numberofOrbitsEntry.grid(row=2, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultNumberOfOrbits = len(orbitFamily.familyData)
        self.numberofOrbitsEntry.insert(END, self.defaultNumberOfOrbits)

        Label(plotOptions, text="Objects:").grid(row=3, column=0, columnspan=2, padx=5, sticky='W')
        self.nameFP = IntVar(value=0)
        self.nameSP = IntVar(value=0)
        self.L1 = IntVar(value=0)
        self.L2 = IntVar(value=0)
        self.L1CB = Checkbutton(plotOptions, text="L1", variable=self.L1, onvalue=1, offvalue=0)
        self.L1CB.grid(row=3, column=2, padx=5)
        self.L2CB = Checkbutton(plotOptions, text="L2", variable=self.L2, onvalue=1, offvalue=0)
        self.L2CB.grid(row=3, column=3, padx=5)
        self.nameFPCB = Checkbutton(plotOptions, text=dynamicalSystem.nameFP, variable=self.nameFP, onvalue=1, offvalue=0)
        self.nameFPCB.grid(row=3, column=4, padx=5)
        self.nameSPCB = Checkbutton(plotOptions, text=dynamicalSystem.nameSP, variable=self.nameSP, onvalue=1, offvalue=0)
        self.nameSPCB.grid(row=3, column=5, padx=5)


        # button for updating plots
        Button(plotOptions, text="Update Plots", command=lambda: OrbitFamilyPage.plot(self, plotFrame, orbitFamily, dynamicalSystem)).grid(row=5, column=0, columnspan=3, padx=5, sticky='W')
        # button for saving data
        Button(plotOptions, text="Save to file").grid(row=5, column=3, columnspan=3, padx=5, sticky='E')
        # commands
        actions = LabelFrame(frameRight, text=" COMMANDS ")
        actions.pack(fill=BOTH, expand=True, pady=(0, 10))
        # button to search for closest NRHO
        getNRHOButton = Button(actions, text="Get HALO Family")
        getNRHOButton.pack(side=LEFT, padx=5)
        # button to plan transfer
        planTransferButton = Button(actions, text="Get NRHO Family")
        planTransferButton.pack(side=RIGHT, padx=5)
        # status frame
        statusFrame = Frame(frameRight)
        statusFrame.pack(fill=X, side=BOTTOM)
        statusLabel = Label(statusFrame, text="Status:")
        statusLabel.config(font=("TkDefaultFont", 11))
        statusLabel.pack(side=LEFT, padx=(5,2))
        self.status = Label(statusFrame)
        self.status.config(font=("TkDefaultFont", 11))
        self.status.pack(side=LEFT, fill=X)
        Button(statusFrame, text="Home", command=lambda: [Window.Window.centerWindow(controller, 700, 600), controller.show_frame("StartPage"), self.destroy()]).pack(side=RIGHT, padx=5)
        # plots
        plotFrame = LabelFrame(frameLeft, text=" ORBIT PLOTS ")
        plotFrame.pack(fill=BOTH, expand=True)
        Window.Window.fullScreen(controller)
        self.plot(plotFrame, orbitFamily, dynamicalSystem)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets the specification entry of the manifolds to Active or Passive.
    # ------------------------------------------------------------------------------------------------------------------
    def activateOrbitProperties(self):
        self.numberofOrbitsEntry.config(highlightbackground='white')
        if self.actualPlot.get() == 0:
            self.L1CB.config(state=NORMAL)
            self.L2CB.config(state=NORMAL)
            self.nameFPCB.config(state=NORMAL)
            self.nameSPCB.config(state=NORMAL)
            self.numberofOrbitsEntry.config(state=NORMAL)
        else:
            self.L1CB.config(state=DISABLED)
            self.L2CB.config(state=DISABLED)
            self.nameFPCB.config(state=DISABLED)
            self.nameSPCB.config(state=DISABLED)
            self.numberofOrbitsEntry.config(state=DISABLED)
    # ------------------------------------------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------------------------------------------
    def saveButton(self, dynamicalSystem):
        pass
    # ------------------------------------------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------------------------------------------
    def saveData(self, orbit, dynamicalSystem):
        pass
    # ------------------------------------------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------------------------------------------
    def updateAttributes(self, orbitFamily, dynamicalSystem, unit):
        numDisplayedOrbits = int(self.numberofOrbitsEntry.get())
        disDistance = len(orbitFamily.familyData)/numDisplayedOrbits * orbitFamily.orbitDistance
        if unit == 0:
            self.calculatedOrbitNumber.config(text=len(orbitFamily.familyData))
            self.displayedOrbitNumber.config(text=numDisplayedOrbits)
            self.calDistanceLabel.config(text=orbitFamily.orbitDistance)
            self.disDistanceLabel.config(text=disDistance)
        elif unit == 1:
            self.calculatedOrbitNumber.config(text=len(orbitFamily.familyData))
            self.displayedOrbitNumber.config(text=numDisplayedOrbits)
            self.calDistanceLabel.config(text="%.2f km" % (orbitFamily.orbitDistance * dynamicalSystem.distance/1.0e3))
            self.disDistanceLabel.config(text="%.2f km" % (disDistance * dynamicalSystem.distance/1.0e3))
    # ------------------------------------------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------------------------------------------
    def plot(self, plotFrame, orbitFamily, dynamicalSystem):
        self.status.config(text="Plotting Data ...")
        self.status.update()
        # closes old figures
        plt.close('all')
        for widget in plotFrame.winfo_children():
            widget.destroy()
        mainFrame = Frame(plotFrame)
        mainFrame.pack(fill=BOTH, expand=True)




        if self.actualPlot.get() == 0:

            if int(self.numberofOrbitsEntry.get()) != self.defaultNumberOfOrbits:
                step = len(orbitFamily.familyData) / (int(self.numberofOrbitsEntry.get()) - 1)
                reducedData = orbitFamily.familyData[0, :]
                for i in range(int(self.numberofOrbitsEntry.get()) - 2):
                    reducedData = np.vstack([reducedData, orbitFamily.familyData[round(step * (i + 1)), :]])
                reducedData = np.vstack([reducedData, orbitFamily.familyData[-1, :]])
                data = reducedData
            else:
                data = orbitFamily.familyData

            upperFrame = Frame(mainFrame)
            upperFrame.pack(side=TOP, fill=BOTH, expand=True)
            downFrame = Frame(mainFrame)
            downFrame.pack(side=BOTTOM, fill=BOTH, expand=True)

            framefig1 = Frame(upperFrame)
            framefig1.pack(side=LEFT, fill=BOTH, expand=True)
            framefig2 = Frame(upperFrame)
            framefig2.pack(side=RIGHT, fill=BOTH, expand=True)
            framefig3 = Frame(downFrame)
            framefig3.pack(side=LEFT, fill=BOTH, expand=True)
            framefig4 = Frame(downFrame)
            framefig4.pack(side=RIGHT, fill=BOTH, expand=True)


            lagrangianPosition = Utility.lagrangianPosition(dynamicalSystem.mu)
            jacobiMax = max(orbitFamily.familyData[:, 0])
            jacobiMin = min(orbitFamily.familyData[:, 0])
            normalize = mpl.colors.Normalize(vmin=jacobiMin, vmax=jacobiMax)
            colormap = plt.cm.hsv_r
            scalarmappaple = plt.cm.ScalarMappable(norm=normalize, cmap=colormap)
            scalarmappaple.set_array(5)

            # sets xz-projection
            self.fig1 = plt.figure(figsize=(3.5, 3.5))
            plt.axis("equal")
            #plt.grid(True)
            self.canvas = FigureCanvasTkAgg(self.fig1, master=framefig1)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("$xz$-Projection")
            plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
            if self.nameFP.get():
                plt.scatter((-dynamicalSystem.mu) * dynamicalSystem.distance, 0, color='grey', s=3)
            if self.nameSP.get():
                plt.scatter((1 - dynamicalSystem.mu) * dynamicalSystem.distance, 0, color='grey', s=3)
            if self.L1.get():
                plt.scatter(lagrangianPosition[0] * dynamicalSystem.distance, 0, color='blue', s=1)
            if self.L2.get():
                plt.scatter(lagrangianPosition[1] * dynamicalSystem.distance, 0, color='blue', s=1)

            # sets yz-projection
            self.fig2 = plt.figure(figsize=(3.5, 3.5))
            plt.axis("equal")
            #plt.grid(True)
            self.canvas = FigureCanvasTkAgg(self.fig2, master=framefig2)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("$yz$-Projection")
            plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
            if self.nameFP.get():
                plt.scatter(0, 0, color='grey', s=3)
            if self.nameSP.get():
                plt.scatter(0, 0, color='grey', s=3)
            if self.L1.get():
                plt.scatter(0, 0, color='blue', s=1)
            if self.L2.get():
                plt.scatter(0, 0, color='blue', s=1)

            # sets xy-projection
            self.fig3 = plt.figure(figsize=(3.5, 3.5))
            plt.axis("equal")
            #plt.grid(True)
            self.canvas = FigureCanvasTkAgg(self.fig3, master=framefig3)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("$xy$-Projection")
            plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
            if self.nameFP.get():
                plt.scatter((-dynamicalSystem.mu) * dynamicalSystem.distance, 0, color='grey', s=3)
            if self.nameSP.get():
                plt.scatter((1 - dynamicalSystem.mu) * dynamicalSystem.distance, 0, color='grey', s=3)
            if self.L1.get():
                plt.scatter(lagrangianPosition[0] * dynamicalSystem.distance, 0, color='blue', s=1)
            if self.L2.get():
                plt.scatter(lagrangianPosition[1] * dynamicalSystem.distance, 0, color='blue', s=1)

            # sets 3d figure
            self.fig4 = plt.figure(figsize=(3.5, 3.5))
            self.canvas = FigureCanvasTkAgg(self.fig4, master=framefig4)
            self.canvas._tkcanvas.pack(expand=False)
            ax = Axes3D(self.fig4)
            ax.tick_params(axis='both', which='both', labelsize=7, direction='out')
            if self.nameFP.get():
                ax.scatter((-dynamicalSystem.mu) * dynamicalSystem.distance, 0, 0, color='grey', s=3)
            if self.nameSP.get():
                ax.scatter((1 - dynamicalSystem.mu) * dynamicalSystem.distance, 0, 0, color='grey', s=3)
            if self.L1.get():
                ax.scatter(lagrangianPosition[0] * dynamicalSystem.distance, 0, 0, color='blue', s=1)
            if self.L2.get():
                ax.scatter(lagrangianPosition[1] * dynamicalSystem.distance, 0, 0, color='blue', s=1)

            for i in range(len(data)):
                norm = (data[i, 0] - jacobiMin) / (jacobiMax - jacobiMin)
                color = plt.cm.hsv_r(norm)
                t = np.linspace(0, data[i, 1], num=300)
                halo = odeint(Utility.sysEquations, data[i, 2:8], t, args=(dynamicalSystem.mu,))#, rtol=2.5e-8, atol=1e-12)
                x = halo[:, 0] * dynamicalSystem.distance
                y = halo[:, 1] * dynamicalSystem.distance
                z = halo[:, 2] * dynamicalSystem.distance
                plt.figure(1)
                plt.plot(x, z, color=color, linewidth=0.25)
                plt.figure(2)
                plt.plot(y, z, color=color, linewidth=0.25)
                plt.figure(3)
                plt.plot(x, y, color=color, linewidth=0.25)
                ax.plot(x, y, z, color=color, linewidth=0.25)
            Plot.setAxesEqual(ax)

        if self.actualPlot.get() == 1:
            self.fig = plt.figure(figsize=(8, 8))
            self.canvas = FigureCanvasTkAgg(self.fig, master=mainFrame)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("%s - %s %s" % (dynamicalSystem.nameFP, dynamicalSystem.nameSP, orbitFamily.lagrangian))
            plt.xlabel("$x$ [m]")
            plt.ylabel("Period [Days]")
            for element in orbitFamily.familyData:
                plt.scatter(element[2] * dynamicalSystem.distance, (element[1] * np.sqrt(dynamicalSystem.distance ** 3 / (dynamicalSystem.G * (dynamicalSystem.massFP + dynamicalSystem.massSP)))) / (60 * 60 * 24), s=0.5, color='blue')

        if self.actualPlot.get() == 2:
            self.fig = plt.figure(figsize=(8, 8))
            self.canvas = FigureCanvasTkAgg(self.fig, master=mainFrame)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("%s - %s %s" % (dynamicalSystem.nameFP, dynamicalSystem.nameSP, orbitFamily.lagrangian))
            plt.xlabel("$x$ [m]")
            plt.ylabel("Jacobi Constant")
            for element in orbitFamily.familyData:
                plt.scatter(element[2] * dynamicalSystem.distance, element[0], s=0.5, color='blue')

        if self.actualPlot.get() == 3:
            self.fig = plt.figure(figsize=(8, 8))
            self.canvas = FigureCanvasTkAgg(self.fig, master=mainFrame)
            self.canvas._tkcanvas.pack(expand=False)
            plt.title("%s - %s %s" % (dynamicalSystem.nameFP, dynamicalSystem.nameSP, orbitFamily.lagrangian))
            plt.xlabel("$x$ [m]")
            plt.ylabel("Stability Index")
            for element in orbitFamily.familyData:
                monodromy = Utility.stm(element[2:8], element[1], dynamicalSystem.mu)
                eigenvalues, eigenvectors = np.linalg.eig(monodromy)
                maximum = max(abs(eigenvalues))
                stability = 1 / 2 * (maximum + 1 / maximum)
                plt.scatter(element[2] * dynamicalSystem.distance, stability, s=0.5, color='blue')


        self.status.config(text="")
        self.status.update()
        self.updateAttributes(orbitFamily, dynamicalSystem, 0)
