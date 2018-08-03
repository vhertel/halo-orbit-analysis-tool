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
# This Orbit main page contains all functionalities for single orbits.
# ----------------------------------------------------------------------------------------------------------------------
class OrbitPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # During initialization, the structure of the GUI is loaded, a three-dimensional view and all three projections
    # of the orbit are displayed.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller, orbit, dynamicalSystem):
        Frame.__init__(self, parent)
        frameLeft = Frame(self)
        frameLeft.pack(fill=BOTH, expand=True, side=LEFT, padx=10, pady=10)
        frameRight = Frame(self)
        frameRight.pack(fill=BOTH, expand=False, side=RIGHT, padx=10, pady=10)
        # orbit attributes
        orbitAttributes = LabelFrame(frameRight, text=" ORBIT ATTRIBUTES ")
        orbitAttributes.pack(fill=BOTH, expand=True, pady=(0, 7))
        orbitAttributes.columnconfigure(0, weight=1)
        # initial state
        Label(orbitAttributes, text="Initial State:", width=16, anchor=W).grid(row=0, column=0, padx=5, sticky='W', ipadx=5)
        self.state0Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[0]))
        self.state0Label.grid(row=0, column=1, columnspan=2, padx=5, sticky='E')
        self.state1Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[1]))
        self.state1Label.grid(row=1, column=1, columnspan=2, padx=5, sticky='E')
        self.state2Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[2]))
        self.state2Label.grid(row=2, column=1, columnspan=2, padx=5, sticky='E')
        self.state3Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[3]))
        self.state3Label.grid(row=3, column=1, columnspan=2, padx=5, sticky='E')
        self.state4Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[4]))
        self.state4Label.grid(row=4, column=1, columnspan=2, padx=5, sticky='E')
        self.state5Label = Label(orbitAttributes, text="%.10f" % (orbit.x0[5]))
        self.state5Label.grid(row=5, column=1, columnspan=2, padx=5, pady=(0,15), sticky='E')
        # period
        Label(orbitAttributes, text="Period:").grid(row=7, column=0, padx=5, sticky='W')
        self.periodLabel = Label(orbitAttributes, text="%.10f" % orbit.period)
        self.periodLabel.grid(row=7, column=1, columnspan=2, padx=5, sticky='E')
        # jacobi constant
        Label(orbitAttributes, text="Jacobi Constant:").grid(row=8, column=0, padx=5, sticky='W')
        self.jacobiLabel = Label(orbitAttributes, text=orbit.jacobi.round(4))
        self.jacobiLabel.grid(row=8, column=1, columnspan=2, padx=5, sticky='E')
        # stability index
        if orbit.stability is None:
            orbit.getStability()
        Label(orbitAttributes, text="Stability Index:").grid(row=9, column=0, padx=5, sticky='W')
        self.stabilityLabel = Label(orbitAttributes, text=orbit.stability.round(2))
        self.stabilityLabel.grid(row=9, column=1, columnspan=2, padx=5, sticky='E')
        # NRHO
        if orbit.NRHO:
            text = "Yes"
        else:
            text = "No"
        Label(orbitAttributes, text="NRHO:").grid(row=10, column=0, padx=5, sticky='W')
        self.NRHOLabel = Label(orbitAttributes, text=text)
        self.NRHOLabel.grid(row=10, column=1, columnspan=2, padx=5, sticky='E')
        # lagrangian
        Label(orbitAttributes, text="Lagrangian:").grid(row=11, column=0, padx=5, sticky='W')
        Label(orbitAttributes, text=orbit.lagrangian).grid(row=11, column=1, columnspan=2, padx=5, sticky='E')

        Label(orbitAttributes, text="Unit:", width=16, anchor=W, font=("TkDefaultFont", 11)).grid(row=12, column=0, padx=5, sticky='W', ipadx=5)
        self.unit = IntVar(value=0)
        Radiobutton(orbitAttributes, text="Normalized", variable=self.unit, value=0, command=lambda: OrbitPage.updateAttributes(self, orbit, dynamicalSystem, 0), font=("TkDefaultFont", 11)).grid(row=12, column=1, sticky='W', padx=5)
        Radiobutton(orbitAttributes, text="Physical", variable=self.unit, value=1, command=lambda: OrbitPage.updateAttributes(self, orbit, dynamicalSystem, 1), font=("TkDefaultFont", 11)).grid(row=12, column=2, sticky='W', padx=5)

        # system attributes
        systemAttributes = LabelFrame(frameRight, text=" SYSTEM ATTRIBUTES ")
        systemAttributes.pack(fill=BOTH, expand=True, pady=(0, 7))
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
        plotOptions.pack(fill=BOTH, expand=True, pady=(0, 7))
        Label(plotOptions, text="Manifold Direction:").grid(row=2, column=0, columnspan=4, padx=5, sticky='W')
        Label(plotOptions, text=dynamicalSystem.nameFP).grid(row=2, column=4, padx=5, pady=2)
        Label(plotOptions, text=dynamicalSystem.nameSP).grid(row=2, column=5, padx=5, pady=2)
        Label(plotOptions, text="Duration Factor:").grid(row=3, column=0, columnspan=3, padx=5, sticky='W')
        self.durationFactorEntryFP = Entry(plotOptions, width=5, justify=CENTER)
        self.durationFactorEntryFP.grid(row=3, column=2, columnspan=3, padx=5, sticky='E')
        self.durationFactorEntrySP = Entry(plotOptions, width=5, justify=CENTER)
        self.durationFactorEntrySP.grid(row=3, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultDurationFactorFP = 1.2
        self.defaultDurationFactorSP = 0.8
        self.durationFactorEntryFP.insert(END, self.defaultDurationFactorFP)
        self.durationFactorEntryFP.config(state=DISABLED)
        self.durationFactorEntrySP.insert(END, self.defaultDurationFactorSP)
        self.durationFactorEntrySP.config(state=DISABLED)
        Label(plotOptions, text="Number of Manifolds:").grid(row=4, column=0, columnspan=3, padx=5, sticky='W')
        self.numberOfManifoldEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.numberOfManifoldEntry.grid(row=4, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultNumberOfManifolds = 30
        self.numberOfManifoldEntry.insert(END, self.defaultNumberOfManifolds)
        self.numberOfManifoldEntry.config(state=DISABLED)
        Label(plotOptions, text="Invariant Manifolds:").grid(row=0, column=0, columnspan=4, padx=5, sticky='W')
        self.stable = IntVar(value=0)
        self.unstable = IntVar(value=0)
        Checkbutton(plotOptions, text="Stable Manifold", variable=self.stable, onvalue=1, offvalue=0, command=self.activateManifolds).grid(row=0, column=4, columnspan=2, padx=5, sticky='W')
        Checkbutton(plotOptions, text="Unstable Manifold", variable=self.unstable, onvalue=1, offvalue=0, command=self.activateManifolds).grid(row=1, column=4, columnspan=2, padx=5, sticky='W')
        Label(plotOptions, text="Objects:").grid(row=5, column=0, columnspan=2, padx=5, sticky='W')
        self.nameFP = IntVar(value=0)
        self.nameSP = IntVar(value=0)
        self.L1 = IntVar(value=0)
        self.L2 = IntVar(value=0)
        Checkbutton(plotOptions, text="L1", variable=self.L1, onvalue=1, offvalue=0).grid(row=5, column=2, padx=5)
        Checkbutton(plotOptions, text="L2", variable=self.L2, onvalue=1, offvalue=0).grid(row=5, column=3, padx=5)
        Checkbutton(plotOptions, text=dynamicalSystem.nameFP, variable=self.nameFP, onvalue=1, offvalue=0).grid(row=5, column=4, padx=5)
        Checkbutton(plotOptions, text=dynamicalSystem.nameSP, variable=self.nameSP, onvalue=1, offvalue=0).grid(row=5, column=5, padx=5)
        # button for updating plots
        Button(plotOptions, text="Update Plots",
               command=lambda: OrbitPage.plot(self, plotFrame, orbit, dynamicalSystem)).grid(row=6, column=0, columnspan=3, padx=5, sticky='W')
        # button for saving data
        Button(plotOptions, text="Save to file",
               command=lambda: OrbitPage.saveButton(self, orbit, dynamicalSystem)).grid(row=6, column=3, columnspan=3, padx=5, sticky='E')
        # commands
        actions = LabelFrame(frameRight, text=" COMMANDS ")
        actions.pack(fill=BOTH, expand=True, pady=(0, 7))
        # button to search for closest NRHO
        getNRHOButton = Button(actions, text="Get closest NRHO", command=lambda: OrbitPage.getClosestNRHO(self, plotFrame, orbit, dynamicalSystem))
        getNRHOButton.pack(side=LEFT, padx=5)
        # button to plan transfer
        planTransferButton = Button(actions, text="Plan Transfer")
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
        self.plot(plotFrame, orbit, dynamicalSystem)
        Window.Window.fullScreen(controller)
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets the specification entry of the manifolds to Active or Passive.
    # ------------------------------------------------------------------------------------------------------------------
    def activateManifolds(self):
        self.numberOfManifoldEntry.config(highlightbackground='white')
        self.durationFactorEntryFP.config(highlightbackground='white')
        self.durationFactorEntrySP.config(highlightbackground='white')
        if self.stable.get() or self.unstable.get():
            self.numberOfManifoldEntry.config(state=NORMAL)
            self.durationFactorEntryFP.config(state=NORMAL)
            self.durationFactorEntrySP.config(state=NORMAL)
        else:
            self.numberOfManifoldEntry.config(state=DISABLED)
            self.durationFactorEntryFP.config(state=DISABLED)
            self.durationFactorEntrySP.config(state=DISABLED)
    # ------------------------------------------------------------------------------------------------------------------
    # This method calls a window for specifying the data to be saved.
    # ------------------------------------------------------------------------------------------------------------------
    def saveButton(self, orbit, dynamicalSystem):
        self.saveLevel = Toplevel()
        self.dataFile = IntVar(value=0)
        self.figure = IntVar(value=0)
        self.xz = IntVar(value=0)
        self.yz = IntVar(value=0)
        self.xy = IntVar(value=0)
        Checkbutton(self.saveLevel, text="Data File", variable=self.dataFile, onvalue=1, offvalue=0).grid(row=1, sticky='W', padx=5, pady=2)
        Checkbutton(self.saveLevel, text="3D Figure", variable=self.figure, onvalue=1, offvalue=0).grid( row=2, sticky='W', padx=5, pady=2)
        Checkbutton(self.saveLevel, text="xz-Projection", variable=self.xz, onvalue=1, offvalue=0).grid(row=3, sticky='W', padx=5, pady=2)
        Checkbutton(self.saveLevel, text="xyz-Projection", variable=self.yz, onvalue=1, offvalue=0).grid(row=4, sticky='W', padx=5, pady=2)
        Checkbutton(self.saveLevel, text="xy-Projection", variable=self.xy, onvalue=1, offvalue=0).grid(row=5, sticky='W', padx=5, pady=2)
        Button(self.saveLevel, text="Confirm", command=lambda: OrbitPage.saveData(self, orbit, dynamicalSystem)).grid(row=19, column=0, sticky='W')
    # ------------------------------------------------------------------------------------------------------------------
    # This method saves the appropriate data.
    # ------------------------------------------------------------------------------------------------------------------
    def saveData(self, orbit, dynamicalSystem):
        timeStamp = time.strftime("%Y-%m-%dT%H.%M.%S")
        fileName = dynamicalSystem.nameFP + dynamicalSystem.nameSP + "_" + timeStamp
        defaultPath = os.path.dirname(os.path.abspath(__file__)) + "/Output/"
        browsedPath = filedialog.askdirectory(initialdir=defaultPath)
        path = browsedPath + "/" + fileName
        os.makedirs(path)
        if self.dataFile.get():
            output = open(path + "/data.txt", "w")
            # writes header data
            output.write("CREATION_DATE            =      " + timeStamp + "\n"
                                                                          "ORIGINATOR               =      ASTOS SOLUTIONS GMBH\n\n")
            # writes meta data
            output.write("META_START\n"
                         "TYPE                     =      ORBIT\n"
                         "NAME FIRST PRIMARY       =      %s\n"
                         "MASS FIRST PRIMARY       =      %e\n"
                         "NAME SECOND PRIMARY      =      %s\n"
                         "MASS SECOND PRIMARY      =      %e\n"
                         "PRIMARY DISTANCE         =      %e\n"
                         "MASS RATIO               =      %11.10f\n"
                         "LAGRANGIAN               =      %s\n"
                         "META_STOP\n\n" % (
                             dynamicalSystem.nameFP, dynamicalSystem.massFP, dynamicalSystem.nameSP,
                             dynamicalSystem.massSP,
                             dynamicalSystem.distance, dynamicalSystem.mu, orbit.lagrangian))
            # writes orbit data
            output.write("DATA_START\n")
            output.write("        JC           Period           x              z            dy/dt\n")
            output.write('{0:15.10f}'.format(orbit.jacobi))
            output.write('{0:15.10f}'.format(orbit.period))
            output.write('{0:15.10f}'.format(orbit.x0[0]))
            output.write('{0:15.10f}'.format(orbit.x0[2]))
            output.write('{0:15.10f}\n'.format(orbit.x0[4]))
            output.write("DATA_STOP")
            # closes file
            output.close()
        self.saveLevel.destroy()
    # ------------------------------------------------------------------------------------------------------------------
    # This method updates the orbit attributes.
    # ------------------------------------------------------------------------------------------------------------------
    def updateAttributes(self, orbit, dynamicalSystem, unit):
        if unit == 0:
            self.state0Label.config(text="%.10f" % (orbit.x0[0]))
            self.state1Label.config(text="%.10f" % (orbit.x0[1]))
            self.state2Label.config(text="%.10f" % (orbit.x0[2]))
            self.state3Label.config(text="%.10f" % (orbit.x0[3]))
            self.state4Label.config(text="%.10f" % (orbit.x0[4]))
            self.state5Label.config(text="%.10f" % (orbit.x0[5]))
            self.periodLabel.config(text="%.10f" % (orbit.period))
        elif unit == 1:
            vel = orbit.x0[3:6] * dynamicalSystem.distance / (np.sqrt(dynamicalSystem.distance ** 3 / (dynamicalSystem.G * (dynamicalSystem.massFP + dynamicalSystem.massSP))))
            self.state0Label.config(text="%.2f km" % (orbit.x0[0] * dynamicalSystem.distance/1.0e3))
            self.state1Label.config(text="%.2f km" % (orbit.x0[1] * dynamicalSystem.distance/1.0e3))
            self.state2Label.config(text="%.2f km" % (orbit.x0[2] * dynamicalSystem.distance/1.0e3))
            self.state3Label.config(text="%.2f m/s" % (vel[0]))
            self.state4Label.config(text="%.2f m/s" % (vel[1]))
            self.state5Label.config(text="%.2f m/s" % (vel[2]))
            period = (orbit.period * np.sqrt(dynamicalSystem.distance ** 3 /
                     (dynamicalSystem.G * (dynamicalSystem.massFP + dynamicalSystem.massSP)))) / (60 * 60 * 24)
            self.periodLabel.config(text="%.2f Days" % period)
        self.jacobiLabel.config(text=orbit.jacobi.round(4))
        self.stabilityLabel.config(text=orbit.stability.round(2))
        if orbit.NRHO:
            text = "Yes"
        else:
            text = "No"
        self.NRHOLabel.config(text=text)
    # ------------------------------------------------------------------------------------------------------------------
    # This method searches for the closest NRHO.
    # ------------------------------------------------------------------------------------------------------------------
    def getClosestNRHO(self, plotFrame, orbit, dynamicalSystem):
        self.status.config(text="Searching for closest NRHO ...")
        self.status.update()
        var = IntVar()
        if abs(orbit.x0[2]) < 1.0e-4:
            self.NRHOLevel = Toplevel()
            Radiobutton(self.NRHOLevel, text="Up", variable=var, value=0).grid(row=0, column=0)
            Radiobutton(self.NRHOLevel, text="Down", variable=var, value=1).grid(row=1, column=0)
            test = Button(self.NRHOLevel, text="Confirm", command=lambda:[self.NRHOLevel.destroy(), orbit.getClosestNRHO(var.get()), self.plot(plotFrame, orbit, dynamicalSystem), self.updateAttributes(orbit, dynamicalSystem, unit=0)])
            test.grid(row=3, column=0)
        else:
            orbit.getClosestNRHO()
            self.plot(plotFrame, orbit, dynamicalSystem)
            self.updateAttributes(orbit, dynamicalSystem, unit=0)
    # ------------------------------------------------------------------------------------------------------------------
    # This method creates and updates the plots of the orbit.
    # ------------------------------------------------------------------------------------------------------------------
    def plot(self, plotFrame, orbit, dynamicalSystem):
        # closes old figures
        plt.close('all')
        for widget in plotFrame.winfo_children():
            widget.destroy()
        mainFrame = Frame(plotFrame)
        mainFrame.pack(fill=BOTH, expand=True)
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
        # calculates position of L1 and L2
        lagrangianPosition = Utility.lagrangianPosition(dynamicalSystem.mu)




        if self.stable.get() or self.unstable.get():
            try:
                # stores input number of manifolds
                numberOfManifolds = int(self.numberOfManifoldEntry.get())
                # stores factor of integration
                durationFactorFP = float(self.durationFactorEntryFP.get())
                durationFactorSP = float(self.durationFactorEntrySP.get())
            except ValueError:
                self.status.config(text="Missing Input in Manifold Options")
                self.status.update()
                self.stable.set(0)
                self.unstable.set(0)
                if self.numberOfManifoldEntry.get() == "":
                    self.numberOfManifoldEntry.config(highlightbackground='red')
                if self.durationFactorEntryFP.get() == "":
                    self.durationFactorEntryFP.config(highlightbackground='red')
                if self.durationFactorEntrySP.get() == "":
                    self.durationFactorEntrySP.config(highlightbackground='red')
                self.numberOfManifoldEntry.config(state=DISABLED)
                self.durationFactorEntryFP.config(state=DISABLED)
                self.durationFactorEntrySP.config(state=DISABLED)
            else:
                # checks whether input of manifold direction, number or duration factor has been changed
                if self.defaultNumberOfManifolds == int(self.numberOfManifoldEntry.get()) and self.defaultDurationFactorFP == float(self.durationFactorEntryFP.get()) and self.defaultDurationFactorSP == float(self.durationFactorEntrySP.get()):
                    inputChange = False
                else:
                    inputChange = True
                # updates manifold number and duration factor
                self.defaultNumberOfManifolds = numberOfManifolds
                self.defaultDurationFactorFP = durationFactorFP
                self.defaultDurationFactorSP = durationFactorSP
                # checks if manifolds need to be recalculated
                if ((self.stable.get() or self.unstable.get()) and orbit.stableManifolds is None) or ((self.stable.get() or self.unstable.get()) and inputChange):
                    self.status.config(text="Calculating invariant Manifolds ...")
                    self.status.update()


                    # calculates initial states of (un)stable manifolds
                    orbit.invariantManifolds(numberOfPoints=numberOfManifolds, direction=0)
                    # declares arrays
                    self.xStableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    self.yStableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    self.zStableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    self.xUnstableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    self.yUnstableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    self.zUnstableManifoldFP = np.zeros((numberOfManifolds, int(durationFactorFP*50)))
                    # integrates all states of (un)stable manifolds and stores data
                    t = np.linspace(0, durationFactorFP * orbit.period, num=int(durationFactorFP*50))
                    for i in range(0, numberOfManifolds):
                        stableManifold = odeint(Utility.backwards, orbit.stableManifolds[i], t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
                        self.xStableManifoldFP[i, :] = stableManifold[:, 0] * dynamicalSystem.distance
                        self.yStableManifoldFP[i, :] = stableManifold[:, 1] * dynamicalSystem.distance
                        self.zStableManifoldFP[i, :] = stableManifold[:, 2] * dynamicalSystem.distance
                        unstableManifold = odeint(Utility.sysEquations, orbit.unstableManifolds[i], t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
                        self.xUnstableManifoldFP[i, :] = unstableManifold[:, 0] * dynamicalSystem.distance
                        self.yUnstableManifoldFP[i, :] = unstableManifold[:, 1] * dynamicalSystem.distance
                        self.zUnstableManifoldFP[i, :] = unstableManifold[:, 2] * dynamicalSystem.distance

                    orbit.invariantManifolds(numberOfPoints=numberOfManifolds, direction=1)
                    # declares arrays
                    self.xStableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    self.yStableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    self.zStableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    self.xUnstableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    self.yUnstableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    self.zUnstableManifoldSP = np.zeros((numberOfManifolds, int(durationFactorSP*50)))
                    # integrates all states of (un)stable manifolds and stores data
                    t = np.linspace(0, durationFactorSP * orbit.period, num=int(durationFactorSP*50))
                    for i in range(0, numberOfManifolds):
                        stableManifold = odeint(Utility.backwards, orbit.stableManifolds[i], t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
                        self.xStableManifoldSP[i, :] = stableManifold[:, 0] * dynamicalSystem.distance
                        self.yStableManifoldSP[i, :] = stableManifold[:, 1] * dynamicalSystem.distance
                        self.zStableManifoldSP[i, :] = stableManifold[:, 2] * dynamicalSystem.distance
                        unstableManifold = odeint(Utility.sysEquations, orbit.unstableManifolds[i], t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
                        self.xUnstableManifoldSP[i, :] = unstableManifold[:, 0] * dynamicalSystem.distance
                        self.yUnstableManifoldSP[i, :] = unstableManifold[:, 1] * dynamicalSystem.distance
                        self.zUnstableManifoldSP[i, :] = unstableManifold[:, 2] * dynamicalSystem.distance




                    self.status.config(text="")
                    self.status.update()




        self.status.config(text="Plotting Data ...")
        self.status.update()
        # integrates initial state of orbit and stores data
        t = np.linspace(0, orbit.period, num=300)
        orbitStates = odeint(Utility.sysEquations, orbit.x0, t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
        self.xOrbit = orbitStates[:, 0] * dynamicalSystem.distance
        self.yOrbit = orbitStates[:, 1] * dynamicalSystem.distance
        self.zOrbit = orbitStates[:, 2] * dynamicalSystem.distance

        # sets xz-projection
        self.fig1 = plt.figure(figsize=(3.5, 3.5))
        plt.axis("equal")
        #plt.grid(True)
        self.canvas = FigureCanvasTkAgg(self.fig1, master=framefig1)
        self.canvas._tkcanvas.pack(expand=False)
        plt.title("$xz$-Projection")
        plt.xlabel('[m]', fontsize=8)
        plt.ylabel('[m]', fontsize=8)
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xStableManifoldFP[i, :], self.zStableManifoldFP[i, :], color='green', linewidth=0.5)
                plt.plot(self.xStableManifoldSP[i, :], self.zStableManifoldSP[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xUnstableManifoldFP[i, :], self.zUnstableManifoldFP[i, :], color='red', linewidth=0.5)
                plt.plot(self.xUnstableManifoldSP[i, :], self.zUnstableManifoldSP[i, :], color='red', linewidth=0.5)
        plt.plot(self.xOrbit, self.zOrbit, color='black', linewidth=0.75)
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
        plt.xlabel('[m]', fontsize=8)
        plt.ylabel('[m]', fontsize=8)
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.yStableManifoldFP[i, :], self.zStableManifoldFP[i, :], color='green', linewidth=0.5)
                plt.plot(self.yStableManifoldSP[i, :], self.zStableManifoldSP[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.yUnstableManifoldFP[i, :], self.zUnstableManifoldFP[i, :], color='red', linewidth=0.5)
                plt.plot(self.yUnstableManifoldSP[i, :], self.zUnstableManifoldSP[i, :], color='red', linewidth=0.5)
        plt.plot(self.yOrbit, self.zOrbit, color='black', linewidth=0.75)
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
        plt.xlabel('[m]', fontsize=8)
        plt.ylabel('[m]', fontsize=8)
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xStableManifoldFP[i, :], self.yStableManifoldFP[i, :], color='green', linewidth=0.5)
                plt.plot(self.xStableManifoldSP[i, :], self.yStableManifoldSP[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xUnstableManifoldFP[i, :], self.yUnstableManifoldFP[i, :], color='red', linewidth=0.5)
                plt.plot(self.xUnstableManifoldSP[i, :], self.yUnstableManifoldSP[i, :], color='red', linewidth=0.5)
        plt.plot(self.xOrbit, self.yOrbit, color='black', linewidth=0.75)
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
        ax.set_xlabel('[m]', fontsize=8)
        ax.set_ylabel('[m]', fontsize=8)
        ax.set_zlabel('[m]', fontsize=8)
        if self.stable.get():
            for i in range(numberOfManifolds):
                ax.plot(self.xStableManifoldFP[i, :], self.yStableManifoldFP[i, :], self.zStableManifoldFP[i, :], color='green', linewidth=0.5)
                ax.plot(self.xStableManifoldSP[i, :], self.yStableManifoldSP[i, :], self.zStableManifoldSP[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                ax.plot(self.xUnstableManifoldFP[i, :], self.yUnstableManifoldFP[i, :], self.zUnstableManifoldFP[i, :], color='red', linewidth=0.5)
                ax.plot(self.xUnstableManifoldSP[i, :], self.yUnstableManifoldSP[i, :], self.zUnstableManifoldSP[i, :], color='red', linewidth=0.5)
        ax.plot(self.xOrbit, self.yOrbit, self.zOrbit, color='black', linewidth=0.75)
        if self.nameFP.get():
            ax.scatter((-dynamicalSystem.mu) * dynamicalSystem.distance, 0, 0, color='grey', s=3)
        if self.nameSP.get():
            ax.scatter((1 - dynamicalSystem.mu) * dynamicalSystem.distance, 0, 0, color='grey', s=3)
        if self.L1.get():
            ax.scatter(lagrangianPosition[0] * dynamicalSystem.distance, 0, 0, color='blue', s=1)
        if self.L2.get():
            ax.scatter(lagrangianPosition[1] * dynamicalSystem.distance, 0, 0, color='blue', s=1)
        Plot.setAxesEqual(ax)
        self.status.config(text="")
        self.status.update()
