"""
File    : main.py
Author  : Victor Hertel
Date    : 20.07.2018

Graphical User Interface of the Application
"""



# imports
from GUI import OrbitPage

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
# This class contains a number of configuration options for an orbit.
# ----------------------------------------------------------------------------------------------------------------------
class OrbitInputPage(Frame):
    # ------------------------------------------------------------------------------------------------------------------
    # The GUI structure is loaded during initialization.
    # ------------------------------------------------------------------------------------------------------------------
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        # initial state
        self.guessState = np.zeros(6)
        stateFrame = LabelFrame(self, text=" INITIAL STATE ")
        stateFrame.grid(row=0, columnspan=2, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        Label(stateFrame, text="Input of Initial Guess:").grid(row=0, column=0, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="x").grid(row=1, column=0, sticky='W', padx=5, pady=2)
        self.xEntry = Entry(stateFrame)
        self.xEntry.grid(row=1, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="y").grid(row=2, column=0, sticky='W', padx=5, pady=2)
        self.yEntry = Entry(stateFrame)
        self.yEntry.grid(row=2, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="z").grid(row=3, column=0, sticky='W', padx=5, pady=2)
        self.zEntry = Entry(stateFrame)
        self.zEntry.grid(row=3, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="dx/dt").grid(row=4, column=0, sticky='W', padx=5, pady=2)
        self.dxEntry = Entry(stateFrame)
        self.dxEntry.grid(row=4, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="dy/dt").grid(row=5, column=0, sticky='W', padx=5, pady=2)
        self.dyEntry = Entry(stateFrame)
        self.dyEntry.grid(row=5, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="dz/dt").grid(row=6, column=0, sticky='W', padx=5, pady=2)
        self.dzEntry = Entry(stateFrame)
        self.dzEntry.grid(row=6, column=1, sticky='W', padx=5, pady=2)
        Label(stateFrame, text="Fixed value:").grid(row=7, column=0, columnspan=2, sticky='W', padx=5, pady=(15, 0))
        self.fixedValue = ('x', 'z', 'dy/dt', 'Period', 'None')
        self.fixedValueCB = ttk.Combobox(stateFrame, state='readonly', values=self.fixedValue, justify='center', width=10)
        self.fixedValueCB.set('x')
        self.fixedValueCB.grid(row=7, column=1, sticky='E', padx=5, pady=(15, 0))

        # initial guess generation
        guessGeneration = LabelFrame(self, text=" INITIAL GUESS GENERATION ")
        guessGeneration.grid(row=0, column=2, columnspan=2, sticky='NS', padx=5, pady=5, ipadx=5, ipady=5)

        guessOption = IntVar()
        Radiobutton(guessGeneration, text="L1 Planar Lyapunov Orbit", variable=guessOption, value=1, command=lambda: [OrbitInputPage.customInput(self, "DISABLED"), OrbitInputPage.guessGenerator(self, dynamicalSystem, lyapunov=True, lagrangian="L1")]).grid(row=0, columnspan=2, sticky='W', padx=5, pady=2)
        Radiobutton(guessGeneration, text="L2 Planar Lyapunov Orbit", variable=guessOption, value=2, command=lambda: [OrbitInputPage.customInput(self, "DISABLED"), OrbitInputPage.guessGenerator(self, dynamicalSystem, lyapunov=True, lagrangian="L2")]).grid(row=1, columnspan=2, sticky='W', padx=5, pady=2)
        Radiobutton(guessGeneration, text="Custom", variable=guessOption, value=3, command=lambda: OrbitInputPage.customInput(self, "NORMAL")).grid(row=2, columnspan=2, sticky='W', padx=5, pady=(7,2))

        Label(guessGeneration, text="Lagrangian:").grid(row=3, sticky='W', padx=(27,5))
        self.lagrangian = ('L1', 'L2')
        self.lagrangianCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.lagrangian, justify='center', width=10)
        self.lagrangianCB.set('L1')
        self.lagrangianCB.grid(row=3, column=1, sticky='W', padx=5, pady=3)

        Label(guessGeneration, text="Family:").grid(row=4, sticky='W', padx=(27,5))
        self.family = ('Northern', 'Southern')
        self.familyCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.family, justify='center', width=10)
        self.familyCB.set("Northern")
        self.familyCB.grid(row=4, column=1, sticky='W', padx=5, pady=3)

        Label(guessGeneration, text="Input value:").grid(row=5, sticky='W', padx=(27,5))
        self.inputValues = ('x', 'z', 'Period')
        self.inputValueCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.inputValues, justify='center', width=10)
        self.inputValueCB.set('x')
        self.inputValueCB.grid(row=5, column=1, sticky='W', padx=5, pady=3)

        self.guessValue = Entry(guessGeneration, width=12)
        self.guessValue.insert(0, 0.0)
        self.guessValue.config(state=DISABLED)
        self.guessValue.grid(row=6, column=1, sticky='W', padx=5)
        self.guessButton = Button(guessGeneration, text="Get Initial Guess", state=DISABLED, command=lambda: OrbitInputPage.guessGenerator(self, dynamicalSystem))
        self.guessButton.grid(row=7, sticky='W', padx=(27,5), pady=(15,0))

        orbitOptions = LabelFrame(self, text=" ORBIT OPTIONS ")
        orbitOptions.grid(row=1, column=0, columnspan=2, rowspan=4, sticky='WE', padx=5, pady=5, ipadx=5, ipady=5)
        Label(orbitOptions, text="Accuracy:").grid(row=0, sticky='W', padx=5)
        self.accuracyEntry = Entry(orbitOptions, width=7)
        self.accuracyEntry.grid(row=0, column=1, sticky='E', padx=5, pady=2)
        self.accuracy = 1.0e-8
        self.accuracyEntry.insert(END, self.accuracy)

        Label(orbitOptions, text="Max. Iterations:").grid(row=1, column=0, sticky='W', padx=5)
        self.maxIterEntry = Entry(orbitOptions, width=7)
        self.maxIterEntry.grid(row=1, column=1, sticky='E', padx=5, pady=2)
        self.maxIter = 10
        self.maxIterEntry.insert(END, self.maxIter)

        Label(orbitOptions, text="NRHO Stability Criteria:").grid(row=2, column=0, sticky='W', padx=5)
        self.stabilityCritEntry = Entry(orbitOptions, width=7)
        self.stabilityCritEntry.grid(row=2, column=1, sticky='E', padx=5, pady=2)
        self.stabilityCrit = 2
        self.stabilityCritEntry.insert(END, self.stabilityCrit)

        # button for calculation of entry data
        Button(self, text="Calculate", width=15, command=lambda: OrbitInputPage.calculate(self, dynamicalSystem)).grid(row=1, column=2, columnspan=2, sticky='SE', padx=5)
        self.nextButton = Button(self, text="Next", width=15, state=DISABLED, command=lambda: OrbitInputPage.nextPage(self, parent, controller, dynamicalSystem))
        self.nextButton.grid(row=2, column=2, columnspan=2, sticky='E', padx=5)
        Button(self, text="Clear Prompt", width=15, command=lambda: [self.statusBar.delete(1.0,END), self.statusBar.insert(INSERT, time.strftime("%Y-%m-%dT%H.%M.%S:\n\n>>>"))]).grid(row=3, column=2, columnspan=2, sticky='NE', padx=5)
        Button(self, text="Home", width=15, command=lambda: [controller.show_frame("StartPage"), self.destroy()]).grid(row=4, column=2, columnspan=2, sticky='NE', padx=5)
        Label(self, font=("TkDefaultFont", 11), text="Status:",).grid(row=5, sticky='NW', padx=5, pady=(15,5))
        self.statusBar = Text(self, height=7, relief=RIDGE, bd=3, spacing1=3)
        self.statusBar.grid(row=5, column=1, columnspan=3, sticky='W', padx=5, pady=(15,5))
        self.statusBar.insert(INSERT, time.strftime("%Y-%m-%dT%H.%M.%S:\n\n>>>"))
    # ------------------------------------------------------------------------------------------------------------------
    # This method sets the input fields for an individually generated initial guest to Active or Passive.
    # ------------------------------------------------------------------------------------------------------------------
    def customInput(self, modus):
        if modus == "NORMAL":
            self.lagrangianCB.config(state='readonly')
            self.familyCB.config(state='readonly')
            self.inputValueCB.config(state='readonly')
            self.guessValue.config(state=NORMAL)
            self.guessButton.config(state=NORMAL)
        elif modus == "DISABLED":
            self.lagrangianCB.config(state=DISABLED)
            self.familyCB.config(state=DISABLED)
            self.inputValueCB.config(state=DISABLED)
            self.guessValue.config(state=DISABLED)
            self.guessButton.config(state=DISABLED)
    # ------------------------------------------------------------------------------------------------------------------
    # This method calculates the initial Guess.
    # ------------------------------------------------------------------------------------------------------------------
    def guessGenerator(self, dynamicalSystem, lyapunov=None, lagrangian=None):
        try:
            if lyapunov:
                guess = InitialGuess(dynamicalSystem, lagrangian, "Northern", "z", 0)
            else:
                guess = InitialGuess(dynamicalSystem, self.lagrangianCB.get(), self.familyCB.get(), self.inputValueCB.get(), float(self.guessValue.get()))
            self.statusBar.insert(INSERT, "   Succesfully calculated Initial Guess.\n>>>")
            self.statusBar.see(END)
        except:
            self.statusBar.insert(INSERT, "   Initial guess could not be calculated. Please check input parameters.\n>>>")
            self.statusBar.see(END)
            return
        else:
            self.guessState = guess.x0
            self.guessPeriod = guess.tau
            self.xEntry.delete(0, "end")
            self.yEntry.delete(0, "end")
            self.zEntry.delete(0, "end")
            self.dxEntry.delete(0, "end")
            self.dyEntry.delete(0, "end")
            self.dzEntry.delete(0, "end")
            self.xEntry.insert(0, guess.x0[0])
            self.yEntry.insert(0, guess.x0[1])
            self.zEntry.insert(0, guess.x0[2])
            self.dxEntry.insert(0, guess.x0[3])
            self.dyEntry.insert(0, guess.x0[4])
            self.dzEntry.insert(0, guess.x0[5])
    # ------------------------------------------------------------------------------------------------------------------
    # This method instantiates a periodic orbit with the given initial guess.
    # ------------------------------------------------------------------------------------------------------------------
    def calculate(self, dynamicalSystem):
        try:
            x0 = np.array([float(self.xEntry.get()),
                           float(self.yEntry.get()),
                           float(self.zEntry.get()),
                           float(self.dxEntry.get()),
                           float(self.dyEntry.get()),
                           float(self.dzEntry.get())])

            Orbit.setAccuracy(float(self.accuracyEntry.get()))
            Orbit.setStabilityCriteria(float(self.stabilityCritEntry.get()))
            NumericalMethods.setMaxIterations(float(self.maxIterEntry.get()))

            if np.array_equal(x0, self.guessState):
                self.orbit = Orbit(x0, self.fixedValueCB.get(), dynamicalSystem, statusBar=self.statusBar, tau=self.guessPeriod)
            else:
                self.orbit = Orbit(x0, self.fixedValueCB.get(), dynamicalSystem, statusBar=self.statusBar)
            self.nextButton.config(state=NORMAL)
        except ValueError:
            self.statusBar.insert(INSERT, "   Input parameters are not complete.\n>>>")
            self.statusBar.see(END)
        except OverflowError:
            self.statusBar.insert(INSERT, "          Orbit could not be calculated. Integrating Initial State did not lead\n"
                                          "          to periodic orbit.\n>>>")
            self.statusBar.see(END)
        except StopIteration:
            self.statusBar.insert(INSERT, "          Orbit could not be calculated. Differential Corrections Method did\n"
                                          "          not converge.\n>>>")
            self.statusBar.see(END)
        except np.linalg.linalg.LinAlgError:
            self.statusBar.insert(INSERT, "          Orbit could not be calculated. Transformation matrix could not be\n"
                                          "          calculated because of singular matrix. Please fix another value.\n>>>")
            self.statusBar.see(END)
        else:
            pass
    # ------------------------------------------------------------------------------------------------------------------
    # This method calls the detail page for the orbit found.
    # ------------------------------------------------------------------------------------------------------------------
    def nextPage(self, parent, controller, dynamicalSystem):
            controller.frames["OrbitPage"] = OrbitPage.OrbitPage(parent=parent, controller=controller, orbit=self.orbit, dynamicalSystem=dynamicalSystem)
            controller.frames["OrbitPage"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("OrbitPage")
            self.update()
            time.sleep(0.5)
            self.destroy()
