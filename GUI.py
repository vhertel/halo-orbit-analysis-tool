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

LARGE_FONT = ("Verdana", 12)


class Window():

    @staticmethod
    def centerWindow(window, width, height):
        # get screen width and height
        screen_width = window.winfo_screenwidth()
        screen_height = window.winfo_screenheight()
        # calculate position x and y coordinates
        x = (screen_width/2) - (width/2)
        y = (screen_height/2) - (height/2)
        window.geometry('%dx%d+%d+%d' % (width, height, x, y))

    @staticmethod
    def fullScreen(window):
        # get screen width and height
        screen_width = window.winfo_screenwidth()
        screen_height = window.winfo_screenheight()
        Window.centerWindow(window, screen_width, screen_height)


# main class
class HaloTool(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        Tk.wm_title(self, "HALO TOOL")
        Window.centerWindow(self, 700, 600)
        container = Frame(self)
        container.pack(side="left", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}
        self.frames["StartPage"] = StartPage(parent=container, controller=self)
        self.frames["StartPage"].grid(row=0, column=0, sticky="nsew")
        self.setMenu()
        self.show_frame("StartPage")

    # sets menu
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

    # command for menu items
    def donothing(self):
        print("DO NOTHING")

    # shows new frames
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


# startpage
class StartPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        # main frame
        frameStartPage = Frame(self)
        frameStartPage.pack(expand=True, fill=BOTH, padx=10, pady=10)
        # describtion label

        fig = Label(frameStartPage, text="BLA")
        fig.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        # button for creating new instance of orbit or orbit family
        instanceButton = Button(frameStartPage, text="Create new instance",
                                command=lambda: StartPage.createNewInstance(self, parent, controller))
        instanceButton.grid(row=1, column=0, sticky='W', padx=5, pady=2)
        # button for creating new instance by loading data
        loadButton = Button(frameStartPage, text="Load from data",
                            command=lambda: StartPage.loadFromData(self, parent, controller))
        loadButton.grid(row=1, column=1, sticky='E', padx=5, pady=2)
        # button to quit program
        quitButton = Button(frameStartPage, text="Exit", command=self.client_exit)
        quitButton.grid(row=2, column=0, sticky='W', padx=5, pady=2)

    # loads next frame
    def createNewInstance(self, parent, controller):
        controller.frames["SystemPage"] = SystemPage(parent=parent, controller=controller)
        controller.frames["SystemPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("SystemPage")

    # loads next frame and chosen data
    def loadFromData(self, parent, controller):
        defaultPath = os.path.dirname(os.path.abspath(__file__)) + "/Output/"
        filePath = filedialog.askopenfilename(initialdir=defaultPath, title="Select file",
                                              filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
        file = open(filePath, "r")
        for lines in file:
            line = lines.split()
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
            if "DATA_START" in lines:
                for i in range(2):
                    line = file.readline()
                words = line.split()
        file.close()
        x0 = np.array([float(words[2]), 0, float(words[3]), 0, float(words[4]), 0])
        dynamicalSystem = System(nameFP, massFP, nameSP, massSP, distance)
        orbit = Orbit(x0, 1, dynamicalSystem)
        controller.frames["OrbitPage"] = OrbitPage(parent=parent, controller=controller, orbit=orbit,
                                                   dynamicalSystem=dynamicalSystem)
        controller.frames["OrbitPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitPage")

    # quits the program
    def client_exit(self):
        exit()


# page for configuration of the dynamical system
class SystemPage(Frame):
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
        earthMoon = Radiobutton(defaultSystems, text="Earth - Moon System", variable=var, value=0,
                                command=lambda: SystemPage.earthMoon(self))
        earthMoon.grid(row=0, sticky='W', padx=5, pady=2)
        sunEarth = Radiobutton(defaultSystems, text="Sun - Earth System", variable=var, value=1,
                               command=lambda: SystemPage.sunEarth(self))
        sunEarth.grid(row=1, sticky='W', padx=5, pady=2)
        sunMars = Radiobutton(defaultSystems, text="Sun - Mars System", variable=var, value=2,
                              command=lambda: SystemPage.sunMars(self))
        sunMars.grid(row=2, sticky='W', padx=5, pady=2)
        sunJupiter = Radiobutton(defaultSystems, text="Sun - Jupiter System", variable=var, value=3,
                                 command=lambda: SystemPage.sunJupiter(self))
        sunJupiter.grid(row=3, sticky='W', padx=5, pady=2)
        # control
        confirmButton = Button(self, text="Confirm", command=lambda: SystemPage.confirm(self, parent, controller))
        confirmButton.grid(row=7, column=1, sticky='E', padx=5, pady=2)
        homeButton = Button(self, text="Home", command=lambda: controller.show_frame("StartPage"))
        homeButton.grid(row=7, column=0, sticky='W', padx=5, pady=2)

    # fills entry with earth-moon data
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

    # fills entry with sun-earth data
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

    # fills entry with sun-mars data
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

    # fills entry with sun-jupiter data
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

    # confirms and stores data from entry
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
            controller.frames["NewInstance"] = NewInstance(parent=parent, controller=controller,
                                                           dynamicalSystem=dynamicalSystem)
            controller.frames["NewInstance"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("NewInstance")


# page for creating new instance
class NewInstance(Frame):
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        # single orbit
        button = Button(self, text="Single Orbit",
                        command=lambda: NewInstance.singleOrbit(self, parent, controller, dynamicalSystem))
        button.pack()
        # orbit family
        button2 = Button(self, text="Orbit Family",
                         command=lambda: NewInstance.orbitFamily(self, parent, controller, dynamicalSystem))
        button2.pack()
        # back to home
        homeButton = Button(self, text="Home", command=lambda: controller.show_frame("StartPage"))
        homeButton.pack()

    # loads next frame
    def singleOrbit(self, parent, controller, dynamicalSystem):
        controller.frames["OrbitInputPage"] = OrbitInputPage(parent=parent, controller=controller,
                                                             dynamicalSystem=dynamicalSystem)
        controller.frames["OrbitInputPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitInputPage")

    # loads next frame
    def orbitFamily(self, parent, controller, dynamicalSystem):
        controller.frames["FamilyInputPage"] = FamilyInputPage(parent=parent, controller=controller,
                                                     dynamicalSystem=dynamicalSystem)
        controller.frames["FamilyInputPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("FamilyInputPage")


# page for orbit details
class OrbitInputPage(Frame):

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
        Radiobutton(guessGeneration, text="Custom", variable=guessOption, value=3, command=lambda: OrbitInputPage.customInput(self, "NORMAL")).grid(row=2, columnspan=2, sticky='W', padx=5, pady=2)

        Label(guessGeneration, text="Lagrangian:").grid(row=3, sticky='W', padx=(27,5))
        self.lagrangian = ('L1', 'L2')
        self.lagrangianCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.lagrangian, justify='center', width=10)
        self.lagrangianCB.set('L1')
        self.lagrangianCB.grid(row=3, column=1, sticky='W', padx=5, pady=2)

        Label(guessGeneration, text="Family:").grid(row=4, sticky='W', padx=(27,5))
        self.family = ('Northern', 'Southern')
        self.familyCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.family, justify='center', width=10)
        self.familyCB.set("Northern")
        self.familyCB.grid(row=4, column=1, sticky='W', padx=5, pady=2)

        Label(guessGeneration, text="Input value:").grid(row=5, sticky='W', padx=(27,5))
        self.inputValues = ('x', 'z', 'Period')
        self.inputValueCB = ttk.Combobox(guessGeneration, state=DISABLED, values=self.inputValues, justify='center', width=10)
        self.inputValueCB.set('x')
        self.inputValueCB.grid(row=5, column=1, sticky='W', padx=5, pady=2)

        self.guessValue = Entry(guessGeneration, width=12)
        self.guessValue.insert(0, 0.0)
        self.guessValue.config(state=DISABLED)
        self.guessValue.grid(row=6, column=1, sticky='W', padx=5)
        self.guessButton = Button(guessGeneration, text="Get Initial Guess", state=DISABLED, command=lambda: OrbitInputPage.guessGenerator(self, dynamicalSystem))
        self.guessButton.grid(row=7, sticky='W', padx=(27,5))

        orbitOptions = LabelFrame(self, text=" ORBIT OPTIONS ")
        orbitOptions.grid(row=1, column=0, columnspan=2, rowspan=4, sticky='WE', padx=5, pady=5, ipadx=5, ipady=5)
        Label(orbitOptions, text="Accuracy:").grid(row=0, sticky='W', padx=5)
        self.accuracyEntry = Entry(orbitOptions, width=7)
        self.accuracyEntry.grid(row=0, column=1, sticky='W', padx=5, pady=2)
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
        self.stabilityCrit = 1
        self.stabilityCritEntry.insert(END, self.stabilityCrit)

        # button for calculation of entry data
        Button(self, text="Calculate", width=15, command=lambda: OrbitInputPage.calculate(self, parent, controller, dynamicalSystem)).grid(row=1, column=2, columnspan=2, sticky='SE', padx=5)
        self.nextButton = Button(self, text="Next", width=15, state=DISABLED, command=lambda: OrbitInputPage.nextPage(self, parent, controller, dynamicalSystem))
        self.nextButton.grid(row=2, column=2, columnspan=2, sticky='E', padx=5)
        Button(self, text="Clear Prompt", width=15, command=lambda: [self.statusBar.delete(1.0,END), self.statusBar.insert(INSERT, time.strftime("%Y-%m-%dT%H.%M.%S:\n\n>>>"))]).grid(row=3, column=2, columnspan=2, sticky='NE', padx=5)
        Button(self, text="Home", width=15, command=lambda: controller.show_frame("StartPage")).grid(row=4, column=2, columnspan=2, sticky='NE', padx=5)
        Label(self, font=("TkDefaultFont", 11), text="Status:").grid(row=5, sticky='NW', padx=5, pady=(15,5))
        self.statusBar = Text(self, height=7, relief=RIDGE, bd=3, spacing1=3)
        self.statusBar.grid(row=5, column=1, columnspan=3, sticky='W', padx=5, pady=(15,5))
        self.statusBar.insert(INSERT, time.strftime("%Y-%m-%dT%H.%M.%S:\n\n>>>"))

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

    # creates object with input state
    def calculate(self, parent, controller, dynamicalSystem):
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
                self.orbit = Orbit(x0, self.fixedValueCB.get(), dynamicalSystem, self.statusBar, tau=self.guessPeriod)
            else:
                self.orbit = Orbit(x0, self.fixedValueCB.get(), dynamicalSystem, self.statusBar)
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

    def nextPage(self, parent, controller, dynamicalSystem):
            controller.frames["OrbitPage"] = OrbitPage(parent=parent, controller=controller, orbit=self.orbit, dynamicalSystem=dynamicalSystem)
            controller.frames["OrbitPage"].grid(row=0, column=0, sticky="nsew")
            controller.show_frame("OrbitPage")


# page for orbit familiy
class FamilyInputPage(Frame):
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

        # button for calculation of entry data
        Button(self, text="Calculate", command=lambda: FamilyInputPage.calculate(self, parent, controller, dynamicalSystem)).grid(row=7, column=1, sticky='E', padx=5, pady=2)
        # back to home
        Button(self, text="Home", command=lambda: controller.show_frame("StartPage")).grid(row=7, column=0, sticky='W', padx=5, pady=2)

    # creates object with input state
    def calculate(self, parent, controller, dynamicalSystem):
        guess = InitialGuess(dynamicalSystem, self.lagrangian.get(), self.family.get(), "z", 0)
        orbitFamily = OrbitFamily(guess.x0, self.family.get(), self.lagrangian.get(), dynamicalSystem)
        controller.frames["OrbitFamilyPage"] = OrbitFamilyPage(parent=parent, controller=controller, orbitFamily=orbitFamily, initialState=guess.x0, dynamicalSystem=dynamicalSystem)
        controller.frames["OrbitFamilyPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitFamilyPage")


# detailed orbit page
class OrbitPage(Frame):
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
        Label(systemAttributes, text="Mass of " + dynamicalSystem.nameSP + ":").grid(row=3, column=0, padx=5,
                                                                                     sticky='W')
        Label(systemAttributes, text="%e kg" % (dynamicalSystem.massSP)).grid(row=3, column=1, padx=5, sticky='E')
        # distance between primaries
        Label(systemAttributes, text="Distance:").grid(row=4, column=0, padx=5, sticky='W')
        Label(systemAttributes, text="%.0f km" % (dynamicalSystem.distance/1.0e3)).grid(row=4, column=1, padx=5, sticky='E')
        # plot options
        plotOptions = LabelFrame(frameRight, text=" PLOT CONFIGURATION ")
        plotOptions.pack(fill=BOTH, expand=True, pady=(0, 7))
        Label(plotOptions, text="Invariant Manifolds:").grid(row=0, column=0, columnspan=4, padx=5, sticky='W')
        self.stable = IntVar(value=0)
        self.unstable = IntVar(value=0)
        Checkbutton(plotOptions, text="Stable Manifold", variable=self.stable, onvalue=1, offvalue=0).grid(row=0,
                                                        column=4, columnspan=2, padx=5, sticky='W')
        Checkbutton(plotOptions, text="Untable Manifold", variable=self.unstable, onvalue=1, offvalue=0).grid(row=1,
                                                        column=4, columnspan=2, padx=5, sticky='W')
        Label(plotOptions, text="Manifold Direction:").grid(row=2, column=0, columnspan=4, padx=5, sticky='W')
        self.manifoldDirection = IntVar(value=0)
        Radiobutton(plotOptions, text=dynamicalSystem.nameFP, variable=self.manifoldDirection, value=0).grid(row=2, column=4, padx=5, pady=2)
        Radiobutton(plotOptions, text=dynamicalSystem.nameSP, variable=self.manifoldDirection, value=1).grid(row=2, column=5, padx=5, pady=2)
        Label(plotOptions, text="Number of Manifolds:").grid(row=3, column=0, columnspan=3, padx=5, sticky='W')
        self.numberOfManifoldEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.numberOfManifoldEntry.grid(row=3, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultNumberOfManifolds = 30
        self.numberOfManifoldEntry.insert(END, self.defaultNumberOfManifolds)
        Label(plotOptions, text="Duration Factor:").grid(row=4, column=0, columnspan=3, padx=5, sticky='W')
        self.durationFactorEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.durationFactorEntry.grid(row=4, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultDurationFactor = 1.2
        self.durationFactorEntry.insert(END, self.defaultDurationFactor)
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
        Button(statusFrame, text="Home", command=lambda: [Window.centerWindow(controller, 700, 600), controller.show_frame("StartPage")]).pack(side=RIGHT, padx=5)
        # plots
        plotFrame = LabelFrame(frameLeft, text=" ORBIT PLOTS ")
        plotFrame.pack(fill=BOTH, expand=True)
        self.plot(plotFrame, orbit, dynamicalSystem)
        Window.fullScreen(controller)


    # opens input frame for configuration of what should be saved
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

    # saves data
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

    # plots data
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
        # stores input number of manifolds
        numberOfManifolds = int(self.numberOfManifoldEntry.get())
        # stores factor of integration
        durationFactor = float(self.durationFactorEntry.get())
        # calculates position of L1 and L2
        lagrangianPosition = Utility.lagrangianPosition(dynamicalSystem.mu)
        # checks whether input of manifold direction, number or duration factor has been changed
        if self.defaultNumberOfManifolds == int(self.numberOfManifoldEntry.get()) and self.defaultDurationFactor == float(self.durationFactorEntry.get()):
            inputChange = False
        else:
            inputChange = True
        # updates manifold number and duration factor
        self.defaultNumberOfManifolds = numberOfManifolds
        self.defaultDurationFactor = durationFactor
        # checks if manifolds need to be recalculated
        if ((self.stable.get() or self.unstable.get()) and orbit.stableManifolds is None) or (
                (self.stable.get() or self.unstable.get()) and inputChange):
            self.status.config(text="Calculating invariant Manifolds ...")
            self.status.update()
            # calculates initial states of (un)stable manifolds
            orbit.invariantManifolds(numberOfPoints=numberOfManifolds, direction=self.manifoldDirection.get())
            # declares arrays
            self.xStableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            self.yStableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            self.zStableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            self.xUnstableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            self.yUnstableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            self.zUnstableManifold = np.zeros((numberOfManifolds, int(durationFactor*50)))
            # integrates all states of (un)stable manifolds and stores data
            t = np.linspace(0, durationFactor * orbit.period, num=int(durationFactor*50))
            for i in range(numberOfManifolds):
                stableManifold = odeint(Utility.backwards, orbit.stableManifolds[i], t, args=(dynamicalSystem.mu,),
                                        rtol=2.5e-13, atol=1e-22)
                self.xStableManifold[i, :] = stableManifold[:, 0] * dynamicalSystem.distance
                self.yStableManifold[i, :] = stableManifold[:, 1] * dynamicalSystem.distance
                self.zStableManifold[i, :] = stableManifold[:, 2] * dynamicalSystem.distance
                unstableManifold = odeint(Utility.sysEquations, orbit.unstableManifolds[i], t,
                                          args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
                self.xUnstableManifold[i, :] = unstableManifold[:, 0] * dynamicalSystem.distance
                self.yUnstableManifold[i, :] = unstableManifold[:, 1] * dynamicalSystem.distance
                self.zUnstableManifold[i, :] = unstableManifold[:, 2] * dynamicalSystem.distance
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
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xStableManifold[i, :], self.zStableManifold[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xUnstableManifold[i, :], self.zUnstableManifold[i, :], color='red', linewidth=0.5)
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
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.yStableManifold[i, :], self.zStableManifold[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.yUnstableManifold[i, :], self.zUnstableManifold[i, :], color='red', linewidth=0.5)
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
        plt.tick_params(axis='both', which='both', labelsize=7, direction='out')
        if self.stable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xStableManifold[i, :], self.yStableManifold[i, :], color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                plt.plot(self.xUnstableManifold[i, :], self.yUnstableManifold[i, :], color='red', linewidth=0.5)
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
        if self.stable.get():
            for i in range(numberOfManifolds):
                ax.plot(self.xStableManifold[i, :], self.yStableManifold[i, :], self.zStableManifold[i, :],
                        color='green', linewidth=0.5)
        if self.unstable.get():
            for i in range(numberOfManifolds):
                ax.plot(self.xUnstableManifold[i, :], self.yUnstableManifold[i, :], self.zUnstableManifold[i, :],
                        color='red', linewidth=0.5)
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




# detailed orbit page
class OrbitFamilyPage(Frame):
    def __init__(self, parent, controller, orbitFamily, initialState, dynamicalSystem):
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
        Label(orbitAttributes, text="Initial State:", width=16, anchor=W).grid(row=0, column=0, padx=5, sticky='W', ipadx=5)
        self.state0Label = Label(orbitAttributes)
        self.state0Label.grid(row=0, column=1, padx=5, sticky='E')
        self.state1Label = Label(orbitAttributes)
        self.state1Label.grid(row=1, column=1, padx=5, sticky='E')
        self.state2Label = Label(orbitAttributes)
        self.state2Label.grid(row=2, column=1, padx=5, sticky='E')
        self.state3Label = Label(orbitAttributes)
        self.state3Label.grid(row=3, column=1, padx=5, sticky='E')
        self.state4Label = Label(orbitAttributes)
        self.state4Label.grid(row=4, column=1, padx=5, sticky='E')
        self.state5Label = Label(orbitAttributes)
        self.state5Label.grid(row=5, column=1, padx=5, pady=(0,15), sticky='E')
        # period
        Label(orbitAttributes, text="Period:").grid(row=7, column=0, padx=5, sticky='W')
        self.periodLabel = Label(orbitAttributes)
        self.periodLabel.grid(row=7, column=1, padx=5, sticky='E')
        # jacobi constant
        Label(orbitAttributes, text="Jacobi Constant:").grid(row=8, column=0, padx=5, sticky='W')
        self.jacobiLabel = Label(orbitAttributes)
        self.jacobiLabel.grid(row=8, column=1, padx=5, sticky='E')
        # stability index
        Label(orbitAttributes, text="Stability Index:").grid(row=9, column=0, padx=5, sticky='W')
        self.stabilityLabel = Label(orbitAttributes)
        self.stabilityLabel.grid(row=9, column=1, padx=5, sticky='E')
        # number of orbits
        Label(orbitAttributes, text="Number of Orbits:").grid(row=10, column=0, padx=5, sticky='W')
        self.NRHOLabel = Label(orbitAttributes)
        self.NRHOLabel.grid(row=10, column=1, padx=5, sticky='E')
        # lagrangian
        Label(orbitAttributes, text="Lagrangian:").grid(row=11, column=0, padx=5, sticky='W')
        self.lagrangianLabel = Label(orbitAttributes)
        self.lagrangianLabel.grid(row=11, column=1, padx=5, sticky='E')

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
        Label(plotOptions, text="Invariant Manifolds:").grid(row=0, column=0, columnspan=3, padx=5, sticky='W')
        self.stable = IntVar(value=0)
        self.unstable = IntVar(value=0)
        Checkbutton(plotOptions, text="Stable Manifold", variable=self.stable, onvalue=1, offvalue=0).grid(row=0, column=3, columnspan=3, padx=5, sticky='W')
        Checkbutton(plotOptions, text="Untable Manifold", variable=self.unstable, onvalue=1, offvalue=0).grid(row=1, column=3, columnspan=3, padx=5, sticky='W')
        Label(plotOptions, text="Number of Manifolds:").grid(row=2, column=0, columnspan=3, padx=5, sticky='W')
        self.numberOfManifoldEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.numberOfManifoldEntry.grid(row=2, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultNumberOfManifolds = 30
        self.numberOfManifoldEntry.insert(END, self.defaultNumberOfManifolds)
        Label(plotOptions, text="Duration Factor:").grid(row=3, column=0, columnspan=3, padx=5, sticky='W')
        self.durationFactorEntry = Entry(plotOptions, width=5, justify=CENTER)
        self.durationFactorEntry.grid(row=3, column=3, columnspan=3, padx=5, sticky='E')
        self.defaultDurationFactor = 1.2
        self.durationFactorEntry.insert(END, self.defaultDurationFactor)
        Label(plotOptions, text="Objects:").grid(row=4, column=0, columnspan=2, padx=5, sticky='W')
        self.nameFP = IntVar(value=0)
        self.nameSP = IntVar(value=0)
        self.L1 = IntVar(value=0)
        self.L2 = IntVar(value=0)
        Checkbutton(plotOptions, text=dynamicalSystem.nameFP, variable=self.nameFP, onvalue=1, offvalue=0).grid(row=4, column=2, padx=5)
        Checkbutton(plotOptions, text=dynamicalSystem.nameSP, variable=self.nameSP, onvalue=1, offvalue=0).grid(row=4, column=3, padx=5)
        Checkbutton(plotOptions, text="L1", variable=self.L1, onvalue=1, offvalue=0).grid(row=4, column=4, padx=5)
        Checkbutton(plotOptions, text="L2", variable=self.L2, onvalue=1, offvalue=0).grid(row=4, column=5, padx=5)
        # button for updating plots
        Button(plotOptions, text="Update Plots", command=lambda: OrbitPage.plot(self, plotFrame, orbit, dynamicalSystem)).grid(row=5, column=0, columnspan=3, padx=5, sticky='W')
        # button for saving data
        Button(plotOptions, text="Save to file", command=lambda: OrbitPage.saveButton(self, orbit, dynamicalSystem)).grid(row=5, column=3, columnspan=3, padx=5, sticky='E')
        # commands
        actions = LabelFrame(frameRight, text=" COMMANDS ")
        actions.pack(fill=BOTH, expand=True, pady=(0, 10))
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
        Button(statusFrame, text="Home", command=lambda: [Window.centerWindow(controller, 700, 600), controller.show_frame("StartPage")]).pack(side=RIGHT, padx=5)
        # plots
        plotFrame = LabelFrame(frameLeft, text=" ORBIT PLOTS ")
        plotFrame.pack(fill=BOTH, expand=True)
        Window.fullScreen(controller)

        orbit = Orbit(initialState, "x", dynamicalSystem)
        self.updateAttributes(orbit, dynamicalSystem)
        #orbitFamily = OrbitFamily(initialState, dynamicalSystem)
        self.plot(plotFrame, orbitFamily, dynamicalSystem)


    # opens input frame for configuration of what should be saved
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

    # saves data
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

    def updateAttributes(self, orbit, dynamicalSystem):
        self.state0Label.config(text="%.2f km" % (orbit.x0[0] * dynamicalSystem.distance / 1.0e3))
        self.state1Label.config(text="%.2f km" % (orbit.x0[1] * dynamicalSystem.distance / 1.0e3))
        self.state2Label.config(text="%.2f km" % (orbit.x0[2] * dynamicalSystem.distance / 1.0e3))
        self.state3Label.config(text="%.2f km/s" % (orbit.x0[3] * dynamicalSystem.distance / 1.0e3))
        self.state4Label.config(text="%.2f km/s" % (orbit.x0[4] * dynamicalSystem.distance / 1.0e3))
        self.state5Label.config(text="%.2f km/s" % (orbit.x0[5] * dynamicalSystem.distance / 1.0e3))
        period = (orbit.period * np.sqrt(dynamicalSystem.distance ** 3 /
                 (dynamicalSystem.G * (dynamicalSystem.massFP + dynamicalSystem.massSP)))) / (60 * 60 * 24)
        self.periodLabel.config(text="%.2f Days" % period)
        self.jacobiLabel.config(text=orbit.jacobi.round(2))
        self.lagrangianLabel.config(text=orbit.lagrangian)
        if orbit.stability is None:
            orbit.getStability()
        if orbit.NRHO:
            text = "Yes"
        else:
            text = "No"
        self.stabilityLabel.config(text=orbit.stability.round(2))
        self.NRHOLabel.config(text=text)

    def getClosestNRHO(self, plotFrame, orbit, dynamicalSystem):
        self.status.config(text="Searching for closest NRHO ...")
        self.status.update()
        var = IntVar()
        if abs(orbit.x0[2]) < 1.0e-4:
            self.NRHOLevel = Toplevel()
            Radiobutton(self.NRHOLevel, text="Up", variable=var, value=0).grid(row=0, column=0)
            Radiobutton(self.NRHOLevel, text="Down", variable=var, value=1).grid(row=1, column=0)
            test = Button(self.NRHOLevel, text="Confirm", command=lambda:[self.NRHOLevel.destroy(), orbit.getClosestNRHO(var.get()), self.plot(plotFrame, orbit, dynamicalSystem), self.updateAttributes(orbit, dynamicalSystem)])
            test.grid(row=3, column=0)
        else:
            orbit.getClosestNRHO()
            self.plot(plotFrame, orbit, dynamicalSystem)
            self.updateAttributes(orbit, dynamicalSystem)

    # plots data
    def plot(self, plotFrame, orbitFamily, dynamicalSystem):
        self.status.config(text="Plotting Data ...")
        self.status.update()
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


        lagrangianPosition = Utility.lagrangianPosition(dynamicalSystem.mu)
        jacobiMax = max(orbitFamily.familyData[:, 0])
        jacobiMin = min(orbitFamily.familyData[:, 0])
        normalize = mpl.colors.Normalize(vmin=jacobiMin, vmax=jacobiMax)
        colormap = plt.cm.hsv_r
        scalarmappaple = plt.cm.ScalarMappable(norm=normalize, cmap=colormap)
        scalarmappaple.set_array(5)




        # stores input number of manifolds
        numberOfManifolds = int(self.numberOfManifoldEntry.get())
        # stores factor of integration
        durationFactor = float(self.durationFactorEntry.get())
        # calculates position of L1 and L2
        lagrangianPosition = Utility.lagrangianPosition(dynamicalSystem.mu)
        # checks whether input of manifold number or duration factor has been changed
        if self.defaultNumberOfManifolds == int(
                self.numberOfManifoldEntry.get()) and self.defaultDurationFactor == float(
            self.durationFactorEntry.get()):
            inputChange = False
        else:
            inputChange = True


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

        for i in range(len(orbitFamily.familyData)):
            norm = (orbitFamily.familyData[i, 0] - jacobiMin) / (jacobiMax - jacobiMin)
            color = plt.cm.hsv_r(norm)
            t = np.linspace(0, orbitFamily.familyData[i, 1], num=300)
            halo = odeint(Utility.sysEquations, orbitFamily.familyData[i, 2:8], t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
            x = halo[:, 0] * dynamicalSystem.distance
            y = halo[:, 1] * dynamicalSystem.distance
            z = halo[:, 2] * dynamicalSystem.distance
            plt.figure(1)
            plt.plot(x, y, color=color, linewidth=0.25)
            plt.figure(2)
            plt.plot(x, z, color=color, linewidth=0.25)
            plt.figure(3)
            plt.plot(y, z, color=color, linewidth=0.25)
            ax.plot(x, y, z, color=color, linewidth=0.25)
        Plot.setAxesEqual(ax)

        self.status.config(text="")
        self.status.update()




class test(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        self.init_window(controller)

    def init_window(self, controller):
        stepOne = LabelFrame(self, text=" 1. Enter File Details: ")
        stepOne.grid(row=0, columnspan=7, sticky='W', \
                     padx=5, pady=5, ipadx=5, ipady=5)
        helpLf = LabelFrame(self, text=" Quick Help ")
        helpLf.grid(row=0, column=9, columnspan=2, rowspan=8, \
                    sticky='NS', padx=5, pady=5)
        helpLbl = Label(helpLf, text="Help will come - ask for it.")
        helpLbl.grid(row=0)
        stepTwo = LabelFrame(self, text=" 2. Enter Table Details: ")
        stepTwo.grid(row=2, columnspan=7, sticky='W', \
                     padx=5, pady=5, ipadx=5, ipady=5)
        stepThree = LabelFrame(self, text=" 3. Configure: ")
        stepThree.grid(row=3, columnspan=7, sticky='W', \
                       padx=5, pady=5, ipadx=5, ipady=5)
        inFileLbl = Label(stepOne, text="Select the File:")
        inFileLbl.grid(row=0, column=0, sticky='E', padx=5, pady=2)
        inFileTxt = Entry(stepOne)
        inFileTxt.grid(row=0, column=1, columnspan=7, sticky="WE", pady=3)
        inFileBtn = Button(stepOne, text="Browse ...")
        inFileBtn.grid(row=0, column=8, sticky='W', padx=5, pady=2)
        outFileLbl = Label(stepOne, text="Save File to:")
        outFileLbl.grid(row=1, column=0, sticky='E', padx=5, pady=2)
        outFileTxt = Entry(stepOne)
        outFileTxt.grid(row=1, column=1, columnspan=7, sticky="WE", pady=2)
        outFileBtn = Button(stepOne, text="Browse ...")
        outFileBtn.grid(row=1, column=8, sticky='W', padx=5, pady=2)
        inEncLbl = Label(stepOne, text="Input File Encoding:")
        inEncLbl.grid(row=2, column=0, sticky='E', padx=5, pady=2)
        inEncTxt = Entry(stepOne)
        inEncTxt.grid(row=2, column=1, sticky='E', pady=2)
        outEncLbl = Label(stepOne, text="Output File Encoding:")
        outEncLbl.grid(row=2, column=5, padx=5, pady=2)
        outEncTxt = Entry(stepOne)
        outEncTxt.grid(row=2, column=7, pady=2)
        outTblLbl = Label(stepTwo, \
                          text="Enter the name of the table to be used in the statements:")
        outTblLbl.grid(row=3, column=0, sticky='W', padx=5, pady=2)
        outTblTxt = Entry(stepTwo)
        outTblTxt.grid(row=3, column=1, columnspan=3, pady=2, sticky='WE')
        fldLbl = Label(stepTwo, \
                       text="Enter the field (column) names of the table:")
        fldLbl.grid(row=4, column=0, padx=5, pady=2, sticky='W')
        getFldChk = Checkbutton(stepTwo, \
                                text="Get fields automatically from input file", \
                                onvalue=1, offvalue=0)
        getFldChk.grid(row=4, column=1, columnspan=3, pady=2, sticky='WE')
        fldRowTxt = Entry(stepTwo)
        fldRowTxt.grid(row=5, columnspan=5, padx=5, pady=2, sticky='WE')
        transChk = Checkbutton(stepThree, \
                               text="Enable Transaction", onvalue=1, offvalue=0)
        transChk.grid(row=6, sticky='W', padx=5, pady=2)
        transRwLbl = Label(stepThree, \
                           text=" => Specify number of rows per transaction:")
        transRwLbl.grid(row=6, column=2, columnspan=2, \
                        sticky='W', padx=5, pady=2)
        transRwTxt = Entry(stepThree)
        transRwTxt.grid(row=6, column=4, sticky='WE')
