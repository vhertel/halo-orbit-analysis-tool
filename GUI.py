

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from tkinter import *
from Utility import Utility, System, Plot
from Orbit import Orbit
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


LARGE_FONT= ("Verdana", 12)

class HaloTool(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)
        Tk.wm_title(self, "HALO TOOL")
        Tk.geometry(self, "700x300")
        container = Frame(self)
        container.pack(side="left", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}
        self.frames["StartPage"] = StartPage(parent=container, controller=self)
        self.frames["StartPage"].grid(row=0, column=0, sticky="nsew")
        self.setMenu()
        self.show_frame("StartPage")
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
    def donothing(self):
       print("DO NOTHING")
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

class StartPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        frameStartPage = LabelFrame(self, text="HALO DETERMINATION CONTINUATION")
        frameStartPage.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        fig = Label(frameStartPage, text="BLA")
        fig.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        button = Button(frameStartPage, text="Create new instance", command=lambda: StartPage.createNewInstance(self, parent, controller))
        button.grid(row=1, column=0, sticky='W', padx=5, pady=2)
        button2 = Button(frameStartPage, text="Load from data", command=lambda: StartPage.loadFromData(self, parent, controller))
        button2.grid(row=1, column=1, sticky='E', padx=5, pady=2)
        quitButton = Button(frameStartPage, text="Exit", command=self.client_exit)
        quitButton.grid(row=2, column=0, sticky='W', padx=5, pady=2)
    def createNewInstance(self, parent, controller):
        controller.frames["SystemPage"] = SystemPage(parent=parent, controller=controller)
        controller.frames["SystemPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("SystemPage")
    def loadFromData(self, parent, controller):
        controller.frames["LoadFromData"] = LoadFromData(parent=parent, controller=controller)
        controller.frames["LoadFromData"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("LoadFromData")

    def client_exit(self):
        exit()


class SystemPage(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        stepOne = LabelFrame(self, text=" DYNAMICAL SYSTEM ")
        stepOne.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        describtion = Label(stepOne, text="Describtion")
        describtion.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        self.nameFPLabel = Label(stepOne, text="Name of first Primary:")
        self.nameFPLabel.grid(row=1, column=0, sticky='W', padx=5, pady=2)
        self.nameFPEntry = Entry(stepOne)
        self.nameFPEntry.grid(row=1, column=1, sticky='W', padx=5, pady=2)
        self.massFPLabel = Label(stepOne, text="Mass of first Primary:")
        self.massFPLabel.grid(row=2, column=0, sticky='W', padx=5, pady=2)
        self.massFPEntry = Entry(stepOne)
        self.massFPEntry.grid(row=2, column=1, sticky='W', padx=5, pady=2)
        unit = Label(stepOne, text="kg")
        unit.grid(row=2, column=2, sticky='W')
        self.nameSPLabel = Label(stepOne, text="Name of second Primary:")
        self.nameSPLabel.grid(row=3, column=0, sticky='W', padx=5, pady=2)
        self.nameSPEntry = Entry(stepOne)
        self.nameSPEntry.grid(row=3, column=1, sticky='W', padx=5, pady=2)
        self.massSPLabel = Label(stepOne, text="Mass of second Primary:")
        self.massSPLabel.grid(row=4, column=0, sticky='W',padx=5, pady=2)
        self.massSPEntry = Entry(stepOne)
        self.massSPEntry.grid(row=4, column=1, sticky='W',padx=5, pady=2)
        unit = Label(stepOne, text="kg")
        unit.grid(row=4, column=2, sticky='W')
        self.distanceLabel = Label(stepOne, text="Distance of Primaries:")
        self.distanceLabel.grid(row=5, column=0, sticky='W',padx=5, pady=2)
        self.distanceEntry = Entry(stepOne)
        self.distanceEntry.grid(row=5, column=1, sticky='W',padx=5, pady=2)
        unit = Label(stepOne, text="m")
        unit.grid(row=5, column=2, sticky='W')
        defaultSystems = LabelFrame(self, text=" DEFAULT SYSTEMS ")
        defaultSystems.grid(row=0, column=9, columnspan=2, rowspan=7, sticky='NS', padx=5, pady=5)
        var = IntVar()
        R1 = Radiobutton(defaultSystems, text="Earth - Moon System", variable=var, value=0, command=lambda: SystemPage.earthMoon(self))
        R1.grid(row=0, sticky='W', padx=5, pady=2)
        R2 = Radiobutton(defaultSystems, text="Sun - Earth System", variable=var, value=1, command=lambda: SystemPage.sunEarth(self))
        R2.grid(row=1, sticky='W', padx=5, pady=2)
        R3 = Radiobutton(defaultSystems, text="Sun - Mars System", variable=var, value=2, command=lambda: SystemPage.sunMars(self))
        R3.grid(row=2, sticky='W', padx=5, pady=2)
        R4 = Radiobutton(defaultSystems, text="Sun - Jupiter System", variable=var, value=3, command=lambda: SystemPage.sunJupiter(self))
        R4.grid(row=3, sticky='W', padx=5, pady=2)
        confirmButton = Button(self, text="Confirm", command=lambda: SystemPage.confirm(self, parent, controller))
        confirmButton.grid(row=7, column=1, sticky='E', padx=5, pady=2)
        backButton = Button(self, text="Back", command=lambda: controller.show_frame("StartPage"))
        backButton.grid(row=7, column=0, sticky='W',padx=5, pady=2)
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
        self.distanceEntry.insert(0, 384402*1.0e3)
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
        self.distanceEntry.insert(0, 149597870.7*1.0e3)
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
        self.distanceEntry.insert(0, 227900000*1.0e3)
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
        self.distanceEntry.insert(0, 778547200*1.0e3)
    def confirm(self, parent, controller):
        nameFP = self.nameFPEntry.get()
        massFP = float(self.massFPEntry.get())
        nameSP = self.nameSPEntry.get()
        massSP = float(self.massSPEntry.get())
        distance = float(self.distanceEntry.get())
        dynamicalSystem = System(nameFP=nameFP, massFP=massFP, nameSP=nameSP, massSP=massSP, distance=distance)
        controller.frames["NewInstance"] = NewInstance(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem)
        controller.frames["NewInstance"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("NewInstance")

class LoadFromData(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        backButton = Button(self, text="Back", command=lambda: controller.show_frame("StartPage"))
        backButton.pack()

class NewInstance(Frame):
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        button = Button(self, text="Single Orbit", command=lambda: NewInstance.singleOrbit(self, parent, controller, dynamicalSystem))
        button.pack()
        button2 = Button(self, text="Orbit Family", command=lambda: NewInstance.orbitFamily(self, parent, controller, dynamicalSystem))
        button2.pack()
        backButton = Button(self, text="Back", command=lambda: controller.show_frame("SystemPage"))
        backButton.pack()
    def singleOrbit(self, parent, controller, dynamicalSystem):
        controller.frames["OrbitInputPage"] = OrbitInputPage(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem)
        controller.frames["OrbitInputPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitInputPage")
    def orbitFamily(self, parent, controller, dynamicalSystem):
        controller.frames["FamilyPage"] = FamilyPage(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem)
        controller.frames["FamilyPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("FamilyPage")


class OrbitInputPage(Frame):
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        stateFrame = LabelFrame(self, text=" INITIAL STATE ")
        stateFrame.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        describtion = Label(stateFrame, text="Input of Initial State")
        describtion.grid(row=0, column=0, sticky='W', padx=5, pady=2)
        self.xLabel = Label(stateFrame, text="X")
        self.xLabel.grid(row=1, column=0, sticky='W', padx=5, pady=2)
        self.xEntry = Entry(stateFrame)
        self.xEntry.grid(row=1, column=1, sticky='W', padx=5, pady=2)
        self.yLabel = Label(stateFrame, text="Y")
        self.yLabel.grid(row=2, column=0, sticky='W', padx=5, pady=2)
        self.yEntry = Entry(stateFrame)
        self.yEntry.grid(row=2, column=1, sticky='W', padx=5, pady=2)
        self.zLabel = Label(stateFrame, text="Z")
        self.zLabel.grid(row=3, column=0, sticky='W', padx=5, pady=2)
        self.zEntry = Entry(stateFrame)
        self.zEntry.grid(row=3, column=1, sticky='W', padx=5, pady=2)
        self.dxLabel = Label(stateFrame, text="XDOT")
        self.dxLabel.grid(row=4, column=0, sticky='W',padx=5, pady=2)
        self.dxEntry = Entry(stateFrame)
        self.dxEntry.grid(row=4, column=1, sticky='W',padx=5, pady=2)
        self.dyLabel = Label(stateFrame, text="YDOT")
        self.dyLabel.grid(row=5, column=0, sticky='W',padx=5, pady=2)
        self.dyEntry = Entry(stateFrame)
        self.dyEntry.grid(row=5, column=1, sticky='W',padx=5, pady=2)
        self.dzLabel = Label(stateFrame, text="ZDOT")
        self.dzLabel.grid(row=6, column=0, sticky='W',padx=5, pady=2)
        self.dzEntry = Entry(stateFrame)
        self.dzEntry.grid(row=6, column=1, sticky='W',padx=5, pady=2)

        guessGeneration = LabelFrame(self, text=" INITIAL GUESS GENERATION ")
        guessGeneration.grid(row=0, column=9, columnspan=2, rowspan=7, sticky='NS', padx=5, pady=5)
        var = IntVar()
        R1 = Radiobutton(guessGeneration, text="L1", variable=var, value=0, command=lambda: OrbitInputPage.l1(self))
        R1.grid(row=0, sticky='W', padx=5, pady=2)
        R2 = Radiobutton(guessGeneration, text="L1 Middle", variable=var, value=1, command=lambda: OrbitInputPage.l1Middle(self))
        R2.grid(row=1, sticky='W', padx=5, pady=2)
        R3 = Radiobutton(guessGeneration, text="L2", variable=var, value=2, command=lambda: OrbitInputPage.l2(self))
        R3.grid(row=2, sticky='W', padx=5, pady=2)
        R4 = Radiobutton(guessGeneration, text="L2 Middle", variable=var, value=3, command=lambda: OrbitInputPage.l2Middle(self))
        R4.grid(row=3, sticky='W', padx=5, pady=2)

        calculate = Button(self, text="Calculate", command=lambda: OrbitInputPage.calculate(self, parent, controller, dynamicalSystem))
        calculate.grid(row=7, column=1, sticky='E', padx=5, pady=2)
        backButton = Button(self, text="Back", command=lambda: controller.show_frame("NewInstance"))
        backButton.grid(row=7, column=0, sticky='W',padx=5, pady=2)



    def l1(self):
        self.xEntry.delete(0, "end")
        self.yEntry.delete(0, "end")
        self.zEntry.delete(0, "end")
        self.dxEntry.delete(0, "end")
        self.dyEntry.delete(0, "end")
        self.dzEntry.delete(0, "end")
        self.xEntry.insert(0, 0.8233901862)
        self.yEntry.insert(0, 0)
        self.zEntry.insert(0, -0.0029876370)
        self.dxEntry.insert(0, 0)
        self.dyEntry.insert(0, 0.1264751431)
        self.dzEntry.insert(0, 0)
    def l1Middle(self):
        self.xEntry.delete(0, "end")
        self.yEntry.delete(0, "end")
        self.zEntry.delete(0, "end")
        self.dxEntry.delete(0, "end")
        self.dyEntry.delete(0, "end")
        self.dzEntry.delete(0, "end")
        self.xEntry.insert(0, 0.8235990912)
        self.yEntry.insert(0, 0)
        self.zEntry.insert(0, -0.0399866715)
        self.dxEntry.insert(0, 0)
        self.dyEntry.insert(0, 0.1492106867)
        self.dzEntry.insert(0, 0)
    def l2(self):
        self.xEntry.delete(0, "end")
        self.yEntry.delete(0, "end")
        self.zEntry.delete(0, "end")
        self.dxEntry.delete(0, "end")
        self.dyEntry.delete(0, "end")
        self.dzEntry.delete(0, "end")
        self.xEntry.insert(0, 1.1808881373)
        self.yEntry.insert(0, 0)
        self.zEntry.insert(0, -0.0032736457)
        self.dxEntry.insert(0, 0)
        self.dyEntry.insert(0, -0.1559184478)
        self.dzEntry.insert(0, 0)
    def l2Middle(self):
        self.xEntry.delete(0, "end")
        self.yEntry.delete(0, "end")
        self.zEntry.delete(0, "end")
        self.dxEntry.delete(0, "end")
        self.dyEntry.delete(0, "end")
        self.dzEntry.delete(0, "end")
        self.xEntry.insert(0, 1.1542349115)
        self.yEntry.insert(0, 0)
        self.zEntry.insert(0, -0.1379744940)
        self.dxEntry.insert(0, 0)
        self.dyEntry.insert(0, -0.2147411949)
        self.dzEntry.insert(0, 0)

    def calculate(self, parent, controller, dynamicalSystem):
        x0 = np.array([float(self.xEntry.get()),
                       float(self.yEntry.get()),
                       float(self.zEntry.get()),
                       float(self.dxEntry.get()),
                       float(self.dyEntry.get()),
                       float(self.dzEntry.get())])
        orbit = Orbit(x0, "x", dynamicalSystem)
        controller.frames["OrbitPage"] = OrbitPage(parent=parent, controller=controller, dynamicalSystem=dynamicalSystem, orbit=orbit)
        controller.frames["OrbitPage"].grid(row=0, column=0, sticky="nsew")
        controller.show_frame("OrbitPage")




class FamilyPage(Frame):
    def __init__(self, parent, controller, dynamicalSystem):
        Frame.__init__(self, parent)
        backButton = Button(self, text="Back", command=lambda: controller.show_frame("NewInstance"))
        backButton.pack()

class OrbitPage(Frame):
    def __init__(self, parent, controller, dynamicalSystem, orbit):
        Frame.__init__(self, parent)
        Tk.geometry(controller, "1200x700")
        plotFrame = LabelFrame(self, text=" ORBIT PLOTS ",)
        plotFrame.grid(row=0, columnspan=7, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
        attributes = LabelFrame(self, text=" ORBIT ATTRIBUTES ")
        attributes.grid(row=0, column=9, columnspan=2, rowspan=8, sticky='NS', padx=5, pady=5)
        stateLabel = Label(attributes, text="Initial State:")
        stateLabel.grid(row=0, column=0, sticky='W')
        state1 = Label(attributes, text=orbit.x0[0])
        state1.grid(row=0, column=1, sticky='W')
        state2 = Label(attributes, text=orbit.x0[1])
        state2.grid(row=1, column=1, sticky='W')
        state3 = Label(attributes, text=orbit.x0[2])
        state3.grid(row=2, column=1, sticky='W')
        state4 = Label(attributes, text=orbit.x0[3])
        state4.grid(row=3, column=1, sticky='W')
        state5 = Label(attributes, text=orbit.x0[4])
        state5.grid(row=4, column=1, sticky='W')
        state6 = Label(attributes, text=orbit.x0[5])
        state6.grid(row=5, column=1, sticky='W')

        periodLabel = Label(attributes, text="Period:")
        periodLabel.grid(row=7, column=0, sticky='W')
        period = Label(attributes, text=orbit.period)
        period.grid(row=7, column=1, sticky='W')

        jacobiLabel = Label(attributes, text="Jacobi Constant:")
        jacobiLabel.grid(row=8, column=0, sticky='W')
        jacobi = Label(attributes, text=orbit.jacobi)
        jacobi.grid(row=8, column=1, sticky='W')

        stabilityLabel = Label(attributes, text="Stability Index:")
        stabilityLabel.grid(row=8, column=0, sticky='W')
        stability = Label(attributes, text=orbit.stability)
        stability.grid(row=8, column=1, sticky='W')


        backButton = Button(attributes, text="Back", command=lambda: controller.show_frame("OrbitInputPage"))
        backButton.grid(row=9, column=0, sticky='W')



        t = np.linspace(0, orbit.period, num=10000)
        orbitStates = odeint(Utility.sysEquations, orbit.x0, t, args=(dynamicalSystem.mu,), rtol=2.5e-13, atol=1e-22)
        x = orbitStates[:, 0] * dynamicalSystem.distance
        y = orbitStates[:, 1] * dynamicalSystem.distance
        z = orbitStates[:, 2] * dynamicalSystem.distance

        self.fig1 = plt.figure(figsize=(3,3))
        plt.axis("equal")
        self.canvas = FigureCanvasTkAgg(self.fig1, master=plotFrame)
        self.canvas.get_tk_widget().pack(side='top', fill='both')
        self.canvas._tkcanvas.grid(row=0, column=0, sticky="nsew")
        plt.title("$xz$-Projection")
        plt.xlabel("$x$ [m]")
        plt.ylabel("$z$ [m]")
        plt.plot(x, z, color='blue', linewidth=0.5)

        self.fig2 = plt.figure(figsize=(3,3))
        plt.axis("equal")
        self.canvas = FigureCanvasTkAgg(self.fig2, master=plotFrame)
        self.canvas.get_tk_widget().pack(side='top', fill='both')
        self.canvas._tkcanvas.grid(row=0, column=1, sticky="nsew")
        plt.title("$yz$-Projection")
        plt.xlabel("$y$ [m]")
        plt.ylabel("$z$ [m]")
        plt.plot(y, z, color='blue', linewidth=0.5)

        self.fig3 = plt.figure(figsize=(3,3))
        plt.axis("equal")
        self.canvas = FigureCanvasTkAgg(self.fig3, master=plotFrame)
        self.canvas.get_tk_widget().pack(side='top', fill='both')
        self.canvas._tkcanvas.grid(row=1, column=0, sticky="nsew")
        plt.title("$xy$-Projection")
        plt.xlabel("$x$ [m]")
        plt.ylabel("$y$ [m]")
        plt.plot(x, y, color='blue', linewidth=0.5)

        self.fig4 = plt.figure(figsize=(3,3))
        self.canvas = FigureCanvasTkAgg(self.fig4, master=plotFrame)
        self.canvas.get_tk_widget().pack(side='top', fill='both')
        self.canvas._tkcanvas.grid(row=1, column=1, sticky="nsew")
        ax = Axes3D(self.fig4)
        ax.plot(x, y, z, color='blue', linewidth=0.5)
        Plot.setAxesEqual(ax)









class test(Frame):
    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        print("INTI")
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
                               text="Get fields automatically from input file",\
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




