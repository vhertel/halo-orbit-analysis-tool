"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 24.04.2018

OrbitFamily Class and L1Family/L2Family Subclasses
"""

# Imports
from Orbit import Orbit
from Utility import Plot, NumericalMethods
import numpy as np
import math
import time
import os


# orbitfamily class
class OrbitFamily:

    dict = "Output/" + time.strftime("%Y-%m-%dT%H.%M.%S") + "/"

    def __init__(self, x0, orbitDistance, orbitNumber, mu, familyData=None):
        self.x0 = x0
        self.mu = mu
        self.orbitDistance = orbitDistance
        self.orbitNumber = orbitNumber
        self.familyData = familyData
        if x0[0] < 1:
            self.lagrangian = "L1"
        elif x0[0] > 1:
            self.lagrangian = "L2"
        else:
            print("Lagrangian type could not be determined.")
            self.lagrangian = None

    def getHaloFamily(self):
        # prints status update
        print("STATUS: Generation of a family of %d Halo Orbits around %s...\n" % (self.orbitNumber, self.lagrangian))
        OrbitContinuation.natParaConti(self)
        print("DONE")

    def getNRHOFamily(self):
        orbit = Orbit(self.x0, "z", self.mu)
        orbit.getClosestNRHO()
        self.x0 = orbit.x0
        OrbitContinuation.natParaConti(self)

    def writeData(self):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        if not os.path.exists(OrbitFamily.dict):
            os.makedirs(OrbitFamily.dict)
        output = open(OrbitFamily.dict + "data.txt", "w")
        output.write("CREATION_DATE        =      " + time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"
                     "ORIGINATOR           =      ASTOS SOLUTIONS GMBH\n\n")
        output.write("META_START\n"
                     "LAGRANGIAN           =      %s\n"
                     "ORBIT NUMBER         =      %d\n"
                     "ORBIT DISTANCE       =      %f\n"
                     "MASS RATIO           =      %11.10f\n"
                     "META_STOP\n\n" % (self.lagrangian, self.orbitNumber, self.orbitDistance, self.mu))
        output.write("DATA_START\n\n")
        output.write("        JC           Period           x              z            dy/dt\n\n")
        for orbit in self.familyData:
            output.write('{0:15.10f}'.format(orbit[0]))
            output.write('{0:15.10f}'.format(orbit[1]))
            output.write('{0:15.10f}'.format(orbit[2]))
            output.write('{0:15.10f}'.format(orbit[4]))
            output.write('{0:15.10f}\n'.format(orbit[6]))
        output.write("\nDATA_STOP")
        output.close()

    def plot(self, haloFamily="both", background="off", save=False):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        Plot.plot(self.familyData, self.mu, OrbitFamily.dict, haloFamily, background, save)

    def plotJacobi(self):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        Plot.plotJacobi(self.familyData, OrbitFamily.dict)

    def plotPeriod(self):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        Plot.plotPeriod(self.familyData, OrbitFamily.dict)




class OrbitContinuation:
    # --------------------------------------------------------------------------
    # NATURAL PARAMETER CONTINUATION
    #
    # DESCRIPTION:      Creates a family of halo orbits using the Natural Parameter Continuation,
    #                   a simple strategy based on a single converged solution to find and construct
    #                   related solutions. One parameter associated with the single converged solution
    #                   is incremented by a small, specific amount. This modified solution is now
    #                   employed as an initial guess for a new trajectory. Depending on the gradient of
    #                   the halo shape, either the x- or z-value is incremented by a specific stepsize,
    #                   that is calculated for each orbit based on the last two solutions for reaching
    #                   a uniformed plot with constant distances between the orbits.
    #
    #                   The distance between the orbits (orbitDistance) has been successfully tested for the range of 0.005 - 0.0075
    #
    # OUTPUT:           No output is generated. The orbits are plotted in the
    #                   figure that has been created before the function call
    #
    # SYNTAX:           natParaConti(x0, mu, epsilon, orbitNumber, orbitDistance, lagrangian, haloFamily, ax)
    #
    # x0            =   Initial State Guess
    # mu            =   Mass ratio of Primaries
    # epsilon       =   Error Tolerance of Constraints at T/2
    # orbitNumber   =   Number of Orbits to search for
    # orbitDistance =   Distance between Orbits
    # lagrangian    =   Lagrangian Point
    #                   Input possibilities: {"L1", "L2"}
    # haloFamily    =   Desired family of Halo Orbits:
    #                   Input possibilities: {"northern", "southern", "both"}
    # --------------------------------------------------------------------------
    @staticmethod
    def natParaConti(family):
        # calculates first two orbits
        x_n = family.x0
        print("        Orbit Number: 1    (fixed z-value)\n")
        for i in range(2):
            orbit = Orbit(x_n, "z", family.mu, comment=False)
            if Orbit.error == True:
                family.orbitNumber = len(output)
                family.familyData = output
                return
            if i == 0:
                outX = orbit.x0
                orbit.getJacobi()
                output = orbit.data
                x_n = outX - np.array([0, 0, 0.0001, 0, 0, 0])
            else:
                lastX = orbit.x0

        # loops through number of additional orbits
        for i in range(family.orbitNumber-1):

            # checks whether x- or z-value has changed more since last iteration
            if abs(lastX[0] - outX[0]) > abs(lastX[2] - outX[2]):
                # x-value changed more than z-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed x-value)\n" % (i + 2))
                dz = abs(outX[2] - lastX[2])
                stepSize = np.sqrt(family.orbitDistance**2 - dz**2)
                if math.isnan(stepSize):
                    stepSize = family.orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if family.lagrangian == "L1":
                    x_n = outX + np.array([stepSize, 0, 0, 0, 0, 0])   # +
                elif family.lagrangian == "L2":
                    x_n = outX - np.array([stepSize, 0, 0, 0, 0, 0])   # -
                else:
                    print("Lagrangian type not supported.")
                    exit()
                # calculates initial state
                orbit = Orbit(x_n, "x", family.mu, comment=False)
                if Orbit.error == True:
                    family.orbitNumber = len(output)
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                output = np.vstack([output, orbit.data])

            else:
                # z-value changed more than x-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed z-value)\n" % (i + 2))
                dx = abs(outX[0] - lastX[0])
                stepSize = np.sqrt(family.orbitDistance ** 2 - dx ** 2)
                if math.isnan(stepSize):
                    stepSize = family.orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if family.lagrangian == "L1":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])   # -
                elif family.lagrangian == "L2":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])   # -
                else:
                    print("Lagrangian type not supported.")
                    exit()
                # calculates initial state
                orbit = Orbit(x_n, "z", family.mu, comment=False)
                if Orbit.error == True:
                    family.orbitNumber = len(output)
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                output = np.vstack([output, orbit.data])

        family.familyData = output
        return
