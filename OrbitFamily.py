"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 28.05.2018

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

    def __init__(self, x0, orbitDistance, system, orbitNumber=None):
        self.x0 = x0
        self.orbitDistance = orbitDistance
        self.system = system
        self.orbitNumber = orbitNumber
        self.familyData = None
        if x0[0] < 1:
            self.lagrangian = "L1"
        elif x0[0] > 1:
            self.lagrangian = "L2"
        else:
            print("Lagrangian type could not be determined.")
            self.lagrangian = None

    def getHaloFamily(self):
        if self.orbitNumber is None:
            print("Orbit Number is not given.")
            exit()
        # prints status update
        print("STATUS: Generation of a family of %d Halo Orbits around %s...\n" % (self.orbitNumber, self.lagrangian))
        OrbitContinuation.natParaConti(self)
        print("DONE")

    def getNRHOFamily(self):
        orbit = Orbit(self.x0, "z", self.system)
        orbit.getClosestNRHO()
        self.x0 = orbit.x0
        OrbitContinuation.natParaConti(self, NRHOFamily=True)

    def writeData(self):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        if not os.path.exists(OrbitFamily.dict):
            os.makedirs(OrbitFamily.dict)
        output = open(OrbitFamily.dict + "data.txt", "w")
        output.write("CREATION_DATE            =      " + time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"
                     "ORIGINATOR               =      ASTOS SOLUTIONS GMBH\n\n")
        output.write("META_START\n"
                     "NAME FIRST PRIMARY       =      %s\n"
                     "MASS FIRST PRIMARY       =      %e\n"
                     "NAME SECOND PRIMARY      =      %s\n"
                     "MASS SECOND PRIMARY      =      %e\n"
                     "PRIMARY DISTANCE         =      %e\n"
                     "MASS RATIO               =      %11.10f\n"
                     "LAGRANGIAN               =      %s\n"
                     "ORBIT NUMBER             =      %d\n"
                     "ORBIT DISTANCE           =      %f\n"
                     "META_STOP\n\n" % (self.system.nameFP, self.system.massFP, self.system.nameSP, self.system.massSP,
                                        self.system.distance, self.system.mu, self.lagrangian, self.orbitNumber, self.orbitDistance))
        output.write("DATA_START\n")
        output.write("        JC           Period           x              z            dy/dt\n")
        for orbit in self.familyData:
            output.write('{0:15.10f}'.format(orbit[0]))
            output.write('{0:15.10f}'.format(orbit[1]))
            output.write('{0:15.10f}'.format(orbit[2]))
            output.write('{0:15.10f}'.format(orbit[4]))
            output.write('{0:15.10f}\n'.format(orbit[6]))
        output.write("DATA_STOP")
        output.close()

    def plot(self):
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        Plot.plot(self.familyData, self.system, OrbitFamily.dict, self.lagrangian)




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
    def natParaConti(family, NRHOFamily=False):
        # calculates first two orbits
        x_n = family.x0
        print("        Orbit Number: 1    (fixed z-value)")
        for i in range(2):
            orbit = Orbit(x_n, "z", family.system, comment=False)
            if Orbit.error is True:
                family.orbitNumber = len(output)
                family.familyData = output
                return
            if i == 0:
                outX = orbit.x0
                orbit.getJacobi()
                output = orbit.data
                x_n = outX - np.array([0, 0, 0.0001, 0, 0, 0])
                print("        Reference Orbit:")
            else:
                lastX = orbit.x0

        stopLoop = False
        i = 0

        while not stopLoop:
        # loops through number of additional orbits
            # checks whether x- or z-value has changed more since last iteration
            if abs(lastX[0] - outX[0]) > abs(lastX[2] - outX[2]):
                # x-value changed more than z-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed x-value)" % (i + 2))
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
                orbit = Orbit(x_n, "x", family.system, comment=False)
                if Orbit.error is True:
                    family.orbitNumber = len(output)
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                output = np.vstack([output, orbit.data])

            else:
                # z-value changed more than x-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed z-value)" % (i + 2))
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
                orbit = Orbit(x_n, "z", family.system, comment=False)
                if Orbit.error is True:
                    family.orbitNumber = len(output)
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                output = np.vstack([output, orbit.data])

            if NRHOFamily is True:
                orbit.getStability()
                print(orbit.stability)
                if not orbit.NRHO:
                    stopLoop = True
            else:
                if i == (family.orbitNumber-2):
                    stopLoop=True

            i += 1

        family.familyData = output
        return
