"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 24.04.2018

OrbitFamily Class and L1Family/L2Family Subclasses
"""

# Imports
from Orbit import Orbit
import math
from Utility import NumericalMethods
#from Utility import OrbitContinuation
from Utility import Plot
import time
import numpy as np


# orbitfamily class
class OrbitFamily:

    def __init__(self, x0, orbitDistance, orbitNumber, mu, familyData=None):
        self.x0 = x0
        self.mu = mu
        self.orbitDistance = orbitDistance
        self.orbitNumber = orbitNumber
        self.familyData = familyData

    def getHaloFamily(self, lagrangian):
        # prints status update
        print("STATUS: Generation of a family of %2d Halo Orbits around %s...\n" % (self.orbitNumber, lagrangian))
        outData = OrbitContinuation.natParaConti(self.x0, self.orbitDistance, lagrangian, self.orbitNumber, self.mu)
        self.familyData = outData
        print("DONE")

    def getNRHOFamily(self, lagrangian):
        pass

    def writeData(self):
        output = open("Output/" + time.strftime("%Y-%m-%dT%H:%M:%S") + ".txt", "w")
        output.write("CREATION_DATE = " + time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"
                     "ORIGINATOR = Astos Solutions GmbH\n\n")
        output.write("META_START\n"
                     "META_STOP\n\n")
        output.write("DATA_START\n")
        output.write("        JC           Period           x              z            dy/dt\n\n")
        for orbit in self.familyData:
            output.write('{0:15.10f}'.format(orbit[6]))
            output.write('{0:15.10f}'.format(orbit[7]))
            output.write('{0:15.10f}'.format(orbit[0]))
            output.write('{0:15.10f}'.format(orbit[2]))
            output.write('{0:15.10f}\n'.format(orbit[4]))
        output.write("\nDATA_STOP")
        output.close()

    def plot(self, lagrangian, haloFamily="both", background="off"):
        if self.familyData is None:
            OrbitFamily.getHaloFamily(self, lagrangian)
        Plot.plot(self.familyData, self.mu, haloFamily, background)


class L1Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu, familyData=None):
        super().__init__(x0, orbitDistance, mu, familyData)
        self.lagrangian = "L1"

    def getHaloFamily(self):
        OrbitFamily.getHaloFamily(self, self.lagrangian)

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)


class L2Family(OrbitFamily):

    def __init__(self, x0, orbitDistance, mu, familyData=None):
        super().__init__(x0, orbitDistance, mu, familyData)
        self.lagrangian = "L2"

    def getHaloFamily(self):
        OrbitFamily.getHaloFamily(self, self.lagrangian)

    def plot(self, haloFamily="both", background="off"):
        OrbitFamily.plot(self, self.lagrangian, haloFamily, background)




class OrbitContinuation:
    accuracy = 1.0e-6

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
    def natParaConti(x0, orbitDistance, lagrangian, orbitNumber, mu):

        output = np.zeros((orbitNumber, 8))
        # calculates first two orbits
        x_n = x0
        for i in range(2):
            print("        Orbit Number: %2d    (fixed z-value)\n" % (i+1))
            orbit = Orbit(x_n, "z", mu)
            orbit.getJacobi()
            if i == 0:
                lastX = orbit.x0
                output[0][0:6] = orbit.x0
                output[0][6] = orbit.period
                output[0][7] = orbit.jacobi
                x_n = lastX - np.array([0, 0, orbitDistance, 0, 0, 0])
            else:
                outX = orbit.x0
                output[1][0:6] = orbit.x0
                output[1][6] = orbit.period
                output[1][7] = orbit.jacobi
        # loops through number of additional orbits
        for i in range(orbitNumber-2):
            # checks whether x- or z-value has changed more since last iteration
            if abs(lastX[0] - outX[0]) > abs(lastX[2] - outX[2]):
                # x-value changed more than z-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed x-value)\n" % (i + 3))
                dz = abs(outX[2] - lastX[2])
                stepSize = np.sqrt(orbitDistance ** 2 - dz ** 2)
                if math.isnan(stepSize):
                    stepSize = orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if lagrangian == "L1":
                    x_n = outX + np.array([stepSize, 0, 0, 0, 0, 0])
                elif lagrangian == "L2":
                    x_n = outX - np.array([stepSize, 0, 0, 0, 0, 0])
                else:
                    print("Lagrangian type not supported.")
                    exit()
                # saves last initial state for comparison of next iteration
                lastX = outX
                # calculates initial state
                orbit = Orbit(x_n, "x", mu)
                orbit.getJacobi()
                outX = orbit.x0
                output[i+2][0:6] = orbit.x0
                output[i+2][6] = orbit.period
                output[i+2][7] = orbit.jacobi

            else:
                # z-value changed more than x-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed z-value)\n" % (i + 3))
                dx = abs(outX[0] - lastX[0])
                stepSize = np.sqrt(orbitDistance ** 2 - dx ** 2)
                if math.isnan(stepSize):
                    stepSize = orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if lagrangian == "L1":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])
                elif lagrangian == "L2":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])
                else:
                    print("Lagrangian type not supported.")
                    exit()
                # saves last initial state for comparison of next iteration
                lastX = outX
                # calculates initial state
                orbit = Orbit(x_n, "z", mu)
                orbit.getJacobi()
                outX = orbit.x0
                output[i+2][0:6] = orbit.x0
                output[i+2][6] = orbit.period
                output[i+2][7] = orbit.jacobi


        return output


    @classmethod
    def setAccuracy(cls, accuracy):
        cls.accuracy = accuracy
