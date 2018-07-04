"""
File    : OrbitFamily.py
Author  : Victor Hertel
Date    : 28.05.2018

Includes the OrbitFamily Class and OrbitContinuation Class
"""

import math
import os
import time
import numpy as np
from Orbit import Orbit
from Utility import Plot


class OrbitFamily:
    # path to the output folder
    dict = "Output/" + time.strftime("%Y-%m-%dT%H.%M.%S") + "/"

    # initializes by setting attributes and checking for the lagrangian
    def __init__(self, x0, system):
        # dynamical system
        self.system = system
        # initial guess of input orbit
        self.x0 = x0
        # distance between orbits
        self.orbitDistance = 0.007
        # sets data for writing and plotting
        self.familyData = None
        # checks for lagrangian
        if x0[0] < 1:
            self.lagrangian = "L1"
        elif x0[0] > 1:
            self.lagrangian = "L2"
        else:
            print("Lagrangian type could not be determined.")
            self.lagrangian = None
        # prints status update
        print("STATUS: Generation of a family of Halo Orbits around %s...\n" % (self.lagrangian))
        # natural parameter continuation
        OrbitContinuation.natParaConti(self)
        # prints status update
        print("DONE")

    # searches for closest NRHO and stores NRHO family
    def getNRHOFamily(self):
        orbit = Orbit(self.x0, "z", self.system)
        orbit.getClosestNRHO()
        self.x0 = orbit.x0
        OrbitContinuation.natParaConti(self, NRHOFamily=True)

    # writes family data including jacobi constant, period and initial states into file
    def writeData(self):
        # checks if the data has been calculated
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        # checks if the folder to store in already exists
        if not os.path.exists(OrbitFamily.dict):
            os.makedirs(OrbitFamily.dict)
        # creates file
        output = open(OrbitFamily.dict + "data.txt", "w")
        # writes header data
        output.write("CREATION_DATE            =      " + time.strftime("%Y-%m-%dT%H:%M:%S") + "\n"
                                                                                               "ORIGINATOR               =      ASTOS SOLUTIONS GMBH\n\n")
        # writes meta data
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
                                        self.system.distance, self.system.mu, self.lagrangian, len(self.familyData),
                                        self.orbitDistance))
        # writes orbit data
        output.write("DATA_START\n")
        output.write("        JC           Period           x              z            dy/dt\n")
        for orbit in self.familyData:
            output.write('{0:15.10f}'.format(orbit[0]))
            output.write('{0:15.10f}'.format(orbit[1]))
            output.write('{0:15.10f}'.format(orbit[2]))
            output.write('{0:15.10f}'.format(orbit[4]))
            output.write('{0:15.10f}\n'.format(orbit[6]))
        output.write("DATA_STOP")
        # closes file
        output.close()

    # plots data depending on user input
    def plot(self):
        # checks if the data has been calculated
        if self.familyData is None:
            print("        No data has been calculated yet\n"
                  "DONE")
            return
        Plot.plot(self.familyData, self.system, OrbitFamily.dict, self.lagrangian)


class OrbitContinuation:

    # natural parameter continuation method
    @staticmethod
    def natParaConti(family, NRHOFamily=False):
        # initial guess of input orbit
        x_n = family.x0
        # prints status update
        print("        Orbit Number: 1    (fixed z-value)")
        # calculates first two orbits
        for i in range(2):
            # creates orbit object
            orbit = Orbit(x_n, "z", family.system, comment=False)
            # checks if error occurred during orbit instantiation
            if orbit.error is True:
                # sets output data
                family.familyData = output
                return
            if i == 0:
                # stores initial state of first orbit
                outX = orbit.x0
                # calculates jacobi constant of first orbit
                orbit.getJacobi()
                # sets output data
                output = orbit.data
                # modifies initial state to initial guess of second orbit
                x_n = outX - np.array([0, 0, 0.0001, 0, 0, 0])
                # prints status update
                print("        Reference Orbit:")
            else:
                # stores initial state of second orbit for comparison
                lastX = orbit.x0

        # stop criteria
        stopLoop = False
        # counter
        i = 0

        # loops through orbits
        while not stopLoop:
            # checks whether x- or z-value has changed more since last iteration
            if abs(lastX[0] - outX[0]) > abs(lastX[2] - outX[2]):
                # x-value changed more than z-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed x-value)" % (i + 2))
                # calculates the change of z
                dz = abs(outX[2] - lastX[2])
                # calculates new stepsize for next orbit to reach constant distances in between
                stepSize = np.sqrt(family.orbitDistance ** 2 - dz ** 2)
                # reduces distance if mathematical error occurred
                if math.isnan(stepSize):
                    stepSize = family.orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if family.lagrangian == "L1":
                    x_n = outX + np.array([stepSize, 0, 0, 0, 0, 0])
                elif family.lagrangian == "L2":
                    x_n = outX - np.array([stepSize, 0, 0, 0, 0, 0])
                # calculates initial state of next orbit
                orbit = Orbit(x_n, "x", family.system, comment=False)
                # checks if error occurred during orbit instantiation
                if orbit.error is True:
                    stopLoop = True
                    # sets output data
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                # includes orbit data to output data
                output = np.vstack([output, orbit.data])

            else:
                # z-value changed more than x-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed z-value)" % (i + 2))
                # calculates the change of x
                dx = abs(outX[0] - lastX[0])
                # calculates new stepsize for next orbit to reach constant distances in between
                stepSize = np.sqrt(family.orbitDistance ** 2 - dx ** 2)
                # reduces distance if mathematical error occurred
                if math.isnan(stepSize):
                    stepSize = family.orbitDistance / 2
                # generates next initial guess depending on continuation direction
                if family.lagrangian == "L1":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])
                elif family.lagrangian == "L2":
                    x_n = outX - np.array([0, 0, stepSize, 0, 0, 0])
                # calculates initial state of next orbit
                orbit = Orbit(x_n, "z", family.system, comment=False)
                # checks if error occurred during orbit instantiation
                if orbit.error is True:
                    stopLoop = True
                    # sets output data
                    family.familyData = output
                    return
                # saves last initial state for comparison of next iteration
                lastX = outX
                outX = orbit.x0
                # includes orbit data to output data
                output = np.vstack([output, orbit.data])

            # when seachring for NRHO family stability index is checked for every orbit
            if NRHOFamily is True:
                # calculates stability index
                orbit.getStability()
                if not orbit.NRHO:
                    # sets stop criteria when stability index out of range
                    stopLoop = True
            # increments counter
            i += 1

        # sets output data
        family.familyData = output
        return
