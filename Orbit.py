"""
File    : Orbit.py
Author  : Victor Hertel
Date    : 28.05.2018

Orbit Class and NRHO Subclass
"""

# Imports
import numpy as np
from Utility import NumericalMethods
from Utility import Plot
from Utility import Utility
import time




# orbit class
class Orbit:

    accuracy = 1.0e-8
    stabilityCriteria = 5
    error = False
    dict = "Output/" + time.strftime("%Y-%m-%dT%H.%M.%S") + "/"


    # initializes by adapting input state to initial state of periodic halo orbits
    def __init__(self, x0, fixedValue, system, comment="True"):
        self.system = system
        self.stability = None
        self.NRHO = None
        # prints status update
        if comment:
            print("STATUS: Adapting input state to periodic Halo Orbit...\n")
        # stores initial state, period and constraints of halo orbit in 1x3 vector outData
        try:
            outData = NumericalMethods.diffCorrections(x0, self.system.mu, Orbit.accuracy, fixedValue=fixedValue)
        except ValueError:
            Orbit.error = True
            return
        # sets attributes of halo orbit
        self.x0 = outData[0:6]
        self.period = outData[6]
        Orbit.getJacobi(self)
        self.data = np.array([self.jacobi, self.period, self.x0[0], self.x0[1], self.x0[2], self.x0[3], self.x0[4], self.x0[5]])
        if comment:
            print("DONE")


    def getClosestNRHO(self):
        if Orbit.error is True:
            return
        # prints status update
        print("STATUS: Adapting input state to initial state of nearest NRHO...\n")
        # calculates highest stability index when attribute not given
        if self.stability is None:
            Orbit.getStability(self)
        # checks whether orbit is aready NRHO
        if self.NRHO is True:
            print("        Halo Orbit is already NRHO.")
        else:
            # calculates T/2 by integrating until y changes sign
            try:
                tau_n = Utility.halfPeriod(self.x0, self.system.mu, 1.0e-11)
            except ValueError:
                return
            stepSize = 0.01
            while abs(self.stability - Orbit.stabilityCriteria) > 1.0e-2:
                # searches for next NRHO using pseudo-arclength continuation method
                while self.NRHO is False:
                    print("SCHLEIFE")
                    lastX = self.x0
                    lastPeriod = self.period
                    print(self.x0)
                    outData = NumericalMethods.diffCorrections(self.x0, self.system.mu, Orbit.accuracy, tau=tau_n)
                    x_n = outData[0]
                    print(x_n)
                    tau_n = outData[1]
                    phi = outData[2]
                    xRef = outData[3]
                    xdot = outData[4]
                    freeVariables = outData[5]
                    DF = outData[6]
                    # calculates the null space vector of Jacobian matrix DF
                    nullSpace = Utility.nullspace(DF)
                    # declares and initializes the augmented free variable vector
                    contiFreeVariables = np.array([x_n[0], x_n[2], x_n[4], tau_n])
                    # declares and initializes the augmented constraints vector
                    contiConstraints = np.array([xRef[-1, 1], xRef[-1, 3], xRef[-1, 5], ((contiFreeVariables - freeVariables).T).dot(nullSpace) + stepSize])
                    # calculates corrections to the initial state to meet a defined margin of error
                    GF = np.array([[phi[1, 0], phi[1, 2], phi[1, 4], xdot[1]],
                                   [phi[3, 0], phi[3, 2], phi[3, 4], xdot[3]],
                                   [phi[5, 0], phi[5, 2], phi[5, 4], xdot[5]],
                                   [nullSpace[0], nullSpace[1], nullSpace[2], nullSpace[3]]])
                    xIter = contiFreeVariables - (np.linalg.inv(GF)).dot(contiConstraints)
                    # sets the updated initial condition vector
                    x_n = np.array([xIter[0], 0, xIter[1], 0, xIter[2], 0])
                    # sets T/2 of updated initial conditions
                    tau_n = xIter[3]
                    # sets attributes
                    self.x0 = x_n
                    self.period = 2 * tau_n
                    Orbit.getJacobi(self)
                    self.data = np.array([self.jacobi, self.period, x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5]])
                    # updates stability index of orbit
                    Orbit.getStability(self)
                    print(self.x0)
                    print("        Stability index: %8.4f" % (self.stability))
                self.x0 = lastX
                print(self.x0)
                self.period = lastPeriod
                Orbit.getStability(self)
                stepSize = stepSize * 1.0e-1
                print("FERTIG")


            print("\n        Initial state:                                    -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]" % (
            x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5]))
            print("        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]\n" % (
            tau_n, contiConstraints[0], contiConstraints[1], contiConstraints[2]))
        print("DONE")


    # calculates jacobi constant of orbit and sets result as attribute
    def getJacobi(self):
        if Orbit.error is True:
            return
        r1 = np.sqrt((self.x0[0] + self.system.mu) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        r2 = np.sqrt((self.x0[0] - (1 - self.system.mu)) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        self.jacobi = -1 / 2 * (self.x0[3] ** 2 + self.x0[4] ** 2 + self.x0[5] ** 2) + 2 * (
                1 / 2 * (self.x0[0] ** 2 + self.x0[1] ** 2)
                + (1 - self.system.mu) / r1 + self.system.mu / r2)


    # calculates highest stability index of orbit and sets result as attribute
    def getStability(self):
        if Orbit.error is True:
            return
        # calculates monodromy matrix
        monodromy = Utility.stm(self.x0, self.period, self.system.mu)
        # calculates eigenvalues of monodromy matrix
        eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        maximum = max(abs(eigenvalues))
        self.stability = 1 / 2 * (maximum + 1 / maximum)
        if self.stability < Orbit.stabilityCriteria:
            self.NRHO = True
        else:
            self.NRHO = False


    # plots orbit
    def plot(self):
        if Orbit.error is True:
            return
        if self.x0[0] < 1:
            self.lagrangian = "L1"
        elif self.x0[0] > 1:
            self.lagrangian = "L2"
        else:
            print("Lagrangian type could not be determined.")
            self.lagrangian = None
        Plot.plot(np.array([self.data]), self.system, Orbit.dict, self.lagrangian)

