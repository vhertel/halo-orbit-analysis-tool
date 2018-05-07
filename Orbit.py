"""
File    : Orbit.py
Author  : Victor Hertel
Date    : 24.04.2018

Orbit Class and NRHO Subclass
"""

# Imports
import numpy as np

from Utility import NumericalMethods
from Utility import Plot
from Utility import Utility


# orbit class
class Orbit:
    # global variable accuracy
    accuracy = 1.0e-6

    # initializes by adapting input state to initial state of periodic halo orbits
    def __init__(self, x0, fixedValue, mu, jacobi=None, stability=None, NRHO=None):
        # prints status update
#        print("STATUS: Adapting input state to periodic Halo Orbit...\n")
        self.mu = mu
        # stores initial state, period and constraints of halo orbit in 1x3 vector outData
        try:
            outData = NumericalMethods.diffCorrections(x0, self.mu, Orbit.accuracy, fixedValue=fixedValue)
        except ValueError:
            return
        # sets attributes of halo orbit
        self.x0 = outData[0:6]
        self.period = outData[6]
        self.jacobi = jacobi
        self.stability = stability
        self.NRHO = NRHO
        self.data = outData
#        print("DONE")

    def getNearestNRHO(self):
        # prints status update
        print("STATUS: Adapting input state to initial state of nearest NRHO...")
        # calculates highest stability index when attribute not given
        if self.stability is None:
            Orbit.getStability()
        # checks whether orbit is aready NRHO
        if self.NRHO is True:
            print("        Halo Orbit is already NRHO.")
        else:
            # calculates T/2 by integrating until y changes sign
            try:
                tau_n = Utility.halfPeriod(self.x0, self.mu, 1.0e-11)
            except ValueError:
                return
            # searches for next NRHO using pseudo-arclength continuation method
            while self.NRHO is False:
                outData = NumericalMethods.diffCorrections(self.x0, self.mu, Orbit.accuracy, tau=tau_n)
                x_n = outData[0]
                tau_n = outData[1]
                phi = outData[2]
                xRef = outData[3]
                xdot = outData[4]
                freeVariables = outData[5]
                DF = outData[6]
                # calculates the null space vector of Jacobian matrix DF
                nullSpace = Utility.nullspace(DF)
                # declares and initializes the augmented free variable vector
                contiFreeVariables = np.array([[x_n[0]],
                                               [x_n[2]],
                                               [x_n[4]],
                                               [tau_n]])
                # declares and initializes the augmented constraints vector
                contiConstraints = np.array([[xRef[2000, 1]],
                                             [xRef[2000, 3]],
                                             [xRef[2000, 5]],
                                             [((contiFreeVariables - freeVariables).T).dot(nullSpace) + 0.02]])
                # calculates corrections to the initial state to meet a defined margin of error
                GF = np.array([[phi[1, 0], phi[1, 2], phi[1, 4], xdot[1]],
                               [phi[3, 0], phi[3, 2], phi[3, 4], xdot[3]],
                               [phi[5, 0], phi[5, 2], phi[5, 4], xdot[5]],
                               [nullSpace[0], nullSpace[1], nullSpace[2], nullSpace[3]]])
                xIter = contiFreeVariables - (np.linalg.inv(GF)).dot(contiConstraints)
                # sets the updated initial condition vector
                x_n = np.array([xIter[0, 0], 0, xIter[1, 0], 0, xIter[2, 0], 0])
                # sets T/2 of updated initial conditions
                tau_n = xIter[3, 0]
                # sets attributes
                self.x0 = x_n
                self.period = 2 * tau_n
                self.data = np.array([x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5], 2 * tau_n])
                # updates stability index of orbit
                Orbit.getStability(self)
                print(self.stability)
                print(self.x0)
            print("        Initial state:                                    -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]" % (
            x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5]))
            print("        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]" % (
            tau_n, contiConstraints[0], contiConstraints[1], contiConstraints[2]))
        print("DONE")

    # calculates jacobi constant of orbit and sets result as attribute
    def getJacobi(self):
        r1 = np.sqrt((self.x0[0] + self.mu) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        r2 = np.sqrt((self.x0[0] - (1 - self.mu)) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        self.jacobi = -1 / 2 * (self.x0[3] ** 2 + self.x0[4] ** 2 + self.x0[5] ** 2) + 2 * (
                1 / 2 * (self.x0[0] ** 2 + self.x0[1] ** 2)
                + (1 - self.mu) / r1 + self.mu / r2)

    # calculates highest stability index of orbit and sets result as attribute
    def getStability(self):
        # calculates monodromy matrix
        monodromy = Utility.stm(self.x0, self.period, self.mu)
        # calculates eigenvalues of monodromy matrix
        eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        maximum = max(abs(eigenvalues))
        self.stability = 1 / 2 * (maximum + 1 / maximum)
        if self.stability < 5:
            self.NRHO = True
        else:
            self.NRHO = False

    # plots orbit
    def plot(self, haloFamily="both", background="off"):
        if self.jacobi is None:
            Orbit.getJacobi(self)
        Plot.plot(np.array([self.data]), self.mu, haloFamily, background)

    # sets accuracy
    @classmethod
    def setAccuracy(cls, accuracy):
        cls.accuracy = accuracy
