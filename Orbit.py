"""
File    : Orbit.py
Author  : Victor Hertel
Date    : 28.05.2018

Includes the Orbit Class
"""

import time

# Imports
import numpy as np

from Utility import NumericalMethods
from Utility import Plot
from Utility import Utility

import matplotlib as mpl

mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import misc
from Utility import Utility, Plot
from scipy.integrate import odeint


class Orbit:
    # accuracy of differential corrections method
    accuracy = 1.0e-8
    # criteria for NRHOs
    stabilityCriteria = 5
    # indicates error during initialization of an orbit
    error = False
    # path to the output folder
    dict = "Output/" + time.strftime("%Y-%m-%dT%H.%M.%S") + "/"

    # initializes by adapting input state to initial state of periodic halo orbits and setting attributes
    def __init__(self, x0, fixedValue, system, comment="True"):
        # dynamical system
        self.system = system
        # stability index
        self.stability = None
        # bool if NRHO or not
        self.NRHO = None
        # prints status updated if requested
        if comment:
            print("STATUS: Adapting input state to periodic Halo Orbit...\n")
        # stores initial state, period and constraints of halo orbit in 1x8 vector outData
        try:
            outData = NumericalMethods.diffCorrections(x0, self.system.mu, Orbit.accuracy, fixedValue=fixedValue)
        except ValueError:
            Orbit.error = True
            return
        # initial state
        self.x0 = outData[0:6]
        # period
        self.period = outData[6]
        # calculates and sets jacobi constant
        Orbit.getJacobi(self)
        # sets data for plot
        self.data = np.array(
            [self.jacobi, self.period, self.x0[0], self.x0[1], self.x0[2], self.x0[3], self.x0[4], self.x0[5]])
        # checks for lagrangian
        if self.x0[0] < 1:
            self.lagrangian = "L1"
        elif self.x0[0] > 1:
            self.lagrangian = "L2"
        self.unstableManifolds = None
        self.stableManifolds = None
        # prints status updated if requested
        if comment:
            print("DONE")

    # searches for closest NRHO by using the pseudo-arclength continuation method
    def getClosestNRHO(self, direction=None):
        # method is canceled when error occured
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
        # starts pseudo-arclength continuation method
        else:
            # adjusts z-value for continuation theme
            print(direction)
            if direction == 0:
                self.x0[2] = 5.0e-4
            elif direction == 1:
                self.x0[2] = -5.0e-4
            # calculates T/2 by integrating until y changes sign
            try:
                tau_n = Utility.halfPeriod(self.x0, self.system.mu, 1.0e-11)
            except ValueError:
                return
            # stepsize
            stepSize = 0.01
            while self.stability > Orbit.stabilityCriteria:
                # searches for next NRHO using pseudo-arclength continuation method
                outData = NumericalMethods.diffCorrections(self.x0, self.system.mu, Orbit.accuracy, tau=tau_n)
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
                contiFreeVariables = np.array([x_n[0], x_n[2], x_n[4], tau_n])
                # declares and initializes the augmented constraints vector
                contiConstraints = np.array([xRef[-1, 1], xRef[-1, 3], xRef[-1, 5],
                                             (contiFreeVariables - freeVariables).T.dot(nullSpace) + stepSize])
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
                lastStability = self.stability
                Orbit.getStability(self)
                print("        Stability index: %8.4f" % self.stability)
            self.stableManifolds = None
            self.unstableManifolds = None
            print(
                "\n        Initial state:                                    -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]" % (
                    x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5]))
            print(
                "        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]\n" % (
                    tau_n, contiConstraints[0], contiConstraints[1], contiConstraints[2]))
        print("DONE")

    # calculates jacobi constant of orbit and sets result as attribute
    def getJacobi(self):
        # method is canceled when error occured
        if Orbit.error is True:
            return
        r1 = np.sqrt((self.x0[0] + self.system.mu) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        r2 = np.sqrt((self.x0[0] - (1 - self.system.mu)) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        self.jacobi = -1 / 2 * (self.x0[3] ** 2 + self.x0[4] ** 2 + self.x0[5] ** 2) + 2 * (
                1 / 2 * (self.x0[0] ** 2 + self.x0[1] ** 2)
                + (1 - self.system.mu) / r1 + self.system.mu / r2)

    # calculates highest stability index of orbit and sets result as attribute
    def getStability(self):
        # method is canceled when error occured
        if Orbit.error is True:
            return
        # calculates monodromy matrix
        monodromy = Utility.stm(self.x0, self.period, self.system.mu)
        # calculates eigenvalues of monodromy matrix
        eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        # takes maximum of eigenvalues and calculates stability index
        maximum = max(abs(eigenvalues))
        self.stability = 1 / 2 * (maximum + 1 / maximum)
        # sets attribute NRHO
        if self.stability < Orbit.stabilityCriteria:
            self.NRHO = True
        else:
            self.NRHO = False

    def invariantManifolds(self, numberOfPoints):
        # perturbation of state in stable/unstable eigenvector direction
        epsilon = 0.00013007216403660752
        # number of points to split the orbit
        if numberOfPoints <= 1000:
            num = 1000
        else:
            num = numberOfPoints
        # integrates orbit in CR3BP
        t = np.linspace(0, self.period, num=num)
        orbitStates = odeint(Utility.sysEquations, self.x0, t, args=(self.system.mu,), rtol=2.5e-13, atol=1e-22)
        # specifies distance to get uniformed manifolds around the orbit
        if numberOfPoints == 0:
            return
        else:
            orbitTags = len(orbitStates) / numberOfPoints
        # declares state vector with desired number of manifolds
        reducedOrbitStates = np.zeros((numberOfPoints, 6))
        # declares time vector with desired number of manifolds
        orbitTimes = np.zeros(numberOfPoints)
        # declares manifold states with desired number of manifolds
        stableManifoldStates = np.zeros((numberOfPoints, 6))
        unstableManifoldStates = np.zeros((numberOfPoints, 6))
        # gets states and time of points on the orbit
        for i in range(numberOfPoints):
            reducedOrbitStates[i, :] = orbitStates[round(i * orbitTags), :]
            orbitTimes[i] = t[round(i * orbitTags)]
        # calculates monodromy matrix of orbit
        monodromy = Utility.stm(self.x0, self.period, self.system.mu)
        # gets eigenvalues and eigenvectors of monodromy matrix
        eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        # stores stable and unstable eigenvector
        i = 0
        for element in eigenvalues:
            if element == min(eigenvalues):
                stableEigenvector = np.real(eigenvectors[:, i])
            if element == max(eigenvalues):
                unstableEigenvector = np.real(eigenvectors[:, i])
            i += 1
        # perturbates states of orbit in direction of stable/unstable eigenvector
        stableManifoldStates[0, :] = self.x0 + epsilon * stableEigenvector
        unstableManifoldStates[0, :] = self.x0 + epsilon * unstableEigenvector
        # perturbates the other points of orbit
        for i in range(1, numberOfPoints):
            # calculates perturbation vector by using STM from t0=0
            stablePerturbationVector = Utility.stm(self.x0, orbitTimes[i], self.system.mu).dot(stableEigenvector)
            unstablePerturbationVector = Utility.stm(self.x0, orbitTimes[i], self.system.mu).dot(unstableEigenvector)
            # normalizing perturbation vectors
            stablePerturbationVector = stablePerturbationVector/np.sqrt(stablePerturbationVector[0]** 2 + stablePerturbationVector[1]** 2 + stablePerturbationVector[2]** 2)
            unstablePerturbationVector = unstablePerturbationVector/np.sqrt(unstablePerturbationVector[0]** 2 + unstablePerturbationVector[1]** 2 + unstablePerturbationVector[2]** 2)
            # perturbates points of orbit
            stableManifoldStates[i, :] = reducedOrbitStates[i, :] + epsilon * stablePerturbationVector
            unstableManifoldStates[i, :] = reducedOrbitStates[i, :] + epsilon * unstablePerturbationVector
        # sets attributes of orbit
        self.stableManifolds = stableManifoldStates
        self.unstableManifolds = unstableManifoldStates

    # plots orbit
    def plot(self):
        # method is canceled when error occured
        if Orbit.error is True:
            return
        Plot.plot(np.array([self.data]), self.system, Orbit.dict)
