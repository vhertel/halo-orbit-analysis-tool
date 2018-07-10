"""
File    : Orbit.py
Author  : Victor Hertel
Date    : 28.05.2018

Includes the Orbit Class
"""


# Imports
import numpy as np
from scipy import misc
from scipy.integrate import odeint
import time
import tkinter as tk
from Utility import NumericalMethods, Utility


class Orbit:
    # accuracy of differential corrections method
    ACCURACY = 1.0e-8
    # criteria for NRHOs
    STABILITY_CRITERIA = 5
    # path to the output folder
    dict = "Output/" + time.strftime("%Y-%m-%dT%H.%M.%S") + "/"

    # initializes by adapting input state to initial state of periodic halo orbits and setting attributes
    def __init__(self, x0, fixedValue, system, statusBar, tau=None):
        # dynamical system
        self.system = system
        # stability index
        self.stability = None
        # bool if NRHO or not
        self.NRHO = None
        # prints status updated
        statusBar.insert(tk.INSERT, "   Adapting input state to periodic Halo Orbit...\n")
        statusBar.see(tk.END)

        # stores initial state, period and constraints of halo orbit in 1x8 vector outData
        try:
            outData = NumericalMethods.diffCorrections(x0, self.system.mu, Orbit.ACCURACY, fixedValue, tau=tau, statusBar=statusBar)
        except OverflowError:
            raise OverflowError
        except StopIteration:
            raise StopIteration
        except np.linalg.linalg.LinAlgError:
            raise np.linalg.linalg.LinAlgError
        else:
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
            # prints status updated
        statusBar.insert(tk.INSERT, "      Done\n>>>")
        statusBar.see(tk.END)

    # searches for closest NRHO by using the pseudo-arclength continuation method
    def getClosestNRHO(self, direction=None):
        # method is canceled when error occured
        # calculates highest stability index when attribute not given
        if self.stability is None:
            Orbit.getStability(self)
        # checks whether orbit is aready NRHO
        if self.NRHO is True:
            print("        Halo Orbit is already NRHO.")
        # starts pseudo-arclength continuation method
        else:
            # adjusts z-value for continuation theme
            if direction == 0:
                self.x0[2] = 5.0e-4
            elif direction == 1:
                self.x0[2] = -5.0e-4
            # stepsize
            stepSize = 0.01
            while self.stability > Orbit.STABILITY_CRITERIA:
                # searches for next NRHO using pseudo-arclength continuation method
                outData = NumericalMethods.diffCorrections(self.x0, self.system.mu, Orbit.ACCURACY, fixedValue="None", returnData=True)
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
        r1 = np.sqrt((self.x0[0] + self.system.mu) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        r2 = np.sqrt((self.x0[0] - (1 - self.system.mu)) ** 2 + self.x0[1] ** 2 + self.x0[2] ** 2)
        self.jacobi = -1 / 2 * (self.x0[3] ** 2 + self.x0[4] ** 2 + self.x0[5] ** 2) + 2 * (
                1 / 2 * (self.x0[0] ** 2 + self.x0[1] ** 2)
                + (1 - self.system.mu) / r1 + self.system.mu / r2)

    # calculates highest stability index of orbit and sets result as attribute
    def getStability(self):

        # calculates monodromy matrix
        monodromy = Utility.stm(self.x0, self.period, self.system.mu)
        # calculates eigenvalues of monodromy matrix
        eigenvalues, eigenvectors = np.linalg.eig(monodromy)
        # takes maximum of eigenvalues and calculates stability index
        maximum = max(abs(eigenvalues))
        self.stability = 1 / 2 * (maximum + 1 / maximum)
        # sets attribute NRHO
        if self.stability < Orbit.STABILITY_CRITERIA:
            self.NRHO = True
        else:
            self.NRHO = False

    # calculates invariant stable and unstable manifolds
    def invariantManifolds(self, numberOfPoints, direction):
        # perturbation of state in stable/unstable eigenvector direction
        epsilon = 0.00013007216403660752
        # epsilon = 0.005
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
        if direction == 0:
            if self.lagrangian == "L1" and stableEigenvector[0] > 0:
                stableEigenvector = -stableEigenvector
            if self.lagrangian == "L2" and stableEigenvector[0] < 0:
                stableEigenvector = -stableEigenvector
            if self.lagrangian == "L1" and unstableEigenvector[0] > 0:
                unstableEigenvector = -unstableEigenvector
            if self.lagrangian == "L2" and unstableEigenvector[0] < 0:
                unstableEigenvector = -unstableEigenvector
        elif direction == 1:
            if self.lagrangian == "L1" and stableEigenvector[0] < 0:
                stableEigenvector = -stableEigenvector
            if self.lagrangian == "L2" and stableEigenvector[0] > 0:
                stableEigenvector = -stableEigenvector
            if self.lagrangian == "L1" and unstableEigenvector[0] < 0:
                unstableEigenvector = -unstableEigenvector
            if self.lagrangian == "L2" and unstableEigenvector[0] > 0:
                unstableEigenvector = -unstableEigenvector
        # perturbates states of orbit in direction of stable/unstable eigenvector
        stableManifoldStates[0, :] = self.x0 + epsilon * stableEigenvector
        unstableManifoldStates[0, :] = self.x0 + epsilon * unstableEigenvector
        # perturbates the other points of orbit
        for i in range(1, numberOfPoints):
            # calculates perturbation vector by using STM from t0=0
            stablePerturbationVector = Utility.stm(self.x0, orbitTimes[i], self.system.mu).dot(stableEigenvector)
            unstablePerturbationVector = Utility.stm(self.x0, orbitTimes[i], self.system.mu).dot(unstableEigenvector)
            # normalizing perturbation vectors
            stablePerturbationVector = stablePerturbationVector / np.sqrt(
                stablePerturbationVector[0] ** 2 + stablePerturbationVector[1] ** 2 + stablePerturbationVector[2] ** 2)
            unstablePerturbationVector = unstablePerturbationVector / np.sqrt(
                unstablePerturbationVector[0] ** 2 + unstablePerturbationVector[1] ** 2 + unstablePerturbationVector[
                    2] ** 2)
            # perturbates points of orbit
            stableManifoldStates[i, :] = reducedOrbitStates[i, :] + epsilon * stablePerturbationVector
            unstableManifoldStates[i, :] = reducedOrbitStates[i, :] + epsilon * unstablePerturbationVector
        # sets attributes of orbit
        self.stableManifolds = stableManifoldStates
        self.unstableManifolds = unstableManifoldStates


    @classmethod
    def setAccuracy(cls, accuracy):
        cls.ACCURACY = accuracy
    @classmethod
    def setStabilityCriteria(cls, stabilityCriteria):
        cls.STABILITY_CRITERIA = stabilityCriteria




class InitialGuess:

    def __init__(self, dynamicalSystem, lagrangian, family, fixedValue, value):
        np.seterr(all='raise')

        try:
            # BERECHNUNGEN LAGRANGE POSITION
            l = 1-dynamicalSystem.mu
            if lagrangian == "L1":
                L1 = np.zeros(3)
                p_L1 = np.array([1, 2 * (dynamicalSystem.mu - l), l ** 2 - 4 * l * dynamicalSystem.mu + dynamicalSystem.mu ** 2, 2 * dynamicalSystem.mu * l * (l - dynamicalSystem.mu) + dynamicalSystem.mu - l,
                                 dynamicalSystem.mu ** 2 * l ** 2 + 2 * (l ** 2 + dynamicalSystem.mu ** 2), dynamicalSystem.mu ** 3 - l ** 3])
                L1roots = np.roots(p_L1)
                for i in range(5):
                    if - dynamicalSystem.mu < L1roots[i] < l:
                        L1[0] = np.real(L1roots[i])

                gamma = abs(1-dynamicalSystem.mu-L1[0])
                c2 = 1/gamma**3 * (dynamicalSystem.mu + ((1-dynamicalSystem.mu)*gamma**3)/(1-gamma)**3)
                c3 = 1/gamma**3 * (dynamicalSystem.mu - ((1-dynamicalSystem.mu)*gamma**4)/(1-gamma)**4)
                c4 = 1/gamma**3 * (dynamicalSystem.mu + ((1-dynamicalSystem.mu)*gamma**5)/(1-gamma)**5)

                F = np.array([[1-dynamicalSystem.mu-L1[0], 0, 0],
                  [0, 1-dynamicalSystem.mu-L1[0], 0],
                  [0, 0, 1-dynamicalSystem.mu-L1[0]]])
                P = np.array([L1[0], L1[1], L1[2]])

            elif lagrangian == "L2":
                L2 = np.zeros(3)
                p_L2 = np.array([1, 2 * (dynamicalSystem.mu - l), l ** 2 - 4 * l * dynamicalSystem.mu + dynamicalSystem.mu ** 2, 2 * dynamicalSystem.mu * l * (l - dynamicalSystem.mu) - (dynamicalSystem.mu + l),
                                 dynamicalSystem.mu ** 2 * l ** 2 + 2 * (l ** 2 - dynamicalSystem.mu ** 2), -(dynamicalSystem.mu ** 3 + l ** 3)])
                L2roots = np.roots(p_L2)
                for i in range(5):
                    if L2roots[i] > - dynamicalSystem.mu and L2roots[i] > l:
                        L2[0] = np.real(L2roots[i])

                gamma = abs(L2[0]-(1-dynamicalSystem.mu))
                c2 = 1/gamma**3 * (dynamicalSystem.mu + ((1-dynamicalSystem.mu)*gamma**3)/(1+gamma)**3)
                c3 = 1/gamma**3 * (-dynamicalSystem.mu - ((1-dynamicalSystem.mu)*gamma**4)/(1+gamma)**4)
                c4 = 1/gamma**3 * (dynamicalSystem.mu + ((1-dynamicalSystem.mu)*gamma**5)/(1+gamma)**5)

                F = np.array([[L2[0]-(1-dynamicalSystem.mu), 0, 0],
                  [0, L2[0]-(1-dynamicalSystem.mu), 0],
                  [0, 0, L2[0]-(1-dynamicalSystem.mu)]])
                P = np.array([L2[0], L2[1], L2[2]])


            # ANALYTISCHE BERECHNUNGEN
            coeff = [1, 0, (c2-2), 0, -(c2-1)*(1+2*c2)]
            res = np.roots(coeff)
            lam = abs(res[0])
            k = 1/(2*lam) * (lam**2 + 1 + 2*c2)
            d1 = (3*lam**2)/(k) * (k*(6*lam**2-1) - 2*lam)
            d2 = (8*lam**2)/(k) * (k*(11*lam**2-1) - 2*lam)
            a21 = (3*c3*(k**2-2))/(4*(1+2*c2))
            a22 = (3*c3)/(4*(1+2*c2))
            a23 = - (3*c3*lam)/(4*k*d1) * (3*k**3*lam - 6*k*(k-lam) + 4)
            a24 = - (3*c3*lam)/(4*k*d1) * (2 + 3*k*lam)
            b21 = - (3*c3*lam)/(2*d1) * (3*k*lam - 4)
            b22 = (3*c3*lam)/(d1)
            d21 = - (c3)/(2*lam**2)
            d31 = 3/(64*lam**2) * (4*c3*a24 + c4)
            d32 = 3/(64*lam**2) * (4*c3*(a23-d21) + c4*(4+k**2))
            a31 = - (9*lam)/(4*d2) * (4*c3*(k*a23 - b21) + k*c4*(4 + k**2)) + (9*lam**2 + 1 - c2)/(2*d2) * (3*c3*(2*a23 - k*b21) + c4*(2 + 3*k**2))
            a32 = - (9*lam)/(4*d2) * (4*c3*(k*a24 - b22) + k*c4) - (3*(9*lam**2 + 1 - c2))/(2*d2) * (c3*(k*b22 + d21 - 2*a24) - c4)
            b31 = (3*lam)/(d2) * (3*c3*(k*b21 - 2*a23) - c4*(2+3*k**2)) + (3*(9*lam**2+1+2*c2))/(8*d2) * (4*c3*(k*a23 - b21) + k*c4*(4+k**2))
            b32 = (9*lam)/(d2) * (c3*(k*b22 + d21 - 2*a24) - c4) + (3*(9*lam**2+1+2*c2))/(8*d2) * (4*c3*(k*a24 - b22) + k*c4)
            s1 = (3/2*c3*(2*a21*(k**2-2) - a23*(k**2+2) - 2*k*b21) - 3/8*c4*(3*k**4 - 8*k**2 + 8))/(2*lam*(lam*(1+k**2) - 2*k))
            s2 = (3/2*c3*(2*a22*(k**2-2) + a24*(k**2+2) + 2*k*b22 + 5*d21) + 3/8*c4*(12-k**2))/(2*lam*(lam*(1+k**2) - 2*k))
            l1 = 2*lam**2*s1 - 3/2*c3*(2*a21 + a23 + 5*d21) - 3/8*c4*(12-k**2)
            l2 = 2*lam**2*s2 + 3/2*c3*(a24 - 2*a22) + 9/8*c4
            delta = lam**2 - c2


            if fixedValue == "x":
                ax = value
                az = None
            elif fixedValue == "z":
                ax = None
                az = value
            elif fixedValue == "Period":
                period = value / (np.sqrt(dynamicalSystem.distance ** 3 / (dynamicalSystem.G * (dynamicalSystem.massFP + dynamicalSystem.massSP))) / (60 * 60 * 24))
                ax = None
                az = np.sqrt((delta*s1-l1) * lam*period + 2*l1*np.pi)/np.sqrt((l1*s2 - l2*s1) * lam * period)


            # EINGABE ORBIT
            if ax is None:
                ax = np.sqrt(-1/l1 * (l2*az**2 + delta))
            elif az is None:
                az = np.sqrt(-1/l2 * (l1*ax**2 + delta))
            elif l1*ax**2 + l2*az**2 + delta > 0.1:
                print("Combination of amplitudes is not possible.")
                exit()
            tau = 0
            phi = 0
            if (family == "Northern" and lagrangian == "L1") or (family == "Southern" and lagrangian == "L2"):
                deltan = 1
            else:
                deltan = -1
            omega = 1 + s1*ax**2 + s2*az**2

            # BERECHNUNGEN
            def xi(tau):
                return a21*ax**2 + a22*az**2 - ax*np.cos(lam*tau + phi) + (a23*ax**2 - a24*az**2)*np.cos(2*lam*tau + 2*phi) + (a31*ax**3 - a32*ax*az**2)*np.cos(3*lam*tau + 3*phi)
            def eta(tau):
                return k*ax*np.sin(lam*tau + phi) + (b21*ax**2 - b22*az**2)*np.sin(2*lam*tau + 2*phi) + (b31*ax**3 - b32*ax*az**2)*np.sin(3*lam*tau + 3*phi)
            def zeta(tau):
                return deltan*az*np.cos(lam*tau + phi) + deltan*d21*ax*az*(np.cos(2*lam*tau + 2*phi) - 3) + deltan*(d32*az*ax**2 - d31*az**3)*np.cos(3*lam*tau + 3*phi)

            # AUSWERTUNG
            x = xi(tau)
            y = eta(tau)
            z = zeta(tau)
            xdot = misc.derivative(xi, tau, dx=1e-15)
            ydot = misc.derivative(eta, tau, dx=1e-15)
            zdot = misc.derivative(zeta, tau, dx=1e-15)
            x0 = np.array([x, y, z, xdot, ydot, zdot])

            x0[0:3] = F.dot(x0[0:3]) + P
            x0[3:6] = F.dot(x0[3:6])

            self.x0 = x0
            self.tau = (2*np.pi)/(lam*omega)/2

        except:
            raise
