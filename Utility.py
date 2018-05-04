"""
File    : Utility.py
Author  : Victor Hertel
Date    : 24.04.2018

Utility classes
"""

# Imports
import math
import matplotlib as mpl
import numpy as np
from numpy.linalg import svd
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

        output = np.zeros((orbitNumber, 7))

        # sets the iteratively adjusted initial condition
        x_n = x0
        lastX = x_n
        outX = x_n
        # loops through number of orbits
        for i in range(orbitNumber):
            # checks whether x- or z-value has changed more since last iteration
            if abs(lastX[0] - outX[0]) > abs(lastX[2] - outX[2]):
                # x-value changed more than z-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed x-value)\n" % (i + 1))
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
                try:
                    outData = NumericalMethods.diffCorrections(x_n, mu, OrbitContinuation.accuracy, fixedValue="x")
                except ValueError:
                    print("\nDONE")
                    return
                outX = outData[0:6]
                output[i][:] = outData
                # output[i][0:6] = outData[0]
                # output[i][6] = outData[1]

            else:
                # z-value changed more than x-value and needs to be fixed
                print("        Orbit Number: %2d    (fixed z-value)\n" % (i + 1))
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
                try:
                    outData = NumericalMethods.diffCorrections(x_n, mu, OrbitContinuation.accuracy, fixedValue="z")
                except ValueError:
                    print("\nDONE")
                    return
                outX = outData[0:6]
                output[i][:] = outData
                # output[i][0:6] = outData[0]
                # output[i][6] = outData[1]

        return output

    # --------------------------------------------------------------------------
    # PSEUDO-ARCLENGTH CONTINUATION
    # DESCRIPTION:      The Pseudo-Arclength Continuation method is used for finding a family of
    #                   orbits e.g. a family of related solution. After determining one solution,
    #                   the method steps in a tangent direction of the free variables vector to
    #                   generate an initial guess for the next orbit, which is then corrected
    #                   using the differential corrections method.
    #
    # OUTPUT:           No output is generated. The orbit is plotted in the
    #                   figure that has been created before the function call
    #
    # SYNTAX:           pseudoArcLenConti(x0, t0, mu, epsilon, orbitNumber, familyStep, direction, haloFamily, ax)
    # x0            =   Initial State Guess
    # t0            =   Time Start
    # mu            =   Mass ratio of Primaries
    # epsilon       =   Error Tolerance of Constraints at T/2
    # orbitNumber   =   Number of Orbits to search for
    # familyStep    =   Stepsize in x- or z-Direction
    # haloFamily    =   Desired family of Halo Orbits:
    #                   Input possibilities: {"northern", "southern", "both"}
    # ax            =   Allows access to the figure
    #   x0, orbitDistance, lagrangian, orbitNumber, mu

    # --------------------------------------------------------------------------

    @classmethod
    def setAccuracy(cls, accuracy):
        cls.accuracy = accuracy


class NumericalMethods:

    # --------------------------------------------------------------------------
    # FOURTH ODER RUNGE-KUTTA METHOD
    #
    # DESCRIPTION:      The function rk4System() approximates the solutions of a system of m
    #                   differential equations that are written in the form
    #
    #                   dy1/dt = f1(t,y1,y2,...,ym)
    #                   dy2/dt = f2(t,y1,y2,...,ym)
    #                   ...
    #                   dym/dt = fm(t,y1,y2,...,ym)
    #
    #                   with t in the interval [0; tf] and the m-dimensional initial conditions
    #                   vector y0. The accuracy depends on the step size n.
    #
    # OUTPUT:           A matrix consisting of a row for each time and 42 columns including the elements
    #                   of the State Transition Matrix and the state vector dx/dt
    #
    # SYNTAX:           rk4System(y, tf, mu)
    # y             =   Position values of state vector
    # tf            =   Final time
    # mu            =   Mass ratio of Primaries
    # --------------------------------------------------------------------------

    @staticmethod
    def rk4System(y, tf, mu):

        # defines step size depending on time interval [0; tf] and n
        n = 2000
        h = tf / n
        # declares and initializes vector t with dimensions (n+1)
        t = np.zeros(n + 1)
        t[0] = 0
        # declares and initializes matrix w with a row for each time step and 42 columns
        R = np.zeros((len(y), n + 1))
        for i in range(len(y)):
            R[i, 0] = y[i]

        # 4th order Runge-Kutta method algorithm
        if len(y) == 42:

            for i in range(n):
                k1 = Utility.sysEquations(R[:, i], mu, t[i])
                k2 = Utility.sysEquations(R[:, i] + h / 2 * k1[:, 0], mu, t[i] + h / 2)
                k3 = Utility.sysEquations(R[:, i] + h / 2 * k2[:, 0], mu, t[i] + h / 2)
                k4 = Utility.sysEquations(R[:, i] + h * k3[:, 0], mu, t[i] + h)
                R[:, i + 1] = R[:, i] + h / 6 * (k1[:, 0] + 2 * k2[:, 0] + 2 * k3[:, 0] + k4[:, 0])
                t[i + 1] = t[i] + h
            R = R.T

        elif len(y) == 6:

            for i in range(n):
                k1 = Utility.sysEquations(R[:, i], mu, t[i])
                k2 = Utility.sysEquations(R[:, i] + h / 2 * k1, mu, t[i] + h / 2)
                k3 = Utility.sysEquations(R[:, i] + h / 2 * k2, mu, t[i] + h / 2)
                k4 = Utility.sysEquations(R[:, i] + h * k3, mu, t[i] + h)
                R[:, i + 1] = R[:, i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                t[i + 1] = t[i] + h
            R = R.T

        else:
            print("Dimension of input parameter y is not supported.")

        return R

    # --------------------------------------------------------------------------
    # DIFFERENTIAL CORRECTIONS METHOD
    #
    #
    # --------------------------------------------------------------------------

    @staticmethod
    def diffCorrections(x, mu, epsilon, tau=None, fixedValue=None):

        if fixedValue is None and tau is not None:
            # declares and initializes constraint vector
            constraints = np.ones((3, 1))
            # prints status update
            #            print("        Differential Corrections Method for adaption of the initial state...")
            # declares and initializes counter for counting updates of initial state
            counter = 1
            # constraints are corrected until they meet a defined margin of error
            while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon or abs(constraints[2]) > epsilon:
                if counter > 15:
                    print("        Differential Corrections Method did not converge.")
                    raise ValueError
                else:
                    counter = counter + 1
                # calculates the state transition matrix
                phi = Utility.stm(x, tau, mu)
                # integrates initial state from t0 to tau_n
                xRef = NumericalMethods.rk4System(x, tau, mu)
                # calculates the derivation of the state at T/2
                xdot = Utility.sysEquations(xRef[2000, :], mu)
                # declares and initializes free variable vector with x, z, ydot and tau
                freeVariables = np.array([[x[0]],
                                          [x[2]],
                                          [x[4]],
                                          [tau]])
                # declares and initializes the constraint vector with y, xdot and zdot
                constraints = np.array([[xRef[2000, 1]],
                                        [xRef[2000, 3]],
                                        [xRef[2000, 5]]])
                # calculates corrections to the initial state to meet a defined margin of error
                DF = np.array([[phi[1, 0], phi[1, 2], phi[1, 4], xdot[1]],
                               [phi[3, 0], phi[3, 2], phi[3, 4], xdot[3]],
                               [phi[5, 0], phi[5, 2], phi[5, 4], xdot[5]]])
                xIter = freeVariables - ((DF.T).dot(np.linalg.inv(DF.dot(DF.T)))).dot(constraints)
                # sets the updated initial condition vector
                x = np.array([xIter[0, 0], 0, xIter[1, 0], 0, xIter[2, 0], 0])
                # sets T/2 of updated initial conditions
                tau = xIter[3, 0]
                counter = counter + 1
            # prints status update
            #            print("        Initial state has been adapted for %d times:       -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]" % (counter, x[0], x[1], x[2], x[3], x[4], x[5]))
            #            print("        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]" % (tau, constraints[0], constraints[1], constraints[2]))
            outX = x
            outTime = tau
            outData = np.array([outX, outTime, phi, xRef, xdot, freeVariables, DF])

            return outData


        elif fixedValue is not None and tau is None:

            print("        Differential Corrections Method for adaption of the initial state...")
            # calculates T/2 by integrating until y changes sign
            try:
                tau = Utility.halfPeriod(x, mu, 1.0e-11)
            except ValueError:
                return

            # declares and initializes constraint vector
            constraints = np.ones((2, 1))
            # declares and initializes counter
            counter = 0

            # case of x-amplitude being the fixed variable
            if fixedValue == "x":

                # constraints are corrected until they meet a defined margin of error
                while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon:
                    if counter > 15:
                        print("        Differential Corrections Method did not converge.")
                        raise ValueError
                    else:
                        counter = counter + 1
                    # calculates the state transition matrix
                    phi = Utility.stm(x, tau, mu)
                    # integrates initial state from 0 to tau_n
                    xRef = NumericalMethods.rk4System(x, tau, mu)
                    # calculates the derivation of the state at T/2
                    xdot = Utility.sysEquations(xRef[2000, :], mu)
                    # declares and initializes free variable vector with z, ydot and time
                    freeVariables = np.array([[x[2]],
                                              [x[4]],
                                              [tau]])
                    # declares and initializes the constraint vector with y, xdot and zdot
                    constraints = np.array([[xRef[2000, 1]],
                                            [xRef[2000, 3]],
                                            [xRef[2000, 5]]])
                    # calculates corrections to the initial state to meet a defined margin of error
                    D = np.array([[phi[1, 2], phi[1, 4], xdot[1]],
                                  [phi[3, 2], phi[3, 4], xdot[3]],
                                  [phi[5, 2], phi[5, 4], xdot[5]]])
                    DF = np.linalg.inv(D)
                    xIter = freeVariables - DF.dot(constraints)
                    # sets the updated initial condition vector
                    x = np.array([x[0], 0, xIter[0, 0], 0, xIter[1, 0], 0])
                    # calculates T/2 with updated initial conditions
                    tau = xIter[2, 0]
                    outX = x
                    outTime = tau

            # case of z-amplitude being the fixed variable
            elif fixedValue == "z":

                # constraints are corrected until they meet a defined margin of error
                while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon:
                    if counter > 15:
                        print("        Differential Corrections Method did not converge.")
                        raise ValueError
                    else:
                        counter = counter + 1
                    # calculates the state transition matrix
                    phi = Utility.stm(x, tau, mu)
                    # integrates initial state from 0 to tau_n
                    xRef = NumericalMethods.rk4System(x, tau, mu)
                    # calculates the derivation of the state at T/2
                    xdot = Utility.sysEquations(xRef[2000, :], mu)
                    # declares and initializes free variable vector with x, ydot and time
                    freeVariables = np.array([[x[0]],
                                              [x[4]],
                                              [tau]])
                    # declares and initializes the constraint vector with y, xdot and zdot
                    constraints = np.array([[xRef[2000, 1]],
                                            [xRef[2000, 3]],
                                            [xRef[2000, 5]]])
                    # calculates corrections to the initial state to meet a defined margin of error
                    D = np.array([[phi[1, 0], phi[1, 4], xdot[1]],
                                  [phi[3, 0], phi[3, 4], xdot[3]],
                                  [phi[5, 0], phi[5, 4], xdot[5]]])
                    DF = np.linalg.inv(D)
                    xIter = freeVariables - DF.dot(constraints)
                    # sets the updated initial condition vector
                    x = np.array([xIter[0, 0], 0, x[2], 0, xIter[1, 0], 0])
                    # calculates T/2 with updated initial conditions
                    tau = xIter[2, 0]
                    outX = x
                    outTime = tau

            # stores the initial condition and T/2 in output vector
            outData = np.array([outX[0], outX[1], outX[2], outX[3], outX[4], outX[5], 2 * outTime])
            # prints status updates
            print("        Initial state has been adapted for %d times:       -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]"
                  % (counter, outData[0], outData[1], outData[2], outData[3], outData[4], outData[5]))
            print(
                "        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]\n"
                % (1 / 2 * outData[6], constraints[0], constraints[1], constraints[2]))

            return outData

        else:
            print("Input parameters not correct.")


class Utility:

    # --------------------------------------------------------------------------
    # CALCULATING MATRIX A
    #
    # DESCRIPTION:      A Taylor series expansion retaining only first-order terms is used to
    #                   linearize the nonlinear system equations of motion. With the six-dimensional
    #                   state vector x = [x, y, z, vx, vy, vz]^T, this system of three second-order
    #                   differential equations can be written in state space form as
    #
    #                        dx/dt = A(t) * x(t)
    #
    # OUTPUT:           Output is the 6x6 matrix A(t)
    #
    # SYNTAX:           AMatrix(x, mu)
    # x             =   Position values of state vector
    # mu            =   Mass ratio of Primaries
    # --------------------------------------------------------------------------

    @staticmethod
    def AMatrix(x, mu):

        # declares and initializes 3x3 submatrix O with zeros
        O = np.zeros((3, 3))
        # declares and initializes 3x3 unit submatrix I
        I = np.eye(3)
        # calculates r1 and r2
        r1 = np.sqrt((x[0] + mu) ** 2 + x[1] ** 2 + x[2] ** 2)
        r2 = np.sqrt((x[0] - (1 - mu)) ** 2 + x[1] ** 2 + x[2] ** 2)
        # calculates the elements of the second partial derivatives of the three-body pseudo-potential U
        Uxx = 1 - (1 - mu) / (r1 ** 3) - (mu) / (r2 ** 3) + (3 * (1 - mu) * (x[0] + mu) ** 2) / (r1 ** 5) + (
                3 * mu * (x[0] - (1 - mu)) ** 2) / (r2 ** 5)
        Uxy = (3 * (1 - mu) * (x[0] + mu) * x[1]) / (r1 ** 5) + (3 * mu * (x[0] - (1 - mu)) * x[1]) / (r2 ** 5)
        Uxz = (3 * (1 - mu) * (x[0] + mu) * x[2]) / (r1 ** 5) + (3 * mu * (x[0] - (1 - mu)) * x[2]) / (r2 ** 5)
        Uyx = Uxy
        Uyy = 1 - (1 - mu) / (r1 ** 3) - (mu) / (r2 ** 3) + (3 * (1 - mu) * x[1] ** 2) / (r1 ** 5) + (
                3 * mu * x[1] ** 2) / (r2 ** 5)
        Uyz = (3 * (1 - mu) * x[1] * x[2]) / (r1 ** 5) + (3 * mu * x[1] * x[2]) / (r2 ** 5)
        Uzx = Uxz
        Uzy = Uyz
        Uzz = - (1 - mu) / (r1 ** 3) - (mu) / (r2 ** 3) + (3 * (1 - mu) * x[2] ** 2) / (r1 ** 5) + (
                3 * mu * x[2] ** 2) / (r2 ** 5)
        # declares and fills 3x3 submatrix U
        U = np.array([[Uxx, Uxy, Uxz],
                      [Uyx, Uyy, Uyz],
                      [Uzx, Uzy, Uzz]])
        # declares and initialized 3x3 submatrix C
        C = np.array([[0, 2, 0],
                      [-2, 0, 0],
                      [0, 0, 0]])
        # assembles 6x6 augmented matrix A with submatrices O, I, U and C
        A = np.bmat("O, I;"
                    "U, C")

        return A

    # --------------------------------------------------------------------------
    # CALCULATING SYSTEMS OF EQUATIONS
    #
    # DESCRIPTION:      Function deals with the equations of motion comprising the dynamical model of the CR3BP.
    #                       case a):    input variable y has dimension of 6:
    #                                   The following state vector is calculated
    #                                   dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
    #                       case b):    input variable y has dimension of 42:
    #                                   The matrix multiplication dPhi(t,0)/dt = A(t) * Phi(t,0) for
    #                                   calculating the State Transition Matrix STM is done at t = 0 with
    #                                   Phi(0,0) = I6. Additionally the state vector
    #                                   dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T is
    #                                   calculated and the values stored in the vector APhi with 42 elements.
    #
    # OUTPUT:           case a):    input variable y has dimension of 6:
    #                               Six dimensional state vector dx/dt
    #                   case b):    input variable y has dimension of 42:
    #                               Vector with 42 elements including the result of the matrix multiplication
    #                               dPhi(t,0)/dt = A(t) * Phi(t,0) as the first 36 values and the six
    #                               dimensional state vector dx/dt as the last 6 elements
    #
    # SYNTAX:           sysEquations(y, mu, t="0")
    # y             =   Position values of state vector
    # mu            =   Mass ratio of Primaries
    # t="0"         =   Parameter is only used for numerical integration with rk4System()
    # --------------------------------------------------------------------------

    @staticmethod
    def sysEquations(y, mu, t="0"):

        if len(y) == 42:
            # declares and initializes vector x including the initial position
            x = np.zeros(3)
            # stores initial position in vector x
            x = y[36:39]
            # matrix A is calculated
            A = Utility.AMatrix(x, mu)
            # declares and initializes 6x6 matrix phi and fills it with components of vector y
            phi = np.zeros((6, 6))
            for i in range(6):
                for j in range(6):
                    phi[i, j] = y[6 * i + j]
            # matrix caluclation of 6x6 matrices A and phi
            res = A.dot(phi)
            # declares and initializes vector a and fills it with components of 6x6 matrix APhi
            a = np.zeros(36)
            for i in range(6):
                for j in range(6):
                    a[6 * i + j] = res[i, j]
            # calculates r1, r2 and the state vector dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
            r1 = np.sqrt((y[36] + mu) ** 2 + y[37] ** 2 + y[38] ** 2)
            r2 = np.sqrt((y[36] - (1 - mu)) ** 2 + y[37] ** 2 + y[38] ** 2)
            ydot = np.array([y[39],
                             y[40],
                             y[41],
                             y[36] + 2 * y[40] - ((1 - mu) * (y[36] + mu)) / (r1 ** 3) - (mu * (y[36] - (1 - mu))) / (
                                     r2 ** 3),
                             y[37] - 2 * y[39] - ((1 - mu) * y[37]) / (r1 ** 3) - (mu * y[37]) / (r2 ** 3),
                             - ((1 - mu) * y[38]) / (r1 ** 3) - (mu * y[38]) / (r2 ** 3)])
            # transposes horizontal vectors a and c to vertical vectors
            a = np.vstack(a)
            ydot = np.vstack(ydot)
            # declares and initializes vector APhi and fills it with components of matrix APhi
            APhi = np.zeros(42)
            APhi = np.vstack((a, ydot))

            return APhi

        elif len(y) == 6:
            # calculates r1, r2 and the state vector dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
            r1 = np.sqrt((y[0] + mu) ** 2 + y[1] ** 2 + y[2] ** 2)
            r2 = np.sqrt((y[0] - (1 - mu)) ** 2 + y[1] ** 2 + y[2] ** 2)
            ydot = np.array([y[3],
                             y[4],
                             y[5],
                             y[0] + 2 * y[4] - ((1 - mu) * (y[0] + mu)) / (r1 ** 3) - (mu * (y[0] - (1 - mu))) / (
                                     r2 ** 3),
                             y[1] - 2 * y[3] - ((1 - mu) * y[1]) / (r1 ** 3) - (mu * y[1]) / (r2 ** 3),
                             - ((1 - mu) * y[2]) / (r1 ** 3) - (mu * y[2]) / (r2 ** 3)])

            return ydot

        else:
            print("Dimension of input parameter y is not supported.")

    # --------------------------------------------------------------------------
    # STATE TRANSITION MATRIX
    #
    # DESCRIPTION:      The State Transition Matrix (STM) is a linear map from the initial state
    #                   at the initial time t = 0 to a state at some later time t and therefore
    #                   offers a tool to approximate the impact of variations in the initial
    #                   state on the evolution of the trajectory
    #                   The STM is defined by the following equations:
    #
    #                             x(t) = Phi(t,0) * x(0)
    #                    dPhi(t,0)/dt = A(t) * Phi(t,0)
    #                        Phi(t,0) = I6
    #
    # OUTPUT:           Output is the 6x6 State Transition Matrix
    #
    # SYNTAX:           stm(x, tf, mu)
    # x             =   Initial State Guess
    # tf            =   Final time
    # mu            =   Mass ratio of Primaries
    # --------------------------------------------------------------------------

    @staticmethod
    def stm(x, tf, mu):

        # declares and initializes initial state vector including STM data
        y0 = np.zeros(42)
        # declares and initializes the initial STM as a 6x6 unit matrix
        I = np.eye(6)
        # fills vector y with STM components and the initial state x
        for i in range(6):
            y0[36 + i] = x[i]
            for j in range(6):
                y0[6 * i + j] = I[i, j]
        # numerically integrates the system of ODEs
        Y = NumericalMethods.rk4System(y0, tf, mu)
        # gets dimension of the matrix Y and stores number of rows in r
        d = np.shape(Y)
        r = d[0]
        # declares and initializes vector including the STM data of the last time step
        y = np.zeros((r, 36))
        y = Y[(r - 1), 0:36]
        # declares and initiazlizes matrix phi and fills it with the elements of q
        phi = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                phi[i, j] = y[6 * i + j]

        return phi

    # --------------------------------------------------------------------------
    # HALF PERIOD
    #
    # DESCRIPTION:      The equations of motion are integrated until y changes sign. Then the step
    #                   size is reduced and the integration goes forward again starting at the last
    #                   point before the change of the sign. This is repeated until abs(y) < epsilon.
    #
    # OUTPUT:           Value of half period time
    #
    # SYNTAX:           halfPeriod(x0, mu, epsilon)
    # x0            =   Initial State
    # mu            =   Mass ratio of Primaries
    # epsilon       =   Error Tolerance
    # --------------------------------------------------------------------------

    @staticmethod
    def halfPeriod(x0, mu, epsilon):
        # declares and initializes actual state x and state at half period xHalfPeriod
        if x0[4] >= 0:
            x = np.ones(6)
            xHalfPeriod = np.ones(6)
        else:
            x = - np.ones(6)
            xHalfPeriod = - np.ones(6)
        # step size
        stepSize = 0.5
        timeStep = 0
        # number of steps for numerical integration
        n = 2000
        # counter for step size reductions
        counter = 0

        while abs(xHalfPeriod[1]) > epsilon:

            if x0[4] >= 0:
                while x[1] > 0:
                    if timeStep > 2:
                        print("        Half Period could not be calculated.")
                        raise ValueError
                    xHalfPeriod = x
                    timeStep = timeStep + stepSize
                    Y = NumericalMethods.rk4System(x0, timeStep, mu)
                    x = Y[n, :]
                    # print("y = %10.8e at t = %10.8e" % (xHalfPeriod[1], timeStep))
            else:
                while x[1] < 0:
                    if timeStep > 2:
                        print("        Half Period could not be calculated.")
                        raise ValueError
                    xHalfPeriod = x
                    timeStep = timeStep + stepSize
                    Y = NumericalMethods.rk4System(x0, timeStep, mu)
                    x = Y[n, :]
                    # print("y = %10.8e at t = %10.8e" % (xHalfPeriod[1], timeStep))

            # last point before change of sign is defined as starting point for next iteration
            timeStep = timeStep - stepSize
            # last state before change of sign is defined as starting point for next iteration
            x = xHalfPeriod
            # step size is reduced
            stepSize = stepSize * 1.0e-1
            # counter is incremented
            counter = counter + 1

        tHalfPeriod = timeStep

        return tHalfPeriod

    @staticmethod
    def nullspace(A, atol=1e-13, rtol=0):
        A = np.atleast_2d(A)
        u, s, vh = svd(A)
        tol = max(atol, rtol * s[0])
        nnz = (s >= tol).sum()
        ns = vh[nnz:].conj().T
        return ns


class Plot:

    @staticmethod
    def setAxesEqual(ax):
        '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

        Input
          ax: a matplotlib axis, e.g., as output from plt.gca().
        '''

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5 * max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    @staticmethod
    def plot(data, mu, haloFamily, background):

        # Prepares figure for plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # norm = (jacobi - 2.7) / (3.2 - 2.7)  #!!! Hard coded
        # color = plt.cm.jet(norm)
        color = 'blue'

        #ax.scatter(-mu, 0, 0, color='black', s = 15)
        ax.scatter(1-mu, 0, 0, color='white', s = 8, label='Moon')

        l = 1-mu
        p_L1= np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)+mu-l, mu**2*l**2+2*(l**2+mu**2), mu**3-l**3])
        L1roots = np.roots(p_L1)
        for i in range(5):
            if L1roots[i] > - mu and L1roots[i] < l:
                L1 = np.real(L1roots[i])
        ax.scatter(L1, 0, 0, color='red', s=3, label='L1/L2')


        p_L2 = np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)-(mu+l), mu**2*l**2+2*(l**2-mu**2), -(mu**3+l**3)])
        L2roots=np.roots(p_L2)
        for i in range(5):
            if L2roots[i] > - mu and L2roots[i] > l:
                L2 = np.real(L2roots[i])
        ax.scatter(L2, 0, 0, color='red', s=3)



        for i in range(len(data)):
            halo = NumericalMethods.rk4System(data[i][0:6], data[i][6], mu)
            x = halo[:, 0]
            y = halo[:, 1]
            z = halo[:, 2]
            if haloFamily == "northern":
                ax.plot(x, y, -z, color=color, linewidth=1)
            elif haloFamily == "southern":
                ax.plot(x, y, z, color=color, linewidth=0.7)
            elif haloFamily == "both":
                ax.plot(x, y, z, color=color, linewidth=1)
                ax.plot(x, y, -z, color=color, linewidth=1)
            else:
                print("Input of haloFamily is not supported.")
                exit()
            if background == "on":
                ax.patch.set_facecolor('black')
                ax.set_axis_off()
            elif background == "off":
                ax.set_xlabel("x Axis")
                ax.set_ylabel("y Axis")
                ax.set_zlabel("z Axis")
            else:
                print("Input of background is not supported.")

        #ax.legend(loc="center right", markerscale=1., scatterpoints=1, fontsize=10)
        ax.view_init(elev = 0, azim = -90)
        fig.savefig("fig.pdf", format='pdf', dpi=500, bbox_inches = 'tight')

        Plot.setAxesEqual(ax)
        plt.show()

    # @staticmethod
    # def setAccuracy(accuracy):
    #    Orbit.setAccuracy(accuracy)
    #    OrbitContinuation.setAccuracy(accuracy)
