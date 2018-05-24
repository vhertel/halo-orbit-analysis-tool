"""
File    : Utility.py
Author  : Victor Hertel
Date    : 24.04.2018

Utility classes
"""

# Imports
import numpy as np
from scipy.integrate import odeint
import os
from numpy.linalg import svd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D





class NumericalMethods:

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
                if counter > 10:
                    print("        Differential Corrections Method did not converge.")
                    raise ValueError
                else:
                    counter = counter + 1
                # calculates the state transition matrix
                phi = Utility.stm(x, tau, mu)
                # integrates initial state from t0 to tau_n
                t = np.linspace(0, tau, num=200000)
                xRef = odeint(Utility.sysEquations, x, t, args=(mu,))
                # calculates the derivation of the state at T/2
                xdot = Utility.sysEquations(xRef[-1, :], 0, mu)
                # declares and initializes free variable vector with x, z, ydot and tau
                freeVariables = np.array([x[0], x[2], x[4], tau])
                # declares and initializes the constraint vector with y, xdot and zdot
                constraints = np.array([xRef[-1, 1], xRef[-1, 3], xRef[-1, 5]])
                # calculates corrections to the initial state to meet a defined margin of error
                DF = np.array([[phi[1, 0], phi[1, 2], phi[1, 4], xdot[1]],
                               [phi[3, 0], phi[3, 2], phi[3, 4], xdot[3]],
                               [phi[5, 0], phi[5, 2], phi[5, 4], xdot[5]]])
                xIter = freeVariables - ((DF.T).dot(np.linalg.inv(DF.dot(DF.T)))).dot(constraints)
                # sets the updated initial condition vector
                x = np.array([xIter[0], 0, xIter[1], 0, xIter[2], 0])
                # sets T/2 of updated initial conditions
                tau = xIter[3]
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
                raise ValueError

            # declares and initializes constraint vector
            constraints = np.ones((2, 1))
            # declares and initializes counter
            counter = 0

            # case of x-amplitude being the fixed variable
            if fixedValue == "x":

                # constraints are corrected until they meet a defined margin of error
                while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon:
                    if counter > 10:
                        print("        Differential Corrections Method did not converge.")
                        raise ValueError
                    else:
                        counter = counter + 1
                    # calculates the state transition matrix
                    phi = Utility.stm(x, tau, mu)
                    # integrates initial state from 0 to tau
                    t = np.linspace(0, tau, num=200000)
                    xRef = odeint(Utility.sysEquations, x, t, args=(mu,))
                    # calculates the derivation of the state at T/2
                    xdot = Utility.sysEquations(xRef[-1, :], 0, mu)
                    # declares and initializes free variable vector with z, ydot and time
                    freeVariables = np.array([x[2], x[4], tau])
                    # declares and initializes the constraint vector with y, xdot and zdot
                    constraints = np.array([xRef[-1, 1], xRef[-1, 3], xRef[-1, 5]])
                    # calculates corrections to the initial state to meet a defined margin of error
                    D = np.array([[phi[1, 2], phi[1, 4], xdot[1]],
                                  [phi[3, 2], phi[3, 4], xdot[3]],
                                  [phi[5, 2], phi[5, 4], xdot[5]]])
                    DF = np.linalg.inv(D)
                    xIter = freeVariables - DF.dot(constraints)
                    # sets the updated initial condition vector
                    x = np.array([x[0], 0, xIter[0], 0, xIter[1], 0])
                    # calculates T/2 with updated initial conditions
                    tau = xIter[2]
                    outX = x
                    outTime = tau

            # case of z-amplitude being the fixed variable
            elif fixedValue == "z":

                # constraints are corrected until they meet a defined margin of error
                while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon:
                    if counter > 10:
                        print("        Differential Corrections Method did not converge.")
                        raise ValueError
                    else:
                        counter = counter + 1
                    # calculates the state transition matrix
                    phi = Utility.stm(x, tau, mu)
                    # integrates initial state from 0 to tau
                    t = np.linspace(0, tau, num=200000)
                    xRef = odeint(Utility.sysEquations, x, t, args=(mu,))
                    # calculates the derivation of the state at T/2
                    xdot = Utility.sysEquations(xRef[-1, :], 0, mu)
                    # declares and initializes free variable vector with x, ydot and time
                    freeVariables = np.array([x[0], x[4], tau])
                    # declares and initializes the constraint vector with y, xdot and zdot
                    constraints = np.array([xRef[-1, 1], xRef[-1, 3], xRef[-1, 5]])
                    # calculates corrections to the initial state to meet a defined margin of error
                    D = np.array([[phi[1, 0], phi[1, 4], xdot[1]],
                                  [phi[3, 0], phi[3, 4], xdot[3]],
                                  [phi[5, 0], phi[5, 4], xdot[5]]])
                    DF = np.linalg.inv(D)
                    xIter = freeVariables - DF.dot(constraints)
                    # sets the updated initial condition vector
                    x = np.array([xIter[0], 0, x[2], 0, xIter[1], 0])
                    # calculates T/2 with updated initial conditions
                    tau = xIter[2]
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
    def sysEquations(y, t, mu):

        if len(y) == 42:
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

            APhi = np.concatenate((a,ydot))

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
        t = np.linspace(0, tf, num=10000)
        Y = odeint(Utility.sysEquations, y0, t, args=(mu,))
        # declares and initializes vector including the STM data of the last time step
        y = Y[-1, 0:36]
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
                    t = np.linspace(0, timeStep, num=10000)
                    Y = odeint(Utility.sysEquations, x0, t, args=(mu,))
                    x = Y[-1, :]
                    #print("y = %10.8e at t = %10.8e" % (xHalfPeriod[1], timeStep))
            else:
                while x[1] < 0:
                    if timeStep > 2:
                        print("        Half Period could not be calculated.")
                        raise ValueError
                    xHalfPeriod = x
                    timeStep = timeStep + stepSize
                    t = np.linspace(0, timeStep, num=10000)
                    Y = odeint(Utility.sysEquations, x0, t, args=(mu,))
                    x = Y[-1, :]
                    #print("y = %10.8e at t = %10.8e" % (xHalfPeriod[1], timeStep))

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
    def lagrangianPosition(mu):
        # L1
        l = 1-mu
        p_L1= np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)+mu-l, mu**2*l**2+2*(l**2+mu**2), mu**3-l**3])
        L1roots = np.roots(p_L1)
        for i in range(5):
            if L1roots[i] > - mu and L1roots[i] < l:
                L1 = np.real(L1roots[i])
        # L2
        p_L2 = np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)-(mu+l), mu**2*l**2+2*(l**2-mu**2), -(mu**3+l**3)])
        L2roots=np.roots(p_L2)
        for i in range(5):
            if L2roots[i] > - mu and L2roots[i] > l:
                L2 = np.real(L2roots[i])

        return np.array([L1, L2])

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
        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()
        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)
        plot_radius = 0.5 * max([x_range, y_range, z_range])
        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


    @staticmethod
    def plot(data, mu, dict, haloFamily, background, save):

        # prepares figure for plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # plots lagrangian points
        lagrangian = Utility.lagrangianPosition(mu)
        ax.scatter(lagrangian[0], 0, 0, color='blue', s=0.3, label='L1/L2')
        ax.scatter(lagrangian[1], 0, 0, color='blue', s=0.3, label='L1/L2')

        # plots second primary
        #ax.scatter(-mu, 0, 0, color='blue', s = 3)
        ax.scatter(1-mu, 0, 0, color='grey', s = 1, label='Moon')

        jacobiMax = max(data[:, 0])
        jacobiMin = min(data[:, 0])

        # plots orbits
        for i in range(len(data)):

            norm = (data[i,0] - jacobiMin) / (jacobiMax - jacobiMin)
            color = plt.cm.hsv_r(norm)
            t = np.linspace(0, data[i][1], num=10000)
            halo = odeint(Utility.sysEquations, data[i][2:8], t, args=(mu,))
            x = halo[:, 0]
            y = halo[:, 1]
            z = halo[:, 2]
            if haloFamily == "northern":
                ax.plot(x, y, -z, color=color, linewidth=0.25)
            elif haloFamily == "southern":
                ax.plot(x, y, z, color=color, linewidth=0.25)
            elif haloFamily == "both":
                ax.plot(x, y, z, color=color, linewidth=0.25)
                ax.plot(x, y, -z, color=color, linewidth=0.25)
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

        Plot.setAxesEqual(ax)
        plt.show()

        if save:
            if not os.path.exists(dict + "plots"):
                os.makedirs(dict + "plots")
            print("\nSaving plots..")
            numOfFigures = 100
            theta = np.linspace(0, 2*np.pi, numOfFigures)
            for i in range(0, numOfFigures, 1):
                 ax.view_init(0, i*(360/numOfFigures) + 270)
                 fig.savefig(dict + "plots/fig%d.pdf" % (i), format='pdf', dpi=400, bbox_inches = 'tight')
            print("DONE")


    @staticmethod
    def plotJacobi(data, dict, save):
        plt.title("Earth - Moon L1/L2")
        plt.xlabel("x-Axis")
        plt.ylabel("Jacobi Constant")
        for element in data:
            if element[2] > 1:
                color = 'blue'
            else:
                color = 'red'
            plt.scatter(element[2], element[0], s=0.5, color=color)
        if save:
            if not os.path.exists(dict + "plots"):
                os.makedirs(dict + "plots")
            plt.savefig(dict + "plots/jacobi.pdf", format='pdf', dpi=500, bbox_inches = 'tight')
        plt.show()


    @staticmethod
    def plotPeriod(data, dict, save):
        plt.title("Earth - Moon L1/L2")
        plt.xlabel("x-Axis")
        plt.ylabel("Period")
        for element in data:
            if element[2] > 1:
                color = 'blue'
            else:
                color = 'red'
            plt.scatter(element[2], element[1], s=0.5, color=color)
        if save:
            if not os.path.exists(dict + "plots"):
                os.makedirs(dict + "plots")
            plt.savefig(dict + "plots/period.pdf", format='pdf', dpi=500, bbox_inches = 'tight')
        plt.show()


    @staticmethod
    def plotStability(data, mu, dict, save):
        plt.title("Earth - Moon L1/L2")
        plt.xlabel("x-Axis")
        plt.ylabel("Stability Index")
        for element in data:
            monodromy = Utility.stm(element[2:8], element[1], mu)
            eigenvalues, eigenvectors = np.linalg.eig(monodromy)
            maximum = max(abs(eigenvalues))
            stability = 1 / 2 * (maximum + 1 / maximum)
            if element[2] > 1:
                color = 'blue'
            else:
                color = 'red'
            plt.scatter(element[2], stability, s=0.5, color=color)
        if save:
            if not os.path.exists(dict + "plots"):
                os.makedirs(dict + "plots")
            plt.savefig(dict + "plots/stability.pdf", format='pdf', dpi=500, bbox_inches = 'tight')
        plt.show()


    @staticmethod
    def plotFromTable(folder, jacobi=False, period=False, stability=False, orbits=False, orbitNumber=None, haloFamily="both", background="off", save=False):

        table = open("Output/" + folder + "/data.txt", "r")
        for lines in table:

            if "MASS RATIO" in lines:
                line = lines.split()
                mu = float(line[3])
            if "ORBIT NUMBER" in lines:
                line = lines.split()
                number = int(line[3])
            if "DATA_START" in lines:
                for i in range(4):
                    line = table.readline()
                words = line.split()
                for i in range(3):
                    words.insert(2*i + 3, 0)
                data = words
                while line != "\n":
                    words = line.split()
                    for i in range(3):
                        words.insert(2*i + 3, 0)
                    data = np.vstack([data, words])
                    line = table.readline()

                data = data.astype(np.float)
        table.close()

        dict = "Output/" + folder + "/"

        if jacobi:
            Plot.plotJacobi(data, dict, save)
        if period:
            Plot.plotPeriod(data, dict, save)
        if stability:
            Plot.plotStability(data, mu, dict, save)
        if orbits:
            if orbitNumber is None:
                Plot.plot(data, mu, dict=dict, haloFamily=haloFamily, background=background, save=save)
            else:
                step = number/(orbitNumber-1)
                reducedData = data[0, :]
                for i in range(orbitNumber-2):
                    reducedData = np.vstack([reducedData, data[round(step*(i+1)), :]])
                reducedData = np.vstack([reducedData, data[-1, :]])
                Plot.plot(reducedData, mu, dict=dict, haloFamily=haloFamily, background=background, save=save)
