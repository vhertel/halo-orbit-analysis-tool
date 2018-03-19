
"""
File    : numMethods.py
Author  : Victor Hertel
Date    : 18.03.2018

Numerical methods for solving differential equations.
"""

# Imports
import numpy as np
from library import haloCalculation




#--------------------------------------------------------------------------
# FOURTH ODER RUNGE-KUTTA METHOD

# DESCRIPTION:      The function rk4System() approximates the solutions of a system of m
#                   differential equations that are written in the form
#
#                   dy1/dt = f1(t,y1,y2,...,ym)
#                   dy2/dt = f2(t,y1,y2,...,ym)
#                   ...
#                   dym/dt = fm(t,y1,y2,...,ym)
#
#                   with t in the interval [t0; tf] and the m-dimensional initial conditions
#                   vector y0. The accuracy depends on the step size n.
#
# OUTPUT:           A matrix consisting of a row for each time and 42 columns including the elements
#                   of the State Transition Matrix and the state vector dx/dt
#
# SYNTAX:           rk4System(y, t0, tf, mu)
# y             =   Position values of state vector
# t0            =   Start time
# tf            =   Final time
# mu            =   Mass ratio of Primaries
#--------------------------------------------------------------------------

def rk4System(y, t0, tf, mu):

    # defines step size depending on time interval [t0; tf] and n
    n = 2000
    h = (tf-t0)/n
    # declares and initializes vector t with dimensions (n+1)
    t = np.zeros(n+1)
    t[0] = t0
    # declares and initializes matrix w with a row for each time step and 42 columns
    R = np.zeros((len(y),n+1))
    for i in range(len(y)):
        R[i,0] = y[i]

    # 4th order Runge-Kutta method algorithm
    if len(y) == 42:

        for i in range(n):
            k1 = haloCalculation.sysEquations(R[:,i], mu, t[i])
            k2 = haloCalculation.sysEquations(R[:,i] + h/2*k1[:,0], mu, t[i] + h/2)
            k3 = haloCalculation.sysEquations(R[:,i] + h/2*k2[:,0], mu, t[i] + h/2)
            k4 = haloCalculation.sysEquations(R[:,i] + h*k3[:,0], mu, t[i] + h)
            R[:,i+1] = R[:,i] + h/6*(k1[:,0] + 2*k2[:,0] + 2*k3[:,0] + k4[:,0])
            t[i+1] = t[i] + h
        R = R.T

    elif len(y) == 6:

        for i in range(n):
            k1 = haloCalculation.sysEquations(R[:,i], mu, t[i])
            k2 = haloCalculation.sysEquations(R[:,i] + h/2*k1, mu, t[i] + h/2)
            k3 = haloCalculation.sysEquations(R[:,i] + h/2*k2, mu, t[i] + h/2)
            k4 = haloCalculation.sysEquations(R[:,i] + h*k3, mu, t[i] + h)
            R[:,i+1] = R[:,i] + h/6*(k1 + 2*k2 + 2*k3 + k4)
            t[i+1] = t[i] + h
        R = R.T

    else:
        print("Dimension of input parameter y is not supported.")

    return R




#--------------------------------------------------------------------------
# DIFFERENTIAL CORRECTIONS METHOD
#
#
#--------------------------------------------------------------------------

def diffCorrections(x0, t0, mu, epsilon, fixedValue):

    print("        Differential Corrections Method for adaption of the initial state...")
    # sets the iteratively adjusted initial condition
    x_n = x0
    # calculates T/2 by integrating until y changes sign
    try:
        tau_n = haloCalculation.halfPeriod(x_n, t0, mu, 1.0e-11)
    except ValueError:
        return

    # declares and initializes constraint vector
    constraints = np.ones((2,1))
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
            phi = haloCalculation.stm(x_n, t0, tau_n, mu)
            # integrates initial state from t0 to tau_n
            xRef = rk4System(x_n, t0, tau_n, mu)
            # calculates the derivation of the state at T/2
            xdot = haloCalculation.sysEquations(xRef[2000, :], mu)
            # declares and initializes free variable vector with z, ydot and time
            freeVariables = np.array([[x_n[2]],
                                      [x_n[4]],
                                      [tau_n]])
            # declares and initializes the constraint vector with y, xdot and zdot
            constraints = np.array([[xRef[2000,1]],
                                    [xRef[2000,3]],
                                    [xRef[2000,5]]])
            # calculates corrections to the initial state to meet a defined margin of error
            D = np.array([[phi[1,2], phi[1,4], xdot[1]],
                          [phi[3,2], phi[3,4], xdot[3]],
                          [phi[5,2], phi[5,4], xdot[5]]])
            DF = np.linalg.inv(D)
            xIter = freeVariables - DF.dot(constraints)
            # sets the updated initial condition vector
            x_n = np.array([x0[0], 0, xIter[0,0], 0, xIter[1,0], 0])
            # calculates T/2 with updated initial conditions
            tau_n = xIter[2,0]

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
            phi = haloCalculation.stm(x_n, t0, tau_n, mu)
            # integrates initial state from t0 to tau_n
            xRef = rk4System(x_n, t0, tau_n, mu)
            # calculates the derivation of the state at T/2
            xdot = haloCalculation.sysEquations(xRef[2000, :], mu)
            # declares and initializes free variable vector with x, ydot and time
            freeVariables = np.array([[x_n[0]],
                                      [x_n[4]],
                                      [tau_n]])
            # declares and initializes the constraint vector with y, xdot and zdot
            constraints = np.array([[xRef[2000,1]],
                                    [xRef[2000,3]],
                                    [xRef[2000,5]]])
            # calculates corrections to the initial state to meet a defined margin of error
            D = np.array([[phi[1,0], phi[1,4], xdot[1]],
                          [phi[3,0], phi[3,4], xdot[3]],
                          [phi[5,0], phi[5,4], xdot[5]]])
            DF = np.linalg.inv(D)
            xIter = freeVariables - DF.dot(constraints)
            # sets the updated initial condition vector
            x_n = np.array([xIter[0,0], 0, x0[2], 0, xIter[1,0], 0])
            # calculates T/2 with updated initial conditions
            tau_n = xIter[2,0]

    # stores the initial condition and T/2 in output vector
    outData = np.array([x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5], tau_n])
    # prints status updates
    print("        Initial state has been adapted for %d times:       -> x0 = [%0.8f, %d, %8.8f, %d, %8.8f, %d]" % (counter, outData[0], outData[1], outData[2], outData[3], outData[4], outData[5]))
    print("        Constraints at T/2 = %6.5f:                     -> [y, dx/dt, dz/dt] = [%6.5e, %6.5e, %6.5e]" % (outData[6], constraints[0], constraints[1], constraints[2]))

    return outData
