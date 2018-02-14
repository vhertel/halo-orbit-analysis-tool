
"""
File    : calculation.py
Author  : Victor Hertel
Date    : 08.02.2018

Calculates individual orbits or families of halo orbits.
"""

# Imports
import numpy as np
import matplotlib as mpl
# necessary to use matplotlib for Mac
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from library import haloCalculation, numMethods, utilities
import math




#--------------------------------------------------------------------------
# SINGLE HALO ORBIT
#
# DESCRIPTION:      Searches for one halo orbit with fixing either the
#                   x- or z- amplitude of the initial guess
#
# OUTPUT:           No output is generated. The orbit is plotted in the
#                   figure that has been created before the function call
#
# SYNTAX:           singleHalo(x0, t0, mu, epsilon, fixedValue, haloFamily, ax)
#
# x0            =   Initial State Guess
# t0            =   Time Start
# mu            =   Mass ratio of Primaries
# epsilon       =   Error Tolerance of Constraints at T/2
# fixedValue    =   The value of the initial guess that should be fixed:
#                   Input possibilities: {"x", "z"}
# haloFamily    =   Desired family of Halo Orbits:
#                   Input possibilities: {"northern", "southern", "both"}
# ax            =   Allows access to the figure
#--------------------------------------------------------------------------
def singleHalo(x0, t0, mu, epsilon, fixedValue, haloFamily, ax):
    # prints status update
    print("STATUS: Single Halo Orbit Computation...")
    # stores initial state and half period of halo orbit in 1x7 vector outData
    outData = numMethods.diffCorrections(x0, t0, mu, epsilon, fixedValue)
    # stores initial state
    outX = outData[0:6]
    # stores half period
    outTime = outData[6]

    # calculates monodromy matrix
    #monodromy = haloCalculation.stm(outX, t0, 2*outTime, mu)
    # calculates eigenvalues of monodromy matrix
    #eigenvalues, eigenvectors = np.linalg.eig(monodromy)
    #print(eigenvalues)
    #stability = haloCalculation.stability(eigenvalues)
    #print(stability)

    # calculates Jacobi constant
    J = haloCalculation.jacobiConst(outX, mu)

    # plots halo orbit
    halo = numMethods.rk4System(outX, t0, 2*outTime, mu)
    x = halo[:, 0]
    y = halo[:, 1]
    z = halo[:, 2]
    # defines colormap in range of Jacobi constant
    norm = (J - 3.02) / (3.2 - 3.02)                                #!!! Hard coded
    color = plt.cm.jet(norm)
    if haloFamily == "northern" :
        ax.plot(x, y, -z, color=color, linewidth=1)
    elif haloFamily == "southern" :
        ax.plot(x, y, z, color=color, linewidth=1)
    elif haloFamily == "both":
        ax.plot(x, y, z, color=color, linewidth=1)
        ax.plot(x, y, -z, color=color, linewidth=1)
    else:
        print("Input of haloFamily is not supported.")
        exit()
    print("DONE")




#--------------------------------------------------------------------------
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
# OUTPUT:           No output is generated. The orbits are plotted in the
#                   figure that has been created before the function call
#
# SYNTAX:           natParaConti(x0, t0, mu, epsilon, orbitNumber, familyStep, lagrangian, haloFamily, ax)
#
# x0            =   Initial State Guess
# t0            =   Time Start
# mu            =   Mass ratio of Primaries
# epsilon       =   Error Tolerance of Constraints at T/2
# orbitNumber   =   Number of Orbits to search for
# familyStep    =   Distance between Orbits
# lagrangian    =   Lagrangian Point
#                   Input possibilities: {"L1", "L2"}
# haloFamily    =   Desired family of Halo Orbits:
#                   Input possibilities: {"northern", "southern", "both"}
# ax            =   Allows access to the figure
#--------------------------------------------------------------------------
def natParaConti(x0, t0, mu, epsilon, orbitNumber, familyStep, lagrangian, haloFamily, ax):

    print("STATUS: Natural Parameter Continuation method for generating a family of %2d orbits...\n" % (orbitNumber))
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
            stepSize = np.sqrt(familyStep**2 - dz**2)
            if math.isnan(stepSize):
                stepSize = familyStep/2
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
            outData = numMethods.diffCorrections(x_n, t0, mu, epsilon, "x")
            outX = outData[0:6]
            outTime = outData[6]

        else:
            # z-value changed more than x-value and needs to be fixed
            print("        Orbit Number: %2d    (fixed z-value)\n" % (i + 1))
            dx = abs(outX[0] - lastX[0])
            stepSize = np.sqrt(familyStep**2 - dx**2)
            if math.isnan(stepSize):
                stepSize = familyStep/2
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
            outData = numMethods.diffCorrections(x_n, t0, mu, epsilon, "z")
            outX = outData[0:6]
            outTime = outData[6]

        if outData[7] > 35:
            break

        J = haloCalculation.jacobiConst(outX, mu)
        halo = numMethods.rk4System(outX, t0, 2*outTime, mu)
        x = halo[:, 0]
        y = halo[:, 1]
        z = halo[:, 2]
        # defines colormap in range of Jacobi constant
        norm = (J - 2.9) / (3.2 - 2.9)                                #!!! Hard coded
        color = plt.cm.jet(norm)
        if haloFamily == "northern" :
            ax.plot(x, y, -z, color=color, linewidth=1)
        elif haloFamily == "southern" :
            ax.plot(x, y, z, color=color, linewidth=1)
        elif haloFamily == "both":
            ax.plot(x, y, z, color=color, linewidth=1)
            ax.plot(x, y, -z, color=color, linewidth=1)
        else:
            print("Input of halo family is not supported.")
            exit()
    print("DONE")




#--------------------------------------------------------------------------
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
# direction     =   Direction of continuation
#                   Input possibilities: {"planar", "vertical"}
# haloFamily    =   Desired family of Halo Orbits:
#                   Input possibilities: {"northern", "southern", "both"}
# ax            =   Allows access to the figure
#--------------------------------------------------------------------------
def pseudoArcLenConti(x0, t0, mu, epsilon, orbitNumber, familyStep, direction, haloFamily, ax):

    print("STATUS: Pseudo Arc-Length Continuation method for generating a family of %2d orbits...\n" % (orbitNumber))
    # sets the iteratively adjusted initial condition
    x_n = x0
    # calculates T/2 by integrating until y changes sign
    tau_n = haloCalculation.halfPeriod(x_n, t0, mu, 1.0e-11)
    # declares and initializes constraint vector
    constraints = np.ones((3,1))
    # calculation of halo orbits
    for i in range(orbitNumber):
        # prints status update
        print("        Orbit Number: %2d" % (i + 1))
        print("        Differential corrections...")
        # declares and initializes counter for counting updates of initial state
        counter = 1
        # constraints are corrected until they meet a defined margin of error
        while abs(constraints[0]) > epsilon or abs(constraints[1]) > epsilon or abs(constraints[2]) > epsilon:
            # calculates the state transition matrix
            phi = haloCalculation.stm(x_n, t0, tau_n, mu)
            # integrates initial state from t0 to tau_n
            xRef = numMethods.rk4System(x_n, t0, tau_n, mu)
            # calculates the derivation of the state at T/2
            xdot = haloCalculation.sysEquations(xRef[2000, :], mu)
            # declares and initializes free variable vector with x, z, ydot and tau
            freeVariables = np.array([[x_n[0]],
                                      [x_n[2]],
                                      [x_n[4]],
                                      [tau_n]])
            # declares and initializes the constraint vector with y, xdot and zdot
            constraints = np.array([[xRef[2000,1]],
                                    [xRef[2000,3]],
                                    [xRef[2000,5]]])
            # calculates corrections to the initial state to meet a defined margin of error
            DF = np.array([[phi[1,0], phi[1,2], phi[1,4], xdot[1]],
                           [phi[3,0], phi[3,2], phi[3,4], xdot[3]],
                           [phi[5,0], phi[5,2], phi[5,4], xdot[5]]])
            xIter = freeVariables - ((DF.T).dot(np.linalg.inv(DF.dot(DF.T)))).dot(constraints)
            # sets the updated initial condition vector
            x_n = np.array([xIter[0,0], 0, xIter[1,0], 0, xIter[2,0], 0])
            # sets T/2 of updated initial conditions
            tau_n = xIter[3,0]
            counter = counter + 1
        # prints status update
        print("        Initial state has been adapted for %2d times:" % (counter))
        print("        -> x0 = [%10.8e, %2d, %10.8e, %2d, %10.8e, %2d]" % (x_n[0], x_n[1], x_n[2], x_n[3], x_n[4], x_n[5]))
        print("        Constraints at T/2 = %10.8e:" % (tau_n))
        print("        -> [y, dx/dt, dz/dt] = [%10.8e, %10.8e, %10.8e]\n" % (constraints[0], constraints[1], constraints[2]))
        outX = x_n
        outTime = tau_n

        # calculates Jacobi constant
        J = haloCalculation.jacobiConst(outX, mu)

        # plots halo orbit
        halo = numMethods.rk4System(outX, t0, 2*outTime, mu)
        x = halo[:, 0]
        y = halo[:, 1]
        z = halo[:, 2]
        # defines colormap in range of Jacobi constant
        norm = (J - 3.02) / (3.2 - 3.02)                                #!!! Hard coded
        color = plt.cm.jet(norm)
        if haloFamily == "northern" :
            ax.plot(x, y, -z, color=color, linewidth=1)
        elif haloFamily == "southern" :
            ax.plot(x, y, z, color=color, linewidth=1)
        elif haloFamily == "both":
            ax.plot(x, y, z, color=color, linewidth=1)
            ax.plot(x, y, -z, color=color, linewidth=1)
        else:
            print("Input of haloFamily is not supported.")
            exit()

        # calculates the null space vector of Jacobian matrix DF
        nullSpace = utilities.nullspace(DF)
        # declares and initializes the augmented free variable vector
        contiFreeVariables = np.array([[x_n[0]],
                                  [x_n[2]],
                                  [x_n[4]],
                                  [tau_n]])
        # declares and initializes the augmented constraints vector
        if direction == "planar":
            contiConstraints = np.array([[xRef[2000,1]],
                                         [xRef[2000,3]],
                                         [xRef[2000,5]],
                                         [((contiFreeVariables - freeVariables).T).dot(nullSpace) - familyStep]])
        elif direction == "vertical":
            contiConstraints = np.array([[xRef[2000,1]],
                                         [xRef[2000,3]],
                                         [xRef[2000,5]],
                                         [((contiFreeVariables - freeVariables).T).dot(nullSpace) + familyStep]])
        else:
            print("Input of continuation direction is not supported.")
            exit()        # calculates corrections to the initial state to meet a defined margin of error
        GF = np.array([[phi[1,0], phi[1,2], phi[1,4], xdot[1]],
                       [phi[3,0], phi[3,2], phi[3,4], xdot[3]],
                       [phi[5,0], phi[5,2], phi[5,4], xdot[5]],
                       [nullSpace[0], nullSpace[1], nullSpace[2], nullSpace[3]]])
        xIter = contiFreeVariables - (np.linalg.inv(GF)).dot(contiConstraints)
        # sets the updated initial condition vector
        x_n = np.array([xIter[0,0], 0, xIter[1,0], 0, xIter[2,0], 0])
        # sets T/2 of updated initial conditions
        tau_n = xIter[3,0]
        # resets constraints for next iteration
        constraints = np.ones((3,1))
    print("DONE")




#--------------------------------------------------------------------------
# PLOT OF PRIMARIES

# SYNTAX:           primaries(mu, primaries, ax)
# mu            =   Mass ratio of Primaries
# primaries     =   Defines which primaries will be plotted
#                   Input possibilities: {"first", "second", "both"}
# ax            =   Allows access to the figure
#--------------------------------------------------------------------------

def primaries(mu, primaries, ax):
    if primaries == "first":
        ax.scatter(-mu, 0, 0, color='black', s = 15)
    elif primaries == "second":
        ax.scatter(1-mu, 0, 0, color='black', s = 8)
    elif primaries == "both":
        ax.scatter(-mu, 0, 0, color='black', s = 15)
        ax.scatter(1-mu, 0, 0, color='black', s = 15)
    else:
        print("Input of primaries is not supported.")
        exit()




#--------------------------------------------------------------------------
# PLOT OF LAGRANGIAN POINTS L1 and L2

# SYNTAX:           lagrangianPoints(mu, ax)
# mu            =   Mass ratio of Primaries
# ax            =   Allows access to the figure
#--------------------------------------------------------------------------

def lagrangianPoints(mu, ax):
    l = 1-mu
    p_L1= np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)+mu-l, mu**2*l**2+2*(l**2+mu**2), mu**3-l**3])
    L1roots = np.roots(p_L1)
    for i in range(5):
        if L1roots[i] > - mu and L1roots[i] < l:
            L1 = np.real(L1roots[i])
    ax.scatter(L1, 0, 0, color='blue', s=2)


    p_L2 = np.array([1, 2*(mu-l), l**2-4*l*mu+mu**2, 2*mu*l*(l-mu)-(mu+l), mu**2*l**2+2*(l**2-mu**2), -(mu**3+l**3)])
    L2roots=np.roots(p_L2)
    for i in range(5):
        if L2roots[i] > - mu and L2roots[i] > l:
            L2 = np.real(L2roots[i])
    ax.scatter(L2, 0, 0, color='blue', s=2)
