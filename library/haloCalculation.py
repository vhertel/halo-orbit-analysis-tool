
"""
File    : haloCalculation.py
Author  : Victor Hertel
Date    : 08.02.2018

A library of specific functions necessary to find periodic solutions
to the Circular Restricted Three Body Problem (CR3BP).
"""

# Imports
import numpy as np
from library import numMethods




#--------------------------------------------------------------------------
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
#--------------------------------------------------------------------------

def AMatrix(x, mu):

    # declares and initializes 3x3 submatrix O with zeros
    O = np.zeros((3,3))
    # declares and initializes 3x3 unit submatrix I
    I = np.eye(3)
    # calculates r1 and r2
    r1 = np.sqrt((x[0]+mu)**2 + x[1]**2 + x[2]**2)
    r2 = np.sqrt((x[0]-(1-mu))**2 + x[1]**2 + x[2]**2)
    # calculates the elements of the second partial derivatives of the three-body pseudo-potential U
    Uxx = 1 - (1-mu)/(r1**3) - (mu)/(r2**3) + (3*(1-mu)*(x[0]+mu)**2)/(r1**5) + (3*mu*(x[0]-(1-mu))**2)/(r2**5)
    Uxy = (3*(1-mu)*(x[0]+mu)*x[1])/(r1**5) + (3*mu*(x[0]-(1-mu))*x[1])/(r2**5)
    Uxz = (3*(1-mu)*(x[0]+mu)*x[2])/(r1**5) + (3*mu*(x[0]-(1-mu))*x[2])/(r2**5)
    Uyx = Uxy
    Uyy = 1 - (1-mu)/(r1**3) - (mu)/(r2**3) + (3*(1-mu)*x[1]**2)/(r1**5) + (3*mu*x[1]**2)/(r2**5)
    Uyz = (3*(1-mu)*x[1]*x[2])/(r1**5) + (3*mu*x[1]*x[2])/(r2**5)
    Uzx = Uxz
    Uzy = Uyz
    Uzz = - (1-mu)/(r1**3) - (mu)/(r2**3) + (3*(1-mu)*x[2]**2)/(r1**5) + (3*mu*x[2]**2)/(r2**5)
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



#--------------------------------------------------------------------------
# CALCULATING SYSTEMS OF EQUATIONS
#
# DESCRIPTION:      Function deals with the equations of motion comprising the dynamical model of the CR3BP.
#                       case a):    input variable y has dimension of 6:
#                                   The following state vector is calculated
#                                   dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
#                       case b):    input variable y has dimension of 42:
#                                   The matrix multiplication dPhi(t,t0)/dt = A(t) * Phi(t,t0) for
#                                   calculating the State Transition Matrix STM is done at t0 with
#                                   Phi(t0,t0) = I6. Additionally the state vector
#                                   dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T is
#                                   calculated and the values stored in the vector APhi with 42 elements.
#
# OUTPUT:           case a):    input variable y has dimension of 6:
#                               Six dimensional state vector dx/dt
#                   case b):    input variable y has dimension of 42:
#                               Vector with 42 elements including the result of the matrix multiplication
#                               dPhi(t,t0)/dt = A(t) * Phi(t,t0) as the first 36 values and the six
#                               dimensional state vector dx/dt as the last 6 elements
#
# SYNTAX:           sysEquations(y, mu, t="0")
# y             =   Position values of state vector
# mu            =   Mass ratio of Primaries
# t="0"         =   Parameter is only used for numerical integration with rk4System()
#--------------------------------------------------------------------------

def sysEquations(y, mu, t="0"):

    if len(y) == 42:
        # declares and initializes vector x including the initial position
        x = np.zeros(3)
        # stores initial position in vector x
        x = y[36:39]
        # matrix A is calculated
        A = AMatrix(x, mu)
        # declares and initializes 6x6 matrix phi and fills it with components of vector y
        phi = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                phi[i,j] = y[6*i+j]
        # matrix caluclation of 6x6 matrices A and phi
        res = A.dot(phi)
        # declares and initializes vector a and fills it with components of 6x6 matrix APhi
        a = np.zeros(36)
        for i in range(6):
            for j in range(6):
                a[6*i+j] = res[i,j]
        # calculates r1, r2 and the state vector dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
        r1 = np.sqrt((y[36]+mu)**2 + y[37]**2 + y[38]**2)
        r2 = np.sqrt((y[36]-(1-mu))**2 + y[37]**2 + y[38]**2)
        ydot = np.array([y[39],
                         y[40],
                         y[41],
                         y[36] + 2*y[40] - ((1-mu)*(y[36]+mu))/(r1**3) - (mu*(y[36]-(1-mu)))/(r2**3),
                         y[37] - 2*y[39] - ((1-mu)*y[37])/(r1**3) - (mu*y[37])/(r2**3),
                         - ((1-mu)*y[38])/(r1**3) - (mu*y[38])/(r2**3)])
        # transposes horizontal vectors a and c to vertical vectors
        a = np.vstack(a)
        ydot = np.vstack(ydot)
        # declares and initializes vector APhi and fills it with components of matrix APhi
        APhi = np.zeros(42)
        APhi = np.vstack((a,ydot))

        return APhi

    elif len(y) == 6:
        # calculates r1, r2 and the state vector dx/dt = [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt]^T
        r1 = np.sqrt((y[0]+mu)**2 + y[1]**2 + y[2]**2)
        r2 = np.sqrt((y[0]-(1-mu))**2 + y[1]**2 + y[2]**2)
        ydot = np.array([y[3],
                         y[4],
                         y[5],
                         y[0] + 2*y[4] - ((1-mu)*(y[0]+mu))/(r1**3) - (mu*(y[0]-(1-mu)))/(r2**3),
                         y[1] - 2*y[3] - ((1-mu)*y[1])/(r1**3) - (mu*y[1])/(r2**3),
                         - ((1-mu)*y[2])/(r1**3) - (mu*y[2])/(r2**3)])

        return ydot

    else:
        print("Dimension of input parameter y is not supported.")




#--------------------------------------------------------------------------
# STATE TRANSITION MATRIX
#
# DESCRIPTION:      The State Transition Matrix (STM) is a linear map from the initial state
#                   at the initial time t0 to a state at some later time t and therefore
#                   offers a tool to approximate the impact of variations in the initial
#                   state on the evolution of the trajectory
#                   The STM is defined by the following equations:
#
#                             x(t) = Phi(t,t0) * x(t0)
#                    dPhi(t,t0)/dt = A(t) * Phi(t,t0)
#                        Phi(t,t0) = I6
#
# OUTPUT:           Output is the 6x6 State Transition Matrix
#
# SYNTAX:           stm(x, t0, tf, mu)
# x             =   Initial State Guess
# t0            =   Start time
# tf            =   Final time
# mu            =   Mass ratio of Primaries
#--------------------------------------------------------------------------

def stm(x, t0, tf, mu):

    # defines number of steps for the numerical ODE solver
    n = 100
    # declares and initializes initial state vector including STM data
    y0 = np.zeros(42)
    # declares and initializes the initial STM as a 6x6 unit matrix
    I = np.eye(6)
    # fills vector y with STM components and the initial state x
    for i in range(6):
        y0[36 + i] = x[i]
        for j in range(6):
            y0[6*i+j] = I[i,j]
    # numerically integrates the system of ODEs
    Y = numMethods.rk4System(y0, t0, tf, mu)
    # gets dimension of the matrix Y and stores number of rows in r
    d = np.shape(Y)
    r = d[0]
    # declares and initializes vector including the STM data of the last time step
    y = np.zeros((r,36))
    y = Y[(r-1),0:36]
    # declares and initiazlizes matrix phi and fills it with the elements of q
    phi = np.zeros((6,6))
    for i in range(6):
        for j in range(6):
            phi[i,j] = y[6*i+j]

    return phi




#--------------------------------------------------------------------------
# HALF PERIOD
#
# DESCRIPTION:      The equations of motion are integrated until y changes sign. Then the step
#                   size is reduced and the integration goes forward again starting at the last
#                   point before the change of the sign. This is repeated until abs(y) < epsilon.
#
# OUTPUT:           Value of half period time
#
# SYNTAX:           halfPeriod(x0, t0, mu, epsilon)
# x0            =   Initial State
# t0            =   Start time
# mu            =   Mass ratio of Primaries
# epsilon       =   Error Tolerance
#--------------------------------------------------------------------------

def halfPeriod(x0, t0, mu, epsilon):
    print("        Calculation of T/2...")
    # declares and initializes actual state x and state at half period xHalfPeriod
    if x0[4] >= 0:
        x = np.ones(6)
        xHalfPeriod = np.ones(6)
    else:
        x = - np.ones(6)
        xHalfPeriod = - np.ones(6)
    # step size
    stepSize = 0.5
    timeStep = t0 + stepSize
    # number of steps for numerical integration
    n = 2000
    # counter for step size reductions
    counter = 0

    while abs(xHalfPeriod[1]) > epsilon:

        if x0[4] >= 0:
            while x[1] > 0:
                xHalfPeriod = x
                timeStep = timeStep + stepSize
                Y = numMethods.rk4System(x0, t0, timeStep, mu)
                x = Y[n,:]
                #print("y = %10.8e at t = %10.8e" % (xHalfPeriod[1], timeStep))
        else:
            while x[1] < 0:
                xHalfPeriod = x
                timeStep = timeStep + stepSize
                Y = numMethods.rk4System(x0, t0, timeStep, mu)
                x = Y[n,:]
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
    print("        Step size has been reduced for %2d times:" % (counter))
    print("        -> y = %10.8e at T/2 = %10.8e\n" % (xHalfPeriod[1], tHalfPeriod))

    return tHalfPeriod




#--------------------------------------------------------------------------
# JACOBI CONSTANT
#
# DESCRIPTION:      Calculates the Jacobi Constant of a state x. If x is a solultion to the
#                   CR3BP that leads to a periodic halo orbit, the Jacobi Constant of x can
#                   be seen as the Jacobi constant related to the orbit.
#
# OUTPUT:           Value of Jacobi constant
#
# SYNTAX:           jacobiConst(x, mu)
# x             =   State vector
# mu            =   Mass ratio of Primaries
#--------------------------------------------------------------------------

def jacobiConst(x, mu):

    r1 = np.sqrt((x[0]+mu)**2 + x[1]**2 + x[2]**2)
    r2 = np.sqrt((x[0]-(1-mu))**2 + x[1]**2 + x[2]**2)
    C = -1/2 * (x[3]**2 + x[4]**2 + x[5]**2) + 2 * (1/2 * (x[0]**2 + x[1]**2) + (1-mu)/r1 + mu/r2)

    return C




# NOT FINISHED YET
def stability(eigenvalues):

    # array for reciprocal pairs of eigenvalues
    reciprocal = np.zeros((2,2), dtype=np.complex)
    # array for trivial pair of eigenvalues
    trivial = np.zeros((1,2), dtype=np.complex)
    # array for complex stability indices
    complexStabilityIndices = np.zeros((1,2), dtype=np.complex)
    # array for stability indices in modulus
    stabilityIndices = np.zeros(2)
    # counter for reciprocal array
    counter = 0
    # loops through eigenvalues and compares each possible pair one single time
    for i in range(6):
        for j in range(5-i):
            # checks if acutal pair if eigenvalues is reciprocal
            if eigenvalues[i].round(3) == (1/eigenvalues[5-j]).round(3):
                # sorts out trivial pair of eigenvalues and stores them
                if  0.99 <= np.sqrt(eigenvalues[i].real**2 + eigenvalues[i].imag**2)  <= 1.01:
                    trivial[0, 0] = eigenvalues[i]
                    trivial[0, 1] = eigenvalues[5-j]
                # sorts out reciprocal pair of eigenvalues and stores them
                else:
                    reciprocal[counter, 0] = eigenvalues[i]
                    reciprocal[counter, 1] = eigenvalues[5-j]
                    counter = counter + 1
    print(reciprocal)
    # calculates stability indices
    for i in range(2):
        complexStabilityIndices[0, i] = 1/2 * (reciprocal[i, 0] + 1/reciprocal[i, 1])
        stabilityIndices[i] = np.sqrt(complexStabilityIndices[0, i].real**2 + complexStabilityIndices[0, i].imag**2)
    print(complexStabilityIndices)
    print(stabilityIndices)

    if stabilityIndices[0] <= 1 and stabilityIndices[1] <= 1:
        return True
    else:
        return False
