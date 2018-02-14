
"""
File    : utilities.py
Author  : Victor Hertel
Date    : 12.01.2018


"""

# Imports
import numpy as np
import matplotlib as mpl
# necessary to use matplotlib for Mac
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from numpy.linalg import svd



def plotTraj(ax):
    #ax.legend()
    ax.set_xlabel("x Axis")
    ax.set_ylabel("y Axis")
    ax.set_zlabel("z Axis")
    setAxesEqual(ax)
    plt.show()

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
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

def colormapLegend():
    normalize = mpl.colors.Normalize(vmin=3.02, vmax=3.2)
    colormap = plt.cm.jet
    scalarmappaple = plt.cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(5)
    plt.colorbar(scalarmappaple)

def saveFig(fig, ax, numOfFigures, format):
    theta = np.linspace(0, 2*np.pi, numOfFigures)
    if format == "png":
        for i in range(0, numOfFigures, 1):
             ax.view_init(elev = 0 + 20 * np.sin(theta[i]), azim = (i*(360/numOfFigures) + 270))
             fig.savefig("fig%d.png" % (i), format='png', dpi=500, bbox_inches = 'tight')
             print("Figure %2d has been saved." % (i))

    elif format == "pdf":
        for i in range(0, numOfFigures, 1):
             ax.view_init(elev = 0 + 20 * np.sin(theta[i]), azim = (i*(360/numOfFigures) + 270))
             fig.savefig("fig%d.pdf" % (i), format='pdf', dpi=500, bbox_inches = 'tight')
             print("Figure %2d has been saved." % (i))

    else:
        print("Format is not supported.")
