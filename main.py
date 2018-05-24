"""
File    : main.py
Author  : Victor Hertel
Date    : 24.04.2018

Computation of halo orbits
"""

# Imports
import numpy as np
from Orbit import Orbit
from OrbitFamily import OrbitFamily
from Utility import Plot



# Initial Conditions
l1 = np.array([0.8233901862, 0, -0.0029876370, 0, 0.1264751431, 0])
l2 = np.array([1.1808881373, 0, -0.0032736457, 0, -0.1559184478, 0])
l2_middle = np.array([1.1542349115, 0, -0.1379744940, 0, -0.2147411949, 0])
mFirstPrimary = 5.97237e+24  # Earth
mSecondPrimary = 7.342e+22  # Moon
mu = mSecondPrimary / (mSecondPrimary + mFirstPrimary)


nro = OrbitFamily(l2, 0.007, 10, mu)
nro.getNRHOFamily()
nro.plot(haloFamily="southern")


#Plot.plotFromTable("L1&L2", orbits=True, haloFamily="southern", orbitNumber=100)


























#Writer = animation.writers["ffmpeg"]
#writer = Writer(fps=15, bitrate=1800)

# def randrange(n, vmin, vmax):
#    return (vmax - vmax) * np.random.rand(n) + vmin
# n = 100
# xx = randrange(n, 23, 32)
# yy = randrange(n, 0, 100)
# zz = randrange(n, -50, -25)
#
# fig = plt.figure()
# ax = Axes3D(fig)
#
# def init():
#     ax.scatter(xx, yy, zz, marker="o", s=20, c="goldenrod", alpha=0.6)
#     return fig,
# def animate(i):
#     ax.view_init(elev=10, azim=i)
#     return fig,
#
# anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
# anim.save("basic.mp4")

# def make_views(ax, angles, elevation=None, width=4, height=3, prefix="tmprot_", **kwargs):
#     files = []
#     ax.figure.set_size_inches(width, height)
#
#     for i, angle in enumerate(angles):
#         ax.view_init(elev=elevation, azim=angle)
#         fname = "%s%03d.png" % (prefix,i)
#         ax.figure.savefig(fname)
#         files.append(fname)
#     return files
#
#
# def make_movie(files, output, fps=10, bitrate=1800, **kwargs):
#     output_name, output_ext = os.path.splitext(output)
#     print(files)
#     command = {'.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\ -lavcopts vcodec=msmpeg4v2:vbitrate=%d' % (",".join(files), fps, output_name,bitrate)}
#     print(command[output_ext])
#     output_ext = os.path.splitext(output)[1]
#     os.system(command[output_ext])
#
#
#
# def rotanimate(ax, angles, output, **kwargs):
#     output_ext = os.path.splitext(output)[1]
#
#     files = make_views(ax, angles, **kwargs)
#
#     D = {'.mp4' : make_movie}
#     D[output_ext](files, output, **kwargs)
#     for f in files:
#         os.remove(f)
#
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# X, Y, Z = axes3d.get_test_data(0.05)
# s = ax.plot_surface(X, Y, Z, cmap=cm.jet)
# plt.axis("off")
# angles = np.linspace(0, 360, 21)[:-1]
# rotanimate(ax, angles, "movie.mp4", fps=10, bitrate=2000)
