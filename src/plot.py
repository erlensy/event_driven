from matplotlib import pyplot as plt, collections
import numpy as np

f = np.loadtxt("../data/skrt.txt")
x, y, r = f[:, 0], f[:, 1], f[:, 4]

fig , ax = plt.figure(), plt.axes()
ax.add_collection(collections.EllipseCollection(
    widths=r, heights=r, angles=0, units='xy', color = "black", offsets=list(zip(x, y)), transOffset=ax.transData))
ax.grid(); ax.set_xlim(0, 1.0); ax.set_ylim(0, 1.0); plt.show()

