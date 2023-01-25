from matplotlib import pyplot as plt, collections
import numpy as np
import os

os.system("../build/event_diven")

for i in range(10):
    filename = f"../data/{i}.txt"
    f = np.loadtxt(filename)
    x, y, r = f[:, 0], f[:, 1], f[:, 4]
    ax = plt.axes() 
    ax.add_collection(collections.EllipseCollection(
        widths=r*2, heights=r*2, angles=0, units='xy', color = "black", offsets=list(zip(x, y)), transOffset=ax.transData))
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.grid()
    plt.show()
