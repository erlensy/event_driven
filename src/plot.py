from matplotlib import pyplot as plt, collections
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import numpy as np

filename = f"../data/skrt.txt"
f = np.loadtxt(filename, skiprows = 6)
N = int(np.loadtxt(filename, skiprows = 2, max_rows = 1)[0])
N_frames = int(len(f) / N)
frames = np.zeros((N_frames, N, 6))
color = np.random.rand(N)

for i in range(N_frames):
    for j in range(N):
        frames[i, j, :] = f[i*N + j, :]

for i in range(N_frames):
    fig, ax = plt.subplots()
    x, y, r = frames[i, :, 0], frames[i, :, 1], frames[i, :, 4]
    ax.add_collection(collections.EllipseCollection(
        widths=r*2, heights=r*2, angles=0, units='xy', offsets=list(zip(x, y)), transOffset=ax.transData, 
        facecolors=plt.cm.brg(color)))
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.grid()
    plt.savefig(f"../figures/{i}.pdf")
    plt.close()

#def animate(i):
#    ax.clear()
#    i_ = round(N_frames / nof * i)
#    
#    x, y, r = frames[i_, :, 0], frames[i_, :, 1], frames[i_, :, 4]
#    ax.add_collection(collections.EllipseCollection(
##        widths=r*2, heights=r*2, angles=0, units='xy', color = "black", offsets=list(zip(x, y)), transOffset=ax.transData))
#
#nof = 500
#fig, ax = plt.subplots()
#ax.set_xlim(0, 1.0)
#ax.set_ylim(0, 1.0)
##ax.grid()
#anim = FuncAnimation(fig, animate, frames = nof, interval = 25)
#anim.save(r"../figures/test.gif", writer=PillowWriter())
