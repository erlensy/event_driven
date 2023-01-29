import numpy as np
import multiprocessing
from multiprocessing import cpu_count, Process

def plot_frame(i_, frames, colors):
    from matplotlib import pyplot as plt, collections
    fig = plt.figure() 
    ax = plt.axes()
    x, y, r = frames[i_, :, 1], frames[i_, :, 2], frames[i_, :, 5]
    ax.add_collection(collections.EllipseCollection(
        widths=r*2, heights=r*2, angles=0, units='xy', offsets=np.dstack((x, y))[0], transOffset=ax.transData, 
        facecolors=plt.cm.brg(colors)))
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    plt.savefig(f"../figures/{i_}.jpeg", dpi = 300)
    plt.close()

def read_file(filename):
    print("reading files...")
    f = np.loadtxt(filename, dtype = float, skiprows = 6)
    N = int(np.loadtxt(filename, skiprows = 2, max_rows = 1)[0])
    N_frames = int(len(f) / N)
    frames = np.zeros((N_frames, N, 7))
    for i in range(N_frames):
        for j in range(N):
            frames[i, j, :] = f[i*N + j, :]
    return frames, N, N_frames

if __name__ == "__main__":
    filename = f"../data/data.txt"
    frames, N, N_frames = read_file(filename)
    colors = np.linspace(1.0/N, 1.0, N)
    print("plotting files...")
    cpus = cpu_count()

    i = 0
    while i < N_frames:
        processes = []
        for j in range(cpus):
            process = Process(target = plot_frame, args = (i, frames, colors))
            processes.append(process)
            process.start()
            i += 1
            if i >= N_frames:
                break
    
        for process in processes:
            process.join()
            process.close()
