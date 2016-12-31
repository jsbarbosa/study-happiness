import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('simulated.dat')

size = int(np.sqrt(data.shape[0]/2))

first = data[:size**2]
second = data[size**2:]

magnets_names = ["single.png", "double.png"]
magnets = [first, second]

for (name, magnet) in zip(magnets_names, magnets):
    potential = magnet.reshape(size, size)

    xi, yi = np.linspace(0, size-1, size), np.linspace(0, size-1, size)
    xi, yi = np.meshgrid(xi, yi)

    dx, dy = np.gradient(-potential)

    plt.imshow(potential)
    plt.streamplot(xi, yi, dy, dx, color='black', arrowsize=3)

    plt.xlim(0, size-1)
    plt.ylim(size-1, 0)

    labels = np.arange(6, dtype=int)*10/5
    ticks = np.arange(6)*size/5
    ticks[-1] += -1
    
    plt.xticks(ticks, labels)
    plt.yticks(ticks[::-1], labels)

    plt.xlabel('$x$ (cm)')
    plt.ylabel('$y$ (cm)')
    plt.grid()

    plt.savefig(name, dpi=300)
    plt.close()

