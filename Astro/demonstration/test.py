from temp import box, solver
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


N = 100
p = np.zeros((3, N))
speeds = np.zeros_like(p)
tree = []
for i in range(3):
    p[i] = np.random.random(N)

#p[0] = 0.5, 1.5
#p[1] = 0, 0
cb = np.zeros(3)
r = 2*np.ones(3)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
ax.plot(p[0], p[1], p[2], "o", c="r", ms=1)

n = solver(p, speeds, N, 10, 0.01, 1, filename='data/')
