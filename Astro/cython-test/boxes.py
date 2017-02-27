import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random, randint, normal
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from core import box, limits, solver

import sys

def galaxy(N, min_v, max_v, arms = 5):
    a = 0.4
    b = 0.7
    N = int(N/arms)
    t = np.random.randn(N)
    
    x = np.array([0])
    y = np.array([0])

    for i in range(arms):
        xt = a*np.exp(b*t)*np.cos(t+2*i*np.pi/arms) + np.random.normal(0, a*0.2, N)
        yt = a*np.exp(b*t)*np.sin(t+2*i*np.pi/arms) + np.random.normal(0, a*0.2, N)
        x = np.hstack((x, xt))
        y = np.hstack((y, yt))
        
    rank = max_v - min_v
    return x[1:]*rank + min_v, y[1:]*rank + min_v

def update(i):
    ax.view_init(30, i*360/frames)

def minimum(value):
    if value > 0:
        return 0.9*value
    else:
        return value*1.1
    
def random_generator(min_value, max_value, N = 10):
    return (max_value-min_value)*random(N) + min_value

def speeds_generator(distances):
    r = np.sqrt(np.sum(distances**2, axis=0))
    n = len(r)
    s = np.zeros((3, n))
    for i in range(n):
        pos = np.where(r < r[i])
        M = sum(masses[pos])
        s[:, i] = np.sqrt(G*M/r[i])
#    s = np.sqrt(G/r)
    e1 = np.zeros((3, n))
    e1[:-1] = -distances[1,:], distances[0,:]
    e2 = np.array([np.cross(distances[:,i], e1[:,i]) for i in range(n)]).T
    n1 = e1/np.sqrt(np.sum(e1**2, axis=0))
    n2 = e2/np.sqrt(np.sum(e2**2, axis=0))
    omega = 2*np.pi*random()
    vel = np.cos(omega)*n1 + np.sin(omega)*n2
    vel *= s
    return vel

N = 1000
N_half = int(N/2)
sys.setrecursionlimit(N*N)

G = 4.302e-3 #pc, solar masses, km/s
density = (N/(0.14*10**(np.log10(N)-10)))**(1/3)
masses = abs(normal(loc = 10, size=(N+2)))
center_mass = 3*max(masses)
masses[-2:] = center_mass
min_value, max_value = -density*0.5, density*0.5
cos45 = np.cos(np.pi/4)
rotx = np.array([[1, 0, 0], [0, cos45, -cos45], [0, cos45, cos45]])

x1, y1 = galaxy(N_half, min_value, max_value) - min_value
x2, y2 = galaxy(N_half, min_value, max_value) + min_value
z1 = random_generator(min_value/10, max_value/10, N_half) - 5*max_value
z2 = random_generator(min_value/10, max_value/10, N_half) + 5*max_value
matrix2 = np.zeros((3, N_half))
matrix2 = x2, y2, z2
center2 = np.mean(matrix2, axis = 1)
distances = np.array([matrix2[i] - center2[i] for i in range(3)])
speeds2 = speeds_generator(distances)

matrix2 = np.dot(rotx, matrix2)
speeds2 = np.dot(rotx, speeds2)

speeds = np.zeros((3, N+2))
matrix = np.zeros((3, N+2))
matrix[:, :N_half] = x1, y1, z1
matrix[:, N_half:-2] = matrix2

center1 = np.mean(matrix[:, :N_half], axis = 1)
center2 = np.mean(matrix[:, N_half:], axis = 1)
matrix[:, -2] = center1
matrix[:, -1] = center2
distances = np.array([matrix[i, :N_half] - center1[i] for i in range(3)])
speeds_ = speeds_generator(distances)
speeds[:2, :N_half] = speeds_[:2]
speeds[:, N_half:-2] = speeds2
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
#ax.plot([center1[0], center2[0]], [center1[1], center2[1]], [center1[2], center2[2]], "o", ms=20)
#ax.set_axis_off()
plot = ax.plot(matrix[0], matrix[1], matrix[2], "o", ms=0.4, c="black")[0]
#plot2 = ax.plot(x2, y2, z2, "o", ms=0.4, c="r")[0]
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")

np.savetxt("Data/0_instant.dat", matrix)

n = solver(matrix, speeds, masses, N, 2e6, 1e4, G)

frames = 100
ratio = int(n/frames)
temp = np.arange(0, n, ratio).astype(int)

def update(i):
    temp = np.genfromtxt("Data/%d_instant.dat"%i)
    plot.set_data(temp[:,0], temp[:,1])
    plot.set_3d_properties(temp[:, 2])

min_value *= 5
max_value *= 5
ax.set_xlim(min_value, max_value)
ax.set_ylim(min_value, max_value)
ax.set_zlim(min_value, max_value)

ani = FuncAnimation(fig, update, frames=frames, interval=frames)
ani.save("temp.gif", writer='imagemagick')
#plt.show()
