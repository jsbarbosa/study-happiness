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
        xt = a*np.exp(b*t)*np.cos(t+2*i*np.pi/arms)# + np.random.normal(0, a*0.2, N)
        yt = a*np.exp(b*t)*np.sin(t+2*i*np.pi/arms)# + np.random.normal(0, a*0.2, N)
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
        M = n
        if M != 0:
            cm = np.sum(distances[:, pos], axis=-1)/M
            temp = r[i] - cm
            temp = np.sqrt(np.sum(temp**2, axis=0))
            s[:, i] = np.sqrt(G*M/temp)
        else:
            s[:, i] = 0
    e1 = np.zeros((3, n))
    e1[:-1] = -distances[1,:], distances[0,:]
    e2 = np.array([np.cross(distances[:,i], e1[:,i]) for i in range(n)]).T
    n1 = e1/np.sqrt(np.sum(e1**2, axis=0))
    n2 = e2/np.sqrt(np.sum(e2**2, axis=0))
    omega = np.pi#2*np.pi*random()
    vel = np.cos(omega)*n1 + np.sin(omega)*n2
    vel *= s
    return vel

N = 1000
N_half = int(N/2)
sys.setrecursionlimit(N*N)

G = 44.97 #kpc, e7 solar masses, Gy
mass = 1
cos45 = np.cos(np.pi/4)
rotx = np.array([[1, 0, 0], [0, cos45, -cos45], [0, cos45, cos45]])
min_value, max_value = -50, 50

x1, y1 = galaxy(N_half, min_value, max_value)
x2, y2 = galaxy(N_half, min_value, max_value)
z1 = random_generator(min_value/10, max_value/10, N_half) - 5*max_value
z2 = random_generator(min_value/10, max_value/10, N_half) + 5*max_value
   
galaxy1 = np.zeros((3, N_half))
galaxy2 = np.zeros((3, N_half))
galaxy1[:] = x1, y1, z1
galaxy2[:] = x2, y2, z2       
#galaxy1 *= 0.25
#galaxy2 *= 0.25

center1 = np.mean(galaxy2, axis = 1)
distances = np.array([-galaxy1[i] + center1[i] for i in range(3)])
speeds1 = speeds_generator(distances)

center2 = np.mean(galaxy2, axis = 1)
distances = np.array([-galaxy2[i] + center2[i] for i in range(3)])
speeds2 = speeds_generator(distances)

galaxy2 = np.dot(rotx, galaxy2)
speeds2 = np.dot(rotx, speeds2)

speeds = np.zeros((3, N))
system = np.zeros((3, N))
system[:, :N_half] = galaxy1
system[:, N_half:] = galaxy2
#speeds[:, :N_half] = speeds1
#speeds[:, N_half:] = speeds2

np.savetxt("Data/0_instant.dat", system)

n = solver(system, speeds, N, 100, 1, G)