import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random, randint, normal
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

import sys

class box():
    def __init__(self, points, r0, rs, masses, root = False, parent = False, tree = None):
        global ax
        self.points = points
        self.masses = masses
        self.parent = parent
        self.root = root
        self.tolerance = 1
        self.r0 = r0
        self.r1 = [r0[i] + rs[i] for i in range(3)]
        self.steps = [0.5*rs[i] for i in range(3)]
        self.halfs = [r0[i] + self.steps[i] for i in range(3)]        
        self.tx = self.r0[0], self.halfs[0]
        self.ty = self.r0[1], self.halfs[1]
        self.tz = self.r0[2], self.halfs[2]
        
        self.mass = sum(self.masses)
        self.center_mass = self.center_mass_calc()
        self.length = np.sqrt(sum([rs[i]**2 for i in range(3)]))
        if tree == None:
            self.tree = []
        else:
            self.tree = tree
        
        self.numParticles = len(points[0])
        self.children = []
        self.sons_points()
        self.tree.append(self)

    def force(self, point):
        l = self.length
        if self.center_mass[0] != point[0]:
            d = np.sqrt(sum([(self.center_mass[i] - point[i])**2 for i in range(3)]))
            theta = l/d
            if theta > self.tolerance:
                cm, m = [], []
                i = 0
                for child in self.children:
                    ans = child.force(point)
                    if ans != None:
                        if ans[0]:
                            cm += ans[1]
                        else:
                            cm.append(ans[1])
                        m += ans[2]
                        i += 1
                if i > 0:
                    if not self.root:
                        return True, cm, m
                    else:
                        if len(cm) == 1:
                            return cm[0], np.array(m)
                        return cm, np.array(m)
            else:
                if not self.root:
                    return False, self.center_mass - point, [self.mass]
                else:
                    return self.center_mass - point, [self.mass]
            
            
    def center_mass_calc(self):
        return [sum(self.masses*self.points[i])/self.mass for i in range(3)]
    
        
    def plot(self, a = 1, c = "b", w = 0.5, child_only = True):
        if child_only:
            if not self.parent:
                xy, yx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[1], self.r1[1]])
                xz, zx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[2], self.r1[2]])
                yz, zy = np.meshgrid([self.r0[1], self.r1[1]], [self.r0[2], self.r1[2]])
                self.surface(xy, yx, self.r1[2], a, c, w)
                self.surface(xy, yx, self.r0[2], a, c, w)
                self.surface(xz, self.r0[1], zx, a, c, w)
                self.surface(xz, self.r1[1], zx, a, c, w)
                self.surface(self.r1[0], yz, zy, a, c, w)
                self.surface(self.r0[0], yz, zy, a, c, w)
        else:
            xy, yx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[1], self.r1[1]])
            xz, zx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[2], self.r1[2]])
            yz, zy = np.meshgrid([self.r0[1], self.r1[1]], [self.r0[2], self.r1[2]])
            self.surface(xy, yx, self.r1[2], a, c, w)
            self.surface(xy, yx, self.r0[2], a, c, w)
            self.surface(xz, self.r0[1], zx, a, c, w)
            self.surface(xz, self.r1[1], zx, a, c, w)
            self.surface(self.r1[0], yz, zy, a, c, w)
            self.surface(self.r0[0], yz, zy, a, c, w)
            
            
    def surface(self, x, y, z, a, c, w):
        ax.plot_wireframe(x, y, z, alpha=a, color=c, linewidth = w)
        
    def sons_points(self):
        data = []
        for i in range(3):
            half = self.halfs[i]            
            low = (self.points[i] <= half) & (self.points[i] >= self.r0[i])
            high = (self.points[i] > half) & (self.points[i] <= self.r1[i])
            data.append(low)
            data.append(high)
            
        for i in range(8):
            j = i%2
            k = int(i/2)%2
            l = int(i/4)
            inter = data[j] & data[k + 2] & data[l + 4]
            inter = np.where(inter)[0]      
            num = len(inter)  
            points = self.points[:, inter]
            if num > 0 and self.numParticles > 1:
                if num == 1:
                    parent = False
                else:
                    parent = True
                self.children.append(box(points, [self.tx[j], self.ty[k],
                         self.tz[l]], self.steps, self.masses[inter], parent=parent, tree = self.tree))

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

def magnitude(pos):
    r = np.sqrt(np.dot(pos, pos))
    if r < 1000:
        r = 1000
    return r**3

def iterator(particle, tree):
    ans = tree.force(particle)
    if ans is not None:
        distance, masses = ans
        n = len(masses)
        range_n = range(n)
        if n == 1:
            mag = magnitude(distance)
            ratio = masses/mag
            acceleration = G*ratio*distance
        else:            
            mag = np.array(list(map(magnitude, distance)))
            ratio = masses/mag
            acceleration = np.array([ratio[i]*np.array(distance[i]) for i in range_n])
            acceleration = G*np.sum(acceleration, axis = 0)
        return acceleration
    else:
        return [0,0,0]

def minimum(value):
    if value > 0:
        return 0.9*value
    else:
        return value*1.1

def limits(matrix):
    x_min = [min(matrix[:,i])*0.9 if min(matrix[:,i]) > 0 else min(matrix[:,i])*1.1 for i in range(3)]
    x_max = [1.1*max(matrix[i, :]) - 1.1*x_min[i] for i in range(3)]
    rank = [x_max[i] - x_min[i] for i in range(3)]
    return np.array(x_min), np.array(rank)
    
def leapfrog(x, v, dt):
    x0, r0 = limits(x)
    tree = box(x.T, x0, r0, masses, parent=True, root=True)        
    accelerations = np.array(list(map(lambda pos: iterator(pos, tree), x)))
    v_half = v + 0.5*dt*accelerations
    x += dt*v_half
    x0, r0 = limits(x)
    tree = box(x.T, x0, r0, masses, parent=True, root=True)
    accelerations = np.array(list(map(lambda pos: iterator(pos, tree), x)))
    v = v_half + 0.5*dt*accelerations
    return x, tree

def solver(positions, speeds, t, dt):
    n = int(t/dt)
    positions_in_time = np.zeros((n, N+2, 3))
    x = positions.T.copy()
    v = speeds.T.copy()
    
    positions_in_time[0] = x
        
    for i in range(n-1):
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)        
        accelerations = np.array(list(map(lambda pos: iterator(pos, tree), x)))
        v_half = v + 0.5*dt*accelerations
        x += dt*v_half
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)
        accelerations = np.array(list(map(lambda pos: iterator(pos, tree), x)))
        v = v_half + 0.5*dt*accelerations
        positions_in_time[i+1] = x
        
    return positions_in_time, n

def random_generator(min_value, max_value, N = 10):
    return (max_value-min_value)*random(N) + min_value

def speeds_generator(distances):
    r = np.sqrt(np.sum(distances**2, axis=0))
    s = np.sqrt(G*center_mass/r)
    n = len(r)
    e1 = np.zeros((3, n))
    e1[:-1] = -distances[1,:], distances[0,:]
    e2 = np.array([np.cross(distances[:,i], e1[:,i]) for i in range(n)]).T
    n1 = e1/np.sqrt(np.sum(e1**2, axis=0))
    n2 = e2/np.sqrt(np.sum(e2**2, axis=0))
    omega = 2*np.pi*random()
    vel = np.cos(omega)*n1 + np.sin(omega)*n2
    vel *= s
    return vel

N = 100
N_half = int(N/2)
sys.setrecursionlimit(N*N)

G = 4.302e-3 #pc, solar masses, km/s
density = (N/(0.14*10**(np.log10(N)-10)))**(1/3)
masses = abs(normal(loc = 10, size=(N+2)))
center_mass = 100*max(masses)
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

np.savetxt("Data/0_data.dat", matrix)
pos, n = solver(matrix, speeds, 2e6, 1e4)
for i, data in enumerate(pos):
    np.savetxt("Data/%d_data.dat"%(i+1), data)

frames = 100
ratio = int(n/frames)
temp = np.arange(0, n, ratio).astype(int)
pos = pos[temp]

def update(i):
    temp = pos[i]
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
