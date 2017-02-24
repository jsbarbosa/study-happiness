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
        self.r0 = r0
        self.r1 = [r0[i] + rs[i] for i in range(3)]
        self.steps = [0.5*rs[i] for i in range(3)]
        self.halfs = [r0[i] + self.steps[i] for i in range(3)]        
        self.tx = self.r0[0], self.halfs[0]
        self.ty = self.r0[1], self.halfs[1]
        self.tz = self.r0[2], self.halfs[2]
        
        self.mass = sum(self.masses)
        self.center_mass = self.center_mass_calc()
        self.length = sum(rs)
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
            if theta > 1:
                cm, m = ['a'], []
                i = 0
                for child in self.children:
                    ans = child.force(point)
                    if ans != None:
                        if 'a' in ans[0]:
                            cm += ans[0][1:]
                        else:
                            cm.append(ans[0])
                        m += ans[1]
                        i += 1
                if i > 0:
                    if not self.root:
                        return cm, m
                    else:
                        return cm[1:], np.array(m)
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
            low = self.points[i] <= half
            high = self.points[i] > half
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
N = 25
sys.setrecursionlimit(N*N)
def random_generator(min_value, max_value, N = 10):
    return (max_value-min_value)*random(N) + min_value

G = 4.302e-3 #pc, solar masses, km/s
min_value, max_value = -20, 20
x = random_generator(min_value, max_value, N)
y = random_generator(min_value, max_value, N)
z = random_generator(min_value, max_value, N)
vx = random_generator(-1, 1, N)/10
vy = random_generator(-1, 1, N)/10
vz = random_generator(-1, 1, N)/10
matrix = np.zeros((3, N))
speeds = np.zeros((3, N))
matrix[0] = x
matrix[1] = y
matrix[2] = z
speeds[0] = vx
speeds[1] = vy
speeds[2] = vz
masses = abs(normal(loc = 10, size=N))
data = [x, y, z]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
#ax.set_axis_off()
plot = ax.plot(x, y, z, "o", ms=1, c="black")[0]
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
def update(i):
    ax.view_init(30, i*360/frames)

def magnitude(pos):
    return np.dot(pos, pos)**(3/2)

def iterator(particle):
    global tree
    ans = tree.force(particle)
    if ans is not None:
        distance, masses = ans
        range_n = range(len(masses))
        mag = np.array(list(map(magnitude, distance)))
        ratio = masses/mag
        acceleration = np.array([ratio[i]*np.array(distance[i]) for i in range_n])
    #    acceleration = np.dot(ratio, distance)
        acceleration = G*np.sum(acceleration, axis = 0)
        maximum = max(abs(acceleration))
        if maximum > 25:
            return 5*acceleration/maximum
        return acceleration
    else:
        return [0, 0, 0]

def minimum(value):
    if value > 0:
        return 0.9*value
    else:
        return value*1.1

def limits(matrix):
    x_min = [min(matrix[:,i])*0.9 if min(matrix[:,i]) > 0 else min(matrix[:,i])*1.1 for i in range(3)]
    x_max = [1.1*max(matrix[i, :]) - 1.1*x_min[i] for i in range(3)]
    return x_min, x_max
    
def solver(positions, speeds, t, dt):
    global tree
    n = int(t/dt)
    positions_in_time = np.zeros((n, N, 3))
    x = positions.T.copy()
    v = speeds.T.copy()
    positions_in_time[0] = x
        
    for i in range(n-1):
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)        
        accelerations = np.array(list(map(iterator, x)))
        v_half = v + 0.5*dt*accelerations
        x += dt*v_half
        
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)
        accelerations = np.array(list(map(iterator, x)))
        v = v_half + 0.5*dt*accelerations
        positions_in_time[i+1] = x
        
    return positions_in_time, n

pos, n = solver(matrix, speeds, 1000, 1)
frames = 100
ratio = int(n/frames)
temp = np.arange(0, n, ratio).astype(int)
pos = pos[temp]

def update(i):
    temp = pos[i]
    plot.set_data(temp[:,0], temp[:,1])
    plot.set_3d_properties(temp[:, 2])


ax.set_xlim(1.5*min_value, 1.5*max_value)
ax.set_ylim(1.5*min_value, 1.5*max_value)
ax.set_zlim(1.5*min_value, 1.5*max_value)

ani = FuncAnimation(fig, update, frames=frames, interval=frames)
ani.save("concept.gif", writer='imagemagick')
plt.show()