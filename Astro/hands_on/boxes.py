import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random, randint
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import mpl_toolkits.mplot3d.art3d as art3d

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as patches

from time import sleep


class box():
    def __init__(self, points, r0, rs):
        global ax
        self.points = points
        print(len(points[0]))
        self.r0 = r0
        self.r1 = [r0[i] + rs[i] for i in range(3)]
        self.steps = [0.5*rs[i] for i in range(3)]
        self.halfs = [r0[i] + self.steps[i] for i in range(3)]
        
        self.tx = self.r0[0], self.halfs[0]
        self.ty = self.r0[1], self.halfs[1]
        self.tz = self.r0[2], self.halfs[2]
        self.children = [None for i in range(8)]
#        self.transitions = [[self.r0[i], self.halfs[i]] for i in range(3)]
        self.sons_points()
        self.plot()
        
    def plot(self):        
        xy, yx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[1], self.r1[1]])
        xz, zx = np.meshgrid([self.r0[0], self.r1[0]], [self.r0[2], self.r1[2]])
        yz, zy = np.meshgrid([self.r0[1], self.r1[1]], [self.r0[2], self.r1[2]])
        self.surface(xy, yx, self.r1[2])
        self.surface(xy, yx, self.r0[2])
        self.surface(xz, self.r0[1], zx)
        self.surface(xz, self.r1[1], zx)
        self.surface(self.r1[0], yz, zy)
        self.surface(self.r0[0], yz, zy)
            
    def surface(self, x, y, z):
        a, c, w = 0.05, "b", 0.5
        ax.plot_surface(x, y, z, alpha = a, color = c, linewidth = w, edgecolor=c)
        
    def sons_points(self):
        data = []
        for i in range(3):
            r1 = self.r1[i]
            r0 = self.r0[i]
            half = self.halfs[i]
            
            low = (self.points[i] < half) & (self.points[i] >= r0)
            high = (self.points[i] >= half) & (self.points[i] <= r1)
            data.append(low)
            data.append(high)
        for i in range(8):
            j = i%2
            k = int(i/2)%2
            l = int(i/4)
            inter = data[j] & data[k + 2] & data[l + 4]
            inter = np.where(inter)[0]                    
            x = self.points[0][inter]
            y = self.points[1][inter]
            z = self.points[2][inter]
            if len(inter) > 0 and len(self.points[0]) > 1:
                self.children[i] = box([x, y, z], [self.tx[j], self.ty[k],
                         self.tz[l]], self.steps) 

N = 10

def random_generator(min_value, max_value, N = 10):
    return (max_value-min_value)*random(N) + min_value

min_x, max_x = 1, 2
min_y, max_y = 1, 2
min_z, max_z = 1, 2
x = random_generator(min_x, max_x, N)
y = random_generator(min_y, max_y, N)
z = random_generator(min_z, max_z, N)

data = [x, y, z]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
ax.plot(x, y, z, "o", c="r")
a = box(data, [min_x, min_y, min_z], [max_x-min_x, max_y-min_y, max_z-min_z])

#temp = box([x, y, z], [min(x), min(y), min(z)], [2, 2, 2])

plt.show()
#class box():
#    def __init__(self, points, x0, y0, z0, xs, ys, zs, buffer = None):
#        global count, ax
#        self.points = points
#        self.stop = False
#        self.empty = False
#        self.status()
#        self.x0 = x0
#        self.y0 = y0
#        self.z0 = z0
#        self.x1 = x0 + xs
#        self.y1 = y0 + ys
#        self.z1 = z0 + zs
#        self.xsize = xs
#        self.ysize = ys
#        self.zsize = zs
#        self.xstep = 0.5*xs
#        self.ystep = 0.5*ys
#        self.zstep = 0.5*zs
#        self.xhalf = x0 + self.xstep
#        self.yhalf = y0 + self.ystep
#        self.zhalf = z0 + self.zstep
#        print(self.xhalf, self.yhalf, self.zhalf)
#        sleep(0.1)
#        self.tx = x0, self.xhalf
#        self.ty = y0, self.yhalf
##        self.tz = self.zhalf, z0
#        self.tz = z0, self.zhalf
##        self.children = np.zeros((2, 2, 2), dtype=box)
#        self.inter = []
#        self.children = [0 for i in range(8)]
#        
#        if not self.stop and not self.empty:
#            self.plot()
#            count += 1
#            self.checksons()
#            self.assign()
#            
#    def plot(self):
#        xy, yx = np.meshgrid([self.x0, self.x1], [self.y0, self.y1])
#        xz, zx = np.meshgrid([self.x0, self.x1], [self.z0, self.z1])
#        yz, zy = np.meshgrid([self.y0, self.y1], [self.z0, self.z1])
#        self.surface(xy, yx, self.z1)
#        self.surface(xy, yx, self.z0)
#        self.surface(xz, self.y0, zx)
#        self.surface(xz, self.y1, zx)
#        self.surface(self.x1, yz, zy)
#        self.surface(self.x0, yz, zy)
#            
#    def surface(self, x, y, z):
#        a, c, w = 0.1, "b", 0.5
#        ax.plot_surface(x, y, z, alpha = a, color = c, linewidth = w, edgecolor=c)
#        
#    def status(self):
#        stops = sum([len(item) for item in self.points])
#        if stops == 3:
#            self.stop = True            
#        if len(self.points[0]) == 0:
#            self.empty = True
#
#    def assign(self):
#        for i in range(8):
#            j, k, l = i%2, int(i/2)%2, int(4/(i+1))%2 - 1
#            l = int(i/4)
#            pos = self.inter[i]
#            x = self.points[0][pos]
#            y = self.points[1][pos]
#            z = self.points[2][pos]
#            print(x, y, z)
#            self.children[i] = box([x, y, z], self.tx[j], self.ty[k],
#                         self.tz[l], self.xstep, self.ystep, self.zstep) 
##            self.children[j, k, l] = box([x, y, z], self.tx[j], self.ty[k],
##                         self.tz[l], self.xstep, self.ystep, self.zstep) 
#        
#    def checksons(self):
#        pos = []
#        values = [self.xhalf, self.yhalf, self.zhalf]
#        low = [self.x0, self.y0, self.z0]
#        high = [self.x1, self.y1, self.z1]
#        for i in range(3):
#            lower = np.where((self.points[i] <= values[i]) & (self.points[i] > low[i]))[0].astype(int)
#            higher = np.where((self.points[i] > values[i]) & (self.points < high[i]))[0].astype(int)
#            print(lower, higher)
#            pos.append(lower)
#            pos.append(higher)
#
#        
##        self.inter = [np.intersect1d(np.intersect1d(pos[i%2], pos[int(i/2)%2+2]), pos[int(4/(i+1))%2-1+4]).astype(int) for i in range(8)]
#        self.inter = [np.intersect1d(np.intersect1d(pos[i%2], pos[int(i/2)%2+2]), pos[int(i/4)+3]).astype(int) for i in range(8)]
#
##

#    
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_aspect("equal")
#N = 20
#count = 0
#x = np.linspace(1, 2, N)
#y = x
#z = x

#

#ax.plot(x, y, z, "o", ms=2)
#ax.set_xlabel("$x$")
#ax.set_ylabel("$y$")
#ax.set_zlabel("$z$")
##plt.savefig("only.png")
#plt.show()

#