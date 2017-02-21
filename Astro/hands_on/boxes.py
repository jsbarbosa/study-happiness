import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random, randint
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import mpl_toolkits.mplot3d.art3d as art3d

from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as patches

class box():
    def __init__(self, points, x0, y0, z0, xs, ys, zs, buffer = None):
        global count, ax
        self.points = points
        self.stop = False
        self.status()
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.x1 = x0 + xs
        self.y1 = y0 + ys
        self.z1 = z0 + zs
        self.xsize = xs
        self.ysize = ys
        self.zsize = zs
        self.xstep = 0.5*xs
        self.ystep = 0.5*ys
        self.zstep = 0.5*zs
        self.xhalf = x0 + self.xstep
        self.yhalf = y0 + self.ystep
        self.zhalf = z0 + self.zstep
        self.tx = x0, self.xhalf
        self.ty = y0, self.yhalf
        self.tz = z0, self.zhalf
        self.children = np.zeros((2, 2, 2), dtype=box)
#        self.children = [0 for i in range(8)]
        
        if not self.stop:
            self.plot()
            count += 1
            self.checksons()
            self.assign()
            
    def plot(self):
        xy, yx = np.meshgrid([self.x0, self.x1], [self.y0, self.y1])
        xz, zx = np.meshgrid([self.x0, self.x1], [self.z0, self.z1])
        yz, zy = np.meshgrid([self.y0, self.y1], [self.z0, self.z1])
        self.surface(xy, yx, self.z1)
        self.surface(xy, yx, self.z0)
        self.surface(xz, self.y0, zx)
        self.surface(xz, self.y1, zx)
        self.surface(self.x1, yz, zy)
        self.surface(self.x0, yz, zy)
            
    def surface(self, x, y, z):
        a, c, w = 0.1, "b", 0.5
        ax.plot_surface(x, y, z, alpha = a, color = c, linewidth = w, edgecolor=c)
        
    def status(self):
        stops = sum([len(item) for item in self.points])
        if stops == 3:
            self.stop = True

    def assign(self):
        for i in range(8):
            j, k, l = i%2, int(i/2)%2, int(i/4)
            pos = self.inter[i]
            if len(pos) > 0:
                x = self.points[0][pos]
                y = self.points[1][pos]
                z = self.points[2][pos]
                self.children[j, k, l] = box([x, y, z], self.tx[j], self.ty[k],
                             self.tz[l], self.xstep, self.ystep, self.zstep) 
        
    def checksons(self):
        pos = []
        values = [self.xhalf, self.yhalf, self.zhalf]
        for i in range(3):
            higher = np.where(values[i] >= self.points[i])[0].astype(int)
            lower = np.where(values[i] < self.points[i])[0].astype(int)
            pos.append(higher)
            pos.append(lower)
        
        self.inter = [np.intersect1d(np.intersect1d(pos[i%2], pos[int(i/2)%2+2]), pos[int(i/4)+3]).astype(int) for i in range(8)]
#        for i in range(8):
#            j, k, l = i%2, int(i/2)%2, int(i/4)
#            inter = np.intersect1d(np.intersect1d(pos[j], pos[k+2]), pos[l+3]).astype(int)
#            self.inter.append(inter)
#
def random_generator(min_value, max_value, N = 10):
    return (max_value-min_value)*random(N) + min_value
    
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
N = 10
count = 0
x = random_generator(1, 2, N)
y = random_generator(1, 2, N)
z = random_generator(1, 2, N)
data = [x, y, z]
a = box(data, min(x), min(y), min(z),
        max(x) - min(x), max(y) - min(y),
           max(z) - min(z))
ax.plot(x, y, z, "o", ms=2)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")
#plt.savefig("only.png")
plt.show()

