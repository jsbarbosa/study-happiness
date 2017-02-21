import numpy as np
from numpy.random import random, randint
import matplotlib.pyplot as plt
import time

class box():
    def __init__(self, points, x0, y0, z0, xs, ys, zs, buffer = None):
        global count
        self.points = points
        self.empty = False
        self.stop = False
        self.status()
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.xstep = 0.5*xs
        self.ystep = 0.5*ys
        self.zstep = 0.5*zs
        self.xhalf = x0 + self.xstep
        self.yhalf = y0 + self.ystep
        self.zhalf = z0 + self.zstep
        self.tx = x0, self.xhalf
        self.ty = y0, self.yhalf
        self.tz = z0, self.zhalf
        self.inter = []
        self.children = np.zeros((2, 2, 2), dtype=box)
        
        plt.plot([x0, x0, x0+xs, x0+xs], [y0, y0+ys, y0+ys, y0])
        if not self.empty and not self.stop:
            count += 1
            self.checksons()
            self.assign()
            
    def status(self):
        empties = 0
        stops = 0
        for item in self.points:
            l = len(item)
            if l == 0:
                empties += 1
            elif l == 1:
                stops += 1
        if empties == 3:
            self.empty = True
        if stops == 3:
            self.stop = True

    def assign(self):
        for i in range(8):
            j, k, l = i%2, int(i/2)%2, int(i/4)
            pos = self.inter[i]
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
            
        for i in range(8):
            j, k, l = i%2, int(i/2)%2, int(i/4)
            inter = np.intersect1d(np.intersect1d(pos[j], pos[k+2]), pos[l+3]).astype(int)
            self.inter.append(inter)

N = 10
count = 0
#x = np.array([randint(N*2) for i in range(N)])
#y = np.array([randint(N*2) for i in range(N)])
#z = np.array([randint(N*2) for i in range(N)])
x = 2 - random(N)
y = 2 - random(N)
z = 2 - random(N)
data = [x, y, z]
a = box(data, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
print(count)
plt.show()
#a.assign()
#plt.plot(x[[0, 1, 2, 3, 4]], y[[0, 1, 2, 3, 4]], "o")


