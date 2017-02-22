import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random, randint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class box():
    def __init__(self, points, r0, rs, parent = False):
        global ax
        self.points = points
        self.r0 = r0
        self.r1 = [r0[i] + rs[i] for i in range(3)]
        self.steps = [0.5*rs[i] for i in range(3)]
        self.halfs = [r0[i] + self.steps[i] for i in range(3)]        
        self.tx = self.r0[0], self.halfs[0]
        self.ty = self.r0[1], self.halfs[1]
        self.tz = self.r0[2], self.halfs[2]
        
        self.numParticles = len(points[0])
        self.children = [None for i in range(8)]
        self.sons_points()
        if self.numParticles == 1:
            self.plot()
        if parent:
            self.plot(c = "r", w = 0.25)
        
    def plot(self, a = 1, c = "b", w = 0.5):        
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
            low = self.points[i] < half
            high = self.points[i] >= half
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
            if len(inter) > 0 and self.numParticles > 1:
                self.children[i] = box([x, y, z], [self.tx[j], self.ty[k],
                         self.tz[l]], self.steps) 

N = 20

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
ax.set_axis_off()
ax.plot(x, y, z, "o", c="r", ms=2)
a = box(data, [min_x, min_y, min_z], [max_x-min_x, max_y-min_y, max_z-min_z], parent=True)
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")

def update(i):
    ax.view_init(30, i*360/frames)
    
frames = 50
ani = FuncAnimation(fig, update, frames=frames, interval=2*frames)
ani.save("concept.gif", writer='imagemagick')
#plt.show()