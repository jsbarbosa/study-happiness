#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:03:14 2017

@author: juan
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib.animation import FuncAnimation

N = 10000

costheta = (np.random.random(N) - 0.5)*2

theta = np.arccos(costheta) 
phi = np.random.random(N)*2*np.pi
                      
xy = np.zeros((2, N))
xy[0] = np.sin(theta)*np.cos(phi)
xy[1] = np.sin(theta)*np.sin(phi)
z = np.cos(theta)

fig = plt.figure()
ax = p3.Axes3D(fig)

max_R = 10

ax.set_xlim(-max_R, max_R)
ax.set_ylim(-max_R, max_R)
ax.set_zlim(-max_R, max_R)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')

#ax.set_axis_off()

points = ax.plot([],[],[], "o", ms=0.5)[0]
                      
def update(i):
    j = i/max_R
    points.set_data(j*xy)
    points.set_3d_properties(j*z)
    ax.view_init(30, i/100*360)
    
ani = FuncAnimation(fig, update, 100, interval = 50,  blit=False)
#ani.save("Increasing_sphere.gif", writer = "imagemagick")

plt.show()