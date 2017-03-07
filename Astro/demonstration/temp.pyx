#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:45:40 2017

@author: juan
"""

import numpy as np
#import matplotlib.pyplot as plt
cimport numpy as np
import copy

#points 3xn array

cdef class box:
    cdef public np.ndarray points, lower_bound, coordinate_size, bound_data 
    cdef public np.ndarray box_half_size, upper_bound, box_half, tmatrix, center_mass
    cdef public int number_of_points
    cdef public list children, tree
    cdef public bint parent
    cdef double length, mass, tolerance
    cdef Py_ssize_t i, j, k, l
    
    def __init__(self, points, lb, cs, tr, root = False,
                                     parent = False, mass = 1, tolerance = 1):
        cdef np.ndarray inter
        
        self.points = points
        self.parent = parent
        self.coordinate_size = cs
        self.lower_bound = lb
        self.number_of_points = points.shape[1]
        self.box_half_size = 0.5*self.coordinate_size
        self.upper_bound = self.lower_bound + self.coordinate_size
        self.box_half = self.lower_bound + self.box_half_size
        self.bound_data = np.zeros((6, self.number_of_points))
        self.mass = mass*self.number_of_points
        self.tolerance = tolerance
        self.children = []
        self.tree = tr
        if self.number_of_points > 0:
            self.tree.append(self)
        self.tmatrix = np.zeros((3, 2))
        self.tmatrix[0] = self.lower_bound[0], self.box_half[0]
        self.tmatrix[1] = self.lower_bound[1], self.box_half[1]
        self.tmatrix[2] = self.lower_bound[2], self.box_half[2]
        
        self.center_mass = np.sum(self.points, axis=-1)/self.mass
        
        self.length = np.sqrt(np.sum(self.coordinate_size**2))
        
        if self.number_of_points > 1:
            for i in range(3):
                half = self.box_half[i]
                low = self.points[i] <= half
                high = self.points[i] > half
                self.bound_data[2*i] = low
                self.bound_data[2*i+1] = high
                
            for i in range(8):
                j = i%2
                k = int(i/2)%2
                l = int(i/4)
                inter = (self.bound_data[j]==1) & (self.bound_data[k+2]==1) & (self.bound_data[l + 4]==1)
                num = sum(inter)
                points = self.points[:, inter]
                if num > 0:
                    if num == 1:
                        parent = False
                    else:
                        parent = True
                    son_box = box(points, np.array([self.tmatrix[0, j], self.tmatrix[1, k],
                             self.tmatrix[2, l]]), 0.5*self.coordinate_size, self.tree, parent=parent)
    
                    self.children.append(son_box)
                    
    cpdef force(self, point):
        cdef double l, d, theta
        cdef box child
        cdef tuple ans
        cdef list m, cm
        cdef Py_ssize_t i
        l = self.length
        if self.center_mass[0] != point[0]:
            r = self.center_mass - point
            d = np.sqrt(sum(r**2) + 0.05)
            theta = l/d
            if theta < self.tolerance:
                a = np.zeros(3)
                for child in self.children:
                    a += child.force(point)
                return a
            else:
                return r*self.mass/(d**3)
        else:
            return np.zeros(3)

    def plot(self, ax, a = 1, c = "b", w = 0.5, child_only = True):
        if child_only:
            if not self.parent:
                xy, yx = np.meshgrid([self.lower_bound[0], self.upper_bound[0]], [self.lower_bound[1], self.upper_bound[1]])
                xz, zx = np.meshgrid([self.lower_bound[0], self.upper_bound[0]], [self.lower_bound[2], self.upper_bound[2]])
                yz, zy = np.meshgrid([self.lower_bound[1], self.upper_bound[1]], [self.lower_bound[2], self.upper_bound[2]])
                self.surface(xy, yx, self.upper_bound[2], a, c, w, ax)
                self.surface(xy, yx, self.lower_bound[2], a, c, w, ax)
                self.surface(xz, self.lower_bound[1], zx, a, c, w, ax)
                self.surface(xz, self.upper_bound[1], zx, a, c, w, ax)
                self.surface(self.upper_bound[0], yz, zy, a, c, w, ax)
                self.surface(self.lower_bound[0], yz, zy, a, c, w, ax)
        else:
            xy, yx = np.meshgrid([self.lower_bound[0], self.upper_bound[0]], [self.lower_bound[1], self.upper_bound[1]])
            xz, zx = np.meshgrid([self.lower_bound[0], self.upper_bound[0]], [self.lower_bound[2], self.upper_bound[2]])
            yz, zy = np.meshgrid([self.lower_bound[1], self.upper_bound[1]], [self.lower_bound[2], self.upper_bound[2]])
            self.surface(xy, yx, self.upper_bound[2], a, c, w, ax)
            self.surface(xy, yx, self.lower_bound[2], a, c, w, ax)
            self.surface(xz, self.lower_bound[1], zx, a, c, w, ax)
            self.surface(xz, self.upper_bound[1], zx, a, c, w, ax)
            self.surface(self.upper_bound[0], yz, zy, a, c, w, ax)
            self.surface(self.lower_bound[0], yz, zy, a, c, w, ax)
            
            
    def surface(self, x, y, z, a, c, w, ax):
        ax.plot_wireframe(x, y, z, alpha=a, color=c, linewidth = w)
    
cpdef limits(matrix):
    cdef Py_ssize_t i
    cdef list x_min, x_max, rank
    x_min = [min(matrix[i])*0.9 if min(matrix[i]) > 0 else min(matrix[i])*1.1 for i in range(3)]
    x_max = [max(matrix[i]) for i in range(3)]
    rank = [1.1*(x_max[i] - x_min[i]) for i in range(3)]
    return np.array(x_min), np.array(rank)

        
def solver(np.ndarray[np.float64_t, ndim=2] positions, np.ndarray[np.float64_t, ndim=2] speeds,
           int N, double t, double dt , double G, filename='Data/'):
    cdef int n
    cdef np.ndarray[np.float64_t, ndim=2] x
    cdef np.ndarray[np.float64_t, ndim=2] v
    cdef np.ndarray[np.float64_t, ndim=1] x0, r0
    cdef box tree
    n = int(t/dt)
    x = positions
    v = speeds
    accelerations = np.zeros((3, N))
    np.savetxt("%s%d_instant.dat"%(filename, 0), x.T) 
    pos = []
    trees = []
    for i in range(n-1):
        print(i)
        x0, r0 = limits(x)
        tree = box(x, x0, r0, [], parent=True, root=True)
        for j in range(N):
            accelerations[:, j] = tree.force(x[:,j])
        v_half = v + 0.5*dt*accelerations
        x += dt*v_half
        x0, r0 = limits(x)
        tree = box(x, x0, r0, [], parent=True, root=True)
        for j in range(N):
            accelerations[:, j] = tree.force(x[:,j])
        v = v_half + 0.5*dt*accelerations
        np.savetxt("%s%d_instant.dat"%(filename, i+1), x.T)        
    return n
