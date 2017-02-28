#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 19:30:32 2017

@author: juan
"""

import numpy as np
cimport numpy as np
from cython.parallel cimport parallel, prange
from libc.stdlib cimport abort, malloc, free
cimport openmp

cdef class box:
    cdef public np.ndarray points, masses, r0, r1, steps
    cdef int tolerance, numParticles
    cdef public bint parent, root
    cdef tuple tx, ty, tz
    cdef double mass, length
    cdef list halfs, center_mass, children
    cdef box tree
    def __init__(self, np.ndarray[np.float64_t, ndim=2] ppoints, 
                 np.ndarray[np.float64_t, ndim=1] pr0, np.ndarray[np.float64_t, ndim=1] prs,
                           np.ndarray[np.float64_t, ndim=1] pmasses, root = False,
                                     parent = False, tree = None):
        cdef Py_ssize_t i                       
        
        self.points = ppoints #array
        self.masses = pmasses #array
        self.parent = parent #boolean
        self.root = root #boolean
        self.tolerance = 1
        self.r0 = pr0 #array
        self.tree = tree
#        self.r1 = np.sum(pr0, prs, axis=0)
        self.r1 = np.array([pr0[i] + prs[i] for i in range(3)])
        self.steps = 0.5*prs
        self.halfs = [pr0[i] + self.steps[i] for i in range(3)]
        self.length = np.sum(prs**2)**0.5
        self.tx = self.r0[0], self.halfs[0]
        self.ty = self.r0[1], self.halfs[1]
        self.tz = self.r0[2], self.halfs[2]
        
        self.mass = sum(self.masses)
        self.center_mass = self.center_mass_calc()
        
        self.numParticles = self.points.shape[1]#len(ppoints[0])
        self.children = []
        self.sons_points()

    def force(self, np.ndarray[np.float64_t, ndim=1] point):
        cdef double l, d, theta
        cdef box child
        cdef tuple ans
        cdef list m, cm
        cdef Py_ssize_t i
        l = self.length
        if self.center_mass[0] != point[0]:
            d = sum([(self.center_mass[i] - point[i])**2 for i in range(3)])**0.5
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
                    return list(self.center_mass - point), np.array([self.mass])
            
    cdef center_mass_calc(self):
        return [sum(self.masses*self.points[i])/self.mass for i in range(3)]            
        
    cdef sons_points(self):
        cdef box son_box
        cdef double half
        cdef Py_ssize_t i, j, k, l
        cdef int num
        cdef bint parent
        cdef np.ndarray[np.float64_t, ndim=2] points
#        cdef np.ndarray[np.bool_t, ndim=1] low, high
        cdef list data
        
        data = []
        for i in range(3):
            half = self.halfs[i]            
            low = self.points[i] <= half# & (self.points[i] >= self.r0[i])
            high = self.points[i] > half# & (self.points[i] <= self.r1[i])
            data.append(low)
            data.append(high)
            
        for i in range(8):
            j = i%2
            k = int(i/2)%2
            l = int(i/4)
            inter = data[j] & data[k + 2] & data[l + 4]
            num = sum(inter)  
            points = self.points[:, inter]
            if num > 0 and self.numParticles > 1:
                if num == 1:
                    parent = False
                else:
                    parent = True
                son_box = box(points, np.array([self.tx[j], self.ty[k],
                         self.tz[l]]), self.steps, self.masses[inter], 
                        parent=parent, tree = self.tree)

                self.children.append(son_box)
                
cpdef limits(matrix):
    cdef Py_ssize_t i
    cdef list x_min, x_max, rank
    x_min = [min(matrix[:,i])*0.9 if min(matrix[:,i]) > 0 else min(matrix[:,i])*1.1 for i in range(3)]
    x_max = [max(matrix[i, :]) for i in range(3)]
#    x_max = [1.1*max(matrix[i, :]) - 1.1*x_min[i] for i in range(3)]
    rank = [1.1*(x_max[i] - x_min[i]) for i in range(3)]
    return np.array(x_min), np.array(rank)

#cpdef limits(matrix):
#    x_maxs = np.zeros(3)
#    x_mins = np.zeros(3)
#    for i in range(3):
#        x_min = min(matrix[:, i])
#        x_maxs[i] = max(matrix[:, i])
#        
#        if x_min > 0:
#            x_min *= 0.9
#        else:
#            x_min *= 1.1
#        if x_max > 0:
#            x_max *= 1.1
#        else:
#            x_max *= 0.9
#        x_mins[i] = x_min
#        x_maxs[i] = x_max
#    
#    rank = (x_maxs - x_mins)*1.5
#    return x_mins, rank


cdef magnitude(pos):
    cdef double r
    r = np.sqrt(np.dot(pos, pos))
    if r < 1000:
        r = 1000
    return r**3

cdef iterator(np.ndarray[np.float64_t, ndim=1] particle, box tree, double G):
    cdef Py_ssize_t i
    cdef tuple ans
#    cdef np.ndarray[np.float64_t, ndim=2] distance
#    cdef np.ndarray[np.float64_t, ndim=2] distance
    cdef list distance
    cdef np.ndarray[np.float64_t, ndim=1] masses
    cdef np.ndarray[np.float64_t, ndim=2] acceleration_array
    cdef np.ndarray[np.float64_t, ndim=1] acceleration_3d
    cdef int n
    ans = tree.force(particle)
    if ans is not None:
        distance = ans[0]
        masses = ans[1]
        n = len(masses)
        range_n = range(n)
        if n == 1:
            mag = magnitude(distance)
            ratio = masses/mag
            acceleration = G*ratio*distance
        else:
            acceleration_array = np.zeros((n, 3))
            for i in range(n):
                acceleration_array[i] = masses[i]*distance[i]/magnitude(distance[i]) 
#            mag_array = np.array(list(map(magnitude, distance)))
#            ratio = masses/mag_array
#            acceleration_array = np.array([ratio[i]*np.array(distance[i]) for i in range_n])
            acceleration_3d = G*np.sum(acceleration_array, axis = 0)
        return acceleration_3d
    else:
        return [0,0,0]
    
#cdef parrallel_acceleration(np.ndarray[np.float64_t, ndim=2] x, box tree, double G):
#    cdef int n, numthreads
#    cdef Py_ssize_t i
#    cdef np.ndarray[np.float64_t, ndim=2] temp
#    cdef double position[3]
#    n = x.shape[0]
#    temp = np.zeros((n, 3))
#    with nogil, parallel(num_threads=4):
#        numthreads = openmp.omp_get_num_threads()
#        for i in prange(n, schedule="dynamic"):
#            with gil:            
#                temp[i] = iterator(x[i], tree, G)
#    return temp
    
def solver(np.ndarray[np.float64_t, ndim=2] positions, np.ndarray[np.float64_t, ndim=2] speeds,
           np.ndarray[np.float64_t, ndim=1] masses, int N, double t,
                     double dt , double G, filename='Data/'):
    cdef int n
    cdef np.ndarray[np.float64_t, ndim=2] x
    cdef np.ndarray[np.float64_t, ndim=2] v
    cdef np.ndarray[np.float64_t, ndim=1] x0, r0
    cdef box tree
    n = int(t/dt)
#    positions_in_time = np.zeros((n, N+2, 3))
    x = positions.T#positions.T.copy()
    v = speeds.T#speeds.T.copy()
    
#    positions_in_time[0] = x
        
    for i in range(n-1):
        print(i)
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)        
        accelerations = np.array(list(map(lambda pos: iterator(pos, tree, G), x)))
        v_half = v + 0.5*dt*accelerations
        x += dt*v_half
        x0, r0 = limits(x)
        tree = box(x.T, x0, r0, masses, parent=True, root=True)
        accelerations = np.array(list(map(lambda pos: iterator(pos, tree, G), x)))#parrallel_acceleration(x, tree, G)
        v = v_half + 0.5*dt*accelerations
#        positions_in_time[i+1] = x
        np.savetxt("%s%d_instant.dat"%(filename, i+1), x)
        
    return n    