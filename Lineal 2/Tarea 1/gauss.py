# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

clase = 11

def exchanger(num, denom):
    for x in range(1, clase):
        val = (num - denom*x)%clase
        if val == 0:
            break
    return x
    
def solver(numerator, denom):
    numerator = numerator%clase
    denom = denom%clase
    try:
        n = len(numerator)
    except:
        n = 0
    
    if n != 0:
        ans = np.zeros_like(numerator)
        for i in range(n):
            ans[i] = exchanger(numerator[i], denom)
    else:
        ans = exchanger(numerator, denom)
        
    return ans

def eliminacion_gaussiana(A, b):    
    A = np.float_(A)
    b = np.float_(b)
    n = len(b)
    for i in range(n):
        #unos en la diagonal
        a = A[i,i]
        A[i,:] = solver(A[i,:], a)
        b[i] = solver(b[i], a)

        # resta bajando
        for ii in range(i+1,n):
            a = A[ii,i]
            A[ii,:] = A[ii,:] - ((a*A[i,:])%clase)
            b[ii] = b[ii] - ((a*b[i])%clase)

    # reemplaza hacia arriba
    for i in range(n-1,-1,-1):
        for ii in range(i+1,n):
            b[i] = b[i] - ((A[i,ii]*b[ii])%clase)
    return b

A1 = np.array([[3, 4, -5],
              [-2, 3, 9],
              [1, 2, -3]])

A2 = np.array([[3, 4, -5],
              [-2, 3, 9],
              [2, 6, 0]])
    
B = np.array([1, 0, 1])

x1 = eliminacion_gaussiana(A1, B)
x2 = eliminacion_gaussiana(A2, B)

A = np.linalg.solve(A1, B)
print(A)
print("A1 = %s, A2 = %s"%(x1, x2))
