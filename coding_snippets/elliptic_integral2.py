# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:02:31 2025

@author: Alexey Silakov
"""
import numpy as np
import matplotlib.pyplot as plt
a, b, c = 1.9, 2.0, 2.1

nPoints= 100
tF = np.linspace(2.00001, 2.1, 1000)
tt = np.linspace(0, 2*np.pi, nPoints)

F, t = np.meshgrid(tF, tt)

A = np.sqrt(1/a**2 - 1/c**2)
B = np.sqrt(1/b**2 - 1/c**2)
C = np.sqrt(1 - F**2/c**2)

x = C/A*np.cos(t)
y = C/B*np.sin(t)
z = np.sqrt(F**2-x**2-y**2)

# A = np.sqrt(F**2/a**2-1)
# B = np.sqrt(1/a**2 - 1/b**2)
# C = np.sqrt(1/a**2 - 1/c**2)

# y = A/B*np.cos(t)
# z =- A/C*np.sin(t)
# x = np.sqrt(F**2-z**2-y**2)

dx = np.diff(x, axis=0)
dy = np.diff(y, axis=0)
dz = np.diff(z, axis=0)
amps = np.sum(np.sqrt(dx**2+dy**2+dz**2), axis=0)


#### other solution ####
alpha = a**2-F**2
beta  = b**2-F**2
gamma = c**2-F**2
tA2 = -gamma/alpha; tB2 = -gamma/beta
M = 1/2.0*(tA2+tB2)-tA2*tB2
N = 1-1/2.0*(tA2+tB2)
ds = np.sqrt(np.abs((M-np.cos(2*t))/(N-np.cos(2*t))))

amps2 = np.sum(ds, axis=0)


plt.plot(71.448*9.5/tF, amps,'r')
plt.plot(71.448*9.5/tF, amps2,'b')
plt.show()