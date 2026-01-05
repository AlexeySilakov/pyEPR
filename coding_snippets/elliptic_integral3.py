# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:02:31 2025

@author: Alexey Silakov
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipe, elliprf, elliprj

Bx, By, Bz = 310, 350, 375

nPoints= 1000
tF = np.linspace(By+0.001, Bz-0.0001, 1000)
tt = np.linspace(0, np.pi/2, nPoints)
gw = 2
fga = np.fft.ifft(np.fft.fftshift(np.exp( -(tF-np.mean(tF))**2/gw**2 )))

dt = np.max(tt)/(nPoints-1)
F, t = np.meshgrid(tF, tt)

A = Bx/F*np.sqrt((Bz**2-F**2)/(Bz**2-Bx**2))
B = By/F*np.sqrt((Bz**2-F**2)/(Bz**2-By**2))

x = A*np.cos(t)
y = B*np.sin(t)
z = np.sqrt(1-x**2-y**2)

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
zz = np.sum(np.sqrt(dz**2), axis=0)

#### other solution ####

zeta2 = (A**2-B**2)/(1-A**2)*np.sin(t)**2
ds = np.sqrt(B**2*(1+zeta2/B**2)/(1+zeta2) )

amps2 = np.sum(ds*dt, axis=0)

# fa = np.fft.ifft(np.fft.fftshift(amps2))


# tA = Bx/tF*np.sqrt((Bz**2-tF**2)/(Bz**2-Bx**2))
# tB = By/tF*np.sqrt((Bz**2-tF**2)/(Bz**2-By**2))
# D = (tA**2-tB**2)/(1-tA**2)
# k2 = D/((tB**2+D)*(1+D))
# n = 1-tB**2
# E_val = ellipe(k2) 
# Pi_val = elliprf(0,1-k2,1)-n/2*elliprj(0,1-k2,1,1+n)

# amps3 = 4*tB*(E_val/tB-n*Pi_val)

plt.plot(tF[1:], np.diff(amps)/np.diff(zz),'r')
# plt.plot(tF, amps2,'b')
# plt.plot(tF, np.fft.fftshift(np.fft.fft(fa*fga)),'m')
#plt.plot(tF, amps3,'g')
plt.show()