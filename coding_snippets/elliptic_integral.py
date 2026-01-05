# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 20:56:40 2025

@author: Alexey Silakov
"""

import numpy as np
import matplotlib.pyplot as plt

def amps(BB, Bx, By, Bz):
    pi = np.pi
    
    if BB <= By:
        S1 = 2.0/pi * Bx * By * Bz / BB**2 / np.sqrt((Bx**2 - By**2)*(BB**2 - Bz**2))
        k = (Bx**2 - BB**2)*(By**2 - Bz**2) / ((Bx**2 - By**2)*(BB**2 - Bz**2))
    else:
        S1 = 2.0/pi * Bx * By * Bz / BB**2 / np.sqrt((By**2 - Bz**2)*(Bx**2 - BB**2))
        k = (Bx**2 - By**2)*(BB**2 - Bz**2) / ((Bx**2 - BB**2)*(By**2 - Bz**2))
    
    # Initialize AGM (arithmetic-geometric mean) variables
    a0 = 1.0
    b0 = np.sqrt(1.0 - k)
    c0 = np.sqrt(k)
    
    an = a0
    bn = b0
    cn = c0
    
    cnt = 0
    KK = 0.0
    
    # AGM iteration
    while cn > 1e-15:
        an = (a0 + b0)/2.0
        bn = np.sqrt(a0*b0)
        cn = (a0 - b0)/2.0
        
        a0 = an
        b0 = bn
        c0 = cn
        
        cnt += 1
        if cnt > 1000:  # safety break
            break
    print(cnt)
    KK = pi / (2.0 * an)
    
    return S1 * KK


# Example parameters
Bx = 2
By = 2.1
Bz = 2.2

# R range: slightly above Bz to slightly below Bx
BB_values = np.linspace(Bx, Bz, 200)
amps_values = np.array([amps(BB, Bx, By, Bz) for BB in BB_values])

# Plot the result
plt.figure(figsize=(8,5))
plt.plot(BB_values, amps_values, label='amps(BB)')
plt.xlabel('BB')
plt.ylabel('amps')
plt.title('amps vs BB')
plt.grid(True)
plt.legend()
plt.show()
