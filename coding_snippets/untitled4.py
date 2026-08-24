# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 07:48:03 2026

@author: Alexey
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

nP = 6
true_sig = 2
x0 = -0.5
x = np.linspace(-5, 5, nP)
p = np.exp( - (x-x0)**2/2/true_sig**2)

log_post = np.log(p)

ii, = np.unravel_index(np.argmax(log_post), log_post.shape)

mini = max(ii-1, 0) if ii==log_post.shape[0] else max(ii-2, 0)

cx = np.polyfit(x[mini:(mini+3)], log_post[mini:(mini+3)], 2)

# minj = max(jj-1, 0) if jj==log_post.shape[1] else max(jj-2, 0)

x0_fit = -cx[1] / (2 * cx[0])
sigmaX_fit = np.sqrt(1 / (2 * abs(cx[0])))


nInterp = 10

finex = np.linspace(x[0], x[-1], nP*nInterp)
spline = CubicSpline(x, log_post)
fine_vals = spline(finex)

fig = plt.figure(num='aaa')
fig.clf()
plt.plot(x, p, 'o')
plt.plot(finex, np.exp(fine_vals), '-b')
plt.plot(finex, np.exp( -(finex-x0_fit)**2/2/sigmaX_fit**2), ':k')

plt.show()