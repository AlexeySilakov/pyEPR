import numpy as np
from scipy.optimize import curve_fit, minimize_scalar# -*- coding: utf-8 -*-
from scipy import special as sp
import matplotlib.pyplot as plt

def model(x, amp, x0, sigma, shape):
    
    d = shape/np.sqrt(1+shape**2)
    mz=d*np.sqrt(2/np.pi)
    sz = np.sqrt(1-mz**2)
    g1 = (4-np.pi)/2*mz**3/((1-mz**2)**(3/2))
    if shape==0:
        zm = mz-g1*sz/2
    else:
        zm = mz-g1*sz/2-np.sign(shape)/2*np.exp(-2*np.pi/np.abs(shape))
    z = (x - (x0-sigma*zm)) / (sigma)
    return amp * np.exp(-z**2 / 2) * (1 + sp.erf(shape * z / np.sqrt(2)))

def fitModel(x, yexp):
    amp0 = np.max(yexp)
    x0 = x[np.argmax(yexp)]
    sigma0 = np.mean(x*yexp)
    shape0 = 0.0
    prop, _ = curve_fit(model, x, yexp, p0=[amp0, x0, sigma0, shape0],
        bounds=([0.1,     x.min() - 5 * sigma0, x[1]-x[0], -3],
                [10*amp0, x.max() + 5 * sigma0,x[-1]-x[0],  3]),
        maxfev=400000)
    return prop


    
    
xe0= 50
sigma = 10
shape = 5
nPoints = 100
rand_Amp = 0.1

np.random.seed(0)


yrand =  rand_Amp*np.random.normal(0, 1.0, nPoints)

xe = np.linspace(1, 100, 100);
y0 = np.exp(-(xe-xe0)**2/2/sigma**2)
y1 = model(xe, 1, xe0, sigma, shape)

yexp = y1 + yrand
############## fit ##############

prop = fitModel(xe, yexp)
amp, x0, sigma, shape = prop
x = np.linspace(max(x0-sigma*10, 0), x0+sigma*10, 1000);
yfit = model(x, amp, x0, sigma, shape)

ww = np.sort(yfit)[::-1]
wi = np.argsort(yfit)[::-1]
res = np.cumsum(ww)
res/=np.max(res)
ii = np.argmin(np.abs(res-0.95))
int_width = float(ii)/(x[1]-x[0])

midI = np.argmin(np.abs(x-x0))
pt1 = np.argmin(np.abs(yfit[0:midI]-ww[ii]))
pt2 = np.argmin(np.abs(yfit[midI:-1]-ww[ii]))+midI

fig, ax = plt.subplots(1, 1, figsize=(13.5, 5), num='test', clear=True)

ax.plot(xe, y0/np.max(y0))
ax.plot(xe, yexp/np.max(y1))
ax.plot(x, yfit/np.max(y1))
ax.plot(x, res)
ax.plot(x[ii], res[ii], '*')
ax.plot(x, ww/np.max(y1), ':')
ax.axhline(y = ww[ii]/np.max(y1), linestyle='--', color = 'gray')
ax.plot(x[pt1], yfit[pt1]/np.max(y1), 'd')
ax.plot(x[pt2], yfit[pt2]/np.max(y1), 'd')
ax.set_title(f'amp={amp:0.3f}(1), x0={x0:0.3f}({xe0:0.3f})')
plt.show()