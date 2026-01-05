import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk


gx, gy,gz= 2.0, 2.04, 2.1

mwFreq = 9.5
Bx, By,Bz= mwFreq*71.448/gx,  mwFreq*71.448/gy,  mwFreq*71.448/gz
gx2, gy2,gz2= gx**2, gy**2, gz**2

nPoints= 1000

tF = np.linspace(Bx+0.001, Bz-0.0001, 1000)

g2 = (mwFreq*71.448/tF)**2

r1 = g2-gx2
r2 = g2-gy2
d32= gz2-gy2
d31= gz2-gx2
k2=(d32*r1-r2*d31)/(r1*d32)

sig = 1.0/np.sqrt(r1*d32)*ellipk(k2)

plt.plot(tF, sig,'r')
plt.show()
