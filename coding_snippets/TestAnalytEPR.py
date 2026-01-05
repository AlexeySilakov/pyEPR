import numpy as np
import matplotlib.pyplot as plt

xField = np.linspace(300, 400, 1000)
PrincipleFields = np.array([310, 350, 370])

a, b, c = np.sort(PrincipleFields)

nPoints = 100
tt = np.linspace(0, np.pi/2, nPoints, endpoint=False)

R, t = np.meshgrid(xField, tt)
gamma = 1.0 - (R**2)/c**2
gamma[gamma<0]=0
alpha = np.ones_like(gamma)*(1.0/a**2 - 1.0/c**2)
beta  = np.ones_like(gamma)*(1.0/b**2 - 1.0/c**2)

X = np.sqrt(gamma / alpha)
Y = np.sqrt(gamma / beta)

xx = X * np.cos(t)
yy = Y * np.sin(t)
RR1 = R**2 - xx**2 - yy**2
xx[RR1<0]=0
yy[RR1<0]=0
RR1[RR1<0]=0
zz = np.sqrt(RR1)


dt = np.pi/2 / (nPoints-1)
x1 = X * np.cos(t+dt)
y1 = Y * np.sin(t+dt)
RR2 = R**2 - x1**2 - y1**2
x1[RR2<0]=0
y1[RR2<0]=0
RR2[RR2<0]=0
z1 = np.sqrt(RR2)

xx-=x1
yy-=y1
zz-=z1

RR3 = xx**2+yy**2+zz**2

Result =  np.sum( np.sqrt(RR3), axis=0)


plt.figure(1, figsize=(8, 4))
plt.plot(xField, Result)
plt.xlabel('Field')
plt.title('1‑D Plot Example')
plt.legend()
plt.grid(True)
plt.show()