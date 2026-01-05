import sys
sys.path.append('../')
#import brukerread as Mybr
import numpy as np

import matplotlib
#matplotlib.use('WXAgg')                 # Force the WXAgg backend
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg 
import matplotlib.tri as mtri
from matplotlib import cm
import matplotlib.pyplot as plt

from SysPar import sysPar, expPar
from hyscore_sim import HYSCOREsim
from eigenfield_my import SpinOp, fibonacci

from scipy.linalg import eig, eigh

nKnots = 20
Sys = sysPar(S=[1/2], g=[np.array([2, 2.05, 2.1])], Nucs='15N', A=[np.array([50, 50, 100])])

mu_B = 9.274009994e-24          # J/T  (Bohr magneton)
hplank = 6.62607015e-34          # J·s

S = Sys.S[0]
I = Sys.I[0]
sx = SpinOp([S, I], 'xe')
sy = SpinOp([S, I], 'ye')
sz = SpinOp([S, I], 'ze')

ix = SpinOp([S, I], 'ex')
iy = SpinOp([S, I], 'ey')
iz = SpinOp([S, I], 'ez')

ndim = sz.shape[0]

omega=9.5*1e12*hplank/mu_B # GHz to mT
Ident = np.eye(sz.shape[0])
# phi = np.array([0.0])
# theta = np.array([1])
Tphi, Ttheta, Tww = fibonacci(nKnots)

HFA = Sys.A[0]*1e9*hplank/mu_B
gn = Sys.gn[0]

ZfHAM = HFA[0]*sx@ix+HFA[1]*sy@iy+HFA[2]*sz@iz

maA = np.kron(Ident, Ident)*omega - np.kron(ZfHAM, Ident) + np.kron(Ident, np.conj(ZfHAM))  # HF coupling goes here

g = Sys.g[0]
#### not sure if this is cosher... I did not think this too hard about:
BHAMx = np.kron(sx*g[0] + ix*gn, Ident) - np.kron(Ident, np.conj(sx*g[0] + ix*gn) )
BHAMy = np.kron(sy*g[1] + iy*gn, Ident) - np.kron(Ident, np.conj(sy*g[1] + iy*gn) )
BHAMz = np.kron(sz*g[2] + iz*gn, Ident) - np.kron(Ident, np.conj(sz*g[2] + iz*gn) )

Detx = -sz
Dety = -sz
Detz = +sy

AB = np.linalg.solve(maA, BHAMx)
eigvals, eigvecs = eigh(AB)
eigvals = np.real(1/eigvals)

for ii in range(ndim**2):
    vv = eigvecs[:, 0] 
    if eigvals[ii]>0 and eigvals[ii]<400:
        sexp = eigvecs[:,5]@sz.ravel()
        iexp = eigvecs[:,5]@iz.ravel()
        print(f'{eigvals[ii]} - {sexp} : {iexp}')
dd = np.diagonal(vv.reshape((ndim,ndim)))
