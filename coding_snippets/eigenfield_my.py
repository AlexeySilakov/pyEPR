import numpy as np
from scipy.linalg import eig, eigh
import matplotlib.pyplot as plt
import time

def SpinOp(spins, corrs):
    """
    Build the tensor‑product operator

        totalout = kron( S( s1 , corr1 ), S( s2 , corr2 ), … )

    Parameters
    ----------
    spins : sequence of floats or ints
        The total spin of each subsystem, e.g. [1/2, 1, 3/2, ...].
    corrs : sequence of str or a single str
        For each spin one of the characters 'x', 'y', 'z' or 'e'.
        'e' denotes the identity operator.

    Returns
    -------
    totalout : ndarray, shape (∏(2*s_i+1), ∏(2*s_i+1))
        The Kronecker product of all local operators.
    """
    # ------------------------------------------------------------------
    # Start with the 1×1 identity – this is the neutral element for
    # Kronecker products.
    totalout = np.array([[1.0]], dtype=complex)

    # ------------------------------------------------------------------
    # Iterate over all subsystems
    for s, corr in zip(spins, corrs):
        # dimension of the local Hilbert space
        dim = int(2 * s + 1)

        if corr == 'x':
            # ----- Sx -------------------------------------------------
            sp = np.zeros((dim, dim), dtype=complex)
            sm = np.zeros((dim, dim), dtype=complex)

            for ii in range(dim - 1):
                jj = ii + 1
                ms = -s + ii
                val = np.sqrt(s * (s + 1) - ms * (ms + 1))
                sp[ii, jj] = val

            for ii in range(1, dim):
                jj = ii - 1
                ms = -s + ii
                val = np.sqrt(s * (s + 1) - ms * (ms - 1))
                sm[ii, jj] = val

            S_s = (sp + sm) / 2.0

        elif corr == 'y':
            # ----- Sy -------------------------------------------------
            sp = np.zeros((dim, dim), dtype=complex)
            sm = np.zeros((dim, dim), dtype=complex)

            for ii in range(dim - 1):
                jj = ii + 1
                ms = -s + ii
                val = np.sqrt(s * (s + 1) - ms * (ms + 1))
                sp[ii, jj] = val

            for ii in range(1, dim):
                jj = ii - 1
                ms = -s + ii
                val = np.sqrt(s * (s + 1) - ms * (ms - 1))
                sm[ii, jj] = val

            S_s = (sp - sm) / (2j)

        elif corr == 'z':
            # ----- Sz -------------------------------------------------
            # diagonal entries are m = -s, -s+1, … , +s
            ms_vals = np.arange(-s, s + 1)
            S_s = np.diag(ms_vals)

        elif corr == 'e':
            # ----- Identity --------------------------------------------
            S_s = np.eye(dim, dtype=complex)

        else:
            raise ValueError(f"Unknown correlation type '{corr}'. "
                             "Must be one of 'x', 'y', 'z', or 'e'.")

        # Kronecker product with the current operator
        totalout = np.kron(totalout, S_s)

    return totalout
def fibonacci(nKnots, Symmetry=None):
    epsilon=0.33
    n   = int(nKnots**1.5)
    thLim=np.pi
    if not Symmetry:
        phLim=2*np.pi
    # Golden angle
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))   # ≈ 2.39996322973
    idx = np.arange(0, n, dtype='float64')
    cos_theta = (idx + epsilon)/ (float(n - 1)+2*epsilon)
    ttheta = np.acos(cos_theta)   # polar angle from +z axis
    tphi = np.remainder(golden_angle*idx, 2*np.pi)
    mask = (ttheta <= thLim) & (tphi <= phLim)

    phi = tphi[mask]
    theta = ttheta[mask] 
    ww = np.ones_like(phi)
    return phi, theta, ww


# if __name__ == "__main__":

nPoints = 1000
nKnots = 100
omega= 9.5 # GHz
lw = 0.1 # in mT
g = np.array([1.99, 2, 2.01])
A = np.array([500., 500., 0]) # MHz
S = 1/2
I = 1/2
# --------------- start ----------------
mu_B = 9.274009994e-24          # J/T  (Bohr magneton)
hplank = 6.62607015e-34          # J·s

sx = SpinOp([S, I], 'xe')
sy = SpinOp([S, I], 'ye')
sz = SpinOp([S, I], 'ze')

ix = SpinOp([S, I], 'ex')
iy = SpinOp([S, I], 'ey')
iz = SpinOp([S, I], 'ez')

A*=1e9*hplank/mu_B # to mT
omega*=1e12*hplank/mu_B # to mT
# phi = np.array([0.0])
# theta = np.array([1])
Tphi, Ttheta, Tww = fibonacci(nKnots)

Ident = np.eye(sz.shape[0])

AHAM = sx@ix*A[0]+A[1]*sy@iy+A[2]*sz@iz

maA = np.kron(Ident, Ident)*omega - np.kron(AHAM, Ident) + np.kron(Ident, np.conj(AHAM))  # HF coupling goes here

Xaxis = np.linspace(300, 400, nPoints)
dx = Xaxis[1]-Xaxis[0]

Spectrum = np.zeros_like(Xaxis)
threshold = 1e-3
timest = time.time()

for ii in range(len(Tphi)):
    lx = np.sin(Ttheta[ii])*np.cos(Tphi[ii])
    ly = np.sin(Ttheta[ii])*np.sin(Tphi[ii])
    lz = np.cos(Ttheta[ii])

    lxd = -np.cos(Ttheta[ii])*np.sin(Tphi[ii])
    lyd = +np.cos(Ttheta[ii])*np.cos(Tphi[ii])
    lzd = -np.sin(Ttheta[ii])
    
    ZHAM = sx*lx*g[0]+sy*ly*g[1]+sz*lz*g[2]
    DetHAM = sx*lxd+sy*lyd+sz*lzd
    # ZHAM*=mu_B/hbar*1e-9

    B = np.kron(ZHAM, Ident) - np.kron(Ident, np.conj(ZHAM))
    # eigvals, eigvecs = eig(maA, B)
    detVec = DetHAM.ravel()
     
    # AB = maA@np.linalg.inv(B)
    AB = np.linalg.solve(maA, B)
    eigvals, eigvecs = eigh(AB)
    eigvals = np.real(1/eigvals)
    
    ###  Performance test. (1) runs as a C-code, so, 
    ### it it is faster than a loop in python below in (2) even though it does more stuff
    # (1)
    C = np.tile(detVec, (B.shape[0], 1))
    prob = np.abs(np.diagonal(C @ eigvecs))**2 
    # (2)
    # prob = np.zeros_like(eigvals)
    # for nn in range(len(prob)):
    #     prob[nn]= np.abs(detVec@eigvecs[:, nn])**2
        
    idx = np.floor((eigvals-Xaxis[0])/dx).astype(np.int64)  
    mask = (idx >= 0) & (idx <= nPoints) & (prob>threshold)
    np.add.at(Spectrum, idx[mask], prob[mask]) 

    # # unoptimized code for now, just to make sure math is right
    # for jj in range(len(eigvals)):
    #     prob = np.abs(detVec@eigvecs[:, jj])**2
    #     if prob>threshold and eigvals[jj].real<np.inf:
    #         # print(f'{np.abs(eigvals[jj])}, {prob}')
    #         B = eigvals[jj].real
    #         idx = int(np.floor((B-Xaxis[0])/dx))
    #         if idx>=0 and idx<nPoints:
    #             Spectrum[idx]+=prob*Tww[ii]/B
         
print(f'{time.time()-timest}')
sig = lw/2.354820045
GA = np.exp(-0.5*((Xaxis-np.mean(Xaxis))/sig)**2)

ispec = np.fft.ifft(np.fft.ifftshift(Spectrum))
iGA = np.fft.ifft(np.fft.ifftshift(GA))

LWSpectrum = np.abs(np.fft.fftshift(np.fft.fft(ispec*iGA)))

plt.figure(figsize=(8, 5))
plt.plot(Xaxis, LWSpectrum)
# print(eigvals)
# print(eigvecs)
plt.show()