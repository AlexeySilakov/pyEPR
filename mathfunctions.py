# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 14:31:25 2025

@author: Alexey Silakov
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def RotMatrix(angles):
    ###  perform no checks for speed
    # Euler transformation matrix
    # Precalculate trigonometric functions of angles
    sa = np.sin(angles[0]);
    ca = np.cos(angles[0]);
    sb = np.sin(angles[1]);
    cb = np.cos(angles[1]);
    sg = np.sin(angles[2]);
    cg = np.cos(angles[2]);

   # Compute passive rotation matrix
    R = np.array([[cg*cb*ca-sg*sa,   cg*cb*sa+sg*ca,  -cg*sb],
                  [-sg*cb*ca-cg*sa,  -sg*cb*sa+cg*ca,   sg*sb],
                  [sb*ca,            sb*sa,            cb]])
    return R
def RotTransDiag(A, angles):
    ###  perform no checks for speed
    R=RotMatrix(angles)
    AA = np.diag(A)
    return R@AA@R.transpose()

def RotTrans3x3(A, angles):
    ###  perform no checks for speed
    R=RotMatrix(angles)
    return R@A@R.transpose()

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


def adaptSim(Data_x, Data_y, Sim, SimX, SimY):
    """Resample a simulated spectrum onto an experimental region's axes."""
    if Sim is None:
        return None
    interp = RegularGridInterpolator((SimX, SimY), Sim,
                                      bounds_error=False, fill_value=0.0)
    gx, gy = np.meshgrid(Data_x, Data_y, indexing='ij')
    return interp(np.stack([gx.ravel(), gy.ravel()], axis=-1)).reshape(gx.shape)

def calcSSR(Exp, Sim, Mask, offset):
    """Sum of squared residuals between an experimental region and a
    simulation scaled to best fit it. Exp and Sim must share the same XY."""
    denom = np.sum(Sim[Mask] ** 2)
    if denom > 1e-30:
        c = max(np.sum((Exp[Mask]-offset) * Sim[Mask]) / denom, 0.0)
    else:
        c = 0.0
    scaled_sim = c * Sim + offset
    resid = Exp[Mask] - scaled_sim[Mask]
    SumSqRes = float(np.sum(resid ** 2))
    N = int(resid.size)

    return SumSqRes, N, scaled_sim

if __name__ == "__main__":
    RR=RotMatrix(np.array([np.pi/3, np.pi/2, np.pi]))
    print(RR)
    
    AA=RotTransDiag(np.array([1., 2., 3.]), np.array([np.pi/3, np.pi/2, np.pi]))
    print(AA)
    
    BB=RotTrans3x3(AA, np.array([-np.pi/3, -np.pi/2, -np.pi]))
    print(BB)
    print(RR.transpose()@AA@RR)
    
    SX1 = SpinOp([1/2, 1/2], ['x','e'])
    print(SX1)
    
    SZ2 = SpinOp([1/2, 1/2], ['e','z'])
    print(SZ2)