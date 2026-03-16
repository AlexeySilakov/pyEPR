import numpy as np

#from eprcsvread import csvread
#import mathfunctions as mf

from SysPar import sysPar, expPar

from scipy.linalg import eigh #, eig,ishermitian, pinvh 
from scipy.fft import fft, ifft, fftshift
from scipy.special import ellipk

import matplotlib.pyplot as plt

### using optimization hints from ChatGPT:
# time no optimization: 0.9194717407226562
# 9)- no difference ... "silent killer" my heinie
# 1) - slight improvement for the provided suggestion, but nothing to write home about: 0.9081077575683594
# 1+4) - 0.887643575668335
# 1+4+2) - time: 0.4007081985473633 @ using option 1
# 4B) BS
# 5) makes no difference
# 6) actually makes it slower
# 7) no difference
# see option 0 in the run function - allows to change from eig to eigh. Works great so far. time: 0.2119410037994384
    # 1) Use eigh instead of eig whenever possibl
    #     eigvals, eigvecs = eigh(self.BHAM[ii], self.AHAM, check_finite=False)
        
    # 2) Pre-diagonalize AHAM once (very important)
    #     Do once:
    #     evalsA, U = eigh(self.AHAM)
    #     Ainv_sqrt = U @ np.diag(1/np.sqrt(evalsA)) @ U.conj().T
    #     Then:
    #     Btilde = Ainv_sqrt @ self.BHAM[ii] @ Ainv_sqrt
    #     eigvals, eigvecs = eigh(Btilde)
    #     eigvecs = Ainv_sqrt @ eigvecs
    #     Removes repeated factorization of AHAM
    #     ✔ Often 2–4× faster overall
    
    # 3) Massive Kronecker products (🔥 memory + CPU)
    #     Where it hurts
    #     self.BHAM[ii] = np.kron(self.Zee[ii], Ident) - np.kron(Ident, self.Zee[ii].T)
    #     self.AHAM = np.kron(Ident, Ident)*self.Omega - ...
    #     Use scipy.sparse for Liouville space
    #     from scipy.sparse import kron, eye, csc_matrix
    #     Ident = eye(n, format='csc')
    #     kron(..., format='csc')
    #     from scipy.sparse.linalg import eigsh
    #     (B) Avoid building BHAM explicitly if possible
    #     You only ever apply BHAM @ v.
    #     from scipy.sparse.linalg import LinearOperator
    #     def B_action(v):
    #         V = v.reshape((n,n))
    #         return (Z @ V - V @ Z.T).ravel()
    #     Bop = LinearOperator((n*n, n*n), matvec=B_action)
    #     This avoids allocating an n*n, n*n matrix entirely.
    
    # 4) findTransition() diagonalizes Hamiltonian repeatedly (🔥)
    #     A Cache eigenpairs per (LL, field)
    #    
    #     B Avoid full eigen-decomposition
    #     Instead:
    #         idx = np.argmin(np.abs(d[:,None]-d[None,:]-Omega))
    #     Replace with:
    #         dd = d[:,None] - d
    #         mask = np.triu(np.ones_like(dd), 1)
    #         idx = np.argmin(np.abs((dd - Omega)*mask + 1e9*(1-mask)))
    
    # 5) FFT pipeline inefficiencies
    #     Where it hurts
    #     gg = np.fft.ifft(np.fft.fftshift(...))
    #     tf = np.fft.ifft(np.fft.fftshift(self.rawY[:, ii])) * gg
    #     Precompute FFT plan once
    #     from numpy.fft import rfft, irfft
    #     Use real FFTs:
    #     gg = irfft(rfft(kernel))
    #     tf = irfft(rfft(self.rawY[:, ii]) * rfft(gg))
        
    # 6) Batch FFTs
    #     Instead of looping over transitions:
    #     fftY = np.fft.fft(self.rawY, axis=0)
    #     fftY *= gg[:,None]
    #     self.rawY = np.real(np.fft.ifft(fftY, axis=0))
    
    # 7) generateSpectrum() vector math (medium)
    #     Avoid repeated divisions
    #         invX3 = 1.0 / self.rawX**3
    #     Avoid np.argmin(abs(...))
    #         idx = np.searchsorted(g2[::-1], gxyz2[0])
    #         This is O(log N) instead of O(N).
    
    # 8) Spin matrices generation (minor but easy)
    #     Where it hurts
    #     Nested Python loops:
    #     for i, mi in enumerate(m):
    #         for j, mj in enumerate(m):
    #
    #     Replace with vectorized version
    #     m = np.arange(S, -S-1, -1)
    #     mp = np.sqrt(S*(S+1) - m[:-1]*m[1:])
    #     Splus = np.diag(mp, -1)
    #     Sminus = np.diag(mp, 1)
    
    # 9) Printing inside loops (🚨 silent killer)

class optEPR():
    def __init__(self, **kwargs):
        self.Symmetry = 'Ci' #### not sure need that
        self.Treshold = 1e-3
        self.Verbosity = False
        self.nKnots = 20
        
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"Unknown parameter '{key}' for optEPR")
    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"Unknown parameter '{key}' for SysPar")
    def setFromCtrl(self, Dict):
        if type(Dict)==dict:
            for key in Dict.keys():
                if hasattr(self, key):
                    setattr(self, key, Dict[key])
                else:
                    raise AttributeError(f"Unknown parameter '{key}' for Opt")
    def getDefaultDict(self):
        return {'Symmetry': 'Ci', 
        'nKnots': 20, 
        'Verbosity': False, 
        'Treshold': 1e-3,
        }

class EPRsim():
    def __init__(self, Sys=None, Exp=None, Opt=None):
        if not isinstance(Sys, sysPar):
            self.Sys=sysPar()
        else:
            self.Sys=Sys
        
        if not isinstance(Exp, expPar):
            self.Exp=expPar()
        else:
            self.Exp=Exp
        
        if not isinstance(Opt, optEPR):
            self.Opt=optEPR()
        else:
            self.Opt=Opt
            
        self.Zee = [np.array([])]*3
        self.ZfHAM = np.array([])
        self.DET = [np.array([])]*3
        self.BHAM = [np.array([])]*3
        self.AHAM = np.array([])
        self.Ix_tot = np.array([])
        self.Iy_tot = np.array([])
        self.Iz_tot = np.array([])
        
        self.Sx_tot = np.array([])
        self.Sy_tot = np.array([])
        self.Sz_tot = np.array([])
        
        self.Omega = 0.0
        self.X = np.array([])
        self.Y = np.array([])        
        self.rawY = np.array([])    
        self.rawX = np.array([])    
        # Planck constant (J·s)
        self.planck = 6.62607015e-34
        # Bohr magneton (J/T)
        self.bmagn = 9.2740100783e-24 #J/T.
        self.nmagn = 5.050783699e-27 #J/T.
        self.RegenHAM = True

    def run(self):
        
        self.Omega =self.Exp.mwFreq*1e12*self.planck/self.bmagn # GHz to mT
        
        #g1F = self.Omega
        if self.RegenHAM:
            self.generateHam()
            

        evalsA, UHF = eigh(self.AHAM)
        #AHAMinv = UHF @ np.diag(1.0/evalsA) @ UHF.conj().T
        AHAMsq = UHF @ np.diag(np.sqrt(evalsA)) @ UHF.conj().T
        AHAMsq_inv = UHF @ np.diag(1.0 / np.sqrt(evalsA)) @ UHF.conj().T
        
        transitions = np.zeros((self.BHAM[0].shape[0], 6))
        #timest = time.time()
        for ii in range(3):
            ### option 0 (ChatGPT):
            C = AHAMsq_inv @ self.BHAM[ii] @ AHAMsq_inv
            eigvals, eigvecs_C = eigh(C, check_finite=False)
            eigvecs = AHAMsq @ eigvecs_C
            eigvecs /= np.linalg.norm(eigvecs, axis=0, keepdims=True)
            
            # eigvals, eigvecs = eigh(self.BHAM[ii], self.AHAM, check_finite=False)
            
            ### option 1 :
            # AB = self.BHAM[ii]@AHAMinv # np.linalg.solve(self.AHAM, self.BHAM[ii])
            # eigvals, eigvecs = eig(AB, check_finite=False) 
            # note that AB is not hermitian unless AHAMinv is diagonal. 
            # usign eigh can be beneficial for performance, so this may need to be researched further
            
            ### option 2 :
            # eigvals, eigvecs = eig(self.BHAM[ii], self.AHAM,  check_finite=False)
            ### option 3: 
            ## type=1 =>     a @ v = w @ b @ v
            ## type=2 => a @ b @ v = w @ v
            ## type=3 => b @ a @ v = w @ v
            # eigvals, eigvecs = eigh(self.BHAM[ii], AHAMinv,  type=3, check_finite=False)
            # Btilde = UHF.conj().T @ self.BHAM[ii] @ UHF
            # AB = np.linalg.solve(Btilde, np.diag(evalsA))
            # eigvals, eigvecs = eigh(AB)
            #
            ######################################################
            
            prob = self.DET[ii].ravel() @ eigvecs
            prob *= np.conj(prob)
            prob = np.abs(prob)
            
            ori = np.zeros(3)
            ori[ii]=1.0
            HAMZ  = self.Zee[0]*ori[0]
            HAMZ += self.Zee[1]*ori[1]
            HAMZ += self.Zee[2]*ori[2]
            
            Iop = self.Ix_tot*ori[0]+self.Iy_tot*ori[1]+self.Iz_tot*ori[2]
            Sop = self.Sx_tot*ori[0]+self.Sy_tot*ori[1]+self.Sz_tot*ori[2]
            
            for jj, ff in enumerate(np.real(eigvals)):
                #transitions[jj][0]+=prob[jj]
                if 0<ff<np.inf:
                    #fields = np.real(ff)  # option 2
                    fields = np.real(1.0/ff) # option 2
                    #C = np.tile(self.DET[ii].ravel(), (self.BHAM[ii].shape[0], 1))
                    #prob = np.diagonal(C @ eigvecs).copy() #### I freaking hate python... if no copy(), it shows "read only" error
                    
                    if prob[jj]>0.1:
                        Szexp, Izexp = self.findTransition(HAMZ*fields+self.ZfHAM, Sop, Iop)
                        id1 = self.Sys.S[0]+Szexp[0] 
                        id2 = self.Sys.I[0]+np.mean(Izexp)
                        kk = int(np.round(id1*(2*self.Sys.S[0]+1)+id2))
                        transitions[kk][ii]=np.sqrt(prob[jj])
                        transitions[kk][ii+3]=fields
                        # print(f'{ii} # {kk}# {id1:2.3f},{id2:2.3f}: B={fields:05.2f}: {prob[jj]:.3f}')
        #print(f'{time.time()-timest}')
            # if eigvals[ii]>0 and eigvals[ii]<g1F:
            #     sexp = eigvecs[:,5]@sz.ravel()
            #     iexp = eigvecs[:,5]@iz.ravel()
            #     print(f'{eigvals[ii]} - {sexp} : {iexp}')
        #print(transitions)
        nSample = 10
        self.rawY = np.zeros((self.Exp.nPoints*nSample, transitions.shape[0]))
        self.rawX = np.linspace(self.Exp.BMin, self.Exp.BMax, self.Exp.nPoints*10)
        
        if self.Sys.lw_G is not None: ## Gaussian lineshape
            gw = self.Sys.lw_G/2.355  
        elif self.Sys.lw is not None:
            gw = self.Sys.lw[0]/2.355
            # gg =  np.fft.ifft(np.fft.fftshift(np.exp(-(self.rawX-np.mean(self.rawX) )**2/gw**2/2.0).astype(complex)))
        else:
            gw = 0
        
        ifftline = None
        if gw>0:
            if self.Exp.Harmonic==1:
                ga=-(self.rawX-np.mean(self.rawX))/gw**2
                ifftline =  ifft(fftshift(ga*np.exp(-(self.rawX-np.mean(self.rawX) )**2/gw**2/2.0).astype(complex)))
            else:
                ifftline =  ifft(fftshift(np.exp(-(self.rawX-np.mean(self.rawX) )**2/gw**2/2.0).astype(complex)))
            
        if self.Sys.lw_L is not None: ## Lorentzian lineshape
            llw = self.Sys.lw_L
            
            if self.Exp.Harmonic==1:
                ll = ifft(fftshift((-(self.rawX-np.mean(self.rawX))/( (self.rawX-np.mean(self.rawX))**2 
                                        + llw**2/4.0      )**2).astype(complex)))
            else:
                ll = ifft(fftshift(1.0/( (self.rawX-np.mean(self.rawX))**2 
                                        + llw**2/4.0      ).astype(complex)))
            if ifftline is not None:
                ifftline+=ll
            else:
                ifftline = ll
                
        if ifftline is None: ifftline= 1.0
        
        self.X = np.linspace(self.Exp.BMin, self.Exp.BMax, self.Exp.nPoints) # in mT
        self.Y = np.zeros_like(self.X) 
        
        ra1 = int(self.Exp.nPoints/2)
        ra2 = int(self.Exp.nPoints*nSample-ra1)
        
        for ii in range(transitions.shape[0]):
            if np.sum(transitions[ii][:3])>0.01:
                amp = np.mean(transitions[ii][:3])
                fields = transitions[ii][3:]
                self.generateSpectrum(ii, amp, fields)
                
                tf = ifft(fftshift(self.rawY[:, ii]))*ifftline
                self.rawY[:, ii] = np.real(fftshift(fft( tf )))
                self.Y += np.real(fftshift(fft( np.delete(tf, np.s_[ra1:ra2]) )))
                
                
        # tf = ifft(fftshift(self.rawY, axes=(0,)), axis=0)
        # tf *=gg[:,None]
        # self.rawY = np.real(fftshift(fft( tf , axis=0), axes=(0,)))
        # self.Y = np.real(np.sum( fftshift(fft( np.delete(tf, np.s_[ra1:ra2], axis=0) ,axis=0),  axes=(0,)), axis=1) )
    def generateSpectrum(self, ii, amp, fields):
        h2 = np.sort(fields**2)[::-1]
        h2x, h2y, h2z = h2[0], h2[1], h2[2]
        # h2x > h2y > h2z
        Xh2 = self.rawX**2
        
        idx1 = np.argmin(np.abs(Xh2-h2x))
        idx2 = np.argmin(np.abs(Xh2-h2y))
        idx3 = np.argmin(np.abs(Xh2-h2z))
        
        if Xh2[idx1]>h2x: idx1-=1
        if Xh2[idx2]<h2y: idx2+=1
        idx21 = idx2
        if Xh2[idx3]<h2z: idx3+=1
        
        k2_1 = (h2y-h2z)*(h2x-Xh2[idx2:idx1])/(Xh2[idx2:idx1]-h2z)/(h2x-h2y)
        pre1 =  np.sqrt(h2x*h2y*h2z)/Xh2[idx2:idx1]/np.sqrt((Xh2[idx2:idx1]-h2z)*(h2x-h2y))
        
        k2_2 = (Xh2[idx3:idx21]-h2z)*(h2x-h2y)/(h2y-h2z)/(h2x-Xh2[idx3:idx21])
        pre2 =  np.sqrt(h2x*h2y*h2z)/Xh2[idx3:idx21]/np.sqrt((h2y-h2z)*(h2x-Xh2[idx3:idx21]))
        
        self.rawY[idx2:idx1,  ii] += amp*pre1*ellipk(k2_1) #self.Omega**2
        self.rawY[idx3:idx21, ii] += amp*pre2*ellipk(k2_2) #self.Omega**2
        
        #self.rawY[:, ii]*=self.rawX        
    # def generateSpectrum(self, ii, amp, fields):
    #     gxyz2 = np.sort(self.Omega/fields)**2 # don't do fields.sort()... it is by pointer, so it will overwrite the source (stupid python crap)
        
    #     g2 = (self.Omega/self.rawX)**2 ### note that the index will run backwards (high g^2 first)! 
    #     idx1 = np.argmin(np.abs(g2-gxyz2[0]))
    #     # if g2[idx1]<gxyz2[0]: idx1-=1
        
    #     idx2 = np.argmin(np.abs(g2-gxyz2[1]))
    #     if g2[idx2]>gxyz2[1]: idx2+=1
    #     idx21 = idx2  ### because stupid numpy does not include the last bit
    #     idx3 = np.argmin(np.abs(g2-gxyz2[2]))
        
    #     if g2[idx3]>gxyz2[2]: idx3+=1
    #     k2_1 = (gxyz2[2]-gxyz2[1])*(g2[idx2:idx1]-gxyz2[0])/(gxyz2[2]-g2[idx2:idx1])/(gxyz2[1]-gxyz2[0])
    #     k2_2 = (gxyz2[2]-g2[idx3:idx21])*(gxyz2[1]-gxyz2[0])/(gxyz2[2]-gxyz2[1])/(g2[idx3:idx21]-gxyz2[0])
        
    #     pre1 = (gxyz2[2]-gxyz2[0])*(gxyz2[2]-gxyz2[1])/(gxyz2[2]-g2[idx2:idx1])/(gxyz2[1]-gxyz2[0])
    #     pre2 = (gxyz2[2]-gxyz2[0])/(g2[idx3:idx21]-gxyz2[0])
        
    #     self.rawY[idx2:idx1,  ii] += amp*self.Omega**2*pre1*ellipk(k2_1)
    #     self.rawY[idx3:idx21, ii] += amp*self.Omega**2*pre2*ellipk(k2_2)
        
    #     self.rawY[:, ii]*=self.rawX
        
    def findTransition(self, HAM, Sop, Iop):
        d, v = eigh(HAM, check_finite=False)
        
        # #### FOLLOWNG CODE IS FASTER FOR LARGE MATRIXES BUT SUCKS FOR EPR
        # idx = np.argsort(d)
        # ds = d[idx]
        # targets = ds - self.Omega
        # j = np.searchsorted(ds, targets)
        # j0 = np.clip(j - 1, 0, len(ds) - 1)
        # j1 = np.clip(j,     0, len(ds) - 1)
        # err0 = np.abs(ds - ds[j0] - self.Omega)
        # err1 = np.abs(ds - ds[j1] - self.Omega)
        # use_j1 = err1 < err0
        # best_j = np.where(use_j1, j1, j0)
        # best_err = np.where(use_j1, err1, err0)
        # i_best = np.argmin(best_err)
        # rr = idx[i_best]
        # cc = idx[best_j[i_best]]
        
        dd = d[:, None] - d[None, :] #
        aa = np.argmin(np.abs(dd-self.Omega))
        rr,cc=np.unravel_index(aa, dd.shape)
        
        iizz = np.real(v.conj().T@(Iop) @ v)
        sszz = np.real(v.conj().T@(Sop) @ v)
        Izexp = [iizz[rr,rr], iizz[cc, cc]]
        Szexp = [sszz[rr,rr], sszz[cc, cc]]    
        
        # mm = np.max(np.abs(np.diag(iizz)))
        # ss = np.sum(np.abs(np.diag(iizz)))
        # np.vdot(v, (self.Ix_tot*LL[0]+self.Iy_tot*LL[1]+self.Iz_tot*LL[2]) @ v)
        # v1 = v[:, rr]
        # v2 = v[:, cc]
        # Izexp = [np.real(v1.conj().T@(Iop) @ v1), np.real(v2.conj().T@(Iop) @ v2)]
        # Szexp = [np.real(v1.conj().T@(Sop) @ v1), np.real(v2.conj().T@(Sop) @ v2)]
        #print(f'{dd[rr,cc]-self.Omega:.4f}')
        # Izexp = [ np.real(np.vdot(v1, (self.Ix_tot*LL[0]+self.Iy_tot*LL[1]+self.Iz_tot*LL[2]) @ v1)),
        #           np.real(np.vdot(v2, (self.Ix_tot*LL[0]+self.Iy_tot*LL[1]+self.Iz_tot*LL[2]) @ v2))]
   
        
        return Szexp, Izexp
                  #exp = v1@(self.Ix_tot*LL[0]+self.Iy_tot*LL[1]+self.Iz_tot*LL[2])@v2
        #
    def generateHam(self):
        # I want to separate this one in case I don't need to regenerate it, e.g. if we only change Exp or Opt parameters

        
        if len(self.Sys.S)>1:
            raise('Currently, only one e-spin supported')
        if len(self.Sys.I)>1:
            raise('Currently, only one nuc-spin supported')
        
        Sx, Sy, Sz = self.spin_matrices(self.Sys.S[0])
        g = self.Sys.g[0]

        if self.Sys.I[0]>0:
            Ix, Iy, Iz = self.spin_matrices(self.Sys.I[0])
            A = self.Sys.A[0]*1e9*self.planck/self.bmagn # MHz to mT
            gn = self.Sys.gn[0]
            doHF=True
        else:
            Ix = np.ones(1.0)
            Iy = np.ones(1.0)
            Iz = np.ones(1.0)
            doHF=False

        sE = np.eye(int(self.Sys.S[0]*2.0+1.0))
        iE = np.eye(int(self.Sys.I[0]*2.0+1.0))

        self.Sx_tot = np.kron(Sx, iE)
        self.Sy_tot = np.kron(Sy, iE)
        self.Sz_tot = np.kron(Sz, iE)

        self.Ix_tot = np.kron(sE, Ix)
        self.Iy_tot = np.kron(sE, Iy)
        self.Iz_tot = np.kron(sE, Iz)

        self.Zee[0] = self.Sx_tot*g[0]+self.Ix_tot*gn*self.nmagn/self.bmagn
        self.Zee[1] = self.Sy_tot*g[1]+self.Iy_tot*gn*self.nmagn/self.bmagn
        self.Zee[2] = self.Sz_tot*g[2]+self.Iz_tot*gn*self.nmagn/self.bmagn

        self.DET[0] = -self.Sz_tot
        self.DET[1] = -self.Sz_tot
        self.DET[2] = +self.Sx_tot

        if doHF:
            self.ZfHAM = (self.Sx_tot.T@self.Ix_tot*A[0] 
                        + self.Sy_tot.conj().T@self.Iy_tot*A[1]
                        + self.Sz_tot@self.Iz_tot*A[2])
        else:
            self.ZfHAM = np.zeros_like(self.Sz_tot)
        
        Ident = np.eye(self.Sz_tot.shape[0])

        self.BHAM[0] = np.kron(self.Zee[0], Ident) - np.kron(Ident, self.Zee[0].T )
        self.BHAM[1] = np.kron(self.Zee[1], Ident) - np.kron(Ident, self.Zee[1].T )
        self.BHAM[2] = np.kron(self.Zee[2], Ident) - np.kron(Ident, self.Zee[2].T )

        self.AHAM = np.kron(Ident, Ident)*self.Omega - np.kron(self.ZfHAM, Ident) + np.kron(Ident, self.ZfHAM.T)  # HF coupling goes here

    def spin_matrices(self, S):
        """Return Sx, Sy, Sz for spin S."""
        dim = int(2*S+1)
        m = np.arange(S, -S-1, -1)  # m = S, S-1, ..., -S
        Splus = np.zeros((dim,dim), dtype=complex)
        Sminus = np.zeros((dim,dim), dtype=complex)

        for i, mi in enumerate(m):
            for j, mj in enumerate(m):
                if mj == mi-1:  # lowering from mi to mj
                    Sminus[i,j] = np.sqrt(S*(S+1) - mj*(mj+1))
                if mj == mi+1:  # raising
                    Splus[i,j] = np.sqrt(S*(S+1) - mj*(mj-1))

        Sx = 0.5*(Splus + Sminus)
        Sy = -0.5j*(Splus - Sminus)
        Sz = np.diag(m).astype(complex)
        return Sx, Sy, Sz

if __name__ == "__main__":
    import time
    es = EPRsim()
    es.Sys.set( S = [1/2],
            g = [np.array([1.981, 1.979, 1.944])], 
            Nucs = ['51V'],
            A = [np.array([519, 185, 192])],
            lw = [10],
            )
    es.Exp.set(mwFreq=9.5, BMin=200, BMax=500, nPoints=1000, Harmonic=1)
    timest = time.time()
    es.run()
    timeed = time.time()-timest
    print(f'time: {timeed}')
    plt.figure(2)
    plt.clf()
    plt.plot(es.rawX, es.rawY,'r')
    plt.plot(es.X, es.Y,'b')
    plt.show()