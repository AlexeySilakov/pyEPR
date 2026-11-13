import time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from spin_system.SysPar import expPar
from spin_system.hyscore_sim import optHYSCORE, HYSCOREsim
from spin_system import mathfunctions as mf


class expENDOR(expPar):
    """Experimental settings of an ENDOR simulation: expPar plus the radio-frequency window and the ENDOR sequences."""
    def __init__(self, **kwargs):
        super().__init__()
        self.MinFreq = 0.0       # MHz, start of the radio-frequency axis
        self.MaxFreq = 30.0      # MHz, end of the radio-frequency axis
        self.nPoints = 1024      # points of the radio-frequency axis
        self.Sequence = ['data', ['data', 'ideal', 'Mims', 'Davies']]   # 'data' -> per-dataset choice
        self.set(**kwargs)

    def getToolTips(self, param=None):
        ss = 'If "-1", take from data'
        dc = {'mwFreq': 'Microwave frequency in [GHz]. '+ss,
              'Field': 'Magnetic field in [mT]. '+ss,
              'MinFreq': 'Start of the simulated radio-frequency axis in [MHz]. '+ss,
              'MaxFreq': 'End of the simulated radio-frequency axis in [MHz]. '+ss,
              'nPoints': 'Number of points of the simulated radio-frequency axis. '+ss,
              'tau': 'Mims ENDOR: delay between the first two microwave pulses [ns]. '+ss,
              'Tinv': 'Davies ENDOR: length of the inversion pulse [ns]. '+ss,
              'Sequence': ('ENDOR pulse sequence. \n'+
                           ' "data" = use the sequence set for each dataset \n'+
                           ' "ideal" = no sequence factor \n'+
                           ' "Mims" = every line times sin(pi*a*tau)**2 \n'+
                           ' "Davies" = every line times 1.4*eta/(0.49+eta**2), eta = a*Tinv \n'+
                           'a = hyperfine shift of the electron resonance for that nuclear transition'),
              'ExciteWidth': 'Excitation width in [MHz]; 0 = no orientation selection',
              }
        if param is None:
            return dc
        return dc.get(param, '')

    def getDefaultDict(self):
        return {'mwFreq': float(-1.0),
                'Field': float(-1.0),
                'Sequence': ['data', ['data', 'ideal', 'Mims', 'Davies']],
                'tau': float(-1.0),
                'Tinv': float(-1.0),
                'MinFreq': float(-1.0),
                'MaxFreq': float(-1.0),
                'nPoints': int(-1),
                'ExciteWidth': float(20.0),
                }


class optENDOR(optHYSCORE):
    """Options of ENDORsim: the orientation grid and selection of optHYSCORE plus the enhancement switch."""
    def __init__(self, **kwargs):
        self.Enhancement = False      # hyperfine enhancement of the radio-frequency field
        super().__init__(**kwargs)

    def getDefaultDict(self):
        return {'nKnots': 20,
                'Grid': ['fibonacci', ['fibonacci', 'sphgrid', 'spiral']],
                'OriSelType': ['g_eff', ['g_eff', 'g_eff+HFC', 'brute force', 'precalculated']],
                'Enhancement': False,
                }

    def getToolTips(self, param=None):
        dc = super().getToolTips()
        dc['Enhancement'] = ("Hyperfine enhancement of the radio-frequency field.\n"+
                             " off = transition moment of b.I, b perpendicular to B0 \n"+
                             " on  = b.I replaced by b.(1 - mS*g.A/(g_eff*nu_n)).I")
        if param is None:
            return dc
        return dc.get(param, '')


class ENDORsim(HYSCOREsim):
    """Powder ENDOR spectrum of an S = 1/2 system, on the orientation grid and selection weights of HYSCOREsim."""
    def __init__(self, Sys=None, Exp=None, Opt=None, errorFunc=AttributeError):
        super().__init__(Sys=Sys,
                         Exp=Exp if isinstance(Exp, expENDOR) else expENDOR(),
                         Opt=Opt if isinstance(Opt, optENDOR) else optENDOR(),
                         errorFunc=errorFunc)
        self.Parts = []      # spectrum of each nucleus, indexed as Sys.Nucs (None for 'ori.sel.only'); Spectrum is their sum

    def preCompute(self):
        """Nothing to prepare: the ENDOR operators are small and built inside compute_transitions."""

    def sequence_name(self):
        """'ideal', 'Mims' or 'Davies'; Exp.Sequence is the resolved string or the [value, [choices]] list of the Exp tab."""
        s = getattr(self.Exp, 'Sequence', 'ideal')
        if isinstance(s, (list, tuple)):
            s = s[0] if len(s) else 'ideal'
        return s if s in ('Mims', 'Davies') else 'ideal'

    def reRun(self):
        sta = time.perf_counter()
        if not self.Exp.Field or self.Exp.Field <= 0:
            self.errorFunc('ENDOR_sim_reRun: Exp.Field is not set')
            return
        if not (self.Exp.MaxFreq > self.Exp.MinFreq) or int(self.Exp.nPoints) < 2:
            self.errorFunc(f'ENDOR_sim_reRun: empty radio-frequency axis, MinFreq={self.Exp.MinFreq}, '
                           f'MaxFreq={self.Exp.MaxFreq}, nPoints={self.Exp.nPoints}')
            return
        seq = self.sequence_name()
        if seq == 'Mims' and not self.Exp.tau > 0:
            self.errorFunc('ENDOR_sim_reRun: Mims ENDOR needs Exp.tau > 0')
            return
        if seq == 'Davies' and not self.Exp.Tinv > 0:
            self.errorFunc('ENDOR_sim_reRun: Davies ENDOR needs Exp.Tinv > 0')
            return
        if self.Opt.OriSelType[0] == 'precalculated':
            inp = self.Opt.OriSelInp
            if type(inp) is not np.ndarray or inp.ndim != 2 or inp.shape[1] != 3:
                self.errorFunc(f'ENDOR_sim_reRun: expected a 3-column Opt.OriSelInp, got {np.shape(inp)}. \n - Opt.OriSelType={self.Opt.OriSelType}')
                return
            self.phi = inp[:, 0]
            self.theta = inp[:, 1]
            self.ak = inp[:, 2]
            self.tri = None
            self.triareas = None
            self.gridSig = None       # a loaded grid is not one we can re-derive
        else:
            self.ensure_grid()        # the grid cache of HYSCOREsim: rebuilt only when grid_signature() changes
        self.compute_transitions()
        if self.verbose:
            print(f"ENDOR total in {(time.perf_counter()-sta)*1000:.2f} ms")

    def compute_transitions(self):
        Sys, Exp = self.Sys, self.Exp
        seq = self.sequence_name()
        npts = int(Exp.nPoints)
        self.X = np.linspace(Exp.MinFreq, Exp.MaxFreq, npts)
        dx = self.X[1]-self.X[0]
        self.Spectrum = np.zeros(npts)
        self.Parts = []

        sel = self.ak != 0           # orientations outside the excitation window cost nothing
        ak = self.ak[sel]
        st, ct = np.sin(self.theta[sel]), np.cos(self.theta[sel])
        sp, cp = np.sin(self.phi[sel]), np.cos(self.phi[sel])
        lxyz = np.column_stack((st*cp, st*sp, ct))            # direction of B0 in the g frame
        g = np.asarray(Sys.g[0], dtype=float)
        geff = np.linalg.norm(lxyz*g, axis=1)
        u = lxyz*g/geff[:, None]                               # quantisation axis of the electron spin, as in local_sham
        # two radio-frequency field directions perpendicular to B0; the intensity is their average
        perp = (np.column_stack((ct*cp, ct*sp, -st)), np.column_stack((-sp, cp, np.zeros_like(sp))))

        for n_nuc in range(len(Sys.I)):
            if Sys.useFor[n_nuc] == 'ori.sel.only':
                self.Parts.append(None)
                continue
            spin = Sys.I[n_nuc]
            mult = int(2*spin+1)
            Iop = np.array([mf.SpinOp([spin], [cc]) for cc in 'xyz'])

            if Sys.A[n_nuc].ndim == 1:
                A = mf.RotTransDiag(Sys.A[n_nuc], Sys.Apa[n_nuc]*np.pi/180.0)
            else:
                A = mf.RotTrans3x3(Sys.A[n_nuc], Sys.Apa[n_nuc]*np.pi/180.0)
            quad = np.zeros((mult, mult), dtype=complex)
            if Sys.Q and mult > 2 and np.sum(np.abs(Sys.Q[n_nuc])) != 0:
                if Sys.Q[n_nuc].ndim == 1:
                    Q = mf.RotTransDiag(Sys.Q[n_nuc], Sys.Qpa[n_nuc]*np.pi/180.0)
                else:
                    Q = mf.RotTrans3x3(Sys.Q[n_nuc], Sys.Qpa[n_nuc]*np.pi/180.0)
                quad = np.einsum('ab,aij,bjk->ik', Q, Iop, Iop)

            # the two nuclear sub-Hamiltonians of local_sham, for all orientations at once
            nu_n = Sys.gamma[n_nuc]*Exp.Field
            zeeman = -nu_n*np.einsum('nc,cij->nij', lxyz, Iop)
            hfi = np.einsum('nc,cij->nij', u@A.T, Iop)        # I.A.u

            lw = Sys.lwEndor[n_nuc] if Sys.lwEndor else 0.0
            sigma = lw/2.354820045/dx                          # 2.354820045 = 2sqrt(2ln2); in points
            pad = int(np.ceil(4*sigma))                        # lines just outside the window still reach into it
            stick = np.zeros(npts+2*pad)
            iu, ju = np.triu_indices(mult, 1)
            for mS in (-0.5, 0.5):
                D, V = np.linalg.eigh(zeeman + mS*hfi + quad)
                Vh = np.conj(np.transpose(V, (0, 2, 1)))
                moment = 0.0
                for b in perp:
                    c = b
                    if self.Opt.Enhancement:
                        # second-order mixing of the electron spin by the radio-frequency field: b.(1 - mS*g.A/(g_eff*nu_n)).I
                        w = b*g
                        w = w - np.sum(w*u, axis=1)[:, None]*u
                        c = b - (mS/(nu_n*geff))[:, None]*(w@A.T)
                    T = Vh@np.einsum('nc,cij->nij', c, Iop)@V
                    moment = moment + 0.5*np.abs(T)**2
                amp = moment[:, iu, ju]*ak[:, None]
                freq = D[:, ju]-D[:, iu]
                # hyperfine shift of the electron resonance when the nucleus goes from level i to level j
                shift = np.real(np.einsum('nij,njk,nki->ni', Vh, hfi, V))
                a = shift[:, ju]-shift[:, iu]
                if seq == 'Mims':
                    amp = amp*np.sin(np.pi*a*Exp.tau*1e-3)**2
                elif seq == 'Davies':
                    eta = np.abs(a)*Exp.Tinv*1e-3
                    amp = amp*1.4*eta/(0.49+eta**2)

                # each line is split between its two neighbouring points, so positions are not rounded to the axis
                pos = (freq.ravel()-self.X[0])/dx + pad
                lo = np.floor(pos).astype(np.int64)
                ok = (lo >= 0) & (lo < stick.size-1)
                frac = (pos-lo)[ok]
                stick += np.bincount(lo[ok], weights=amp.ravel()[ok]*(1-frac), minlength=stick.size)
                stick += np.bincount(lo[ok]+1, weights=amp.ravel()[ok]*frac, minlength=stick.size)

            if sigma > 0:
                stick = gaussian_filter1d(stick, sigma, mode='constant')
            nrep = max(1, int(Sys.nNucs[n_nuc])) if Sys.nNucs is not None else 1
            self.Parts.append(nrep*stick[pad:pad+npts]/dx)   # per MHz, so the scale does not depend on nPoints
            self.Spectrum += self.Parts[-1]

    def computeEPR(self, fields):
        """EPR as the orientation selection sees it: for every field [mT] in `fields`, the sum over the grid
        of the selection weights at Exp.mwFreq, with the current Grid, nKnots, OriSelType and ExciteWidth."""
        self.make_grid(grid=self.grid_name(), nKnots=self.Opt.nKnots)
        self.gridSig = None           # the weights are left at the last field, not at Exp.Field
        out = np.zeros(len(fields))
        if self.Opt.OriSelType[0] in ('g_eff', 'g_eff+HFC'):
            # compute_weights itself, field by field; away from resonance its 'no resonances' is the answer, not a fault
            errorFunc, field = self.errorFunc, self.Exp.Field
            self.errorFunc = lambda message: None
            try:
                for kk, B in enumerate(fields):
                    self.Exp.Field = B
                    self.compute_weights()
                    out[kk] = np.sum(self.ak)
            finally:
                self.errorFunc, self.Exp.Field = errorFunc, field
            return out

        # 'brute force' (and 'precalculated', which stores its results): computeOrisel_eig for all orientations at once
        sig = self.Exp.ExciteWidth/2.35482
        mwFreqMHz = self.Exp.mwFreq*1e3
        Ssys, HFc, gamma = [self.Sys.S[0]], [], []
        for n_nuc in range(len(self.Sys.nNucs)):
            if self.Sys.useFor[n_nuc] == 'sim.only': continue
            for ii in range(int(self.Sys.nNucs[n_nuc])):
                if self.Sys.A[n_nuc].ndim == 1:
                    HFc.append(mf.RotTransDiag(self.Sys.A[n_nuc], self.Sys.Apa[n_nuc]*np.pi/180.0))
                else:
                    HFc.append(mf.RotTrans3x3(self.Sys.A[n_nuc], self.Sys.Apa[n_nuc]*np.pi/180.0))
                Ssys.append(self.Sys.I[n_nuc])
                gamma.append(self.Sys.gamma[n_nuc])
        # x, y, z operators of spin number pos in the product space of the electron and those nuclei
        ops = lambda pos: np.array([mf.SpinOp(Ssys, ['e']*pos + [cc] + ['e']*(len(Ssys)-pos-1)) for cc in 'xyz'])
        st, ct = np.sin(self.theta), np.cos(self.theta)
        sp, cp = np.sin(self.phi), np.cos(self.phi)
        lxyz = np.column_stack((st*cp, st*sp, ct))
        S = ops(0)
        zeeman = np.einsum('nc,cij->nij', lxyz*np.asarray(self.Sys.g[0], dtype=float), S)/self.h_ov_beta_mTMHz   # per mT
        HFHam = np.zeros(S.shape[1:], dtype=complex)
        for kk in range(len(gamma)):
            I = ops(kk+1)
            zeeman += gamma[kk]*np.einsum('nc,cij->nij', lxyz, I)
            HFHam += np.einsum('ab,aij,bjk->ik', HFc[kk], I, S)
        # microwave field along the direction at theta + 90 degrees, as there
        DHam = np.einsum('nc,cij->nij', np.column_stack((ct*cp, ct*sp, -st)), S)
        iu, ju = np.triu_indices(S.shape[1], 1)
        for kk, B in enumerate(fields):
            D, V = np.linalg.eigh(B*zeeman + HFHam)
            moment = np.abs(np.conj(np.transpose(V, (0, 2, 1)))@DHam@V)[:, iu, ju]
            df = np.abs(np.abs(D[:, iu]-D[:, ju]) - mwFreqMHz)
            ak = np.sum(np.where(df < 5.0*sig, moment*np.exp(-0.5*(df/sig)**2), 0.0), axis=1)
            out[kk] = np.sum(ak*self.ww)
        return out


# ----------------------------------------------------------------------
# Example usage
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    es = ENDORsim()
    es.Sys.set(g=[np.array([2.2, 2.1, 2.0])],
               S=[0.5],
               A=[np.array([2.0, 2.0, 5.0]), np.array([4.0, 4.0, 6.0])],
               Apa=[np.array([0, 0, 0]), np.array([0, 0, 0])],
               Q_K=[0.0, 0.8],
               Q_eta=[0.0, 0.3],
               lwEndor=[0.1, 0.2],
               useFor=['all', 'all'],
               Nucs=['1H', '14N'])      # last: setNucs builds Q from Q_K and Q_eta
    es.Exp.set(mwFreq=9.5, Field=325, ExciteWidth=50, nPoints=2048, MinFreq=0, MaxFreq=25)
    es.Opt.nKnots = 60

    fig = plt.figure(1, figsize=(6, 4))
    ax = fig.add_subplot(111)
    for seq in ('ideal', 'Mims', 'Davies'):
        es.Exp.set(Sequence=seq, tau=200, Tinv=100)
        es.run()
        ax.plot(es.X, es.Spectrum/np.max(es.Spectrum), label=seq)
    ax.set_xlabel('RF, MHz')
    ax.legend()
    plt.show()
