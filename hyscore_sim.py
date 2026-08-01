import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from SysPar import sysPar, expPar
import mathfunctions as mf
import time



class optHYSCORE():
    def __init__(self, **kwargs):
        self.Symmetry = 'Ci'
        self.nKnots = 20
        self.Grid = ['fibonacci', ['fibonacci', 'sphgrid', 'spiral']]
        self.OriSelInp = None #'orisel.mat'; alternatively, one can load pregenerated grid
                              # it should contain three columns : "phi", "theta" and "weights"
        self.Verbosity = False
        self.KillNeg = True           # flag for suppression of negative time in timedomain spectrum
        self.ProdRule = True      # use product rule, producing combination frequencies (memory/time consuming)
        self.AmpRatios = None # relative amplitudes of the peaks from different HF couplings
        self.Treshold = 1e-3
        self.OriSelType = ['g_eff', ['g_eff', 'g_eff+HFC', 'brute force', 'precalculated']]
        self.errorFunc=AttributeError
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.errorFunc(f"Unknown parameter '{key}' for optHYSCORE")
        
    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.errorFunc(f"Unknown parameter '{key}' for SysPar")
    def setFromCtrl(self, Dict):
        if type(Dict)==dict:
            for key in Dict.keys():
                if hasattr(self, key):
                    setattr(self, key, Dict[key])
                else:
                    self.errorFunc(f"Unknown parameter '{key}' for Opt")
    def getDefaultDict(self):
        return {'Symmetry': 'Ci',
        'nKnots': 20,
        'Grid': ['fibonacci', ['fibonacci', 'sphgrid', 'spiral']],
        'Verbosity': False,
        'KillNeg': True,           # flag for suppression of negative time in timedomain spectrum
        'ProdRule': True,      # use product rule, producing combination frequencies (memory/time consuming)
        'AmpRatios': 1.0, # relative amplitudes of the peaks from different HF couplings
        'Treshold': 1e-3,
        'OriSelType': ['g_eff', ['g_eff', 'g_eff+HFC', 'brute force', 'precalculated']]
        }
    def getToolTips(self, param=None):
        dc = {'nKnots':'number of Knots in the orientation grid',
              'Grid': ("Orientation grid construction.\n"+
                       " 'fibonacci' = golden-angle spiral, equal-area, uniform weights \n"+
                       " 'sphgrid' = SOPHE/triangular grid (EasySpin Ci), non-uniform weights \n"+
                       " 'spiral' = latitude-ring spiral over the full sphere \n"+
                       "nKnots is rescaled per method so the knot count stays comparable"),
              'OriSelInp': 'Input matrix of orientations. Must have three columns [phi, theta, weights]',
              'AmpRatios': 'Modify relative absolute amplitudes of crosspeaks from different nuclei',
              'OriSelType': ("Orientation selection algorithm.\n"+
                             " 'g_eff' = only based on g-anisotropy \n"+
                             " 'g_eff+HFC' = also account for HFC \n"+
                             " 'brute force' = fully diagonalize spin Hamiltonian for all grid points \n"+
                             " 'precalculated '= load from Opt.OriSelInp (e.g. generated via computeOrisel_eig)")
              }
        if param is None:
            return dc
        else:
            if param in dc.keys():
                return dc[param]
            else:
                return ''
class HYSCOREsim():
    def __init__(self, Sys=None, Exp=None, Opt=None, errorFunc=AttributeError):
        if not isinstance(Sys, sysPar):
            self.Sys=sysPar()
        else:
            self.Sys=Sys
        
        if not isinstance(Exp, expPar):
            self.Exp=expPar()
        else:
            self.Exp=Exp
        
        if not isinstance(Opt, optHYSCORE):
            self.Opt=optHYSCORE()
        else:
            self.Opt=Opt
        self.errorFunc=errorFunc
        # Planck constant (J·s)
        self.planck = 6.62607015e-34

        # Bohr magneton (J/T)
        self.bmagn = 9.2740100783e-24
        self.nmagn = 5.050783699e-27 #J/T.
        self.h_ov_beta_mTMHz = self.planck/self.bmagn*1e9 # mT/MHz
        self.h_ov_nmagn_mTMHz = self.planck/self.nmagn*1e9 # mT/MHz
        
        self.sIx = [None]     # index = dimension (2*I+1)
        self.sIy = [None]
        self.sIz = [None]
        
        self.sIxx = [None]
        self.sIxy = [None]
        self.sIxz = [None]

        self.sIyx = [None]
        self.sIyy = [None]
        self.sIyz = [None]

        self.sIzx = [None]
        self.sIzy = [None]
        self.sIzz = [None]
    
        self.HFIx = []
        self.HFIy = []
        self.HFIz = []
    
        self.Quad = []
    
        self.cii = [None]
        self.ckk = [None]
        self.cll = [None]
        self.cnn = [None]
        self.cjj = [None]
        self.cmm = [None]
        self.geff= None
        self.orisel = None
        self.phi = None
        self.theta = None
        self.ak = None
        self.ww = None
        self.tri = None            # (nTri,3) knot indices, or None if not built
        self.triareas = None       # solid angle per triangle, sums to 4*pi
        self.lxyz = None
        
        self.rawSpec = None
        self.X = None
        self.Y = None
        self.verbose = True
        self.Spectrum = None
        # Grid construction. Left as None so Opt.Grid (set from the Opt tab)
        # decides; set it to a name to override the Opt setting.
        self.Grid = None
        self.gridSig = None        # signature of the grid currently in self.phi

    def grid_name(self):
        """Selected grid construction: 'fibonacci', 'sphgrid' or 'spiral'.

        self.Grid wins if set; otherwise Opt.Grid, which the property grid
        stores in the [value, [choices]] form used by OriSelType.
        """
        g = self.Grid if self.Grid is not None else getattr(self.Opt, 'Grid', 'fibonacci')
        if isinstance(g, (list, tuple)):
            g = g[0] if len(g) else 'fibonacci'
        return g

    def get_orisel(self, file):
        pass
    def preCompute(self):
        self.build_spin_operators()

    # =====================================================================
    #  >>> GRID CACHE <<<   grid_signature / ensure_grid
    #
    #  Look here if the orientation grid is rebuilt too often, or -- worse --
    #  not rebuilt when it should be. grid_signature() is the single list of
    #  everything the grid and its ak weights depend on; ensure_grid() only
    #  rebuilds when that list changes.
    #
    #  If you add a parameter that must force a grid rebuild, add it to
    #  grid_signature(). If you add one that must NOT, keep it out.
    # =====================================================================
    def grid_signature(self):
        """Everything the orientation grid and its ak weights depend on.

        Excluded on purpose:
          - parameters that only shape the spectrum once the orientations are
            known (tau, nPoints, MaxFreq, ProdRule, KillNeg, Treshold, ...)
          - nuclei flagged useFor='sim.only', which contribute to the spectrum
            but are skipped by computeOrisel_eig and so cannot move the grid.
        """
        def tag(v):
            # canonical, exactly comparable form -- no rounding, so a change of
            # any size counts as a change
            if isinstance(v, np.ndarray):
                return (v.shape, v.dtype.str, v.tobytes())
            if isinstance(v, (list, tuple)):
                return tuple(tag(x) for x in v)
            if isinstance(v, dict):
                return tuple((k, tag(v[k])) for k in sorted(v, key=repr))
            return v

        Sys, Exp, Opt = self.Sys, self.Exp, self.Opt

        # Which nuclei may influence the grid: everything except 'sim.only'
        nucs = getattr(Sys, 'Nucs', None)
        useFor = getattr(Sys, 'useFor', None)
        keep = []
        for ii in range(len(nucs) if nucs is not None else 0):
            flag = useFor[ii] if (useFor is not None and ii < len(useFor)) else 'all'
            if flag != 'sim.only':
                keep.append(ii)

        sig = [('grid', self.grid_name()),
               ('nKnots', tag(getattr(Opt, 'nKnots', None))),
               ('Symmetry', tag(getattr(Opt, 'Symmetry', None))),
               ('useFor', tag(useFor)),
               ('OriSelType', tag(getattr(Opt, 'OriSelType', [None])[0])),
               ('OriSelInp', tag(getattr(Opt, 'OriSelInp', None)))]

        # Experimental settings that enter the orientation weighting
        for name in ('mwFreq', 'Field', 'gField', 'ExciteWidth'):
            sig.append((name, tag(getattr(Exp, name, None))))

        # Electron part of the spin system
        for name in ('S', 'g', 'gFrame', 'D', 'DFrame'):
            sig.append((name, tag(getattr(Sys, name, None))))

        # Per-nucleus parameters, restricted to the nuclei kept above
        for name in ('Nucs', 'nNucs', 'I', 'gn', 'gamma',
                     'A', 'Apa', 'AFrame', 'Q', 'Qpa', 'QFrame'):
            v = getattr(Sys, name, None)
            if v is None:
                sig.append((name, None))
                continue
            try:
                sig.append((name, tuple(tag(v[ii]) for ii in keep)))
            except (TypeError, IndexError, KeyError):
                sig.append((name, tag(v)))    # not per-nucleus after all

        return tuple(sig)

    def ensure_grid(self):
        """Build the orientation grid only when grid_signature() has changed.

        Returns True if the grid was (re)built. The signature is stored after
        the build, not before, because compute_weights writes Exp.Field back
        when Exp.gField is set -- capturing it first would cost a spurious
        rebuild on the very next call.
        """
        if (self.gridSig is not None and self.phi is not None
                and self.grid_signature() == self.gridSig):
            return False

        if self.Opt.OriSelType[0] in ['g_eff', 'g_eff+HFC']:
            self.make_grid(grid=self.grid_name(), nKnots=self.Opt.nKnots)
            self.compute_weights()
        elif self.Opt.OriSelType[0] == 'brute force':
            self.computeOrisel_eig(grid=self.grid_name(), nKnots=self.Opt.nKnots,
                                   epsilon=0.33, returnMatrix=False, setActive=True)
        self.gridSig = self.grid_signature()
        return True

    def reRun(self):
        if self.verbose: print('----start calculaitons ----');
        sta = time.perf_counter()
        if self.Opt.OriSelType[0] == 'precalculated':
            if type(self.Opt.OriSelInp) is not np.ndarray:
                self.errorFunc(f'HYSCORE_sim_reRun: expected a non-empty Opt.OriSelInp, got None. \n - Opt.OriSelType={self.Opt.OriSelType}')
                return
            if len(self.Opt.OriSelInp.shape)!=2:
                self.errorFunc(f'HYSCORE_sim_reRun: expected a 3-column Opt.OriSelInp, got shape={self.Opt.OriSelInp.shape}. \n - Opt.OriSelType={self.Opt.OriSelType}')
                return
            if self.Opt.OriSelInp.shape[1]!=3:
                self.errorFunc(f'HYSCORE_sim_reRun: expected a 3-column Opt.OriSelInp, got shape={self.Opt.OriSelInp.shape}. \n - Opt.OriSelType={self.Opt.OriSelType}')
                return
            self.phi = self.Opt.OriSelInp[:, 0]
            self.theta = self.Opt.OriSelInp[:, 1]
            self.ak = self.Opt.OriSelInp[:, 2]
            self.tri = None
            self.triareas = None
            self.gridSig = None       # a loaded grid is not one we can re-derive
        else:
            # see the GRID CACHE banner above: this is a no-op unless something
            # the grid depends on has actually changed
            self.ensure_grid()


        #elapsed = time.perf_counter() - sta
        #if self.verbose: print(f"HYSCORE grid in {elapsed*1000:.2f} ms")
        
        
        #elapsed2 = time.perf_counter() - elapsed
        #if self.verbose: print(f"HYSCORE spin_op in {elapsed2*1000:.2f} ms")
        
        self.rawSpec = np.zeros((self.Exp.nPoints, self.Exp.nPoints), dtype=complex)
        self.X = np.linspace(-self.Exp.MaxFreq, self.Exp.MaxFreq, self.Exp.nPoints, endpoint=False)
        self.Y = np.linspace(-self.Exp.MaxFreq, self.Exp.MaxFreq, self.Exp.nPoints, endpoint=False)
        
        self.compute_transitions()
        self.makeFFT()
        end = time.perf_counter()
        #rel =  end - elapsed2
        tot = end-sta
        if self.verbose: 
            #print(f"HYSCORE compute_trans in {rel*1000:.2f} ms")
            print(f"HYSCORE total in {tot*1000:.2f} ms")
    def run(self):
        self.preCompute()
        self.reRun()
    def makeFFT(self):
        if self.Sys.lw:
            ispec = np.fft.ifft2(np.fft.ifftshift(self.rawSpec))
            X, Y = np.meshgrid(self.X, self.Y)
            sig = self.Sys.lw[0]/2.354820045 # 2.354820045 = 2sqrt(2ln2)


            GA = np.exp(-0.5*X**2/sig**2)*np.exp(-0.5*Y**2/sig**2)
            fGA = np.fft.ifft2(np.fft.ifftshift(GA))
            
            if self.Opt.KillNeg:
                np2 = int(self.Exp.nPoints/2)
                fGA[:, np2:] = 0.0
                fGA[np2:, :] = 0.0
            self.Spectrum=np.abs(np.fft.fftshift(np.fft.fft2(ispec*fGA)))
            
        else:
            self.Spectrum = np.abs(self.rawSpec)
        
    def make_grid(self, grid='fibonacci', nKnots=20,
                  Symmetry=None, phiDeg=None, thetaDeg=None, weights=None,
                  epsilon=0.33):
        # 1.  Exp contains both 'phi' and 'theta'
        phi=None
        theta=None
        ww=None
        tri=None
        triareas=None
        if grid=='input':
            if phiDeg and thetaDeg:
                phiT   = np.asarray(phiDeg)*np.pi/180 # inputs are assumed in MHz
                thetaT = np.asarray(thetaDeg)*np.pi/180
                phi, theta = np.meshgrid(phiT, thetaT, indexing='ij')
                ww = np.ones_like(phi).ravel()
            else:
                self.errorFunc('either phi or theta are not provided')
        elif grid=='fibonacci':
            n   = int(nKnots**1.5)
            thLim=np.pi
            phLim=2*np.pi
            if Symmetry:
                # Symmetry-restricted azimuthal ranges are not implemented here,
                # so the full azimuth is used regardless. Bound unconditionally
                # above because errorFunc does not necessarily abort (the GUI
                # handler only shows a dialog), and the mask below needs phLim.
                # Use grid='sphgrid' for a symmetry-aware (Ci) EasySpin grid.
                self.errorFunc(f"make_grid: Symmetry={Symmetry!r} is not "
                               "implemented for grid='fibonacci'; "
                               "covering the full azimuth")
            # Golden angle
            golden_angle = np.pi * (3.0 - np.sqrt(5.0))   # ≈ 2.39996322973
            idx = np.arange(0, n, dtype='float64')
            cos_theta = (idx + epsilon)/ (float(n - 1)+2*epsilon)
            ttheta = np.arccos(cos_theta)   # polar angle from +z axis
            tphi = np.remainder(golden_angle*idx, 2*np.pi)
            mask = (ttheta <= thLim) & (tphi <= phLim)
            
            phi = tphi[mask]
            theta = ttheta[mask]
            # Equal-area grid, so the knots carry equal weight; normalised by
            # the sum so that the intensity scale does not depend on nKnots and
            # is directly comparable with the 'sphgrid' branch below.
            ww = np.ones_like(phi)
            ww = ww/np.sum(ww)
            tri, triareas = self.grid_triangulation(phi, theta)
        elif grid=='spiral':
            # Curtesy of Ed Reijerse /MAGRES.. currently not working well
            # Rescaled like the 'sphgrid' branch so nKnots means a comparable
            # cost in every method. This grid holds about
            #   sum_i floor(sin(i*pi/K)*4K) ~ 4K*cot(pi/2K) ~ 8K^2/pi
            # knots, so invert that against what 'fibonacci' would produce.
            krid   = max(2, int(round(np.sqrt(int(nKnots**1.5)*np.pi/8))))
            itheta, iphi = krid, 0
            step = np.pi/krid
            phi_l, theta_l = [], []
            thetaa=0.0
            dthe=0.0
            dphi=0.0
            nphi=0            # 0 so the first pass opens a ring straight away
            while True:
                if itheta <= 0:
                    self.errorFunc('pdrpeal too many calls')
                    break                        # errorFunc may not abort
                # A ring is full after its own nphi knots. Testing against
                # 2*krid instead put every ring's surplus knots on top of each
                # other -- with nKnots=20 the first 40 knots were all the north
                # pole, which is what made this grid "not work well".
                if iphi == nphi:                 # new latitude ring
                    iphi   = 0
                    itheta -= 1
                    thetam  = itheta * step
                    nphi    = int(np.floor(np.sin(thetam)*(4*krid)))
                    if nphi <= 0:
                        self.errorFunc('BUG pdrpeal 1')
                    thetaa = thetam + 0.5*step
                    dthe   = -step/nphi
                    # radians: theta is built in radians from step = pi/krid, so
                    # a degree increment here left phi on a different scale
                    dphi   = 2*np.pi/nphi
                thetap = thetaa + iphi*dthe
                phip   = iphi*dphi
                phi_l.append(phip)
                theta_l.append(thetap)
                iphi += 1
                if iphi == nphi and itheta == 1:
                    break
            phi = np.remainder(np.array(phi_l), 2*np.pi)
            theta = np.array(theta_l)
            # Rings sit at constant steps in theta but hold nphi ~ sin(theta)
            # knots each, so a knot covers (2*pi*sin(th)*dth)/nphi = const solid
            # angle: the grid is equal-area and the weights stay uniform.
            # Normalised by the sum, as in the other branches.
            ww = np.ones_like(theta)
            ww = ww/np.sum(ww)
            tri, triareas = self.grid_triangulation(phi, theta)

        elif grid=='sphgrid':
            # SOPHE / triangular grid, Ci branch (upper hemisphere, four
            # octants, open phi). Port of sphgrid_ in EasySpin's sphgrid.m.
            # The grid is Eq. (15) of Stoll & Schweiger, J.Magn.Reson. 178
            # (2006) 42-55, p.47 -- which prints only the one-octant case --
            # generalised by sphgrid.m to
            #    theta_k = (pi/2)*k/(K-1), nPhi = nOct*k+1, phi = 0..maxPhi
            # Original grid: Wang & Hanson, J.Magn.Reson. A 117, 1-8 (1995).
            #
            # nKnots is rescaled rather than passed through: EasySpin's GridSize
            # counts knots along a quarter meridian, so GridSize=20 would give
            # 761 knots against the 89 that 'fibonacci' produces for nKnots=20.
            # Inverting the knot count 1+2*K*(K-1) keeps the two branches
            # interchangeable at equal cost (one diagonalisation per knot).
            nTarget = int(nKnots**1.5)
            K = max(2, int(round((1.0+np.sqrt(max(2.0*nTarget-1.0, 1.0)))/2.0)))
            nOct = 4
            maxPhi = 2*np.pi
            dtheta = (np.pi/2)/(K-1)

            n = K + nOct*K*(K-1)//2
            phi = np.zeros(n)
            theta = np.zeros(n)
            ww = np.zeros(n)
            sindth2 = np.sin(dtheta/2)
            w0 = 1.0                          # open phi -> full weight at phi=0

            # North pole. Handled outside the loop, which is also why Eq. (15)'s
            # phi = (pi/2)*(q/k) never has to evaluate 0/0.
            ww[0] = maxPhi*(1.0-np.cos(dtheta/2))

            start = 1
            for iSlice in range(2, K):        # every slice but the equatorial one
                nPhi = nOct*(iSlice-1)+1
                dPhi = maxPhi/(nPhi-1)
                idx = start+np.arange(nPhi)
                ww[idx] = (2*np.sin((iSlice-1)*dtheta)*sindth2*dPhi
                           * np.r_[w0, np.ones(nPhi-2), 0.5])
                phi[idx] = np.linspace(0.0, maxPhi, nPhi)
                theta[idx] = (iSlice-1)*dtheta
                start += nPhi

            # Equatorial slice. The missing 2*sin(theta) factor is deliberate:
            # on a hemisphere the equator spans only half a band.
            nPhi = nOct*(K-1)+1
            dPhi = maxPhi/(nPhi-1)
            idx = start+np.arange(nPhi)
            phi[idx] = np.linspace(0.0, maxPhi, nPhi)
            theta[idx] = np.pi/2
            ww[idx] = sindth2*dPhi*np.r_[w0, np.ones(nPhi-2), 0.5]

            # Border removal: drop the phi=2*pi duplicate ending every ring
            keep = np.ones(n, dtype=bool)
            keep[np.cumsum(nOct*np.arange(1, K)+1)] = False
            phi = phi[keep]
            theta = theta[keep]
            ww = 2*(2*np.pi/maxPhi)*ww[keep]              # EasySpin: sum = 4*pi

            # This grid is NOT equal-area (nPhi grows as k, i.e. as theta,
            # whereas equal-area sampling needs sin(theta)), so the weights are
            # essential -- uniform weights would bias the powder average toward
            # the poles. Normalised by the sum as for 'fibonacci'; that is a
            # pure global scale, so relative weights stay exact vs EasySpin.
            ww = ww/np.sum(ww)
            tri, triareas = self.grid_triangulation(phi, theta)
        # 3.  Default case – latitude rings
        else:
            nknots = int(nKnots)
            t_theta = np.linspace(0, np.pi, nknots)
    
            phi_list   = []
            theta_list = []
    
            for th in t_theta:
                num = int(np.floor(nknots * abs(np.sin(th))))
                if num == 0:
                    continue
                t_phi = np.linspace(0, 2*np.pi, num+1)[:-1]   # exclude last point
                phi_list.extend(t_phi)
                theta_list.extend([th]*num)
    
            phi   = np.array(phi_list)
            theta = np.array(theta_list)
            ww    = np.ones_like(phi)
        self.phi = phi
        self.theta = theta
        self.ww = ww
        # TODO(pyHYSCORE): 'fibonacci' and 'sphgrid' now both publish a knot
        # triangulation here, so the orientation-selection plot in pyHYSCORE can
        # consume self.tri/self.triareas instead of building its own Delaunay.
        # Not wired up yet -- ask before modifying pyHYSCORE.py.
        self.tri = tri
        self.triareas = triareas

    def grid_triangulation(self, phi, theta):
        """Delaunay triangulation of an orientation grid, with exact solid angles.

        Hemisphere grids follow the nOctants=4 branch of sphgrid.m: triangulate
        in the azimuthal-equidistant projection (r = theta), which is a
        bijection there and keeps triangles well conditioned out to the equator.
        Full-sphere grids fall back to the 3D convex hull instead, since that
        projection is no longer injective past the equator.
        Returns (tri, areas), areas summing to 4*pi.
        """
        phi = np.asarray(phi, dtype=float).ravel()
        theta = np.asarray(theta, dtype=float).ravel()
        if phi.size < 4:
            return np.zeros((0, 3), dtype=int), np.zeros(0)

        st = np.sin(theta)
        vecs = np.vstack((st*np.cos(phi), st*np.sin(phi), np.cos(theta)))

        if theta.max() <= np.pi/2 + 1e-9:
            # Hemisphere ('fibonacci', 'sphgrid'): the azimuthal-equidistant
            # projection is a bijection there, so a planar Delaunay suffices.
            tri = Delaunay(np.column_stack((theta*np.cos(phi),
                                            theta*np.sin(phi)))).simplices
        else:
            # Full sphere ('spiral'): that projection folds the two hemispheres
            # onto each other, so triangulate on the sphere itself. For points
            # on a sphere the convex hull IS the spherical Delaunay.
            tri = ConvexHull(vecs.T).simplices

        # Triangle solid angles by l'Huilier's formula, as triangleareas() in
        # sphgrid.m: edge arc lengths, then the spherical excess.
        x1 = vecs[:, tri[:, 0]]
        x2 = vecs[:, tri[:, 1]]
        x3 = vecs[:, tri[:, 2]]
        a1 = np.arccos(np.clip(np.sum(x2*x3, axis=0), -1.0, 1.0))
        a2 = np.arccos(np.clip(np.sum(x3*x1, axis=0), -1.0, 1.0))
        a3 = np.arccos(np.clip(np.sum(x1*x2, axis=0), -1.0, 1.0))
        s = (a1+a2+a3)/2                      # half the perimeter
        t = np.tan(s/2)*np.tan((s-a1)/2)*np.tan((s-a2)/2)*np.tan((s-a3)/2)
        # MATLAB takes real(sqrt(t)); for t<0 that is 0, which np.maximum matches
        areas = 4*np.arctan(np.sqrt(np.maximum(t, 0.0)))

        # Discard slivers thrown off by the Delaunay call, as sphgrid.m does
        if areas.size:
            keep = areas >= areas.mean()*1e-3
            tri, areas = tri[keep], areas[keep]

        tri = np.sort(tri, axis=1)
        order = np.lexsort((tri[:, 2], tri[:, 1], tri[:, 0]))
        tri, areas = tri[order], areas[order]
        total = areas.sum()
        if total > 0:
            areas = areas/total*4*np.pi
        return tri, areas

    def compute_weights(self):
        if (self.Exp.ExciteWidth > 0) and (self.Exp.mwFreq > 0):
            ctheta = np.cos(self.theta)
            stheta = np.sin(self.theta)
            cphi   = np.cos(self.phi)
            sphi   = np.sin(self.phi)
            nangles = len(cphi)
            angle = np.column_stack((stheta * cphi,
                             stheta * sphi,
                             ctheta))          # shape (n_angles, 3)
            self.lxyz = angle
            # ---------- effective g‑tensor ---------------------------------------- #
            Gx = angle[:, 0] * self.Sys.g[0][0]
            Gy = angle[:, 1] * self.Sys.g[0][1]
            Gz = angle[:, 2] * self.Sys.g[0][2]
            self.geff = np.sqrt(Gx**2 + Gy**2 + Gz**2)
            
    
            #Exp.ExciteWidth = safeget(Exp, 'ExciteWidth', 0.0)
            #Exp.mwFreq      = safeget(Exp, 'mwFreq', 0.0)   # in GHz
            
            # If the experiment is defined via gField, compute the field
            if self.Exp.gField:
                if self.Exp.mwFreq == 0:
                    self.errorFunc('HYSCORE_compute_weights: Exp.mwFreq is not set to calculate B0')
                # B0 in milliTesla
                self.Exp.Field = self.Exp.mwFreq /self.Exp.gField *self.h_ov_beta_mTMHz*1e3
            if not self.Exp.Field:
                self.errorFunc('HYSCORE_compute_weights: Exp.Field is not set')
            
            #  Excitation profile (ak)

            # Effective frequency (MHz)
            feff = self.Exp.Field * self.geff / self.h_ov_beta_mTMHz   # shape (nangles,)
            sig = self.Exp.ExciteWidth/2.35482
            if self.Opt.OriSelType[0]=='g_eff':
                # Simple Gaussian line shape
                ak = np.exp(-0.5 * ((feff - self.Exp.mwFreq*1e3) / sig)**2)
            elif self.Opt.OriSelType[0]=='g_eff+HFC':
                # ---------- optional HF selection ------------------------------------- #
                # NEEDS CHECKING.
                
                total_nuclei = int(np.sum(self.Sys.nNucs))
                projA = np.zeros((nangles, total_nuclei), dtype=float)
                k = 0  # column counter in projA
                mult = []
                for ii in range(len(self.Sys.nNucs)):
                    for jj in range(int(self.Sys.nNucs[ii])):
                        k += 1
                        # Hyperfine tensor for nucleus ii
                        if np.sum(self.Sys.Apa[ii]) != 0:
                            At = mf.RotTransDiag(self.Sys.A[ii], self.Sys.Apa[ii]*np.pi/180.0)
                            
                        else:
                            At = np.diag(self.Sys.A[ii])
                
                        # Projected components (vectorised)
                        LA = angle @ At          # shape (nangles, 3)
                        projA[:, k-1] = np.linalg.norm(LA, axis=1)  # sqrt(sum of squares)
                mult.append(int(2 * self.Sys.I[ii] + 1))

                # Build the MIs matrix – all possible nuclear spin projections
                mult_arr = np.array(mult, dtype=int)          # multiplicities
                nTrans    = int(np.prod(mult_arr))           # total number of hyperfine transitions
        
                # Generate all combinations of M_I values
                # Example: for I=3/2 (mult=4) -> [-1.5, -0.5, 0.5, 1.5]
                vals = [np.arange(-(m-1)/2, (m-1)/2 + 1, 1) for m in mult_arr]
                grids = np.meshgrid(*vals, indexing='ij')
                MIs = np.stack([g.ravel() for g in grids], axis=1)  # shape (nTrans, len(mult))
        
                ak = np.zeros(nangles, dtype=float)
        
                for ii in range(nTrans):
                    # Hyperfine shift for this transition
                    ffre = np.sum(projA[:, :len(mult)] * MIs[ii, :], axis=1)
                    ak += np.exp(-0.5 * ((feff + ffre - self.Exp.mwFreq * 1e3) /sig)**2)
        
        else:
            ak = np.ones_like(self.theta)

        if np.max(ak) < 1e-4:
            self.errorFunc('HYSCORE_compute_weights: no resonances')
        self.ak = ak*self.ww
        
    def build_spin_operators(self):
        dims = [int(2 * s + 1) for s in self.Sys.I]
        max_dim = max(dims)

        total_nuclei = int(np.sum(self.Sys.nNucs))
        ndiffnucs= len(self.Sys.I)
        # --- containers --------------------------------------------------------
        # start with a clean slate for global variables
        self.sIx = [None] * (max_dim + 1)      # index = dimension (2*I+1)
        self.sIy = [None] * (max_dim + 1)
        self.sIz = [None] * (max_dim + 1)
        
        self.sIxx = [None] * (max_dim + 1)
        self.sIxy = [None] * (max_dim + 1)
        self.sIxz = [None] * (max_dim + 1)

        self.sIyx = [None] * (max_dim + 1)
        self.sIyy = [None] * (max_dim + 1)
        self.sIyz = [None] * (max_dim + 1)

        self.sIzx = [None] * (max_dim + 1)
        self.sIzy = [None] * (max_dim + 1)
        self.sIzz = [None] * (max_dim + 1)
    
        self.HFIx = []
        self.HFIy = []
        self.HFIz = []
        
        self.HamA = []
        self.HamB = []
        
    
        self.Quad = []
    
        self.cii = [None] * (max_dim + 1)
        self.ckk = [None] * (max_dim + 1)
        self.cll = [None] * (max_dim + 1)
        self.cnn = [None] * (max_dim + 1)
        self.cjj = [None] * (max_dim + 1)
        self.cmm = [None] * (max_dim + 1)
    

        
        if len(self.Sys.A)<len(self.Sys.I):
            self.errorFunc('HYSCORE: len(Sys.A) is different than len(Sys.I)')
        # -----------------------------------------------------------------------
        #  loop over nuclei. There is no need for multiple entries 
        #  for the same I, so we will sort record them into lists by multiplicity
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv    
        for ii in range(ndiffnucs):
            # Spin quantum number of this nucleus (MATLAB index is 1‑based)
            if self.Sys.useFor[ii]=='ori.sel.only': continue

            s = self.Sys.I[ii]                     # Sys.I is expected to be a list/array
            mult = int(2 * s + 1)              # dimension of the local space
            if not self.sIx[mult]:
                # 1.  Spin operators ------------------------------------------------
                self.sIx[mult] = mf.SpinOp([s], ['x'])
                self.sIy[mult] = mf.SpinOp([s], ['y'])
                self.sIz[mult] = mf.SpinOp([s], ['z'])
        
                ### what I need to do:
                # idx = np.arange(mult, dtype=int)
                # The 6 arrays are simply the broadcasted index vectors:
                # cnt = 0
                # tii=np.zeros(mult**6)
                # tkk=np.zeros(mult**6)
                # tll=np.zeros(mult**6)
                # tnn=np.zeros(mult**6)
                # tjj=np.zeros(mult**6)
                # tmm=np.zeros(mult**6)
                # for ii in range(mult):
                #     for kk in range(mult):
                #         for ll in range(mult):
                #             for nn in range(mult):
                #                 for jj in range(mult):
                #                     for mm in range(mult):
                #                         tii[cnt]=ii
                #                         tkk[cnt]=kk
                #                         tll[cnt]=ll
                #                         tnn[cnt]=nn
                #                         tjj[cnt]=jj
                #                         tmm[cnt]=mm
                #                         cnt+=1
                ### That is highly inefficient way
                ### or thanks to ChatGPT:
                
                flat = np.arange(mult ** 6 , dtype=np.int32)    # 0, 1, 2, … , size-1
                # (tii, tkk, tll, tnn, tjj, tmm) are each shape (size,)
                tii, tkk, tll, tnn, tjj, tmm = np.unravel_index(flat, (mult,)*6)
                
                self.cii[mult]=tii
                self.ckk[mult]=tkk
                self.cll[mult]=tll
                self.cnn[mult]=tnn
                self.cjj[mult]=tjj
                self.cmm[mult]=tmm
        
            # 2.  Hyperfine tensor A

            if self.Sys.A[ii].ndim == 1:
                if not (self.Sys.A[ii].shape[0]==3):
                    self.errorFunc(f'Sys.A[{ii}] has an unexpected shape')
                A = mf.RotTransDiag(self.Sys.A[ii], self.Sys.Apa[ii]*np.pi/180.0)
            elif self.Sys.A[ii].ndim == 2:
                if not (self.Sys.A[ii].shape[0]==3)&(self.Sys.A[ii].shape[1]==3):
                    self.errorFunc(f'Sys.A[{ii}] has an unexpected shape')
                A = mf.RotTrans3x3(self.Sys.A[ii], self.Sys.Apa[ii]*np.pi/180.0)
            else:
                self.errorFunc('Sys.A has an unexpected shape')
    
            # Build the three components of the hyperfine operator
            self.HFIx.append( self.sIx[mult] * A[0, 0] + self.sIy[mult] * A[1, 0] + self.sIz[mult] * A[2, 0] )
            self.HFIy.append( self.sIx[mult] * A[0, 1] + self.sIy[mult] * A[1, 1] + self.sIz[mult] * A[2, 1] )
            self.HFIz.append( self.sIx[mult] * A[0, 2] + self.sIy[mult] * A[1, 2] + self.sIz[mult] * A[2, 2] )
    
            # 3.  Quadrupole coupling
            
            quad_mat = np.zeros((mult, mult), dtype=complex)   # default 0
            if self.Sys.Q:
                if (mult>2)&(np.sum(np.abs(self.Sys.Q[ii]))!=0):
                    if not self.sIxx[mult]:
                        # products
                        self.sIxx[mult] = self.sIx[mult] @ self.sIx[mult]
                        self.sIxy[mult] = self.sIx[mult] @ self.sIy[mult]
                        self.sIxz[mult] = self.sIx[mult] @ self.sIz[mult]
                
                        self.sIyx[mult] = self.sIy[mult] @ self.sIx[mult]
                        self.sIyy[mult] = self.sIy[mult] @ self.sIy[mult]
                        self.sIyz[mult] = self.sIy[mult] @ self.sIz[mult]
                
                        self.sIzx[mult] = self.sIz[mult] @ self.sIx[mult]
                        self.sIzy[mult] = self.sIz[mult] @ self.sIy[mult]
                        self.sIzz[mult] = self.sIz[mult] @ self.sIz[mult]
                        
                    if self.Sys.Q[ii].ndim == 1:
                        if not (self.Sys.Q[ii].shape[0]==3):
                            self.errorFunc(f'Sys.Q[{ii}] has an unexpected shape. Need three values')
                        Q = mf.RotTransDiag(self.Sys.Q[ii], self.Sys.Qpa[ii]*np.pi/180.0)
                    elif self.Sys.Q[ii].ndim == 2:
                        if not (self.Sys.Q[ii].shape[0]==3)&(self.Sys.Q[ii].shape[1]==3):
                            self.errorFunc(f'Sys.Q[{ii}] has an unexpected shape. Expected 3x3')
                        Q = mf.RotTrans3x3(self.Sys.Q[ii], self.Sys.Qpa[ii].np.pi/180.0)
                    else:
                        self.errorFunc('Sys.Q has an unexpected shape')
                        
                    quad_mat = (
                        Q[0, 0] * self.sIxx[mult] + Q[0, 1] * self.sIxy[mult] + Q[0, 2] * self.sIxz[mult] +
                        Q[1, 0] * self.sIyx[mult] + Q[1, 1] * self.sIyy[mult] + Q[1, 2] * self.sIyz[mult] +
                        Q[2, 0] * self.sIzx[mult] + Q[2, 1] * self.sIzy[mult] + Q[2, 2] * self.sIzz[mult]
                    )
            self.Quad.append(quad_mat)
            
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #  loop over nuclei
    def local_sham(self, n_ang, n_nuc, B0):
        costheta = np.cos(self.theta[n_ang])
        sintheta = np.sin(self.theta[n_ang])

        cosphi = np.cos(self.phi[n_ang])
        sinphi = np.sin(self.phi[n_ang])
     
        lx = sintheta*cosphi
        ly = sintheta*sinphi
        lz = costheta
    
        nu_n = self.Sys.gamma[n_nuc]*B0
        mult = int(2*self.Sys.I[n_nuc]+1)
        
        Zeeman_A = -nu_n*(lx*self.sIx[mult] + ly*self.sIy[mult] + lz*self.sIz[mult])
        
        Gx=lx*self.Sys.g[0][0]
        Gy=ly*self.Sys.g[0][1]
        Gz=lz*self.Sys.g[0][2]
        geff=np.sqrt(Gx**2+Gy**2+Gz**2)

        HFI_A = self.HFIx[n_nuc]*Gx/geff*1/2 + self.HFIy[n_nuc]*Gy/geff*1/2 + self.HFIz[n_nuc]*Gz/geff*1/2
        
        Ham_A = Zeeman_A - HFI_A + self.Quad[n_nuc]
        Ham_B = Zeeman_A + HFI_A + self.Quad[n_nuc]
        return Ham_A, Ham_B
    
    def compute_transitions(self):
        ### Do all the checks here. TBD
        for n_ang in range(len(self.phi)):

            Amp = None
            W = None
            for n_nuc in range(len(self.Sys.I)):
                if self.Sys.useFor[n_nuc]=='ori.sel.only': continue
                mult = int(2*self.Sys.I[n_nuc]+1)
                nidx = len(self.cii[mult])
                
                tAmp = np.zeros((nidx, 2), dtype=complex)
                tW = np.zeros((nidx, 2), dtype=float)

                Ham_A, Ham_B = self.local_sham(n_ang, n_nuc, self.Exp.Field)
                # For Hermitian matrix
                DA, VA = np.linalg.eigh(Ham_A)       # eigenvalues + eigenvectors
                DB, VB = np.linalg.eigh(Ham_B)       # eigenvalues + eigenvectors
                # np.linalg.eigh returns V whose columns are eigenvectors, so VA[:, i] = |α_i> and VB[:, l] = |β_l>
                M = VA.conj().T@VB  # check if needed a transpose !
                # cM=np.linalg.inv(M)
                cM = M.conj().T # Because VA and VB are unitary (eigenvectors of Hermitian Hamiltonians from eigh), their product M = V_A^H V_B is also unitary. So
                
                wa = np.real(DA[:, None] - DA[None, :])   # shape (n, n)
                wb = np.real(DB[:, None] - DB[None, :])   # shape (n, n)
                
                ii, ll, jj, nn, kk, mm = (self.cii[mult],
                                          self.cll[mult],
                                          self.cjj[mult],
                                          self.cnn[mult],
                                          self.ckk[mult],
                                          self.cmm[mult])
                                         
                ampar = np.array(
                         M[ii, ll]*
                         cM[ll, jj]*
                          M[jj, nn]*
                         cM[nn, kk]*
                          M[kk, mm]*
                         cM[mm, ii]
                         )
                
                # Blind sponts
                exppar =np.array(np.exp(-2j * np.pi * (wa[ii, jj] + wb[ll, mm]) * self.Exp.tau*1e-3) +
                          np.exp(+2j * np.pi * (wa[kk, jj] + wb[nn, mm]) * self.Exp.tau*1e-3))*self.ak[n_ang]
                if self.Exp.tDead:
                    phasedr = np.exp(-2j * np.pi * (wa[ii, kk] + wb[ll, nn]) * self.Exp.tDead)
                    tAmp[:, 0] = ampar*exppar*phasedr
                    tAmp[:, 1] = np.conj(ampar)*exppar*phasedr                    
                else:
                    tAmp[:, 0] = ampar*exppar
                    tAmp[:, 1] = np.conj(ampar)*exppar
        
                tW[:, 0]=wa[ii, kk]
                tW[:, 1]=wb[ll, nn]
                
                if not Amp:
                    Amp = np.array(tAmp, dtype=complex)
                else:
                    Amp = np.concatenate((Amp, tAmp), axis=0)
                    
                if not W:
                    W = np.array(tW, dtype=float)
                else:
                    W = np.concatenate((W, tW), axis=0)
                    
            # if self.Opt.ProdRule:
            #     print('tbd')
            
            self.bin_hyscore(W[:, 0], W[:, 1], Amp[:, 0], self.Exp.MaxFreq)
            self.bin_hyscore(W[:, 1], W[:, 0], Amp[:, 1], self.Exp.MaxFreq)
            
  # def bin_hyscore(self, omega_a, omega_b, amp, max_freq):
  #
  #     # ---------- sanity checks ----------
  #     n_points = self.Exp.nPoints
  #     
  #     # Ensure all inputs are NumPy arrays of the right dtype
  #     #omega_a = np.asarray(omega_a, dtype=np.float64)
  #     #omega_b = np.asarray(omega_b, dtype=np.float64)
  #     # ---------- compute mapping ----------
  #     n_points_d = float(n_points)
  #     dx = (n_points_d - 1.0) / (2.0 * max_freq)
  # 
  #     # Convert frequencies to integer indices
  #     idx1 = np.floor((omega_a + max_freq) * dx).astype(np.int64)
  #     idx2 = np.floor((omega_b + max_freq) * dx).astype(np.int64)
  # 
  #     # Keep only peaks that fall inside the array bounds
  #     valid = (idx1 >= 0) & (idx1 < n_points) & (idx2 >= 0) & (idx2 < n_points) & ( ((omega_a==0.0)|(omega_b==0.0))==False)
  #     if not np.any(valid):
  #         return  # nothing to add
  # 
  #     idx1 = idx1[valid]
  #     idx2 = idx2[valid]
  #     amp_valid = amp[valid]
  #     
  #     #self.rawSpec[idx1, idx2]=amp_valid
  #     #self.rawSpec[idx2, idx1]=amp_valid
  #     
  #     # Linear index in column‑major order (as MATLAB does)
  #     linear_idx = idx1 + idx2 * n_points
  #     # Symmetric counterpart
  #     sym_linear_idx = idx2 + idx1 * n_points
  # 
  #     # ---------- add contributions ----------
  #     np.add.at(self.rawSpec.ravel(), linear_idx, amp_valid)
  #     np.add.at(self.rawSpec.ravel(), sym_linear_idx, amp_valid)
    def bin_hyscore(self, omega_a, omega_b, amp, max_freq):

        # ---------- sanity checks ----------
        n_points = self.Exp.nPoints
        
        # Ensure all inputs are NumPy arrays of the right dtype
        #omega_a = np.asarray(omega_a, dtype=np.float64)
        #omega_b = np.asarray(omega_b, dtype=np.float64)
        # ---------- compute mapping ----------
        n_points_d = float(n_points)
        dx = n_points_d / (2.0 * max_freq)   # was (n_points_d - 1.0) / (2*max_freq)
    
        # Convert frequencies to integer indices ## in the future seek ways to make subbin.
        idx1 = np.round((omega_a + max_freq) * dx).astype(np.int64)
        idx2 = np.round((omega_b + max_freq) * dx).astype(np.int64)
    
        # Keep only peaks that fall inside the array bounds
        valid = (idx1 >= 0) & ((idx1+1) < n_points) & (idx2 >= 0) & ((idx2+1) < n_points) & ( ((omega_a==0.0)|(omega_b==0.0))==False)
        if not np.any(valid):
            return  # nothing to add
    
        idx1 = idx1[valid]
        idx2 = idx2[valid]
        amp_valid = amp[valid]
        
        #self.rawSpec[idx1, idx2]=amp_valid
        #self.rawSpec[idx2, idx1]=amp_valid
        
        # Linear index in column‑major order (as MATLAB does)
        linear_idx = idx1 + idx2 * n_points
        # Symmetric counterpart
        sym_linear_idx = idx2 + idx1 * n_points
    
        # ---------- add contributions ----------
        np.add.at(self.rawSpec.ravel(), linear_idx, amp_valid)
        np.add.at(self.rawSpec.ravel(), sym_linear_idx, amp_valid)        
    def computeOrisel_eig(self, grid='fibonacci', nKnots=20, epsilon=0.33, returnMatrix=False, setActive=False):
        # Brut force calculation of orientation selection using diagonalization of a complete spin Hamiltonian. For complex cases
        self.preCompute()

        self.make_grid(grid=grid, nKnots=nKnots, epsilon=epsilon)

        OriSelOut = np.zeros((len(self.phi), 3))
        OriSelOut[:, 0] =self.phi
        OriSelOut[:, 1] =self.theta

        if self.Exp.gField:
            if self.Exp.mwFreq == 0:
                self.errorFunc('computeOrisel_eig: Exp.mwFreq is not set to calculate B0')
            # B0 in milliTesla
            self.Exp.Field = self.Exp.mwFreq /self.Exp.gField *self.h_ov_beta_mTMHz*1e3
        if not self.Exp.Field:
            self.errorFunc('computeOrisel_eig: Exp.Field is not set')

        if not ((self.Exp.ExciteWidth > 0) and (self.Exp.mwFreq > 0)):
            self.errorFunc('computeOrisel_eig: Exp.ExciteWidth and Exp.mwFreq are not set')

        sig = self.Exp.ExciteWidth/2.35482
        mwFreqMHz = self.Exp.mwFreq * 1e3

        if self.Sys.S[0] != 0.5:
            self.errorFunc('computeOrisel_eig: for now only S=1/2 is allowed')

        smult = int(2*self.Sys.S[0]+1)

        stee = ['e']
        Ssys = [self.Sys.S[0]]
        HFc = []
        #Qc = [None]
        gamma  = []
        mult = []
        for n_nuc in range(len(self.Sys.nNucs)):
            if self.Sys.useFor[n_nuc] == 'sim.only': continue
            for ii in range(int(self.Sys.nNucs[n_nuc])):
                if self.Sys.A[n_nuc].ndim == 1:
                    if not (self.Sys.A[n_nuc].shape[0]==3):
                        self.errorFunc(f'Sys.A[{n_nuc}] has an unexpected shape')
                    A = mf.RotTransDiag(self.Sys.A[n_nuc], self.Sys.Apa[n_nuc]*np.pi/180.0)
                elif self.Sys.A[n_nuc].ndim == 2:
                    if not (self.Sys.A[n_nuc].shape[0]==3)&(self.Sys.A[n_nuc].shape[1]==3):
                        self.errorFunc(f'Sys.A[{n_nuc}] has an unexpected shape')
                    A = mf.RotTrans3x3(self.Sys.A[n_nuc], self.Sys.Apa[n_nuc]*np.pi/180.0)
                else:
                    self.errorFunc('Sys.A has an unexpected shape')
                #### quadrupole coupling is not needed for orientation selection... at least in most cases. Revisit at some point
                # if self.Sys.Q and self.Sys.I[n_nuc]>0.5:

                #     if self.Sys.Q[ii].ndim == 1:
                #         if not (self.Sys.Q[ii].shape[0]==3):
                #             raise ValueError(f'Sys.Q[{ii}] has an unexpected shape. Need three values')
                #         Q = mf.RotTransDiag(self.Sys.Q[ii], self.Sys.Qpa[ii]*np.pi/180.0)
                #     elif self.Sys.Q[ii].ndim == 2:
                #         if not (self.Sys.Q[ii].shape[0]==3)&(self.Sys.Q[ii].shape[1]==3):
                #             raise ValueError(f'Sys.Q[{ii}] has an unexpected shape. Expected 3x3')
                #         Q = mf.RotTrans3x3(self.Sys.Q[ii], self.Sys.Qpa[ii].np.pi/180.0)
                #     else:
                #         raise ValueError('Sys.Q has an unexpected shape')
                # else:
                #
                stee.append('e')
                Ssys.append(self.Sys.I[n_nuc])
                HFc.append(A)
                gamma.append(self.Sys.gamma[n_nuc])
                mult.append(int(2*self.Sys.I[n_nuc]+1))



        stee[0]='x'
        sSx = mf.SpinOp(Ssys, stee)
        stee[0]='y'
        sSy = mf.SpinOp(Ssys, stee)
        stee[0]='z'
        sSz = mf.SpinOp(Ssys, stee)
        stee[0]='e'

        HFHam = np.zeros_like(sSx)
        Ix = [None]*len(gamma)
        Iy = [None]*len(gamma)
        Iz = [None]*len(gamma)

        for n_nuc in range(len(gamma)):

            nu_n = gamma[n_nuc]*self.Exp.Field
            #mult = int(2*self.Sys.I[n_nuc]+1)

            stee[n_nuc+1]='x'
            Ix[n_nuc] = mf.SpinOp(Ssys, stee)
            stee[n_nuc+1]='y'
            Iy[n_nuc] = mf.SpinOp(Ssys, stee)
            stee[n_nuc+1]='z'
            Iz[n_nuc] = mf.SpinOp(Ssys, stee)
            stee[n_nuc+1]='e'

            #tHam.append(-nu_n*(lx*self.sIx[mult] + ly*self.sIy[mult] + lz*self.sIz[mult]) )
            HFHam += (
                     HFc[n_nuc][0, 0] * Ix[n_nuc]@sSx + HFc[n_nuc][0, 1] * Ix[n_nuc]@sSy + HFc[n_nuc][0, 2] * Ix[n_nuc]@sSz +
                     HFc[n_nuc][1, 0] * Iy[n_nuc]@sSx + HFc[n_nuc][1, 1] * Iy[n_nuc]@sSy + HFc[n_nuc][1, 2] * Iy[n_nuc]@sSz +
                     HFc[n_nuc][2, 0] * Iz[n_nuc]@sSx + HFc[n_nuc][2, 1] * Iz[n_nuc]@sSy + HFc[n_nuc][2, 2] * Iz[n_nuc]@sSz
                     )
        self.ak = np.zeros_like(self.theta)
        for n_ang in range(len(self.phi)):

            costheta = np.cos(self.theta[n_ang])
            sintheta = np.sin(self.theta[n_ang])
    
            cosphi = np.cos(self.phi[n_ang])
            sinphi = np.sin(self.phi[n_ang])

            lx = sintheta*cosphi
            ly = sintheta*sinphi
            lz = costheta

            ############### Build full hamiltonian ############################
            Gx=lx*self.Sys.g[0][0]
            Gy=ly*self.Sys.g[0][1]
            Gz=lz*self.Sys.g[0][2]

            Ham = (sSx*Gx + sSy*Gy + sSz*Gz)*self.Exp.Field/self.h_ov_beta_mTMHz


            for n_nuc in range(len(gamma)):
                Ham += (lx*Ix[n_nuc] + ly*Iy[n_nuc] + lz*Iz[n_nuc])*self.Exp.Field*gamma[n_nuc]
            Ham+=HFHam

            D, V = np.linalg.eigh(Ham)

            ctd = np.cos(self.theta[n_ang]+np.pi/2.0)
            std = np.sin(self.theta[n_ang]+np.pi/2.0)

            dlx = std*cosphi
            dly = std*sinphi
            dlz = ctd

            DHam = (sSx*dlx + sSy*dly + sSz*dlz)

            #tF = D[None, :]-D[:, None]

            for ii in range(Ham.shape[0]):
                for  jj in range(ii+1,Ham.shape[0]):
                    df =abs(abs(D[ii]-D[jj])-mwFreqMHz)
                    if df<5.0*sig:
                        self.ak[n_ang] += np.abs(V[:, ii].conj()@DHam@V[:, jj].T) * np.exp(-0.5 * (df / sig)**2)

        if np.max(self.ak) < 1e-4:
            self.errorFunc(f'HYSCORE_computeOrisel_eig: no orientations were found with P(θ,ϕ)>1e-4. Aborting the run. \n\t mwFreq = {self.Exp.mwFreq :.8} [GHz] \n\t Field = {self.Exp.Field :.8} [mT]')

        # Apply the grid weights, as compute_weights does. Grids like 'sphgrid'
        # are not equal-area, so skipping this would bias the powder average
        # toward the poles. Deliberately after the threshold test above, which
        # has to see the unweighted ak: normalised weights are O(1/N) and would
        # trip it on every run.
        if self.ww is not None:
            self.ak = self.ak*self.ww

        OriSelOut[:, 2] =self.ak
    
    
        if setActive:
            self.Opt.OriSelInp = OriSelOut
        if returnMatrix:
            return OriSelOut
# ----------------------------------------------------------------------
# Example usage
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib import cm


    hs = HYSCOREsim()

    hs.Sys.set( g = [np.array([2.2, 2.1, 2.0])],
                S = [0.5],
                Nucs = ['1H'],
                A = [np.array([5, 4, 6])*100.0],
                Apa = [np.array([0, 0, 0])],
                Q = [np.array([0, 0, 0])*0.0],
                Qpa = [np.array([0, 0, 0])],
                lw = [0.1],
                useFor = ['all']
                )
    
    hs.Exp.set(mwFreq = 9.5,
               Field = 325,
               ExciteWidth=25,
               nPoints=512,
               MaxFreq = 10,
               )
    
    #hs.Opt.useHFsel = True
    hs.Opt.OriSelType = ['brute force', ['g_eff', 'g_eff+HFC', 'brute force', 'precalculated']]
    hs.Opt.nKnots = 200

    #matr = hs.computeOrisel_eig(nKnots=hs.Opt.nKnots, returnMatrix=True, setActive=True)
    # phi = matr[:, 0]
    # theta = matr[:, 1]
    # ak = matr[:, 2]


    hs.run()
    phi, theta, ak = hs.phi,hs.theta,hs.ak
    #hs.make_grid(grid='fibonacci', nKnots=150)

    
    #hs.compute_weights(phi, theta, ww)
        
    # Convert spherical to Cartesian for plotting
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)

    fig = plt.figure(1, figsize=(6,5))
    fig.clf()
    ax  = fig.add_subplot(111, projection='3d')
    
    norm = plt.Normalize(vmin=0, vmax=np.max(ak))

    tri = mtri.Triangulation(x, y)    
    triangle_ak = ak[tri.triangles].max(axis=1)          # one value per triangle
    facecolors = cm.jet(norm(triangle_ak))

    surf = ax.plot_trisurf(
        x, y, z,
        triangles=tri.triangles,
        linewidth=0, antialiased=True,
        shade=False,  # we supply colour via facecolors
        )
        #
    
    surf.set_facecolor(facecolors)
    #surf.set_edgecolor(cm.jet(norm(ak)))
    #surf.set_array(ak/np.max(ak))
    #surf.autoscale()
    #ax.scatter(x, y, z, s=1, c='b', alpha=0.3)
    #ax.set_box_aspect([1,1,1])
    ax.set_aspect('equal') 
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.title.set_text('Generated grid')

    # fig1 = plt.figure(2, figsize=(6,5))
    # fig1.clf()
    # ax1  = fig1.add_subplot(121)
    # ax2  = fig1.add_subplot(122)
    # ax1.imshow(
    #     np.abs(hs.rawSpec),
    #     extent=[-hs.Exp.MaxFreq, hs.Exp.MaxFreq, -hs.Exp.MaxFreq, hs.Exp.MaxFreq],
    #     origin="lower",
    #     cmap="jet",
    #     interpolation="nearest",
    # )
    # ax1.title.set_text("HYSCORE")
    # ax1.set_xlabel(r"$\nu_a, MHz$")
    # ax1.set_ylabel(r"$\nu_b, MHz$")
    
    # ax2.imshow(
    #     np.abs(hs.Spectrum),
    #     extent=[-hs.Exp.MaxFreq, hs.Exp.MaxFreq, -hs.Exp.MaxFreq, hs.Exp.MaxFreq],
    #     origin="lower",
    #     cmap="jet",
    #     interpolation="nearest",
    # )
    # ax2.title.set_text("HYSCORE")
    # ax2.set_xlabel(r"$\nu_a, MHz$")
    # ax2.set_ylabel(r"$\nu_b, MHz$")
    
    # plt.tight_layout()
    
    plt.show()
    del(hs)

    SZ2 = mf.SpinOp([1/2, 1/2], ['e','z'])
    print(SZ2)