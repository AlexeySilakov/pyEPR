# pyHYSCORE — Build Specification for an AI Coding Agent

**Target:** a from-scratch re-implementation of *pyHYSCORE* (v0.9 in this folder),
a desktop application that reads Bruker 2D HYSCORE data, processes it (baseline,
apodization, zero-fill, 2D FFT), simulates HYSCORE spectra from an EasySpin-style
spin system, overlays experiment and simulation, and fits two spin-Hamiltonian
parameters at a time with a 2D grid search that reports RMSD maps, Bayesian
posteriors and credible intervals.

This document is the single source of truth for the agent. It describes *what*
to build, *which algorithms and equations* to use, *how the modules fit together*,
and *which defects of v0.9 must not be carried over*. Where the reference code
and this document disagree, this document wins; where a choice is genuinely the
author's, it is marked **[AUTHOR DECISION]** and must be asked, not guessed.

---

## 0. Non-negotiable priorities

Three priorities override every other consideration, in this order of tie-break:

### 0.1 Code simplicity (highest)

The author debugs by stepping through the code. A shallow call stack and few
lines matter more than "clean architecture".

1. **If code is called from exactly one place, it is not a function.** Write it
   inline at the call site. No single-use helpers, no "extract method" refactors.
2. **Algorithm variants chosen by `if/elif` live inline in their branch**, however
   long (a 60-line `elif grid == 'sphgrid':` branch is correct style).
3. **Only logic with ≥ 2 real callers becomes a function**, and then it is a
   **public method on the relevant class**. No module-level helper functions
   (exception: the worker entry points that `multiprocessing` must pickle, §11.4),
   **no leading-underscore names** — every method is public.
4. **Keep line count down.** Write the shortest correct expression of a numeric
   routine. No defensive guards for inputs that cannot occur, no speculative
   options, no configuration nobody asked for. Docstrings say what the next reader
   needs (units, shapes, the one non-obvious fact) — not a derivation.
5. Comments explain *why*, not *what*. Match the density of the surrounding code.
6. No new dependencies beyond §1.2 without the author's approval.

### 0.2 Cross-platform compatibility

Must run unmodified on **Windows 10/11, macOS (Intel + Apple Silicon) and Linux**
from a plain `pip install`. Rules in §13.

### 0.3 Computational performance

Simulation and grid search are the product. They must be vectorised NumPy with no
per-orientation Python loops in the hot path, no `np.add.at`, no repeated
allocation inside loops, and multiprocessing for the grid search. Rules and
targets in §12.

When simplicity and performance conflict, prefer the simpler code **unless** the
hot path (§12.1) is affected by more than ~2×.

---

## 1. Technology stack and repository layout

### 1.1 Language

Python **≥ 3.10**. Type hints optional; f-strings; `pathlib` or `os.path` (never
string-concatenated paths).

### 1.2 Dependencies (complete list)

| Package | Use |
|---|---|
| `numpy` | everything numeric |
| `scipy` | `spatial.Delaunay`, `spatial.ConvexHull`, `interpolate.RegularGridInterpolator`, `interpolate.RectBivariateSpline`, `special.ndtri`, `stats` |
| `matplotlib` | all plotting (`backend_wxagg`), colormaps |
| `wxPython` ≥ 4.2 | GUI (kept for continuity with the sibling app *pyEPR*, which shares `SysPar`, the Bruker reader and the property-grid widgets) |

Standard library only otherwise (`xml.etree`, `base64`, `ast`, `multiprocessing`,
`multiprocessing.shared_memory`, `threading`, `math`, `re`, `copy`, `time`).

### 1.3 Module map (keep it this flat — one file per concern)

| File | Contents | wx allowed? |
|---|---|---|
| `pyHYSCORE.py` | `MainFrame`, `FileBrowserPanel`, `MatplotlibPanel`, `TabulatedPanel`, session XML, FFT processing, orchestration | yes |
| `hyscore_sim.py` | `optHYSCORE` (options), `HYSCOREsim` (grids, operators, transitions, binning, broadening, orientation selection, analysis-box geometry) | **no** |
| `SysPar.py` | `sysPar`, `expPar`, isotope table | **no** |
| `mathfunctions.py` | `RotMatrix`, `RotTransDiag`, `RotTrans3x3`, `SpinOp`, `adaptSim`, `calcSSR` | **no** |
| `sys_functions.py` | expression engine for Function Modifiers | **no** |
| `brukerread.py` | BES3T (`.DSC/.DTA`) and ESP (`.par/.spc`) reader | **no** |
| `grid_search.py` | `GridSearchFrame`, `HeatmapPanel`, `ResidualPanel` | yes |
| `grid_search_worker.py` | everything one grid point needs; `SimStore`; pool entry points | **no** |
| `classPropGridPanel.py` | `PropGridPanel`, `SpinCtrlDoubleEditor` (spin box with right-click step menu) | yes |
| `classFunctionModPG.py` | `FunctionModPanel` (wx front end of `sys_functions`) | yes |
| `theme.py` | light/dark palettes, pill widgets, icon buttons, themed toolbar, redraw debouncer, resize freeze guard | yes |
| `colormap_editor.py` | `ColormapEditorDialog`, `colormap_to_stops`, `stops_to_cmap` | yes |
| `icons/` | image assets, located relative to `__file__` (never the CWD) | — |
| `tests/` | headless pytest checks of the non-wx modules (§17) | — |

**Hard rule:** the four "no-wx" modules plus `grid_search_worker.py` must import
and run with wxPython absent. This is what makes the grid-search worker processes
cheap to spawn and lets the engine be tested headlessly.

---

## 2. Units and conventions (use everywhere, no exceptions)

| Quantity | Unit | Notes |
|---|---|---|
| Hyperfine `A`, quadrupole `Q`, `Q_K`, linewidth `lw`, `ExciteWidth`, all spectral frequencies | **MHz** | |
| Microwave frequency `mwFreq` | **GHz** | ×10³ when compared with MHz |
| Field `Field`, `BMin`, `BMax` | **mT** | `gField` is a g-value; if set, `Field` is derived (§6.5) |
| `tau` | **ns** | ×10⁻³ → µs when multiplied by MHz |
| `tDead` | **ns** | ×10⁻³ → µs (v0.9 used it unconverted — fix) |
| Time axes of loaded data | **µs** | converted from ns on load |
| Euler angles `Apa`, `Qpa`, `AFrame`, `QFrame`, `gFrame` | **degrees** | radians internally |
| Nuclear `gamma` | **MHz/mT** | `nu350mT/350` from the isotope table |
| Orientation angles `phi`, `theta` | radians | θ from +z, φ from +x |

Physical constants (CODATA 2018):
`h = 6.62607015e-34 J s`, `μB = 9.2740100783e-24 J/T`, `μN = 5.050783699e-27 J/T`.
Derived: `h_ov_beta_mTMHz = h/μB·1e9 = 0.0714477 mT/MHz`, so
$f_e\,[\mathrm{MHz}] = B_0\,[\mathrm{mT}]\cdot g_\mathrm{eff} / 0.0714477$.

"−1" in any `Exp` field means **"take it from the dataset"** (§9.3).

Array index convention for 2D maps: **axis 0 ↔ `x` axis (ν₁ / t₁), axis 1 ↔ `y`
axis (ν₂ / t₂)**. Every function that takes two axes takes them in that order.
For `imshow` the array is transposed or the axes are swapped consistently — the
plotting code alone decides; the engine never does.

Frequency axes are always
`np.linspace(-F, F, N, endpoint=False)` ≡ `np.fft.fftshift(np.fft.fftfreq(N, dt))`
with `F = 1/(2·dt)`. **Never** `endpoint=True` (it stretches every frequency by
N/(N−1) and makes peak positions depend on the zero-fill).

---

## 3. Parameter containers (`SysPar.py`, `hyscore_sim.optHYSCORE`)

### 3.1 `sysPar` — spin system (EasySpin-inspired)

Attributes default to `None`. All per-nucleus attributes are **lists indexed by
nucleus type** (not by equivalent nucleus).

| Attribute | Shape per entry | Meaning |
|---|---|---|
| `S` | float | electron spin (simulation supports S = 1/2; §6.1) |
| `g` | ndarray(3) | principal g values; assumed collinear with the molecular frame |
| `gFrame` | ndarray(3) | reserved (ignored in v0.9 — keep ignored but listed) |
| `lw` | float | Gaussian FWHM broadening of the HYSCORE map, MHz |
| `Nucs` | str, e.g. `'14N'` | isotope label; sets `I`, `gn`, `gamma` from the table |
| `nNucs` | int | number of equivalent nuclei of this type |
| `I`, `gn`, `gamma` | float | from the table (`gamma = nu350mT/350`) |
| `A` | ndarray(3) principal values **or** ndarray(3,3) full tensor | MHz |
| `Apa` | ndarray(3) | Euler angles (deg) rotating A into the g frame (old EasySpin) |
| `AFrame` | ndarray(3) | reserved (new EasySpin convention) |
| `Q_K`, `Q_eta` | float | quadrupole K (MHz) and asymmetry η; converted to `Q` |
| `Q` | ndarray(3) or (3,3) | quadrupole tensor, MHz |
| `Qpa` | ndarray(3) | Euler angles of Q (deg) |
| `useFor` | `'all'` \| `'sim.only'` \| `'ori.sel.only'` | whether the nucleus contributes to the spectrum, to orientation selection, or both |

Methods (public): `set(**kw)`, `setDict(d)`, `setNucs(Nucs)`, `setFromCtrl(paramDict)`,
`getDefaultDict()`, `getAll()`, `set_All_None()`, `isotopes(Nucs=None)`,
`getToolTips(param=None)`. Unknown keys raise `AttributeError`.

**Quadrupole from (K, η)** — use the EasySpin convention:
$$Q = K\,\mathrm{diag}\big(-(1-\eta),\; -(1+\eta),\; 2\big),\qquad K=\frac{e^2qQ}{4I(2I-1)h}.$$
v0.9 builds `[(1+η)K, (1−η)K, −2K]` (= EasySpin with K→−K, η→−η) and also has a
logic bug that zeroes η when it is given and raises `NameError` when it is not.
**[AUTHOR DECISION]** confirm EasySpin convention before implementing; do not
replicate the bug either way. `Q` is only built for `I ≥ 1`.

**GUI parameter dictionary** (what `PropGridPanel.parameters` holds and what
`setFromCtrl` consumes, also what is saved in sessions):

```python
{
 'spin(1)': {'S': 0.5, 'g': np.array([2., 2., 2.]), 'lw': 0.1},
 'nuc(1)':  {'Nucs': ['1H', [<all isotope labels>]], 'nNucs': 1,
             'A': np.array([1., 1., -2.]), 'Apa': np.zeros(3),
             'Q_K': 0.0, 'Q_eta': 0.0, 'Qpa': np.zeros(3),
             'useFor': ['all', ['sim.only', 'ori.sel.only', 'all']]},
 'nuc(2)':  {...},
 'functions': {'func': [['nuc(1).A(3)', '59/r**3 + a'], ...],
               'var':  {'r': 2.5, 'a': 0.3}},
}
```

Rules: keys are `name(number)`, 1-based; `[value, [choices]]` is an enum;
`setFromCtrl` resets every `getAll()` attribute (and `useFor`), then fills list
index `number−1`; the `'functions'` key is skipped by `setFromCtrl` (it is applied
earlier, §10). Default `useFor` is `'all'` (v0.9 mixes `'both'` and `'all'` —
use `'all'` only).

**Isotope table:** copy the table in `SysPar.isotopes()` verbatim (≈115 isotopes,
fields `abundance %, I, gn, Q (barn), nu350mT (MHz)`). Its comment notes it was
partly AI-generated; keep the comment and the source list. Do not alter values.

### 3.2 `expPar` — experiment

Fields and GUI defaults (`getDefaultDict`):

| Field | Default | Meaning |
|---|---|---|
| `mwFreq` | −1 | GHz |
| `Field` | −1 | mT |
| `gField` | None | if set: `Field = mwFreq/gField · h_ov_beta_mTMHz · 1e3` |
| `tau` | −1 | ns |
| `nPoints` | 256 | simulation grid size per axis (−1 → dataset FFT length) |
| `MaxFreq` | −1 | simulated half-range, MHz (−1 → dataset `fmax`) |
| `ExciteWidth` | 20 | FWHM of the excitation profile, MHz |
| `tDead` | 0 | ns |

Also present (non-GUI): `BMin`, `BMax`, `Tinv`, `phi`, `theta`, `Harmonic`.

### 3.3 `optHYSCORE` — simulation options

| Option | Default | Choices / meaning |
|---|---|---|
| `Grid` | `'fibonacci'` | `fibonacci`, `sphgrid`, `spiral` (stored as `[value, choices]`) |
| `nKnots` | 20 | grid resolution in each method's own units (§5) |
| `OriSelType` | `'g_eff'` | `g_eff`, `g_eff+HFC`, `brute force`, `precalculated` |
| `OriSelInp` | None | N×3 array `[phi, theta, weight]` for `precalculated` |
| `Symmetry` | `'Ci'` | only Ci is implemented (upper hemisphere) |
| `KillNeg` | True | suppress negative-time half after broadening (§6.4) |
| `ProdRule` | True | product rule for >1 nucleus (§6.3) |
| `AmpRatios` | 1.0 | per-nucleus amplitude multipliers (scalar or list) |
| `Treshold` | 1e-3 | drop orientations with normalised weight below this (§12.1 item 6) |
| `Verbosity` | False | print timings |

Methods (public): `set`, `setFromCtrl`, `getDefaultDict`, `getToolTips`,
`nGridPoints(grid, nKnots)`, `nKnotsForGridPoints(target, grid)`,
`diagonalBox(axis0, axis1, fmin, fmax, antidiagspan)` (§8.1).

v0.9 declares `ProdRule`, `AmpRatios`, `Treshold`, `Symmetry` but never uses them;
the rebuild must either implement them as specified here or remove them from the
GUI. Implementing is preferred.

---

## 4. Mathematical foundations of the simulation

### 4.1 Rotations and tensors (`mathfunctions.py`)

Passive z-y-z Euler rotation with angles (α, β, γ) in radians:

$$
R=\begin{pmatrix}
c_\gamma c_\beta c_\alpha - s_\gamma s_\alpha & c_\gamma c_\beta s_\alpha + s_\gamma c_\alpha & -c_\gamma s_\beta\\
-s_\gamma c_\beta c_\alpha - c_\gamma s_\alpha & -s_\gamma c_\beta s_\alpha + c_\gamma c_\alpha & s_\gamma s_\beta\\
s_\beta c_\alpha & s_\beta s_\alpha & c_\beta
\end{pmatrix}
$$

- `RotTransDiag(a3, ang) = R · diag(a3) · Rᵀ`
- `RotTrans3x3(A, ang) = R · A · Rᵀ`

A tensor given as 3 values goes through `RotTransDiag(values, Apa·π/180)`; a 3×3
tensor through `RotTrans3x3`. (v0.9 has a typo `Qpa[ii].np.pi` in the 3×3 Q
branch — do not copy.)

### 4.2 Spin operators

`SpinOp(spins, corrs)` returns $\bigotimes_k O_k$ where for spin $s$, dimension
$d=2s+1$, basis ordered $m=-s,\dots,+s$:

$$
(S_+)_{m+1,m}=\sqrt{s(s+1)-m(m+1)},\quad S_x=\tfrac12(S_++S_-),\quad S_y=\tfrac1{2i}(S_+-S_-),\quad S_z=\mathrm{diag}(m)
$$

and `'e'` is the identity. Build the single-spin matrices with vectorised NumPy
(`np.diag(vals, 1)`), not Python loops. Cache per dimension `d` in the
`HYSCOREsim` instance: $I_x, I_y, I_z$ and, for $d>2$, the nine products
$I_aI_b$ needed for the quadrupole term.

### 4.3 Nuclear sub-Hamiltonians (high-field, S = 1/2)

For orientation unit vector $\mathbf l=(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$
of $B_0$ in the molecular (g) frame, the electron spin is quantised along the
g-weighted direction

$$
\mathbf n = \frac{g\,\mathbf l}{|g\,\mathbf l|},\qquad g_\mathrm{eff}(\mathbf l)=|g\,\mathbf l|=\sqrt{\textstyle\sum_i (g_i l_i)^2}.
$$

For each nucleus type $k$ (treated independently, §6.3), with
$\nu_I = \gamma_k B_0$:

$$
H_{\alpha} = -\nu_I\,(\mathbf l\cdot\mathbf I) - \tfrac12\,\mathbf n^{\mathsf T} A^{\mathsf T}\mathbf I + \mathbf I^{\mathsf T}Q\,\mathbf I,\qquad
H_{\beta} = -\nu_I\,(\mathbf l\cdot\mathbf I) + \tfrac12\,\mathbf n^{\mathsf T} A^{\mathsf T}\mathbf I + \mathbf I^{\mathsf T}Q\,\mathbf I
$$

($\alpha \leftrightarrow m_S=-\tfrac12$, $\beta \leftrightarrow m_S=+\tfrac12$).
Implementation: precompute per nucleus the three operators
$\mathrm{HFI}_j=\sum_i I_i A_{ij}$ ($j=x,y,z$) and the constant $\mathbf I^\mathsf T Q\mathbf I$;
then per orientation $H_{\alpha,\beta}=Z \mp \tfrac12\sum_j n_j\mathrm{HFI}_j + Q_\mathrm{op}$.

Sanity identity (use as a unit test, §17): for $I=\tfrac12$, $Q=0$ the
eigenfrequencies are exactly
$\nu_{\alpha/\beta}=\big|\,\mp\tfrac12A^{\mathsf T}\mathbf n-\nu_I\mathbf l\,\big|$.

### 4.4 Diagonalisation and overlap

$H_\alpha=V_\alpha D_\alpha V_\alpha^\dagger$, $H_\beta=V_\beta D_\beta V_\beta^\dagger$
(`np.linalg.eigh`, columns = eigenvectors). Nuclear frequencies
$\omega^\alpha_{ij}=D_{\alpha,i}-D_{\alpha,j}$, $\omega^\beta_{ln}=D_{\beta,l}-D_{\beta,n}$.
Overlap $M=V_\alpha^\dagger V_\beta$ is unitary, so $M^{-1}=M^\dagger$ (never invert).

### 4.5 HYSCORE amplitudes (4-pulse, ideal pulses)

Sequence π/2–τ–π/2–t₁–π–t₂–π/2–τ–echo. Each contributing term has amplitude

$$
a_{ijklmn}= M_{il}\,M^*_{jl}\,M_{jn}\,M^*_{kn}\,M_{km}\,M^*_{im}
$$

(v0.9: `M[i,l]*cM[l,j]*M[j,n]*cM[n,k]*M[k,m]*cM[m,i]` with `cM = M.conj().T`),
the τ-dependent ("blind-spot") factor

$$
b_{ijklmn}(\tau)= e^{-2\pi i(\omega^\alpha_{ij}+\omega^\beta_{lm})\tau}+e^{+2\pi i(\omega^\alpha_{kj}+\omega^\beta_{nm})\tau},
$$

the optional dead-time phase $p=e^{-2\pi i(\omega^\alpha_{ik}+\omega^\beta_{ln})t_\mathrm{dead}}$,
and the orientation weight $w(\mathbf l)$ (§7). It contributes two peaks:

| peak at (ν₁, ν₂) | complex amplitude |
|---|---|
| $(\omega^\alpha_{ik},\ \omega^\beta_{ln})$ | $a\,b\,p\,w$ |
| $(\omega^\beta_{ln},\ \omega^\alpha_{ik})$ | $a^*\,b\,p\,w$ |

#### 4.5.1 Mandatory contraction (performance-critical)

The naive sum runs over $d^6$ index tuples per orientation ($d=2I+1$:
64, 729, 4096, 46 656, 262 144 for I = ½…7/2). Because $b$ is a sum of two
products that each separate in $j$ and in $m$, contract $j$ and $m$ first:

$$
\begin{aligned}
X^{(1)}_{iln}&=\textstyle\sum_j M^*_{jl}M_{jn}\,e^{-2\pi i\omega^\alpha_{ij}\tau}, &
Y^{(1)}_{ikl}&=\textstyle\sum_m M_{km}M^*_{im}\,e^{-2\pi i\omega^\beta_{lm}\tau},\\
X^{(2)}_{kln}&=\textstyle\sum_j M^*_{jl}M_{jn}\,e^{+2\pi i\omega^\alpha_{kj}\tau}, &
Y^{(2)}_{ikn}&=\textstyle\sum_m M_{km}M^*_{im}\,e^{+2\pi i\omega^\beta_{nm}\tau},
\end{aligned}
$$

$$
\mathcal A_{ikln}= M_{il}\,M^*_{kn}\,\big(X^{(1)}_{iln}Y^{(1)}_{ikl}+X^{(2)}_{kln}Y^{(2)}_{ikn}\big).
$$

This gives the identical spectrum with $d^4$ peaks and $O(d^5)$ work. Implement
with `np.einsum(..., optimize=True)` **batched over orientations** (leading axis
`o`). Verify against the brute $d^6$ sum in a test for I = ½ and I = 1 (§17).

### 4.6 Stick spectrum → binned map

Collect all peaks of all orientations (and nuclei, §6.3) as flat arrays
`(f1, f2, amp)` and bin onto an N×N complex grid:

- `dx = N / (2·MaxFreq)`; `i1 = round((f1 + MaxFreq)·dx)`, `i2 = round((f2 + MaxFreq)·dx)`.
- keep `0 ≤ i < N−1` on both axes and drop peaks with `f1 == 0` or `f2 == 0`
  exactly (zero-frequency combinations within a manifold).
- accumulate with **`np.bincount`** on the flattened index `i1·N + i2`, once for
  `amp.real` and once for `amp.imag` (`minlength = N²`) — 5–20× faster than
  `np.add.at`.
- symmetrise: add the same contributions at the transposed index (`i2·N + i1`).

### 4.7 Line broadening and negative-time suppression

With `rawSpec` the binned complex map and σ = `lw / 2.354820045`:

1. `s = ifft2(ifftshift(rawSpec))` (to the time domain).
2. Gaussian convolution kernel in the time domain — use the analytic transform
   instead of FFT'ing a sampled Gaussian:
   $g(t_1,t_2)=e^{-2\pi^2\sigma^2(t_1^2+t_2^2)}$ on the FFT time grid
   (`t = fftfreq(N, d=1/(2·MaxFreq))`, i.e. ordered as `ifft2` output).
   (v0.9 builds the frequency-domain Gaussian and FFTs it; same result, one more
   FFT.)
3. `KillNeg`: zero the second half of the time indices on both axes
   (`g[N//2:, :] = 0; g[:, N//2:] = 0`), which removes the mirror (negative-time)
   artefacts of the circular transform.
4. `Spectrum = |fftshift(fft2(s · g))|`.

If `lw == 0`: `Spectrum = |rawSpec|`. Simulated axes `X = Y = linspace(-MaxFreq, MaxFreq, N, endpoint=False)`.

---

## 5. Orientation grids (`HYSCOREsim.make_grid`)

All grids cover the **upper hemisphere** (Ci symmetry) except `spiral` per its own
definition below. Every grid returns `phi`, `theta`, weights `ww` **normalised to
Σ = 1** (so intensity is independent of grid size and comparable across
methods), and a triangulation `tri`, `triareas` (§5.4). Each method is one inline
`if/elif` branch.

`nKnots` keeps each method's native meaning (deliberately *not* rescaled):

| Grid | Points for `nKnots = K` | K = 20 → |
|---|---|---|
| `fibonacci` | `int(K**1.5)` | 89 |
| `sphgrid` | `1 + 2K(K−1)` (K ≥ 2) | 761 |
| `spiral` | $\sum_{r=1}^{K-1}\lfloor\sin(r\pi/2K)\cdot4K\rfloor$ (K ≥ 2) | 968 |

`optHYSCORE.nGridPoints(grid, nKnots)` returns these counts without building the
grid; a test asserts it equals `len(make_grid(...).phi)` for K = 2…60 for every
method.

### 5.1 `fibonacci` (golden-angle, equal area)

$n=\lfloor K^{1.5}\rfloor$, $\varepsilon=0.33$, golden angle $\psi=\pi(3-\sqrt5)$:

$$
\cos\theta_i=\frac{i+\varepsilon}{n-1+2\varepsilon},\qquad \phi_i=(i\,\psi)\bmod 2\pi,\qquad i=0,\dots,n-1,\qquad w_i=1/n .
$$

If `Symmetry` other than Ci is requested, report "not implemented, full azimuth
used" through the error callback but continue.

### 5.2 `sphgrid` (SOPHE / triangular, EasySpin Ci)

Port of `sphgrid_` from EasySpin `sphgrid.m` (Stoll & Schweiger, *J. Magn. Reson.*
178 (2006) 42, Eq. 15; Wang & Hanson, *J. Magn. Reson. A* 117 (1995) 1). With
K = max(2, nKnots), nOct = 4, maxPhi = 2π, Δθ = (π/2)/(K−1):

1. North pole: θ = 0, φ = 0, $w_0 = \text{maxPhi}\,(1-\cos(\Delta\theta/2))$.
2. Slices s = 2…K−1: θ = (s−1)Δθ, nPhi = nOct(s−1)+1 points
   φ = linspace(0, maxPhi, nPhi), ΔΦ = maxPhi/(nPhi−1),
   $w = 2\sin\theta\,\sin(\Delta\theta/2)\,\Delta\Phi\cdot[1,1,\dots,1,\tfrac12]$.
3. Equator: θ = π/2, nPhi = nOct(K−1)+1, $w=\sin(\Delta\theta/2)\,\Delta\Phi\cdot[1,\dots,1,\tfrac12]$
   (no factor 2 sinθ: half a band on a hemisphere).
4. Drop the φ = 2π duplicate closing every ring (indices `cumsum(nOct·arange(1,K)+1)`).
5. Weights ×2·(2π/maxPhi) (EasySpin: Σ = 4π), then normalise to Σ = 1.

This grid is **not** equal-area; the weights are essential.

### 5.3 `spiral` (latitude rings, after E. Reijerse)

krid = max(2, K), step = (π/2)/krid. Walk rings from the equator to the pole:
ring index r = krid−1 … 1, ring centre θ_m = r·step, nPhi = ⌊sin(θ_m)·4·krid⌋,
start θ = θ_m + step/2, then for q = 0…nPhi−1:
θ = θ_m + step/2 − q·step/nPhi, φ = q·2π/nPhi. A ring is full after **its own**
nPhi knots (v0.9's earlier bug stacked surplus knots on the pole). Equal area ⇒
uniform weights, normalised.

### 5.4 Triangulation and solid angles (`grid_triangulation(phi, theta)`, public — two callers)

- Hemisphere (`max θ ≤ π/2`): Delaunay in the azimuthal-equidistant projection
  $(\theta\cos\phi,\ \theta\sin\phi)$ — a bijection on the hemisphere.
- Otherwise: `ConvexHull` of the unit vectors (for points on a sphere the hull
  *is* the spherical Delaunay).
- Triangle solid angle by l'Huilier: edge arcs $a_k=\arccos(\mathbf x_i\cdot\mathbf x_j)$,
  $s=\tfrac12\sum a_k$,
  $E = 4\arctan\sqrt{\max(0,\ \tan\tfrac s2\tan\tfrac{s-a_1}2\tan\tfrac{s-a_2}2\tan\tfrac{s-a_3}2)}$.
- Drop slivers with area < 10⁻³ × mean; sort each triangle's indices, then rows
  lexicographically; rescale areas to Σ = 4π.

Used by the orientation-selection 3D plot (§14.4).

### 5.5 Target grid size (`nKnotsForGridPoints`)

Since only a ladder of totals exists, find the `nKnots` whose count is closest to
a requested total: bracket by doubling from 2, then bisect (counts are monotone),
then compare with the rung below; ties → smaller grid. The GUI shows "actual: N
points (nKnots = K)" next to the requested target.

### 5.6 Grid cache (`grid_signature`, `ensure_grid`)

The grid and its weights are rebuilt **only** when something they depend on
changes. `grid_signature()` returns a hashable tuple of exact (unrounded) values:
grid name, nKnots, Symmetry, useFor, OriSelType, OriSelInp, `Exp.mwFreq/Field/
gField/ExciteWidth`, `Sys.S/g/gFrame/D/DFrame`, and per-nucleus
`Nucs, nNucs, I, gn, gamma, A, Apa, AFrame, Q, Qpa, QFrame` **restricted to nuclei
not flagged `'sim.only'`** (they cannot move the grid). Arrays are tagged as
`(shape, dtype, bytes)`. Excluded on purpose: tau, nPoints, MaxFreq, lw,
ProdRule, KillNeg, Treshold, AmpRatios.

`ensure_grid()` compares, rebuilds if needed (`make_grid` + `compute_weights`, or
`computeOrisel_eig` for brute force), and stores the signature **after** the
build (because `compute_weights` may write `Exp.Field` from `gField`). Put a
banner comment above these two methods saying "add grid dependencies here".

---

## 6. The simulation engine (`HYSCOREsim`)

### 6.1 Public interface

```python
hs = HYSCOREsim(Sys=sysPar, Exp=expPar, Opt=optHYSCORE, errorFunc=callable)
hs.preCompute()   # spin operators, tensors (depends on Sys only)
hs.reRun()        # grid (cached), transitions, binning, broadening
hs.run()          # preCompute + reRun
# results:
hs.Spectrum (N×N float), hs.X, hs.Y, hs.rawSpec,
hs.phi, hs.theta, hs.ak (orientation weights incl. grid weights), hs.tri, hs.triareas
```

`errorFunc(message)` is the only error channel; the engine calls it and then
**raises** a plain exception class defined in `hyscore_sim.py` (`HYSCOREError`)
so callers never depend on the callback aborting (v0.9's GUI callback showed a
modal dialog *and* raised a bare `Exception`; the grid search passed a callback
that raised its own class — keep one exception class, and let the GUI show the
dialog when it catches it).

Only S = 1/2 is simulated. If `Sys.S != 0.5`, raise with a clear message (v0.9's
README claims arbitrary S; the engine never supported it).

### 6.2 `reRun()` algorithm

1. Orientation weights:
   - `precalculated`: take `phi, theta, ak` from `Opt.OriSelInp` (validate N×3),
     `tri = None`, clear the grid signature.
   - otherwise `ensure_grid()`.
2. Apply `Treshold`: keep orientations with `ak ≥ Treshold·max(ak)` (cheap and
   often removes 50–90 % of the work at narrow excitation widths).
3. For every nucleus type with `useFor != 'ori.sel.only'`, **batched over all
   kept orientations at once**:
   1. Build stacked Hamiltonians `Hα, Hβ` of shape `(n_ori, d, d)` by broadcasting
      `l`, `n` against the cached operators.
   2. `Dα, Vα = np.linalg.eigh(Hα)`; same for β (eigh is batched natively).
   3. `M = Vα.conj().transpose(0,2,1) @ Vβ`; `ωα = Dα[:, :, None] − Dα[:, None, :]`,
      likewise `ωβ`.
   4. Contracted amplitudes `𝒜[o,i,k,l,n]` (§4.5.1), times `ak[o]`, times the
      dead-time phase, times `AmpRatios[k]`.
   5. Peak frequencies `ωα[o,i,k]`, `ωβ[o,l,n]` broadcast to the same shape.
4. Combine nuclei (§6.3), bin (§4.6), broaden (§4.7).

Chunk the orientation axis if `n_ori·d⁴` would exceed ~5·10⁷ complex numbers, so
memory stays bounded for I ≥ 5/2.

### 6.3 Several nuclei: sum vs product rule

- **Sum (`ProdRule = False`)** — what v0.9 actually does: concatenate the peak
  lists of all nucleus types. Fast; no combination peaks.
- **Product rule (`ProdRule = True`)** — exact for nuclei coupled only through the
  electron: the time-domain modulation is the product of single-nucleus
  modulations within each manifold pairing. In the stick picture, for two nuclei
  with peak sets $\{(f^{(1)}_p, g^{(1)}_p, a^{(1)}_p)\}$ and $\{(f^{(2)}_q, g^{(2)}_q, a^{(2)}_q)\}$
  (each set including its zero-frequency "1" term, i.e. `i==k` / `l==n` entries
  that the binner normally discards), the combined set is
  $\{(f^{(1)}_p+f^{(2)}_q,\ g^{(1)}_p+g^{(2)}_q,\ a^{(1)}_p a^{(2)}_q)\}$ per
  orientation, applied separately to the (α→ν₁, β→ν₂) and (β→ν₁, α→ν₂) branches.
  Equivalent nuclei (`nNucs > 1`) are combined with themselves the same way.
  Peak counts multiply, so prune amplitudes below `Treshold·max|a|` after each
  combination. **[AUTHOR DECISION]** confirm the exact normalisation of the
  product-rule expression (½ prefactor per branch) against a reference
  (EasySpin `saffron`) before enabling by default; ship with `ProdRule = False`
  as default until validated.
- `nNucs` is ignored by v0.9's simulation; the rebuild applies it through the
  product rule, or (sum rule) multiplies that nucleus's amplitudes by `nNucs`.

### 6.4 Notes on what is *not* in the model

g-strain, A-strain, finite pulse lengths, relaxation, S > 1/2, ZFS (`D`) and
non-collinear g (`gFrame`) are not modelled. State this in the README; do not
add them unasked.

### 6.5 Field from gField

Wherever `Exp.gField` is set (orientation weights and brute-force selection):
`Exp.Field = Exp.mwFreq / Exp.gField * h_ov_beta_mTMHz * 1e3` (mT), error if
`mwFreq == 0`.

---

## 7. Orientation selection (`compute_weights`, `computeOrisel_eig`)

Common: σ = `ExciteWidth / 2.35482`, ν_mw = `mwFreq·10³` MHz. If
`ExciteWidth ≤ 0` or `mwFreq ≤ 0`: all weights 1. After computing `a(l)`: if
`max(a) < 1e-4` → error "no resonances" (test on the **unweighted** values), then
`ak = a · ww`.

### 7.1 `g_eff`

$f_e(\mathbf l)=B_0\,g_\mathrm{eff}(\mathbf l)/0.0714477$,
$a(\mathbf l)=\exp\!\big[-\tfrac12\big((f_e-\nu_\mathrm{mw})/\sigma\big)^2\big]$.

### 7.2 `g_eff+HFC`

For nuclei with `useFor != 'sim.only'`, each equivalent copy counted, projected
coupling $p_k(\mathbf l)=|A_k^{\mathsf T}\mathbf l|$ (MHz). Sum over all
combinations of nuclear projections $m_k\in\{-I_k,\dots,I_k\}$:

$$
a(\mathbf l)=\sum_{\{m_k\}}\exp\!\Big[-\tfrac12\Big(\frac{f_e(\mathbf l)+\sum_k m_k\,p_k(\mathbf l)-\nu_\mathrm{mw}}{\sigma}\Big)^2\Big].
$$

Vectorise: build the `(n_comb, n_nuc)` matrix of m-values with `np.meshgrid`,
then one `(n_ori, n_comb)` exponent. (v0.9 mis-indents the multiplicity
bookkeeping and uses only the first few columns — do not copy.)

### 7.3 `brute force` (`computeOrisel_eig`)

Full Hilbert space of the electron plus all non-`sim.only` nuclei (each
equivalent copy explicit). Per orientation:

$$
H = \frac{B_0}{0.0714477}\,(g\mathbf l)\cdot\mathbf S \;-\; \sum_k \gamma_kB_0\,\mathbf l\cdot\mathbf I_k \;+\; \sum_k \mathbf I_k^{\mathsf T}A_k\,\mathbf S .
$$

(Use the same nuclear-Zeeman sign as §4.3; v0.9 uses the opposite sign here.)
Diagonalise; for every level pair (i < j) with
$\delta=\big||E_i-E_j|-\nu_\mathrm{mw}\big| < 5\sigma$, add
$|\langle i|\mathbf S\cdot\mathbf b_1|j\rangle|^2\,e^{-\delta^2/2\sigma^2}$ where
$\mathbf b_1$ is the unit vector at (θ+π/2, φ) (⟂ to B₀). v0.9 uses the
magnitude, not its square. **[AUTHOR DECISION]** squared (transition probability,
recommended) vs. magnitude (v0.9).

Performance: batch `eigh` over orientations, compute `V† (S·b₁) V` batched with
`einsum`, and evaluate all pairs as an `(n_ori, D, D)` mask — no double Python
loop. Only S = 1/2 allowed.

Return `N×3 [phi, theta, ak]`; with `setActive` store it in `Opt.OriSelInp`.

### 7.4 `precalculated`

The Opt tab's **"Update Ori.Sel. Grid"** button runs `computeOrisel_eig` once per
dataset (each dataset's own field/frequency), stores the matrix in
`dataset['orisel']`, and from then on that dataset simulates with
`OriSelType = 'precalculated'` using it (saved in the session).

---

## 8. Analysis region geometry

### 8.1 Tilted analysis box (`optHYSCORE.diagonalBox(axis0, axis1, fmin, fmax, antidiagspan)`)

The region every comparison uses (diagonal projection, z-scaling, grid-search
fit) is the 45°-tilted rectangle inscribed in the square window
`[fmin, fmax]²`. Copy this algorithm exactly:

1. `sel0 = (axis0 ≥ fmin) & (axis0 ≤ fmax)`, same for `sel1`; offsets `off0, off1`,
   lengths `n0, n1`. Return `None` if either is empty.
2. Samples on diagonal `d = i − j` inside the window:
   `count(d) = min(n0−1, n1−1+d) − max(0, d) + 1` (non-increasing in |d|).
3. Half-width in diagonals: `nwin = min(n0, n1)`, `step = |axis0[1]−axis0[0]|`.
   If `antidiagspan ≥ 0`: `D = round(antidiagspan/step)`; if `−1` (the only
   sentinel): `D = (nwin−1)//2` (square box, corners at edge midpoints). Clamp
   `0 ≤ D ≤ nwin−1`. (0 means the bare diagonal.)
4. Run length along the diagonal: `L = min(count(−D), count(D))`; `None` if `L < 1`.
5. For `d ∈ [−D, D]`: centred run `i = max(0,d) + (count(d)−L)//2 + arange(L)`,
   `j = i − d`. Every diagonal has exactly `L` samples.
6. `mask[i, j] = True` on an `(n0, n1)` array.
7. Display outline: `s = i + j`; `s_lo = ½(min s[:,0] + max s[:,0])`,
   `s_hi = ½(min s[:,−1] + max s[:,−1])`; corners in polygon order
   `(s_lo,−D), (s_hi,−D), (s_hi,D), (s_lo,D)` mapped to
   `(axis0[off0] + ½(s+d)·step0, axis1[off1] + ½(s−d)·step1)`.

Return `{'sel0','sel1','mask','idx0': off0+i,'idx1': off1+j,'corners',
'ndiag': 2D+1,'nalong': L,'nwindow': nwin}`.

Never rotate/resample the image (a 45° rotation of a lattice leaves half the
pixels empty). Keep this method in one place (v0.9 duplicates it in
`MainFrame`; the rebuild calls `self.Opt.diagonalBox`).

### 8.2 Dataset region for fitting (`MainFrame.dataset_region(dd, fit_noise)`)

From the experimental spectrum (FFT result, or raw data for already-FFT'd files)
and the dataset's `fmin, fmax, antidiagspan`:

```
box   = diagonalBox(expX, expY, fmin, fmax, antidiagspan)
exp   = spectrum[np.ix_(box.sel0, box.sel1)]
null  = Σ exp[mask]²                  # no-simulation SSE
sigma = dd['noise'] if dd['noiselevel'] < 0 else max(spectrum)·dd['noiselevel']
offset= sigma if (fit_noise and sigma) else 0
region = {exp, mask, n = mask.sum(), corners, x = expX[sel0], y = expY[sel1],
          null, null_rms = sqrt(null/n), sigma, offset, field}
```

The region is cut once per scan and is the only experimental data a grid-search
worker receives.

---

## 9. Data handling

### 9.1 Bruker reader (`brukerread.py`)

`brukerread(filename, return_ax=False, return_dsc=False)` → `data` or
`(ax, data[, dsc])`.

- **Pairing:** `.dsc/.dta` (BES3T) or `.par/.spc` (ESP); extension case-insensitive.
  On case-sensitive file systems try both `.DTA/.dta` spellings.
- **DSC/PAR parsing:** one `KEY value` per line, split at the first whitespace;
  skip empty lines and lines starting with any of ``~!@#$%^&*().//\``; a line
  whose first token starts with a digit continues the previous key; strip `'`.
  Read with `encoding='latin-1'` (Bruker files are not UTF-8).
- **Axes (BES3T):** `XPTS, XMIN, XWID` → `x = XMIN + XWID/(XPTS−1)·arange(XPTS)`,
  same for Y; label from `XNAM`/`XUNI`. Keep the ESP `JSS` branches of v0.9
  (2 = CW sweep, 32 = ESP380, 4128/4144 = 2D) for completeness.
- **Complex:** `IKKF == CPLX` (BES3T) or `XQAD == ON` (ESP) → interleaved re/im.
- **Binary:** `BSEQ` `BIG` → `'>'` else `'<'`; `IRFMT` `I` → `int32`, else `float64`;
  `np.fromfile(..., count=N)`; reshape `order='F'` to `(nx, ny)`. Return 2D for 2D
  data (v0.9 returns a trailing singleton axis and every caller slices `[:,:,0]` —
  drop the third axis at the source).
- **Other fields:** `freq1 = MWFQ` (Hz), `cf = CenterField`, `title = TITL`.
- v0.9 bugs not to copy: undefined `file_path` in the not-found error; undefined
  `name` when reading `.XGF/.YGF` axis files (use `ppath + '.XGF'`).

### 9.2 Loading a dataset (`MainFrame.loadData(path, BG=False)`)

1. Read with `brukerread(..., return_ax=True, return_dsc=True)`.
2. `freq = MWFQ·1e-9` (GHz).
3. `tau`: from `FTEzDelay1` (`"136 ns"` → 136) if present, else from a title token
   containing `tau` or `ns` (digits only). `field`: title token with `mT`
   (digits) or `G` (digits/10). Default 0. Use one regex per quantity, e.g.
   `r'(\d+(?:\.\d+)?)\s*mT'`, rather than v0.9's digit filtering (which turns
   "345.6mT" into 3456).
4. Time axes in ns (`'ns' in xlabel`) → ×10⁻³ to µs.
5. `isfft` if `'fft'` in file name or title (case-insensitive), or
   `DSRC == 'MAN' and PROCESS == 'prAbs'` (v0.9 checks the misspelt key `DSCR`, so
   this never fires — fix). For FFT'd data: `fmin, fmax = x[0], x[-1]`; for time
   data: `fmax = 1/(2·dt)`, `fmin = −fmax`.
6. Append the dataset dict (§9.4) to `self.Data` (or `self.BGData`), refresh the
   list and plots. Datasets are dicts — keep them dicts (they are serialised
   wholesale into the session).

### 9.3 `Exp` fields taken from the dataset (`actuateSimMethod`)

Before simulating dataset `dd`, copy the Exp parameter dict and replace every
negative value: `tau ← dd['tau']`, `Field ← dd['field']`, `mwFreq ← dd['freq']`,
`MaxFreq ← dd['fmax']`, `nPoints ← len(dd['fftax']['x'])` (or the raw axis for
FFT'd data). If `dd['orisel']` exists, set `Opt.OriSelInp`. Store the resolved
dicts in `dd['simmethod']`.

### 9.4 Dataset dictionary (complete key list)

```python
{'ax', 'data', 'dsc', 'fname', 'fullpath', 'title',
 'bgdata': 'none' | <BG title>, 'bgscale': 1.0,
 'freq', 'tau', 'field',
 'show': True, 'zmax': 1.0, 'zmin': 0.0,
 'fmin', 'fmax', 'antidiagspan': -1.0,
 'isfft', 'fftdata', 'fftax', 'fftmethod', 'fftactual',
 'noise': 0.0, 'noisesky': 0.0, 'noiselevel': -1.0,
 'simdata', 'simax', 'simmethod', 'simactual',
 'orisel', 'region'}
```

`simax = {'x','y','xlabel','ylabel','orisel': array([phi, theta, ak]), 'oriseltri'}`.

### 9.5 Time-domain processing (`MainFrame.update_FFT(fftmethod=None)`)

Per dataset (FFT parameters are **per dataset**: explicit argument → the
dataset's stored `fftmethod` → the UI controls). Skip transform for `isfft`
datasets but still measure their noise (step 9).

FFT parameters (Data tab, bottom):

| key | label | default | type |
|---|---|---|---|
| `polynomial` | Polynom | 2 | int |
| `apodization` | Apodization | `Hamming` | `Hamming` / `Lor-Gau` / `Gaussian` / `none` |
| `awidth` | A.Width | 1.0 | float (step 0.05) |
| `aalpha` | Lor.Width | 0.6 | float |
| `ashift` | A.Shift | 0.0 | float |
| `zfill` | Zero fill | 2 | int |
| `imag` | Use 'imaginary' | False | bool |

Pipeline on `D` (complex or real, shape nx×ny, axes `X`, `Y` in µs):

1. **Background:** if `bgdata != 'none'`, `D −= bgscale · BG` (shapes must match,
   else error).
2. If not `imag`: `D = D.real`.
3. **Polynomial baseline, both axes:** fit order `p` along `Y` for every row, subtract;
   then along `X` for every column, subtract. **Vectorise**: one
   `np.polynomial.polynomial.polyfit(Y, D.T, p)` call fits all rows at once
   (`polyfit` accepts a 2D `y`), `polyval` evaluates all; same for columns. (v0.9
   loops `np.polyfit` per row — O(n) Python calls.)
4. **Apodization** with `a_x = max(X)·awidth`, `a_y = max(Y)·awidth`,
   `σ² = (a_x/2.354820045)²`:
   - Hamming: $w=(0.54+0.46\cos(\pi t_x/a_x))(0.54+0.46\cos(\pi t_y/a_y))$
   - Gaussian: $w=\exp(-(t_x^2+t_y^2)/\sigma^2)$
   - Lor-Gau: $w=\exp\!\big(-((t_x-s)^2+(t_y-s)^2)/\sigma^2 + t_x/(a_x\alpha) + t_y/(a_y\alpha)\big)$,
     s = `ashift`, α = `aalpha`.
   Build as an outer product of two 1D windows where separable (Hamming, Gaussian)
   — never a full `meshgrid` when not needed.
5. **Zero fill:** pad both axes at the end to $2^{\lceil\log_2 n\rceil+z}$; error if
   the exponent > 13.
6. `F = |fftshift(fft2(D))|`.
7. Axes: `nyq = 1/(2·(X[1]−X[0]))` (MHz), `linspace(−nyq, nyq, N, endpoint=False)`.
   Warn once if Δν₁ ≠ Δν₂ (relative 10⁻⁶).
8. Store `fftdata, fftax, fftactual = True`.
9. **Noise level** (the far corner of the map holds no HYSCORE signal):
   `sliver = F[0:N//4, 0]`; `noise = mean(sliver)`, `noisesky = max(sliver)`.
   These ride along into the session and are the grid search's χ² denominator.

The **auto FFT** toolbar toggle re-runs this on every FFT-parameter change.

---

## 10. Function modifiers (`sys_functions.py`, `classFunctionModPG.py`)

Purpose: drive one Sys parameter from an expression of custom scalars and/or
other Sys parameters (e.g. a dipolar `nuc(1).A(3) = 59/r**3 + a`). Driven
parameters become read-only in the property grid and are honoured inside the
grid search.

### 10.1 Engine (UI-free, thread-safe)

- `SAFE_FUNCS`: `sqrt exp log log10 sin cos tan asin acos atan atan2 pi e abs min max pow`.
- `extract_refs_and_freevars(expr, known_labels)`: take the right-hand side of `=`,
  replace `^` with `**`; replace each literal occurrence of a known label (labels
  look like `nuc(1).A(3)`, `spin(1).S`) by a placeholder token `__PGREF{i}__`,
  **longest label first**; free variables = identifiers (`[A-Za-z_]\w*`) minus
  placeholders minus `SAFE_FUNCS`.
- `topo_sort(deps)`: Kahn's algorithm; returns `(order, cyclic_nodes)`.
- `evaluate_functions(rows, variable_values, prop_value_lookup, known_labels)`:
  1. reject duplicate targets (both rows error);
  2. parse every row; dependency edge when a row references another row's target;
  3. topo-sort; cyclic rows error "circular reference";
  4. evaluate in order with `eval(work, {"__builtins__": {}}, env)` where `env`
     holds `SAFE_FUNCS`, custom variables and resolved references (other driven
     targets from this pass, else `prop_value_lookup`);
  5. return `{'values': {row: float}, 'errors': {row: str}, 'label_values', 'order'}`.
- `UnresolvedReference(LookupError)` is what `prop_value_lookup` raises.

These are module-level functions because both the wx panel and the worker
process call them (≥2 callers in different processes) — the one sanctioned
exception to "methods on classes" besides the pool entry points.

### 10.2 Panel

Rows of `[target choice] = [expression text]`; a variables block with one
spin control per free variable (auto-added/removed as expressions change; step
via right-click menu); on any change `recompute_all()`: evaluate, write values
into the property grid without triggering a simulation per row, mark errors in a
red row colour and print them, set driven properties read-only, store
`parameters['functions'] = {'func': [[label, expr]...], 'var': {...}}`, and fire
one simulation if anything changed.

---

## 11. 2D grid search (`grid_search.py`, `grid_search_worker.py`)

### 11.1 Parameters that can be scanned

`flatten_sys_params(sys_dict)` lists every plain float/int and every element of a
1D ndarray in the Sys dict as
`{'kind': 'sysleaf', 'cat': 'nuc(1)', 'key': 'A', 'idx': 2, 'value': ...}`
(bools, strings, enums skipped). Label: `nuc(1).A(3)` (1-based index in labels).
Function-modifier variables are scannable as `{'kind': 'funcvar', 'key': 'r'}`.
Parameters driven by a function are shown locked (not scannable).

UI (left pane): "Refresh Parameters from Main Window" button; a property grid of
checkboxes (max two ticked); per ticked parameter a row `min | max | points`
(remembered per parameter during the session); options: *subtract exp. noise*
(fit_noise), *Keep simulations in memory*, *parallel, workers: N* (default
`os.cpu_count()−1`), Run/Cancel button, progress gauge with
"count / total, elapsed, ETA". Display controls: **Display**
(RMSD / Fit Quality norm(RMSD) / Probability / log(probability)),
**Smooth** (None / Linear / Cubic), **Oversampling** F (default `4**zfill`, §11.5),
*As Contours* + line thickness, colormap button + *Invert*, contour settings
dialog, "Show/Hide Residual Maps".

### 11.2 Scan set-up (`_start_run` — keep as one method)

1. Exactly two parameters, `max > min`, at least one dataset with `show=True`.
2. `grid1 = linspace(min1, max1, n1)`, `grid2 = ...`.
3. Cut every active dataset's region once (§8.2); abort if any is empty.
4. If storing simulations: estimate `Σ region.size · n1 · n2 · 4 bytes`; above
   0.5 GiB ask Yes (store) / No (don't store) / Cancel.
5. Payload (deep copies, picklable, **no full 2D maps**):
   `sys_template, exp_template, opt_template, func_rows, base_variable_values,
   param1, param2, records, store_layout` where
   `records[k] = {title, tau, field, freq, fmax, nPoints, orisel, region}`.
6. Tasks: `[(i, j, v1, v2) for i, v1 in enumerate(grid1) for j, v2 in enumerate(grid2)]`.
7. Run `worker_run` on a background `threading.Thread` (daemon) so the GUI stays
   live; it drives either a process pool or an in-thread loop. All GUI updates
   from it go through `wx.CallAfter`.

### 11.3 One grid point (`run_point(task)` → `run_gridpoint(...)`)

```
sys_dict = deepcopy(sys_template); vars = dict(base_variable_values)
apply v1, v2 (sysleaf → into sys_dict; funcvar → into vars)
ok = apply_functions_to_sys_dict(sys_dict, func_rows, vars)   # False → NaN point
Sys = sysPar(); Sys.setFromCtrl(sys_dict)
hs.Sys = Sys; hs.preCompute()
for each record k:
    Exp from exp_template with −1 replaced from the record (§9.3)
    Opt from opt_template (+ record.orisel if any)
    hs.reRun()                         # grid cache survives between points
    s  = adaptSim(region.x, region.y, hs.Spectrum, hs.X, hs.Y)
    SSE_k, N_k, scaled = calcSSR(region.exp, s, region.mask, region.offset)
    accumulate ΣSSE, ΣN, Σ N_k·SSE_k/null_k ; keep (SSE_k, N_k) ; store scaled
return (i, j, RMSD, ΣN, RMSD_norm, [(SSE_k, N_k)...], present_flags, error_log)
```

- `adaptSim`: `RegularGridInterpolator((SimX, SimY), Sim, bounds_error=False,
  fill_value=0)` evaluated on the region's mesh (`indexing='ij'`). If the
  simulated axes coincide with the region's (same step and aligned origin), skip
  interpolation and slice directly — the common case when `nPoints`/`MaxFreq`
  follow the data; check with `np.allclose`.
- `calcSSR(Exp, Sim, Mask, offset)`:
  $c=\max\!\big(\langle e-o, s\rangle_\mathrm{mask}/\lVert s\rVert^2_\mathrm{mask},\,0\big)$
  (0 if $\lVert s\rVert^2<10^{-30}$), scaled = $c\,s+o$,
  SSE = $\sum_\mathrm{mask}(e - \text{scaled})^2$, N = mask count.
- Errors never pop dialogs in workers: collect `{message: count}` and return it;
  the coordinator merges, prints, and shows the five most common at the end.

### 11.4 Pool plumbing and shared memory

- Context: **`multiprocessing.get_context('spawn')` on every OS** (the only method
  on Windows; using it everywhere keeps failure modes identical).
- `Pool(processes=N, initializer=init_worker, initargs=(payload,))`;
  `init_worker` builds one `HYSCOREsim` per process (so its grid cache persists
  across points) and attaches to the `SimStore`. Results via
  `imap_unordered(run_point, tasks, chunksize=1)` (point costs vary wildly).
- In each worker set `OMP_NUM_THREADS = MKL_NUM_THREADS = OPENBLAS_NUM_THREADS = 1`
  **before importing NumPy** (set in the parent's environment prior to creating the
  pool) to avoid BLAS over-subscription; in single-process mode leave BLAS threads
  alone.
- Single-process mode (`workers = 1`) runs the *same* `init_worker/run_point`
  in the background thread, then `release_worker()`. A class constant
  `RUN_IN_THREAD = True` can be flipped to run synchronously for debugging.
- **`SimStore`**: one `shared_memory.SharedMemory` block of
  `n1·n2·Σ region.size` float32, each (i, j, k) owning a fixed disjoint slice
  (`base = (i·n2 + j)·stride + offset_k`). Workers write directly; nothing is
  pickled. The owner records which slices were written (`note(i, j, present)`),
  exposes a dict-like `get((i, j))` returning copies, and `dispose()` (drop views,
  `close()`, `unlink()` by the owner — required on POSIX). All stores are disposed
  when the window closes.
- **Cancel:** set an `Event`, `pool.terminate()`, discard the whole scan (a
  partially filled map is misleading); the previous results stay on display.

Coordinator arrays: `rmsd, rmsd_norm (n1,n2)`, `total_n (n1,n2) int`,
`SSE (n1,n2,K)`, `N (n1,n2,K)`.

### 11.5 Statistics (cite `grid_search_statistics.tex` in the docstring)

For dataset k at θ: $e_k$, $s_k$, $N_k$, $\mathrm{SSE}_k$, $\mathrm{null}_k=\lVert e_k\rVert^2$.

1. **RMSD**: $\sqrt{\sum_k\mathrm{SSE}_k/\sum_kN_k}$ (dominated by bright datasets).
2. **Fit Quality, norm(RMSD)**: $\sqrt{\sum_kN_k(\mathrm{SSE}_k/\mathrm{null}_k)/\sum_kN_k}\in(0,1]$;
   equals $\sqrt{1-r_k^2}$ with $r_k$ the uncentred correlation (single dataset).
3. **Posterior** (Gaussian residuals, unknown per-dataset σ_k, Jeffreys prior
   $1/\sigma_k$, marginalised):
   $$p(\theta\mid D)\propto\prod_k \mathrm{SSE}_k(\theta)^{-N^\mathrm{eff}_k/2},\qquad N^\mathrm{eff}_k=N_k/F,\quad F=4^{z}$$
   with z the zero-fill exponent (each axis gets 2^z points, so the 2D count
   gains 4^z). Compute in log space:
   `logp = −(N/F/2)·log(SSE)` summed over k, minus its max.
   Brightness-invariant and normalisation-invariant by construction.
4. **Displays 3/4**: interpolate **log p** (smooth, near-quadratic) — not p — onto
   the fine grid, then `p = exp(logp)` normalised to Σ = 1 (failed points → 0).
   - *Probability (credible level)*: $L_{ij}=\sum_{p_{ab}\ge p_{ij}}p_{ab}$
     (`credible_level_map`: argsort descending, cumsum, ties share the level of
     the end of their block via `searchsorted`). Contours at 0.68 / 0.95 *are*
     the HPD regions.
   - *Probability (density)*: $p_{ij}$ itself; contour at
     `density_threshold(p, conf) = interp(conf, levels_sorted, p_sorted)`.
   - *log(probability)*: `logp − max`, contour at the same density threshold in
     log units.
5. **Intervals**: marginals $p^{(1)}_i=\sum_jp_{ij}$, $p^{(2)}_j=\sum_ip_{ij}$,
   95 % highest-density interval (`marginal_interval`: sort descending, take the
   smallest set reaching 0.95, report its span). Also report a Gaussian
   approximation: fit a parabola to log p through ±5 fine-grid cells around the
   maximum along each axis, $x_0=-b/2a$, $\sigma=\sqrt{1/(2|a|)}$, half-width
   $z_{0.975}\sigma$ (`scipy.special.ndtri`). Show "label = x₀ ± hw (95 %)" for both
   parameters.
6. **Edge warning**: if the 95 % region touches the scanned box edge, say so (the
   box is the prior).
7. For RMSD displays, **F-test** contour (unknown σ): threshold
   $\mathrm{SSE}\le\mathrm{SSE}_\mathrm{min}(1+\tfrac{p}{N-p}F_{p,N-p}(\mathrm{conf}))$,
   p = 2, drawn in RMSD units ($\mathrm{RMSD}_\mathrm{min}\sqrt{1+\Delta}$); and,
   if every dataset has a measured σ_k, the known-noise χ² variant
   $\Delta\chi^2_2(\mathrm{conf})$ with pooled
   $\sigma^2=\sum_kN_k\sigma_k^2/\sum_kN_k$ (normalised mode: σ_k/null_rms_k). Report
   reduced χ² of the best point = $\mathrm{SSE}_\mathrm{min}/(\sigma^2(N^\mathrm{eff}-2))$
   so the user sees whether the residual reached the noise.
8. **Interpolation** (`_interpolate_grid` → name it `interpolate_grid`, it has
   several callers): `none` (raw grid), `linear` (`RegularGridInterpolator`),
   `cubic` (`RectBivariateSpline`, only if ≥4×4 and no failed points, else
   linear). Fine grid ≥ 160 per axis. Failed (NaN) points are filled with
   2·max+ε for interpolation only.

### 11.6 Heat-map interaction

- `imshow` of the map with cell-centred extent, colorbar (remove the previous
  colorbar **before** `ax.clear()` and restore the stored axes position, or the
  axes shrink on every redraw), dotted rules + dots at the actually-scanned
  points, optional background contours, labelled confidence contours, ★ at the
  best fit, distinct markers for "residuals shown here" and "main window's
  current Sys" (read live from the main window).
- Hover → status bar `p1 = …, p2 = …, value = …`.
- Right-click → "Load these values into main Sys grid" (sysleaf → parameter;
  funcvar → `functions.var`, then rebuild the panel so driven values recompute).
- Double-click (when residual row visible and sims stored) → nearest grid point:
  per dataset show residual `exp − scaled_sim` (NaN outside the box), normalised by
  `null_rms`, plus sim contours, titled `B₀ = … mT`, with each map's RMSD.

---

## 12. Performance requirements

### 12.1 Hot path rules (simulation and grid point)

1. **No Python loop over orientations.** Batch `eigh`, matrix products and
   `einsum` along a leading orientation axis.
2. **Use the d⁴ contraction (§4.5.1)**, never the d⁶ index expansion.
3. **Bin with `np.bincount`** (real and imaginary separately), never `np.add.at`,
   never `np.concatenate` inside a loop — preallocate or concatenate once.
4. **Cache everything that depends only on Sys** (`preCompute`): operators,
   rotated tensors, quadrupole operator. Cache the orientation grid by signature
   (§5.6). Cache the broadening kernel by `(N, MaxFreq, lw, KillNeg)`.
5. **Analytic Gaussian kernel** in the time domain (saves one N² FFT).
6. **Threshold** orientations by weight (`Opt.Treshold`).
7. Use `complex128` for amplitudes (phases accumulate error), `float32` only for
   stored display data (`SimStore`).
8. Avoid `np.meshgrid` where an outer product or broadcasting does.
9. Skip `adaptSim` interpolation when grids coincide (§11.3).
10. Optional, behind a single constant, off by default: `scipy.fft` with
    `workers=-1` for the 2D FFTs (fallback `numpy.fft`). No numba/Cython.

### 12.2 GUI responsiveness

- Run the main-window simulation on a worker thread (v0.9 runs it on the GUI
  thread and the window freezes); disable SIM/FFT buttons while running; progress
  via `wx.CallAfter`; ignore re-entrant requests (queue at most one pending
  re-run, coalescing rapid parameter edits — "debounce auto SIM" ~250 ms).
- Redraw: rebuild subplots only when the layout (number of datasets × rows)
  changes; otherwise update `AxesImage.set_data/set_extent/set_clim`. Contours
  cannot be updated in place: remove and redraw them, and contour only the
  visible window (+1 cell) — a full 2048² contour takes seconds.
- Debounce `canvas.draw` during window resize / sash drags
  (`theme.CanvasRedrawDebouncer`), freeze the frame during live resize
  (`theme.ResizeFreezeGuard`).

### 12.3 Targets (to be measured and reported by the agent in `tests/bench.py`; not hard gates)

| Case | Target (one modern laptop core) |
|---|---|
| 1 × ¹H, `sphgrid` K = 20 (761 ori), N = 256, g_eff | < 0.15 s |
| 1 × ¹⁴N (I = 1) + 1 × ¹H, 761 ori, N = 512 | < 0.5 s |
| 1 × ⁶³Cu (I = 3/2), 761 ori, N = 512 | < 1.5 s |
| Grid search 21 × 21, 3 datasets, the ¹H case, 8 workers | < 30 s |

Report v0.9 vs new timings for the same inputs.

---

## 13. Cross-platform requirements

1. Paths: `os.path.join` / `pathlib` only. Icons and other assets resolved from
   `os.path.dirname(os.path.abspath(__file__))`, never the CWD (v0.9 requires
   starting from the project folder).
2. Per-user settings (last directory, theme, window size) in a single INI under
   `wx.StandardPaths.Get().GetUserConfigDir()` (falls back to `~`), not a hidden
   file sprinkled in `~`.
3. Text files: always pass `encoding=` — sessions `'utf-8'`, Bruker DSC
   `'latin-1'`. (v0.9 writes session XML with the platform default encoding:
   breaks on Windows for non-ASCII titles.)
4. `multiprocessing`: `spawn` context, entry point guarded by
   `if __name__ == "__main__": mp.freeze_support(); main()`; worker functions are
   top-level in a wx-free module; payload picklable; `SharedMemory.unlink()` by the
   owner.
5. No `os.startfile`, no `\\` literals, no Windows-only APIs, no shell calls.
6. macOS: set `wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES` before file dialogs; use
   `Cmd` for accelerators via `wx.ACCEL_CMD` (maps to Ctrl elsewhere). Theme
   toggle `Ctrl/Cmd+D`.
7. HiDPI: draw icons via `wx.BitmapBundle` or scale by `GetContentScaleFactor()`;
   never hard-code pixel sizes for fonts.
8. Owner-drawn widgets must not rely on MSW-specific colour behaviour (v0.9 had to
   replace `wx.StatusBar` with a custom strip because MSW ignored text colours in
   dark mode — keep the custom `StatusStrip`).
9. Case-sensitive file systems: Bruker pairs may be `.DSC/.DTA` or `.dsc/.dta`.
10. CI matrix (if a CI is added): `windows-latest`, `macos-latest`, `ubuntu-latest`
    × Python 3.10/3.12, running only the wx-free tests.

---

## 14. GUI specification (`pyHYSCORE.py`)

### 14.1 Main window layout

```
┌ toolbar: [folder][up] | [save session] | [FFT][auto FFT] | [SIM][auto SIM] | [grid search] ……… [grid][map/contour][colormap][theme] ┐
│◂│ FileBrowserPanel │ MatplotlibPanel (checkboxes row, canvas, nav toolbar) │ TabulatedPanel (Data | Sys | Exp | Opt) │▸│
└ status strip: [message] [x, y readout] [progress gauge] ┘
```

- Two nested `wx.SplitterWindow`s; collapse chevrons outside the splitters for
  the left and right panes.
- Toolbar icons are drawn in code (`theme.IconButton`, `IconToggleButton`) — no
  asset dependency.

### 14.2 File browser

Editable path box (Enter/blur navigates), filter `All | DSC | XML`, list with
columns *name* and *title* (title read from `TITL` in each DSC — read lazily and
cache by mtime), sortable columns, double-click: DSC → load dataset, XML → load
session, folder → enter. "Load as background" button. Remembers last directory.

### 14.3 Tabs

- **Data**: per-dataset category in a property grid —
  *Show, f.name, MW Freq [GHz], B0 [mT], tau [ns], z.max, z.min, f.max, f.min,
  adiag span [MHz], noise/max(y), BG Data (enum), BG Scale*; buttons delete,
  show DSC (key/value dialog), move up/down, expand all; FFT parameter block.
  Editing *noise/max(y)*: −1 = measured; 0…1 = fraction of the peak in the box.
- **Sys**: property grid from the Sys dict (§3.1), "Add Nuc" (isotope chooser
  listing label and I), delete nucleus (v0.9 stub — implement), Function
  Modifiers panel below.
- **Exp**: property grid of §3.2 fields; tooltips explain −1.
- **Opt**: property grid of §3.3; *Target Grid Points* spin + read-only "actual"
  box (§5.5); "Update Ori.Sel. Grid" (§7.4); contour levels / line width /
  overlay colour; *View pad, MHz* + *Fit box*; quadrant radio buttons
  (All 4 / ±ν₁ / ±ν₂).

All float fields use the spin editor with a right-click step menu (0.001 … 100)
that commits on Enter, arrows and focus loss.

### 14.4 Plot panel (`update_graph`)

One column per shown dataset; rows (top → bottom) when enabled:
**data**, **simulation**, **skyline**, **orientation map (3D)**,
**diagonal projection**. Checkboxes: FFT(y), Background, Skyline, Orientat maps,
Diag.Projection, Overlay Sim.

- Data row: FFT spectrum (or raw time data if FFT is off/not done). Analysis
  window `[fmin, fmax]²`, cut by quadrant mode; view window = analysis window
  widened by `viewPad` (clamped to data extent, never crossing a quadrant cut).
- Colour scale: `ref = max over the tilted box` (fallback: window, then whole
  map); `scale = ref·|zmax|` (zmax ≠ 0, else ref); `vmax = scale`,
  `vmin = scale·zmin`. Everything (colour, projection normalisation, noise floor)
  uses the analysis window only — `viewPad` only adds surroundings.
- Fit box: grey wash (α = 0.18) over the view **outside** the tilted box (a path
  with outer rectangle + inner polygon of opposite winding — check the shoelace
  sign and reverse if needed), dashed white outline.
- Density mode: `imshow(origin='lower', interpolation='nearest')`; contour mode:
  colormap contours with `contourLevels` linear levels between vmin and vmax,
  `set_aspect('equal')`.
- Sim row: simulation normalised to max = 1, scale `[zmin/|zmax|, 1]`; "no
  simulation"/"simulation failed" annotation. Overlay Sim: contours of the
  normalised sim at `levels/vmax` in `overlayColor` on both the data and sim
  panels.
- Noise handling (right-click on any panel): *none* | *lift* (sim → sim·(1−u)+u,
  u = noise/vmax) | *subtract* (exp projection − noise floor, clipped at 0); toggle
  "Show noise floor" on projection panels. Noise floor for a max over n samples
  drawn at `sqrt(ln n)·u`.
- Skyline: max of the box-masked map along ν₂ (projected on ν₁) — or along ν₁
  in "±ν₂" mode; shares the x axis with the data panel when it projects onto it;
  exp in blue, sim dashed red; grey wash over the view-pad margins.
- Diagonal projection: max of `data[idx0, idx1]` along each diagonal (axis 1 of
  the box index arrays), plotted against the centred anti-diagonal offset
  `linspace(xw[0], xw[-1], ndiag) − mid`; sim likewise.
- Orientation map: `plot_trisurf` of the unit sphere points with the engine's
  triangulation (`simax['oriseltri']`; fallback `mtri.Triangulation(x, y)` for
  precalculated grids), faces and edges coloured by per-triangle max weight in
  the current colormap.
- Title `B₀ = {field} mT`; mouse-move → `x = …, y = …` in status field 1.

### 14.5 Simulation orchestration (`runSim`)

1. `Sys.setFromCtrl(Sys_param.parameters)`; one `HYSCOREsim`; `preCompute()` once.
2. For each dataset: `actuateSimMethod` (§9.3), `reRun()`, store `simdata`,
   `simax`, `simactual = True`; on `HYSCOREError` store `simdata = None`, record
   the title, continue. Progress gauge steps per dataset regardless of outcome.
3. Status: "Ready. It took X s" or "Simulation failed for n/m dataset(s)".
4. Auto SIM: run on Sys changes (debounced, §12.2). Changing `OriSelType`
   recomputes precalculated grids if selected. Changing Grid/nKnots updates the
   grid-size info.

### 14.6 Theme and colormaps

Light/dark palettes as token dicts in `theme.py` (`bg, ctrl_bg, text, accent,
pill_load, pill_warn, pill_calm, tab_inactive, ax_bg, ax_fg, propgrid_caption_fg,
…`), applied recursively to windows and to matplotlib figures
(`theme.style_figure`). Colormap = list of `(position, '#rrggbb')` stops →
`LinearSegmentedColormap`; editor dialog with draggable stops, undo, presets
from matplotlib names; default `jet`, 9 stops.

---

## 15. Session files (XML)

### 15.1 Content

```python
{'Data': [dataset dicts], 'BGData': [...], 'fftmethod': {...},
 'Sys': Sys_param.parameters, 'Exp': ..., 'Opt': ...,
 'cmap_name': str, 'cmap_stops': [(pos, hex), ...],
 'view': {mapMode, contourLevels, contourWidth, overlayColor, showGrid, ...}}
```

### 15.2 Encoding (`save_mixed_dict_xml` / `load_mixed_dict_xml`)

Root `<dictionary>`; every value is `<item key="…">` containing exactly one of:

- `<numpy_array dtype="float64" shape="(256, 256)"><data>BASE64(tobytes())</data></numpy_array>`
- `<dict>` of items, `<list type="list|tuple">` of items,
- `<primitive type="int|float|float64|bool|str|NoneType">text</primitive>`.

Rules for the rebuild:
- Parse `shape` with **`ast.literal_eval`** (v0.9 uses `eval` on file content).
- Write with `encoding='utf-8'` and an XML declaration.
- `wx.Colour` values are stored as `#rrggbbaa` strings and restored.
- Loading lays saved Opt values **over current defaults** (`withDefaults`): enum
  entries keep the saved selection but take the choice list from the defaults; a
  selection no longer offered falls back to the default. Do the same for Exp and
  Sys leaves so older sessions gain new parameters.
- Loading replaces the current state: **ask for confirmation** if datasets are
  loaded (v0.9 does not). Disable auto SIM while restoring; afterwards offer to
  re-run the simulation if the session contained one.
- Unknown keys are ignored with a console note; a failed array is `None` with a
  console note (never abort the whole load).

---

## 16. Known v0.9 defects — do not replicate

| # | Location (v0.9) | Defect | Required behaviour |
|---|---|---|---|
| 1 | `SysPar.setNucs` | Q_eta zeroed when given / NameError when absent; non-EasySpin sign | §3.1 |
| 2 | `hyscore_sim.compute_weights` (`g_eff+HFC`) | multiplicity list built outside the per-nucleus loop | §7.2 |
| 3 | `build_spin_operators` | `Qpa[ii].np.pi` typo in 3×3 Q branch | §4.1 |
| 4 | `compute_transitions` | per-orientation Python loop, d⁶ expansion, `concatenate` in loop, `np.add.at` | §6.2, §12.1 |
| 5 | `optHYSCORE` | `ProdRule`, `AmpRatios`, `Treshold`, `Symmetry` unused; `nNucs` ignored in simulation | §3.3, §6.3 |
| 6 | `computeOrisel_eig` | nuclear Zeeman sign opposite to simulation; |m.e.| not squared; double Python loop | §7.3 |
| 7 | `tDead` | used in µs while documented nowhere | ns, §2 |
| 8 | `brukerread` | undefined `file_path`, undefined `name`; 3D return for 2D data | §9.1 |
| 9 | `loadData` | `DSCR` vs `DSRC` key; digit-filter parsing of field/tau | §9.2 |
| 10 | `update_FFT` | per-row `polyfit` loops; full meshgrid windows; axis step from `X[2]-X[1]` | §9.5 |
| 11 | `MainFrame.calcRMSD` | broken stub referencing undefined names | remove |
| 12 | `MainFrame.diagonalBox` | duplicate of `optHYSCORE.diagonalBox` | single copy |
| 13 | `runSim` | runs on GUI thread, window frozen | §12.2 |
| 14 | `RaiseError` | shows modal dialog then raises bare `Exception`, used as engine callback | §6.1 |
| 15 | session XML | `eval` on file content; no encoding on write; no overwrite warning | §15 |
| 16 | README | claims arbitrary S; icons need CWD = project folder | §6.1, §13 |
| 17 | `useFor` | `'both'` vs `'all'` | `'all'` only |
| 18 | `on_delete_nuc` | stub | implement |

---

## 17. Testing (`tests/`, pytest, wx-free modules only)

The author tests the GUI personally; automated tests cover the engine.

1. `SpinOp`: commutators $[S_x,S_y]=iS_z$, $S^2=s(s+1)\mathbb 1$ for s = ½…9/2.
2. `RotMatrix` orthogonal, det = 1; `RotTransDiag(a, 0) = diag(a)`.
3. Grids: `nGridPoints == len(phi)` for K = 2…60, all methods; weights sum to 1;
   triangle areas sum to 4π; `sphgrid` relative weights match EasySpin's for K = 5
   (hard-code the reference vector).
4. `nKnotsForGridPoints` returns the closest rung (brute-force comparison).
5. I = ½ frequencies equal $|\mp\tfrac12A^{\mathsf T}\mathbf n-\nu_I\mathbf l|$ for random
   orientations/tensors.
6. Contracted amplitudes (§4.5.1) equal the brute d⁶ sum to 1e-10 for I = ½ and
   I = 1 with random Hamiltonians.
7. Single ¹H, isotropic A < 2ν_I, no orientation selection: peaks at
   (ν_α, ν_β) = (ν_I + A/2, |ν_I − A/2|) and its mirror in the (+,+) quadrant, within one bin.
8. `bincount` binning equals the `np.add.at` reference.
9. Grid cache: changing `tau`, `lw`, `nPoints` does **not** rebuild; changing `g`,
   `nKnots`, a non-`sim.only` `A` does.
10. `diagonalBox`: every diagonal has exactly `nalong` samples; all selected
    samples inside `[fmin, fmax]²`; `antidiagspan = 0` gives one diagonal;
    `= fmax−fmin` gives `nalong = 1`.
11. `calcSSR`: SSE ≤ null; c ≥ 0; scaling `e` by a scales SSE by a².
12. Posterior invariance: multiplying one dataset's `e` by 250 leaves the
    normalised posterior unchanged (to 1e-12).
13. `credible_level_map`: monotone in density; ties share a level; max = 1.
14. Function engine: dependency order, duplicate target, cycle, undefined
    variable, unresolved reference, `^` → `**`, label containing `(`/`.`.
15. Session round trip: save → load reproduces arrays (dtype, shape), nested
    dicts, tuples vs lists, None, bools.
16. `run_point` in-process equals pool result for a 3 × 3 scan (determinism).
17. `tests/bench.py` prints the §12.3 timings (not asserted).

---

## 18. Implementation order for the agent

Each phase ends with its tests green and a short report to the author
(what was built, what deviates from v0.9 and why, open **[AUTHOR DECISION]**
items). Do not start the next phase before the author confirms.

1. `SysPar.py`, `mathfunctions.py`, `sys_functions.py`, `brukerread.py` + tests 1, 2, 14.
2. `hyscore_sim.py`: grids, triangulation, cache, weights (§5, §7.1–7.2) + tests 3, 4, 9.
3. `hyscore_sim.py`: operators, batched transitions, contraction, binning,
   broadening (§4, §6) + tests 5–8; benchmark vs v0.9.
4. `computeOrisel_eig` batched (§7.3); product rule (§6.3) behind `ProdRule`.
5. Region geometry and fitting maths (§8, §11.3, §11.5) + tests 10–13.
6. `grid_search_worker.py` + `SimStore` + test 16; benchmark.
7. GUI: main frame, file browser, Data/Sys/Exp/Opt tabs, FFT pipeline, plot panel,
   theme, colormap editor (§9, §14).
8. Grid-search window (§11.1–11.6).
9. Session XML (§15) + test 15; README with the module map, units table and
   model limitations (§6.4).

---

## 19. Glossary

- **HYSCORE** — Hyperfine Sublevel Correlation spectroscopy: 2D four-pulse ESEEM
  correlating nuclear frequencies of the α and β electron-spin manifolds.
- **Orientation selection** — only molecules whose EPR line falls within the
  microwave excitation bandwidth contribute; modelled by the weights `ak`.
- **Blind spots** — τ-dependent suppression of cross peaks (the factor *b*, §4.5).
- **Skyline projection** — maximum (not sum) along one axis of the 2D map.
- **Tilted box / analysis box** — §8.1; the only region ever fitted.
- **Credible level / HPD** — §11.5 item 4.
- **Oversampling F** — spectrum points per independent measurement, 4^zfill.
