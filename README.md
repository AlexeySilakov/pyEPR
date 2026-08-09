# pyHYSCORE

A desktop application for simulating and fitting **HYSCORE** (Hyperfine Sublevel
Correlation) spectra, and for comparing those simulations against experimental
Bruker data.

pyHYSCORE reads Bruker 2D datasets, does the time-domain processing (baseline,
apodization, zero-filling, FFT) in the app, simulates the spectrum from an
EasySpin-style spin-system description, and overlays the two. A built-in 2D grid
search scans any pair of spin-Hamiltonian parameters and maps the fit quality,
with confidence intervals and Bayesian credible regions.

Written in Python with wxPython + matplotlib; version **0.9**.

---

## Features

**Data handling**

- Loads Bruker 2D datasets (`.DSC`/`.DTA`) directly from a built-in file browser
  with a live directory listing and dataset titles.
- Recognizes already-FFT'd files (by name, title, or DSC processing flags) and
  skips re-transforming them.
- Pulls τ, field, and microwave frequency out of the DSC/title where possible, so
  `Exp` parameters left at `-1` are filled in per dataset.
- Separate background datasets with a per-dataset scale factor.
- Any number of datasets loaded at once; each gets its own column of plots, and
  each can be shown/hidden, reordered, or deleted.

**Processing**

- Polynomial baseline correction, apodization (Hamming, Lorentz–Gauss, Gaussian)
  with width/shift/α controls, zero-filling, optional use of the imaginary part.
- FFT parameters are per-dataset — change them for one dataset without touching
  the others — with an "auto FFT" toggle that re-transforms on every change.
- A noise level is measured from every spectrum (far corner of the ±/± quadrant,
  where a HYSCORE map holds no signal) and stored with the dataset, so it
  survives a session save/load and can be used later as a χ² denominator.

**Simulation**

- Arbitrary electron spin `S`, any number of nuclei from a built-in isotope
  table, hyperfine tensors with Euler angles, nuclear quadrupole (`K`/η form),
  g-anisotropy, linewidth.
- Orientation grids: `fibonacci` (golden-angle, equal-area), `sphgrid`
  (SOPHE/triangular, EasySpin `Ci`), `spiral` (latitude rings). Since `nKnots`
  means something different in each, a **Target Grid Points** box lets the
  resolution be set by the one figure that is comparable across methods.
- Orientation selection: `g_eff`, `g_eff+HFC`, `brute force` (full diagonalization
  at every grid point), or a `precalculated` grid loaded per dataset.
- Product rule for combination frequencies, suppression of negative-time signal,
  amplitude ratios between nuclei, per-nucleus `useFor` flags
  (`sim.only` / `ori.sel.only` / `all`).
- The orientation grid is cached and only rebuilt when something it actually
  depends on changes (see `HYSCOREsim.grid_signature`).

**Display**

- 2D contour maps of data and simulation, per dataset, with quadrant selection
  (all four, ±ν₁, ±ν₂).
- Simulation overlay on the experimental map, skyline (max) projections,
  diagonal projections, orientation-selection maps on the sphere.
- A "view pad" that widens what is *drawn* around the f.min…f.max analysis window
  without changing the analysis window itself — colour scale, projection
  normalization, noise floor and fit region all stay on f.min…f.max.
- Projection panels can lift or subtract the experimental noise floor before
  comparing with a noiseless simulation (right-click menu).
- Light/dark theme (Ctrl+D) and a full colormap editor.

**Fitting**

- **Sys 2D Grid Search**: scan two scalar `Sys` parameters over a user-defined
  grid, re-simulating every active dataset at every point. Fitness is the pooled
  RMSD with a per-dataset least-squares amplitude, and the map can be displayed
  as raw RMSD, normalized (0–1) RMSD, Bayesian credible level, or posterior
  density.
- Confidence intervals at 68% and 95% from an F-test, a naive χ², or the
  measured noise level; reduced χ² of the best fit is reported so it is obvious
  whether the residual has reached the noise or stalled above it.
- Optionally keeps every grid point's simulation in memory (interpolated onto the
  ROI only) so double-clicking any point on the finished map shows that point's
  per-dataset residual maps without re-simulating.
- **Function modifiers**: drive one `Sys` parameter from an expression of custom
  scalars and/or other parameters (e.g. a dipolar `A` from `r`), evaluated in
  dependency order with cycle detection. Driven parameters become read-only, and
  the expressions are honoured inside the grid search too.

**Sessions**

- The whole session — datasets, FFT settings, spin system, display options — is
  saved to and restored from XML.

---

## Requirements

- Python 3.9+
- [wxPython](https://wxpython.org/)
- NumPy
- SciPy
- matplotlib

```bash
pip install wxPython numpy scipy matplotlib
```

## Running

Run from the project directory — the file browser loads its icons from `./icons`
relative to the working directory:

```bash
python pyHYSCORE.py
```

## Quick start

1. **Pick a folder.** Use the folder button in the toolbar (or type a path into
   the box at the top of the file browser). The `Filter` box narrows the listing
   to `DSC`, `XML`, or `All`.
2. **Load data.** Double-click a dataset to load it. Loaded datasets appear in
   the **Data** tab and get a column of plots.
3. **Transform.** Set the FFT parameters at the bottom of the **Data** tab and
   press the `FFT` toolbar button — or leave `auto FFT` on and they re-transform
   as you type.
4. **Describe the spin system.** On the **Sys** tab, set `S`, `g` and `lw`, then
   `Add Nuc` to add nuclei and set `A`, `Apa`, `Q_K`, `Q_eta`.
5. **Set the experiment.** On the **Exp** tab, anything left at `-1` is taken
   from the loaded dataset (field, τ, mwFreq, MaxFreq, nPoints).
6. **Choose the grid.** On the **Opt** tab, pick a `Grid` method and set either
   `nKnots` or the **Target Grid Points** box. Start coarse.
7. **Simulate.** Press the `SIM` toolbar button (or leave `auto SIM` on). Tick
   `Overlay Sim` to draw the simulation as contours over the data.
8. **Fit.** Open the grid-search window from the toolbar, tick two parameters,
   set ranges, and run.
9. **Save.** The save button writes the whole session to XML; loading that XML
   from the file browser restores it.

---

## The tabs

| Tab | What lives there |
| --- | --- |
| **Data** | Loaded datasets — show/hide, reorder, delete, background file and scale, per-dataset f.min/f.max, DSC viewer, and the FFT parameter block. |
| **Sys** | Spin system: electron spin(s) and nuclei, plus the Function Modifiers panel. |
| **Exp** | Experimental settings: `mwFreq`, `Field`/`gField`, `tau`, `nPoints`, `MaxFreq`, `ExciteWidth`, `tDead`. |
| **Opt** | Simulation options: `Grid`, `nKnots`, Target Grid Points, `OriSelType`, `Symmetry`, `ProdRule`, `KillNeg`, `Treshold`, `AmpRatios`. |

Every parameter has a tooltip explaining its units and meaning; the `nKnots` and
`Grid` tooltips in particular spell out how the three grid methods differ.

## FFT parameters

| Parameter | Meaning |
| --- | --- |
| `Polynom` | Order of the polynomial baseline subtracted before apodization. |
| `Apodization` | `Hamming`, `Lor-Gau`, or `Gaussian`. |
| `A.Width` | Apodization window width. |
| `Lor.Width` | Lorentzian α for the Lorentz–Gauss window. |
| `A.Shift` | Shift of the apodization window. |
| `Zero fill` | Zero-filling exponent — the transform is 2ⁿ× longer per axis. |
| `Use 'imaginary'` | Include the imaginary channel in the transform. |

Note that zero-filling correlates adjacent spectrum points, which is why the grid
search's **Oversampling** divisor defaults to `4**zfill` (2ⁿ per axis, squared for
a 2D point count).

## Grid search statistics

The four **Display** modes are not cosmetic variants of each other:

- **RMSD** — pooled residual over all active datasets, each with its own
  non-negative best-fit amplitude.
- **Fit Quality, norm(RMSD)** — each dataset's SSE divided by its own null SSE,
  giving an absolute 0–1 score (0 = perfect, 1 = the simulation explains
  nothing). Equivalent to √(1−r²) with r the uncentred correlation.
- **Probability (credible level)** — every cell labelled with the smallest
  credible region containing it, so the 0.95 contour *is* the 95% region.
- **Probability (density)** — the posterior density itself.

The posterior assumes Gaussian residuals with a separate unknown noise level per
dataset and integrates those levels out, which makes it immune to dataset
brightness and to the normalization choice. The full derivation, with references,
is in [`grid_search_statistics.tex`](grid_search_statistics.tex); the module
docstring at the top of [`grid_search.py`](grid_search.py) is the working
summary.

---

## Module map

| File | Role |
| --- | --- |
| `pyHYSCORE.py` | Main application: file browser, plot panel, tabs, session XML, simulation orchestration. |
| `hyscore_sim.py` | Simulation engine — `optHYSCORE` (options) and `HYSCOREsim` (grids, spin operators, transitions, binning, orientation selection). |
| `SysPar.py` | `sysPar` / `expPar` — EasySpin-inspired spin-system and experiment parameter containers, plus the isotope table. |
| `grid_search.py` | The 2D grid-search window: scan worker, fitness map, residual maps, confidence intervals. |
| `sys_functions.py` | UI-free expression engine behind the function modifiers (safe `eval`, dependency ordering, cycle detection). |
| `classFunctionModPG.py` | wx front end for the above, embedded in the Sys tab. |
| `classPropGridPanel.py` | Property-grid panel and the spin-with-step-menu float editor used throughout. |
| `theme.py` | Light/dark theme, pill widgets, themed matplotlib toolbar, redraw debouncing. |
| `colormap_editor.py` | Colormap stop editor. |
| `brukerread.py` | Bruker reader — BES3T (`.DSC`/`.DTA`) and ESP (`.par`/`.spc`). |
| `mathfunctions.py` | Rotation matrices, tensor transforms, spin operators. |
| `pyEPR.py`, `EPR_sim.py` | Sibling CW-EPR application sharing `SysPar`, `brukerread` and the property-grid widgets. |
| `xeprpar.py`, `eprcsvread.py`, `periodic_dialog.py`, `binHYSCORE.py` | Standalone helpers/reference implementations, not imported by the main app. |

## Notes and limitations

- Windows is the primary development platform; nothing in the code is
  Windows-specific, but the theming has only been exercised on MSW.
- The simulation runs on the main thread — the window is unresponsive during a
  run, with progress shown in the status bar. The grid search does run on a
  worker thread and can be cancelled.
- `loadData` currently accepts `.DSC`/`.DTA` only, though `brukerread` can also
  read `.par`/`.spc`.
- The file browser remembers its last directory in `~/.filebrowser.ini`.
- Loading a session XML replaces the current state without warning.

Further items are listed in the `TO DO` block at the top of `pyHYSCORE.py`.

## License

BSD 3-Clause — see [`LICENSE`](LICENSE). The same license used by NumPy, SciPy
and matplotlib: use, modify and redistribute freely, including commercially,
provided the copyright notice is kept and the author's name is not used to
endorse derived products.

## Authors

- Alexey Silakov
- Megan Lavigne

If pyHYSCORE contributes to work you publish, a citation is appreciated.
