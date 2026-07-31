"""
2D grid-search / fitness-mapping window for pyHYSCORE.

Scans two scalar Sys parameters over a user-defined grid. At every grid
point, a full HYSCORE simulation is re-run for every "active" dataset
(active == the dataset's "Show" checkbox in the main window's Data tab)
and compared against that dataset's experimental spectrum, restricted to
the [fmin, fmax] x [fmin, fmax] region shown in the main window's 2D
plots (i.e. the same square region regardless of the "Quadrants" radio
button). Comparison is always done in the frequency domain (FFT'd data
for time-domain files, raw data for files already loaded as FFT spectra),
independent of the main window's own "FFT(y)" display toggle.

Fitness = pooled RMSD:
  For each active dataset k, the simulated spectrum is amplitude-scaled
  by a non-negative factor c_k that minimizes ||exp_k - c_k*sim_k||^2
  over its ROI (closed-form least squares). The residuals of all active
  datasets are then pooled together:

      RMSD = sqrt( sum_k sum_roi (exp_k - c_k*sim_k)^2  /  sum_k N_k )

  i.e. each dataset gets its own best-fit amplitude, but the reported
  number is the RMSD of the combined residual pool, not an average of
  per-dataset RMSDs.

Once a scan finishes, the coarse RMSD grid is interpolated (linear or
cubic, selectable live without re-scanning -- see _interpolate_grid) onto
a finer grid for display, and used to estimate a confidence interval for
each scanned parameter (F-test or naive chi-square, selectable -- see
_rmsd_threshold_factor) at both 68% and 95%, reported as the bounding box
of the joint 2-parameter confidence region projected onto each axis.
"""
import copy
import time
import threading

import numpy as np
import wx
import wx.propgrid as wxpg

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
import matplotlib.cm as cm

from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
from scipy import stats

from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim
from colormap_editor import ColormapEditorDialog
import sys_functions as sysfun
import theme
from theme import PillButton


class _SilentSimError(Exception):
    """Raised by the grid-search's own errorFunc instead of popping a modal
    dialog from a background thread. Caught per grid point / per dataset so
    one bad combination (e.g. "no resonances") doesn't abort the scan."""
    pass


def _silent_error_func(message):
    raise _SilentSimError(message)


# --------------------------------------------------------------------------- #
# Sys-dict flattening helpers
# --------------------------------------------------------------------------- #
def flatten_sys_params(sys_dict):
    """
    Flatten a Sys parameter dict (as stored in PropGridPanel.parameters,
    i.e. {'spin(1)': {...}, 'nuc(1)': {...}, ...}) into scannable scalar
    leaves. Only plain float/int values and elements of 1-D np.ndarray
    values are included; bools, strings and [value, choices] entries
    (e.g. 'Nucs', 'useFor') are skipped since a numeric grid scan over
    them isn't meaningful. The 'functions' entry (Function Modifiers state)
    is a dict too but contains no float/int/ndarray leaves at this level,
    so it's naturally skipped here -- see flatten_func_vars() for scanning
    the custom parameters (r, a, ...) defined there instead.
    """
    out = []
    for cat, sub in sys_dict.items():
        if not isinstance(sub, dict):
            continue
        for key, val in sub.items():
            if isinstance(val, bool):
                continue
            if isinstance(val, (int, float, np.floating, np.integer)):
                out.append({'kind': 'sysleaf', 'cat': cat, 'key': key, 'idx': None, 'value': float(val)})
            elif isinstance(val, np.ndarray) and val.ndim == 1:
                for ii in range(val.size):
                    out.append({'kind': 'sysleaf', 'cat': cat, 'key': key, 'idx': ii, 'value': float(val[ii])})
    return out


# Qualified label used by both this file and classFunctionModPG.py to name a
# property leaf inside a function's expression, e.g. "nuc(1).A(1)".
def qualified_label(p):
    return f"{p['cat']}.{param_leaf_label(p)}"


def flatten_func_vars(variable_values):
    """
    Flatten a {name: value} dict of Function-Modifier custom parameters
    (e.g. {'r': 3.0, 'a': 1.5}) into scannable leaves, grouped under a
    pseudo-category so they show up as their own section in the grid-search
    parameter picker.
    """
    out = []
    for name, val in variable_values.items():
        out.append({'kind': 'funcvar', 'cat': 'Functions', 'key': name, 'idx': None, 'value': float(val)})
    return out


def param_name(p):
    kind = p.get('kind', 'sysleaf')
    return f"{kind}:{p['cat']}:{p['key']}:{-1 if p['idx'] is None else p['idx']}"


def param_leaf_label(p):
    return p['key'] if p['idx'] is None else f"{p['key']}({p['idx'] + 1})"


def param_full_label(p):
    return f"{p['cat']} {param_leaf_label(p)}"


def label_to_syspath_map(sys_dict):
    """qualified 'cat.leaf' label (e.g. 'nuc(1).A(1)') -> sysleaf param
    descriptor, for patching Function-Modifier results back into a plain
    sys_dict via apply_param_value()."""
    return {qualified_label(p): p for p in flatten_sys_params(sys_dict)}


def apply_label_value(sys_dict, syspath_map, label, value):
    p = syspath_map.get(label)
    if p is None:
        return False
    apply_param_value(sys_dict, p, value)
    return True


def apply_param_value(sys_dict, param, value):
    if param.get('kind', 'sysleaf') != 'sysleaf':
        raise ValueError(f"apply_param_value only applies to 'sysleaf' params, got {param}")
    sub = sys_dict[param['cat']]
    if param['idx'] is None:
        sub[param['key']] = float(value)
    else:
        arr = sub[param['key']]
        if not isinstance(arr, np.ndarray):
            arr = np.array(arr, dtype=float)
            sub[param['key']] = arr
        arr[param['idx']] = float(value)


# --------------------------------------------------------------------------- #
# Fitness (pooled, amplitude-scaled RMSD)
# --------------------------------------------------------------------------- #
def _dataset_sse(dd, simdata, simX, simY):
    """Amplitude-scale simdata to dd's experimental data (minimizing SSE
    over the [fmin, fmax] region) and return (sum-of-squared-residuals,
    n_points) for this one dataset."""
    if dd['isfft']:
        raw = np.asarray(dd['data'])
        expdata = np.real(raw[:, :, 0]) if raw.ndim > 2 else np.real(raw)
        expX = np.asarray(dd['ax']['x'])
        expY = np.asarray(dd['ax']['y'])
    else:
        if not dd['fftactual'] or dd['fftdata'] is None:
            return 0.0, 0
        expdata = dd['fftdata']
        expX = np.asarray(dd['fftax']['x'])
        expY = np.asarray(dd['fftax']['y'])

    fmin, fmax = dd['fmin'], dd['fmax']
    xsel = (expX >= fmin) & (expX <= fmax)
    ysel = (expY >= fmin) & (expY <= fmax)
    if not xsel.any() or not ysel.any():
        return 0.0, 0

    # Array axis 0 tracks the 'x' coordinate array and axis 1 tracks 'y'
    # throughout this codebase (see MainFrame.update_FFT, which builds
    # fftX from fftData.shape[0] and fftY from fftData.shape[1]) -- index
    # accordingly rather than assuming the usual imshow row=y convention.
    exp_roi = expdata[np.ix_(xsel, ysel)]
    sub_x = expX[xsel]
    sub_y = expY[ysel]

    interp = RegularGridInterpolator((simX, simY), simdata,
                                      bounds_error=False, fill_value=0.0)
    gx, gy = np.meshgrid(sub_x, sub_y, indexing='ij')
    sim_on_exp = interp(np.stack([gx.ravel(), gy.ravel()], axis=-1)).reshape(gx.shape)

    denom = np.sum(sim_on_exp ** 2)
    c = max(np.sum(exp_roi * sim_on_exp) / denom, 0.0) if denom > 1e-30 else 0.0
    resid = exp_roi - c * sim_on_exp
    return float(np.sum(resid ** 2)), int(resid.size)


def _log_gridpoint_error(error_log, where, exc):
    """Print each distinct failure once (grid scans can hit the same
    condition, e.g. 'no resonances', at hundreds of points; printing every
    occurrence would just flood the console) and keep a count for the
    end-of-scan summary shown to the user."""
    msg = f"{where}: {type(exc).__name__}: {exc}"
    if error_log is None:
        return
    if msg not in error_log:
        print(f"[GridSearch] {msg}")
        error_log[msg] = 0
    error_log[msg] += 1


def run_gridpoint(hs, Sys, exp_template, opt_template, active_data, error_log=None):
    """Run one HYSCORE simulation per active dataset for the given Sys and
    return (pooled RMSD, pooled N) over all datasets (NaN, 0 if none
    succeeded). N (the total number of pooled ROI residual points) is
    returned alongside the RMSD so confidence-interval estimation can later
    convert RMSD back to a sum-of-squares without re-running anything.
    Failures are not fatal to the scan (a single bad grid point/dataset
    combination, e.g. "no resonances", shouldn't abort hundreds of other
    points) but are logged via _log_gridpoint_error so they are never
    silently invisible."""
    hs.Sys = Sys
    try:
        hs.preCompute()
    except Exception as e:
        _log_gridpoint_error(error_log, "preCompute", e)
        return np.nan

    total_sse = 0.0
    total_n = 0
    for dd in active_data:
        Exp = expPar()
        Exp.setDict(copy.deepcopy(exp_template))
        if Exp.tau < 0:
            Exp.tau = dd['tau']
        if Exp.Field < 0:
            Exp.Field = dd['field']
        if Exp.mwFreq < 0:
            Exp.mwFreq = dd['freq']
        if Exp.MaxFreq < 0:
            Exp.MaxFreq = dd['fmax']
        if Exp.nPoints < 0:
            Exp.nPoints = (len(dd['fftax']['x']) if not dd['isfft']
                            else len(dd['ax']['x']))

        Opt = optHYSCORE()
        Opt.setFromCtrl(copy.deepcopy(opt_template))
        if dd.get('orisel') is not None:
            Opt.OriSelInp = dd['orisel']

        hs.Exp = Exp
        hs.Opt = Opt
        try:
            hs.reRun()
        except Exception as e:
            _log_gridpoint_error(error_log, f"reRun ({dd.get('title', dd.get('fname', '?'))})", e)
            continue
        if hs.Spectrum is None:
            continue

        sse, n = _dataset_sse(dd, hs.Spectrum, hs.X, hs.Y)
        if n > 0:
            total_sse += sse
            total_n += n

    if total_n == 0:
        return np.nan, 0
    return float(np.sqrt(total_sse / total_n)), total_n


def apply_functions_to_sys_dict(sys_dict, func_rows, variable_values, error_log=None):
    """
    If Function Modifiers are defined (func_rows non-empty), evaluate them
    -- using sys_functions.evaluate_functions(), the same UI-free engine
    classFunctionModPG.py uses -- and patch the results into `sys_dict` in
    place, resolving any "other property" references directly against
    `sys_dict` itself (so e.g. a formula referencing "spin(1).S" picks up
    whatever value a scanned sysleaf parameter already patched into
    sys_dict for this grid point).

    Returns True on success. Returns False (sys_dict left only partially
    patched) if any row failed to evaluate -- callers should treat that
    grid point as a failure (NaN) rather than silently simulating with
    stale/wrong values for the properties that were supposed to be driven.
    """
    if not func_rows:
        return True

    syspath_map = label_to_syspath_map(sys_dict)
    known_labels = list(syspath_map.keys())

    def lookup(label):
        p = syspath_map.get(label)
        if p is None:
            raise sysfun.UnresolvedReference(label)
        sub = sys_dict[p['cat']]
        return float(sub[p['key']]) if p['idx'] is None else float(sub[p['key']][p['idx']])

    result = sysfun.evaluate_functions(func_rows, variable_values, lookup, known_labels=known_labels)

    if result['errors']:
        for idx, msg in result['errors'].items():
            lbl, expr = func_rows[idx]
            _log_gridpoint_error(error_log, f'function "{lbl} = {expr}"', RuntimeError(msg))
        return False

    for label, value in result['label_values'].items():
        apply_label_value(sys_dict, syspath_map, label, value)
    return True


# --------------------------------------------------------------------------- #
# Interpolation of the (coarse) scanned grid, and confidence intervals
# derived from the interpolated surface.
# --------------------------------------------------------------------------- #
def _interpolate_grid(method, grid1, grid2, results, n_out=160):
    """Interpolate a (possibly NaN-containing) coarse RMSD grid onto a
    finer regular grid, for a smoother heatmap and for locating the
    minimum/confidence contours more precisely than the coarse sampling
    allows. Returns (fine1, fine2, fine_results, used_method) or None if
    there is no finite data at all.

    NaN entries (grid points where every dataset's simulation failed) are
    replaced with a value worse than the worst finite point before
    interpolating -- treating a failed point as "confirmed bad fit" rather
    than "unknown" keeps the interpolators (which don't accept NaN) well
    defined without silently ignoring the failure. Bilinear interpolation
    of a sentinel value stays safely bounded to its local grid cell, but a
    global bicubic spline does not -- a single sentinel outlier can ring
    across the *entire* surface and fabricate a spurious dip well below
    every real sampled value (verified: on an 11x6 test grid, one bad
    point turned a true minimum of 0 into an interpolated minimum of -24).
    So 'cubic' silently falls back to 'linear' whenever the coarse grid
    contains any failed point at all, not only when an axis is too short.

    method='cubic' uses a bicubic spline (RectBivariateSpline, kx=ky=3):
    smooth and best for locating the minimum/contours when every grid
    point succeeded, but needs >=4 points along each scanned axis and can
    still overshoot between very noisy (but all-finite) samples.
    """
    finite = results[np.isfinite(results)]
    if finite.size == 0:
        return None
    has_failures = finite.size < results.size
    sentinel = float(finite.max()) * 2.0 + 1e-9
    filled = np.where(np.isfinite(results), results, sentinel)

    n1, n2 = len(grid1), len(grid2)
    fine1 = np.linspace(grid1[0], grid1[-1], max(n_out, n1))
    fine2 = np.linspace(grid2[0], grid2[-1], max(n_out, n2))

    used_method = method
    if method == 'cubic' and n1 >= 4 and n2 >= 4 and not has_failures:
        spline = RectBivariateSpline(grid1, grid2, filled, kx=3, ky=3)
        fine_vals = spline(fine1, fine2)
    else:
        used_method = 'linear'
        interp = RegularGridInterpolator((grid1, grid2), filled, method='linear',
                                          bounds_error=False, fill_value=None)
        g1, g2 = np.meshgrid(fine1, fine2, indexing='ij')
        fine_vals = interp(np.stack([g1.ravel(), g2.ravel()], axis=-1)).reshape(g1.shape)
    return fine1, fine2, fine_vals, used_method


def _rmsd_threshold_factor(method, conf, N, p=2):
    """sqrt of the SSE-inflation factor that bounds the `conf` (e.g. 0.68,
    0.95) confidence region: any grid point with rmsd <= rmsd_min*factor
    lies inside it. p is the number of jointly-scanned parameters (always
    2 here).

    'fstat' (recommended): Draper & Smith's F-test region for nonlinear
    least squares when the noise variance isn't independently known --
    SSE(theta) <= SSE_min * [1 + (p/(N-p))*F_p,(N-p)(conf)]. This is what
    the "F-test (unknown noise)" UI option uses.

    'chi2': textbook Delta-chi^2 = chi2.ppf(conf, df=p), but since there is
    no independently calibrated noise level for these spectra either, the
    fit's own rmsd_min is (ab)used as sigma -- which makes chi^2_min equal
    to N by construction and gives SSE(theta) <= SSE_min*(1 + chi2/N). This
    is simpler but under-corrects for finite N, so it tends to be more
    optimistic (narrower) than the F-test region; the two converge as
    N -> infinity.
    """
    if N is None or N <= p:
        return None
    if method == 'fstat':
        Fval = stats.f.ppf(conf, p, N - p)
        factor = 1.0 + (p / (N - p)) * Fval
    else:
        chi2val = stats.chi2.ppf(conf, df=p)
        factor = 1.0 + chi2val / N
    return float(np.sqrt(factor))


def _axis_confidence_bounds(axis_vals, profile, thresh):
    """Where a 1D profile (the interpolated surface minimized over the
    *other* scanned axis, so profile[i] = min_j fine_vals[i, j] or
    min_i fine_vals[i, j]) first drops to/below `thresh` and last rises
    back above it -- i.e. the projection of the 2D confidence region onto
    this axis. The crossing points are linearly interpolated between the
    bracketing samples rather than snapped to the nearest one, so the
    reported interval responds continuously to small changes in `thresh`
    (e.g. from switching CI method) instead of only changing once the
    threshold crosses an entire fine-grid cell -- at typical N the F-test
    vs. chi-square factors can differ by 10-25% while still landing in the
    same grid cell, which made the two methods look identical otherwise.
    Returns (lo, hi) or None if the profile never reaches thresh.
    """
    below = profile <= thresh
    if not np.any(below):
        return None
    idx = np.where(below)[0]
    i0, i1 = int(idx[0]), int(idx[-1])

    if i0 == 0:
        lo = axis_vals[0]
    else:
        x0, x1 = axis_vals[i0 - 1], axis_vals[i0]
        y0, y1 = profile[i0 - 1], profile[i0]
        lo = x1 if y1 == y0 else x0 + (thresh - y0) * (x1 - x0) / (y1 - y0)

    if i1 == len(axis_vals) - 1:
        hi = axis_vals[-1]
    else:
        x0, x1 = axis_vals[i1], axis_vals[i1 + 1]
        y0, y1 = profile[i1], profile[i1 + 1]
        hi = x0 if y1 == y0 else x0 + (thresh - y0) * (x1 - x0) / (y1 - y0)

    return float(lo), float(hi)


# --------------------------------------------------------------------------- #
# Heatmap display panel
# --------------------------------------------------------------------------- #
class HeatmapPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.figure = Figure(figsize=(5, 4.5), dpi=100)
        self.ax = self.figure.add_subplot(111)

        # self.figure.patch.set_facecolor(
        #     theme.theme_colors()["ctrl_bg"]
        # )

        # self.ax = self.figure.add_subplot(111)

        # self.ax.set_facecolor(
        #     theme.theme_colors()["ctrl_bg"]
        # )







        # figure.colorbar(im, ax=self.ax) shrinks self.ax's position to make
        # room for the colorbar, and that shrink sticks around after
        # ax.clear() (clear() wipes content, not geometry) -- so redrawing
        # repeatedly (e.g. every time the interpolation/CI dropdown changes)
        # would make the axes shrink a little more each time unless we put
        # it back to this original position first.
        self._ax_home_position = self.ax.get_position().frozen()
        theme.style_figure(self.figure, [self.ax])
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(sizer)

        # Debounce the (expensive) full-figure redraw during a live window
        # resize or splitter-sash drag -- see theme.CanvasRedrawDebouncer.
        self._redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self._redraw_debouncer.on_size)

        self.results = None
        self.grid1 = None
        self.grid2 = None
        self.label1 = ''
        self.label2 = ''
        self.cmap = cm.get_cmap('viridis')
        self.im = None
        self.cbar = None

        self.on_hover_callback = None
        self.on_rightclick_callback = None
        self.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.canvas.mpl_connect('button_press_event', self.on_click)

    def set_colormap(self, cmap):
        self.cmap = cmap
        if self.im is not None:
            self.im.set_cmap(cmap)
            self.canvas.draw_idle()

    def show_results(self, results, grid1, grid2, label1, label2,
                     contour_levels=None, raw_points=None, subtitle="",
                     background_levels=None, half_width=None):
        """results/grid1/grid2 are whatever should currently be displayed
        and hit-tested (hover/right-click) -- the interpolated fine grid
        once a scan has completed, or the raw coarse grid as a fallback.
        contour_levels: optional list of (rmsd_value, label) pairs (e.g.
        the confidence-interval boundary), drawn as a bold labelled dashed
        contour line.
        background_levels: optional list of plain rmsd values (e.g. 10%
        steps from the surface's min to max), drawn as thin unlabelled
        reference contour lines underneath contour_levels.
        raw_points: optional (grid1, grid2) of the original *coarse* scan
        locations, overlaid as small dots for reference."""
        self.results = results
        self.grid1 = grid1
        self.grid2 = grid2
        self.label1 = label1
        self.label2 = label2

        # Remove the *old* colorbar before clearing/replacing the old image
        # (self.im) -- Colorbar.remove() needs its mappable's `.axes` to
        # still be valid to restore the main axes' position/subplotspec.
        # ax.clear() detaches the old image from the axes (setting its
        # .axes to None), so doing this after clear() made remove() raise
        # partway through and silently swallow the error (caught below),
        # leaving the axes position never restored -- which is exactly why
        # the plot kept shrinking a little more on every redraw (e.g. every
        # time the interpolation/CI-method dropdown was changed).
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
        self.ax.set_position(self._ax_home_position)
        self.ax.clear()
        theme.style_figure(self.figure, [self.ax])
        extent = [grid1[0], grid1[-1], grid2[0], grid2[-1]]
        self.im = self.ax.imshow(results.T, extent=extent, origin='lower',
                                  aspect='auto', cmap=self.cmap,
                                  interpolation='nearest')
        self.ax.set_xlabel(label1)
        self.ax.set_ylabel(label2)
        self.ax.set_title('RMSD fitness map' + (f'  ({subtitle})' if subtitle else ''))

        self.cbar = self.figure.colorbar(self.im, ax=self.ax)
        self.cbar.set_label('RMSD')
        ax_fg = theme.theme_colors()["ax_fg"]
        self.cbar.ax.yaxis.label.set_color(ax_fg)
        self.cbar.ax.tick_params(colors=ax_fg)

        if raw_points is not None:
            rg1, rg2 = raw_points
            gx, gy = np.meshgrid(rg1, rg2, indexing='ij')
            self.ax.plot(gx.ravel(), gy.ravel(), '.', color='white',
                         markersize=2, alpha=0.35, zorder=3)

        if background_levels:
            blevels = sorted({lvl for lvl in background_levels if np.isfinite(lvl)})
            if blevels:
                try:
                    self.ax.contour(grid1, grid2, results.T, levels=blevels,
                                    colors='white', linewidths=0.5,
                                    linestyles='-', alpha=0.35, zorder=3)
                except Exception:
                    pass

        if contour_levels:
            levels = sorted({lvl for lvl, _lbl in contour_levels if np.isfinite(lvl)})
            if levels:
                try:
                    cs = self.ax.contour(grid1, grid2, results.T, levels=levels,
                                         colors='white', linewidths=1.6,
                                         linestyles='--', alpha=0.95, zorder=4)

                    def _fmt(v, _pairs=contour_levels):
                        for lvl, lbl in _pairs:
                            if abs(v - lvl) < 1e-9 * max(1.0, abs(lvl)):
                                return lbl
                        return f"{v:.3g}"
                    self.ax.clabel(cs, fmt=_fmt, fontsize='x-small', colors='white')
                except Exception:
                    pass

        if np.any(np.isfinite(results)):
            bi, bj = np.unravel_index(np.nanargmin(results), results.shape)
            self.ax.plot(grid1[bi], grid2[bj], marker='*', markersize=16,
                         markeredgecolor='k', markerfacecolor='white', zorder=5)

            text = (f"Best: ")
            best_ci1 = (f"")
            best_ci2 = (f"")
            if half_width is not None:
                hw1 = half_width[1]
                hw2 = half_width[2]

                if hw1 is not None: 
                    best_ci1 += f"{grid1[bi]:.4g} ± {hw1:.4g} (95%)\n"
                if hw2 is not None:
                    best_ci2 += f"{grid2[bj]:4g} ± {hw2:.4g} (95%)\n"
            text += best_ci1
            text += best_ci2
            self.ax.annotate(
                f"{text}RMSD={results[bi, bj]:.4g}",xy=(grid1[bi], grid2[bj]), 
                textcoords='offset points',
                xytext=(8, 8), fontsize='small', color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))

        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _nearest_index(self, x, y):
        if self.grid1 is None or self.grid2 is None:
            return None
        i = int(np.argmin(np.abs(self.grid1 - x)))
        j = int(np.argmin(np.abs(self.grid2 - y)))
        return i, j

    def on_move(self, event):
        if event.inaxes != self.ax or self.results is None:
            return
        idx = self._nearest_index(event.xdata, event.ydata)
        if idx is None:
            return
        i, j = idx
        if self.on_hover_callback:
            self.on_hover_callback(self.grid1[i], self.grid2[j], self.results[i, j])

    def on_click(self, event):
        if event.button != 3 or event.inaxes != self.ax or self.results is None:
            return
        idx = self._nearest_index(event.xdata, event.ydata)
        if idx is None:
            return
        i, j = idx
        if self.on_rightclick_callback:
            self.on_rightclick_callback(self.grid1[i], self.grid2[j], self.results[i, j])


# --------------------------------------------------------------------------- #
# Main grid-search window
# --------------------------------------------------------------------------- #
class GridSearchFrame(wx.Frame):
    def __init__(self, mainWindow):
        super().__init__(mainWindow, title="Sys 2D Grid Search", size=(1200, 720))
        self.main = mainWindow

        self._selected = []       # ordered list of selected param names (max 2)
        self._param_by_name = {}
        self._bool_props = {}
        self.range_ctrls = {}
        self._worker_thread = None
        self._cancel_event = threading.Event()
        self._running = False
        self._closing = False
        self._last_param1 = None
        self._last_param2 = None

        # results of the most recently completed scan (coarse grid), kept
        # around so the interpolation method / CI method dropdowns can
        # redraw and recompute instantly without re-running anything.
        self._last_results = None
        self._last_grid1 = None
        self._last_grid2 = None
        self._last_N = 0

        # remembers min/max/points the user set for each scannable
        # parameter, keyed by its stable param_name() -- so toggling a
        # parameter off and back on (or leaving one of the two checked
        # parameters untouched while swapping the other) doesn't reset
        # ranges that were already configured.
        self._range_memory = {}

        # colormap, seeded from the main window's current choice
        self.cmap_name = self.main.current_cmap_name
        self.cmap_stops = list(self.main.current_cmap_stops)
        self.cmap_base = self.main.current_cmap
        self.cmap_inverted = False
        self.cmap = self.cmap_base

        self.statusbar = self.CreateStatusBar(1)
        self.statusbar.SetStatusText("Select two Sys parameters, set their ranges, then Run.")

        splitter = wx.SplitterWindow(self, style=wx.SP_3D)
        left = wx.Panel(splitter)
        self.heatmap = HeatmapPanel(splitter)
        self.heatmap.set_colormap(self.cmap)
        self.heatmap.on_hover_callback = self._on_hover
        self.heatmap.on_rightclick_callback = self._on_rightclick

        splitter.SplitVertically(left, self.heatmap, sashPosition=440)
        splitter.SetMinimumPaneSize(260)

        leftsizer = wx.BoxSizer(wx.VERTICAL)

        lbl = wx.StaticText(left, label="Check up to two Sys parameters to scan:")
        leftsizer.Add(lbl, 0, wx.ALL, 4)

        self.param_pg = wxpg.PropertyGrid(
            left, style=wxpg.PG_DEFAULT_STYLE | wxpg.PG_HIDE_MARGIN | wxpg.PG_TOOLTIPS)
        leftsizer.Add(self.param_pg, 1, wx.EXPAND | wx.ALL, 4)
        wx.CallAfter(self.param_pg.SetSplitterPosition, 100)

        self.btn_refresh = PillButton(left, "Refresh Parameters from Main Window",
                                      color_key="pill_load")
        leftsizer.Add(self.btn_refresh, 0, wx.EXPAND | wx.ALL, 4)

        leftsizer.Add(wx.StaticLine(left), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 4)

        self.range_panel = wx.Panel(left)
        self.range_sizer = wx.FlexGridSizer(cols=4, vgap=4, hgap=6)
        self.range_sizer.AddGrowableCol(1)
        self.range_sizer.AddGrowableCol(2)
        self.range_panel.SetSizer(self.range_sizer)
        leftsizer.Add(self.range_panel, 0, wx.EXPAND | wx.ALL, 4)

        cmap_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_cmap = PillButton(left, "Colormap…", color_key="accent")
        self.cmap_swatch = wx.Panel(left, size=(70, 22), style=wx.BORDER_SIMPLE)
        self.chk_invert_cmap = wx.CheckBox(left, "Invert")
        cmap_sizer.Add(self.btn_cmap, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        cmap_sizer.Add(self.cmap_swatch, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        cmap_sizer.Add(self.chk_invert_cmap, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(cmap_sizer, 0, wx.ALL, 4)

        analysis_sizer = wx.BoxSizer(wx.HORIZONTAL)
        analysis_sizer.Add(wx.StaticText(left, label="Interpolation:"),
                           0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.interp_choice = wx.Choice(left, choices=["Linear", "Cubic"])
        self.interp_choice.SetSelection(1)
        self.interp_choice.SetToolTip(
            "Linear: bilinear, always stable, blocky/faceted.\n"
            "Cubic: smooth bicubic spline, better for locating the minimum "
            "and tracing confidence contours, but can overshoot on noisy "
            "data; needs >=4 points per scanned axis and no failed grid "
            "points, and falls back to linear otherwise.")
        analysis_sizer.Add(self.interp_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        analysis_sizer.Add(wx.StaticText(left, label="CI method:"),
                           0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.ci_method_choice = wx.Choice(
            left, choices=["F-test (unknown noise)", "Chi-square (simple)"])
        self.ci_method_choice.SetSelection(0)
        self.ci_method_choice.SetToolTip(
            "F-test: standard when the noise level isn't independently "
            "known (our case) -- more conservative (wider) intervals.\n"
            "Chi-square: textbook Delta-chi^2 using the fit's own RMSD as "
            "the noise estimate -- simpler but can be optimistic/too narrow.")
        analysis_sizer.Add(self.ci_method_choice, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(analysis_sizer, 0, wx.ALL, 4)

        self.lbl_ci = wx.StaticText(left, label="")
        self.lbl_ci.SetFont(wx.Font(9, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL,
                                     wx.FONTWEIGHT_NORMAL))
        leftsizer.Add(self.lbl_ci, 0, wx.EXPAND | wx.ALL, 4)

        note = wx.StaticText(left, label=(
            "Note: every grid point re-runs a full simulation for each active\n"
            "(Show=True) dataset. Large grids x many datasets can take a while."))
        note.SetForegroundColour(wx.Colour(theme.theme_colors()["text_dim"]))
        leftsizer.Add(note, 0, wx.ALL, 4)

        self.gauge = wx.Gauge(left, range=100)
        leftsizer.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 4)

        self.lbl_progress = wx.StaticText(left, label="")
        leftsizer.Add(self.lbl_progress, 0, wx.EXPAND | wx.ALL, 4)

        self.btn_run = PillButton(left, "Run Grid Search", color_key="accent")
        leftsizer.Add(self.btn_run, 0, wx.EXPAND | wx.ALL, 6)

        left.SetSizer(leftsizer)
        self.left_panel = left

        self.btn_refresh.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.btn_cmap.Bind(wx.EVT_BUTTON, self.on_choose_cmap)
        self.chk_invert_cmap.Bind(wx.EVT_CHECKBOX, self.on_invert_cmap_toggle)
        self.interp_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        self.ci_method_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        self.btn_run.Bind(wx.EVT_BUTTON, self.on_run_or_cancel)
        self.param_pg.Bind(wxpg.EVT_PG_CHANGED, self.on_param_toggle)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        self._load_sys_snapshot()
        self._build_param_grid()
        self._paint_swatch()
        self.Centre()

        # Inherit whichever theme the main window currently has active.
        theme.apply_theme_to_window(self, theme.get_theme())
        theme.theme_propgrid(self.param_pg, caption_key="accent")

    # ------------------------------------------------------------------
    # Colormap
    # ------------------------------------------------------------------
    def _paint_swatch(self):
        w, h = 70, 22
        bmp = wx.Bitmap(w, h)
        dc = wx.MemoryDC(bmp)
        for x in range(w):
            rgba = self.cmap(x / (w - 1))
            colour = wx.Colour(int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255))
            dc.SetPen(wx.Pen(colour))
            dc.DrawLine(x, 0, x, h)
        dc.SelectObject(wx.NullBitmap)

        self.cmap_swatch.Unbind(wx.EVT_PAINT)
        self.cmap_swatch.Bind(wx.EVT_PAINT, lambda evt, b=bmp: self._on_swatch_paint(evt, b))
        self.cmap_swatch.Refresh()

    def _on_swatch_paint(self, event, bmp):
        dc = wx.PaintDC(self.cmap_swatch)
        dc.DrawBitmap(bmp, 0, 0)

    def on_choose_cmap(self, event):
        with ColormapEditorDialog(parent=self, preset_stops=self.cmap_stops,
                                   preset_name=self.cmap_name,
                                   title="Heatmap Colormap") as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                cmap = dlg.get_colormap()
                if cmap is not None:
                    self.cmap_base = cmap
                    self.cmap_name = dlg.get_colormap_name()
                    self.cmap_stops = dlg.get_stops()
                    self._apply_effective_cmap()

    def on_invert_cmap_toggle(self, event):
        self.cmap_inverted = self.chk_invert_cmap.GetValue()
        self._apply_effective_cmap()

    def _apply_effective_cmap(self):
        self.cmap = self.cmap_base.reversed() if self.cmap_inverted else self.cmap_base
        self.heatmap.set_colormap(self.cmap)
        self._paint_swatch()

    # ------------------------------------------------------------------
    # Parameter selection grid
    # ------------------------------------------------------------------
    def _load_sys_snapshot(self):
        self.sys_snapshot = copy.deepcopy(self.main.tabulated_panel.Sys_param.parameters)
        func_panel = self.main.tabulated_panel.Sys_param.func_panel
        if func_panel is not None:
            self.func_rows_snapshot = list(func_panel.get_dict()['func'])
            self.func_var_snapshot = {name: info[func_panel.V_SPIN_ID].GetValue()
                                       for name, info in func_panel.variables.items()}
        else:
            self.func_rows_snapshot = []
            self.func_var_snapshot = {}
        # Sys leaves that a Function Modifier row currently targets aren't
        # independently scannable (they're read-only in the main grid for
        # the same reason: a direct scan of them would just be silently
        # overwritten by the formula at every grid point).
        self._driven_sys_labels = {lbl for lbl, _expr in self.func_rows_snapshot}

    def _build_param_grid(self):
        self.param_pg.Clear()
        self._bool_props = {}
        flat = [p for p in flatten_sys_params(self.sys_snapshot)
                if qualified_label(p) not in self._driven_sys_labels]
        flat += flatten_func_vars(self.func_var_snapshot)
        self._param_by_name = {param_name(p): p for p in flat}
        # drop selections whose parameter no longer exists (e.g. nucleus removed)
        self._selected = [n for n in self._selected if n in self._param_by_name]

        cur_cat = None
        cat_prop = None
        for p in flat:
            if p['cat'] != cur_cat:
                cur_cat = p['cat']
                cat_prop = self.param_pg.Append(wxpg.PropertyCategory(cur_cat))
            pname = param_name(p)
            disp = f"{param_leaf_label(p)} = {p['value']:.5g}"
            prop = self.param_pg.AppendIn(
                cat_prop, wxpg.BoolProperty(disp, pname, value=(pname in self._selected)))
            prop.SetAttribute(wxpg.PG_BOOL_USE_CHECKBOX, True)
            self._bool_props[pname] = prop
        self.param_pg.ExpandAll()
        self._sync_param_lock_state()
        self._rebuild_range_controls()

    def on_refresh(self, event):
        if self._running:
            wx.MessageBox("Cannot refresh while a scan is running.", "Grid Search",
                          wx.OK | wx.ICON_WARNING)
            return
        self._load_sys_snapshot()
        self._build_param_grid()

    def on_param_toggle(self, event):
        prop = event.GetProperty()
        pname = prop.GetName()
        if pname not in self._param_by_name:
            event.Skip()
            return
        if prop.GetValue():
            if len(self._selected) >= 2:
                prop.SetValue(False)
                self.statusbar.SetStatusText(
                    "Only two parameters can be scanned at once - deselect one first.")
                return
            self._selected.append(pname)
        else:
            if pname in self._selected:
                self._selected.remove(pname)
        self._sync_param_lock_state()
        self._rebuild_range_controls()

    def _sync_param_lock_state(self):
        lock = len(self._selected) >= 2
        for pname, prop in self._bool_props.items():
            prop.Enable(True if pname in self._selected else not lock)

    def _remember_current_ranges(self):
        """Stash whatever is currently showing in the range controls into
        self._range_memory (keyed by the stable param_name()) before the
        controls get torn down -- so a parameter that stays selected keeps
        its exact min/max/points across a rebuild, and one that gets
        deselected then reselected later comes back exactly as the user
        left it instead of reverting to a freshly-computed default range."""
        for pname, ctrls in self.range_ctrls.items():
            self._range_memory[pname] = (
                ctrls['min'].GetValue(), ctrls['max'].GetValue(), ctrls['points'].GetValue())

    def _relayout_left(self):
        # range_panel's own Layout() only repositions its children; the
        # *parent* (left) also needs relaying out so it re-queries
        # range_panel's now-changed best size and grows/shrinks the space
        # allocated to it (and shifts everything below it accordingly).
        # Skipping this step is why newly-added range controls used to
        # stay invisible until the whole window was manually resized.
        self.range_panel.Layout()
        self.left_panel.Layout()
        self.Layout()

    def _rebuild_range_controls(self):
        self._remember_current_ranges()
        self.range_sizer.Clear(True)
        self.range_ctrls = {}
        if not self._selected:
            self.range_sizer.Add(
                wx.StaticText(self.range_panel, label="(select two parameters above)"))
            self._relayout_left()
            return

        for header in ("Parameter", "Min", "Max", "Points"):
            self.range_sizer.Add(wx.StaticText(self.range_panel, label=header),
                                 0, wx.ALIGN_CENTER_VERTICAL)

        for pname in self._selected:
            p = self._param_by_name[pname]
            remembered = self._range_memory.get(pname)
            if remembered is not None:
                min_str, max_str, npts = remembered
            else:
                v = p['value']
                halfrange = max(abs(v) * 0.5, 0.5)
                min_str, max_str, npts = f"{v - halfrange:.6g}", f"{v + halfrange:.6g}", 11
                self._range_memory[pname] = (min_str, max_str, npts)
            self.range_sizer.Add(
                wx.StaticText(self.range_panel, label=param_full_label(p)),
                0, wx.ALIGN_CENTER_VERTICAL)
            min_ctrl = wx.TextCtrl(self.range_panel, value=min_str)
            max_ctrl = wx.TextCtrl(self.range_panel, value=max_str)
            pts_ctrl = wx.SpinCtrl(self.range_panel, min=2, max=500, initial=npts)
            for ctrl in (min_ctrl, max_ctrl, pts_ctrl):
                theme.theme_control(ctrl)
            self.range_sizer.Add(min_ctrl, 0, wx.EXPAND)
            self.range_sizer.Add(max_ctrl, 0, wx.EXPAND)
            self.range_sizer.Add(pts_ctrl, 0, wx.EXPAND)
            self.range_ctrls[pname] = {'min': min_ctrl, 'max': max_ctrl, 'points': pts_ctrl}

        self._relayout_left()

    # ------------------------------------------------------------------
    # Run / Cancel
    # ------------------------------------------------------------------
    def on_run_or_cancel(self, event):
        if self._running:
            self._cancel_event.set()
            self.btn_run.Disable()
            self.lbl_progress.SetLabel("Cancelling...")
            return
        self._start_run()

    def _start_run(self):
        if len(self._selected) != 2:
            wx.MessageBox("Select exactly two Sys parameters to scan.", "Grid Search",
                          wx.OK | wx.ICON_WARNING)
            return

        active_data = [dd for dd in self.main.Data if dd['show']]
        if not active_data:
            wx.MessageBox("No active (Show=True) datasets are loaded.", "Grid Search",
                          wx.OK | wx.ICON_WARNING)
            return

        try:
            ranges = []
            for pname in self._selected:
                ctrls = self.range_ctrls[pname]
                vmin = float(ctrls['min'].GetValue())
                vmax = float(ctrls['max'].GetValue())
                npts = int(ctrls['points'].GetValue())
                if vmax <= vmin:
                    raise ValueError(
                        f"Max must be greater than Min for "
                        f"{param_full_label(self._param_by_name[pname])}.")
                ranges.append((vmin, vmax, npts))
        except ValueError as e:
            wx.MessageBox(str(e), "Grid Search", wx.OK | wx.ICON_ERROR)
            return

        # Make sure every active dataset's FFT is up to date before comparing.
        try:
            self.main.update_FFT()
        except Exception as e:
            wx.MessageBox(f"Could not update FFT data before scanning:\n{e}",
                          "Grid Search", wx.OK | wx.ICON_ERROR)
            return

        param1 = self._param_by_name[self._selected[0]]
        param2 = self._param_by_name[self._selected[1]]
        grid1 = np.linspace(ranges[0][0], ranges[0][1], ranges[0][2])
        grid2 = np.linspace(ranges[1][0], ranges[1][1], ranges[1][2])
        self._last_param1, self._last_param2 = param1, param2

        sys_template = copy.deepcopy(self.main.tabulated_panel.Sys_param.parameters)
        exp_template = copy.deepcopy(self.main.tabulated_panel.Exp_param.parameters)
        opt_template = copy.deepcopy(self.main.tabulated_panel.Opt_param.parameters)
        func_rows = list(self.func_rows_snapshot)          # [[label, expr], ...]
        base_variable_values = dict(self.func_var_snapshot)  # baseline for un-scanned funcvars

        self._running = True
        self._cancel_event = threading.Event()
        self.btn_run.SetLabel("Cancel")
        self.btn_run.SetColorKey("pill_warn")
        self.param_pg.Enable(False)
        self.btn_refresh.Disable()
        self.range_panel.Enable(False)
        self.btn_cmap.Disable()
        total = len(grid1) * len(grid2)
        self.gauge.SetRange(total)
        self.gauge.SetValue(0)
        self.lbl_progress.SetLabel(f"0 / {total}")
        self.statusbar.SetStatusText("Running grid search...")

        self._worker_thread = threading.Thread(
            target=self._worker_run,
            args=(param1, param2, grid1, grid2, sys_template, exp_template,
                  opt_template, active_data, self._cancel_event,
                  func_rows, base_variable_values),
            daemon=True)
        self._worker_thread.start()

    def _worker_run(self, param1, param2, grid1, grid2, sys_template,
                     exp_template, opt_template, active_data, cancel_event,
                     func_rows, base_variable_values):
        n1, n2 = len(grid1), len(grid2)
        results = np.full((n1, n2), np.nan)
        results_n = np.zeros((n1, n2), dtype=int)
        hs = HYSCOREsim(errorFunc=_silent_error_func)
        hs.verbose = False
        start = time.perf_counter()
        count = 0
        total = n1 * n2
        n_failed_points = 0
        error_log = {}   # unique error message -> occurrence count; also printed live

        for i, v1 in enumerate(grid1):
            if cancel_event.is_set():
                break
            for j, v2 in enumerate(grid2):
                if cancel_event.is_set():
                    break
                sys_dict = copy.deepcopy(sys_template)
                variable_values = dict(base_variable_values)
                for p, v in ((param1, v1), (param2, v2)):
                    if p.get('kind', 'sysleaf') == 'funcvar':
                        variable_values[p['key']] = v
                    else:
                        apply_param_value(sys_dict, p, v)
                try:
                    func_ok = apply_functions_to_sys_dict(
                        sys_dict, func_rows, variable_values, error_log=error_log)
                    if not func_ok:
                        rmsd, n_pts = np.nan, 0
                    else:
                        Sys = sysPar()
                        Sys.setFromCtrl(sys_dict)
                        rmsd, n_pts = run_gridpoint(hs, Sys, exp_template, opt_template,
                                                     active_data, error_log=error_log)
                except Exception as e:
                    _log_gridpoint_error(error_log, "Sys.setFromCtrl", e)
                    rmsd, n_pts = np.nan, 0
                if np.isnan(rmsd):
                    n_failed_points += 1
                results[i, j] = rmsd
                results_n[i, j] = n_pts
                count += 1
                elapsed = time.perf_counter() - start
                wx.CallAfter(self._on_progress, count, total, elapsed)

        wx.CallAfter(self._on_finished, results, results_n, grid1, grid2, param1, param2,
                     cancel_event.is_set(), n_failed_points, total, error_log)

    def _on_progress(self, count, total, elapsed):
        if self._closing:
            return
        self.gauge.SetValue(count)
        rate = elapsed / count if count else 0.0
        remaining = rate * (total - count)
        self.lbl_progress.SetLabel(
            f"{count} / {total}   elapsed {elapsed:.0f}s   ETA {remaining:.0f}s")

    def _on_finished(self, results, results_n, grid1, grid2, param1, param2, cancelled,
                      n_failed, total, error_log=None):
        if self._closing:
            return
        self._running = False
        self.btn_run.SetLabel("Run Grid Search")
        self.btn_run.SetColorKey("accent")
        self.btn_run.Enable()
        self.param_pg.Enable(True)
        self.btn_refresh.Enable()
        self.range_panel.Enable(True)
        self.btn_cmap.Enable()
        self._sync_param_lock_state()

        self._last_results = results
        self._last_grid1 = grid1
        self._last_grid2 = grid2
        self._last_param1 = param1
        self._last_param2 = param2
        if np.any(np.isfinite(results)):
            bi, bj = np.unravel_index(np.nanargmin(results), results.shape)
            self._last_N = int(results_n[bi, bj])
        else:
            self._last_N = 0
        self._refresh_display()

        msg = "Grid search cancelled." if cancelled else "Grid search finished."
        if n_failed:
            msg += f" {n_failed}/{total} grid points had no successful simulation (see console)."
        self.statusbar.SetStatusText(msg)

        if error_log:
            examples = "\n".join(f"  x{count}  {m}" for m, count in
                                 sorted(error_log.items(), key=lambda kv: -kv[1])[:5])
            wx.MessageBox(
                f"{n_failed}/{total} grid points had no successful simulation.\n\n"
                f"Most common reasons:\n{examples}\n\n"
                f"(full detail was printed to the console as it happened)",
                "Grid Search - some points failed", wx.OK | wx.ICON_WARNING)

    # ------------------------------------------------------------------
    # Interpolation / confidence intervals
    # ------------------------------------------------------------------
    def on_analysis_choice_changed(self, event):
        self._refresh_display()

    def _refresh_display(self):
        """(Re)draw the heatmap from the coarse results cached in
        self._last_results using the currently-selected interpolation
        method, and (re)compute the confidence-interval display using the
        currently-selected CI method. Cheap enough to call on every
        dropdown change -- no re-scanning involved."""
        if self._last_results is None:
            return
        grid1, grid2 = self._last_grid1, self._last_grid2
        results = self._last_results
        param1, param2 = self._last_param1, self._last_param2
        label1, label2 = param_full_label(param1), param_full_label(param2)

        method = 'cubic' if self.interp_choice.GetSelection() == 1 else 'linear'
        out = _interpolate_grid(method, grid1, grid2, results)
        if out is None:
            self.heatmap.show_results(results, grid1, grid2, label1, label2)
            self.lbl_ci.SetLabel(
                "No successful grid points -- cannot interpolate or estimate confidence intervals.")
            self.left_panel.Layout()
            self.Layout()
            return
        fine1, fine2, fine_vals, used_method = out

        bi, bj = np.unravel_index(np.nanargmin(fine_vals), fine_vals.shape)
        rmsd_min, best1, best2 = fine_vals[bi, bj], fine1[bi], fine2[bj]

        # Thin reference contours at 10% steps of the interpolated surface's
        # own range -- purely a visual aid for reading the landscape, drawn
        # under the (bold) confidence-interval contour.
        vmin, vmax = float(np.nanmin(fine_vals)), float(np.nanmax(fine_vals))
        background_levels = list(np.linspace(vmin, vmax, 11)[1:-1]) if vmax > vmin else []

        ci_method = 'fstat' if self.ci_method_choice.GetSelection() == 0 else 'chi2'
        N = self._last_N
        CONF = 0.95
        factor = _rmsd_threshold_factor(ci_method, CONF, N, p=2)
        contour_levels = []
        half_width = {1: None, 2: None}
        if factor is not None:
            thresh = rmsd_min * factor
            if np.nanmin(fine_vals) <= thresh:
                contour_levels.append((thresh, "95%"))
                profile1 = np.nanmin(fine_vals, axis=1)  # min over axis2, vs. axis1
                profile2 = np.nanmin(fine_vals, axis=0)  # min over axis1, vs. axis2
                bounds1 = _axis_confidence_bounds(fine1, profile1, thresh)
                bounds2 = _axis_confidence_bounds(fine2, profile2, thresh)
                if bounds1 is not None:
                    half_width[1] = (bounds1[1] - bounds1[0]) / 2.0
                if bounds2 is not None:
                    half_width[2] = (bounds2[1] - bounds2[0]) / 2.0

        def fmt_line(label, best, hw):
            if hw is None:
                return f"{label} = {best:.5g}  (95% CI: n/a)"
            return f"{label} = {best:.5g} ± {hw:.5g}  (95%)"

        lines = [fmt_line(label1, best1, half_width[1]), fmt_line(label2, best2, half_width[2])]
        if used_method != method:
            lines.append("(cubic needs >=4 points/axis and no failed grid points -- used linear instead)")
        if N <= 2:
            lines.append(f"(only N={N} pooled residual points at the best fit -- "
                         f"too few for a meaningful confidence interval)")
        lines.append("± is half the width of the joint 2-parameter 95% region's bounding "
                     "box on each axis (conservative vs. a true 1D profile interval).")
        self.lbl_ci.SetLabel("\n".join(lines))

        subtitle = f"{used_method} interp, N={N}"
        self.heatmap.show_results(fine_vals, fine1, fine2, label1, label2,
                                  contour_levels=contour_levels,
                                  background_levels=background_levels,
                                  raw_points=(grid1, grid2), subtitle=subtitle, half_width=half_width)
        self.left_panel.Layout()
        self.Layout()

    # ------------------------------------------------------------------
    # Heatmap interaction
    # ------------------------------------------------------------------
    def _on_hover(self, x, y, val):
        vtxt = "n/a" if np.isnan(val) else f"{val:.5g}"
        self.statusbar.SetStatusText(
            f"{self.heatmap.label1} = {x:.5g},  {self.heatmap.label2} = {y:.5g},  RMSD = {vtxt}")

    def _on_rightclick(self, x, y, val):
        if self._last_param1 is None or self._last_param2 is None:
            return
        menu = wx.Menu()
        item = menu.Append(wx.ID_ANY, "Load these values into main Sys grid")
        self.Bind(wx.EVT_MENU, lambda evt, xx=x, yy=y: self._load_into_main(xx, yy), item)
        pos = self.heatmap.canvas.ScreenToClient(wx.GetMousePosition())
        self.heatmap.canvas.PopupMenu(menu, pos)
        menu.Destroy()

    def _load_into_main(self, v1, v2):
        if self._running:
            return
        main_params = self.main.tabulated_panel.Sys_param.parameters
        try:
            for p, v in ((self._last_param1, v1), (self._last_param2, v2)):
                if p.get('kind', 'sysleaf') == 'funcvar':
                    # Custom Function-Modifier parameters (r, a, ...) don't
                    # live as a normal Sys leaf -- update the stashed
                    # 'functions'.'var' entry instead; SetFromParClean()
                    # below rebuilds the func panel from it, which in turn
                    # recomputes every driven Sys leaf from the new value.
                    fdict = main_params.get('functions')
                    if fdict is None or 'var' not in fdict or p['key'] not in fdict['var']:
                        raise KeyError(p['key'])
                    fdict['var'][p['key']] = v
                else:
                    apply_param_value(main_params, p, v)
        except KeyError:
            wx.MessageBox(
                "The scanned parameter no longer exists in the main Sys grid "
                "(it may have been removed/renamed). Could not load values.",
                "Grid Search", wx.OK | wx.ICON_WARNING)
            return
        self.main.tabulated_panel.Sys_param.SetFromParClean(main_params)
        self.main.Sys.setFromCtrl(main_params)
        self.main.runSim()
        msg = (f"Loaded {param_full_label(self._last_param1)}={v1:.5g}, "
               f"{param_full_label(self._last_param2)}={v2:.5g} into main window "
               f"and re-ran simulation.")
        # The main Sys grid just changed under us -- refresh this window's
        # own parameter list so the picker's displayed values (and which
        # leaves are Function-Modifier-driven) reflect the new baseline,
        # same as clicking "Refresh Parameters from Main Window".
        self.on_refresh(None)
        self.statusbar.SetStatusText(msg)

    # ------------------------------------------------------------------
    def on_close(self, event):
        self._closing = True
        if self._running:
            self._cancel_event.set()
        event.Skip()
