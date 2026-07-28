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

from scipy.interpolate import RegularGridInterpolator

from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim
from colormap_editor import ColormapEditorDialog


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
    them isn't meaningful.
    """
    out = []
    for cat, sub in sys_dict.items():
        if not isinstance(sub, dict):
            continue
        for key, val in sub.items():
            if isinstance(val, bool):
                continue
            if isinstance(val, (int, float, np.floating, np.integer)):
                out.append({'cat': cat, 'key': key, 'idx': None, 'value': float(val)})
            elif isinstance(val, np.ndarray) and val.ndim == 1:
                for ii in range(val.size):
                    out.append({'cat': cat, 'key': key, 'idx': ii, 'value': float(val[ii])})
    return out


def param_name(p):
    return f"{p['cat']}:{p['key']}:{-1 if p['idx'] is None else p['idx']}"


def param_leaf_label(p):
    return p['key'] if p['idx'] is None else f"{p['key']}({p['idx'] + 1})"


def param_full_label(p):
    return f"{p['cat']} {param_leaf_label(p)}"


def apply_param_value(sys_dict, param, value):
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
    return the pooled RMSD over all datasets (NaN if none succeeded).
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
        return np.nan
    return float(np.sqrt(total_sse / total_n))


# --------------------------------------------------------------------------- #
# Heatmap display panel
# --------------------------------------------------------------------------- #
class HeatmapPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.figure = Figure(figsize=(5, 4.5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(sizer)

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

    def show_results(self, results, grid1, grid2, label1, label2):
        self.results = results
        self.grid1 = grid1
        self.grid2 = grid2
        self.label1 = label1
        self.label2 = label2

        self.ax.clear()
        extent = [grid1[0], grid1[-1], grid2[0], grid2[-1]]
        self.im = self.ax.imshow(results.T, extent=extent, origin='lower',
                                  aspect='auto', cmap=self.cmap,
                                  interpolation='nearest')
        self.ax.set_xlabel(label1)
        self.ax.set_ylabel(label2)
        self.ax.set_title('RMSD fitness map')

        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
        self.cbar = self.figure.colorbar(self.im, ax=self.ax)
        self.cbar.set_label('RMSD')

        if np.any(np.isfinite(results)):
            bi, bj = np.unravel_index(np.nanargmin(results), results.shape)
            self.ax.plot(grid1[bi], grid2[bj], marker='*', markersize=16,
                         markeredgecolor='k', markerfacecolor='white', zorder=5)
            self.ax.annotate(
                f"best: {grid1[bi]:.4g}, {grid2[bj]:.4g}\nRMSD={results[bi, bj]:.4g}",
                (grid1[bi], grid2[bj]), textcoords='offset points',
                xytext=(8, 8), fontsize='small', color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))

        self.figure.tight_layout()
        self.canvas.draw()

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

        # colormap, seeded from the main window's current choice
        self.cmap_name = self.main.current_cmap_name
        self.cmap_stops = list(self.main.current_cmap_stops)
        self.cmap = self.main.current_cmap

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

        self.btn_refresh = wx.Button(left, label="Refresh Parameters from Main Window")
        leftsizer.Add(self.btn_refresh, 0, wx.EXPAND | wx.ALL, 4)

        leftsizer.Add(wx.StaticLine(left), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 4)

        self.range_panel = wx.Panel(left)
        self.range_sizer = wx.FlexGridSizer(cols=4, vgap=4, hgap=6)
        self.range_sizer.AddGrowableCol(1)
        self.range_sizer.AddGrowableCol(2)
        self.range_panel.SetSizer(self.range_sizer)
        leftsizer.Add(self.range_panel, 0, wx.EXPAND | wx.ALL, 4)

        cmap_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_cmap = wx.Button(left, label="Colormap…")
        self.cmap_swatch = wx.Panel(left, size=(70, 22), style=wx.BORDER_SIMPLE)
        cmap_sizer.Add(self.btn_cmap, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        cmap_sizer.Add(self.cmap_swatch, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(cmap_sizer, 0, wx.ALL, 4)

        note = wx.StaticText(left, label=(
            "Note: every grid point re-runs a full simulation for each active\n"
            "(Show=True) dataset. Large grids x many datasets can take a while."))
        note.SetForegroundColour(wx.Colour(120, 120, 120))
        leftsizer.Add(note, 0, wx.ALL, 4)

        self.gauge = wx.Gauge(left, range=100)
        leftsizer.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 4)

        self.lbl_progress = wx.StaticText(left, label="")
        leftsizer.Add(self.lbl_progress, 0, wx.EXPAND | wx.ALL, 4)

        self.btn_run = wx.Button(left, label="Run Grid Search")
        leftsizer.Add(self.btn_run, 0, wx.EXPAND | wx.ALL, 6)

        left.SetSizer(leftsizer)

        self.btn_refresh.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.btn_cmap.Bind(wx.EVT_BUTTON, self.on_choose_cmap)
        self.btn_run.Bind(wx.EVT_BUTTON, self.on_run_or_cancel)
        self.param_pg.Bind(wxpg.EVT_PG_CHANGED, self.on_param_toggle)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        self._load_sys_snapshot()
        self._build_param_grid()
        self._paint_swatch()
        self.Centre()

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
                    self.cmap = cmap
                    self.cmap_name = dlg.get_colormap_name()
                    self.cmap_stops = dlg.get_stops()
                    self.heatmap.set_colormap(self.cmap)
                    self._paint_swatch()

    # ------------------------------------------------------------------
    # Parameter selection grid
    # ------------------------------------------------------------------
    def _load_sys_snapshot(self):
        self.sys_snapshot = copy.deepcopy(self.main.tabulated_panel.Sys_param.parameters)

    def _build_param_grid(self):
        self.param_pg.Clear()
        self._bool_props = {}
        flat = flatten_sys_params(self.sys_snapshot)
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

    def _rebuild_range_controls(self):
        self.range_sizer.Clear(True)
        self.range_ctrls = {}
        if not self._selected:
            self.range_sizer.Add(
                wx.StaticText(self.range_panel, label="(select two parameters above)"))
            self.range_panel.Layout()
            self.Layout()
            return

        for header in ("Parameter", "Min", "Max", "Points"):
            self.range_sizer.Add(wx.StaticText(self.range_panel, label=header),
                                 0, wx.ALIGN_CENTER_VERTICAL)

        for pname in self._selected:
            p = self._param_by_name[pname]
            v = p['value']
            halfrange = max(abs(v) * 0.5, 0.5)
            self.range_sizer.Add(
                wx.StaticText(self.range_panel, label=param_full_label(p)),
                0, wx.ALIGN_CENTER_VERTICAL)
            min_ctrl = wx.TextCtrl(self.range_panel, value=f"{v - halfrange:.6g}")
            max_ctrl = wx.TextCtrl(self.range_panel, value=f"{v + halfrange:.6g}")
            pts_ctrl = wx.SpinCtrl(self.range_panel, min=2, max=500, initial=11)
            self.range_sizer.Add(min_ctrl, 0, wx.EXPAND)
            self.range_sizer.Add(max_ctrl, 0, wx.EXPAND)
            self.range_sizer.Add(pts_ctrl, 0, wx.EXPAND)
            self.range_ctrls[pname] = {'min': min_ctrl, 'max': max_ctrl, 'points': pts_ctrl}

        self.range_panel.Layout()
        self.Layout()

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

        self._running = True
        self._cancel_event = threading.Event()
        self.btn_run.SetLabel("Cancel")
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
                  opt_template, active_data, self._cancel_event),
            daemon=True)
        self._worker_thread.start()

    def _worker_run(self, param1, param2, grid1, grid2, sys_template,
                     exp_template, opt_template, active_data, cancel_event):
        n1, n2 = len(grid1), len(grid2)
        results = np.full((n1, n2), np.nan)
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
                apply_param_value(sys_dict, param1, v1)
                apply_param_value(sys_dict, param2, v2)
                try:
                    Sys = sysPar()
                    Sys.setFromCtrl(sys_dict)
                    rmsd = run_gridpoint(hs, Sys, exp_template, opt_template,
                                          active_data, error_log=error_log)
                except Exception as e:
                    _log_gridpoint_error(error_log, "Sys.setFromCtrl", e)
                    rmsd = np.nan
                if np.isnan(rmsd):
                    n_failed_points += 1
                results[i, j] = rmsd
                count += 1
                elapsed = time.perf_counter() - start
                wx.CallAfter(self._on_progress, count, total, elapsed)

        wx.CallAfter(self._on_finished, results, grid1, grid2, param1, param2,
                     cancel_event.is_set(), n_failed_points, total, error_log)

    def _on_progress(self, count, total, elapsed):
        if self._closing:
            return
        self.gauge.SetValue(count)
        rate = elapsed / count if count else 0.0
        remaining = rate * (total - count)
        self.lbl_progress.SetLabel(
            f"{count} / {total}   elapsed {elapsed:.0f}s   ETA {remaining:.0f}s")

    def _on_finished(self, results, grid1, grid2, param1, param2, cancelled,
                      n_failed, total, error_log=None):
        if self._closing:
            return
        self._running = False
        self.btn_run.SetLabel("Run Grid Search")
        self.btn_run.Enable()
        self.param_pg.Enable(True)
        self.btn_refresh.Enable()
        self.range_panel.Enable(True)
        self.btn_cmap.Enable()
        self._sync_param_lock_state()

        label1 = param_full_label(param1)
        label2 = param_full_label(param2)
        self.heatmap.show_results(results, grid1, grid2, label1, label2)

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
            apply_param_value(main_params, self._last_param1, v1)
            apply_param_value(main_params, self._last_param2, v2)
        except KeyError:
            wx.MessageBox(
                "The scanned parameter no longer exists in the main Sys grid "
                "(it may have been removed/renamed). Could not load values.",
                "Grid Search", wx.OK | wx.ICON_WARNING)
            return
        self.main.tabulated_panel.Sys_param.SetFromParClean(main_params)
        self.main.Sys.setFromCtrl(main_params)
        self.main.runSim()
        self.statusbar.SetStatusText(
            f"Loaded {param_full_label(self._last_param1)}={v1:.5g}, "
            f"{param_full_label(self._last_param2)}={v2:.5g} into main window "
            f"and re-ran simulation.")

    # ------------------------------------------------------------------
    def on_close(self, event):
        self._closing = True
        if self._running:
            self._cancel_event.set()
        event.Skip()
