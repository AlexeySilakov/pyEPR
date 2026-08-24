"""
2D grid-search / fitness-mapping window for pyHYSCORE.

Scans two scalar Sys parameters over a user-defined grid. At every grid
point, a full HYSCORE simulation is re-run for every "active" dataset
(active == the dataset's "Show" checkbox in the main window's Data tab)
and compared against that dataset's experimental spectrum, restricted to
the region the main window's 2D plots outline with a dashed line: the
45-degree tilted box inscribed in [fmin, fmax] x [fmin, fmax], whose
sides run parallel and perpendicular to f2 = f1 and whose width is set by
that dataset's AntiDiagSpan (see optHYSCORE.diagonalBox). That box is also exactly
what the main window's diagonal skyline panel projects, and it is the
same region regardless of the "Quadrants" radio button. Comparison is
always done in the frequency domain (FFT'd data for time-domain files,
raw data for files already loaded as FFT spectra), independent of the
main window's own "FFT(y)" display toggle.

Fitness = pooled RMSD:
  For each active dataset k, the simulated spectrum is amplitude-scaled
  by a non-negative factor c_k that minimizes ||exp_k - c_k*sim_k||^2
  over its ROI (closed-form least squares). The residuals of all active
  datasets are then pooled together:

      RMSD = sqrt( sum_k sum_roi (exp_k - c_k*sim_k)^2  /  sum_k N_k )

  i.e. each dataset gets its own best-fit amplitude, but the reported
  number is the RMSD of the combined residual pool, not an average of
  per-dataset RMSDs.

Normalized RMSD (optional, "Normalized RMSD" checkbox) rescales that into
an absolute 0-1 score by dividing each dataset's SSE by its own *null*
SSE -- sum_roi(exp_k^2), the residual left by no simulation at all:

    RMSD_norm = sqrt( sum_k N_k*(SSE_k/null_k) / sum_k N_k )

Because the amplitude fit clamps c_k >= 0, c_k = 0 reproduces the null
SSE exactly and the fit can never do worse, so each fraction -- and hence
RMSD_norm -- is bounded in (0, 1]: 0 is a perfect fit, 1 means the
simulation explains nothing. Equivalently RMSD_norm = sqrt(1 - r^2) with
r the uncentred correlation between experiment and simulation, so 0.3
means 91% of the experimental power is accounted for.
Normalizing per dataset *before* averaging (rather than pooling both sums)
stops a bright dataset from dominating a multi-dataset fit, but makes this
a genuinely different statistic: with more than one dataset its minimum
can sit somewhere other than the raw pooled RMSD's. With one dataset it is
just the raw RMSD over a constant. Note the normalization is uncentred, so
a large baseline in the experimental spectrum inflates null_k and flatters
the score.

Optionally ("Keep simulations in memory") every grid point's simulation
is retained -- interpolated onto each dataset's ROI only, which is what
makes this affordable at all -- so that double-clicking any point on the
finished fitness map can show that point's per-dataset residual maps
without re-simulating. The predicted cost is
8 bytes x n1 x n2 x sum_k(ROI points of dataset k); the user is warned and
offered a no-storage run before the scan starts if that exceeds 0.5 GB.

The "Display:" pulldown chooses what the finished map shows:

  RMSD                         the raw pooled fitness, above.
  Normalized RMSD              the 0-1 null-model score, above.
  Probability (credible level) each cell labelled with the smallest
                               credible region containing it, so the 0.95
                               contour *is* the 95% region.
  Probability (density)        the posterior density itself; same contours,
                               but a sharp posterior leaves most cells at
                               effectively zero.

The posterior assumes Gaussian residuals with a separate unknown noise
level per dataset and integrates those levels out, giving
p(theta) ~ prod_k SSE_k(theta)**(-N_k/2) -- so what is minimized is
sum_k (N_k/2)*log SSE_k. Marginalizing per dataset makes this immune both
to dataset brightness and to the RMSD normalization choice above (they
only shift it by a constant), which is why the probability map does not
move when the fitness map does. Confidence intervals in these modes come
from the posterior marginals rather than from the F-test.

Because N sits in an exponent there, the "Oversampling" divisor matters:
adjacent spectrum points are correlated by apodization and especially by
zero-filling, so N is divided by the number of points per independent
measurement. It defaults to 4**zfill from the main window's zero-fill
setting (2**zfill per axis, squared for the 2-D point count). Left at 1
the credible regions are optimistically narrow. The panel warns when the
posterior spills off the edge of the scanned box (widen the ranges) or
collapses inside a single cell (narrow them).

Once a scan finishes, the coarse RMSD grid is interpolated (linear or
cubic, selectable live without re-scanning -- see _interpolate_grid) onto
a finer grid for display, and used to estimate a confidence interval for
each scanned parameter (F-test, naive chi-square, or measured-noise
chi-square, selectable -- see _rmsd_threshold_factor) at both 68% and 95%,
reported as the bounding box of the joint 2-parameter confidence region
projected onto each axis.

The measured-noise option is the only one of the three that uses a noise
level from outside the fit. pyHYSCORE measures one per dataset at the end of
every FFT, as the mean of the second half of the +- quadrant's extreme-y
x-trace -- the far corner of the map in both frequencies, where a HYSCORE
spectrum holds no signal but the same noise as everywhere else. Since the
spectrum is a magnitude it never falls to zero, and that mean is within
about 11% of the RMS a chi-square denominator wants. It rides along on the
dataset dict, so it survives a session save/load and needs no rescan
(MainFrame.noise_from_spectrum). Here it is pooled across datasets by point
count (pooled_sigma), in the units of whichever statistic is on display, and
the region becomes SSE <= SSE_min + divisor*Delta-chi^2*sigma^2: a width set
by how noisy the data is rather than by how good the fit happens to be. The
same oversampling divisor as the posterior applies, for the same reason.
Whenever a noise level is available the panel also reports the reduced
chi^2 of the best fit, which says plainly whether the residual has reached
the noise or stalled well above it -- and warns when it has, since a region
scaled to the noise alone is then optimistically narrow.
"""
import copy
import time
import threading

import numpy as np
import wx
import wx.propgrid as wxpg

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
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
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
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
        self.value_label = 'RMSD'
        self.cmap = cm.get_cmap('viridis')
        self.im = None
        self.cbar = None

        self.on_hover_callback = None
        self.on_rightclick_callback = None
        self.on_doubleclick_callback = None
        self.canvas.mpl_connect('motion_notify_event', self.on_move)
        self.canvas.mpl_connect('button_press_event', self.on_click)

    def set_colormap(self, cmap):
        self.cmap = cmap
        if self.im is not None:
            # A ContourSet (contour display mode) is a mappable too, but
            # recolouring it in place does not restyle the already-drawn
            # line collections -- the frame follows this with a full
            # redraw whenever a scan's results exist.
            try:
                self.im.set_cmap(cmap)
            except Exception:
                pass
            self.canvas.draw_idle()

    def show_results(self, results, grid1, grid2, label1, label2,
                     contour_levels=None, raw_points=None, subtitle="",
                     background_levels=None, half_width=None,
                     as_contours=False, rmsd_contour_width=1.2,
                     value_label="RMSD", best_index=None, map_title=None,
                     bg_contour_color='white', bg_contour_width=0.5,
                     grid_lines=True, mark_residual=None, mark_loaded=None):
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
        self.value_label = value_label

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
        ax_fg = theme.theme_colors()["ax_fg"]
        blevels = sorted({lvl for lvl in (background_levels or [])
                          if np.isfinite(lvl)})

        # Contour display mode replaces the filled image with colormap-
        # coloured lines. Falling back to the image when there are no
        # levels to draw (a perfectly flat surface) avoids an empty plot.
        drew_contours = False
        if as_contours and blevels:
            try:
                self.im = self.ax.contour(grid1, grid2, results.T, levels=blevels,
                                           cmap=self.cmap,
                                           linewidths=rmsd_contour_width)
                self.ax.set_xlim(grid1[0], grid1[-1])
                self.ax.set_ylim(grid2[0], grid2[-1])
                drew_contours = True
            except Exception:
                drew_contours = False
        if not drew_contours:
            self.im = self.ax.imshow(results.T, extent=extent, origin='lower',
                                      aspect='auto', cmap=self.cmap,
                                      interpolation='nearest')

        # Overlays are drawn on top of the colormap in image mode (white
        # reads against any colormap) but on top of the plain axes
        # background in contour mode, where white would be invisible.
        overlay_color = ax_fg if drew_contours else 'white'

        self.ax.set_xlabel(label1)
        self.ax.set_ylabel(label2)
        self.ax.set_title((map_title or f'{value_label} fitness map')
                          + (f'  ({subtitle})' if subtitle else ''))

        self.cbar = self.figure.colorbar(self.im, ax=self.ax)
        self.cbar.set_label(value_label)
        self.cbar.ax.yaxis.label.set_color(ax_fg)
        self.cbar.ax.tick_params(colors=ax_fg)

        if raw_points is not None:
            rg1, rg2 = raw_points
            # Rules through every scanned value: their intersections are
            # exactly the grid points that were actually simulated, so it is
            # obvious how coarse the scan underneath the smooth interpolated
            # surface really is.
            if grid_lines:
                self.ax.vlines(rg1, grid2[0], grid2[-1], colors=overlay_color,
                               linewidths=0.5, alpha=0.4, linestyles=':', zorder=2)
                self.ax.hlines(rg2, grid1[0], grid1[-1], colors=overlay_color,
                               linewidths=0.5, alpha=0.4, linestyles=':', zorder=2)
            gx, gy = np.meshgrid(rg1, rg2, indexing='ij')
            self.ax.plot(gx.ravel(), gy.ravel(), '.', color=overlay_color,
                         markersize=2, alpha=0.35, zorder=3)

        # In contour mode these thin reference lines *are* the display, so
        # drawing them again underneath would only double up.
        if blevels and not drew_contours:
            try:
                self.ax.contour(grid1, grid2, results.T, levels=blevels,
                                colors=bg_contour_color,
                                linewidths=bg_contour_width,
                                linestyles='-', alpha=0.35, zorder=3)
            except Exception:
                pass

        if contour_levels:
            levels = sorted({lvl for lvl, _lbl in contour_levels if np.isfinite(lvl)})
            if levels:
                try:
                    cs = self.ax.contour(grid1, grid2, results.T, levels=levels,
                                         colors=overlay_color, linewidths=1.6,
                                         linestyles='--', alpha=0.95, zorder=4)

                    def _fmt(v, _pairs=contour_levels):
                        for lvl, lbl in _pairs:
                            if abs(v - lvl) < 1e-9 * max(1.0, abs(lvl)):
                                return lbl
                        return f"{v:.3g}"
                    self.ax.clabel(cs, fmt=_fmt, fontsize='x-small',
                                   colors=overlay_color)
                except Exception:
                    pass

        if np.any(np.isfinite(results)):
            # best_index is supplied when "best" is not the minimum of the
            # displayed surface -- the posterior-density map peaks at the
            # best fit rather than dipping to it.
            if best_index is not None:
                bi, bj = best_index
            else:
                bi, bj = np.unravel_index(np.nanargmin(results), results.shape)
            self.ax.plot(grid1[bi], grid2[bj], marker='*', markersize=16,
                         markeredgecolor='k', markerfacecolor='white',
                         linestyle='none', label='best fit', zorder=5)

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
                f"{text}{value_label}={results[bi, bj]:.4g}",xy=(grid1[bi], grid2[bj]),
                textcoords='offset points',
                xytext=(8, 8), fontsize='small', color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))

        # Two more places worth finding again: the grid point whose
        # residuals are on screen, and where the main window's Sys currently
        # sits. Deliberately different shapes as well as different colours,
        # so they stay distinguishable against any colormap and for anyone
        # who reads the two colours as similar.
        marked = False
        if mark_residual is not None:
            self.ax.plot(mark_residual[0], mark_residual[1], marker='o',
                         markersize=11, markerfacecolor='none',
                         markeredgecolor='#00e5ff', markeredgewidth=2.0,
                         linestyle='none', label='residuals shown', zorder=6)
            marked = True
        if mark_loaded is not None:
            self.ax.plot(mark_loaded[0], mark_loaded[1], marker='D',
                         markersize=10, markerfacecolor='none',
                         markeredgecolor='#ff9800', markeredgewidth=2.0,
                         linestyle='none', label='current Sys', zorder=6)
            marked = True
        if marked:
            leg = self.ax.legend(loc='upper left', fontsize='xx-small',
                                 framealpha=0.65, handlelength=1.2,
                                 borderpad=0.4, labelspacing=0.3)
            leg.set_zorder(7)
            for text in leg.get_texts():
                text.set_color('black')

        # Pin the view to the scanned range. Everything overlaid above --
        # markers, the coarse-point dots, the grid rules -- goes through
        # ax.plot and so feeds the autoscaler; without this a single stray
        # point anywhere off-grid silently rescales the whole map.
        self.ax.set_xlim(grid1[0], grid1[-1])
        self.ax.set_ylim(grid2[0], grid2[-1])

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
        if event.inaxes != self.ax or self.results is None:
            return
        # Double-click reports the raw data coordinates rather than a
        # snapped grid1/grid2 pair: what is displayed here is the *fine*
        # interpolated grid, but stored simulations only exist at the
        # original coarse scan points, so the frame does its own snapping
        # against the coarse axes it kept.
        if event.button == 1 and event.dblclick:
            if self.on_doubleclick_callback:
                self.on_doubleclick_callback(event.xdata, event.ydata)
            return
        if event.button != 3:
            return
        idx = self._nearest_index(event.xdata, event.ydata)
        if idx is None:
            return
        i, j = idx
        if self.on_rightclick_callback:
            self.on_rightclick_callback(self.grid1[i], self.grid2[j], self.results[i, j])


# --------------------------------------------------------------------------- #
# Residual maps (one per active dataset, for a single grid point)
# --------------------------------------------------------------------------- #
class ContourSettingsDialog(wx.Dialog):
    """Appearance of a set of contour lines: on/off, colour, how many
    levels, how thick. Shared by the simulation contours over the residual
    maps (right-click the row) and the reference contours over the fitness
    map (the "Contours..." button), so the two behave identically.

    Seeded from the caller's current settings so re-opening it always shows
    what is actually on screen."""

    def __init__(self, parent, title, enabled, color, levels, width,
                 enable_label="Draw contours"):
        super().__init__(parent, title=title, style=wx.DEFAULT_DIALOG_STYLE)

        self.chk_on = wx.CheckBox(self, label=enable_label)
        self.chk_on.SetValue(bool(enabled))

        self.colour = wx.ColourPickerCtrl(self, colour=wx.Colour(color))
        self.levels = wx.SpinCtrl(self, min=1, max=40, initial=int(levels))
        self.width = wx.SpinCtrlDouble(self, min=0.1, max=6.0, inc=0.1,
                                        initial=float(width),
                                        style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.width.SetDigits(1)
        for ctrl in (self.levels, self.width):
            theme.theme_control(ctrl)
        # Enter commits the typed number instead of being swallowed (or, in a
        # dialog, triggering the default button with the old value still set).
        theme.bind_spin_enter(self.width)

        grid = wx.FlexGridSizer(cols=2, vgap=6, hgap=8)
        grid.AddGrowableCol(1)
        for label, ctrl in (("Colour:", self.colour),
                            ("Number of levels:", self.levels),
                            ("Line thickness:", self.width)):
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(ctrl, 0, wx.EXPAND)

        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(self.chk_on, 0, wx.ALL, 10)
        outer.Add(grid, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        btns = self.CreateButtonSizer(wx.OK | wx.CANCEL)
        if btns is not None:
            outer.Add(btns, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizerAndFit(outer)

        self.chk_on.Bind(wx.EVT_CHECKBOX, self._on_toggle)
        self._on_toggle(None)
        theme.apply_theme_to_window(self, theme.get_theme())

    def _on_toggle(self, event):
        """The three appearance controls mean nothing while the contours
        are switched off -- grey them out rather than letting the user set
        values that will not be used."""
        on = self.chk_on.GetValue()
        for ctrl in (self.colour, self.levels, self.width):
            ctrl.Enable(on)

    def get_enabled(self):
        return self.chk_on.GetValue()

    def get_color(self):
        return self.colour.GetColour().GetAsString(wx.C2S_HTML_SYNTAX)

    def get_levels(self):
        return int(self.levels.GetValue())

    def get_width(self):
        return float(self.width.GetValue())


class ResidualPanel(wx.Panel):
    """A single row of residual maps -- experiment minus amplitude-scaled
    simulation, over each active dataset's [fmin, fmax] region -- for one
    grid point. Sits above the fitness heatmap, hidden until the user
    toggles it on, and is repopulated by double-clicking the heatmap.

    Residuals are signed (positive = the simulation under-predicts that
    peak, negative = it puts intensity where the experiment has none), so
    these are drawn with a diverging colormap on limits made symmetric
    about zero rather than with the heatmap's own colormap -- otherwise
    zero residual would land at an arbitrary colour and the sign, which is
    the whole diagnostic value of a residual map, would be unreadable.
    """

    # Fixed full scale: +-2 is twice a dataset's typical experimental
    # intensity, since residual maps are always drawn normalized (they are a
    # display, not the statistic -- see GridSearchFrame._render_residuals).
    # Fixed rather than autoscaled so maps stay comparable between datasets
    # and between grid points.
    RESIDUAL_SCALE_MIN = -2.0
    RESIDUAL_SCALE_MAX = 2.0

    def __init__(self, parent):
        super().__init__(parent)
        self.figure = Figure(figsize=(6, 3.0), dpi=100)
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetMinSize((-1, 300))

        self.cmap = cm.get_cmap('RdBu_r')

        # Simulation-contour overlay settings, edited by right-clicking the
        # row (see _on_click). Kept here rather than in the frame because
        # they are purely a property of this display.
        self.contour_on = True
        self.contour_color = '#000000'
        self.contour_n = 6
        self.contour_width = 0.8

        # Last drawn content, so a settings change can redraw without the
        # user having to double-click the same grid point again.
        self._last_entries = []
        self._last_subtitle = ""
        self._last_value_label = "RMSD"

        self._redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self._redraw_debouncer.on_size)
        self.canvas.mpl_connect('button_press_event', self._on_click)
        self.show_message("Double-click a point on the fitness map to show "
                          "its residual maps.\n"
                          "Right-click here for simulation-contour settings.")

    # -- simulation contour settings ---------------------------------------
    def _on_click(self, event):
        if event.button == 3:
            self.edit_contour_settings()

    def edit_contour_settings(self):
        with ContourSettingsDialog(
                self, "Simulation Contours", self.contour_on,
                self.contour_color, self.contour_n, self.contour_width,
                enable_label="Draw simulation contours") as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            self.contour_on = dlg.get_enabled()
            self.contour_color = dlg.get_color()
            self.contour_n = dlg.get_levels()
            self.contour_width = dlg.get_width()
        if self._last_entries:
            self.show_maps(self._last_entries, self._last_subtitle,
                           self._last_value_label)

    # -- drawing ------------------------------------------------------------
    @staticmethod
    def _layout_margins(has_subtitle):
        """Constant subplot margins. Chosen once, by eye, to clear two-line
        titles above and tick labels below at this row's height, with
        wspace wide enough for each map's colorbar and its labels."""
        return {'left': 0.055, 'right': 0.935, 'bottom': 0.17,
                'top': 0.78 if has_subtitle else 0.86, 'wspace': 0.55}

    def show_message(self, text):
        """Blank the row and print a single centred note (no stored data,
        nothing selected yet, ...)."""
        self._last_entries = []
        self._last_subtitle = ""
        self.figure.clear()
        theme.style_figure(self.figure, [])
        self.figure.text(0.5, 0.5, text, ha='center', va='center',
                         fontsize='small',
                         color=theme.theme_colors()["ax_fg"], wrap=True)
        self.canvas.draw_idle()

    def _draw_sim_contours(self, ax, entry):
        """Overlay the amplitude-scaled simulation as contour lines, so the
        residual underneath can be read against where the simulation
        actually put intensity. Levels are evenly spaced across the
        simulation's own range rather than left to matplotlib's automatic
        choice, so the requested number of lines is what gets drawn."""
        sim = entry.get('sim')
        if sim is None or sim.size == 0:
            return
        lo, hi = float(np.nanmin(sim)), float(np.nanmax(sim))
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return
        levels = np.linspace(lo, hi, self.contour_n + 2)[1:-1]
        ax.contour(entry['x'], entry['y'], sim.T, levels=levels,
                   colors=self.contour_color, linewidths=self.contour_width,
                   alpha=0.9, zorder=3)

    def show_maps(self, entries, subtitle="", value_label="RMSD"):
        """entries: list of dicts with 'title' (what to head the map with --
        the field value, matching the main window's own B0 titles), 'resid'
        and 'sim' (2-D arrays, axis 0 = x as everywhere else in this
        codebase), 'x'/'y' (the region's frequency axes), 'extent'
        ([x0, x1, y0, y1]) and 'rmsd' (that dataset's own RMSD here)."""
        if not entries:
            self.show_message("This grid point has no successful simulation "
                              "for any active dataset.")
            return

        self._last_entries = entries
        self._last_subtitle = subtitle
        self._last_value_label = value_label
        self.figure.clear()
        axes = self.figure.subplots(1, len(entries), squeeze=False)[0]
        for ax, entry in zip(axes, entries):
            resid = entry['resid']
            # Fixed colour limits rather than per-map autoscaling, so the
            # maps are directly comparable with each other and from one
            # grid point to the next. Residuals always arrive normalized,
            # so +-2 means "twice the typical experimental intensity of
            # this dataset" for every map regardless of display mode.
            im = ax.imshow(resid.T, extent=entry['extent'], origin='lower',
                           aspect='auto', cmap=self.cmap,
                           interpolation='nearest',
                           vmin=self.RESIDUAL_SCALE_MIN, vmax=self.RESIDUAL_SCALE_MAX)
            if self.contour_on:
                self._draw_sim_contours(ax, entry)
            # Square the axes box explicitly rather than via aspect='equal':
            # both frequency axes span the same [fmin, fmax] here, so the
            # two agree in practice, but set_box_aspect keeps the panel
            # square even if a dataset's x and y sampling differ slightly.
            ax.set_box_aspect(1)
            ax.set_title(f"{entry['title']}\n{value_label}={entry['rmsd']:.4g}",
                         fontsize='x-small')
            ax.tick_params(labelsize='xx-small')
            cbar = self.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize='xx-small',
                                colors=theme.theme_colors()["ax_fg"])
            # Attaching a colorbar anchors its parent east so the two stay
            # together; combined with the square box that pins each map to
            # the right of its column and leaves a lopsided gap (glaring
            # with one dataset). Re-centre it in its own column.
            ax.set_anchor('C')
        theme.style_figure(self.figure, list(axes))
        if subtitle:
            self.figure.suptitle(subtitle, fontsize='small',
                                 color=theme.theme_colors()["ax_fg"])

        # Deliberately NOT tight_layout: it re-measures the tick labels and
        # titles on every draw, so switching units (say "RMSD" to
        # "Normalized RMSD", or 250.0 to 1.2e-4 on a colorbar) changed the
        # reserved text width, which changed each slot, which shifted the
        # squared axes sideways -- the row appeared to creep every time
        # anything was touched. Fixed margins make the geometry depend only
        # on the number of datasets, so a redraw with different numbers in
        # it lands in exactly the same place.
        self.figure.subplots_adjust(**self._layout_margins(bool(subtitle)))
        self.canvas.draw_idle()


# --------------------------------------------------------------------------- #
# Main grid-search window
# --------------------------------------------------------------------------- #
class GridSearchFrame(wx.Frame):

    # Display statistics offered by the "Display:" pulldown, in menu order.
    DISPLAY_RMSD = 0
    DISPLAY_NORMALIZED = 1
    DISPLAY_CREDIBLE = 2
    DISPLAY_DENSITY = 3
    DISPLAY_LABELS = ["RMSD", "Fit Quality, norm(RMSD)",
                      "Probability (credible level)", "Probability (density)"]
    DISPLAY_VALUE_LABELS = ["RMSD", "Fit Quality",
                            "Credible level", "Posterior density"]
    PROBABILITY_MODES = (DISPLAY_CREDIBLE, DISPLAY_DENSITY)

    # Every stored simulation is one 64-bit float per experimental point in
    # the region of interest (the interpolator output is real --
    # HYSCOREsim.Spectrum is an np.abs magnitude spectrum).
    STORE_BYTES_PER_POINT = 8
    # Above this predicted total the user is asked whether to keep storing.
    STORE_WARN_BYTES = 0.5 * 1024 ** 3

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
        self._last_results_norm = None
        self._last_results_n = None
        self._last_results_sse = None
        self._last_results_nk = None
        self._last_n_datasets = 0
        # (sigma, N_k, null_rms_k) per usable dataset, captured when the scan
        # starts. Kept separately from _roi_cache, which is released whenever
        # simulations aren't being stored -- the noise level has to outlive it
        # so switching CI method or display can re-derive the interval without
        # a rescan.
        self._last_noise = None
        self._last_grid1 = None
        self._last_grid2 = None
        self._last_N = 0
        # coarse (i, j) of the grid point whose residuals are on screen, so
        # switching RMSD units can redraw it in the new units
        self._last_click_ij = None

        # Interpolated simulations kept from the last scan, when the user
        # asked for them: {(i, j): [scaled_sim_or_None per active dataset]},
        # indexed by *coarse* grid indices. _roi_cache holds the matching
        # experimental regions (one per dataset, constant across the scan),
        # so a residual map is just roi['exp'] - stored_sim. Both are
        # dropped at the start of every new scan so the old scan's arrays
        # are freed before the new one starts allocating.
        self._sim_store = {}
        self._roi_cache = None
        self._show_residuals = False

        # Reference contours drawn over the fitness map; edited through the
        # same dialog the residual row uses. Defaults reproduce what the
        # map drew before this was configurable.
        self.map_contour_on = True
        self.map_contour_color = '#ffffff'
        self.map_contour_n = 9
        self.map_contour_width = 0.5


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

        # theme.StatusStrip rather than a native wx.StatusBar, whose field
        # text MSW draws in the system colour regardless of the theme.
        self.statusbar = theme.StatusStrip(self)
        self.statusbar.SetStatusText("Select two Sys parameters, set their ranges, then Run.")

        splitter = wx.SplitterWindow(self, style=wx.SP_3D)
        left = wx.Panel(splitter)

        # Right-hand side is a column: the (initially hidden) row of
        # residual maps on top, the fitness heatmap filling the rest.
        right = wx.Panel(splitter)
        self.residuals = ResidualPanel(right)
        self.heatmap = HeatmapPanel(right)
        self.heatmap.set_colormap(self.cmap)
        self.heatmap.on_hover_callback = self._on_hover
        self.heatmap.on_rightclick_callback = self._on_rightclick
        self.heatmap.on_doubleclick_callback = self._on_doubleclick

        rightsizer = wx.BoxSizer(wx.VERTICAL)
        rightsizer.Add(self.residuals, 0, wx.EXPAND)
        rightsizer.Add(self.heatmap, 1, wx.EXPAND)
        right.SetSizer(rightsizer)
        self.residuals.Hide()
        self.right_panel = right

        splitter.SplitVertically(left, right, sashPosition=400)
        splitter.SetMinimumPaneSize(260)
        # Give the left pane a share of any extra width, so the parameter
        # range boxes grow with the window instead of staying pinned at
        # their starting size, and relayout live while the sash is dragged.
        splitter.SetSashGravity(0)
        splitter.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self._on_sash_moved)
        splitter.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self._on_sash_moved)
        self.splitter = splitter

        leftsizer = wx.BoxSizer(wx.VERTICAL)

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
        self.btn_map_contours = PillButton(left, "Contours…", color_key="pill_calm")
        self.btn_map_contours.SetToolTip(
            "Colour, number of levels and line thickness of the reference "
            "contours drawn over the fitness map, and whether to draw them "
            "at all -- the same controls the residual row offers on "
            "right-click.\n"
            "The level count also sets how many lines the \"RMSD as "
            "contours\" display mode draws.")
        self.btn_cmap = PillButton(left, "Colormap…", color_key="accent")
        self.cmap_swatch = wx.Panel(left, size=(70, 22), style=wx.BORDER_SIMPLE)
        self.chk_invert_cmap = wx.CheckBox(left, label="Invert")
        cmap_sizer.Add(self.btn_map_contours, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        cmap_sizer.Add(self.btn_cmap, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        cmap_sizer.Add(self.cmap_swatch, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        cmap_sizer.Add(self.chk_invert_cmap, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(cmap_sizer, 0, wx.ALL, 4)

        analysis_sizer = wx.BoxSizer(wx.HORIZONTAL)
        analysis_sizer.Add(wx.StaticText(left, label="Smooth: "),
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
        analysis_sizer.Add(wx.StaticText(left, label="CI :"),
                           0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.ci_method_choice = wx.Choice(
            left, choices=["F-test", "χ²",
                           "χ²/w noise"])
        self.ci_method_choice.SetSelection(0)
        self.ci_method_choice.SetToolTip(
            "F-test: standard when the noise level isn't independently known\n"
            "χ²: Delta-χ² using the fit's own RMSD as the noise estimate\n"
            "χ² /w noise: the same Delta-χ² , but with the noise level for each spectrum")
        analysis_sizer.Add(self.ci_method_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        self.chk_use_noise = wx.CheckBox(left, label="Use exp. noise")
        self.chk_use_noise.SetValue(True)
        self.chk_use_noise.SetToolTip(
            "Use the noise level pyHYSCORE measured on each spectrum's far "
            "high-frequency corner (where a HYSCORE spectrum carries no "
            "signal) as an independent sigma.\n"
            "On: enables the measured-noise chi-square, and reports the "
            "reduced χ² of the best fit -- whether the residual is at the "
            "noise level or well above it.\n"
            "Off: the noise estimate is ignored entirely and the interval is "
            "whatever the F-test or simple chi-square says, as before.")
        analysis_sizer.Add(self.chk_use_noise, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(analysis_sizer, 0, wx.ALL, 4)

        display_sizer = wx.BoxSizer(wx.HORIZONTAL)
        display_sizer.Add(wx.StaticText(left, label="Display:"),
                          0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.display_choice = wx.Choice(left, choices=self.DISPLAY_LABELS)
        self.display_choice.SetSelection(self.DISPLAY_RMSD)
        self.display_choice.SetToolTip(
            "RMSD: the raw pooled fitness, in experimental intensity units.\n"
            "Normalized RMSD: each dataset scored against its own null model "
            "(the residual left by no simulation at all, which the "
            "non-negative amplitude fit can never do worse than), averaged by "
            "point count. 0 = perfect, 1 = explains nothing; 0.3 means 91% of "
            "the experimental power is accounted for.\n"
            "Probability (credible level): each cell labelled with the "
            "smallest credible region containing it, so 0.95 is the boundary "
            "of the 95% region. Built by marginalizing an unknown noise level "
            "per dataset, which makes it independent of dataset brightness "
            "and of the normalization choice above.\n"
            "Probability (density): the posterior density itself. The primary "
            "object, but a sharp posterior puts nearly every cell at "
            "effectively zero.\n"
            "This chooses the fitness map only -- residual maps are always "
            "drawn in normalized units on a fixed -2 to +2 scale, so they "
            "stay comparable between datasets and between grid points.")
        display_sizer.Add(self.display_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)

        self.lbl_divisor = wx.StaticText(left, label="Oversampling:")
        display_sizer.Add(self.lbl_divisor, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.neff_divisor = wx.SpinCtrlDouble(left, min=1.0, max=100000.0,
                                               inc=1.0, initial=1.0,
                                               style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.neff_divisor.SetDigits(2)
        self.neff_divisor.SetToolTip(
            "Spectrum points per independent measurement.\n"
            "Oversampling = 4^(zero fill), i.e. 2^n per axis, squared.\n"
            "At 1 the credible regions come out too narrow.\n"
            "Applies to the probability maps and to the measured-noise "
            "chi-square, both of which would otherwise count correlated "
            "points as independent ones.")
        self.neff_divisor.SetValue(self._default_oversampling())
        theme.theme_control(self.neff_divisor)
        display_sizer.Add(self.neff_divisor, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(display_sizer, 0, wx.ALL, 4)

        rmsd_contour_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_rmsd_contours = wx.CheckBox(left, label="RMSD as contours")
        self.chk_rmsd_contours.SetToolTip(
            "Draw the fitness surface as contour lines coloured by the "
            "colormap instead of as a filled density plot. Uses the same "
            "10%-step levels the density plot draws as thin reference "
            "lines.")
        rmsd_contour_sizer.Add(self.chk_rmsd_contours, 0,
                               wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        rmsd_contour_sizer.Add(wx.StaticText(left, label="Line thickness:"),
                               0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.rmsd_contour_width = wx.SpinCtrlDouble(
            left, min=0.1, max=6.0, inc=0.1, initial=1.2,
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.rmsd_contour_width.SetDigits(1)
        theme.theme_control(self.rmsd_contour_width)
        rmsd_contour_sizer.Add(self.rmsd_contour_width, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(rmsd_contour_sizer, 0, wx.ALL, 4)

        resid_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_store_sims = wx.CheckBox(
            left, label="Keep simulations in memory")
        self.chk_store_sims.SetValue(True)
        self.chk_store_sims.SetToolTip(
            "Keep every grid point's simulation (interpolated onto each "
            "dataset's [fmin, fmax] region only) so residual maps can be "
            "shown afterwards without re-simulating.\n"
            "You are warned before the scan starts if this would need more "
            f"than {self.format_bytes(self.STORE_WARN_BYTES)}.")
        resid_sizer.Add(self.chk_store_sims, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)
        self.btn_residuals = PillButton(left, "Show Residual Maps",
                                         color_key="pill_calm")
        self.btn_residuals.SetToolTip(
            "Show/hide the row of per-dataset residual maps above the "
            "fitness map. While shown, double-click anywhere on the "
            "fitness map to load that grid point's residuals.")
        resid_sizer.Add(self.btn_residuals, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(resid_sizer, 0, wx.ALL, 4)

        self.lbl_ci = wx.StaticText(left, label="")
        # Same face and size as every other label in this panel, just bold --
        # it used to be a fixed-pitch face that sat oddly against the rest.
        self.lbl_ci.SetFont(
            wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT).Bold())
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

        # The splitter used to be this frame's only child and so filled it
        # automatically; with the status strip below it they need a sizer.
        framesizer = wx.BoxSizer(wx.VERTICAL)
        framesizer.Add(splitter, 1, wx.EXPAND)
        framesizer.Add(self.statusbar, 0, wx.EXPAND)
        self.SetSizer(framesizer)

        self.btn_refresh.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.btn_cmap.Bind(wx.EVT_BUTTON, self.on_choose_cmap)
        self.btn_map_contours.Bind(wx.EVT_BUTTON, self.on_edit_map_contours)
        self.chk_invert_cmap.Bind(wx.EVT_CHECKBOX, self.on_invert_cmap_toggle)
        self.interp_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        self.ci_method_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        self.chk_use_noise.Bind(wx.EVT_CHECKBOX, self.on_analysis_choice_changed)
        self.btn_residuals.Bind(wx.EVT_BUTTON, self.on_toggle_residuals)
        self.chk_rmsd_contours.Bind(wx.EVT_CHECKBOX, self.on_analysis_choice_changed)
        self.display_choice.Bind(wx.EVT_CHOICE, self.on_display_mode_changed)
        self.neff_divisor.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_display_mode_changed)
        self.neff_divisor.Bind(wx.EVT_KILL_FOCUS, self.on_display_mode_changed)
        self.rmsd_contour_width.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_analysis_choice_changed)
        # SpinCtrlDouble only fires EVT_SPINCTRLDOUBLE for the arrows; typing
        # a value and tabbing away emits a plain text event. Enter is handled
        # by theme.bind_spin_enter, since wx.EVT_TEXT_ENTER doesn't reliably
        # reach these compound controls (see that function).
        theme.bind_spin_enter(self.neff_divisor,
                              lambda: self.on_display_mode_changed(None))
        theme.bind_spin_enter(self.rmsd_contour_width,
                              lambda: self.on_analysis_choice_changed(None))
        self.rmsd_contour_width.Bind(wx.EVT_KILL_FOCUS, self.on_analysis_choice_changed)
        self.btn_run.Bind(wx.EVT_BUTTON, self.on_run_or_cancel)
        self.param_pg.Bind(wxpg.EVT_PG_CHANGED, self.on_param_toggle)
        self.Bind(wx.EVT_CLOSE, self.on_close)

        self._load_sys_snapshot()
        self._build_param_grid()
        self._paint_swatch()
        self._sync_display_controls()
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

    def on_edit_map_contours(self, event):
        """Same dialog the residual row uses, applied to the fitness map's
        reference contours."""
        with ContourSettingsDialog(
                self, "Fitness Map Contours", self.map_contour_on,
                self.map_contour_color, self.map_contour_n,
                self.map_contour_width,
                enable_label="Draw contours on the fitness map") as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            self.map_contour_on = dlg.get_enabled()
            self.map_contour_color = dlg.get_color()
            self.map_contour_n = dlg.get_levels()
            self.map_contour_width = dlg.get_width()
        self._refresh_display()

    def _on_sash_moved(self, event):
        """Keep the parameter range boxes tracking the left panel's width
        as the sash is dragged, rather than only once the drag ends."""
        event.Skip()
        self._relayout_left()

    def on_invert_cmap_toggle(self, event):
        self.cmap_inverted = self.chk_invert_cmap.GetValue()
        self._apply_effective_cmap()

    def _apply_effective_cmap(self):
        self.cmap = self.cmap_base.reversed() if self.cmap_inverted else self.cmap_base
        self.heatmap.set_colormap(self.cmap)
        self._paint_swatch()
        # In contour display mode recolouring the mappable in place does not
        # restyle the drawn lines, so redraw from the cached results (cheap
        # -- no re-scanning; see _refresh_display).
        if self._last_results is not None and self.chk_rmsd_contours.GetValue():
            self._refresh_display()

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
        flat = [p for p in self.flatten_sys_params(self.sys_snapshot)
                if self.qualified_label(p) not in self._driven_sys_labels]
        flat += self.flatten_func_vars(self.func_var_snapshot)
        self._param_by_name = {self.param_name(p): p for p in flat}
        # drop selections whose parameter no longer exists (e.g. nucleus removed)
        self._selected = [n for n in self._selected if n in self._param_by_name]

        cur_cat = None
        cat_prop = None
        for p in flat:
            if p['cat'] != cur_cat:
                cur_cat = p['cat']
                cat_prop = self.param_pg.Append(wxpg.PropertyCategory(cur_cat))
            pname = self.param_name(p)
            disp = f"{self.param_leaf_label(p)} = {p['value']:.5g}"
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
        # The Sys marker is read from the main window, so re-pulling the
        # parameters is exactly when it may have moved.
        self._refresh_display()

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
                wx.StaticText(self.range_panel, label=self.param_full_label(p)),
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
                        f"{self.param_full_label(self._param_by_name[pname])}.")
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

        # Slice each dataset's experimental region once, up front: the scan
        # reuses it at every grid point instead of re-slicing, and its size
        # is what determines the cost of keeping simulations.
        roi_cache = [self.dataset_roi(dd) for dd in active_data]
        self._last_noise = [(roi['sigma'], roi['n'], roi['null_rms'])
                            for roi in roi_cache if roi is not None]

        store_sims = self.chk_store_sims.GetValue()
        if store_sims:
            # predicted cost of keeping every grid point's interpolated
            # simulation: one float per region-of-interest point, per dataset
            per_point = sum(roi['exp'].size for roi in roi_cache if roi is not None)
            est = per_point * len(grid1) * len(grid2) * self.STORE_BYTES_PER_POINT
            if est > self.STORE_WARN_BYTES:
                answer = wx.MessageBox(
                    f"Keeping every simulation from this scan would need about "
                    f"{self.format_bytes(est)} of memory\n"
                    f"({len(grid1)} x {len(grid2)} grid points x "
                    f"{len(active_data)} dataset(s), interpolated onto the "
                    f"[fmin, fmax] region only).\n\n"
                    f"Yes  - store anyway (residual maps available)\n"
                    f"No   - run without storing (residual maps unavailable)\n"
                    f"Cancel - don't run at all",
                    "Grid Search - large memory requirement",
                    wx.YES_NO | wx.CANCEL | wx.ICON_WARNING)
                if answer == wx.CANCEL:
                    return
                store_sims = (answer == wx.YES)

        # Release the previous scan's stored arrays *before* the new scan
        # starts allocating, so the two sets never coexist.
        self._sim_store = {}
        self._roi_cache = roi_cache if store_sims else None
        # The residual marker indexes a grid that is about to change; the
        # Sys marker needs no clearing, being read live at draw time.
        self._last_click_ij = None
        # Reset the row even while it is hidden, so toggling it on later
        # can never surface the previous scan's residuals.
        self.residuals.show_message(
            "Scan running -- double-click a grid point when it finishes."
            if store_sims else
            "Simulations are not being stored for this scan.")

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
        self.chk_store_sims.Enable(False)
        total = len(grid1) * len(grid2)
        self.gauge.SetRange(total)
        self.gauge.SetValue(0)
        self.lbl_progress.SetLabel(f"0 / {total}")
        self.statusbar.SetStatusText("Running grid search...")

        self._worker_thread = threading.Thread(
            target=self._worker_run,
            args=(param1, param2, grid1, grid2, sys_template, exp_template,
                  opt_template, active_data, self._cancel_event,
                  func_rows, base_variable_values, roi_cache, store_sims),
            daemon=True)
        self._worker_thread.start()

    def _worker_run(self, param1, param2, grid1, grid2, sys_template,
                     exp_template, opt_template, active_data, cancel_event,
                     func_rows, base_variable_values, roi_cache, store_sims):
        n1, n2 = len(grid1), len(grid2)
        n_ds = len(active_data)
        results = np.full((n1, n2), np.nan)
        results_norm = np.full((n1, n2), np.nan)
        results_n = np.zeros((n1, n2), dtype=int)
        # Per-dataset residuals, needed by the posterior. Two floats per
        # dataset per grid point -- negligible next to anything else here.
        results_sse = np.full((n1, n2, n_ds), np.nan)
        results_nk = np.zeros((n1, n2, n_ds), dtype=int)
        hs = HYSCOREsim(errorFunc=self._silent_error_func)
        hs.verbose = False
        start = time.perf_counter()
        count = 0
        total = n1 * n2
        n_failed_points = 0
        error_log = {}   # unique error message -> occurrence count; also printed live
        sim_store = {} if store_sims else None

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
                        self.apply_param_value(sys_dict, p, v)
                slot = [None] * len(active_data) if store_sims else None
                stats = [None] * len(active_data)
                try:
                    func_ok = self.apply_functions_to_sys_dict(
                        sys_dict, func_rows, variable_values, error_log=error_log)
                    if not func_ok:
                        rmsd, n_pts, rmsd_norm = np.nan, 0, np.nan
                    else:
                        Sys = sysPar()
                        Sys.setFromCtrl(sys_dict)
                        rmsd, n_pts, rmsd_norm = self.run_gridpoint(
                            hs, Sys, exp_template, opt_template,
                            active_data, error_log=error_log,
                            roi_cache=roi_cache, sim_store=slot,
                            stats_out=stats)
                except Exception as e:
                    self._log_gridpoint_error(error_log, "Sys.setFromCtrl", e)
                    rmsd, n_pts, rmsd_norm = np.nan, 0, np.nan
                # Only points that actually produced a simulation are worth
                # a dict entry -- a failed point would otherwise cost a
                # list of Nones for nothing.
                if slot is not None and any(s is not None for s in slot):
                    sim_store[(i, j)] = slot
                if np.isnan(rmsd):
                    n_failed_points += 1
                results[i, j] = rmsd
                results_norm[i, j] = rmsd_norm
                results_n[i, j] = n_pts
                for k, st in enumerate(stats):
                    if st is not None:
                        results_sse[i, j, k], results_nk[i, j, k] = st
                count += 1
                elapsed = time.perf_counter() - start
                wx.CallAfter(self._on_progress, count, total, elapsed)

        wx.CallAfter(self._on_finished, results, results_n, grid1, grid2, param1, param2,
                     cancel_event.is_set(), n_failed_points, total, error_log, sim_store,
                     results_norm, len(active_data), results_sse, results_nk)

    def _on_progress(self, count, total, elapsed):
        if self._closing:
            return
        self.gauge.SetValue(count)
        rate = elapsed / count if count else 0.0
        remaining = rate * (total - count)
        self.lbl_progress.SetLabel(
            f"{count} / {total}   elapsed {elapsed:.0f}s   ETA {remaining:.0f}s")

    def _on_finished(self, results, results_n, grid1, grid2, param1, param2, cancelled,
                      n_failed, total, error_log=None, sim_store=None,
                      results_norm=None, n_datasets=0,
                      results_sse=None, results_nk=None):
        if self._closing:
            # The frame is going away; don't take a reference to the
            # worker's stored arrays, just let them be collected.
            return
        self._running = False
        self.btn_run.SetLabel("Run Grid Search")
        self.btn_run.SetColorKey("accent")
        self.btn_run.Enable()
        self.param_pg.Enable(True)
        self.btn_refresh.Enable()
        self.range_panel.Enable(True)
        self.btn_cmap.Enable()
        self.chk_store_sims.Enable(True)
        self._sync_param_lock_state()

        self._sim_store = sim_store if sim_store else {}
        if not self._sim_store:
            self._roi_cache = None
        self._last_click_ij = None
        self.residuals.show_message(
            "Double-click a point on the fitness map to show its "
            "residual maps." if self._sim_store else
            "No simulations were stored for this scan.")

        self._last_results = results
        self._last_results_norm = results_norm
        self._last_results_n = results_n
        self._last_results_sse = results_sse
        self._last_results_nk = results_nk
        self._last_n_datasets = n_datasets
        self._last_grid1 = grid1
        self._last_grid2 = grid2
        self._last_param1 = param1
        self._last_param2 = param2
        self._refresh_display()

        msg = "Grid search cancelled." if cancelled else "Grid search finished."
        if n_failed:
            msg += f" {n_failed}/{total} grid points had no successful simulation (see console)."
        if self._sim_store:
            held = sum(s.nbytes for slot in self._sim_store.values()
                       for s in slot if s is not None)
            msg += (f" Simulations kept for {len(self._sim_store)} grid points "
                    f"({self.format_bytes(held)}).")
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
        # EVT_KILL_FOCUS is one of the sources here and must keep
        # propagating, or focus handling for the control breaks.
        if event is not None:
            event.Skip()
        # The CI method and the "use measured noise" checkbox both decide
        # whether the oversampling divisor is doing anything, so the control
        # states have to follow them and not only the Display pulldown.
        self._sync_display_controls()
        self._refresh_display()

    def _refresh_display(self):
        """(Re)draw the heatmap from the coarse results cached by the last
        scan, in whichever statistic the Display pulldown selects, and
        recompute the interval readout to match. Cheap enough to call on
        every control change -- no re-scanning involved."""
        if self._last_results is None:
            return
        mode = self.display_mode()
        drew = False
        if mode in self.PROBABILITY_MODES:
            drew = self._display_probability(mode)
        if not drew:
            # Probability needs per-dataset residuals; fall back to the raw
            # fitness map if they are missing or the posterior degenerates.
            self._display_fitness(self.DISPLAY_RMSD if mode in self.PROBABILITY_MODES else mode)
        self.left_panel.Layout()
        self.Layout()

    def _background_levels(self, vmin, vmax):
        """Evenly spaced reference contour levels across the displayed
        surface's range, as many as the contour dialog asks for. Empty when
        the surface is flat, or when the overlay is switched off *and* the
        map is not itself drawn as contours -- in contour display mode
        these levels are the plot, so they are always needed."""
        if not (vmax > vmin):
            return []
        if not self.map_contour_on and not self.chk_rmsd_contours.GetValue():
            return []
        return list(np.linspace(vmin, vmax, self.map_contour_n + 2)[1:-1])

    def _current_sys_point(self):
        """Where the main window's Sys currently sits, in the plane of the
        two scanned parameters -- read live rather than remembered, so it
        follows whatever the main grid holds now, however it got there
        (loaded from this window, typed in by hand, driven by a Function
        Modifier). Returns None if either parameter no longer exists."""
        if self._last_param1 is None or self._last_param2 is None:
            return None
        try:
            params = self.main.tabulated_panel.Sys_param.parameters
        except AttributeError:
            return None
        v1 = self.read_param_value(params, self._last_param1)
        v2 = self.read_param_value(params, self._last_param2)
        if v1 is None or v2 is None:
            return None
        return (v1, v2)

    def _marker_points(self):
        """(residual-shown point, current-Sys point) in data coordinates,
        either possibly None.

        Both are dropped when they fall outside the axes currently being
        drawn. Markers are plotted with ax.plot, which feeds the autoscaler,
        so a point outside the scanned range -- the Sys sitting somewhere
        this scan never visited, say -- would otherwise drag the axis limits
        out to reach it and wreck the scale."""
        if self._last_grid1 is None or self._last_grid2 is None:
            return None, None
        lo1, hi1 = self._last_grid1[0], self._last_grid1[-1]
        lo2, hi2 = self._last_grid2[0], self._last_grid2[-1]

        def in_range(pt):
            if pt is None:
                return None
            x, y = pt
            if min(lo1, hi1) <= x <= max(lo1, hi1) and \
               min(lo2, hi2) <= y <= max(lo2, hi2):
                return (x, y)
            return None

        resid_pt = None
        if self._last_click_ij is not None:
            i, j = self._last_click_ij
            if i < len(self._last_grid1) and j < len(self._last_grid2):
                resid_pt = (self._last_grid1[i], self._last_grid2[j])
        return in_range(resid_pt), in_range(self._current_sys_point())

    def _display_fitness(self, mode):
        """Draw the RMSD or normalized-RMSD map, with the F-test or
        chi-square confidence region."""
        grid1, grid2 = self._last_grid1, self._last_grid2
        param1, param2 = self._last_param1, self._last_param2
        label1, label2 = self.param_full_label(param1), self.param_full_label(param2)

        # Everything downstream -- interpolation, the minimum, the
        # confidence region, the contours -- runs on whichever statistic is
        # selected, so the whole display stays self-consistent.
        normalized = (mode == self.DISPLAY_NORMALIZED)
        results = self._last_results_norm if normalized else self._last_results
        if results is None:
            results = self._last_results
            normalized = False
        value_label = "Normalized RMSD" if normalized else "RMSD"

        as_contours = self.chk_rmsd_contours.GetValue()
        rmsd_cw = float(self.rmsd_contour_width.GetValue())

        method = 'cubic' if self.interp_choice.GetSelection() == 1 else 'linear'
        out = self._interpolate_grid(method, grid1, grid2, results)
        if out is None:
            self.heatmap.show_results(results, grid1, grid2, label1, label2,
                                      as_contours=as_contours,
                                      rmsd_contour_width=rmsd_cw,
                                      value_label=value_label)
            self.lbl_ci.SetLabel(
                "No successful grid points -- cannot interpolate or estimate confidence intervals.")
            return
        fine1, fine2, fine_vals, used_method = out

        bi, bj = np.unravel_index(np.nanargmin(fine_vals), fine_vals.shape)
        rmsd_min, best1, best2 = fine_vals[bi, bj], fine1[bi], fine2[bj]

        # Thin reference contours at 10% steps of the interpolated surface's
        # own range -- purely a visual aid for reading the landscape, drawn
        # under the (bold) confidence-interval contour.
        vmin, vmax = float(np.nanmin(fine_vals)), float(np.nanmax(fine_vals))
        background_levels = self._background_levels(vmin, vmax)

        ci_method = ('fstat', 'chi2', 'noise')[max(self.ci_method_choice.GetSelection(), 0)]
        # Degrees of freedom come from the coarse point that is best under
        # the *displayed* statistic -- with per-dataset normalization the
        # two statistics can put their minimum in different places.
        if np.any(np.isfinite(results)) and self._last_results_n is not None:
            ci, cj = np.unravel_index(np.nanargmin(results), results.shape)
            N = int(self._last_results_n[ci, cj])
        else:
            N = 0
        self._last_N = N

        # The measured noise level, in the same units as the displayed
        # statistic, and the residual it has to be compared against. Both are
        # needed for the measured-noise region and for the reduced chi^2, which
        # is worth reporting whichever method drew the contour -- it is the one
        # number that says whether the fit has reached the noise or stalled
        # well above it.
        divisor = float(self.neff_divisor.GetValue())
        use_noise = self.chk_use_noise.GetValue()
        sigma = self.pooled_sigma(normalized) if use_noise else None
        sse_min = N * rmsd_min ** 2 if N > 0 else None
        chi2_red = None
        if sigma and sse_min and (N / max(divisor, 1e-9)) > 2:
            chi2_red = sse_min / (sigma ** 2 * (N / max(divisor, 1e-9) - 2))
        # Falling back silently would leave the panel claiming a measured-noise
        # interval it never computed, so say which of the two reasons it was.
        noise_note = None
        if ci_method == 'noise' and sigma is None:
            noise_note = ('("Use measured noise" is off -- used the F-test instead)'
                          if not use_noise else
                          "(no measured noise level for every dataset -- re-run the "
                          "FFT, or load data saved with one; used the F-test instead)")
            ci_method = 'fstat'

        CONF = 0.95
        factor = self._rmsd_threshold_factor(ci_method, CONF, N, p=2, sigma=sigma,
                                             sse_min=sse_min, divisor=divisor)
        contour_levels = []
        half_width = {1: None, 2: None}
        if factor is not None:
            thresh = rmsd_min * factor
            if np.nanmin(fine_vals) <= thresh:
                contour_levels.append((thresh, "95%"))
                profile1 = np.nanmin(fine_vals, axis=1)  # min over axis2, vs. axis1
                profile2 = np.nanmin(fine_vals, axis=0)  # min over axis1, vs. axis2
                bounds1 = self._axis_confidence_bounds(fine1, profile1, thresh)
                bounds2 = self._axis_confidence_bounds(fine2, profile2, thresh)
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
        if noise_note is not None:
            lines.append(noise_note)
        if sigma:
            noise_line = (f"Measured noise sigma = {sigma:.4g} "
                          f"(high-frequency corner of each spectrum), "
                          f"best-fit {value_label} = {rmsd_min:.4g}")
            if chi2_red is not None:
                noise_line += f",  reduced chi^2 = {chi2_red:.3g}"
            lines.append(noise_line)
        if chi2_red is not None and chi2_red > 2.0:
            # Worth spelling out rather than silently drawing a narrow contour:
            # at this point the residual is dominated by model error, not by
            # noise, and a region scaled to the noise alone answers the wrong
            # question. (Some excess is expected regardless -- the spectrum is
            # a magnitude, so its noise floor has a positive mean that an
            # amplitude-scaled simulation cannot reproduce.)
            lines.append(f"(residual sits {np.sqrt(chi2_red):.1f}x above the noise -- "
                         f"systematic misfit, so the measured-noise region is "
                         f"optimistically narrow; the F-test is the safer read)")
        if N <= 2:
            lines.append(f"(only N={N} pooled residual points at the best fit -- "
                         f"too few for a meaningful confidence interval)")
        if normalized:
            lines.append("Normalized: 0 = perfect fit, "
                         "1 = no better than no simulation at all.")
            if self._last_n_datasets > 1:
                # Worth stating plainly: this is not the raw statistic with a
                # different label, and the F-test/chi-square factors assume a
                # plain sum-of-squares ratio, which a dataset-reweighted mean
                # only approximates.
                lines.append("(each dataset normalized to its own null model, then "
                             "averaged -- best fit and interval can differ from raw RMSD)")
        lines.append("± is half the width of the joint 2-parameter 95% region's bounding "
                     "box on each axis (conservative vs. a true 1D profile interval).")
        self.lbl_ci.SetLabel("\n".join(lines))

        mark_resid, mark_loaded = self._marker_points()
        subtitle = f"{used_method} interp, N={N}"
        self.heatmap.show_results(fine_vals, fine1, fine2, label1, label2,
                                  contour_levels=contour_levels,
                                  background_levels=background_levels,
                                  raw_points=(grid1, grid2), subtitle=subtitle,
                                  half_width=half_width,
                                  as_contours=as_contours,
                                  rmsd_contour_width=rmsd_cw,
                                  value_label=value_label,
                                  bg_contour_color=self.map_contour_color,
                                  bg_contour_width=self.map_contour_width,
                                  mark_residual=mark_resid, mark_loaded=mark_loaded)

    def _display_probability(self, mode):
        """Draw the posterior over the scanned box, either as the credible
        level of each cell or as the density itself, with marginal credible
        intervals. Returns False if the posterior could not be formed, so
        the caller can fall back to the fitness map."""
        if self._last_results_sse is None or self._last_results_nk is None:
            return False
        grid1, grid2 = self._last_grid1, self._last_grid2
        label1 = self.param_full_label(self._last_param1)
        label2 = self.param_full_label(self._last_param2)

        divisor = float(self.neff_divisor.GetValue())
        nlp = self.neg_log_posterior(self._last_results_sse, self._last_results_nk, divisor)
        if nlp is None:
            return False

        # Interpolate the *negative log* posterior, not the posterior:
        # it is smooth and roughly quadratic near the peak, whereas the
        # probability itself spans many orders of magnitude and would
        # interpolate atrociously. It also keeps the min-is-best convention
        # that _interpolate_grid's failed-point sentinel assumes.
        method = 'cubic' if self.interp_choice.GetSelection() == 1 else 'linear'
        out = self._interpolate_grid(method, grid1, grid2, nlp)
        if out is None:
            return False
        fine1, fine2, nlp_fine, used_method = out

        p = self.posterior_grid(nlp_fine)
        if p is None:
            return False
        level = self.credible_level_map(p)
        bi, bj = np.unravel_index(np.argmax(p), p.shape)
        best1, best2 = fine1[bi], fine2[bj]

        CONF_INNER, CONF_OUTER = 0.68, 0.95
        if mode == self.DISPLAY_CREDIBLE:
            displayed = level
            contour_levels = [(CONF_INNER, "68%"), (CONF_OUTER, "95%")]
        else:
            displayed = p
            contour_levels = [(self.density_threshold(p, CONF_INNER), "68%"),
                              (self.density_threshold(p, CONF_OUTER), "95%")]
        value_label = self.DISPLAY_VALUE_LABELS[mode]

        vmin, vmax = float(np.nanmin(displayed)), float(np.nanmax(displayed))
        background_levels = self._background_levels(vmin, vmax)

        bounds1 = self.marginal_interval(fine1, p.sum(axis=1), CONF_OUTER)
        bounds2 = self.marginal_interval(fine2, p.sum(axis=0), CONF_OUTER)
        half_width = {1: None, 2: None}
        if bounds1 is not None:
            half_width[1] = (bounds1[1] - bounds1[0]) / 2.0
        if bounds2 is not None:
            half_width[2] = (bounds2[1] - bounds2[0]) / 2.0

        # Pooled point count at the best coarse point, and what the divisor
        # leaves of it -- the exponent that actually set the width above.
        if self._last_results_n is not None and np.any(np.isfinite(nlp)):
            ci, cj = np.unravel_index(np.nanargmin(nlp), nlp.shape)
            N = int(self._last_results_n[ci, cj])
        else:
            N = 0
        self._last_N = N
        n_eff = N / max(divisor, 1e-9)

        def fmt_line(label, best, bounds):
            if bounds is None:
                return f"{label} = {best:.5g}  (95% credible: n/a)"
            return (f"{label} = {best:.5g}  "
                    f"(95% credible: {bounds[0]:.5g} to {bounds[1]:.5g})")

        lines = [fmt_line(label1, best1, bounds1), fmt_line(label2, best2, bounds2)]
        if used_method != method:
            lines.append("(cubic needs >=4 points/axis and no failed grid points -- "
                         "used linear instead)")
        lines.append(f"Posterior with each dataset's noise level marginalized out; "
                     f"unaffected by dataset brightness or by the RMSD normalization.")
        lines.append(f"Effective points: {N} / {divisor:g} = {n_eff:.0f} "
                     f"(sets how sharp the posterior is).")
        if divisor <= 1.0:
            lines.append("(divisor is 1 -- if the spectra were zero-filled, these "
                         "regions are optimistically narrow)")

        # If the posterior has not decayed by the edges, normalizing over
        # the box is meaningless and the region is an artifact of where the
        # scan happened to stop.
        edge = float(p[0, :].sum() + p[-1, :].sum() + p[:, 0].sum() + p[:, -1].sum())
        if edge > 0.01:
            lines.append(f"WARNING: {edge*100:.0f}% of the posterior mass sits on the "
                         f"edge of the scanned box -- widen the ranges.")
        # The opposite failure: with a few hundred points in the exponent the
        # posterior can collapse below one cell, which renders as a single
        # dot and tells the user nothing. Say why rather than drawing a blank.
        peak_mass = float(p.max())
        if peak_mass > 0.5:
            lines.append(f"WARNING: {peak_mass*100:.0f}% of the posterior sits in one "
                         f"grid cell -- scan a narrower range around the best fit, or "
                         f"raise the divisor, to resolve it.")
        self.lbl_ci.SetLabel("\n".join(lines))

        mark_resid, mark_loaded = self._marker_points()
        subtitle = f"{used_method} interp, N={N}, eff {n_eff:.0f}"
        self.heatmap.show_results(displayed, fine1, fine2, label1, label2,
                                  contour_levels=contour_levels,
                                  background_levels=background_levels,
                                  raw_points=(grid1, grid2), subtitle=subtitle,
                                  half_width=half_width,
                                  as_contours=self.chk_rmsd_contours.GetValue(),
                                  rmsd_contour_width=float(self.rmsd_contour_width.GetValue()),
                                  value_label=value_label,
                                  best_index=(int(bi), int(bj)),
                                  map_title=f'Posterior: {value_label.lower()}',
                                  bg_contour_color=self.map_contour_color,
                                  bg_contour_width=self.map_contour_width,
                                  mark_residual=mark_resid, mark_loaded=mark_loaded)
        return True

    # ------------------------------------------------------------------
    # Heatmap interaction
    # ------------------------------------------------------------------
    def _on_hover(self, x, y, val):
        vtxt = "n/a" if np.isnan(val) else f"{val:.5g}"
        # One string, built by implicit concatenation across the two lines --
        # note there is deliberately no comma between them. SetStatusText's
        # second argument is the field *index* (an int), so a comma here
        # makes wx read the second piece of text as a field number.
        self.statusbar.SetStatusText(
            f"{self.heatmap.label1} = {x:.5g},  {self.heatmap.label2} = {y:.5g},  "
            f"{self.heatmap.value_label} = {vtxt}")

    def on_toggle_residuals(self, event):
        """Show/hide the residual-map row. While hidden it is not merely
        invisible -- double-clicking the fitness map does nothing at all,
        so the (potentially expensive) redraw never happens unless the row
        is actually on screen."""
        self._show_residuals = not self._show_residuals
        self.residuals.Show(self._show_residuals)
        self.btn_residuals.SetLabel(
            "Hide Residual Maps" if self._show_residuals else "Show Residual Maps")
        self.btn_residuals.SetColorKey(
            "pill_warn" if self._show_residuals else "pill_calm")
        if self._show_residuals and not self._sim_store:
            self.residuals.show_message(
                "No stored simulations. Run a scan with "
                "\"Keep simulations in memory\" enabled to get residual maps.")
        self.right_panel.Layout()
        self.Layout()

    def _default_oversampling(self):
        """Oversampling implied by the main window's zero-fill setting.

        update_FFT pads to 2**(ceil(log2 N) + zfill), so each axis gains a
        factor 2**zfill in points and the 2-D count gains that squared --
        4**zfill. Taken from the datasets that are shown (each carries its
        own fftmethod once it has been transformed), falling back to the
        FFT parameter grid, and using the largest if they disagree.

        Assumes the acquired size is a power of two, as asked. When it is
        not, the true factor is (2**(ceil(log2 N) + zfill) / N)**2, which is
        larger -- so this is a floor, and the credible regions it gives are
        the optimistic ones."""
        zfill = None
        try:
            for dd in self.main.Data:
                if not dd.get('show'):
                    continue
                method = dd.get('fftmethod')
                z = method.get('zfill') if isinstance(method, dict) else None
                if z is None:
                    prop = self.main.tabulated_panel.ffttree.GetProperty('zfill')
                    z = None if prop is None else prop.GetValue()
                if z is not None:
                    zfill = int(z) if zfill is None else max(zfill, int(z))
        except Exception:
            return 1.0
        if not zfill or zfill <= 0:
            return 1.0
        return float(min(4 ** zfill, 100000.0))

    def display_mode(self):
        sel = self.display_choice.GetSelection()
        return self.DISPLAY_RMSD if sel == wx.NOT_FOUND else sel

    def on_display_mode_changed(self, event):
        """Switching statistic redraws the fitness map and, if a grid
        point's residuals are on screen, re-renders those in the new units
        too. The divisor and the confidence-interval method only apply to
        some modes, so their controls follow the selection."""
        if event is not None:
            event.Skip()      # EVT_KILL_FOCUS must keep propagating
        self._sync_display_controls()
        self._refresh_display()
        if self._last_click_ij is not None:
            self._render_residuals(*self._last_click_ij)

    def _sync_display_controls(self):
        probability = self.display_mode() in self.PROBABILITY_MODES
        # The oversampling divisor is no longer the posterior's alone: the
        # measured-noise chi-square divides by it too, for the same reason
        # (neighbouring points of a zero-filled spectrum are not independent
        # measurements), so leaving it greyed out there would hide a control
        # that is setting the width of the region on screen.
        noise_ci = (self.ci_method_choice.GetSelection() == 2
                    and self.chk_use_noise.GetValue())
        self.neff_divisor.Enable(probability or noise_ci)
        self.lbl_divisor.Enable(probability or noise_ci)
        # Credible intervals come from the posterior itself in probability
        # mode, so the F-test/chi-square choice has nothing to act on.
        self.ci_method_choice.Enable(not probability)
        self.chk_use_noise.Enable(not probability)

    def _on_doubleclick(self, x, y):
        """Rebuild the residual row for whichever *coarse* grid point the
        double-click landed nearest. The heatmap shows the interpolated
        fine surface, but simulations only exist where they were actually
        run, so the click is snapped back onto the coarse axes."""
        if not self._show_residuals:
            return
        if self._last_grid1 is None or self._last_grid2 is None:
            return
        if not self._sim_store or self._roi_cache is None:
            self.residuals.show_message(
                "No stored simulations for this scan.")
            return
        if x is None or y is None:
            return
        i = int(np.argmin(np.abs(self._last_grid1 - x)))
        j = int(np.argmin(np.abs(self._last_grid2 - y)))
        self._render_residuals(i, j)
        # Move the "residuals shown" marker to the point just picked.
        self._refresh_display()

    def _render_residuals(self, i, j):
        """Draw the residual row for coarse grid point (i, j).

        These are always drawn in normalized units, whatever the fitness
        map is showing: each dataset's residual is divided by its own null
        root-mean-square (the region's experimental root-mean-square
        intensity), so a displayed residual of 1 is as large as the typical
        experimental intensity there and every map shares the fixed +-2
        scale. They are a display, not the statistic -- raw intensity units
        would differ by orders of magnitude between datasets and make the
        maps incomparable and mostly saturated.

        The simulation contours are scaled by the same factor; that leaves
        the lines exactly where they were (levels are taken from the
        simulation's own range) but keeps their values in the same units as
        the map underneath."""
        if self._last_grid1 is None or not self._sim_store or self._roi_cache is None:
            return
        slot = self._sim_store.get((i, j))
        v1, v2 = self._last_grid1[i], self._last_grid2[j]
        label1 = self.param_full_label(self._last_param1)
        label2 = self.param_full_label(self._last_param2)
        header = f"{label1} = {v1:.5g},  {label2} = {v2:.5g}"
        self._last_click_ij = (i, j)

        if slot is None:
            self.residuals.show_message(
                f"{header}\nNo successful simulation was stored at this grid point.")
            self.statusbar.SetStatusText(f"No stored simulation at {header}.")
            return

        entries = []
        for roi, sim in zip(self._roi_cache, slot):
            if roi is None or sim is None:
                continue
            # Blanked outside the tilted box: those points were never fitted,
            # so showing their residual would invite reading a miss into a
            # region the scan never scored.
            resid = np.where(roi['mask'], roi['exp'] - sim, np.nan)
            if roi['null_rms'] > 0.0:
                scale = 1.0 / roi['null_rms']
                resid = resid * scale
                sim = sim * scale
            # Head each map with the field, matching how the main window
            # titles its spectra (see MainFrame's 'B$_0$={field} mT'),
            # rather than with the file name -- which of several tau/field
            # measurements this is, is the useful identity here.
            field = roi.get('field')
            title = (f"B$_0$={field:g} mT" if isinstance(field, (int, float))
                     else str(roi['label']))
            entries.append({
                'title': title,
                'resid': resid,
                'sim': sim,
                'x': roi['x'],
                'y': roi['y'],
                'extent': [roi['x'][0], roi['x'][-1], roi['y'][0], roi['y'][-1]],
                'rmsd': float(np.sqrt(np.nanmean(resid ** 2))),
            })
        self.residuals.show_maps(
            entries, subtitle=f"Residuals at {header} (normalized)",
            value_label="Normalized RMSD")
        self.statusbar.SetStatusText(f"Residual maps at {header} (normalized).")

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
                    self.apply_param_value(main_params, p, v)
        except KeyError:
            wx.MessageBox(
                "The scanned parameter no longer exists in the main Sys grid "
                "(it may have been removed/renamed). Could not load values.",
                "Grid Search", wx.OK | wx.ICON_WARNING)
            return
        self.main.tabulated_panel.Sys_param.SetFromParClean(main_params)
        self.main.Sys.setFromCtrl(main_params)
        self.main.runSim()
        msg = (f"Loaded {self.param_full_label(self._last_param1)}={v1:.5g}, "
               f"{self.param_full_label(self._last_param2)}={v2:.5g} into main window "
               f"and re-ran simulation.")
        # The main Sys grid just changed under us -- refresh this window's
        # own parameter list so the picker's displayed values (and which
        # leaves are Function-Modifier-driven) reflect the new baseline,
        # same as clicking "Refresh Parameters from Main Window".
        self.on_refresh(None)
        self._refresh_display()
        self.statusbar.SetStatusText(msg)

    # ------------------------------------------------------------------
    def on_close(self, event):
        self._closing = True
        if self._running:
            self._cancel_event.set()
        # Drop the stored simulations now rather than waiting for the frame
        # itself to be collected -- this can be hundreds of megabytes.
        self._sim_store = {}
        self._roi_cache = None
        event.Skip()

    # ------------------------------------------------------------------
    # Parameter descriptors and the Sys dictionary
    # ------------------------------------------------------------------
    def flatten_sys_params(self, sys_dict):
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
    def qualified_label(self, p):
        return f"{p['cat']}.{self.param_leaf_label(p)}"
    def flatten_func_vars(self, variable_values):
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
    def param_name(self, p):
        kind = p.get('kind', 'sysleaf')
        return f"{kind}:{p['cat']}:{p['key']}:{-1 if p['idx'] is None else p['idx']}"
    def param_leaf_label(self, p):
        return p['key'] if p['idx'] is None else f"{p['key']}({p['idx'] + 1})"
    def param_full_label(self, p):
        return f"{p['cat']} {self.param_leaf_label(p)}"
    def read_param_value(self, sys_dict, param):
        """Current value of a scannable parameter in `sys_dict`, or None if it
        is no longer there (a nucleus removed, a Function-Modifier variable
        renamed, ...). The read counterpart of apply_param_value, but it also
        handles 'funcvar' leaves, which live under 'functions'.'var' rather
        than as an ordinary Sys leaf."""
        try:
            if param.get('kind', 'sysleaf') == 'funcvar':
                fdict = sys_dict.get('functions')
                if not isinstance(fdict, dict):
                    return None
                variables = fdict.get('var')
                if not isinstance(variables, dict) or param['key'] not in variables:
                    return None
                return float(variables[param['key']])

            sub = sys_dict.get(param['cat'])
            if not isinstance(sub, dict) or param['key'] not in sub:
                return None
            val = sub[param['key']]
            if param['idx'] is None:
                return float(val)
            arr = np.asarray(val)
            if arr.ndim != 1 or param['idx'] >= arr.size:
                return None
            return float(arr[param['idx']])
        except (TypeError, ValueError):
            return None
    def apply_param_value(self, sys_dict, param, value):
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
    def apply_functions_to_sys_dict(self, sys_dict, func_rows, variable_values, error_log=None):
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

        # qualified 'cat.leaf' label (e.g. 'nuc(1).A(1)') -> sysleaf descriptor,
        # for patching the results back in via apply_param_value below
        syspath_map = {self.qualified_label(p): p for p in self.flatten_sys_params(sys_dict)}
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
                self._log_gridpoint_error(error_log, f'function "{lbl} = {expr}"', RuntimeError(msg))
            return False

        for label, value in result['label_values'].items():
            p = syspath_map.get(label)
            if p is not None:
                self.apply_param_value(sys_dict, p, value)
        return True

    # ------------------------------------------------------------------
    # Fitness: one grid point, all datasets
    # ------------------------------------------------------------------
    def _silent_error_func(self, message):
        raise _SilentSimError(message)
    def dataset_roi(self, dd):
        """Extract one dataset's experimental spectrum restricted to its
        [fmin, fmax] x [fmin, fmax] region, together with the frequency axes
        of that region. Returns a dict with keys 'exp', 'x', 'y', 'label', or
        None when the dataset has no usable frequency-domain data or the
        region selects nothing.

        'exp' stays the full rectangle -- residual maps are drawn from it and
        want a grid -- but what is actually fitted is the 45-degree tilted box
        inscribed in that square (optHYSCORE.diagonalBox, widened or narrowed
        by the dataset's own AntiDiagSpan), and 'mask' marks it. Every sum below runs over the mask,
        so the region scored here is exactly the one the main window outlines
        on its 2D maps and projects in its diagonal skyline panel. At the
        default AntiDiagSpan of -1 that is about half the square's points, and
        less again as the span widens, so 'n' rather than exp.size is the count
        the statistics must use.

        Split out of _dataset_sse so a scan can build this once per dataset up
        front instead of re-slicing the same unchanging experimental data at
        every grid point, and so the region's size is known before the scan
        starts (needed to predict how much memory keeping every simulation
        would cost -- see the estimate in _start_run)."""
        if dd['isfft']:
            raw = np.asarray(dd['data'])
            expdata = np.real(raw[:, :, 0]) if raw.ndim > 2 else np.real(raw)
            expX = np.asarray(dd['ax']['x'])
            expY = np.asarray(dd['ax']['y'])
        else:
            if not dd['fftactual'] or dd['fftdata'] is None:
                return None
            expdata = dd['fftdata']
            expX = np.asarray(dd['fftax']['x'])
            expY = np.asarray(dd['fftax']['y'])

        fmin, fmax = dd['fmin'], dd['fmax']
        # Array axis 0 tracks the 'x' coordinate array and axis 1 tracks 'y'
        # throughout this codebase (see MainFrame.update_FFT, which builds
        # fftX from fftData.shape[0] and fftY from fftData.shape[1]) -- index
        # accordingly rather than assuming the usual imshow row=y convention,
        # and hand diagonalBox the axes in that same order so its masks and
        # indices come back matching.
        box = self.main.Opt.diagonalBox(expX, expY, fmin, fmax,
                                        dd.get('antidiagspan', -1.0))
        if box is None:
            return None
        xsel, ysel, mask = box['sel0'], box['sel1'], box['mask']
        n_fit = int(mask.sum())
        if not n_fit:
            return None

        exp_roi = expdata[np.ix_(xsel, ysel)]
        # 'null' is this dataset's no-simulation sum of squares -- the residual
        # you are left with at c = 0. Since the amplitude fit clamps c to be
        # non-negative it can never do worse than that, so SSE/null is bounded
        # in (0, 1] and is what the normalized RMSD divides by. 'null_rms' is
        # the same thing per point (the region's root-mean-square intensity),
        # used to put residual *maps* into the same normalized units.
        null = float(np.sum(exp_roi[mask] ** 2))
        # 'sigma' is the per-point noise level pyHYSCORE measured on this
        # spectrum's far high-frequency corner when the FFT was built
        # (MainFrame.noise_from_spectrum) -- an *independent* estimate, which
        # is what turns the confidence region into a real chi-square instead
        # of one calibrated on the fit's own residual. None for a session
        # saved before this existed, or when the border was unmeasurable; the
        # CI code falls back to the F-test in that case.
        sigma = dd.get('noise')
        return {'exp': exp_roi,
                'mask': mask,
                'n': n_fit,
                'corners': box['corners'],
                'x': expX[xsel],
                'y': expY[ysel],
                'null': null,
                'null_rms': float(np.sqrt(null / n_fit)),
                'sigma': float(sigma) if sigma else None,
                'label': dd.get('title', dd.get('fname', '?')),
                'field': dd.get('field')}
    def _dataset_sse(self, dd, simdata, simX, simY, roi=None, want_sim=False):
        """Amplitude-scale simdata to dd's experimental data (minimizing SSE
        over the [fmin, fmax] region) and return (sum-of-squared-residuals,
        n_points, scaled_sim) for this one dataset.

        roi: a precomputed dataset_roi(dd) result, to avoid re-slicing the
        experimental data at every grid point; computed here when omitted.

        want_sim: when True the third element is the simulation interpolated
        onto the experimental region *and already multiplied by the best-fit
        amplitude* c, i.e. exactly the array the residual was formed from, so
        a stored copy can be turned back into a residual map later by
        subtracting it from the same region's experimental data. It is None
        otherwise -- callers that only need the fitness number keep the old
        memory profile, since the scaled array is a temporary either way."""
        roi = self.dataset_roi(dd) if roi is None else roi
        if roi is None:
            return 0.0, 0, None

        exp_roi, sub_x, sub_y, mask = roi['exp'], roi['x'], roi['y'], roi['mask']

        interp = RegularGridInterpolator((simX, simY), simdata,
                                          bounds_error=False, fill_value=0.0)
        gx, gy = np.meshgrid(sub_x, sub_y, indexing='ij')
        sim_on_exp = interp(np.stack([gx.ravel(), gy.ravel()], axis=-1)).reshape(gx.shape)

        # Both the amplitude and the residual are taken over the tilted box
        # alone (roi['mask']). The simulation is still interpolated onto the
        # whole rectangle, since `scaled` is what a stored residual map is
        # later rebuilt from and that map is drawn as a grid.
        denom = np.sum(sim_on_exp[mask] ** 2)
        c = max(np.sum(exp_roi[mask] * sim_on_exp[mask]) / denom, 0.0) if denom > 1e-30 else 0.0
        scaled = c * sim_on_exp
        resid = exp_roi[mask] - scaled[mask]
        return float(np.sum(resid ** 2)), int(resid.size), (scaled if want_sim else None)
    def format_bytes(self, nbytes):
        for unit in ('bytes', 'kilobytes', 'megabytes', 'gigabytes'):
            if nbytes < 1024 or unit == 'gigabytes':
                return f"{nbytes:.0f} {unit}" if unit == 'bytes' else f"{nbytes:.2f} {unit}"
            nbytes /= 1024.0
    def _log_gridpoint_error(self, error_log, where, exc):
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
    def run_gridpoint(self, hs, Sys, exp_template, opt_template, active_data, error_log=None,
                       roi_cache=None, sim_store=None, stats_out=None):
        """Run one HYSCORE simulation per active dataset for the given Sys and
        return (pooled RMSD, pooled N, normalized RMSD) over all datasets
        (NaN, 0, NaN if none succeeded). N (the total number of pooled ROI
        residual points) is returned alongside the RMSD so confidence-interval
        estimation can later convert RMSD back to a sum-of-squares without
        re-running anything.

        The normalized RMSD scores each dataset against its own null model
        (SSE_k / null_k, i.e. how much of that dataset's power the simulation
        failed to explain -- see dataset_roi) and then takes the N-weighted
        mean of those fractions before the square root. Normalizing per
        dataset rather than pooling the sums stops one bright dataset from
        dominating a multi-dataset fit, at the cost of being a genuinely
        different statistic: with more than one dataset its landscape (and so
        its best-fit point) can differ from the raw pooled RMSD's. With a
        single dataset it is exactly the raw RMSD divided by a constant.

        Failures are not fatal to the scan (a single bad grid point/dataset
        combination, e.g. "no resonances", shouldn't abort hundreds of other
        points) but are logged via _log_gridpoint_error so they are never
        silently invisible.

        roi_cache: optional list, parallel to active_data, of precomputed
        dataset_roi() results (None entries for unusable datasets).

        sim_store: optional list, also parallel to active_data, into which
        each dataset's amplitude-scaled interpolated simulation is written by
        index -- writing by index rather than appending keeps the slots
        aligned with active_data even when some datasets fail and are skipped.
        Left as None (i.e. not requested) this costs nothing.

        stats_out: optional list, parallel to active_data, filled with each
        dataset's (SSE, N) pair. Two floats per dataset per grid point, which
        is nothing next to the simulations themselves, and it is what the
        posterior needs -- the pooled RMSD alone cannot reconstruct the
        per-dataset terms that noise-level marginalization requires."""
        hs.Sys = Sys
        if roi_cache is None:
            roi_cache = [self.dataset_roi(dd) for dd in active_data]
        try:
            hs.preCompute()
        except Exception as e:
            self._log_gridpoint_error(error_log, "preCompute", e)
            return np.nan, 0, np.nan

        total_sse = 0.0
        total_n = 0
        norm_acc = 0.0    # sum of N_k * (SSE_k / null_k)
        norm_n = 0        # sum of N_k over datasets that could be normalized
        for k, dd in enumerate(active_data):
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
                self._log_gridpoint_error(error_log, f"reRun ({dd.get('title', dd.get('fname', '?'))})", e)
                continue
            if hs.Spectrum is None:
                continue

            roi = roi_cache[k]
            sse, n, scaled = self._dataset_sse(dd, hs.Spectrum, hs.X, hs.Y, roi=roi,
                                           want_sim=sim_store is not None)
            if n > 0:
                total_sse += sse
                total_n += n
                if sim_store is not None:
                    sim_store[k] = scaled
                if stats_out is not None:
                    stats_out[k] = (sse, n)
                # A region that is identically zero has no power to explain and
                # cannot be normalized; it still counts towards the raw RMSD.
                if roi is not None and roi['null'] > 0.0:
                    norm_acc += n * (sse / roi['null'])
                    norm_n += n

        if total_n == 0:
            return np.nan, 0, np.nan
        rmsd_norm = float(np.sqrt(norm_acc / norm_n)) if norm_n > 0 else np.nan
        return float(np.sqrt(total_sse / total_n)), total_n, rmsd_norm

    # ------------------------------------------------------------------
    # Interpolating the scanned grid, and confidence regions
    # ------------------------------------------------------------------
    def _interpolate_grid(self, method, grid1, grid2, results, n_out=160):
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
    def pooled_sigma(self, normalized):
        """The measured noise level in the units of whichever fitness
        statistic is on display, as an N-weighted pooled variance over the
        datasets that carry an estimate:

            raw:        sigma^2     = sum_k N_k sigma_k^2 / sum_k N_k
            normalized: sigma_norm^2 = sum_k N_k (sigma_k/null_rms_k)^2 / sum_k N_k

        The second form is the exact counterpart of what run_gridpoint does to
        the residuals themselves -- the normalized statistic divides each
        dataset's mean square residual by its own mean square intensity
        (null_rms_k^2), so the noise has to be divided by the same thing or
        the threshold and the surface would be in different units.

        Returns None unless *every* contributing dataset has a noise estimate:
        pooling a measured level with a missing one would quietly report a
        confidence region calibrated on only part of the data.
        """
        if not self._last_noise:
            return None
        acc, n_tot = 0.0, 0
        for sigma, n_k, null_rms in self._last_noise:
            if not sigma or n_k <= 0:
                return None
            if normalized:
                if not null_rms:
                    return None
                sigma = sigma / null_rms
            acc += n_k * sigma ** 2
            n_tot += n_k
        return float(np.sqrt(acc / n_tot)) if n_tot else None
    def _rmsd_threshold_factor(self, method, conf, N, p=2, sigma=None, sse_min=None,
                                divisor=1.0):
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

        'noise': the same Delta-chi^2, but with a noise level measured
        independently on the spectrum's own far high-frequency corner rather
        than borrowed from the fit -- which is what a chi-square is supposed
        to be.
        The region is chi^2(theta) <= chi^2_min + chi2.ppf(conf, df=p) with
        chi^2 = SSE/(divisor*sigma^2), i.e.

            SSE(theta) <= SSE_min + divisor * Delta-chi^2 * sigma^2

        so the width no longer depends on how good the fit happens to be, only
        on how noisy the data is. `divisor` is the same oversampling correction
        the posterior uses: adjacent points of a zero-filled, apodized spectrum
        are not independent measurements, and counting them as if they were
        would shrink the region by exactly that factor. Needs `sigma` and
        `sse_min`; returns None without them, so the caller can fall back.
        """
        if N is None or N <= p:
            return None
        if method == 'noise':
            if not sigma or not sse_min or sse_min <= 0:
                return None
            chi2val = stats.chi2.ppf(conf, df=p)
            factor = 1.0 + max(float(divisor), 1e-9) * chi2val * sigma ** 2 / sse_min
        elif method == 'fstat':
            Fval = stats.f.ppf(conf, p, N - p)
            factor = 1.0 + (p / (N - p)) * Fval
        else:
            chi2val = stats.chi2.ppf(conf, df=p)
            factor = 1.0 + chi2val / N
        return float(np.sqrt(factor))
    def _axis_confidence_bounds(self, axis_vals, profile, thresh):
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

    # ------------------------------------------------------------------
    # Posterior probability over the scanned grid
    # ------------------------------------------------------------------
    def neg_log_posterior(self, results_sse, results_nk, divisor=1.0):
        """Negative log posterior (up to an additive constant) on the coarse
        grid, from the per-dataset sums of squared residuals.

        Assuming Gaussian residuals with a noise level that is *unknown and
        separate for each dataset*, and integrating each of those unknown
        levels out under the usual scale-invariant prior, the posterior over
        the scanned parameters comes out as a product of per-dataset terms,

            p(theta) ~ prod_k SSE_k(theta) ** (-N_k/2)

        so the quantity to minimize is sum_k (N_k/2)*log SSE_k. Two things
        follow that are worth knowing. First, it is automatically invariant to
        each dataset's overall brightness: rescaling a dataset changes
        log SSE_k by a constant, which shifts the whole surface without moving
        anything -- so unlike the fitness map, this does not care whether the
        RMSD was normalized. Second, it is the *logarithm* of each dataset's
        residual that gets averaged, not the residual itself, which is what
        makes the noise-level marginalization show up as a plain reweighting.

        divisor: how many raw points make up one independent one. Adjacent
        points in a zero-filled/apodized spectrum are not independent, and N
        sits in an exponent here, so overstating it drives the posterior
        towards a delta function. Dividing N_k by this before use is the
        crude but visible correction.

        Returns a (n1, n2) array, offset so its minimum is 0, with NaN where
        no dataset produced a usable residual -- or None if that is everywhere.
        """
        sse = np.asarray(results_sse, dtype=float)
        nk = np.asarray(results_nk, dtype=float)
        div = max(float(divisor), 1e-9)

        usable = np.isfinite(sse) & (nk > 0)
        # A dataset that fits perfectly would send log SSE to -inf; floor it
        # rather than dropping the dataset, which would silently change the
        # weighting at that one grid point.
        safe = np.where(usable & (sse > 0), sse, 1.0)
        safe = np.maximum(safe, 1e-300)

        contrib = np.where(usable, (nk / (2.0 * div)) * np.log(safe), 0.0)
        any_ok = usable.any(axis=2)
        if not any_ok.any():
            return None
        nlp = np.where(any_ok, contrib.sum(axis=2), np.nan)
        return nlp - np.nanmin(nlp)
    def posterior_grid(self, nlp):
        """exp(-nlp) normalized to sum to 1 over the grid. Failed points get
        exactly zero probability. Returns None if nothing is finite.

        Note this makes the scanned box the prior: the result is only
        meaningful if the posterior has decayed to negligible before the edges,
        which is worth checking by eye on the map itself."""
        nlp = np.asarray(nlp, dtype=float)
        good = np.isfinite(nlp)
        if not good.any():
            return None
        p = np.zeros(nlp.shape, dtype=float)
        p[good] = np.exp(-(nlp[good] - nlp[good].min()))
        total = p.sum()
        if not np.isfinite(total) or total <= 0:
            return None
        return p / total
    def _sorted_levels(self, p):
        """(densities sorted high to low, credible level of each). Cells of
        equal density share one level -- the cumulative mass through the end of
        their tied block. Without that, an arbitrary sort order splits a ring of
        equal-probability cells across the contour, and the credible-level and
        density displays disagree about where the boundary is."""
        flat = np.asarray(p, dtype=float).ravel()
        order = np.argsort(flat)[::-1]
        s = flat[order]
        cum = np.cumsum(s)
        # index of the last cell in each run of equal density
        last = np.searchsorted(-s, -s, side='right') - 1
        return order, s, cum[last]
    def credible_level_map(self, p):
        """Label every cell with the smallest credible region that contains
        it: near 0 at the most probable cell, 1 at the least. Contouring the
        result at 0.68 / 0.95 draws the credible regions, and cells of zero
        probability (failed points) land at 1."""
        order, _s, lvl_sorted = self._sorted_levels(p)
        level = np.empty(order.size, dtype=float)
        level[order] = lvl_sorted
        return level.reshape(np.shape(p))
    def density_threshold(self, p, conf):
        """The density value whose credible level is `conf` -- i.e. the contour
        level that draws the same boundary on a plot of the density as
        contouring the credible-level map at `conf` does. Interpolated between
        adjacent cells so the two displays agree to better than one cell."""
        _order, s, lvl_sorted = self._sorted_levels(p)
        # lvl_sorted is non-decreasing, s non-increasing, so this reads off the
        # density at which the credible level passes conf.
        return float(np.interp(conf, lvl_sorted, s))
    def marginal_interval(self, axis_vals, marg, conf):
        """Highest-density interval of a 1-D marginal, returned as the span of
        the selected cells. As with the existing confidence bounds, a
        multimodal marginal is reported as the range that brackets every
        included mode rather than as disjoint pieces -- conservative, and
        consistent with what the F-test path reports."""
        m = np.asarray(marg, dtype=float)
        total = m.sum()
        if not np.isfinite(total) or total <= 0:
            return None
        m = m / total
        order = np.argsort(m)[::-1]
        cum = np.cumsum(m[order])
        k = min(int(np.searchsorted(cum, conf)) + 1, m.size)
        sel = np.sort(order[:k])
        return float(axis_vals[sel[0]]), float(axis_vals[sel[-1]])
