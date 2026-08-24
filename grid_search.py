"""
2D grid-search / fitness-mapping window for pyHYSCORE.

"""
import copy
import os
import time
import threading
import multiprocessing as mp

import numpy as np
import wx
import wx.propgrid as wxpg

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
import matplotlib.cm as cm

from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
from scipy import stats
from scipy.special import ndtri

from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim
from colormap_editor import ColormapEditorDialog
import sys_functions as sysfun
import grid_search_worker as gsw
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
                     value_label="RMSD", best_value=None, map_title=None,
                     bg_contour_color='white', bg_contour_width=0.5,
                     grid_lines=True, mark_residual=None, mark_loaded=None):
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

        #extent = [grid1[0], grid1[-1], grid2[0], grid2[-1]]
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
            dx = (grid1[-1] - grid1[0]) / (len(grid1) - 1) if len(grid1) > 1 else 1.0
            dy = (grid2[-1] - grid2[0]) / (len(grid2) - 1) if len(grid2) > 1 else 1.0
            extent = [grid1[0] - dx/2, grid1[-1] + dx/2,
                      grid2[0] - dy/2, grid2[-1] + dy/2]
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

        if best_value is not None:
            self.ax.plot(best_value[0], best_value[1], marker='*', markersize=16,
                     markeredgecolor='k', markerfacecolor='white',
                     linestyle='none', label='best fit', zorder=5)

            # text = (f"Best: ")
            # best_ci1 = (f"")
            # best_ci2 = (f"")
            # if half_width is not None:
            #     hw1 = half_width[0]
            #     hw2 = half_width[1]

            #     if hw1 is not None: 
            #         best_ci1 += f"{best_value[0]:.4g} ± {hw1:.4g} (95%)\n"
            #     if hw2 is not None:
            #         best_ci2 += f"{best_value[1]:4g} ± {hw2:.4g} (95%)\n"
            # text += best_ci1
            # text += best_ci2
            # bi = np.argmin( np.abs(grid1-best_value[0]))
            # bj = np.argmin( np.abs(grid2-best_value[1]))
            # self.ax.annotate(
            #     f"{text}{value_label}={results[bi, bj]:.4g}",xy=(grid1[bi], grid2[bj]),
            #     textcoords='offset points',
            #     xytext=(8, 8), fontsize='small', color='white',
            #     bbox=dict(boxstyle='round', fc='black', alpha=0.6))

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
            # Attaching a colorbar anchors its parent east
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
    # DISPLAY_RMSD = 0
    # DISPLAY_NORMALIZED = 1
    # DISPLAY_CREDIBLE = 2
    # DISPLAY_DENSITY = 3
    DISPLAY_LABELS = ["RMSD", "Fit Quality, norm(RMSD)",
                      "Probability", "log(probability)"]
    # DISPLAY_VALUE_LABELS = ["RMSD", "Fit Quality",
    #                         "Credible level", "Posterior density"]
    #PROBABILITY_MODES = (DISPLAY_CREDIBLE, DISPLAY_DENSITY)

    # Every stored simulation is one 64-bit float per experimental point in
    # the region of interest (the interpolator output is real --
    # HYSCOREsim.Spectrum is an np.abs magnitude spectrum).
    # Flip to False to run a scan synchronously on the GUI thread, which is
    # the only way to step through a grid point in a debugger that cannot
    # follow a background thread. Untick "parallel" as well.
    RUN_IN_THREAD = True

    STORE_BYTES_PER_POINT = 4
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
        self._last_rmsd = None
        self._last_rmsd_norm = None
        self._last_rmsd_n = None
        self._last_rmsd_SumSqRes = None
        self._last_rmsd_nk = None
        self._last_n_datasets = 0
        self._last_noise = None
        self._fit_noise = False
        self._last_grid1 = None
        self._last_grid2 = None
        self._last_N = 0
        self._last_click_ij = None
        self._sim_store = {}
        # Every store this panel has created, including cancelled scans' and
        # superseded ones. Nothing is erased until the window closes, and this
        # is the list of what to erase then.
        self._sim_stores = []
        self._pool = None          # live worker pool, so Cancel can kill it
        self._region_cache = None
        self._show_residuals = False

        self.map_contour_on = True
        self.map_contour_color = '#ffffff'
        self.map_contour_n = 9
        self.map_contour_width = 0.5

        self._range_memory = {}

        self.cmap_name = self.main.current_cmap_name
        self.cmap_stops = list(self.main.current_cmap_stops)
        self.cmap_base = self.main.current_cmap
        self.cmap_inverted = False
        self.cmap = self.cmap_base

        self.statusbar = theme.StatusStrip(self)
        self.statusbar.SetStatusText("Select two Sys parameters, set their ranges, then Run.")

        splitter = wx.SplitterWindow(self, style=wx.SP_3D)
        left = wx.Panel(splitter)

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
        splitter.SetSashGravity(0)
        splitter.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self._on_sash_moved)
        splitter.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self._on_sash_moved)
        self.splitter = splitter

        leftsizer = wx.BoxSizer(wx.VERTICAL)
        
        #----------------------------------------------------------------------
        
        self.btn_refresh = PillButton(left, "Refresh Parameters from Main Window",
                                      color_key="pill_load", 
                                      draw_badge=theme.draw_textlines_badge)
        leftsizer.Add(self.btn_refresh, 0, wx.EXPAND | wx.ALL, 4)
        self.param_pg = wxpg.PropertyGrid(
            left, style=wxpg.PG_DEFAULT_STYLE | wxpg.PG_HIDE_MARGIN | wxpg.PG_TOOLTIPS)
        leftsizer.Add(self.param_pg, 1, wx.EXPAND | wx.ALL, 4)
        wx.CallAfter(self.param_pg.SetSplitterPosition, 100)

        leftsizer.Add(wx.StaticLine(left), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 4)
        
        #----------------------------------------------------------------------
        self.range_panel = wx.Panel(left)
        self.range_sizer = wx.FlexGridSizer(cols=4, vgap=4, hgap=6)
        self.range_sizer.AddGrowableCol(1)
        self.range_sizer.AddGrowableCol(2)
        self.range_panel.SetSizer(self.range_sizer)
        leftsizer.Add(self.range_panel, 0, wx.EXPAND | wx.ALL, 4)
        
        self.chk_use_noise = wx.CheckBox(left, label="subtract exp. noise")
        self.chk_use_noise.SetValue(False)
        self.chk_use_noise.SetToolTip(
            "subtract each spectrum's per-point noise level.")
        self.chk_store_sims = wx.CheckBox(
            left, label="Keep simulations in memory")
        self.chk_store_sims.SetValue(True)
        self.chk_store_sims.SetToolTip(
            "Keep every grid point's simulation to display residuals. \n Untick if disk space is limiting.")

        leftsizer.Add(self.chk_store_sims, 0, wx.ALIGN_LEFT | wx.RIGHT, 10)

        leftsizer.Add(self.chk_use_noise, 0, wx.ALIGN_LEFT)

        # ---- Parallel execution ------------------------------------------
        # Grid points are independent, so they can be farmed out to separate
        # processes -- separate, not threads, because the simulation spends
        # its time in small numpy calls with the GIL held and would not
        # overlap at all. One core is left to the GUI by default so the
        # window keeps repainting while a scan runs.
        n_cpu = max(1, os.cpu_count() or 1)
        parsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_parallel = wx.CheckBox(left, label="parallel, workers:")
        self.chk_parallel.SetValue(n_cpu > 1)
        self.chk_parallel.SetToolTip(
            "Run grid points in separate processes.\n"
            "Untick to run the scan in a single background thread, which is "
            "what you want when stepping through a grid point in a debugger.")
        self.spin_workers = wx.SpinCtrl(left, min=1, max=n_cpu,
                                        initial=max(1, n_cpu - 1),
                                        size=(60, -1))
        self.spin_workers.SetToolTip(
            f"How many worker processes to run. This machine reports {n_cpu} "
            f"cores; leaving one free keeps the interface responsive.")
        self.spin_workers.Enable(self.chk_parallel.GetValue())
        self.chk_parallel.Bind(
            wx.EVT_CHECKBOX,
            lambda e: self.spin_workers.Enable(self.chk_parallel.GetValue()))
        parsizer.Add(self.chk_parallel, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        parsizer.Add(self.spin_workers, 0, wx.ALIGN_CENTER_VERTICAL)
        leftsizer.Add(parsizer, 0, wx.ALIGN_LEFT | wx.TOP, 4)
        self.btn_run = PillButton(left, "Run Grid", color_key="pill_load",
                                  draw_badge=theme.draw_play_badge)
        leftsizer.Add(self.btn_run, 0, wx.EXPAND | wx.ALL, 6)
        self.gauge = wx.Gauge(left, range=100)
        leftsizer.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 4)
        self.lbl_progress = wx.StaticText(left, label="")
        leftsizer.Add(self.lbl_progress, 0, wx.EXPAND | wx.ALL, 4)
        leftsizer.Add(wx.StaticLine(left), 0, wx.EXPAND | wx.TOP | wx.BOTTOM, 4)

        #----------------------------------------------------------------------
        
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

        
        
        # analysis_sizer.Add(wx.StaticText(left, label="CI :"),
        #                    0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        # self.ci_method_choice = wx.Choice(
        #     left, choices=["F-test", "χ²", "χ²/w noise", "skew_gau fit"])
        # self.ci_method_choice.SetSelection(0)
        # self.ci_method_choice.SetToolTip(
        #     "F-test: standard when the noise level isn't independently known\n"
        #     "χ²: Delta-χ² using the fit's own RMSD as the noise estimate\n"
        #     "χ² /w noise: the same Delta-χ² , but with the noise level for each spectrum")
        # analysis_sizer.Add(self.ci_method_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        #leftsizer.Add(analysis_sizer, 0, wx.ALL, 4)

        #####################
        display_sizer = wx.BoxSizer(wx.HORIZONTAL)
        display_sizer.Add(wx.StaticText(left, label="Display:"),
                          0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.display_choice = wx.Choice(left, choices=self.DISPLAY_LABELS)
        self.display_choice.SetSelection(2)
        self.display_choice.SetToolTip(
            "RMSD: the raw pooled fitness, in experimental intensity units.\n"
            "Normalized RMSD: each dataset scored against its own null model (0-perferct fit, 1 - no sim)\n"
            "Probability: p(param|data)~[Sum Square Roots]^(-N/2) ")
        display_sizer.Add(self.display_choice, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 10)

        display_sizer.Add(wx.StaticText(left, label="Smooth: "),
                           0, wx.ALIGN_LEFT | wx.RIGHT, 4)
        self.interp_choice = wx.Choice(left, choices=["None", "Linear", "Cubic"])
        self.interp_choice.SetSelection(2)
        self.interp_choice.SetToolTip(
            "Linear: bilinear, always stable, blocky/faceted.\n"
            "Cubic: smooth bicubic spline")
        display_sizer.Add(self.interp_choice, 0, wx.ALIGN_LEFT | wx.RIGHT, 12)
        leftsizer.Add(display_sizer, 0, wx.ALL, 4)
        oversample_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.lbl_divisor = wx.StaticText(left, label="Oversampling:")
        oversample_sizer.Add(self.lbl_divisor, 0, wx.ALIGN_LEFT | wx.RIGHT, 4)
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
        oversample_sizer.Add(self.neff_divisor, 0, wx.ALIGN_LEFT)
        leftsizer.Add(oversample_sizer, 0, wx.ALL, 4)
        #############################################
        rmsd_contour_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_rmsd_contours = wx.CheckBox(left, label="As Contours")
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

        # note = wx.StaticText(left, label=(
        #     "Note: every grid point re-runs a full simulation for each active\n"
        #     "(Show=True) dataset. Large grids x many datasets can take a while."))
        # note.SetForegroundColour(wx.Colour(theme.theme_colors()["text_dim"]))
        # leftsizer.Add(note, 0, wx.ALL, 4)



        left.SetSizer(leftsizer)
        self.left_panel = left

        framesizer = wx.BoxSizer(wx.VERTICAL)
        framesizer.Add(splitter, 1, wx.EXPAND)
        framesizer.Add(self.statusbar, 0, wx.EXPAND)
        self.SetSizer(framesizer)

        self.btn_refresh.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.btn_cmap.Bind(wx.EVT_BUTTON, self.on_choose_cmap)
        self.btn_map_contours.Bind(wx.EVT_BUTTON, self.on_edit_map_contours)
        self.chk_invert_cmap.Bind(wx.EVT_CHECKBOX, self.on_invert_cmap_toggle)
        self.interp_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        #self.ci_method_choice.Bind(wx.EVT_CHOICE, self.on_analysis_choice_changed)
        self.chk_use_noise.Bind(wx.EVT_CHECKBOX, self.on_analysis_choice_changed)
        self.btn_residuals.Bind(wx.EVT_BUTTON, self.on_toggle_residuals)
        self.chk_rmsd_contours.Bind(wx.EVT_CHECKBOX, self.on_analysis_choice_changed)
        self.display_choice.Bind(wx.EVT_CHOICE, self.on_display_mode_changed)
        self.neff_divisor.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_display_mode_changed)
        self.neff_divisor.Bind(wx.EVT_KILL_FOCUS, self.on_display_mode_changed)
        self.rmsd_contour_width.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_analysis_choice_changed)
        # IMPORTANT:
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
        #self._sync_display_controls()
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
        if self._last_rmsd is not None and self.chk_rmsd_contours.GetValue():
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

    def _relayout_left(self):
        self.range_panel.Layout()
        self.left_panel.Layout()
        self.Layout()

    def _rebuild_range_controls(self):
        for pname, ctrls in self.range_ctrls.items():
            self._range_memory[pname] = (
                ctrls['min'].GetValue(), ctrls['max'].GetValue(), ctrls['points'].GetValue())
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
                # TBD: recognize physocal range and implement that (probably from sys class)
                halfrange = max(abs(v) * 0.5, 0.5)
                min_str, max_str, npts = f"{v - halfrange:.6g}", f"{v + halfrange:.6g}", 4
                self._range_memory[pname] = (min_str, max_str, npts)
            self.range_sizer.Add(
                wx.StaticText(self.range_panel, label=self.param_full_label(p)),
                0, wx.ALIGN_CENTER_VERTICAL)
            min_ctrl = wx.TextCtrl(self.range_panel, value=min_str)
            max_ctrl = wx.TextCtrl(self.range_panel, value=max_str)
            pts_ctrl = wx.SpinCtrl(self.range_panel, min=4, max=500, initial=npts)
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
            # Cancel means stop now and keep nothing: the flag stops the
            # collector, and terminating the pool kills whatever the workers
            # are in the middle of rather than waiting out the current grid
            # point, which on a slow Sys can be tens of seconds.
            self._cancel_event.set()
            pool = self._pool
            if pool is not None:
                try:
                    pool.terminate()
                except Exception:
                    pass
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
            wx.MessageBox(str(e), "Grid Search->_start_run", wx.OK | wx.ICON_ERROR)
            return

        # # Make sure every active dataset's FFT is up to date before comparing.
        # try:
        #     self.main.update_FFT()
        # except Exception as e:
        #     wx.MessageBox(f"Could not update FFT data before scanning:\n{e}",
        #                   "Grid Search->_start_run", wx.OK | wx.ICON_ERROR)
        #     return

        param1 = self._param_by_name[self._selected[0]]
        param2 = self._param_by_name[self._selected[1]]
        grid1 = np.linspace(ranges[0][0], ranges[0][1], ranges[0][2])
        grid2 = np.linspace(ranges[1][0], ranges[1][1], ranges[1][2])
        self._last_param1, self._last_param2 = param1, param2

        # Slice each dataset's experimental region once, up front
        self._fit_noise = self.chk_use_noise.GetValue()
        
        self._region_cache = []
        self._last_noise = []
        per_point = 0
        for dd in active_data:
            region = self.main.dataset_region(dd, self._fit_noise)
            if region is not None:
                self._last_noise.append((region['sigma'], region['n'], region['null_rms']))
                per_point += region['exp'].size
            else:
                wx.MessageBox("_start_run: active region is empty \n"
                              "at least for one dataset. Abort run", "Grid Search->_start_run", wx.OK | wx.ICON_ERROR)
                return    
            self._region_cache.append(region)

        store_sims = self.chk_store_sims.GetValue()
        if store_sims:
            # predicted cost of keeping every grid point's interpolated simulation:
            est = per_point * len(grid1) * len(grid2) * self.STORE_BYTES_PER_POINT
            if est > self.STORE_WARN_BYTES:
                answer = wx.MessageBox(
                    f"Keeping every simulation  = {self.format_bytes(est)} on disk\n"
                    f"Yes = store anyway (residual maps available)\n"
                    "Grid Search - large storage requirement",
                    wx.YES_NO | wx.CANCEL | wx.ICON_WARNING)
                if answer == wx.CANCEL:
                    return
                store_sims = (answer == wx.YES)

        self._last_click_ij = None

        self.residuals.show_message(
            "Scan running -- double-click a grid point when it finishes."
            if store_sims else
            "Simulations are not being stored for this scan.")

        sys_template = copy.deepcopy(self.main.tabulated_panel.Sys_param.parameters)
        exp_template = copy.deepcopy(self.main.tabulated_panel.Exp_param.parameters)
        opt_template = copy.deepcopy(self.main.tabulated_panel.Opt_param.parameters)
        func_rows = list(self.func_rows_snapshot)          # [[label, expr], ...]
        base_variable_values = dict(self.func_var_snapshot)  # baseline for un-scanned funcvars

        # Grid points are handed out as (i, j, v1, v2); 
        sim_store = None
        if store_sims:
            try:
                sim_store = gsw.SimStore(
                    [region['exp'].shape for region in self._region_cache],
                    len(grid1), len(grid2))
            except (MemoryError, OSError, ValueError) as e:
                wx.MessageBox(
                    f"Could not reserve {self.format_bytes(est)} of shared "
                    f"memory for the stored simulations:\n\n"
                    f"{type(e).__name__}: {e}\n\n"
                    f"The scan will run without storing them.",
                    "Grid Search - cannot store simulations",
                    wx.OK | wx.ICON_WARNING)
                store_sims = False
        if sim_store is not None:
            self._sim_stores.append(sim_store)

        payload = {
            'sys_template': sys_template,
            'exp_template': exp_template,
            'opt_template': opt_template,
            'func_rows': func_rows,
            'base_variable_values': base_variable_values,
            'param1': param1,
            'param2': param2,
            'records': gsw.dataset_records(active_data, self._region_cache),
            'store_layout': sim_store.layout if sim_store is not None else None,
            }
        tasks = [(i, j, float(v1), float(v2))
                 for i, v1 in enumerate(grid1)
                 for j, v2 in enumerate(grid2)]

        total = len(tasks)
        n_workers = self.spin_workers.GetValue() if self.chk_parallel.GetValue() else 1
        use_pool = n_workers > 1

        self._running = True
        self._cancel_event = threading.Event()
        self.btn_run.SetLabel("Cancel")
        self.btn_run.SetColorKey("pill_warn")
        self.param_pg.Enable(False)
        self.btn_refresh.Disable()
        self.range_panel.Enable(False)
        self.btn_cmap.Disable()
        self.chk_store_sims.Enable(False)
        self.chk_parallel.Enable(False)
        self.spin_workers.Enable(False)
        self.gauge.SetRange(total)
        self.gauge.SetValue(0)
        self.lbl_progress.SetLabel(f"0 / {total}")
        if use_pool:
            self.statusbar.SetStatusText(
            f"Starting grid search on {n_workers} workers...") 


        args = (payload, tasks, grid1, grid2, param1, param2,
                self._cancel_event, sim_store, use_pool, n_workers)
        self.RUN_IN_THREAD = False
        if self.RUN_IN_THREAD:
            self._worker_thread = threading.Thread(target=self.worker_run,
                                                   args=args, daemon=True)
            self._worker_thread.start()
        else:
            self.worker_run(*args)

    def worker_run(self, payload, tasks, grid1, grid2, param1, param2,
                   cancel_event, sim_store, use_pool, n_workers):
        """Drive one scan to completion from a background thread.
        All modes run gsw.run_point
        """
        n1, n2 = len(grid1), len(grid2)
        n_ds = len(payload['records'])
        rmsd = np.full((n1, n2), np.nan)
        rmsd_norm = np.full((n1, n2), np.nan)
        total_n = np.zeros((n1, n2), dtype=int)
        # Per-dataset residuals, needed by the posterior. Two floats per
        # dataset per grid point -- negligible next to anything else here.
        rmsd_SumSqRes = np.full((n1, n2, n_ds), np.nan)
        rmsd_nk = np.zeros((n1, n2, n_ds), dtype=int)
        start = time.perf_counter()
        count = 0
        total = len(tasks)
        n_failed_points = 0
        error_log = {}   # unique error message -> occurrence count

        def collect(result):
            nonlocal count, n_failed_points
            (i, j, t_rmsd, t_total_n, t_rmsd_norm,
             stats, present, errs) = result
            if errs:
                self._merge_error_log(error_log, errs)
            if sim_store is not None:
                sim_store.note(i, j, present)
                
                
            if np.isnan(t_rmsd):
                n_failed_points += 1
            rmsd[i, j] = t_rmsd
            rmsd_norm[i, j] = t_rmsd_norm
            total_n[i, j] = t_total_n
            for k, st in enumerate(stats):
                if st is not None:
                    rmsd_SumSqRes[i, j, k], rmsd_nk[i, j, k] = st
            count += 1
            wx.CallAfter(self._on_progress, count, total,
                         time.perf_counter() - start)

        if use_pool:
            # 'spawn' explicitly: it is the only start method on Windows, and
            # asking for it everywhere keeps the payload requirements (and so
            # the failure modes) identical on every platform.
            ctx = mp.get_context('spawn')
            pool = ctx.Pool(processes=n_workers,
                            initializer=gsw.init_worker, initargs=(payload,))
            self._pool = pool
            self.statusbar.SetStatusText(
                f"Running grid search on {n_workers} workers...") 
            try:
                # chunksize=1: grid points differ wildly in cost (a point
                # with no resonances returns almost immediately), so handing
                # them out one at a time is what keeps the workers level.
                for result in pool.imap_unordered(gsw.run_point, tasks, chunksize=1):
                    if cancel_event.is_set():
                        break
                    collect(result)
            except Exception as e:
                # A hard cancel terminates the workers
                if not cancel_event.is_set():
                    self._log_gridpoint_error(error_log, "worker pool", e)
            finally:
                self._pool = None
                try:
                    pool.terminate()
                    pool.join()
                except Exception:
                    pass
        else:
            gsw.init_worker(payload)
            try:
                for task in tasks:
                    if cancel_event.is_set():
                        break
                    collect(gsw.run_point(task))
            finally:
                gsw.release_worker()

        if sim_store is not None:
            print(f"[GridSearch] has {len(sim_store)}/{total} grid points in memory")
        else:
            print("[GridSearch] simulations were not stored (store_sims off)")

        wx.CallAfter(self._on_finished, rmsd, total_n, grid1, grid2, param1, param2,
                     cancel_event.is_set(), n_failed_points, total, rmsd_norm,
                     error_log, sim_store,
                     n_ds, rmsd_SumSqRes, rmsd_nk)

    def _on_progress(self, count, total, elapsed):
        if self._closing:
            return
        self.gauge.SetValue(count)
        rate = elapsed / count if count else 0.0
        remaining = rate * (total - count)
        self.lbl_progress.SetLabel(
            f"{count} / {total}   elapsed {elapsed:.0f}s   ETA {remaining:.0f}s")

    def _on_finished(self, rmsd, total_n, grid1, grid2, param1, param2, cancelled,
                      n_failed, total, rmsd_norm, error_log=None, sim_store=None,
                      n_datasets=0, rmsd_SumSqRes=None, rmsd_nk=None):
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
        self.chk_parallel.Enable(True)
        self.spin_workers.Enable(self.chk_parallel.GetValue())
        self._sync_param_lock_state()

        if cancelled:
            # A cancelled scan is discarded whole. Half of a fitness map is
            # not a fitness map -- the points that did finish are scattered
            # across it in completion order, not filled in from one corner --
            # so the previous scan's results are left on display untouched.
            # Whatever the workers managed to write stays on disk until the
            # window closes, like every other scan's.
            self.lbl_progress.SetLabel("")
            self.statusbar.SetStatusText("Grid search cancelled - results discarded.")
            return

        self._sim_store = sim_store if sim_store is not None else {}
        if not self._sim_store:
            self._region_cache = None
        self._last_click_ij = None
        self.residuals.show_message(
            "Double-click a point on the fitness map to show its "
            "residual maps." if self._sim_store else
            "No simulations were stored for this scan.")

        self._last_rmsd = rmsd
        self._last_rmsd_norm = rmsd_norm
        self._last_rmsd_n = total_n
        self._last_rmsd_SumSqRes = rmsd_SumSqRes
        self._last_rmsd_nk = rmsd_nk
        self._last_n_datasets = n_datasets
        self._last_grid1 = grid1
        self._last_grid2 = grid2
        self._last_param1 = param1
        self._last_param2 = param2
        self._refresh_display()

        msg = "Grid search finished."
        if n_failed:
            msg += f" {n_failed}/{total} grid points had no successful simulation (see console)."
        if self._sim_store:
            msg += (f" Simulations kept for {len(self._sim_store)} grid points "
                    f"({self.format_bytes(self._sim_store.nbytes)} in memory).")
        self.statusbar.SetStatusText(msg)

        if error_log:
            examples = "\n".join(f"  x{count}  {m}" for m, count in
                                 sorted(error_log.items(), key=lambda kv: -kv[1])[:5])
            # n_failed counts points that produced no simulation; the log can
            # also hold failures that cost no point at all (storing one, say),
            # so the headline follows the count rather than assuming it.
            headline = (f"{n_failed}/{total} grid points had no successful simulation."
                        if n_failed else
                        f"All {total} grid points simulated, but errors were logged.")
            wx.MessageBox(
                f"{headline}\n\n"
                f"Most common reasons:\n{examples}\n\n"
                f"(full detail was printed to the console as it happened)",
                "Grid Search - errors during scan", wx.OK | wx.ICON_WARNING)

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
        #self._sync_display_controls()
        self._refresh_display()

    def _refresh_display(self):
        """(Re)draw the heatmap from the coarse results cached by the last
        scan, in whichever statistic the Display pulldown selects, and
        recompute the interval readout to match. Cheap enough to call on
        every control change -- no re-scanning involved."""
        if self._last_rmsd is None:
            return

        drew = self.displayResults()
        # drew = False
        # if mode in self.PROBABILITY_MODES:
        #     drew = self._display_probability(mode)
        # if not drew:
        #     # Probability needs per-dataset residuals; fall back to the raw
        #     # fitness map if they are missing or the posterior degenerates.
        #     self._display_fitness(self.DISPLAY_RMSD if mode in self.PROBABILITY_MODES else mode)
        self.left_panel.Layout()
        self.Layout()

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
    def displayResults(self):
        
        CONF_LEVEL = 0.95
        if self._last_rmsd_SumSqRes is None or self._last_rmsd_nk is None:
            return False
        
        sel = self.display_choice.GetSelection()
        mode = ('RMSD', 'Normalized RMSD', 'Posterior Probability', 'log Posterior')[max(sel, 0)]
        grid1 = self._last_grid1
        grid2 = self._last_grid2
        param1 = self._last_param1
        param2 = self._last_param2
        label1 = self.param_full_label(param1)
        label2 = self.param_full_label(param2)
        method = ['none', 'linear', 'cubic'][self.interp_choice.GetSelection()]
        
        # CALCULATE CONFIDENCE INTERVALS FROM POSTERIOR
        # Using Jeffrey's prior for noise
        divisor = float(self.neff_divisor.GetValue())
        n_eff = self._last_rmsd_nk/divisor 
        # posterior probability p(theta|Data)~ sum(residuals^2)^(-n_eff/2)
        # log(p(theta|Data)) ~ - n_eff/2*log(sum(residuals^2))
        logJ = -(n_eff/2.0)*np.log(self._last_rmsd_SumSqRes)
        # in log space,posteriors from all datasets add, so we can sum them up to get total posterior
        log_post = logJ.sum(axis=2)
        log_post -= log_post.max()
        
        # ii, jj = np.unravel_index(np.argmax(log_post), log_post.shape)
        # mini = max(ii-1, 0) if ii<(log_post.shape[0]-1) else max(ii-2, 0)
        # minj = max(jj-1, 0) if jj<(log_post.shape[1]-1) else max(jj-2, 0)
        # cx = np.polyfit(grid1[mini:(mini+3)], log_post[mini:(mini+3), jj], 2)
        # cy = np.polyfit(grid2[minj:(minj+3)], log_post[ii, minj:(minj+3)], 2)
        # x0_fit = -cx[1] / (2 * cx[0])
        # sigmaX_fit = np.sqrt(1 / (2 * abs(cx[0])))
        # y0_fit = -cy[1] / (2 * cy[0])
        # sigmaY_fit = np.sqrt(1 / (2 * abs(cy[0])))
        # zp = ndtri((1 + CONF_LEVEL) / 2) # quantile
        # half_width = [zp * sigmaX_fit, zp * sigmaY_fit]
        # A =   (grid1[ii] - x0_fit)**2 / (2 * sigmaX_fit**2) \
        #     + (grid2[jj] - y0_fit)**2 / (2 * sigmaY_fit**2)
        # log_post -= A
        # # techincally, the threshold is -1/2 chi^2_k(CONF) (k-dimentionality, so 2)
        # # but for that case, the threshold is equal to log(1-CONF)
        # thresh_log_post = np.log(1-CONF_LEVEL)
        # thresh_post     = np.exp(thresh_log_post)
        #chi2 = 
        # self.lbl_ci.SetLabel(f"{label1} = {x0_fit:.5g} ± {half_width[0]:.5g}  ({CONF_LEVEL*100:.2f}%)\n"
        #                      f"{label2} = {y0_fit:.5g} ± {half_width[1]:.5g}  ({CONF_LEVEL*100:.2f}%)")
        
        ########### WHAT TO PLOT ##############################################
        contour_levels = None
        
        fine1, fine2, fine_logp, used_method  = self._interpolate_grid(method, grid1, grid2, log_post)
        
        if mode == 'RMSD':
            if self._last_rmsd is None: return False
            out = self._interpolate_grid(method, grid1, grid2, self._last_rmsd)
            fine_vals = out[2]
            value_label = "RMSD"
        elif mode == 'Normalized RMSD':
            if self._last_rmsd_norm is None: return False
            out = self._interpolate_grid(method, grid1, grid2, self._last_rmsd_norm)
            fine_vals = out[2]
            value_label = "norm_RMSD"
        elif mode == 'Posterior Probability':
            fine_vals = np.exp(fine_logp) 
            value_label = "p"
            thresh_post = (1-CONF_LEVEL)/np.sum(fine_vals)
            contour_levels = [(thresh_post, "95%")]
            fine_vals  /=np.sum(fine_vals)
        elif mode == 'log Posterior':
            fine_vals = fine_logp-np.max(fine_logp)
            value_label = "log(p)"
            thresh_post = np.log((1-CONF_LEVEL))
            contour_levels = [(thresh_post, "95%")]
            # fine_vals  /=np.sum(fine_vals)
            
        else:
            return False
        ########### CREDIBLE INTERVALS ########################################
        # do gausian fit over the max of interpolated map using polynomial method
        ii, jj = np.unravel_index(np.argmax(fine_logp), fine_logp.shape)
        
        npts = 5 # how many points to take on each side of the high point
        
        n0 = fine_logp.shape[0]
        n1 = fine_logp.shape[1]
        width = 2 * npts + 1
        
        mini = min(max(ii - npts, 0), max(n0 - width, 0))
        minj = min(max(jj - npts, 0), max(n1 - width, 0))
        
        cx = np.polyfit(fine1[mini:mini+width], fine_logp[mini:mini+width, jj], 2)
        cy = np.polyfit(fine2[minj:minj+width], fine_logp[ii, minj:minj+width], 2)
        
        x0_fit = -cx[1] / (2 * cx[0])
        sigmaX_fit = np.sqrt(1 / (2 * abs(cx[0])))
        y0_fit = -cy[1] / (2 * cy[0])
        sigmaY_fit = np.sqrt(1 / (2 * abs(cy[0])))
        
        zp = ndtri((1 + CONF_LEVEL) / 2) # quantile
        
        half_width = [zp * sigmaX_fit, zp * sigmaY_fit]
        # A =   (fine1[ii] - x0_fit)**2 / (2 * sigmaX_fit**2) \
        #     + (fine2[jj] - y0_fit)**2 / (2 * sigmaY_fit**2)
        
        # techincally, the threshold is -1/2 chi^2_k(CONF) (k-dimentionality, so 2)
        # but for that case, the threshold is equal to log(1-CONF)
        

        self.lbl_ci.SetLabel(f"{label1} = {x0_fit:.5g} ± {half_width[0]:.5g}  ({CONF_LEVEL*100:.2f}%)\n"
                             f"{label2} = {y0_fit:.5g} ± {half_width[1]:.5g}  ({CONF_LEVEL*100:.2f}%)")
        
            
            #bi, bj = np.unravel_index(np.nanargmax(fine_vals), fine_vals.shape)
            
        # else:
        #     #bi, bj = np.unravel_index(np.nanargmin(fine_vals), fine_vals.shape)
        #     #rmsd_min = fine_vals[bi, bj]
        
        background_levels = None
        vmin, vmax = float(np.nanmin(fine_vals)), float(np.nanmax(fine_vals))
        if (vmax > vmin) and self.map_contour_on and self.chk_rmsd_contours.GetValue():
            background_levels = list(np.linspace(vmin, vmax, self.map_contour_n + 2)[1:-1])


        as_contours = self.chk_rmsd_contours.GetValue()
        rmsd_cw = float(self.rmsd_contour_width.GetValue())
        mark_resid, mark_loaded = self._marker_points()
        self.heatmap.show_results(fine_vals, fine1, fine2, label1, label2,
                                  contour_levels=contour_levels,
                                  background_levels=background_levels,
                                  raw_points=(grid1, grid2), subtitle=mode,
                                  half_width=half_width,
                                  as_contours=as_contours,
                                  rmsd_contour_width=rmsd_cw,
                                  value_label=value_label,
                                  best_value=(x0_fit, y0_fit),
                                  bg_contour_color=self.map_contour_color,
                                  bg_contour_width=self.map_contour_width,
                                  mark_residual=mark_resid, 
                                  mark_loaded=mark_loaded)
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
        """
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

    def on_display_mode_changed(self, event):
        if event is not None:
            event.Skip()      # EVT_KILL_FOCUS must keep propagating
        #self._sync_display_controls()
        self._refresh_display()
        if self._last_click_ij is not None:
            self._render_residuals(*self._last_click_ij)

    def _on_doubleclick(self, x, y):
        if not self._show_residuals:
            return
        if self._last_grid1 is None or self._last_grid2 is None:
            return
        if not self._sim_store or self._region_cache is None:
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
        if self._last_grid1 is None or not self._sim_store or self._region_cache is None:
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
        for region, sim in zip(self._region_cache, slot):
            if region is None or sim is None:
                continue
            # Blanked outside the tilted box: those points were never fitted,
            # so showing their residual would invite reading a miss into a
            # region the scan never scored.
            resid = np.where(region['mask'], region['exp'] - sim, np.nan)
            if region['null_rms'] > 0.0:
                scale = 1.0 / region['null_rms']
                resid = resid * scale
                sim = sim * scale
            # Head each map with the field, matching how the main window
            # titles its spectra (see MainFrame's 'B$_0$={field} mT'),
            # rather than with the file name -- which of several tau/field
            # measurements this is, is the useful identity here.
            field = region.get('field')
            if field is not None:
                title = f"B$_0$={field:g} mT"
            else:
                title = ''
            dx = region['x'][1]-region['x'][0]
            dy = region['y'][1]-region['y'][0]
            entries.append({
                'title': title,
                'resid': resid,
                'sim': sim,
                'x': region['x'],
                'y': region['y'],
                'extent': [region['x'][0]-dx/2, region['x'][-1]+dx/2, 
                           region['y'][0]-dy/2, region['y'][-1]+dy/2],
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
        #self.main.Sys.setFromCtrl(main_params)
        #self.main.runSim()
        msg = (f"Loaded {self.param_full_label(self._last_param1)}={v1:.5g}, "
               f"{self.param_full_label(self._last_param2)}={v2:.5g}")
        self.on_refresh(None)
        self._refresh_display()
        self.statusbar.SetStatusText(msg)

    # ------------------------------------------------------------------
    def on_close(self, event):
        self._closing = True
        if self._running:
            self._cancel_event.set()
            pool = self._pool
            if pool is not None:
                # Workers are children of this process; leaving them running
                # after the window is gone would keep burning cores with
                # nothing left to collect the results.
                try:
                    pool.terminate()
                except Exception:
                    pass
        # The one and only place stored simulations are erased. Every scan in
        # this session wrote its own directory and none of them were cleared
        # as they went, so the sweep picks up the ones this panel is no longer
        # holding a handle to -- earlier scans, and any that were cancelled.
        for store in self._sim_stores:
            store.dispose()
        self._sim_stores = []
        self._sim_store = {}
        self._region_cache = None
        event.Skip()

    # ------------------------------------------------------------------
    # Parameter descriptors and the Sys dictionary
    # ------------------------------------------------------------------
    def flatten_sys_params(self, sys_dict):
        """Scannable scalar leaves of a Sys parameter dict. Lives in
        grid_search_worker so a worker process can reach it without wx."""
        return gsw.flatten_sys_params(sys_dict)
    def qualified_label(self, p):
        return gsw.qualified_label(p)
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
        return gsw.param_leaf_label(p)
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
        gsw.apply_param_value(sys_dict, param, value)
    def apply_functions_to_sys_dict(self, sys_dict, func_rows, variable_values, error_log=None):
        """Evaluate the Function Modifiers into `sys_dict` in place. Lives in
        grid_search_worker, which is where the scan actually calls it from."""
        return gsw.apply_functions_to_sys_dict(sys_dict, func_rows, variable_values,
                                               error_log=error_log)

    def _silent_error_func(self, message):
        gsw.silent_error_func(message)

    def format_bytes(self, nbytes):
        for unit in ('bytes', 'kilobytes', 'megabytes', 'gigabytes'):
            if nbytes < 1024 or unit == 'gigabytes':
                return f"{nbytes:.0f} {unit}" if unit == 'bytes' else f"{nbytes:.2f} {unit}"
            nbytes /= 1024.0
    def _log_gridpoint_error(self, error_log, where, exc):
        """Count a failure and print it the first time it is seen. Workers do
        the counting only -- see gsw.log_gridpoint_error -- because a spawned
        process usually has no console to print to."""
        if error_log is None:
            return
        msg = f"{where}: {type(exc).__name__}: {exc}"
        if msg not in error_log:
            print(f"[GridSearch] {msg}")
        error_log[msg] = error_log.get(msg, 0) + 1

    def _merge_error_log(self, error_log, incoming):
        """Fold one grid point's failures into the scan-wide log, printing
        each distinct message the first time it turns up."""
        for msg, count in incoming.items():
            if msg not in error_log:
                print(f"[GridSearch] {msg}")
                error_log[msg] = 0
            error_log[msg] += count

    # ------------------------------------------------------------------
    # Interpolating the scanned grid, and confidence regions
    # ------------------------------------------------------------------
    def _interpolate_grid(self, method, grid1, grid2, results, n_out=160):
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
        elif method == 'linear':
            used_method = 'linear'
            interp = RegularGridInterpolator((grid1, grid2), filled, method='linear',
                                              bounds_error=False, fill_value=None)
            g1, g2 = np.meshgrid(fine1, fine2, indexing='ij')
            fine_vals = interp(np.stack([g1.ravel(), g2.ravel()], axis=-1)).reshape(g1.shape)
        else:
            used_method = 'none'
            fine_vals = results
            fine1, fine2 = grid1, grid2
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

    def _axis_confidence_bounds(self, axis_vals, profile, thresh):
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
        # density at which the credible level paSumSqRess conf.
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
