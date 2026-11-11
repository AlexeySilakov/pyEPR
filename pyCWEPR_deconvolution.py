"""Lorentzian deconvolution of a power series, for pyCWEPR.

A side window opened from the main toolbar, built like the power saturation window: the spectra come
from the very results the main plot draws, so a background, phase, baseline or integration-range
change there moves them here too.

Nothing is fitted. One component measured at several powers gets two knobs for the whole series,
gamma(P) = gamma_0*(1 + P/P_sat)^b with b at 1/2 for the Bloch law, and every spectrum takes its own
power from its file. Turn the knobs until the deconvolved traces collapse onto one curve: that curve
is the inhomogeneous envelope. Traces are divided by the square root of the power, so shapes are
what is compared and the amplitude saturation stays visible as a spread.

Parity decomposition splits each trace about the field its own microwave frequency puts the chosen g
value at, so the even/odd ratio measures the dispersion admixture with no model at all.
"""

import io
import os

import wx
import wx.lib.scrolledpanel as scrolled

import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

from theme import theme
from theme.theme import PillButton, PillCheckBox, SegmentedPill

# the Bohr magneton over the Planck constant, Hz/T, and the free-electron g factor
MU_B_OVER_H = 1.39962449e10
G_FREE = 2.00231930436


class DeconvolutionFrame(wx.Frame):
    """Separate window: Lorentzian deconvolution of a power series by FFT, and the even/odd parity
    decomposition of every spectrum. Display scale, Y shift and the Compare normalisation of the
    main window are display-only and stay out; the vertical separation here is this window's own."""

    KEYS = ('g0', 'Psat', 'b')
    LABELS = ("γ₀, mT", "P_sat, mW", "b")
    LIMITS = {'g0': (0.0, 1e4), 'Psat': (1e-9, 1e9), 'b': (-5.0, 5.0)}
    EQUATION = (r"$S' = \mathrm{IFFT}\left[\,\hat{S}(k)\;/\;\hat{L}(k)\right]$"
                "\n" r"$\gamma(P)=\gamma_0\,(1+P/P_{\mathrm{sat}})^{\,b}$")

    def __init__(self, parent, mainWindow):
        super().__init__(parent, title="Lorentzian deconvolution", size=(1180, 760))
        self.mf = mainWindow
        self.procs = []           # process() result per spectrum, handed over by the main plot
        self.pars = {'g0': 0.05, 'Psat': 1.0, 'b': 0.5}
        self.steps = {}           # parameter -> arrow step chosen by right-clicking its value box
        self.started = False      # the first non-empty spectrum set gets a starting width
        self.trace = 'deriv'      # 'deriv' | 'abs', the trace everything here reads
        self.view = 'spectra'     # 'spectra' | 'parity' | 'width'
        self.logscale = {'x': False, 'y': False}
        self.ax = None

        splitter = wx.SplitterWindow(self, style=wx.SP_3D | wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(140)
        left = scrolled.ScrolledPanel(splitter, style=wx.BORDER_NONE)
        right = wx.Panel(splitter)
        cards = wx.BoxSizer(wx.VERTICAL)

        # Signal card: which trace everything reads
        card = theme.Card(left, "Signal")
        self.seg_trace = SegmentedPill(card, ["Deriv.", "Absorp."])
        self.seg_trace.SetToolTip("The processed trace everything here reads: the derivative, or the\n"
                                  "absorption (1st integral). Differentiating is a factor ik in Fourier,\n"
                                  "which commutes with the deconvolution, so either gives the same widths.\n"
                                  "The derivative needs no integration baseline.")
        self.seg_trace.Bind(wx.EVT_RADIOBOX, self.on_signal)
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(card, label="Trace:", size=(76, -1)), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        row.Add(self.seg_trace, 1, wx.EXPAND)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # Deconvolution card: the equation, then the three knobs
        card = theme.Card(left, "Deconvolution")
        self.bmp_equation = wx.StaticBitmap(card, bitmap=self.equation_bitmap())
        self.bmp_equation.SetToolTip("γ₀ is the homogeneous width at zero power and P_sat the power where the\n"
                                     "saturation factor s reaches 1. b is 1/2 for the Bloch law; free it only\n"
                                     "to let the data argue. Two numbers serve the whole series, so the traces\n"
                                     "collapsing onto one curve is a real test rather than a fitted result.")
        card.body.Add(self.bmp_equation, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 6)
        grid = wx.FlexGridSizer(cols=2, vgap=2, hgap=6)
        grid.AddGrowableCol(1)
        self.param_widgets = {}
        tips = {'g0': "The width of the Lorentzian taken out at zero power, in mT.\n"
                      "Right-click to set the arrow step.",
                'Psat': "The power at which the saturation factor s = P/P_sat reaches 1, which with\n"
                        "b = 1/2 is where the homogeneous width has grown by √2. Not the\n"
                        "half-saturation power of the power saturation window, which is the smaller\n"
                        "P½ = P_sat·(2^(1/ε) − 1) and refers to the amplitude, not the width.\n"
                        "Right-click to set the arrow step.",
                'b': "The exponent of (1 + P/P_sat). 1/2 is the Bloch saturation law;\n"
                     "step it away from 1/2 only to see whether the data prefers something else.\n"
                     "Right-click to set the arrow step."}
        for key, label in zip(self.KEYS, self.LABELS):
            grid.Add(wx.StaticText(card, label=label, size=(76, -1)), 0, wx.ALIGN_CENTER_VERTICAL)
            spin = wx.SpinCtrlDouble(card, min=self.LIMITS[key][0], max=self.LIMITS[key][1],
                                     inc=0.01, initial=self.pars[key], size=(118, -1))
            spin.SetDigits(4)
            spin.SetToolTip(tips[key])
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=key: self.on_param(k))
            theme.bind_spin_enter(spin, lambda k=key: self.on_param(k))
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                spin.Bind(evt, lambda e, s=spin, k=key: self.on_step(k, s))
            grid.Add(spin, 0, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
            self.param_widgets[key] = spin
        card.body.Add(grid, 0, wx.EXPAND | wx.ALL, 3)
        self.lbl_stat = wx.StaticText(card, label="")
        card.body.Add(self.lbl_stat, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # Plot card: which view, how far apart the traces sit, what else is drawn
        card = theme.Card(left, "Plot")
        self.seg_view = SegmentedPill(card, ["Spectra", "Parity", "Width"])
        self.seg_view.SetToolTip("Spectra: the traces, deconvolved, divided by √P.\n"
                                 "Parity: their even and odd parts about the centre g.\n"
                                 "Width: γ(P) over the measured powers.")
        self.seg_view.Bind(wx.EVT_RADIOBOX, self.on_signal)
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(card, label="View:", size=(76, -1)), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        row.Add(self.seg_view, 1, wx.EXPAND)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(card, label="Separation:", size=(76, -1)), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.spin_shift = wx.SpinCtrlDouble(card, min=-100, max=100, inc=0.1, initial=0)
        self.spin_shift.SetDigits(3)
        self.spin_shift.SetToolTip("Vertical gap between consecutive traces, in units of the tallest one,\n"
                                   "as on the main plot. Fan them out to dial the knobs in, then set it\n"
                                   "back to 0 to see whether they really collapse onto each other.\n"
                                   "Right-click to set the arrow step.")
        self.spin_shift.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e: self.update_graph())
        theme.bind_spin_enter(self.spin_shift, lambda: self.update_graph())
        for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
            self.spin_shift.Bind(evt, lambda e: self.mf.tabulated_panel.on_step_menu(
                self.spin_shift, (0.01, 0.05, 0.1, 0.25, 0.5, 1), 3))
        row.Add(self.spin_shift, 1, wx.EXPAND)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        self.chk_grid = PillCheckBox(card, "Grid")
        self.chk_grid.SetValue(True)
        self.chk_grid.SetToolTip("Axes grid, on both scales.")
        self.chk_grid.Bind(wx.EVT_CHECKBOX, lambda e: self.update_graph())
        self.chk_original = PillCheckBox(card, "Original")
        self.chk_original.SetValue(True)
        self.chk_original.SetToolTip("Draw each trace as it came in, faintly, under its deconvolved self.")
        self.chk_original.Bind(wx.EVT_CHECKBOX, lambda e: self.update_graph())
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.chk_grid, 0, wx.RIGHT, 4)
        row.Add(self.chk_original, 0)
        card.body.Add(row, 0, wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # Parity card: the centre, as a g value so each frequency gets its own field
        card = theme.Card(left, "Parity decomposition")
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(card, label="Centre g:", size=(76, -1)), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self.spin_g = wx.SpinCtrlDouble(card, min=0.01, max=100, inc=0.0001, initial=G_FREE)
        self.spin_g.SetDigits(5)
        self.spin_g.SetToolTip("The split is taken about the field this g sits at in each spectrum, using\n"
                               "that spectrum's own microwave frequency, so a frequency drift across the\n"
                               "series does not smear the decomposition.\n"
                               "Right-click to set the arrow step.")
        self.spin_g.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e: self.update_graph())
        theme.bind_spin_enter(self.spin_g, lambda: self.update_graph())
        for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
            self.spin_g.Bind(evt, lambda e: self.mf.tabulated_panel.on_step_menu(
                self.spin_g, (1e-5, 1e-4, 1e-3, 0.01, 0.1), 5))
        row.Add(self.spin_g, 1, wx.EXPAND)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        self.lbl_parity = wx.StaticText(card, label="")
        card.body.Add(self.lbl_parity, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        self.btn_save = PillButton(left, "Save CSV", color_key="pill_calm", draw_badge=theme.draw_textlines_badge)
        self.btn_save.SetToolTip("Save the traces, their deconvolved selves and the even and odd parts as a\n"
                                 "comma-separated file that Excel, Origin, pandas or pyCWEPR can read.")
        self.btn_save.Bind(wx.EVT_BUTTON, self.on_save)
        cards.Add(self.btn_save, 0, wx.EXPAND | wx.ALL, 6)
        left.SetSizer(cards)
        left.SetupScrolling(scroll_x=True, rate_y=10)    # a narrow pane scrolls instead of clipping a card
        left.Scroll(0, 0)

        self.figure = Figure(figsize=(2, 2), dpi=100, layout='constrained')
        self.canvas = FigureCanvasWxAgg(right, -1, self.figure)
        self.canvas.mpl_connect('button_release_event', self.on_axis_menu)
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.canvas, 1, wx.EXPAND)
        box.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        right.SetSizer(box)
        self.toolbar.update()

        self.statusbar = theme.StatusStrip(self, widths=[-1, 360])
        self.statusbar.SetStatusText("Ready", 0)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('figure_leave_event', lambda e: self.statusbar.SetStatusText("", 1))

        splitter.SplitVertically(left, right, sashPosition=400)
        splitter.SetSashGravity(0)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND)
        sizer.Add(self.statusbar, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self.redraw_debouncer.on_size)
        self.Bind(wx.EVT_CLOSE, lambda e: (setattr(self.mf, 'deconv_window', None), self.Destroy()))
        self.apply_theme()

    # -----------------------------------------------------------------------
    def apply_theme(self):
        theme.apply_theme_to_window(self, theme.get_theme())
        self.bmp_equation.SetBitmap(self.equation_bitmap())     # it is drawn in the theme's text colour
        theme.style_figure(self.figure, self.figure.axes)
        self.canvas.draw_idle()
        self.Refresh()

    def equation_bitmap(self):
        """The operation as mathtext in the theme's text colour, at this window's DPI."""
        t = theme.theme_colors()
        scale = self.GetDPIScaleFactor()
        figure = Figure(figsize=(3.4, 0.8), dpi=round(96*scale), facecolor=t["card_bg"])
        figure.text(0.5, 0.5, self.EQUATION, color=t["text"], fontsize=12, ha='center', va='center')
        buf = io.BytesIO()
        figure.savefig(buf, format='png', bbox_inches='tight', pad_inches=0.08, facecolor=t["card_bg"])
        buf.seek(0)
        bitmap = wx.Bitmap(wx.Image(buf, wx.BITMAP_TYPE_PNG))
        bitmap.SetScaleFactor(scale)        # a static bitmap composites its alpha against white on MSW
        return bitmap

    # -----------------------------------------------------------------------
    def gamma(self, P):
        """The width [mT] of the Lorentzian taken out at power P [mW], from the three knobs."""
        return self.pars['g0']*(1 + P/max(self.pars['Psat'], 1e-12))**self.pars['b']

    def deconvolve(self, x, y, P):
        gg = self.gamma(P)
        Lor = 2*gg/np.pi/( 4*(x-np.mean(x))**2+gg**2)

        spec = np.fft.fft(y - np.mean(y))
        fgg  = np.fft.fft(Lor)
        result = np.fft.ifft(spec/fgg)
        return result

    def centre_field(self, freq):
        """The field [mT] the chosen g sits at for a spectrum of this microwave frequency [GHz]."""
        g = self.spin_g.GetValue()
        return 1e12*freq/(MU_B_OVER_H*g) if freq > 0 and g > 0 else np.nan

    def spectra(self):
        """(index, power, field, trace/sqrt(P), deconvolved) of the shown spectra, and the warnings."""
        mf = self.mf
        out, warn, nopower, bad = [], [], 0, 0
        for ii, (dd, p) in enumerate(zip(mf.Data, self.procs)):
            if not dd['show']:
                continue
            if dd['power'] <= 0:
                nopower += 1
                continue
            x, y = ((p['x'], p['y']) if self.trace == 'deriv' else (p['xi'], p['abs'] - p['bl2']))
            x = np.asarray(x, dtype=float)
            P = float(dd['power'])
            y = np.asarray(y, dtype=float)/np.sqrt(P)
            # ifft returns complex; a real trace over a real kernel leaves nothing in the imaginary part
            dec = np.real(self.deconvolve(x, y, P))
            if not np.all(np.isfinite(dec)):
                bad += 1
                print(f"[pyCWEPR] deconvolve returned {np.count_nonzero(~np.isfinite(dec))} of "
                      f"{dec.size} non-finite points for '{dd['title']}' at {P:g} mW")
            out.append((ii, P, x[:dec.size], y[:dec.size], dec))
        if nopower:
            warn.append(f"{nopower} shown spectra state no MW power (0) and stay out")
        if bad:
            warn.append(f"{bad} deconvolved traces hold inf or nan and cannot be drawn")
        out.sort(key=lambda r: r[1])
        return out, warn

    def split_parity(self, x, y, centre):
        """(field, even, odd) of a trace mirrored about `centre`, where it has data on both sides."""
        if not np.isfinite(centre):
            return None
        reach = min(centre - x[0], x[-1] - centre)
        if reach <= 0:
            return None
        d = np.linspace(0, reach, max(int(0.5*x.size), 64))
        right, left = np.interp(centre + d, x, y), np.interp(centre - d, x, y)
        even, odd = 0.5*(right + left), 0.5*(right - left)
        return (np.concatenate([centre - d[::-1], centre + d[1:]]),
                np.concatenate([even[::-1], even[1:]]),
                np.concatenate([-odd[::-1], odd[1:]]))

    def refresh(self, procs):
        """Freshly processed spectra from the main window: new traces, and a redraw."""
        self.procs = procs
        data = self.spectra()[0]
        if not self.started and data:
            self.started = True
            # a hundredth of the sweep starts under-deconvolved rather than ringing
            self.pars['g0'] = round(0.01*float(np.ptp(data[0][2])), 4)
            self.pars['Psat'] = float(np.median([d[1] for d in data]))
            for key, spin in self.param_widgets.items():
                v = self.pars[key]
                step = self.steps.get(key) or (10.0**(np.floor(np.log10(abs(v))) - 2) if v else 0.01)
                # never so few decimals that a value this small could not be typed back in
                spin.SetDigits(max(int(np.clip(4 - np.floor(np.log10(abs(v))) if v else 4, 2, 12)),
                                   int(np.ceil(-np.log10(step)))))
                spin.SetIncrement(float(step))
                spin.SetValue(float(np.clip(v, self.LIMITS[key][0], self.LIMITS[key][1])))
        self.update_graph()

    # -----------------------------------------------------------------------
    def update_graph(self):
        data, warn = self.spectra()
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        self.ax = ax
        colors = theme.theme_colors()

        if not data:
            ax.text(0.5, 0.5, "No shown spectrum states an MW power", ha='center', va='center',
                    transform=ax.transAxes)
            self.lbl_stat.SetLabel("No spectra: show one that states its MW power")
            self.lbl_parity.SetLabel("")
        elif self.view in ('spectra', 'parity'):
            # one step between consecutive traces, in units of the tallest one, as on the main plot
            spans = [float(np.ptp(d[4])) for d in data]
            height = max([s for s in spans if np.isfinite(s) and s > 0], default=1.0)
            step = self.spin_shift.GetValue()*height
            ratios = []
            for kk, (ii, P, x, y, dec) in enumerate(data):
                col = self.mf.trace_color(self.mf.Data[ii])
                off = kk*step
                if self.view == 'spectra':
                    if self.chk_original.GetValue():
                        ax.plot(x, y + off, lw=0.7, color=col, alpha=0.4)
                    ax.plot(x, dec + off, lw=1.1, color=col, label=f"{P:g} mW")
                    continue
                got = self.split_parity(x, dec, self.centre_field(self.mf.Data[ii]['freq']))
                if got is None:
                    warn.append(f"'{self.mf.Data[ii]['title']}' does not reach the centre g")
                    continue
                field, even, odd = got
                ax.plot(field, odd + off, lw=1.0, color=col, label=f"{P:g} mW: odd")
                ax.plot(field, even + off, lw=1.2, ls='--', color=col, label=f"{P:g} mW: even")
                span = float(np.max(np.abs(odd)))
                ratios.append(float(np.max(np.abs(even)))/span if span > 0 else np.nan)
            ax.set_ylabel(("derivative" if self.trace == 'deriv' else "absorption") + r" / $\sqrt{P}$" +
                          ("" if self.view == 'spectra' else ",  even (dashed) and odd (solid)"))
            self.lbl_parity.SetLabel("|even|/|odd| per spectrum:\n" + "  ".join(f"{r:.3f}" for r in ratios) +
                                     "\nflat and non-zero is an intrinsic admixture" if ratios else "")
        else:
            powers = np.array([d[1] for d in data])
            fine = np.logspace(np.log10(powers.min()) - 0.3, np.log10(powers.max()) + 0.3, 400)
            ax.plot(fine, self.gamma(fine), lw=1.3, color=colors["ax_fg"], label=r"$\gamma(P)$")
            ax.plot(powers, self.gamma(powers), 'o', ms=7, ls='none', color=colors["accent"],
                    label="measured powers")
            ax.set_ylabel(r"$\gamma$,  mT")

        ax.set_xlabel("field, mT" if self.view in ('spectra', 'parity') else "P, mW")
        ax.set_xscale('log' if self.logscale['x'] and self.view == 'width' else 'linear')
        ax.set_yscale('log' if self.logscale['y'] else 'linear')
        ax.grid(False)          # matplotlib turns the grid back on if line properties come with False
        if self.chk_grid.GetValue():
            ax.grid(True, which='both', lw=0.5, alpha=0.4)
        if data:
            handles = ax.get_legend_handles_labels()[0]
            if handles and len(handles) <= 14:
                ax.legend(fontsize='small', frameon=False, loc='best')
            lo, hi = (self.gamma(min(d[1] for d in data)), self.gamma(max(d[1] for d in data)))
            self.lbl_stat.SetLabel(f"{len(data)} spectra,  γ from {lo:.4g} to {hi:.4g} mT\n"
                                   f"turn γ₀ and P_sat until the traces collapse;\n"
                                   f"ringing and negative wings mean you have gone past")
        theme.style_figure(self.figure, [ax])
        self.toolbar.update()
        self.canvas.draw_idle()
        self.statusbar.SetStatusText("; ".join(warn) if warn else "Ready", 0)
        self.Layout()

    # -----------------------------------------------------------------------
    def on_mouse_move(self, event):
        if event.inaxes is None or event.xdata is None:
            self.statusbar.SetStatusText("", 1)
            return
        unit = "mT" if self.view in ('spectra', 'parity') else "mW"
        self.statusbar.SetStatusText(f"x = {event.xdata:.6g} {unit}   y = {event.ydata:.6g}", 1)

    def on_step(self, key, spin):
        """Right-click on a parameter's value box: its arrow step, which then survives everything."""
        self.mf.tabulated_panel.on_step_menu(spin, digits=2)
        self.steps[key] = spin.GetIncrement()

    def on_signal(self, event=None):
        self.trace = ('deriv', 'abs')[self.seg_trace.GetSelection()]
        self.view = ('spectra', 'parity', 'width')[self.seg_view.GetSelection()]
        self.update_graph()

    def on_param(self, key):
        self.pars[key] = self.param_widgets[key].GetValue()
        self.update_graph()

    def on_axis_menu(self, event):
        """Right click on an axis: logarithmic or linear scale for it."""
        if event.button != 3 or event.inaxes is not None or str(self.toolbar.mode) or self.ax is None:
            return
        renderer = self.canvas.get_renderer()
        for which, axis in (('x', self.ax.xaxis), ('y', self.ax.yaxis)):
            bb = axis.get_tightbbox(renderer)
            if bb is None or not bb.contains(event.x, event.y):
                continue
            menu = wx.Menu()
            ids = {}
            for name, log in (("Linear", False), ("Logarithmic", True)):
                item = menu.AppendCheckItem(wx.ID_ANY, f"{which.upper()} axis: {name}")
                item.Check(self.logscale[which] == log)
                ids[item.GetId()] = log
            choice = self.canvas.GetPopupMenuSelectionFromUser(menu)
            menu.Destroy()
            if choice in ids and ids[choice] != self.logscale[which]:
                self.logscale[which] = ids[choice]
                self.update_graph()
            return

    def on_save(self, event):
        """Traces, deconvolved traces and parity parts as CSV, the settings as '# key: value' lines."""
        data, _ = self.spectra()
        if not data:
            wx.MessageBox("There are no spectra to save.", "Save CSV", wx.OK | wx.ICON_INFORMATION)
            return
        with wx.FileDialog(self, "Save Lorentzian deconvolution", defaultDir=self.mf.file_browser.path,
                           defaultFile="deconvolution_pyCWEPR.csv",
                           wildcard="CSV (*.csv)|*.csv|All files (*.*)|*.*",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            pathname = dlg.GetPath()
        meta = ["pyCWEPR Lorentzian deconvolution", f"trace: {self.trace}",
                "operation: S' = IFFT[ FFT(S - mean S) / FFT(Lorentzian of width gamma(P)) ]",
                "gamma(P) = gamma_0*(1 + P/P_sat)^b [mT]",
                f"gamma_0 [mT]: {self.pars['g0']:.10g}", f"P_sat [mW]: {self.pars['Psat']:.10g}",
                f"b: {self.pars['b']:.10g}", f"centre g: {self.spin_g.GetValue():.10g}",
                "traces divided by sqrt(P)"]
        meta += [f"gamma at {P:g} mW [mT]: {self.gamma(P):.10g}" for _, P, _, _, _ in data]
        cols = ["title", "P [mW]", "B [mT]", "signal/sqrt(P)", "deconvolved/sqrt(P)",
                "parity field [mT]", "even", "odd"]
        try:
            with open(pathname, 'w', encoding='utf-8', newline='') as f:
                f.write("".join(f"# {line}\n" for line in meta))
                f.write(",".join(cols) + "\n")
                for ii, P, x, y, dec in data:
                    got = self.split_parity(x, dec, self.centre_field(self.mf.Data[ii]['freq']))
                    title = '"' + self.mf.Data[ii]['title'].replace('"', "'") + '"'
                    for jj in range(x.size):
                        cells = [title, f"{P:.10g}", f"{x[jj]:.10g}", f"{y[jj]:.10g}", f"{dec[jj]:.10g}"]
                        cells += ([f"{got[0][jj]:.10g}", f"{got[1][jj]:.10g}", f"{got[2][jj]:.10g}"]
                                  if got is not None and jj < got[0].size else [""]*3)
                        f.write(",".join(cells) + "\n")
            self.mf.statusbar.SetStatusText(f"Lorentzian deconvolution saved to {pathname}", 0)
            self.mf.file_browser.refresh_file_list(select_name=os.path.basename(pathname))
        except OSError as fault:
            wx.MessageBox(f"Could not write '{pathname}':\n{fault}", "Save failed", wx.OK | wx.ICON_ERROR)
