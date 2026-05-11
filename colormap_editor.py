"""
colormap_editor.py  —  Reusable wxPython dialog for matplotlib colormaps.

USAGE
-----
    import wx
    from colormap_editor import ColormapEditorDialog, stops_to_cmap, STANDARD_CMAPS

    with ColormapEditorDialog(parent, theme="dark") as dlg:
        if dlg.ShowModal() == wx.ID_OK:
            cmap  = dlg.get_colormap()
            name  = dlg.get_colormap_name()
            stops = dlg.get_stops()

UI
--
  Interactive matplotlib graph:
    • X axis  = stop index (evenly spaced)
    • Y axis  = colormap position 0-1
    • Markers = colour stops drawn with scatter (screen-space circles, always round)
    • Drag up/down to change position; first/last locked at 0 and 1
    • Right-click  → colour picker
    • Double-click → floating SpinCtrlDouble to type an exact value; Enter confirms
    • Colorbar on the right updates live
    • 🌙 Dark / ☀ Light button in toolbar toggles theme live
    • theme="light"|"dark" constructor arg sets initial theme
"""

import copy
import json
import pathlib
from collections import deque

import wx
import numpy as np
import matplotlib
matplotlib.use("WXAgg")
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

DEFAULT_CMAP_FILE = pathlib.Path(__file__).parent / ".colormap_editor_default.json"
_SENTINEL = object()

# ── colour themes ─────────────────────────────────────────────────────────────
THEMES = {
    "light": dict(
        # wx panel / dialog colours
        dlg_bg       = wx.Colour(240, 240, 240),
        panel_bg     = wx.Colour(240, 240, 240),
        text_fg      = wx.Colour(20,  20,  20),
        label_fg     = wx.Colour(60, 120,  60),   # "saved as default" label
        btn_bg       = wx.Colour(200,  200,  200),
        # matplotlib figure
        fig_bg       = "#efefef",
        ax_bg        = "#f8f8f8",
        ax_fg        = "#222222",   # spines, ticks, labels
        grid_col     = "#cccccc",
        line_col     = "#555555",
        # scatter edge colours (inactive stops)
        edge_dark    = "#333333",
        edge_light   = "#ffffff",
        # active stop highlight
        active_edge  = "#555555",
        active_lw    = 3.5,
        inactive_lw  = 1.8,
    ),
    "dark": dict(
        dlg_bg       = wx.Colour(30,  30,  35),
        panel_bg     = wx.Colour(30,  30,  35),
        text_fg      = wx.Colour(220, 220, 220),
        label_fg     = wx.Colour(100, 200, 120),
        btn_bg       = wx.Colour(60,  50,  50),
        fig_bg       = "#1e1e23",
        ax_bg        = "#28282f",
        ax_fg        = "#cccccc",
        grid_col     = "#3a3a44",
        line_col     = "#aaaaaa",
        edge_dark    = "#dddddd",
        edge_light   = "#111111",
        active_edge  = "#AAAAAA",
        active_lw    = 3.5,
        inactive_lw  = 1.8,
    ),
}

STANDARD_CMAPS = [
    "jet", "viridis", "plasma", "inferno", "magma", "cividis",
    "hot", "cool", "gray", "bone", "copper", "spring", "summer",
    "autumn", "winter", "rainbow", "turbo", "hsv",
    "Spectral", "RdYlBu", "RdBu", "coolwarm", "seismic",
]

def colormap_to_stops(cmap_name: str, n: int = 9) -> list:
    cmap = cm.get_cmap(cmap_name)
    return [(round(float(p), 4), to_hex(cmap(p)[:3]))
            for p in np.linspace(0.0, 1.0, n)]

def stops_to_cmap(stops, name: str = "custom") -> LinearSegmentedColormap:
    stops = sorted(stops, key=lambda s: s[0])
    colors = [(p, to_rgb(c)) for p, c in stops]
    r = [(p, rgb[0], rgb[0]) for p, rgb in colors]
    g = [(p, rgb[1], rgb[1]) for p, rgb in colors]
    b = [(p, rgb[2], rgb[2]) for p, rgb in colors]
    return LinearSegmentedColormap(name, {"red": r, "green": g, "blue": b}, N=512)


class ColormapEditorDialog(wx.Dialog):

    _MARKER_SIZE   = 180   # scatter marker size in points² — always a circle
    _HIT_RADIUS_PT = 14    # hit-test radius in screen points

    def __init__(self, parent=None, initial_cmap=_SENTINEL,
                 initial_stops: int = 9, title: str = "Colormap Editor",
                 theme: str = "light", preset_stops=None, preset_name: str = ""):
        """
        preset_stops : list of (pos, hex) pairs — if supplied the dialog opens
                       showing exactly these stops under preset_name, ignoring
                       initial_cmap entirely.  This is the correct way to
                       re-open the editor after a custom/loaded colormap has
                       already been applied so the user sees their current state.
        preset_name  : name string shown in the Name field when preset_stops
                       is used.
        """
        super().__init__(parent, title=title, size=(680, 560),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self._theme_name = theme if theme in THEMES else "light"
        self._T          = THEMES[self._theme_name]   # shorthand
        self._stops      = []      # [[pos, hex], ...]
        self._drag_idx   = None
        self._spin_popup = None    # floating SpinCtrlDouble while editing
        self._current_stop = None
        self._active_idx = None    # single-click selected stop
        self._history    = deque(maxlen=10)   # undo stack (snapshots of _stops)
        self._build_ui()
        self.SetMinSize((480, 400))

        # Priority 1: caller supplied explicit stops (e.g. reopen after load)
        if preset_stops is not None and len(preset_stops) >= 2:
            self._load_stops([(float(p), str(c)) for p, c in preset_stops])
            self._name_ctrl.SetValue(preset_name or "custom")
            self._n_spin.SetValue(len(preset_stops))
            return

        # Priority 2: no cmap specified → try the saved default file
        if initial_cmap is _SENTINEL or initial_cmap is None:
            try:
                data = json.loads(DEFAULT_CMAP_FILE.read_text())
                self._name_ctrl.SetValue(data["name"])
                self._n_spin.SetValue(len(data["stops"]))
                self._load_stops([(float(p), str(c)) for p, c in data["stops"]])
                self._default_label.SetLabel(f"Default: {data['name']}")
                return
            except Exception:
                initial_cmap = "jet"

        # Priority 3: seed from a named standard colormap
        sel = initial_cmap if initial_cmap in STANDARD_CMAPS else "jet"
        self._cmap_choice.SetStringSelection(sel)
        self._n_spin.SetValue(initial_stops)
        self._load_stops(colormap_to_stops(sel, initial_stops))
        self._name_ctrl.SetValue(f"{sel}_custom")

    # ── public API ────────────────────────────────────────────────────────────

    def get_colormap(self):
        stops = self.get_stops()
        return stops_to_cmap(stops, self.get_colormap_name()) if len(stops) >= 2 else None

    def get_colormap_name(self) -> str:
        return self._name_ctrl.GetValue().strip() or "custom"

    def get_stops(self) -> list:
        return [(p, h) for p, h in self._stops]

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        self._panel = wx.Panel(self)
        panel = self._panel
        main  = wx.BoxSizer(wx.VERTICAL)

        # load bar
        load_box = wx.StaticBoxSizer(wx.HORIZONTAL, panel, "Load standard colormap")
        self._cmap_choice = wx.Choice(panel, choices=STANDARD_CMAPS, style=wx.BORDER_NONE)
        self._cmap_choice.SetStringSelection("jet")
        self._n_spin = wx.SpinCtrl(panel, value="9", min=3, max=24, size=(52, -1))
        self.load_btn = wx.Button(panel, label="Load", style=wx.BORDER_NONE)
        self.load_btn.Bind(wx.EVT_BUTTON, lambda e: self._on_load())
        load_box.Add(self._cmap_choice, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        load_box.Add(wx.StaticText(panel, label="Stops:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)
        load_box.Add(self._n_spin, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        load_box.Add(self.load_btn,    0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        main.Add(load_box, 0, wx.EXPAND | wx.ALL, 5)

        # toolbar
        tb = wx.BoxSizer(wx.HORIZONTAL)
        self.add_btn    = wx.Button(panel, label="+ Add stop", style=wx.BORDER_NONE)
        self.del_btn    = wx.Button(panel, label="− Delete", style=wx.BORDER_NONE)
        self.resort_btn = wx.Button(panel, label="↕ Re-sort", style=wx.BORDER_NONE)
        self._undo_btn  = wx.Button(panel, label="↩ Undo", style=wx.BORDER_NONE)
        self._theme_btn = wx.Button(panel, label="", style=wx.BORDER_NONE)   # label set by _apply_theme
        self.add_btn.SetToolTip("Insert a stop after the selected stop (click a circle to select)")
        self.del_btn.SetToolTip("Delete the selected stop (click a circle to select)")
        self.resort_btn.SetToolTip("Re-order interior stops by position value")
        self._undo_btn.SetToolTip("Undo last change (Ctrl+Z)")
        self._theme_btn.SetToolTip("Switch between light and dark theme")
        self._undo_btn.Enable(False)
        self.add_btn.Bind(wx.EVT_BUTTON,         lambda e: self._on_add_stop())
        self.del_btn.Bind(wx.EVT_BUTTON,         lambda e: self._on_delete_stop())
        self.resort_btn.Bind(wx.EVT_BUTTON,      lambda e: self._on_resort())
        self._undo_btn.Bind(wx.EVT_BUTTON,  lambda e: self._on_undo())
        self._theme_btn.Bind(wx.EVT_BUTTON, lambda e: self._toggle_theme())
        tb.Add(self.add_btn,          0, wx.RIGHT, 4)
        tb.Add(self.del_btn,          0, wx.RIGHT, 4)
        tb.Add(self.resort_btn,       0, wx.RIGHT, 12)
        tb.Add(self._undo_btn,   0, wx.RIGHT, 12)
        tb.AddStretchSpacer()
        tb.Add(self._theme_btn,  0)
        main.Add(tb, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM | wx.EXPAND, 5)

        # Ctrl-Z / Cmd-Z accelerator
        undo_id = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, lambda e: self._on_undo(), id=undo_id)
        self.SetAcceleratorTable(wx.AcceleratorTable([
            (wx.ACCEL_CTRL, ord('Z'), undo_id),
        ]))

        # matplotlib figure
        self._figure = Figure(facecolor="none")
        self._ax     = self._figure.add_axes([0.08, 0.06, 0.74, 0.90])
        self._ax_cb  = self._figure.add_axes([0.87, 0.06, 0.06, 0.90])

        self._ax.set_ylim(-0.05, 1.05)
        self._ax.set_ylabel("Colormap position", fontsize=9)
        self._ax.tick_params(labelbottom=False, bottom=False)
        self._ax_cb.set_axis_off()

        # persistent artists
        self._line, = self._ax.plot([], [], '-', linewidth=1.5, zorder=1)
        self._scat  = self._ax.scatter([], [], s=self._MARKER_SIZE,
                                       zorder=3, linewidths=1.8,
                                       edgecolors='#333333')
        self._cb_img = self._ax_cb.imshow(
            np.linspace(0, 1, 256).reshape(-1, 1),
            aspect='auto', origin='lower',
            cmap=cm.get_cmap("jet"), extent=[0, 1, 0, 1])

        self._canvas = FigureCanvas(panel, -1, self._figure)
        self._canvas.SetMinSize((-1, 180))
        self._canvas.mpl_connect("button_press_event",   self._on_mpl_press)
        self._canvas.mpl_connect("motion_notify_event",  self._on_mpl_motion)
        self._canvas.mpl_connect("button_release_event", self._on_mpl_release)
        main.Add(self._canvas, 1, wx.EXPAND | wx.ALL, 0)

        # bottom bar
        bot = wx.BoxSizer(wx.HORIZONTAL)
        bot.Add(wx.StaticText(panel, label="Name:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self._name_ctrl = wx.TextCtrl(panel, value="my_colormap", size=(140, -1))
        self._name_ctrl.Bind(wx.EVT_TEXT, lambda e: self._refresh_colorbar())
        bot.Add(self._name_ctrl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 8)
        save_btn = wx.Button(panel, label="★ Save as Default", style=wx.BORDER_NONE)
        save_btn.Bind(wx.EVT_BUTTON, lambda e: self._on_save_default())
        bot.Add(save_btn, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        self._default_label = wx.StaticText(panel, label="")
        bot.Add(self._default_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        bot.AddStretchSpacer()
        export_btn = wx.Button(panel, label="Export code…", style=wx.BORDER_NONE)
        export_btn.Bind(wx.EVT_BUTTON, lambda e: self._on_export())
        ok_btn = wx.Button(panel, id=wx.ID_OK, label="OK", style=wx.BORDER_NONE)
        ok_btn.Bind(wx.EVT_BUTTON, lambda e: self._on_ok())
        ok_btn.SetDefault()
        bot.Add(export_btn, 0, wx.RIGHT, 6)
        bot.Add(ok_btn,     0, wx.RIGHT, 6)
        cancel_btn = wx.Button(panel, id=wx.ID_CANCEL, label="Cancel", style=wx.BORDER_NONE)
        cancel_btn.Bind(wx.EVT_BUTTON, lambda e: self._on_cancel())
        bot.Add(cancel_btn)
        main.Add(bot, 0, wx.EXPAND | wx.ALL, 6)

        panel.SetSizer(main)
        root = wx.BoxSizer(wx.VERTICAL)
        root.Add(panel, 1, wx.EXPAND)
        self.SetSizer(root)

        try:
            data = json.loads(DEFAULT_CMAP_FILE.read_text())
            self._default_label.SetLabel(f"Default: {data['name']}")
        except Exception:
            pass

        # Apply initial theme colours
        self._apply_theme()

    # ── undo history ──────────────────────────────────────────────────────────

    def _push_history(self):
        """Snapshot current stops onto the undo stack before a mutation."""
        self._history.append(copy.deepcopy(self._stops))
        self._update_undo_btn()

    def _on_undo(self):
        if not self._history:
            return
        self._stops = self._history.pop()
        # Clamp active_idx to valid range after undo
        if self._active_idx is not None:
            self._active_idx = min(self._active_idx, len(self._stops) - 1)
        self._update_undo_btn()
        self._rebuild_graph()

    def _update_undo_btn(self):
        self._undo_btn.Enable(bool(self._history))

    # ── theming ───────────────────────────────────────────────────────────────

    def _toggle_theme(self):
        self._theme_name = "dark" if self._theme_name == "light" else "light"
        self._T = THEMES[self._theme_name]
        self._apply_theme()

    def _apply_theme(self):
        T = self._T
        # Update toggle button label
        self._theme_btn.SetLabel("☀ Light" if self._theme_name == "dark" else "🌙 Dark")

        # wx widget colours — walk all children of the panel
        self.SetBackgroundColour(T["dlg_bg"])
        self._panel.SetBackgroundColour(T["panel_bg"])
        for w in self._panel.GetChildren():
            if isinstance(w, (wx.StaticText,)):
                w.SetForegroundColour(T["text_fg"])
                w.SetBackgroundColour(T["panel_bg"])
                
            elif isinstance(w, (wx.TextCtrl, wx.SpinCtrl)):
                w.SetForegroundColour(T["text_fg"])
                w.SetBackgroundColour(T["panel_bg"])
            elif isinstance(w, (wx.Button)):
                w.SetForegroundColour(T["text_fg"])
                w.SetBackgroundColour(T["btn_bg"])      
            elif isinstance(w, (wx.Choice)):
                w.SetForegroundColour(T["text_fg"])
                w.SetBackgroundColour(T["btn_bg"])  
                #w.BackgroundColour = T["btn_bg"] 
                
        # Override the "saved" label colour specifically
        self._default_label.SetForegroundColour(T["label_fg"])
        self._panel.Refresh()
        self.Refresh()

        # matplotlib figure colours
        self._figure.set_facecolor(T["fig_bg"])
        self._ax.set_facecolor(T["ax_bg"])
        self._ax.tick_params(colors=T["ax_fg"])
        self._ax.yaxis.label.set_color(T["ax_fg"])
        for spine in self._ax.spines.values():
            spine.set_edgecolor(T["ax_fg"])
        self._ax.grid(True, axis='y', color=T["grid_col"], linewidth=0.7)
        self._line.set_color(T["line_col"])

        # Redraw scatter with updated edge logic (uses _T)
        self._rebuild_graph()

    # ── stop management ───────────────────────────────────────────────────────

    def _load_stops(self, stops: list):
        self._push_history()
        self._stops = [[float(p), str(h)] for p, h in stops]
        self._stops[0][0]  = 0.0
        self._stops[-1][0] = 1.0
        self._stops.sort(key=lambda s: s[0])
        self._drag_idx = None
        self._active_idx = None
        self._dismiss_spin_popup()
        self._rebuild_graph()

    def _rebuild_graph(self):
        """Full scatter rebuild — called when stop count or colours change."""
        T   = self._T
        n   = len(self._stops)
        xs  = list(range(n))
        ys  = [s[0] for s in self._stops]
        fcs = [s[1] for s in self._stops]
        ecs = []
        lws = []
        for i, hx in enumerate(fcs):
            r, g, b = [int(hx.lstrip('#')[k:k+2], 16) / 255 for k in (0, 2, 4)]
            if i == self._active_idx:
                ecs.append(T["active_edge"])
                lws.append(T["active_lw"])
            else:
                lum = 0.299*r + 0.587*g + 0.114*b
                # In dark theme: dark stops get a light border and vice-versa
                if self._theme_name == "dark":
                    ecs.append(T["edge_light"])
                else:
                    ecs.append(T["edge_dark"])
                lws.append(T["inactive_lw"])

        self._scat.set_offsets(np.column_stack([xs, ys]))
        self._scat.set_facecolors(fcs)
        self._scat.set_edgecolors(ecs)
        self._scat.set_linewidths(lws)
        self._line.set_data(xs, ys)
        self._ax.set_xlim(-0.5, n - 0.5)
        self._ax.set_xticks(xs)
        self._ax.set_xticklabels([])
        self._refresh_colorbar()
        

    def _fast_update_marker(self, idx: int, new_y: float):
        """Move one marker and the line during drag — minimal redraw."""
        offsets = self._scat.get_offsets().copy()
        offsets[idx, 1] = new_y
        self._scat.set_offsets(offsets)
        self._line.set_ydata([s[0] for s in self._stops])
        try:
            cmap = stops_to_cmap(self._stops, self.get_colormap_name())
            self._cb_img.set_cmap(cmap)
        except Exception:
            pass
        self._canvas.draw()

    def _refresh_colorbar(self):
        if len(self._stops) < 2:
            return
        try:
            cmap = stops_to_cmap(self._stops, self.get_colormap_name())
            self._cb_img.set_cmap(cmap)
            self._canvas.draw()
        except Exception as e:
            print(f"[ColormapEditor] {e}")

    # ── hit testing (screen-space, matches marker size) ───────────────────────

    def _hit_index(self, event) -> int | None:
        if event.inaxes is not self._ax or event.x is None:
            return None
        # transform each stop's data coords → display coords, compare in pixels
        n  = len(self._stops)
        xs = list(range(n))
        ax = self._ax
        for i, (x_d, (pos, _)) in enumerate(zip(xs, self._stops)):
            x_px, y_px = ax.transData.transform((x_d, pos))
            dist = ((event.x - x_px)**2 + (event.y - y_px)**2) ** 0.5
            if dist <= self._HIT_RADIUS_PT:
                return i
        return None

    # ── matplotlib mouse events ───────────────────────────────────────────────

    def _on_mpl_press(self, event):
        if event.inaxes is not self._ax:
            return
        self._dismiss_spin_popup()
        idx = self._hit_index(event)
        if idx is None:
            # click on empty area clears selection
            if self._active_idx is not None:
                self._active_idx = None
                self._rebuild_graph()
            return
        if event.button == 3:                      # right-click → colour
            self._pick_color(idx)
        elif event.dblclick:                        # double-click → spin popup
            self._active_idx = idx
            self._rebuild_graph()
            self._show_spin_popup(idx)
        else:                                       # single left → activate + maybe drag
            self._active_idx = idx
            self._rebuild_graph()
            if idx not in (0, len(self._stops) - 1):
                self._push_history()
                self._drag_idx = idx

    def _on_mpl_motion(self, event):
        if self._drag_idx is None or event.inaxes is not self._ax or event.ydata is None:
            return
        idx = self._drag_idx
        y   = float(np.clip(event.ydata, 0.0, 1.0))
        self._stops[idx][0] = round(y, 4)
        self._fast_update_marker(idx, y)

    def _on_mpl_release(self, event):
        if self._drag_idx is not None:
            # If the marker wasn't actually moved, discard the snapshot we pushed
            if (self._history and
                    self._history[-1][self._drag_idx][0] == self._stops[self._drag_idx][0]):
                self._history.pop()
                self._update_undo_btn()
        self._drag_idx = None

    # ── colour picker ─────────────────────────────────────────────────────────

    def _pick_color(self, idx: int):
        hx = self._stops[idx][1]
        try:
            r, g, b = [int(hx.lstrip('#')[k:k+2], 16) for k in (0, 2, 4)]
            init_col = wx.Colour(r, g, b)
        except Exception:
            init_col = wx.Colour(128, 128, 128)
        cdlg = wx.ColourDialog(self)
        cdlg.GetColourData().SetColour(init_col)
        cdlg.GetColourData().SetChooseFull(True)
        if cdlg.ShowModal() == wx.ID_OK:
            col = cdlg.GetColourData().GetColour()
            self._push_history()
            self._stops[idx][1] = "#{:02x}{:02x}{:02x}".format(
                col.Red(), col.Green(), col.Blue())
            self._rebuild_graph()
        cdlg.Destroy()

    # ── floating spin popup for exact value entry ─────────────────────────────

    def _show_spin_popup(self, idx: int):
        if idx in (0, len(self._stops) - 1):
            return   # anchors are fixed
        self._dismiss_spin_popup()
        self._current_stop = idx
        # find screen position of the marker in canvas pixel coords
        ax   = self._ax
        x_d  = idx
        y_d  = self._stops[idx][0]
        x_px, y_px = ax.transData.transform((x_d, y_d))

        # canvas is the wx widget; convert matplotlib display coords
        # (origin bottom-left) to wx coords (origin top-left)
        canvas_h = self._canvas.GetSize()[1]
        wx_x = int(x_px)
        wx_y = int(canvas_h - y_px)

        spin = wx.SpinCtrlDouble(self._canvas, value=f"{y_d:.4f}",
                                 min=0.0, max=1.0, inc=0.001,
                                 size=(90, -1),
                                 style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        spin.SetDigits(4)
        spin.SetValue(y_d)
        # centre the popup on the marker
        spin_w, spin_h = spin.GetSize()
        spin.SetPosition(wx.Point(wx_x - spin_w // 2, wx_y - spin_h // 2))

        spin.Bind(wx.EVT_TEXT_ENTER, self.spin_commit)
        spin.Bind(wx.EVT_KILL_FOCUS, self.spin_commit)
        spin.SetFocus()
        self._spin_popup = spin
        
    def spin_commit(self, _evt=None):
        try:
            val = np.clip(float(self._spin_popup.GetTextValue()), 0, 1)
            print(val)
        except Exception:
            self._dismiss_spin_popup()
            self._rebuild_graph()
            return
        idx = self._current_stop
        if idx is None: return
        self._push_history()
        self._stops[idx][0] = round(val, 4)
        self._dismiss_spin_popup()
        self._rebuild_graph()
    def _dismiss_spin_popup(self, _evt=None):
        if self._spin_popup is not None:
            try:
                self._spin_popup.Destroy()
            except Exception:
                pass
            self._spin_popup = None

    # ── toolbar ───────────────────────────────────────────────────────────────

    def _on_load(self):
        name = self._cmap_choice.GetStringSelection()
        self._load_stops(colormap_to_stops(name, self._n_spin.GetValue()))
        self._name_ctrl.SetValue(f"{name}_custom")

    def _on_add_stop(self):
        if len(self._stops) >= 24:
            return
        self._push_history()
        # Determine insertion point: after active stop, or after largest gap
        if self._active_idx is not None:
            idx = self._active_idx
        else:
            gaps = [(abs(self._stops[i+1][0] - self._stops[i][0]), i)
                    for i in range(len(self._stops) - 1)]
            _, idx = max(gaps)

        # Insert after idx; if idx is the last stop, insert before it (midpoint with prev)
        if idx < len(self._stops) - 1:
            new_pos = round((self._stops[idx][0] + self._stops[idx+1][0]) / 2.0, 4)
            insert_at = idx + 1
        else:
            # active stop is the last one — insert between it and its predecessor
            new_pos = round((self._stops[idx-1][0] + self._stops[idx][0]) / 2.0, 4)
            insert_at = idx  # inserts before the last stop, i.e. after second-to-last
        try:
            hex_col = to_hex(stops_to_cmap(self._stops)(new_pos)[:3])
        except Exception:
            hex_col = "#808080"
        self._stops.insert(insert_at, [new_pos, hex_col])
        self._active_idx = insert_at   # select the newly added stop
        self._rebuild_graph()

    def _on_delete_stop(self):
        n = len(self._stops)
        if n <= 3:
            wx.MessageBox("Need at least 3 stops (first and last are locked).",
                          "Cannot delete", wx.ICON_WARNING)
            return
        idx = self._active_idx
        if idx is None:
            wx.MessageBox("Click a stop circle to select it first, then delete.",
                          "No stop selected", wx.ICON_INFORMATION)
            return
        if idx in (0, n - 1):
            wx.MessageBox("The first and last stops are locked and cannot be deleted.",
                          "Cannot delete", wx.ICON_WARNING)
            return
        self._push_history()
        self._stops.pop(idx)
        # Keep selection on a valid neighbour
        self._active_idx = min(idx, len(self._stops) - 1)
        self._rebuild_graph()

    def _on_resort(self):
        self._push_history()
        first, last = self._stops[0], self._stops[-1]
        interior = sorted(self._stops[1:-1], key=lambda s: s[0])
        self._stops = [first] + interior + [last]
        self._rebuild_graph()

    # ── persistence / export ──────────────────────────────────────────────────

    def _on_save_default(self):
        stops = self.get_stops()
        name  = self.get_colormap_name()
        try:
            DEFAULT_CMAP_FILE.write_text(json.dumps({"name": name, "stops": stops}, indent=2))
        except Exception as e:
            wx.MessageBox(f"Could not save:\n{e}", "Error", wx.ICON_ERROR)
            return
        self._default_label.SetLabel(f"Default: {name}")
        self._default_label.GetParent().Layout()
        wx.MessageBox(f'"{name}" saved as default.', "Saved", wx.ICON_INFORMATION)

    def _on_ok(self):
        if self._spin_popup is not None:
            return
        self.EndModal(wx.ID_OK)

    def _on_cancel(self):
        if self._spin_popup is not None:
            self._dismiss_spin_popup()
            return
        self.EndModal(wx.ID_CANCEL)

    def _on_export(self):
        stops = self.get_stops()
        name  = self.get_colormap_name()
        code  = "\n".join([
            "from matplotlib.colors import LinearSegmentedColormap, to_rgb",
            "",
            f"# Colormap: {name}",
            "stops = [",
            *[f"    ({p!r}, {c!r})," for p, c in stops],
            "]",
            "colors = [(p, to_rgb(c)) for p, c in stops]",
            "r = [(p, rgb[0], rgb[0]) for p, rgb in colors]",
            "g = [(p, rgb[1], rgb[1]) for p, rgb in colors]",
            "b = [(p, rgb[2], rgb[2]) for p, rgb in colors]",
            f'{name} = LinearSegmentedColormap("{name}",'
            f' {{"red":r,"green":g,"blue":b}}, N=512)',
        ])
        dlg = wx.Dialog(self, title="Exported Python code", size=(560, 380),
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        vsz = wx.BoxSizer(wx.VERTICAL)
        tc  = wx.TextCtrl(dlg, value=code,
                          style=wx.TE_MULTILINE | wx.TE_READONLY | wx.HSCROLL | wx.TE_RICH2)
        tc.SetFont(wx.Font(9, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        btn_row  = wx.BoxSizer(wx.HORIZONTAL)
        copy_btn = wx.Button(dlg, label="Copy to clipboard")
        close_b  = wx.Button(dlg, id=wx.ID_CLOSE, label="Close")
        close_b.Bind(wx.EVT_BUTTON, lambda _: dlg.EndModal(wx.ID_CLOSE))
        def _copy(_):
            if wx.TheClipboard.Open():
                wx.TheClipboard.SetData(wx.TextDataObject(code))
                wx.TheClipboard.Close()
        copy_btn.Bind(wx.EVT_BUTTON, _copy)
        btn_row.AddStretchSpacer()
        btn_row.Add(copy_btn, 0, wx.RIGHT, 6)
        btn_row.Add(close_b)
        vsz.Add(tc,      1, wx.EXPAND | wx.ALL, 8)
        vsz.Add(btn_row, 0, wx.EXPAND | wx.ALL, 8)
        dlg.SetSizer(vsz)
        dlg.ShowModal()
        dlg.Destroy()


# ── standalone test ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt

    # Accept optional --theme light|dark on the command line
    _theme = "light"
    if "--theme" in sys.argv:
        _i = sys.argv.index("--theme")
        if _i + 1 < len(sys.argv):
            _theme = sys.argv[_i + 1]

    app = wx.App(False)
    with ColormapEditorDialog(theme=_theme) as dlg:
        if dlg.ShowModal() == wx.ID_OK:
            cmap  = dlg.get_colormap()
            name  = dlg.get_colormap_name()
            stops = dlg.get_stops()
            print(f"name : {name}")
            print(f"stops: {stops}")
            x = np.linspace(-np.pi, np.pi, 300)
            X, Y = np.meshgrid(x, x)
            plt.figure(figsize=(6, 5))
            plt.imshow(np.sin(X) * np.cos(Y), cmap=cmap, origin="lower")
            plt.colorbar(label=name)
            plt.title(f"Custom colormap: {name}")
            plt.tight_layout()
            plt.show()
        else:
            print("Cancelled.")
    app.MainLoop()