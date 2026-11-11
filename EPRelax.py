"""Regularized exponential analysis of Bruker relaxation data, in the style of pyCWEPR.

Loads BES3T (.DSC/.DTA) or ESP (.PAR/.SPC) files holding a relaxation measurement -- inversion
recovery, saturation recovery, a 2-pulse or 3-pulse echo decay, anything that is a sum of decaying
exponentials -- and decomposes it into a distribution of relaxation times, the way CONTIN does for
dynamic light scattering.

The model is

    S(t) = c  +  s * sum_j a_j exp(-t/T_j),     a_j >= 0,  c free,  s = +1 or -1

with the T_j a fixed, logarithmically spaced grid, so only the amplitudes a_j and the offset c are
unknown and the problem is linear. s = -1 (the Recovery switch) turns the same solver into
inversion or saturation recovery, where the exponentials come up to the offset from below. The
decomposition is ill-posed, so it is regularized:

    minimize  || A x - y ||^2  +  alpha^2 || L x ||^2      subject to  a_j >= 0

L is either the second difference along the log-T grid (CONTIN's choice: it penalises curvature, so
the answer is the smoothest distribution still consistent with the data, with the distribution taken
as zero outside the grid so that its two end points are priced like the middle) or the identity
(plain Tikhonov, which penalises amplitude instead). Non-negativity is imposed exactly by solving the
augmented system with non-negative least squares; the offset rides along as a pair of non-negative
columns so that it stays free in sign. alpha can be dialled by hand, taken from the corner of the
L-curve, from the minimum of the generalized cross-validation score, or from CONTIN's F-test on the
residual sum of squares.

2-pulse ESEEM carries deep nuclear modulation that no sum of decaying exponentials can represent.
The trace is decomposed exactly as recorded by default, with the residual plotted underneath so the
modulation is visible rather than hidden; a boxcar one modulation period wide is available to take
it out first, which leaves an exponential exponential (convolution with a boxcar only rescales it).

A discrete fit sits alongside: you say how many exponentials, the starting relaxation times are read
off the maxima of the regularized distribution, and a non-linear least-squares fit refines them.

pip install wxPython numpy matplotlib scipy
"""

import os
import ast
import base64
import fnmatch
import warnings
import traceback
import xml.etree.ElementTree as ET
from xml.dom import minidom

import wx
import wx.lib.mixins.listctrl as listmix
import wx.lib.scrolledpanel as scrolled

import numpy as np
from scipy.optimize import nnls, curve_fit, OptimizeWarning
from scipy.ndimage import uniform_filter1d
from scipy.stats import f as fdist, t as tdist

from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

from theme import theme
from theme.theme import PillButton, PillCheckBox, PillTabBar, SegmentedPill
from readers import brukerread as Mybr
from spin_system.sys_functions import extract_refs_and_freevars

VERSION = 0.1
HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = 33                               # points on the alpha sweep behind the L-curve, GCV and the F-test
ROWS = ('data', 'resid', 'dist', 'lcurve')
ROW_HEIGHT = {'data': 2.2, 'resid': 1.0, 'dist': 2.0, 'lcurve': 1.6}
# the ladder a loaded time axis is rescaled along, so the numbers on screen stay between 1 and 1000
TIME_UNITS = (('ps', 1e-12), ('ns', 1e-9), ('µs', 1e-6), ('ms', 1e-3), ('s', 1.0))
TWO_SIDED_95 = 1.959963985           # the normal quantile the Gaussian peak width is turned into a range with
# what a hand-typed fit equation may call: the suite's own names, over arrays rather than scalars
FIT_FUNCS = {'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log, 'log10': np.log10,
             'sin': np.sin, 'cos': np.cos, 'tan': np.tan, 'asin': np.arcsin, 'acos': np.arccos,
             'atan': np.arctan, 'atan2': np.arctan2, 'sinh': np.sinh, 'cosh': np.cosh,
             'tanh': np.tanh, 'pi': np.pi, 'e': np.e, 'abs': np.abs, 'min': np.minimum,
             'max': np.maximum, 'pow': np.power}


class DictPopup(wx.Dialog):
    """Key/value view of a file header, shown by the browser's Header button without loading."""
    def __init__(self, parent, data_dict, title="Header"):
        super().__init__(parent, title=title, size=(560, 640),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.list_ctrl = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.BORDER_NONE)
        self.list_ctrl.InsertColumn(0, "Key", width=220)
        self.list_ctrl.InsertColumn(1, "Value", width=300)
        for key, value in data_dict.items():
            index = self.list_ctrl.InsertItem(self.list_ctrl.GetItemCount(), str(key))
            self.list_ctrl.SetItem(index, 1, str(value))
        vbox.Add(self.list_ctrl, 1, wx.EXPAND | wx.ALL, 8)
        vbox.Add(wx.Button(panel, wx.ID_OK, "OK"), 0, wx.ALIGN_CENTER | wx.BOTTOM, 8)
        panel.SetSizer(vbox)
        theme.apply_theme_to_window(self, theme.get_theme())


# --------------------------------------------------------------------------- #
# File Browser Panel
# --------------------------------------------------------------------------- #
class FileBrowserPanel(wx.Panel, listmix.ColumnSorterMixin):
    """Folder listing with the Bruker pairs shown once, as in pyCWEPR; loaded files get a marker."""

    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.parent = mainWindow
        self.ini_file = os.path.join(wx.StandardPaths.Get().GetUserConfigDir(), "EPRelax.ini")
        self.path = os.getcwd()
        self.filters = {"EPR": ('.dsc', '.par', '.xml'), "Bruker": ('.dsc', '.par'),
                        "Session": ('.xml',), "All": None}
        self.current_filter = "EPR"

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.lbl_path = wx.TextCtrl(self, value="", style=wx.TE_PROCESS_ENTER | wx.TE_RIGHT | wx.BORDER_NONE)
        self.lbl_path.Bind(wx.EVT_TEXT_ENTER, self.on_path_entered)
        self.lbl_path.Bind(wx.EVT_KILL_FOCUS, self.on_path_entered)
        vbox.Add(self.lbl_path, 0, wx.EXPAND)

        grid = wx.FlexGridSizer(cols=2, vgap=4, hgap=5)
        grid.AddGrowableCol(1)
        grid.Add(wx.StaticText(self, label="Filter:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.choice_filter = wx.ComboBox(self, choices=list(self.filters), value=self.current_filter,
                                         style=wx.CB_READONLY | wx.BORDER_NONE)
        self.choice_filter.SetToolTip("Which file extensions the list shows.")
        self.choice_filter.Bind(wx.EVT_COMBOBOX, self.on_filter_changed)
        grid.Add(self.choice_filter, 1, wx.EXPAND)
        grid.Add(wx.StaticText(self, label="Title:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.txt_keyword = wx.TextCtrl(self, value="*", style=wx.TE_PROCESS_ENTER | wx.BORDER_NONE)
        self.txt_keyword.SetToolTip("Keep only the files whose stored title matches this, on top of the\n"
                                    "extension filter. Case does not matter, and folders always stay.\n"
                                    "A plain word matches anywhere in the title, so 'nitroxide' finds\n"
                                    "'Nitroxide 20 K'. Use * ? [ ] to say exactly where: 'T1*' only at the\n"
                                    "start, '*20 K' only at the end. * on its own, the default, keeps\n"
                                    "everything -- including the files whose header carries no title.\n"
                                    "Press Enter or click away to apply.")
        self.txt_keyword.Bind(wx.EVT_TEXT_ENTER, self.on_filter_changed)
        self.txt_keyword.Bind(wx.EVT_KILL_FOCUS, self.on_filter_changed)
        grid.Add(self.txt_keyword, 1, wx.EXPAND)
        vbox.Add(grid, 0, wx.EXPAND | wx.ALL, 5)

        self.list = wx.ListCtrl(self, style=wx.LC_REPORT | wx.LC_VRULES | wx.LC_HRULES | wx.BORDER_NONE)
        self.list.InsertColumn(0, "Name", width=220)
        self.list.InsertColumn(1, "Title", width=240)
        vbox.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)

        img_list = wx.ImageList(16, 16)
        icon = lambda name: wx.Bitmap(os.path.join(HERE, 'icons', name), wx.BITMAP_TYPE_ICO)
        self.folder_img_id = img_list.Add(icon('Folder.ico'))
        self.file_img_id = img_list.Add(icon('New document.ico'))
        self.data_img_id = img_list.Add(icon('List.ico'))
        self.xml_img_id = img_list.Add(icon('Script.ico'))
        self.loaded_img_id = img_list.Add(icon('Green bookmark.ico'))
        self.sort_up_img_id = img_list.Add(wx.ArtProvider.GetBitmap(wx.ART_GO_UP, wx.ART_OTHER, (16, 16)))
        self.sort_down_img_id = img_list.Add(wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN, wx.ART_OTHER, (16, 16)))
        self.list.AssignImageList(img_list, wx.IMAGE_LIST_SMALL)

        btns = wx.BoxSizer(wx.VERTICAL)
        self.btn_load = PillButton(self, "Load", color_key="pill_load", draw_badge=theme.draw_squiggle_badge)
        self.btn_load.SetToolTip("Load the selected Bruker file(s). Ctrl- or Shift-click selects several.")
        self.btn_load.Bind(wx.EVT_BUTTON, self.on_load_clicked)
        btns.Add(self.btn_load, 0, wx.EXPAND | wx.ALL, 2)
        self.btn_header = PillButton(self, "Header", color_key="pill_calm", draw_badge=theme.draw_textlines_badge)
        self.btn_header.SetToolTip("Show the parameters stored in the selected file without loading it.")
        self.btn_header.Bind(wx.EVT_BUTTON, self.on_header_clicked)
        btns.Add(self.btn_header, 0, wx.EXPAND | wx.ALL, 2)
        vbox.Add(btns, 0, wx.EXPAND)
        self.SetSizer(vbox)

        listmix.ColumnSorterMixin.__init__(self, 2)
        self.items = []          # (name, full_path, is_dir); row data id = index here
        self.itemDataMap = {}

        self.parent.tb_change_folder.Bind(wx.EVT_BUTTON, self.on_change_folder)
        self.parent.tb_up.Bind(wx.EVT_BUTTON, self.on_up_clicked)
        self.parent.tb_load_all.Bind(wx.EVT_BUTTON, self.on_load_all)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_item_activated)
        self.Bind(wx.EVT_WINDOW_DESTROY, lambda e: self.save_path())

        try:
            with open(self.ini_file, "r", encoding="utf-8") as f:
                stored = f.readline().strip()
            if os.path.isdir(stored):
                self.path = stored
        except OSError as fault:
            print(f"[EPRelax] Could not read the last folder from {self.ini_file}: {fault}")
        self.lbl_path.SetValue(self.path)
        self.lbl_path.SetInsertionPointEnd()
        self.refresh_file_list()

    def save_path(self):
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.path + "\n")
        except OSError as fault:
            print(f"[EPRelax] Could not save the last folder to {self.ini_file}: {fault}")

    def refresh_file_list(self, select_name=None):
        """Rebuild the list; keeps `select_name` (a name or a list of names) selected if still there."""
        try:
            raw = os.listdir(self.path)
        except OSError as fault:
            self.parent.statusbar.SetStatusText(f"Could not list folder: {fault}", 0)
            raw = []
        exts = self.filters[self.current_filter]
        dirs = sorted([n for n in raw if os.path.isdir(os.path.join(self.path, n))], key=str.lower)
        files = sorted([n for n in raw if not os.path.isdir(os.path.join(self.path, n))
                        and (exts is None or n.lower().endswith(exts))], key=str.lower)
        loaded = {os.path.splitext(dd['fullpath'])[0].lower() for dd in self.parent.Data}
        # a plain word matches anywhere in the title; * ? [ ] say where
        pattern = self.txt_keyword.GetValue().strip().lower() or '*'
        if not any(c in pattern for c in '*?['):
            pattern = f'*{pattern}*'

        self.list.DeleteAllItems()
        self.items = []
        self.itemDataMap = {}
        hidden = 0
        for name in dirs + files:
            full = os.path.join(self.path, name)
            is_dir = name in dirs
            title = ''
            low = name.lower()
            if is_dir:
                image_id = self.folder_img_id
            elif low.endswith(('.dsc', '.par')):
                image_id = self.data_img_id
                try:                             # the title Xepr wrote, so the list reads like the experiment log
                    with open(full, 'r', encoding='latin-1') as f:
                        for line in f:
                            if line.startswith(('TITL', 'JCO')):
                                title = line.split("'")[1] if "'" in line else line.split(None, 1)[-1].strip()
                                break
                except (OSError, IndexError) as fault:
                    print(f"[EPRelax] Could not read the title from '{full}': {fault}")
            elif low.endswith('.xml'):
                image_id = self.xml_img_id
                try:                             # an EPRelax session says so in its first few lines
                    with open(full, 'r', encoding='utf-8', errors='replace') as f:
                        if 'EPRelax' in f.read(4096):
                            title = 'EPRelax session'
                except OSError as fault:
                    print(f"[EPRelax] Could not look inside '{full}': {fault}")
            else:
                image_id = self.file_img_id
            if not is_dir and not fnmatch.fnmatch(title.lower(), pattern):
                hidden += 1
                continue
            if os.path.splitext(full)[0].lower() in loaded:
                image_id = self.loaded_img_id
            idx = self.list.InsertItem(self.list.GetItemCount(), name, image_id)
            self.list.SetItem(idx, 1, title)
            self.items.append((name, full, is_dir))
            self.list.SetItemData(idx, len(self.items) - 1)
            prefix = "0_" if is_dir else "1_"    # folders stay on top whatever the sort
            self.itemDataMap[len(self.items) - 1] = (prefix + low, prefix + (title or name).lower())
        self.list.SetColumnWidth(0, -1)
        if self.GetSortState()[0] >= 0:
            self.SortListItems()
        names = {select_name} if isinstance(select_name, str) else set(select_name or ())
        for ridx in range(self.list.GetItemCount()):
            if self.items[self.list.GetItemData(ridx)][0] in names:
                self.list.Select(ridx)
        if hidden:
            self.parent.statusbar.SetStatusText(
                f"Title '{self.txt_keyword.GetValue().strip()}' hides {hidden} of {len(files)} files", 0)

    # ColumnSorterMixin contract
    def GetListCtrl(self):
        return self.list

    def GetSortImages(self):
        return self.sort_down_img_id, self.sort_up_img_id

    def selected_items(self):
        """[(name, full_path, is_dir), ...] of the selected rows; [] with a warning if none."""
        items, idx = [], self.list.GetFirstSelected()
        while idx != -1:
            items.append(self.items[self.list.GetItemData(idx)])
            idx = self.list.GetNextSelected(idx)
        if not items:
            wx.MessageBox("Please select a file first.", "No file selected", wx.OK | wx.ICON_WARNING)
        return items

    def on_load_clicked(self, event):
        items = self.selected_items()
        if len(items) == 1 and items[0][2]:      # a single folder: open it
            self.set_path(items[0][1])
            return
        for name, full, is_dir in items:
            if not is_dir:
                self.parent.loadData(full)

    def on_load_all(self, event):
        """Load every data file the list shows; folders, sessions and failures are skipped."""
        files = [full for name, full, is_dir in self.items
                 if not is_dir and not name.lower().endswith('.xml')]
        if not files:
            wx.MessageBox(f"This folder holds no '{self.current_filter}' files.", "Load all",
                          wx.OK | wx.ICON_INFORMATION)
            return
        if len(files) > 20 and wx.MessageBox(f"{len(files)} files from '{self.path}' will be loaded.\n\nContinue?",
                                             "Load all", wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
            return
        failed = []
        wait = wx.BusyCursor()
        for full in files:
            fault = self.parent.loadData(full, report=False, refresh=False)
            if fault is not None:
                failed.append(f"{os.path.basename(full)}: {type(fault).__name__}: {fault}")
        del wait
        self.parent.current = len(self.parent.Data) - 1
        self.parent.update_everything()
        self.parent.statusbar.SetStatusText(
            f"Loaded {len(files) - len(failed)} of {len(files)} files" + (f", {len(failed)} failed" if failed else ""), 0)
        if failed:
            wx.MessageBox("These files were skipped:\n\n" + "\n".join(failed[:20])
                          + ("\n..." if len(failed) > 20 else ""), "Load all", wx.OK | wx.ICON_WARNING)

    def on_header_clicked(self, event):
        items = self.selected_items()
        if not items or items[0][2]:
            return
        name, full, _ = items[0]
        try:
            info = Mybr.parse_dsc(os.path.splitext(full)[0] + ('.par' if full.lower().endswith(('.par', '.spc'))
                                                               else '.dsc'))
        except Exception as fault:
            traceback.print_exc()
            wx.MessageBox(f"Could not read the header of '{name}':\n{type(fault).__name__}: {fault}",
                          "Header", wx.OK | wx.ICON_ERROR)
            return
        with DictPopup(self, info, title=name) as dlg:
            dlg.ShowModal()

    def on_filter_changed(self, event):
        """Either half of the filter changed: the extension list, or the title keyword."""
        self.current_filter = self.choice_filter.GetStringSelection()
        self.refresh_file_list()
        event.Skip()                 # the keyword box reaches here on kill-focus, which must carry on

    def on_change_folder(self, event):
        with wx.DirDialog(self, "Select a folder", defaultPath=self.path,
                          style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                self.set_path(dlg.GetPath())

    def on_up_clicked(self, event):
        parent = os.path.dirname(self.path)
        if parent and parent != self.path:
            self.set_path(parent)

    def on_item_activated(self, event):
        name, full, is_dir = self.items[self.list.GetItemData(event.GetIndex())]
        if is_dir:
            self.set_path(full)
        else:
            self.parent.loadData(full)

    def set_path(self, new_path):
        new_path = os.path.abspath(new_path)
        if os.path.isdir(new_path):
            self.path = new_path
            self.lbl_path.SetValue(self.path)
            self.lbl_path.SetInsertionPointEnd()
            self.refresh_file_list()
            self.save_path()

    def on_path_entered(self, event):
        typed = self.lbl_path.GetValue().strip()
        if typed and os.path.isdir(typed):
            if os.path.abspath(typed) != self.path:
                self.set_path(typed)
        else:
            self.lbl_path.SetValue(self.path)
            self.lbl_path.SetInsertionPointEnd()
        event.Skip()


# --------------------------------------------------------------------------- #
# Plot Panel
# --------------------------------------------------------------------------- #
class PlotPanel(wx.Panel):
    """Rows, top to bottom: the trace with its fits, the residual, the distribution of relaxation
    times, and the L-curve. Every shown dataset is overlaid in the first three; the L-curve is a
    diagnostic of one decomposition, so it only ever draws the current dataset."""

    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.parent = mainWindow
        self.show = {'data': True, 'resid': True, 'dist': True, 'lcurve': False}
        self.showGrid = True
        self.norm = 'raw'         # 'raw' | 'area' | 'fraction', the distribution's vertical axis
        self.xmode = 'time'       # 'time' | 'rate', what the distribution is drawn against
        self.sep = 0.0            # gap between distribution traces, in tallest-trace units
        self.stack = False        # the discrete components, stacked, under the trace
        self.stackmode = 'lines'  # 'lines' | 'bands', how those components are drawn
        self.logt = False         # logarithmic time axis on the trace row
        self.snorm = 'raw'        # 'raw' | 'span', the experimental trace's vertical axis
        self.showfit = True       # the regularized fit drawn over the trace
        self.showdisc = True      # the discrete fit drawn over the trace
        self.legend = True        # the dataset names on the trace row
        self.axrows = {}          # axes -> row name of the last draw
        self.picks = {}           # axes -> [(dataset index, x, y)] of the last draw, for clicking
        self.drag = None          # index of the time-window limit being dragged, 0 or 1

        self.figure = Figure(figsize=(2, 2), dpi=100, layout='constrained')
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('figure_leave_event', lambda e: self.parent.statusbar.SetStatusText("", 1))
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()
        self.SetSizer(sizer)

        self.redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self.redraw_debouncer.on_size)

    def on_press(self, event):
        """Left: grab a time-window limit, the dotted lines being draggable, or else make the trace
        under the pointer the current dataset. Right: that trace's colour."""
        mf = self.parent
        if event.inaxes is None or str(self.toolbar.mode):
            return
        if event.button == 1 and mf.current >= 0 and self.axrows.get(event.inaxes) == 'data':
            for j, b in enumerate(mf.Data[mf.current]['win']):
                if abs(event.inaxes.transData.transform((b, 0))[0] - event.x) <= 5:
                    self.drag = j
                    return
        if event.button not in (1, 3):
            return
        # the nearest drawn point within 20 pixels, so clicking an empty part of the axes is not
        # a selection: with traces overlaid, picking the wrong one is worse than picking none
        best, far = None, 20.0
        for ii, x, y in self.picks.get(event.inaxes, ()):
            pts = event.inaxes.transData.transform(np.column_stack((x, y)))
            d = np.hypot(pts[:, 0] - event.x, pts[:, 1] - event.y)
            d = d[np.isfinite(d)]
            if d.size and float(d.min()) < far:
                best, far = ii, float(d.min())
        if best is None:
            return
        if event.button == 1:
            mf.set_current(best)
        else:
            mf.choose_color(best, self)

    def on_release(self, event):
        if self.drag is None:
            return
        self.drag = None
        self.parent.control_panel.update_controls()
        self.update_graph()

    def on_mouse_move(self, event):
        mf = self.parent
        row = self.axrows.get(event.inaxes)
        if self.drag is not None and row == 'data' and event.xdata is not None:
            dd = mf.Data[mf.current]
            dd['win'][self.drag] = float(event.xdata)
            self.update_graph()
            return
        if row is None or event.xdata is None:
            mf.statusbar.SetStatusText("", 1)
            return
        unit = mf.Data[mf.current]['tunit'] if mf.current >= 0 else ''
        if row == 'dist':
            name, u = ("rate", f"1/{unit}") if self.xmode == 'rate' else ("T", unit)
            mf.statusbar.SetStatusText(f"{name} = {event.xdata:.4g} {u}   P = {event.ydata:.4g}", 1)
        elif row == 'lcurve':
            mf.statusbar.SetStatusText(f"residual = {event.xdata:.4g}   ||L a|| = {event.ydata:.4g}", 1)
        else:
            mf.statusbar.SetStatusText(f"t = {event.xdata:.4g} {unit}   y = {event.ydata:.4g}", 1)

    def update_graph(self):
        """Analyse every shown dataset and redraw all visible rows."""
        mf = self.parent
        for ax in self.figure.axes:    # matplotlib warns about the limits when it clears a shared log axis
            ax.set_xscale('linear')
            ax.set_yscale('linear')
        self.figure.clf()
        self.axrows, self.picks = {}, {}
        rows = [k for k in ROWS if self.show[k]]
        if not rows:
            self.canvas.draw_idle()
            return
        gs = self.figure.add_gridspec(len(rows), 1, height_ratios=[ROW_HEIGHT[k] for k in rows])
        axes = {}
        for i, k in enumerate(rows):
            axes[k] = self.figure.add_subplot(gs[i], sharex=axes.get('data') if k == 'resid' else None)
            self.axrows[axes[k]] = k
        # the scales go on before anything is plotted, so autoscaling never has to undo linear margins
        for k in ('data', 'resid') if self.logt else ():
            if k in axes:
                axes[k].set_xscale('log', nonpositive='clip')
        if 'dist' in axes:
            axes['dist'].set_xscale('log')
        if 'lcurve' in axes:
            axes['lcurve'].set_xscale('log')
            axes['lcurve'].set_yscale('log')

        shown = [ii for ii, dd in enumerate(mf.Data) if dd['show']]
        for ii in shown:
            mf.analyse(mf.Data[ii])
        unit = mf.Data[mf.current]['tunit'] if mf.current >= 0 else ''

        # the displayed distribution of every shown dataset, worked out first: the separation is
        # measured in units of the tallest of them, so nothing can be placed until all are known
        amps = {}
        for ii in shown:
            r = mf.Data[ii]['res']
            if r is None:
                continue
            total = max(float(r['a'].sum()), 1e-300)
            if self.norm == 'area':
                amps[ii] = r['a']/(total*float(np.log(r['T'][-1]/r['T'][0]))/max(r['T'].size - 1, 1))
            elif self.norm == 'fraction':
                amps[ii] = r['a']/total
            else:
                amps[ii] = r['a']
        tallest = max((float(a.max()) for a in amps.values() if a.size), default=0.0) or 1.0

        for ii in shown:
            dd = mf.Data[ii]
            r = dd['res']
            if r is None:
                continue
            col = mf.trace_color(dd)
            lw = 1.6 if ii == mf.current else 0.9
            # every trace mapped onto 0..1 over its own window, or left alone; the fits and the
            # stacked components go through the same map, the residual through its scale alone
            lo = float(r['y'].min()) if self.snorm == 'span' else 0.0
            gain = 1.0/max(float(r['y'].max()) - lo, 1e-300) if self.snorm == 'span' else 1.0
            scale = lambda v: (v - lo)*gain
            if 'data' in axes:
                axes['data'].plot(r['t'], scale(r['y']), color=col, lw=lw, ls='none', marker='.',
                                  ms=2.5, label=f"{dd['title']}")
                self.picks.setdefault(axes['data'], []).append((ii, r['t'], scale(r['y'])))
                if self.showfit:
                    axes['data'].plot(r['t'], scale(r['fit']), color=col, lw=lw, alpha=0.9)
                if self.showdisc and dd['disc'] is not None:
                    axes['data'].plot(r['t'], scale(mf.discrete_curve(dd, r['t'])), color=col,
                                      lw=1.0, ls='--', label=f"{dd['title']}, discrete")
                if self.stack and dd['disc'] is not None:
                    # each component laid on the one below it, slowest at the bottom, so every
                    # curve holds its level until its own time comes and then peels away, and the
                    # topmost curve is the whole fit. The horizontal is where that curve starts,
                    # so the gap between two of them is what the component between them is worth
                    d = dd['disc']
                    sign = -1.0 if dd['recovery'] else 1.0
                    level = np.full(r['t'].shape, d['c'])
                    plateau = d['c']
                    for k, i in enumerate(np.argsort(d['T'])[::-1]):
                        top = level + sign*d['a'][i]*np.exp(-r['t']/d['T'][i])
                        plateau += sign*d['a'][i]
                        if self.stackmode == 'bands':
                            axes['data'].fill_between(r['t'], scale(level), scale(top), color=col,
                                                      lw=0, alpha=0.10 + 0.08*k)
                            axes['data'].plot(r['t'], scale(top), color=col, lw=0.7, alpha=0.6)
                        else:
                            axes['data'].plot(r['t'], scale(top), color=col, lw=0.9, alpha=0.85)
                            axes['data'].axhline(scale(plateau), color=col, lw=0.7, alpha=0.5)
                        level = top
            if 'resid' in axes:
                if self.showfit:
                    res = (r['y'] - r['fit'])*gain
                    axes['resid'].plot(r['t'], res, color=col, lw=0.8)
                    self.picks.setdefault(axes['resid'], []).append((ii, r['t'], res))
                if self.showdisc and dd['disc'] is not None:
                    res = (r['y'] - mf.discrete_curve(dd, r['t']))*gain
                    axes['resid'].plot(r['t'], res, color=col, lw=0.7, ls='--')
                    self.picks.setdefault(axes['resid'], []).append((ii, r['t'], res))
            if 'dist' in axes and ii in amps:
                base = self.sep*tallest*shown.index(ii)
                amp = amps[ii] + base
                # a density in ln T is unchanged by T -> 1/T, so the rate view only mirrors it
                x = 1.0/r['T'] if self.xmode == 'rate' else r['T']
                axes['dist'].plot(x, amp, color=col, lw=lw, marker='.', ms=3, label=f"{dd['title']}")
                self.picks.setdefault(axes['dist'], []).append((ii, x, amp))
                axes['dist'].fill_between(x, base, amp, color=col, alpha=0.12)
                if dd['disc'] is not None:
                    for Td in dd['disc']['T']:
                        axes['dist'].axvline(1.0/Td if self.xmode == 'rate' else Td,
                                             color=col, lw=0.9, ls='--', alpha=0.8)

        if 'data' in axes:
            ax = axes['data']
            ax.set_ylabel("Signal, 0 to 1" if self.snorm == 'span' else "Signal")
            ax.set_xlabel(f"t, {unit}" if 'resid' not in axes else "")
            if mf.current >= 0:
                for b in mf.Data[mf.current]['win']:
                    ax.axvline(b, color='0.45', lw=0.9, ls=':')
            if shown and self.legend:
                ax.legend(fontsize='x-small', frameon=False)
        if 'resid' in axes:
            axes['resid'].axhline(0, color='0.5', lw=0.8)
            axes['resid'].set_ylabel("Residual")
            axes['resid'].set_xlabel(f"t, {unit}")
        if 'dist' in axes:
            ax = axes['dist']
            ax.set_xlabel(f"Rate, 1/{unit}" if self.xmode == 'rate' else f"Time constant, {unit}")
            ax.set_ylabel({'area': "Density per unit ln T", 'fraction': "Fraction of total"}
                          .get(self.norm, "Amplitude"))
            ax.set_ylim(bottom=0)
        if 'lcurve' in axes:
            ax = axes['lcurve']
            r = mf.Data[mf.current]['res'] if mf.current >= 0 else None
            if r is not None and r['rho'] is not None:
                ax.plot(r['rho'], r['eta'], color=theme.theme_colors()["accent"], lw=1.2, marker='.', ms=3)
                # the solve that is actually drawn, not the nearest swept alpha: the knob moves the
                # marker along the curve between the points, which is the whole use of watching it
                ax.plot(*r['point'], 'o', ms=9, mfc='none', mew=1.8,
                        color=theme.theme_colors()["ax_fg"], label=f"alpha = {r['alpha']:.3g}")
                ax.legend(fontsize='x-small', frameon=False)
            else:
                ax.text(0.5, 0.5, "no current dataset, or its trace could not be decomposed",
                        ha='center', va='center',
                        transform=ax.transAxes, fontsize='small', color='0.5')
            ax.set_xlabel("Residual norm")
            ax.set_ylabel("||L a||")

        for ax in axes.values():
            ax.grid(True, alpha=0.3) if self.showGrid else ax.grid(False)
        theme.style_figure(self.figure, self.figure.axes)
        self.canvas.draw_idle()
        mf.control_panel.update_stats()
        if mf.rate_window is not None:     # its points are read from the very fits drawn here
            mf.rate_window.update_graph()


# --------------------------------------------------------------------------- #
# Control Panel
# --------------------------------------------------------------------------- #
class ControlPanel(wx.Panel):
    """The dataset list, then four tabbed pages of cards. Most settings belong to the current
    dataset; the grid density, the penalty and where alpha comes from are shared by every dataset,
    because comparing one distribution with another only means something when both were made the
    same way. The manual alpha stays individual, since that is the knob you turn per trace."""

    TABS = ("Signal", "Regularization", "Discrete fit", "Plot")
    DISC_ROWS = ("T", "± 95%", "a", "± 95%", "iLT low", "iLT high", "note")

    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.mf = mainWindow
        # while the list is being filled it ticks and selects rows itself, and every one of those
        # would otherwise read as a click: a redraw per row, and the whole analysis behind it
        self.filling = False
        outer = wx.BoxSizer(wx.VERTICAL)

        # Datasets card, above the tabs: tick to show, click to make current
        card = theme.Card(self, "Datasets")
        self.filelist = wx.ListCtrl(card, size=(-1, 128),
                                    style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.BORDER_NONE)
        self.filelist.InsertColumn(0, "Dataset", width=168)
        self.filelist.InsertColumn(1, self.mf.rates['name'], width=76)
        # a column wide enough to read titles in must not become a width the card may never go
        # below: the list scrolls, where the buttons under it cannot
        self.filelist.SetMinSize(wx.Size(60, 128))
        self.filelist.EnableCheckBoxes(True)
        self.filelist.SetToolTip("Tick a dataset to draw it, click it to make it the one the tabs below edit.\n"
                                 "Each dataset keeps its own channel, window, time limits and manual alpha;\n"
                                 "the grid density, the penalty and where alpha comes from are shared by all.\n"
                                 "Double-click a row, or press Enter on it, to type the parameter it was\n"
                                 "measured at -- a temperature, a concentration -- which the rates window\n"
                                 "plots the fitted relaxation rates against. The column is named there.")
        self.filelist.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda e: self.edit_param(e.GetIndex()))
        self.filelist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_select)
        self.filelist.Bind(wx.EVT_LIST_ITEM_CHECKED, self.on_show)
        self.filelist.Bind(wx.EVT_LIST_ITEM_UNCHECKED, self.on_show)
        self.filelist.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK,
                           lambda e: self.mf.choose_color(e.GetIndex(), self.filelist))
        self.btn_showall = PillButton(card, "Show all", color_key="pill_calm",
                                      draw_badge=theme.draw_plus_badge)
        self.btn_showall.SetToolTip("Tick every dataset, so they are all drawn.")
        self.btn_showall.Bind(wx.EVT_BUTTON, lambda e: self.on_showing(True))
        self.btn_showonly = PillButton(card, "Show only", color_key="pill_calm",
                                       draw_badge=theme.draw_minus_badge)
        self.btn_showonly.SetToolTip("Draw the selected dataset alone: untick all the others.")
        self.btn_showonly.Bind(wx.EVT_BUTTON, lambda e: self.on_showing(False))
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.btn_showall, 0, wx.RIGHT, 4)
        row.Add(self.btn_showonly, 0)
        card.body.Add(row, 0, wx.ALL, 3)
        card.body.Add(self.filelist, 0, wx.EXPAND | wx.ALL, 3)
        self.btn_remove = PillButton(card, "Remove", color_key="pill_warn", draw_badge=theme.draw_cross_badge)
        self.btn_remove.SetToolTip("Drop the selected dataset.")
        self.btn_remove.Bind(wx.EVT_BUTTON, self.on_remove)
        self.btn_removeall = PillButton(card, "Remove all", color_key="pill_warn",
                                        draw_badge=theme.draw_cross_badge)
        self.btn_removeall.SetToolTip("Drop every loaded dataset and start over. Asks first.")
        self.btn_removeall.Bind(wx.EVT_BUTTON, self.on_remove_all)
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.btn_remove, 0, wx.RIGHT, 4)
        row.Add(self.btn_removeall, 0)
        card.body.Add(row, 0, wx.ALL, 3)
        self.fix_buttons((self.btn_showall, self.btn_showonly, self.btn_remove, self.btn_removeall))
        outer.Add(card, 0, wx.EXPAND | wx.ALL, 5)
        # this card sits above the tabs and cannot scroll sideways as they do, so it is what keeps
        # the sash from squashing two buttons into the width of one
        self.floor = card.GetEffectiveMinSize().width + 16

        self.book = wx.Simplebook(self)
        pages = {name: scrolled.ScrolledPanel(self.book, style=wx.BORDER_NONE) for name in self.TABS}
        stacks = {name: wx.BoxSizer(wx.VERTICAL) for name in self.TABS}

        # ---- Signal tab -------------------------------------------------------------------
        page = pages["Signal"]
        card = theme.Card(page, "Signal")
        self.seg_comp = SegmentedPill(card, ["Real", "Imag.", "Magn."])
        self.seg_comp.SetToolTip("Which channel of a complex dataset is decomposed. Real is the file's own real\n"
                                 "channel, untouched -- no phase is applied anywhere, so a mis-set spectrometer\n"
                                 "phase stays visible instead of being corrected behind your back.\n"
                                 "Magnitude is positive-biased in the noise, which shows up as spurious long-T\n"
                                 "amplitude in the tail.")
        card.body.Add(self.row(card, "Channel:", self.seg_comp), 0, wx.EXPAND | wx.ALL, 3)
        self.spin_slice = wx.SpinCtrl(card, min=0, max=0, initial=0, size=(70, -1))
        self.spin_slice.SetToolTip("Which slice of a 2D dataset is decomposed, counted from 0.")
        self.chk_sum = PillCheckBox(card, "Sum")
        self.chk_sum.SetToolTip("Add all slices together before decomposing, for signal-to-noise.\n"
                                "Only meaningful when every slice decays the same way.")
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.label(card, "Slice:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        row.Add(self.spin_slice, 0, wx.RIGHT, 6)
        row.Add(self.chk_sum, 0, wx.ALIGN_CENTER_VERTICAL)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        self.spin_from = wx.SpinCtrlDouble(card, min=-1e12, max=1e12, inc=1, initial=0, size=(104, -1))
        self.spin_to = wx.SpinCtrlDouble(card, min=-1e12, max=1e12, inc=1, initial=0, size=(104, -1))
        for spin, tip in ((self.spin_from, "First time point the decomposition sees. Raise it to cut the dead time,\n"
                                           "the ringdown, or the rising edge of an echo -- and to clear a spike\n"
                                           "at the short end of the distribution, which is where anything faster\n"
                                           "than the grid ends up."),
                          (self.spin_to, "Last time point the decomposition sees. Lower it to cut a noisy tail,\n"
                                         "which otherwise feeds the longest-T end of the distribution.")):
            spin.SetDigits(3)
            spin.SetToolTip(tip + "\nBoth limits are the dotted lines on the trace, and can be dragged there.\n"
                                  "Right-click to set the arrow step.")
        self.btn_full = PillButton(card, "Full", color_key="pill_calm")
        self.btn_full.SetToolTip("Put both limits back to the ends of the recorded trace.")
        self.btn_full.Bind(wx.EVT_BUTTON, self.on_full)
        self.fix_buttons((self.btn_full,))
        self.lbl_window = wx.StaticText(card, label="")
        grid = wx.FlexGridSizer(cols=2, vgap=2, hgap=6)
        grid.AddGrowableCol(1)
        for label, spin in (("From:", self.spin_from), ("To:", self.spin_to)):
            grid.Add(self.label(card, label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(spin, 0, wx.EXPAND)
        card.body.Add(grid, 0, wx.EXPAND | wx.ALL, 3)
        card.body.Add(self.lbl_window, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 3)
        card.body.Add(self.btn_full, 0, wx.ALL, 3)
        self.chk_rec = PillCheckBox(card, "Recovery")
        self.chk_rec.SetToolTip("The exponentials come up to the offset from below: inversion or saturation\n"
                                "recovery rather than a decay. The amplitudes stay non-negative, only the\n"
                                "sign in front of the sum turns over. Set from the trace when the file loads.")
        self.chk_off = PillCheckBox(card, "Fit offset")
        self.chk_off.SetValue(True)
        self.chk_off.SetToolTip("Fit a free constant along with the amplitudes. Leave it on unless the\n"
                                "baseline is already exactly zero: a residual constant otherwise has nowhere\n"
                                "to go but the longest-T end of the distribution.")
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.chk_rec, 0, wx.RIGHT, 4)
        row.Add(self.chk_off, 0)
        card.body.Add(row, 0, wx.ALL, 3)
        stacks["Signal"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        card = theme.Card(page, "ESEEM suppression")
        self.chk_smooth = PillCheckBox(card, "Suppress")
        self.chk_smooth.SetToolTip("Run a boxcar average over the trace before decomposing, to take nuclear\n"
                                   "modulation out of a 2-pulse decay. A boxcar only rescales an exponential\n"
                                   "(by sinh(w/2T)/(w/2T)), so the relaxation times come out unchanged while\n"
                                   "the amplitudes are slightly inflated. Off by default: the envelope is a\n"
                                   "model choice sitting between the file and the result.")
        card.body.Add(self.chk_smooth, 0, wx.ALL, 3)
        self.spin_swin = wx.SpinCtrlDouble(card, min=0, max=1e12, inc=1, initial=0, size=(104, -1))
        self.spin_swin.SetDigits(3)
        self.spin_swin.SetToolTip("Width of the boxcar, in the dataset's own time unit. Set it to one\n"
                                  "modulation period -- 1/nu_I for the nucleus you see -- and the modulation\n"
                                  "cancels. Right-click to set the arrow step.")
        card.body.Add(self.row(card, "Width:", self.spin_swin), 0, wx.EXPAND | wx.ALL, 3)
        stacks["Signal"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        card = theme.Card(page, "Resampling")
        self.chk_resample = PillCheckBox(card, "Log t")
        self.chk_resample.SetToolTip("Average the trace into bins evenly spaced in log t, after the boxcar and\n"
                                     "inside the time window. A linearly sampled trace spends most of its\n"
                                     "points on the tail, where little is happening; binning them hands that\n"
                                     "redundancy back as signal-to-noise and leaves the early decay alone,\n"
                                     "since a bin holding no points is dropped rather than interpolated.\n"
                                     "It also shrinks the kernel, so every solve is quicker.\n"
                                     "The trace row then draws the binned points, which are what is fitted.")
        card.body.Add(self.chk_resample, 0, wx.ALL, 3)
        self.spin_perdec = wx.SpinCtrlDouble(card, min=2, max=500, inc=1, initial=16, size=(104, -1))
        self.spin_perdec.SetDigits(1)
        self.spin_perdec.SetToolTip("Bins per decade of time, so the density is the same for every dataset\n"
                                    "however long it ran and two traces stay comparable. A trace spanning two\n"
                                    "decades at 16 per decade comes out as about 32 points.\n"
                                    "Right-click to set the arrow step.")
        card.body.Add(self.row(card, "Per decade:", self.spin_perdec), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_rweight = SegmentedPill(card, ["By count", "Equal"])
        self.seg_rweight.SetToolTip("How much each binned point counts for in the fit.\n"
                                    "By count: a bin of m points carries weight sqrt(m), which is how much\n"
                                    "more precise its average is. The fit is then statistically the same as\n"
                                    "fitting the raw trace, so resampling changes the picture and the speed\n"
                                    "but not the answer.\n"
                                    "Equal: every bin counts once, so each decade of time has a say in\n"
                                    "proportion to its bins rather than to how densely it was sampled. That\n"
                                    "matches the logarithmic relaxation-time grid and stops a long flat tail\n"
                                    "from dominating -- but it is a different question, and the answer moves.")
        card.body.Add(self.row(card, "Weights:", self.seg_rweight), 0, wx.EXPAND | wx.ALL, 3)
        self.btn_resall = PillButton(card, "Apply to All", color_key="pill_calm",
                                     draw_badge=theme.draw_plus_badge)
        self.btn_resall.SetToolTip("Give every loaded dataset the resampling set here: the switch, the\n"
                                   "density and the weighting. A density per decade means the same thing\n"
                                   "on any time axis, so nothing is rescaled on the way and traces of very\n"
                                   "different length still come out comparable.\n"
                                   "Switching it off here switches it off everywhere, the same way.")
        self.btn_resall.Bind(wx.EVT_BUTTON, self.on_resample_all)
        self.fix_buttons((self.btn_resall,))
        card.body.Add(self.btn_resall, 0, wx.ALL, 3)
        stacks["Signal"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # ---- Regularization tab -----------------------------------------------------------
        page = pages["Regularization"]
        card = theme.Card(page, "Relaxation-time grid")
        self.chk_auto = PillCheckBox(card, "Auto")
        self.chk_auto.SetValue(True)
        self.chk_auto.SetToolTip("Span the grid from a hundredth of the time window to three times it, which\n"
                                 "is roughly what a trace of that length can resolve. Untick to type your own\n"
                                 "range: a grid much wider than the data supports invents amplitude at both\n"
                                 "ends -- a relaxation time far longer than the measurement is the same curve\n"
                                 "as the fitted offset, and the two then simply trade amplitude. The Auto\n"
                                 "limits are left in the boxes to start from.")
        card.body.Add(self.chk_auto, 0, wx.ALL, 3)
        self.spin_tmin = wx.SpinCtrlDouble(card, min=1e-9, max=1e12, inc=1, initial=1, size=(104, -1))
        self.spin_tmax = wx.SpinCtrlDouble(card, min=1e-9, max=1e12, inc=1, initial=1000, size=(104, -1))
        for spin in (self.spin_tmin, self.spin_tmax):
            spin.SetDigits(3)
            spin.SetToolTip("End of the logarithmically spaced grid of relaxation times, in the dataset's\n"
                            "own time unit -- always a time, even when the plot is showing rates.\n"
                            "Right-click to set the arrow step.")
        self.spin_npts = wx.SpinCtrl(card, min=10, max=400, initial=80, size=(104, -1))
        self.spin_npts.SetToolTip("Number of grid points, shared by every dataset. The regularizer, not the\n"
                                  "grid, sets how much structure the answer can have, so more points cost\n"
                                  "time and change little above about 60.")
        grid = wx.FlexGridSizer(cols=2, vgap=2, hgap=6)
        grid.AddGrowableCol(1)
        for label, spin in (("T min:", self.spin_tmin), ("T max:", self.spin_tmax), ("Points:", self.spin_npts)):
            grid.Add(self.label(card, label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(spin, 0, wx.EXPAND)
        card.body.Add(grid, 0, wx.EXPAND | wx.ALL, 3)
        stacks["Regularization"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        card = theme.Card(page, "Regularization")
        self.seg_pen = SegmentedPill(card, ["Curvature", "Amplitude"])
        self.seg_pen.SetToolTip("Shared by every dataset.\n"
                                "Curvature: the second difference along the log-T grid, as in CONTIN. It asks\n"
                                "for the smoothest distribution consistent with the data and leaves broad\n"
                                "features at full height. The distribution is taken as zero outside the grid,\n"
                                "so the two end points are penalised like the middle; left open they would be\n"
                                "2.4 times cheaper and collect a spike.\n"
                                "Amplitude: plain Tikhonov on the amplitudes themselves. Very stable, but it\n"
                                "shrinks a broad distribution instead of smoothing it.")
        card.body.Add(self.row(card, "Penalty:", self.seg_pen), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_amode = SegmentedPill(card, ["Manual", "L-curve", "GCV", "F-test"])
        self.seg_amode.SetToolTip("Shared by every dataset, so all of them are regularized on the same rule.\n"
                                  "Manual: alpha is the knob below, and that one is per dataset.\n"
                                  "L-curve: Hansen's corner, the maximum curvature of the log-log curve of\n"
                                  "||L a|| against the residual norm.\n"
                                  "GCV: the minimum of the generalized cross-validation score.\n"
                                  "F-test: CONTIN's own rule -- the largest alpha whose residual sum of squares\n"
                                  "is still statistically indistinguishable from the unregularized one.")
        card.body.Add(self.row(card, "alpha from:", self.seg_amode), 0, wx.EXPAND | wx.ALL, 3)
        self.chk_autoalpha = PillCheckBox(card, "Auto range")
        self.chk_autoalpha.SetValue(True)
        self.chk_autoalpha.SetToolTip("Find the swept range from each trace instead of fixing it: climb until the\n"
                                      "penalty has flattened the distribution to a single effective component,\n"
                                      "and start a decade below where the solve stops changing. The corner moves\n"
                                      "by decades between traces -- with the noise, the point count and how many\n"
                                      "components there are -- so a fixed ladder can leave the whole informative\n"
                                      "part of the curve outside the window, which shows up as an alpha chosen\n"
                                      "at an end of the L-curve with the residual still at its floor.\n"
                                      "Untick to sweep the fixed range below, the same for every dataset, which\n"
                                      "is what makes two L-curves directly comparable.")
        card.body.Add(self.chk_autoalpha, 0, wx.ALL, 3)
        self.spin_alo = wx.SpinCtrlDouble(card, min=-14, max=14, inc=1, initial=-5, size=(104, -1))
        self.spin_ahi = wx.SpinCtrlDouble(card, min=-14, max=14, inc=1, initial=5, size=(104, -1))
        for spin, tip in ((self.spin_alo, "log10 of the smallest alpha swept."),
                          (self.spin_ahi, "log10 of the largest alpha swept.")):
            spin.SetDigits(1)
            spin.SetToolTip(tip + f"\nThe range is shared by every dataset and carries {SWEEP} points.\n"
                                  "Right-click to set the arrow step.")
        card.body.Add(self.row(card, "Sweep from:", self.spin_alo), 0, wx.EXPAND | wx.ALL, 3)
        card.body.Add(self.row(card, "Sweep to:", self.spin_ahi), 0, wx.EXPAND | wx.ALL, 3)
        self.spin_alpha = wx.SpinCtrlDouble(card, min=-14, max=14, inc=0.1, initial=-2, size=(104, -1))
        self.spin_alpha.SetDigits(2)
        self.spin_alpha.SetToolTip("log10 of alpha, because the useful range covers decades. This one belongs\n"
                                   "to the current dataset. The data is scaled to unit peak before the solve,\n"
                                   "so alpha means the same thing whatever the signal units. Too small and the\n"
                                   "distribution breaks into spikes; too large and everything melts into one\n"
                                   "broad hump -- watch the residual row while you turn it. The automatic\n"
                                   "modes write their choice back into this box.\n"
                                   "Right-click to set the arrow step.")
        card.body.Add(self.row(card, "log10 alpha:", self.spin_alpha), 0, wx.EXPAND | wx.ALL, 3)
        self.spin_prob = wx.SpinCtrlDouble(card, min=0.01, max=0.99, inc=0.05, initial=0.5, size=(104, -1))
        self.spin_prob.SetDigits(2)
        self.spin_prob.SetToolTip("The F-test's probability to reject, shared by every dataset since it is a\n"
                                  "knob of the criterion rather than of a trace. CONTIN uses 0.5: the chosen\n"
                                  "solution is the smoothest one that a coin toss would not call a worse fit\n"
                                  "than the unregularized answer. Raise it to allow more smoothing.")
        card.body.Add(self.row(card, "Probability:", self.spin_prob), 0, wx.EXPAND | wx.ALL, 3)
        self.lbl_stat = wx.StaticText(card, label="")
        card.body.Add(self.lbl_stat, 0, wx.EXPAND | wx.ALL, 3)
        stacks["Regularization"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # ---- Discrete fit tab -------------------------------------------------------------
        page = pages["Discrete fit"]
        card = theme.Card(page, "Discrete fit")
        self.spin_ndisc = wx.SpinCtrl(card, min=1, max=6, initial=2, size=(104, -1))
        self.spin_ndisc.SetToolTip("How many exponentials to fit. Their starting relaxation times are the\n"
                                   "largest maxima of the regularized distribution, in order, so the two\n"
                                   "analyses are asked the same question; a non-linear least-squares fit then\n"
                                   "refines them. Ask for more components than the distribution has peaks and\n"
                                   "the extra ones start from the open grid, where they rarely stay put.")
        # this one keeps its width: a single digit next to a label, with nothing to gain from
        # stretching across the tab or from giving width back when the pane narrows
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(self.label(card, "Components:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        row.Add(self.spin_ndisc, 0)
        card.body.Add(row, 0, wx.ALL, 3)
        self.chk_fixt = PillCheckBox(card, "Hold T at the peaks")
        self.chk_fixt.SetToolTip("Do not fit the time constants: hold each one at the average of a\n"
                                 "distribution peak and fit only the amplitudes and the offset.\n"
                                 "The peaks are taken tallest first, so the tallest becomes component 1,\n"
                                 "and each average is the amplitude-weighted mean of log10 T across that\n"
                                 "peak -- the geometric mean, which is what a logarithmic grid makes of an\n"
                                 "average and the centre the 95 % peak width is measured about.\n"
                                 "Asking for more components than the distribution has peaks fits as many\n"
                                 "as there are peaks: there is nothing to hold the extras at.\n"
                                 "A held time constant is not a fitted parameter, so the interval beside\n"
                                 "it is marked ~ : it is the one T would have had if it had been free at\n"
                                 "that value, which says how tightly the trace pins T down there rather\n"
                                 "than whether that value is the right one.")
        card.body.Add(self.chk_fixt, 0, wx.ALL, 3)
        self.btn_fit = PillButton(card, "Fit", color_key="pill_accent_2", draw_badge=theme.draw_play_badge)
        self.btn_fit.SetToolTip("Run the fit now, from the distribution as it currently stands.")
        self.btn_fit.Bind(wx.EVT_BUTTON, self.on_fit)
        self.btn_clear = PillButton(card, "Clear", color_key="pill_calm", draw_badge=theme.draw_minus_badge)
        self.btn_clear.SetToolTip("Drop the discrete fit of the current dataset.")
        self.btn_clear.Bind(wx.EVT_BUTTON, self.on_clear)
        self.btn_fitall = PillButton(card, "Fit all", color_key="pill_accent_1",
                                     draw_badge=theme.draw_play_badge)
        self.btn_fitall.SetToolTip("Run the same fit on every loaded dataset: the component count and the\n"
                                   "hold-T switch above are handed to all of them first, so the whole\n"
                                   "session is fitted the same way and the table compares like with like.\n"
                                   "Each dataset still starts from its own distribution. Datasets that do\n"
                                   "not converge are listed afterwards rather than stopping the run.")
        self.btn_fitall.Bind(wx.EVT_BUTTON, self.on_fit_all)
        self.btn_savedisc = PillButton(card, "Save", color_key="pill_load",
                                       draw_badge=theme.draw_textlines_badge)
        self.btn_savedisc.SetToolTip("Write every component of every fitted dataset to one CSV table,\n"
                                     "relaxation times in seconds throughout so the column sorts and\n"
                                     "plots across samples whatever unit each trace reads in.")
        self.btn_savedisc.Bind(wx.EVT_BUTTON, lambda e: self.mf.save_discrete())
        # two to a row rather than four, which would set the whole tab's smallest width on its own
        self.fix_buttons((self.btn_fit, self.btn_fitall, self.btn_clear, self.btn_savedisc))
        for first, second in ((self.btn_fit, self.btn_fitall), (self.btn_clear, self.btn_savedisc)):
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(first, 0, wx.RIGHT, 4)
            row.Add(second, 0)
            card.body.Add(row, 0, wx.ALL, 3)
        self.disclist = wx.ListCtrl(card, size=(-1, 176), style=wx.LC_REPORT | wx.BORDER_NONE)
        self.disclist.InsertColumn(0, "Parameter", width=84)
        for row, name in enumerate(self.DISC_ROWS):
            self.disclist.InsertItem(row, name)
        # the component columns arrive with the fit; only the names column is here from the start,
        # so the table never sets the tab's smallest width and the buttons keep theirs
        self.disclist.SetMinSize(wx.Size(60, 176))
        self.disclist.SetToolTip("One column per component.\n"
                                 "T and a with their 95 % intervals come from the non-linear fit's own\n"
                                 "trust region: that is how tightly the trace pins the component down, and it\n"
                                 "goes wide when two components are not really separable. An interval marked\n"
                                 "~ belongs to a held T and is an estimate: the curvature of the residual\n"
                                 "surface along T at the held value, with the amplitudes free to move with\n"
                                 "it, worked out as though that value were the best fit.\n"
                                 "'iLT low' and 'iLT high' are something else entirely: the ends of the\n"
                                 "distribution peak the component sits in, from a Gaussian fitted to that\n"
                                 "peak in log10 T and taken out to 95 %. They are the width of the\n"
                                 "distribution, not an error on T, they are not symmetric about T because a\n"
                                 "width in log T is a factor rather than a difference, and a perfectly\n"
                                 "monoexponential sample still has them because the regularizer smooths.\n"
                                 "They are left empty when the fitted T sits on no peak at all -- on the\n"
                                 "flank of one, or in the valley between two -- since there is then no peak\n"
                                 "there to have a width. A held T is measured on the peak it was held at.\n"
                                 "The note row says when the Gaussian could not be trusted.")
        card.body.Add(self.disclist, 1, wx.EXPAND | wx.ALL, 3)
        self.lbl_disc = wx.StaticText(card, label="")
        card.body.Add(self.lbl_disc, 0, wx.EXPAND | wx.ALL, 3)
        stacks["Discrete fit"].Add(card, 1, wx.EXPAND | wx.ALL, 5)

        # ---- Plot tab ---------------------------------------------------------------------
        page = pages["Plot"]
        card = theme.Card(page, "Plot")
        self.seg_xmode = SegmentedPill(card, ["Time constant", "Rate"])
        self.seg_xmode.SetToolTip("What the distribution is drawn against: the relaxation time T, or the rate\n"
                                  "1/T in the reciprocal of the same unit. Both axes are logarithmic and the\n"
                                  "amplitude is a density in ln T, which is unchanged by T -> 1/T, so the\n"
                                  "curve is simply mirrored and no number has to be rescaled. The grid boxes\n"
                                  "stay in time either way.")
        card.body.Add(self.row(card, "Axis:", self.seg_xmode), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_norm = SegmentedPill(card, ["Raw", "Area", "Fraction"])
        self.seg_norm.SetToolTip("Raw: the fitted amplitudes in the signal's own units, so the height is\n"
                                 "physical and the area under the curve is the total decaying signal.\n"
                                 "Area: a density in ln T -- the integral over ln T is set to 1, which is\n"
                                 "what compares shapes between samples of different concentration. The curve\n"
                                 "itself runs well above 1, since it is a density and not a share.\n"
                                 "Fraction: every point divided by the sum of all of them, so the points add\n"
                                 "up to 1 and each one reads directly as its share of the decaying signal.\n"
                                 "Unlike Area this is not a density, so its height depends on how many grid\n"
                                 "points the session is using.")
        card.body.Add(self.row(card, "Distribution:", self.seg_norm), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_snorm = SegmentedPill(card, ["Raw", "0 to 1"])
        self.seg_snorm.SetToolTip("The experimental trace's vertical axis, the same for every dataset.\n"
                                  "Raw: the signal as recorded, in the file's own units.\n"
                                  "0 to 1: each trace mapped onto 0..1 over its own time window, lowest\n"
                                  "point to 0 and highest to 1, so traces of different intensity lie on top\n"
                                  "of each other and a decay and a recovery share one scale. The fits and\n"
                                  "the stacked components follow the same map and the residual its scale,\n"
                                  "so the rows stay consistent; a partial inversion is no longer visible as\n"
                                  "such, since every trace is stretched to the full height.")
        card.body.Add(self.row(card, "Signal:", self.seg_snorm), 0, wx.EXPAND | wx.ALL, 3)
        self.chk_showfit = PillCheckBox(card, "Regularized fit")
        self.chk_showfit.SetValue(True)
        self.chk_showfit.SetToolTip("Draw the regularized fit over the experimental trace, for every\n"
                                    "dataset, and its residual in the residual row. Untick both fits\n"
                                    "and the residual row is left empty.")
        card.body.Add(self.chk_showfit, 0, wx.ALL, 3)
        self.chk_showdisc = PillCheckBox(card, "Discrete fit")
        self.chk_showdisc.SetValue(True)
        self.chk_showdisc.SetToolTip("Draw the discrete fit, dashed, over the experimental trace, for every\n"
                                     "dataset that has one, and its residual in the residual row.")
        card.body.Add(self.chk_showdisc, 0, wx.ALL, 3)
        self.chk_legend = PillCheckBox(card, "Legend")
        self.chk_legend.SetValue(True)
        self.chk_legend.SetToolTip("The box of dataset names on the trace row. Untick it for a crowded\n"
                                   "overlay, or for a figure whose caption names the traces anyway; the\n"
                                   "colours are still the ones the Datasets list shows.\n"
                                   "The L-curve keeps its own alpha box, which is a readout rather than\n"
                                   "a key to the colours.")
        card.body.Add(self.chk_legend, 0, wx.ALL, 3)
        self.spin_sep = wx.SpinCtrlDouble(card, min=0, max=100, inc=0.1, initial=0, size=(104, -1))
        self.spin_sep.SetDigits(3)
        self.spin_sep.SetToolTip("Vertical gap between consecutive distribution traces, in units of the\n"
                                 "tallest one, in the order the Datasets list has them. Fan them out to see\n"
                                 "each one on its own, then set it back to 0 to compare them against a\n"
                                 "common baseline. Only the distribution row moves; with a gap set, its\n"
                                 "vertical axis no longer reads as amplitude for anything but the first\n"
                                 "trace. Right-click to set the arrow step.")
        card.body.Add(self.row(card, "Separation:", self.spin_sep), 0, wx.EXPAND | wx.ALL, 3)
        self.chk_stack = PillCheckBox(card, "Components")
        self.chk_stack.SetToolTip("Under each fitted trace, the discrete components stacked on one another\n"
                                  "starting at the fitted offset, slowest at the bottom. Each curve holds\n"
                                  "its level until its own time constant comes and then peels away, so the\n"
                                  "gap between two curves is what the component between them is worth and\n"
                                  "the topmost curve is the whole fit.\n"
                                  "Needs a discrete fit; datasets without one are drawn as before.")
        card.body.Add(self.chk_stack, 0, wx.ALL, 3)
        self.seg_stackmode = SegmentedPill(card, ["Lines", "Bands"])
        self.seg_stackmode.SetToolTip("Lines: each cumulative curve as a thin line, with a horizontal where it\n"
                                      "starts. Stays readable with several datasets overlaid and prints well.\n"
                                      "Bands: the gaps filled in, which reads more easily when one dataset is\n"
                                      "shown on its own.")
        card.body.Add(self.row(card, "Drawn as:", self.seg_stackmode), 0, wx.EXPAND | wx.ALL, 3)
        self.chk_logt = PillCheckBox(card, "Log t")
        self.chk_logt.SetToolTip("Logarithmic time axis on the trace, which spreads a recovery curve measured\n"
                                 "on a logarithmic time list out evenly.")
        card.body.Add(self.chk_logt, 0, wx.ALL, 3)
        for w in (self.seg_xmode, self.seg_norm, self.seg_snorm, self.chk_showfit, self.chk_showdisc,
                  self.chk_legend, self.chk_stack, self.seg_stackmode, self.chk_logt):
            w.Bind(wx.EVT_RADIOBOX if isinstance(w, SegmentedPill) else wx.EVT_CHECKBOX, self.on_display)
        self.spin_sep.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_display)
        theme.bind_spin_enter(self.spin_sep, self.on_display)
        for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
            self.spin_sep.Bind(evt, lambda e: self.on_step(self.spin_sep))
        stacks["Plot"].Add(card, 0, wx.EXPAND | wx.ALL, 5)

        for name in self.TABS:
            pages[name].SetSizer(stacks[name])
            pages[name].SetupScrolling(scroll_x=True)
            pages[name].Bind(wx.EVT_MOUSEWHEEL, lambda e, p=pages[name]: self.on_wheel(e, p))
            self.book.AddPage(pages[name], name)
        self.tabs = PillTabBar(self, self.book, list(self.TABS),
                               ["The trace the decomposition is given, and the ESEEM boxcar",
                                "The relaxation-time grid and how hard the answer is smoothed",
                                "The conventional fit, started from the distribution, with its error margins",
                                "What the plot draws and how it is scaled"])
        outer.Add(self.tabs, 0, wx.EXPAND)
        outer.Add(self.book, 1, wx.EXPAND)

        # every control that belongs to a dataset or to the whole session goes through one handler
        for w in (self.seg_comp, self.seg_pen, self.seg_amode, self.seg_rweight):
            w.Bind(wx.EVT_RADIOBOX, self.on_change)
        for w in (self.chk_sum, self.chk_rec, self.chk_off, self.chk_smooth, self.chk_auto,
                  self.chk_resample, self.chk_autoalpha, self.chk_fixt):
            w.Bind(wx.EVT_CHECKBOX, self.on_change)
        for w in (self.spin_from, self.spin_to, self.spin_swin, self.spin_perdec, self.spin_tmin,
                  self.spin_tmax, self.spin_alpha, self.spin_alo, self.spin_ahi, self.spin_prob):
            w.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_change)
            theme.bind_spin_enter(w, self.on_change)
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                w.Bind(evt, lambda e, s=w: self.on_step(s))
        for w in (self.spin_slice, self.spin_npts, self.spin_ndisc):
            w.Bind(wx.EVT_SPINCTRL, self.on_change)
            theme.bind_spin_enter(w, self.on_change)

        self.SetSizer(outer)
        self.allow_shrinking(self)
        self.enable()

    def fix_buttons(self, buttons):
        """One fixed width for a card's buttons: the widest label of the group, so no text is cut,
        none is wider than it has to be, and none stretches when the pane does."""
        widest = max(b.GetBestSize().width for b in buttons)
        for b in buttons:
            b.SetMinSize(wx.Size(widest, b.GetBestSize().height))
            b.SetMaxSize(wx.Size(widest, -1))

    def label(self, card, text):
        """A row label. 76 px lines the controls up down a card, but the label asks for more when
        its own text needs it: a fixed width is also a maximum, and a clipped label is just wrong."""
        lbl = wx.StaticText(card, label=text)
        lbl.SetMinSize(wx.Size(max(76, lbl.GetBestSize().width), -1))
        return lbl

    def row(self, card, label, widget):
        """A label of the fixed card width with its control beside it."""
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(self.label(card, label), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        box.Add(widget, 1, wx.EXPAND)
        return box

    def allow_shrinking(self, window):
        """Let a narrow pane take its width out of the number boxes and nothing else. Every control
        is laid out at a comfortable width, which a sizer then reads as a width it may never go
        below, so a tab runs off the edge rather than shrink. A spin box can give up some of its
        and stay usable; a label, a button or a pill cannot give up any without losing letters, so
        they keep what they need and the pane scrolls sideways instead. The component count keeps
        its width too: one digit is already as narrow as it goes."""
        for child in window.GetChildren():
            if isinstance(child, (wx.SpinCtrl, wx.SpinCtrlDouble)) and child is not self.spin_ndisc:
                child.SetMinSize(wx.Size(44, child.GetMinSize().height))
            self.allow_shrinking(child)

    def on_wheel(self, event, page):
        lines = -event.GetWheelRotation()/event.GetWheelDelta()*event.GetLinesPerAction()
        x, y = page.GetViewStart()
        page.Scroll(x, y + int(round(lines)))

    def on_step(self, spin):
        """Right-click on a spin box: choose its arrow step, as everywhere else in the suite."""
        steps = (1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 1e4, 1e5, 1e6)
        menu = wx.Menu()
        ids = {}
        for step in steps:
            item = menu.AppendCheckItem(wx.ID_ANY, f"Set step to {step:g}")
            item.Check(step == spin.GetIncrement())
            ids[item.GetId()] = step
        choice = spin.GetPopupMenuSelectionFromUser(menu)
        menu.Destroy()
        if choice not in ids:
            return
        value = spin.GetValue()
        spin.SetDigits(max(2, int(np.ceil(-np.log10(ids[choice])))))
        spin.SetIncrement(ids[choice])
        spin.SetValue(value)

    def on_select(self, event):
        if not self.filling:
            self.mf.set_current(event.GetIndex())

    def on_show(self, event):
        ii = event.GetIndex()
        if not self.filling and 0 <= ii < len(self.mf.Data):
            self.mf.Data[ii]['show'] = self.filelist.IsItemChecked(ii)
            self.mf.plot_panel.update_graph()

    def on_showing(self, every):
        """Tick every dataset, or only the selected one: the ticks alone, so the list stays where
        it is and the plot is drawn once rather than once per row."""
        mf = self.mf
        self.filling = True
        for ii, dd in enumerate(mf.Data):
            dd['show'] = every or ii == mf.current
            self.filelist.CheckItem(ii, dd['show'])
        self.filling = False
        mf.plot_panel.update_graph()

    def on_remove(self, event):
        mf = self.mf
        if mf.current < 0:
            return
        mf.Data.pop(mf.current)
        mf.current = min(mf.current, len(mf.Data) - 1)
        mf.harmonize_units()          # the longest trace may have been the one that just went
        mf.update_everything()

    def on_remove_all(self, event):
        mf = self.mf
        if not mf.Data or wx.MessageBox(f"Drop all {len(mf.Data)} datasets?", "Remove all",
                                        wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
            return
        mf.Data, mf.current = [], -1
        mf.update_everything()

    def on_full(self, event):
        mf = self.mf
        if mf.current < 0:
            return
        t = mf.Data[mf.current]['t']
        mf.Data[mf.current]['win'] = [float(t[0]), float(t[-1])]
        self.update_controls()
        mf.plot_panel.update_graph()

    def on_resample_all(self, event):
        """Hand the current dataset's resampling to every other one. A density per decade means the
        same thing on any time axis, so nothing has to be rescaled on the way."""
        mf = self.mf
        if mf.current < 0:
            return
        dd = mf.Data[mf.current]
        for other in mf.Data:
            other['resample'], other['perdec'], other['rweight'] = (dd['resample'], dd['perdec'],
                                                                    dd['rweight'])
        mf.plot_panel.update_graph()
        note = (f"at {dd['perdec']:g} points per decade, weighted {dd['rweight']}" if dd['resample']
                else "switched off")
        mf.statusbar.SetStatusText(f"Resampling {note} for all {len(mf.Data)} datasets", 0)

    def on_display(self, event=None):
        pp = self.mf.plot_panel
        pp.xmode = 'rate' if self.seg_xmode.GetSelection() else 'time'
        pp.norm = ('raw', 'area', 'fraction')[self.seg_norm.GetSelection()]
        pp.snorm = 'span' if self.seg_snorm.GetSelection() else 'raw'
        pp.showfit = self.chk_showfit.GetValue()
        pp.showdisc = self.chk_showdisc.GetValue()
        pp.legend = self.chk_legend.GetValue()
        pp.sep = self.spin_sep.GetValue()
        pp.stack = self.chk_stack.GetValue()
        pp.stackmode = 'bands' if self.seg_stackmode.GetSelection() else 'lines'
        pp.logt = self.chk_logt.GetValue()
        self.seg_stackmode.Enable(pp.stack)
        pp.update_graph()

    def on_change(self, event=None):
        """Any control changed: read the shared settings into the frame, the rest into the current
        dataset, and redraw. The shared ones can be set with nothing loaded."""
        mf = self.mf
        mf.npts = self.spin_npts.GetValue()
        mf.penalty = self.seg_pen.GetSelection()
        mf.amode = self.seg_amode.GetSelection()
        mf.prob = self.spin_prob.GetValue()
        mf.autoalpha = self.chk_autoalpha.GetValue()
        mf.alpha_lo, mf.alpha_hi = self.spin_alo.GetValue(), self.spin_ahi.GetValue()
        if mf.current >= 0:
            dd = mf.Data[mf.current]
            dd['comp'] = self.seg_comp.GetSelection()
            dd['sum'] = self.chk_sum.GetValue()
            dd['slice'] = self.spin_slice.GetValue()
            dd['win'] = [self.spin_from.GetValue(), self.spin_to.GetValue()]
            dd['recovery'] = self.chk_rec.GetValue()
            dd['offset'] = self.chk_off.GetValue()
            dd['smooth'] = self.chk_smooth.GetValue()
            dd['swin'] = self.spin_swin.GetValue()
            dd['resample'] = self.chk_resample.GetValue()
            dd['perdec'] = self.spin_perdec.GetValue()
            dd['rweight'] = 'equal' if self.seg_rweight.GetSelection() else 'count'
            dd['autogrid'] = self.chk_auto.GetValue()
            dd['Tmin'] = self.spin_tmin.GetValue()
            dd['Tmax'] = self.spin_tmax.GetValue()
            dd['logalpha'] = self.spin_alpha.GetValue()
            dd['ndisc'] = self.spin_ndisc.GetValue()
            dd['fixT'] = self.chk_fixt.GetValue()
        self.enable()
        mf.plot_panel.update_graph()

    def on_fit(self, event):
        mf = self.mf
        if mf.current >= 0:
            mf.discrete(mf.Data[mf.current])
            mf.plot_panel.update_graph()

    def on_fit_all(self, event):
        """The same fit on every loaded dataset, failures collected instead of interrupting."""
        mf = self.mf
        if not mf.Data:
            return
        wait = wx.BusyCursor()
        failed = []
        # the whole session is fitted the same way, so the components compare across datasets
        count, hold = self.spin_ndisc.GetValue(), self.chk_fixt.GetValue()
        for dd in mf.Data:
            dd['ndisc'], dd['fixT'] = count, hold
            mf.analyse(dd)
            if not mf.discrete(dd, report=False):
                failed.append(dd['title'])
        del wait
        mf.plot_panel.update_graph()
        mf.statusbar.SetStatusText(f"Fitted {len(mf.Data) - len(failed)} of {len(mf.Data)} datasets "
                                   f"with {count} component{'s' if count != 1 else ''}"
                                   + (", T held" if hold else "")
                                   + (f"; {len(failed)} did not converge" if failed else ""), 0)
        if failed:
            wx.MessageBox("These datasets did not converge:\n\n" + "\n".join(failed[:20])
                          + ("\n..." if len(failed) > 20 else ""), "Fit all", wx.OK | wx.ICON_WARNING)

    def on_clear(self, event):
        mf = self.mf
        if mf.current >= 0:
            mf.Data[mf.current]['disc'] = None
            mf.plot_panel.update_graph()

    def enable(self):
        """Grey out what does not apply: no dataset, 1D data, Auto grid, the mode alpha comes from."""
        mf = self.mf
        dd = mf.Data[mf.current] if mf.current >= 0 else None
        for w in (self.seg_comp, self.spin_from, self.spin_to, self.btn_full, self.chk_rec, self.chk_off,
                  self.chk_smooth, self.chk_resample, self.btn_resall, self.chk_auto, self.spin_ndisc,
                  self.chk_fixt, self.btn_fit,
                  self.btn_clear, self.btn_remove, self.btn_showonly):
            w.Enable(dd is not None)
        for w in (self.btn_fitall, self.btn_showall, self.btn_removeall):
            w.Enable(bool(mf.Data))
        self.btn_savedisc.Enable(any(d['disc'] is not None for d in mf.Data))
        twod = dd is not None and dd['data'].ndim == 2
        self.chk_sum.Enable(twod)
        self.spin_slice.Enable(twod and not dd['sum'])
        self.spin_swin.Enable(dd is not None and dd['smooth'])
        for w in (self.spin_perdec, self.seg_rweight):
            w.Enable(dd is not None and dd['resample'])
        for w in (self.spin_tmin, self.spin_tmax):
            w.Enable(dd is not None and not dd['autogrid'])
        self.seg_stackmode.Enable(self.chk_stack.GetValue())
        self.spin_alpha.Enable(dd is not None and mf.amode == 0)
        self.spin_prob.Enable(mf.amode == 3)
        for w in (self.spin_alo, self.spin_ahi):
            w.Enable(not mf.autoalpha and (mf.amode or mf.plot_panel.show['lcurve']))
        self.chk_autoalpha.Enable(bool(mf.amode) or mf.plot_panel.show['lcurve'])

    def edit_param(self, ii):
        """Type the parameter straight into its cell: a box over the list, Enter or clicking away
        keeps what was typed, Escape leaves the old value. Anything that is not a number clears it,
        and an empty cell simply keeps that dataset out of the rates window."""
        mf = self.mf
        if not 0 <= ii < len(mf.Data):
            return
        value = mf.Data[ii].get('param', float('nan'))
        rect = self.filelist.GetItemRect(ii)
        editor = wx.TextCtrl(self.filelist, value="" if not np.isfinite(value) else f"{value:g}",
                             pos=(rect.x + self.filelist.GetColumnWidth(0), rect.y),
                             size=(self.filelist.GetColumnWidth(1), rect.height),
                             style=wx.TE_PROCESS_ENTER | wx.BORDER_SIMPLE)
        editor.SelectAll()
        editor.SetFocus()
        done = []                 # the kill-focus that follows Enter would otherwise run twice

        def finish(event, keep=True):
            if done:
                return
            done.append(True)
            text = editor.GetValue().strip()
            wx.CallAfter(editor.Destroy)
            if keep:
                try:
                    mf.Data[ii]['param'] = float(text) if text else float('nan')
                except ValueError:
                    print(f"[EPRelax] '{text}' is not a number, so '{mf.Data[ii]['title']}' is left "
                          f"without a {mf.rates['name'].split(',')[0].strip().lower()}")
                    mf.Data[ii]['param'] = float('nan')
                self.show_param(ii)        # that one cell, not the whole list
                if mf.rate_window is not None:
                    mf.rate_window.update_graph()

        editor.Bind(wx.EVT_TEXT_ENTER, finish)
        editor.Bind(wx.EVT_KILL_FOCUS, finish)
        editor.Bind(wx.EVT_CHAR, lambda e: finish(e, keep=False) if e.GetKeyCode() == wx.WXK_ESCAPE
                    else e.Skip())

    def show_param(self, ii):
        """That dataset's parameter into its own cell. An empty cell is one that has none."""
        value = self.mf.Data[ii].get('param', float('nan'))
        self.filelist.SetItem(ii, 1, "" if not np.isfinite(value) else f"{value:g}")

    def param_column(self):
        """The parameter column's heading, which is whatever the rates window calls it."""
        column = self.filelist.GetColumn(1)
        column.SetText(self.mf.rates['name'])
        self.filelist.SetColumn(1, column)

    def update_list(self):
        """Rebuild the dataset list, each row written in that dataset's trace colour. Only for when
        the datasets themselves change: a tick, a colour or a parameter is written into its own row
        instead, so the list neither scrolls back to the top nor costs a redraw per row."""
        mf = self.mf
        top = self.filelist.GetTopItem()
        self.filling = True
        self.filelist.DeleteAllItems()
        self.param_column()
        for ii, dd in enumerate(mf.Data):
            self.filelist.InsertItem(ii, dd['title'])
            self.show_param(ii)
            self.filelist.CheckItem(ii, dd['show'])
            self.filelist.SetItemTextColour(ii, wx.Colour(to_hex(mf.trace_color(dd))))
        if 0 <= mf.current < len(mf.Data):
            self.filelist.Select(mf.current)
        self.filling = False
        if 0 < top < self.filelist.GetItemCount():
            self.filelist.EnsureVisible(self.filelist.GetItemCount() - 1)
            self.filelist.EnsureVisible(top)      # back to the row the list was showing from

    def update_controls(self):
        """Fill every control: the shared ones from the frame, the rest from the current dataset."""
        mf = self.mf
        self.spin_npts.SetValue(mf.npts)
        self.seg_pen.SetSelection(mf.penalty)
        self.seg_amode.SetSelection(mf.amode)
        self.spin_prob.SetValue(mf.prob)
        self.chk_autoalpha.SetValue(mf.autoalpha)
        self.spin_alo.SetValue(mf.alpha_lo)
        self.spin_ahi.SetValue(mf.alpha_hi)
        if mf.current < 0:
            self.lbl_window.SetLabel("")
            self.enable()
            return
        dd = mf.Data[mf.current]
        self.seg_comp.SetSelection(dd['comp'])
        self.chk_sum.SetValue(dd['sum'])
        self.spin_slice.SetRange(0, max(dd['data'].shape[1] - 1, 0) if dd['data'].ndim == 2 else 0)
        self.spin_slice.SetValue(min(dd['slice'], self.spin_slice.GetMax()))
        # enough decimals for the unit the session settled on: a short trace beside a long one
        # reads in thousandths of it, and three places would show every limit as 0.000
        places = lambda v: int(min(max(3.0, np.ceil(-np.log10(abs(v))) + 2 if v else 3.0), 10.0))
        step = max(abs(float(dd['t'][-1] - dd['t'][0]))/100, 1e-12)
        for spin, value in ((self.spin_from, dd['win'][0]), (self.spin_to, dd['win'][1]),
                            (self.spin_swin, dd['swin'])):
            spin.SetDigits(places(step))
            spin.SetIncrement(step)
            spin.SetValue(value)
        self.lbl_window.SetLabel(f"of {dd['t'][0]:.4g} … {dd['t'][-1]:.4g} {dd['tunit']}")
        self.chk_rec.SetValue(dd['recovery'])
        self.chk_off.SetValue(dd['offset'])
        self.chk_smooth.SetValue(dd['smooth'])
        self.chk_resample.SetValue(dd['resample'])
        self.spin_perdec.SetValue(dd['perdec'])
        self.seg_rweight.SetSelection(1 if dd['rweight'] == 'equal' else 0)
        self.chk_auto.SetValue(dd['autogrid'])
        for spin, value in ((self.spin_tmin, dd['Tmin']), (self.spin_tmax, dd['Tmax'])):
            spin.SetDigits(places(dd['Tmin']))
            spin.SetValue(value)
        self.spin_alpha.SetValue(dd['logalpha'])
        self.spin_ndisc.SetValue(dd['ndisc'])
        self.chk_fixt.SetValue(dd['fixT'])
        self.enable()

    def update_stats(self):
        """What the last decomposition and the last discrete fit came to, for the current dataset."""
        mf = self.mf
        self.enable()            # a fit just finishing is what wakes the Save button up
        rate = mf.plot_panel.xmode == 'rate'
        r = mf.Data[mf.current]['res'] if mf.current >= 0 else None
        if r is None:
            self.lbl_stat.SetLabel("")
            self.lbl_disc.SetLabel("")
            self.disclist.DeleteAllItems()
            while self.disclist.GetColumnCount() > 1:
                self.disclist.DeleteColumn(1)
            for row, name in enumerate(self.DISC_ROWS):
                self.disclist.InsertItem(row, name)
            self.Layout()
            return
        dd = mf.Data[mf.current]
        unit = f"1/{dd['tunit']}" if rate else dd['tunit']
        name = "rate" if rate else "T"
        self.spin_alpha.SetValue(float(np.log10(r['alpha'])))      # the automatic modes report back here
        if dd['autogrid']:
            self.spin_tmin.SetDigits(int(min(max(3.0, np.ceil(-np.log10(r['T'][0])) + 2), 10.0)))
            self.spin_tmax.SetDigits(self.spin_tmin.GetDigits())
            self.spin_tmin.SetValue(float(r['T'][0]))
            self.spin_tmax.SetValue(float(r['T'][-1]))
        total = float(r['a'].sum())
        mean = float((r['a']*r['T']).sum()/total) if total > 0 else float('nan')
        peaks = [j for j in range(1, r['T'].size - 1) if r['a'][j] > r['a'][j - 1] and r['a'][j] >= r['a'][j + 1]]
        peaks.sort(key=lambda j: -r['a'][j])
        show = lambda T: 1.0/T if rate else T
        self.lbl_stat.SetLabel(
            f"alpha = {r['alpha']:.4g}   d.o.f. = {r['df']:.1f} of {r['T'].size}\n"
            f"RMS residual = {r['rms']:.4g}  ({100*r['rms']/max(abs(r['y']).max(), 1e-300):.2f} % of peak)\n"
            f"offset = {r['c']:.4g}   amplitude = {total:.4g}\n"
            f"mean {name} = {show(mean):.4g} {unit}\n"
            + "peaks at " + (", ".join(f"{show(r['T'][j]):.4g}" for j in peaks[:4]) + f" {unit}" if peaks else "none"))

        self.disclist.DeleteAllItems()
        while self.disclist.GetColumnCount() > 1:
            self.disclist.DeleteColumn(1)
        for row, text in enumerate((f"{'k' if rate else 'T'}, {unit}", "± 95%", "a", "± 95%",
                                    f"iLT low, {unit}", f"iLT high, {unit}", "note")):
            self.disclist.InsertItem(row, text)
        d = dd['disc']
        if d is None:
            self.lbl_disc.SetLabel("")
        else:
            mark = "~" if d.get('Thold') else ""
            for i in range(len(d['T'])):
                # dT/T is what inverting preserves, so the rate's interval is dT/T^2
                value = show(d['T'][i])
                err = d['dT'][i]/d['T'][i]**2 if rate else d['dT'][i]
                lo, hi = d['T'][i]/d['wide'][i], d['T'][i]*d['wide'][i]
                if rate:
                    lo, hi = 1.0/hi, 1.0/lo        # the slow end of the peak is the small rate
                # a component that sits on no peak has no peak width: the cell stays empty
                edge = lambda v: f"{v:.5g}" if np.isfinite(v) else ""
                self.disclist.InsertColumn(i + 1, f"{i + 1}", width=76)
                for row, text in enumerate((f"{value:.5g}", f"{mark}{err:.3g}", f"{d['a'][i]:.5g}",
                                            f"{d['da'][i]:.3g}", edge(lo),
                                            edge(hi), d['note'][i])):
                    self.disclist.SetItem(row, i + 1, text)
                self.disclist.SetColumnWidth(i + 1, wx.LIST_AUTOSIZE)      # the note row can be long
                self.disclist.SetColumnWidth(i + 1, max(76, self.disclist.GetColumnWidth(i + 1) + 8))
            stale = ("\nThe trace has changed since this fit: the curve follows the parameters, "
                     "but the numbers are the old window's. Fit again." if d.get('sig') != dd['sig']
                     else "")
            if d.get('Thold'):
                stale = ("\nT was held, so its ~ interval is the one it would have had if it were "
                         "free at that value." + stale)
            self.lbl_disc.SetLabel(f"offset = {d['c']:.4g}   RMS residual = {d['rms']:.4g}\n"
                                   f"95 % intervals from the fit's trust region, "
                                   f"{d['dof']} degrees of freedom" + stale)
        self.Layout()


# --------------------------------------------------------------------------- #
# Rates against a parameter
# --------------------------------------------------------------------------- #
class RateParameterFrame(wx.Frame):
    """Separate window: the relaxation rates of the discrete fits against a number given to each
    dataset in the Datasets list -- a temperature, a concentration, a field, whatever the series was
    measured along. Every point carries the 95 % interval its component was fitted with, held time
    constants included, and the set can be fitted with an equation typed in by hand.

    The points are read from the fits as they stand, so refitting in the main window moves them
    here too. Every dataset that has both a discrete fit and a parameter takes part, ticked or not:
    the Datasets ticks say what the trace plot draws, which is a different question from this one."""

    def __init__(self, parent, mainWindow):
        super().__init__(parent, title="Rates against a parameter", size=(1060, 700))
        self.mf = mainWindow
        self.steps = {}           # fit parameter -> arrow step chosen by right-clicking its box
        self.boxes = {}           # fit parameter -> (value box, hold tick)
        self.result = None        # the last fit: its parameters, their intervals, RMSD, R^2 and n
        self.touched = False      # the values have been typed since that fit, so it is not it

        splitter = wx.SplitterWindow(self, style=wx.SP_3D | wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(150)
        left = scrolled.ScrolledPanel(splitter, style=wx.BORDER_NONE)
        right = wx.Panel(splitter)
        cards = wx.BoxSizer(wx.VERTICAL)
        rates = self.mf.rates

        # ---------------- Points card: what is drawn against what
        card = theme.Card(left, "Points")
        self.txt_name = wx.TextCtrl(card, value=rates['name'], style=wx.TE_PROCESS_ENTER)
        self.txt_name.SetToolTip("What the parameter is called: the heading of its column in the Datasets\n"
                                 "list and the label of the axis here. 'Temperature, K' and the like, since\n"
                                 "the unit is yours to state -- the program never converts it.")
        self.txt_name.Bind(wx.EVT_TEXT_ENTER, self.on_name)
        self.txt_name.Bind(wx.EVT_KILL_FOCUS, self.on_name)
        card.body.Add(self.labelled(card, "Name:", self.txt_name), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_series = SegmentedPill(card, ["By T", "By ampl.", "All in one"])
        self.seg_series.SetToolTip("Which component of one dataset belongs on the same series as which\n"
                                   "component of the next.\n"
                                   "By T: component 1 is the fastest everywhere, 2 the next, as the Discrete\n"
                                   "fit table has them. The honest choice when the components keep their\n"
                                   "order along the series.\n"
                                   "By ampl.: series 1 is the tallest component of each dataset, whatever its\n"
                                   "T. Follows the dominant process, but jumps where two amplitudes cross.\n"
                                   "All in one: every component as one set of points, matched up with\n"
                                   "nothing -- which is what to draw when the components are not comparable.")
        card.body.Add(self.labelled(card, "Series:", self.seg_series), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_y = SegmentedPill(card, ["Rate", "Time const."])
        self.seg_y.SetToolTip("The vertical axis: the relaxation rate 1/T, or the time constant itself.\n"
                              "The interval follows the quantity -- dT/T^2 for a rate, dT for a time.\n"
                              "What the fit works in is set in the Fit card and is free of this: a law\n"
                              "fitted to the rates can be looked at on a time-constant axis.")
        card.body.Add(self.labelled(card, "Y axis:", self.seg_y), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_x = SegmentedPill(card, ["Parameter", "1 / parameter"])
        self.seg_x.SetToolTip("The horizontal axis: the parameter as typed, or one over it -- 1/T for an\n"
                              "Arrhenius plot, say. The fit equation always sees the parameter itself,\n"
                              "so flipping the axis re-draws the same fit rather than changing it.")
        card.body.Add(self.labelled(card, "X axis:", self.seg_x), 0, wx.EXPAND | wx.ALL, 3)
        self.chk_logx = PillCheckBox(card, "Log x")
        self.chk_logx.SetToolTip("Logarithmic horizontal axis. A parameter at or below zero is left out of\n"
                                 "the drawing while it is on, and said so in the line under the plot.")
        self.chk_logy = PillCheckBox(card, "Log y")
        self.chk_logy.SetToolTip("Logarithmic vertical axis, which is usually what a rate spanning decades\n"
                                 "wants.")
        self.chk_bars = PillCheckBox(card, "Error bars")
        self.chk_bars.SetToolTip("The 95 % interval of each component, from the discrete fit's own trust\n"
                                 "region. A held time constant has one too -- the interval it would have had\n"
                                 "if it had been free at that value -- and it is drawn like any other.")
        for w in (self.chk_logx, self.chk_logy, self.chk_bars):
            card.body.Add(w, 0, wx.ALL, 3)
        for w in (self.seg_series, self.seg_y, self.seg_x):
            w.Bind(wx.EVT_RADIOBOX, self.on_display)
        for w in (self.chk_logx, self.chk_logy, self.chk_bars):
            w.Bind(wx.EVT_CHECKBOX, self.on_display)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # ---------------- Fit card: an equation typed in, and its parameters
        card = theme.Card(left, "Fit")
        self.txt_eq = wx.TextCtrl(card, value=rates['equation'], style=wx.TE_PROCESS_ENTER)
        self.txt_eq.SetToolTip("y as a function of x, written out: 'a*exp(-b/x) + c', 'a + b*x', 'a*x**n'.\n"
                               "x is the parameter itself, never its inverse, whatever the axis shows.\n"
                               "Every other name is a parameter of the fit and gets a box below.\n"
                               "exp, log, log10, sqrt, sin, cos, tan, sinh, cosh, tanh, abs, pow, pi and e\n"
                               "are available; ^ works as **. Press Enter to read the equation in.")
        self.txt_eq.Bind(wx.EVT_TEXT_ENTER, self.on_equation)
        self.txt_eq.Bind(wx.EVT_KILL_FOCUS, self.on_equation)
        card.body.Add(self.labelled(card, "y =", self.txt_eq), 0, wx.EXPAND | wx.ALL, 3)
        self.seg_fit = SegmentedPill(card, ["Rate", "Time const."])
        self.seg_fit.SetToolTip("What y is to the fit: the rate 1/T, or the time constant itself. This is\n"
                                "the quantity the equation describes, the one the residuals and the\n"
                                "weights are taken in, and it is free of what the Y axis shows -- fit an\n"
                                "Arrhenius law to the rates and look at it in T if that reads better.\n"
                                "Where the two differ the curve is drawn as one over itself, so it still\n"
                                "lies on the points, and the legend says which way round it was fitted.")
        self.seg_fit.Bind(wx.EVT_RADIOBOX, self.on_display)
        card.body.Add(self.labelled(card, "Fit in:", self.seg_fit), 0, wx.EXPAND | wx.ALL, 3)
        self.choice_target = wx.Choice(card, choices=["All points"])
        self.choice_target.SetSelection(0)
        self.choice_target.SetToolTip("Which points the fit goes through: all of them, or one series on its\n"
                                      "own. One fit is kept at a time.")
        self.choice_target.Bind(wx.EVT_CHOICE, self.on_display)
        card.body.Add(self.labelled(card, "Fit to:", self.choice_target), 0, wx.EXPAND | wx.ALL, 3)
        self.parsizer = wx.FlexGridSizer(cols=3, vgap=2, hgap=6)
        self.parsizer.AddGrowableCol(1)
        card.body.Add(self.parsizer, 0, wx.EXPAND | wx.ALL, 3)
        self.chk_weight = PillCheckBox(card, "Weight by error bars")
        self.chk_weight.SetToolTip("Give each point the weight its own interval earns: sigma is the interval\n"
                                   "divided by the Student t the discrete fit used, so a component the trace\n"
                                   "barely pins down pulls less than one it pins down well, and the\n"
                                   "parameter intervals below mean what they say.\n"
                                   "Untick for a straight unweighted fit, where every point counts the same.")
        card.body.Add(self.chk_weight, 0, wx.ALL, 3)
        self.chk_weight.Bind(wx.EVT_CHECKBOX, self.on_display)
        self.chk_draw = PillCheckBox(card, "Draw the equation")
        self.chk_draw.SetToolTip("Draw the equation over the points with the values in the boxes, fitted\n"
                                 "or typed. Turn a value and the curve follows at once, which is how a\n"
                                 "hand-typed equation is brought near enough to fit from.\n"
                                 "The legend says 'fit' while the values are the ones the fit left, and\n"
                                 "'equation' once they have been typed over.")
        self.chk_draw.Bind(wx.EVT_CHECKBOX, self.on_display)
        card.body.Add(self.chk_draw, 0, wx.ALL, 3)
        self.chk_terms = PillCheckBox(card, "Term by term")
        self.chk_terms.SetToolTip("Draw each term of the equation on its own as well, dashed and in its\n"
                                  "own colour, with the same values as the whole curve -- so a sum of two\n"
                                  "processes can be seen as the two it is made of.\n"
                                  "The equation is broken at every + and - that is not inside brackets,\n"
                                  "so 'a*exp(-b/x) + c' draws as 'a*exp(-b/x)' and 'c'. The sign goes with\n"
                                  "the term, and the terms always add back up to the whole curve.\n"
                                  "An equation of one term has nothing to break apart and draws nothing.")
        self.chk_terms.Bind(wx.EVT_CHECKBOX, self.on_display)
        card.body.Add(self.chk_terms, 0, wx.ALL, 3)
        self.btn_fit = PillButton(card, "Fit", color_key="pill_accent_2", draw_badge=theme.draw_play_badge)
        self.btn_fit.SetToolTip("Fit the equation to the chosen points, starting from the values in the\n"
                                "boxes. A hand-typed equation has no automatic starting guess, so put the\n"
                                "boxes within an order of magnitude first.")
        self.btn_fit.Bind(wx.EVT_BUTTON, self.on_fit)
        self.btn_clear = PillButton(card, "Clear", color_key="pill_calm", draw_badge=theme.draw_minus_badge)
        self.btn_clear.SetToolTip("Drop the fitted curve. The typed values stay where they are.")
        self.btn_clear.Bind(wx.EVT_BUTTON, self.on_clear)
        self.btn_save = PillButton(card, "Save", color_key="pill_load", draw_badge=theme.draw_textlines_badge)
        self.btn_save.SetToolTip("Write the points, their intervals and the fit to one CSV table.")
        self.btn_save.Bind(wx.EVT_BUTTON, self.on_save)
        widest = max(b.GetBestSize().width for b in (self.btn_fit, self.btn_clear, self.btn_save))
        row = wx.BoxSizer(wx.HORIZONTAL)
        for b in (self.btn_fit, self.btn_clear, self.btn_save):
            b.SetMinSize(wx.Size(widest, b.GetBestSize().height))
            b.SetMaxSize(wx.Size(widest, -1))
            row.Add(b, 0, wx.RIGHT, 4)
        card.body.Add(row, 0, wx.ALL, 3)
        self.lbl_fit = wx.StaticText(card, label="")
        card.body.Add(self.lbl_fit, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        left.SetSizer(cards)
        left.SetupScrolling(scroll_x=True)
        self.figure = Figure(figsize=(5, 4), dpi=100, layout='constrained')
        self.canvas = FigureCanvasWxAgg(right, -1, self.figure)
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        self.lbl_stat = wx.StaticText(right, label="")
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.canvas, 1, wx.EXPAND)
        box.Add(self.toolbar, 0, wx.EXPAND)
        box.Add(self.lbl_stat, 0, wx.EXPAND | wx.ALL, 6)
        right.SetSizer(box)
        splitter.SplitVertically(left, right, 300)
        splitter.SetSashGravity(0)
        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(splitter, 1, wx.EXPAND)
        self.SetSizer(outer)

        self.redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self.redraw_debouncer.on_size)
        self.Bind(wx.EVT_CLOSE, lambda e: (setattr(self.mf, 'rate_window', None), self.Destroy()))
        self.Centre()
        self.pull()
        self.apply_theme()

    # -----------------------------------------------------------------------
    def labelled(self, card, text, widget):
        """A row label of one width down the card, with its control beside it."""
        lbl = wx.StaticText(card, label=text)
        lbl.SetMinSize(wx.Size(max(62, lbl.GetBestSize().width), -1))
        box = wx.BoxSizer(wx.HORIZONTAL)
        box.Add(lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 4)
        box.Add(widget, 1, wx.EXPAND)
        return box

    def apply_theme(self):
        theme.apply_theme_to_window(self, theme.get_theme())
        self.update_graph()

    def pull(self):
        """Every control from the settings the main window keeps, so a session opens as it was saved."""
        rates = self.mf.rates
        self.txt_name.SetValue(rates['name'])
        self.txt_eq.SetValue(rates['equation'])
        self.seg_series.SetSelection(rates['series'])
        self.seg_y.SetSelection(rates['ymode'])
        self.seg_x.SetSelection(rates['xmode'])
        self.chk_logx.SetValue(rates['logx'])
        self.chk_logy.SetValue(rates['logy'])
        self.chk_bars.SetValue(rates['bars'])
        self.chk_weight.SetValue(rates['weight'])
        self.chk_draw.SetValue(rates['draw'])
        self.chk_terms.SetValue(rates['terms'])
        self.chk_terms.Enable(rates['draw'])
        self.seg_fit.SetSelection(rates['fitmode'])
        self.touched = False
        self.rebuild_params()

    def on_name(self, event):
        """The parameter's name: it heads its column in the Datasets list and labels the axis here."""
        self.mf.rates['name'] = self.txt_name.GetValue().strip() or "Parameter"
        self.mf.control_panel.param_column()
        self.update_graph()
        event.Skip()

    def on_display(self, event=None):
        """A switch moved: read them all back and redraw."""
        rates = self.mf.rates
        rates['series'] = self.seg_series.GetSelection()
        rates['ymode'] = self.seg_y.GetSelection()
        rates['xmode'] = self.seg_x.GetSelection()
        rates['logx'] = self.chk_logx.GetValue()
        rates['logy'] = self.chk_logy.GetValue()
        rates['bars'] = self.chk_bars.GetValue()
        rates['weight'] = self.chk_weight.GetValue()
        rates['draw'] = self.chk_draw.GetValue()
        rates['terms'] = self.chk_terms.GetValue()
        self.chk_terms.Enable(rates['draw'])
        if rates['fitmode'] != self.seg_fit.GetSelection():
            rates['fitmode'] = self.seg_fit.GetSelection()
            self.result, self.touched = None, False   # the old fit was of the other quantity
        rates['target'] = max(self.choice_target.GetSelection(), 0)
        self.update_graph()

    def on_equation(self, event):
        """A new equation: keep the parameters it still has, start the new ones at 1."""
        self.mf.rates['equation'] = self.txt_eq.GetValue().strip()
        self.result, self.touched = None, False
        self.rebuild_params()
        self.update_graph()
        event.Skip()

    def on_step(self, spin, key):
        """Right-click on a value box: choose its arrow step, as everywhere else in the suite."""
        menu = wx.Menu()
        ids = {}
        for step in (1e-6, 1e-4, 1e-2, 0.1, 1, 10, 100, 1e4, 1e6):
            item = menu.AppendCheckItem(wx.ID_ANY, f"Set step to {step:g}")
            item.Check(step == spin.GetIncrement())
            ids[item.GetId()] = step
        choice = spin.GetPopupMenuSelectionFromUser(menu)
        menu.Destroy()
        if choice not in ids:
            return
        value = spin.GetValue()
        self.steps[key] = ids[choice]
        spin.SetDigits(max(2, int(np.ceil(-np.log10(ids[choice])))))
        spin.SetIncrement(ids[choice])
        spin.SetValue(value)

    def equation(self):
        """The equation as eval() can take it, and the names that are its parameters. x is the
        parameter the points are drawn against and is never one of them."""
        work, _tokens, free = extract_refs_and_freevars(self.mf.rates['equation'], (), FIT_FUNCS)
        return work, sorted(n for n in free if n != 'x')

    def rebuild_params(self):
        """One row of boxes per name in the equation, values kept across a re-reading of it."""
        rates = self.mf.rates
        names = self.equation()[1]
        rates['pars'] = {n: float(rates['pars'].get(n, 1.0)) for n in names}
        rates['hold'] = {n: bool(rates['hold'].get(n, False)) for n in names}
        self.parsizer.Clear(delete_windows=True)
        self.boxes = {}
        card = self.txt_eq.GetParent()
        for name in names:
            value = rates['pars'][name]
            step = self.steps.get(name) or (10.0**(np.floor(np.log10(abs(value))) - 1) if value else 0.1)
            spin = wx.SpinCtrlDouble(card, min=-1e15, max=1e15, inc=float(step), initial=value,
                                     size=(110, -1))
            spin.SetDigits(max(2, int(np.ceil(-np.log10(step)))))
            spin.SetValue(value)
            spin.SetToolTip(f"Starting value of '{name}', and where the fit leaves it.\n"
                            "Right-click to set the arrow step.")
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_param)
            theme.bind_spin_enter(spin, self.on_param)
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                spin.Bind(evt, lambda e, s=spin, k=name: self.on_step(s, k))
            hold = PillCheckBox(card, "Hold")
            hold.SetValue(rates['hold'][name])
            hold.SetToolTip(f"Keep '{name}' at the value in the box instead of fitting it.")
            hold.Bind(wx.EVT_CHECKBOX, self.on_param)
            self.parsizer.Add(wx.StaticText(card, label=name), 0, wx.ALIGN_CENTER_VERTICAL)
            self.parsizer.Add(spin, 1, wx.EXPAND)
            self.parsizer.Add(hold, 0, wx.ALIGN_CENTER_VERTICAL)
            self.boxes[name] = (spin, hold)
        self.parsizer.Layout()
        card.Layout()
        self.Layout()

    def on_param(self, event=None):
        """A value or a hold tick typed or turned: the curve follows at once, and the fit the
        numbers below describe is no longer the one being drawn."""
        rates = self.mf.rates
        for name, (spin, hold) in self.boxes.items():
            rates['pars'][name] = float(spin.GetValue())
            rates['hold'][name] = bool(hold.GetValue())
        self.touched = True
        self.update_graph()

    def push_params(self):
        """Fitted values back into their boxes, with decimals that suit the magnitude."""
        for name, (spin, hold) in self.boxes.items():
            value = self.mf.rates['pars'][name]
            step = self.steps.get(name) or (10.0**(np.floor(np.log10(abs(value))) - 2) if value else 0.01)
            spin.SetDigits(max(int(np.clip(4 - np.floor(np.log10(abs(value))) if value else 4, 2, 12)),
                               int(np.ceil(-np.log10(step)))))
            spin.SetIncrement(float(step))
            spin.SetValue(value)

    # -----------------------------------------------------------------------
    def points(self):
        """Every component that has both a discrete fit and a parameter, gathered into series, as
        {series: {'x', 'y', 'dy', 'held', 'dof', 'name'}}. Also the time unit the numbers are in and
        whatever had to be left out, which the line under the plot reports."""
        mf, rates = self.mf, self.mf.rates
        out, unit, mixed, nofit, noparam = {}, '', set(), 0, 0
        for dd in mf.Data:
            if dd['disc'] is None:
                nofit += 1
                continue
            if not np.isfinite(dd.get('param', float('nan'))):
                noparam += 1
                continue
            d = dd['disc']
            unit = unit or dd['tunit']
            if dd['tunit'] != unit:
                mixed.add(dd['tunit'])
            # by T they are already in order; by amplitude the tallest leads
            order = np.argsort(d['a'])[::-1] if rates['series'] == 1 else range(len(d['T']))
            for rank, i in enumerate(order):
                T, dT = float(d['T'][i]), float(d['dT'][i])
                if not np.isfinite(T) or T <= 0:
                    continue
                # the fitted time constant as it stands; the rate, where it is wanted, is made of
                # it by whoever asks -- what is drawn and what is fitted need not be the same
                row = out.setdefault(0 if rates['series'] == 2 else rank,
                                     {'x': [], 'T': [], 'dT': [], 'held': [], 'dof': [], 'name': []})
                row['x'].append(float(dd['param']))
                row['T'].append(T)
                row['dT'].append(dT)
                row['held'].append(bool(d.get('Thold')))
                row['dof'].append(int(d.get('dof', 1)))
                row['name'].append(dd['title'])
        warn = []
        if mixed:
            warn.append(f"{', '.join(sorted(mixed))} mixed with {unit}: the numbers are each file's own")
        if nofit:
            warn.append(f"{nofit} not fitted")
        if noparam:
            warn.append(f"{noparam} without a {rates['name'].split(',')[0].strip().lower()}")
        for key, row in out.items():
            order = np.argsort(np.asarray(row['x'], dtype=float))
            for field in ('x', 'T', 'dT', 'dof'):
                row[field] = np.asarray(row[field], dtype=float)[order]
            row['held'] = [row['held'][j] for j in order]
            row['name'] = [row['name'][j] for j in order]
        return out, unit, warn

    def terms(self):
        """The equation's top-level terms, each carrying the sign in front of it, as text eval()
        can take. A + or a - counts as a break only outside brackets and only where something has
        gone before it: the sign of the first term, the minus of x**-2 and the one inside 1e-5 are
        all part of their term, not a break between two."""
        work = self.equation()[0]
        out, depth, start = [], 0, 0
        for k, ch in enumerate(work):
            if ch in '([':
                depth += 1
            elif ch in ')]':
                depth -= 1
            elif depth == 0 and ch in '+-' and k > start:
                before = work[:k].rstrip()
                exponent = before[-1:] in ('e', 'E') and (before[-2:-1].isdigit() or before[-2:-1] == '.')
                if before and before[-1] not in '+-*/^(,' and not exponent:
                    out.append(work[start:k])
                    start = k
        out.append(work[start:])
        return [t.strip() for t in out if t.strip() not in ('', '+', '-')]

    def curve(self, x, values, work=None):
        """The typed equation at x, with these parameter values, or one term of it. NaN all through
        where it will not evaluate -- a misspelt name, a log of a negative number -- with the reason
        on the console."""
        work = self.equation()[0] if work is None else work
        env = dict(FIT_FUNCS)
        env.update({k: float(v) for k, v in values.items()})
        env['x'] = np.asarray(x, dtype=float)
        try:
            with np.errstate(all='ignore'):
                out = np.asarray(eval(work, {"__builtins__": {}}, env), dtype=float)
        except Exception as fault:
            print(f"[EPRelax] '{self.mf.rates['equation']}' could not be evaluated: "
                  f"{type(fault).__name__}: {fault}")
            return np.full(np.shape(x), np.nan)
        return np.broadcast_to(out, np.shape(x)).astype(float)

    def quantity(self, row, mode):
        """One series as (y, its 95 % interval) in the quantity asked for: the rate 1/T, or T
        itself. dT/T is what inverting preserves, so the rate's interval is dT/T^2."""
        T, dT = row['T'], row['dT']
        return (1.0/T, dT/T**2) if mode == 0 else (T, dT)

    def chosen(self, series, mode):
        """The points 'Fit to' picks out -- all of them, or one series -- as (keys, x, y, dy, dof),
        in the quantity asked for. The fit works on these, and so do its residual numbers."""
        rates = self.mf.rates
        keys = sorted(series) if rates['target'] == 0 else [k for k in series if k + 1 == rates['target']]
        if not keys:
            return [], *(np.array([]) for _ in range(4))
        pieces = [self.quantity(series[k], mode) for k in keys]
        return (keys, np.concatenate([series[k]['x'] for k in keys]),
                np.concatenate([p[0] for p in pieces]), np.concatenate([p[1] for p in pieces]),
                np.concatenate([series[k]['dof'] for k in keys]))

    def on_fit(self, event):
        """Fit the chosen points with the typed equation. x is the parameter itself, whatever the
        axis is showing, so flipping the axis never changes what was fitted."""
        series, unit, _warn = self.points()
        rates = self.mf.rates
        keys, x, y, dy, dof = self.chosen(series, rates['fitmode'])
        if not keys:
            wx.MessageBox("Those points are not there to fit.", "Fit", wx.OK | wx.ICON_INFORMATION)
            return
        names = self.equation()[1]
        free = [n for n in names if not rates['hold'][n]]
        if not free:
            wx.MessageBox("Every parameter is held, so there is nothing to fit.", "Fit",
                          wx.OK | wx.ICON_INFORMATION)
            return
        if x.size <= len(free):
            wx.MessageBox(f"{x.size} point{'s' if x.size != 1 else ''} cannot pin down "
                          f"{len(free)} parameters.", "Fit", wx.OK | wx.ICON_INFORMATION)
            return
        sigma = None
        if rates['weight']:
            # the stored interval is t(0.975, dof) standard errors wide, so this is the error itself
            sigma = dy/np.array([float(tdist.ppf(0.975, max(v, 1))) for v in dof])
            good = np.isfinite(sigma) & (sigma > 0)
            if not good.any():
                print("[EPRelax] no point has a usable interval, so the fit is unweighted")
                sigma = None
            elif not good.all():
                print(f"[EPRelax] {int((~good).sum())} of {x.size} points have no usable interval; "
                      f"they are given the median of the others")
                sigma = np.where(good, sigma, float(np.median(sigma[good])))
        try:
            popt, pcov = curve_fit(lambda xx, *vals: self.curve(xx, {**rates['pars'], **dict(zip(free, vals))}),
                                   x, y, p0=[rates['pars'][n] for n in free], sigma=sigma,
                                   absolute_sigma=sigma is not None, maxfev=40000)
        except (RuntimeError, ValueError, TypeError) as fault:
            self.result = None
            print(f"[EPRelax] the fit of '{rates['equation']}' did not converge: "
                  f"{type(fault).__name__}: {fault}")
            wx.MessageBox(f"The fit did not converge:\n\n{type(fault).__name__}: {fault}\n\n"
                          f"Starting values an order of magnitude out are the usual reason.",
                          "Fit", wx.OK | wx.ICON_WARNING)
            self.update_graph()
            return
        rates['pars'].update({n: float(v) for n, v in zip(free, popt)})
        nu = max(x.size - len(free), 1)
        err = np.sqrt(np.abs(np.diag(pcov)))*float(tdist.ppf(0.975, nu))
        self.result = {'pars': dict(rates['pars']), 'err': dict(zip(free, err.astype(float))),
                       'keys': keys, 'n': int(x.size), 'nu': nu, 'weighted': sigma is not None}
        self.touched = False
        self.push_params()
        self.update_graph()

    def on_clear(self, event):
        """Drop the fit. The equation and its values stay, and so does the curve while it is ticked."""
        self.result, self.touched = None, False
        self.update_graph()

    def on_save(self, event):
        """The points, their intervals and the fit as one CSV table."""
        series, unit, _warn = self.points()
        if not series:
            wx.MessageBox("There are no points to save.", "Save", wx.OK | wx.ICON_INFORMATION)
            return
        rates = self.mf.rates
        with wx.FileDialog(self, "Save the rates", defaultDir=self.mf.file_browser.path,
                           defaultFile="EPRelax_rates.csv", wildcard="CSV files (*.csv)|*.csv",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            path = dlg.GetPath()
        quantity = "rate (1/%s)" % unit if rates['ymode'] == 0 else f"T ({unit})"
        try:
            with open(path, 'w', encoding='utf-8-sig') as f:
                f.write(f"# EPRelax {VERSION} rates against '{rates['name']}'\n")
                f.write(f"# series: {('by T', 'by amplitude', 'all in one')[rates['series']]}, "
                        f"x drawn as {'1/parameter' if rates['xmode'] else 'parameter'}\n")
                if rates['pars']:
                    # the curve as it stands: the fit's own values, or whatever was typed over them
                    fresh = self.result is not None and not self.touched
                    err = self.result['err'] if fresh else {}
                    _keys, px, py, _dy, _dof = self.chosen(series, rates['fitmode'])
                    f.write(f"# {'fit' if fresh else 'equation'}: y = {rates['equation']}  "
                            f"(x is the parameter itself, y is "
                            f"{'the rate 1/T' if rates['fitmode'] == 0 else 'the time constant T'})\n")
                    f.write("# " + ", ".join(f"{n} = {v:.6g}"
                                             + (f" +-{err[n]:.3g}" if n in err else
                                                " (held)" if rates['hold'].get(n) else "")
                                             for n, v in rates['pars'].items()) + "\n")
                    if px.size:
                        resid = py - self.curve(px, rates['pars'])
                        spread = float(np.sum((py - py.mean())**2))
                        r2 = 1 - float(np.nansum(resid**2))/spread if spread > 0 else float('nan')
                        f.write(f"# over {px.size} points, "
                                f"{'weighted' if fresh and self.result['weighted'] else 'unweighted'}"
                                f", RMSD = {float(np.sqrt(np.nanmean(resid**2))):.6g}, R2 = {r2:.6f}\n")
                f.write(f"series,dataset,{rates['name'].replace(',', ' ')},{quantity},+-95%,T held\n")
                for key in sorted(series):
                    row = series[key]
                    y, dy = self.quantity(row, rates['ymode'])     # the column the header names
                    for j in range(row['x'].size):
                        f.write(f"{key + 1},\"{row['name'][j]}\",{row['x'][j]:.8g},{y[j]:.8g},"
                                f"{dy[j]:.6g},{'yes' if row['held'][j] else 'no'}\n")
        except OSError as fault:
            wx.MessageBox(f"Could not write '{path}':\n{fault}", "Save", wx.OK | wx.ICON_ERROR)
            return
        self.mf.statusbar.SetStatusText(f"Rates written to {path}", 0)

    def update_graph(self):
        """Read the fits as they stand and redraw."""
        series, unit, warn = self.points()
        rates = self.mf.rates
        # the Fit to list follows the series that are actually there
        labels = ["All points"] + [f"Component {k + 1}" for k in sorted(series)]
        if labels != [self.choice_target.GetString(i) for i in range(self.choice_target.GetCount())]:
            self.choice_target.Set(labels)
            self.choice_target.SetSelection(min(rates['target'], len(labels) - 1))
            rates['target'] = self.choice_target.GetSelection()
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        inverse = rates['xmode'] == 1
        drawn, held, dropped = 0, 0, 0
        for key in sorted(series):
            row = series[key]
            keep = row['x'] != 0 if inverse else np.ones(row['x'].shape, dtype=bool)
            if rates['logx']:
                keep &= (1.0/row['x'] if inverse else row['x']) > 0
            dropped += int((~keep).sum())
            if not keep.any():
                continue
            x = (1.0/row['x'][keep]) if inverse else row['x'][keep]
            y, dy = self.quantity(row, rates['ymode'])
            label = "components" if rates['series'] == 2 else f"component {key + 1}"
            if rates['bars']:
                ax.errorbar(x, y[keep], yerr=dy[keep], fmt='o', ms=5, elinewidth=1,
                            capsize=3, color=f"C{key % 10}", label=label)
            else:
                ax.plot(x, y[keep], 'o', ms=5, ls='none', color=f"C{key % 10}", label=label)
            drawn += int(keep.sum())
            held += sum(1 for j, h in enumerate(row['held']) if h and keep[j])
        # the curve is drawn from the values in the boxes, fitted or typed, so turning one moves it.
        # It is worked out in the quantity it was fitted in, and turned over where the axis shows
        # the other one, so it keeps lying on the points either way
        flip = rates['fitmode'] != rates['ymode']
        show = lambda v: np.divide(1.0, v, out=np.full_like(v, np.nan), where=v != 0) if flip else v
        keys, px, py, pdy, pdof = self.chosen(series, rates['fitmode'])
        fresh = self.result is not None and not self.touched
        resid = np.array([])
        if rates['draw'] and drawn and keys and self.equation()[1]:
            every = px[np.isfinite(px)]
            if every.size:
                lo, hi = float(every.min()), float(every.max())
                fine = (np.logspace(np.log10(lo), np.log10(hi), 400) if lo > 0 and hi/lo > 20
                        else np.linspace(lo, hi, 400))
                fine = fine[fine != 0] if inverse else fine
                # the curve is worked out in the parameter, then drawn against whatever is shown
                ax.plot(1.0/fine if inverse else fine, show(self.curve(fine, rates['pars'])),
                        lw=1.4, color=theme.theme_colors()["ax_fg"],
                        label=("fit" if fresh else "equation")
                              + (" in 1/T" if flip and rates['fitmode'] == 0 else
                                 " in T" if flip else ""))
                resid = py - self.curve(px, rates['pars'])
                # each term of the sum on its own, in colours the series have not taken. A term is
                # a term of the fitted quantity, so it is turned over with the whole, not summed
                # on the axis -- on the flipped axis the terms no longer add up to the curve
                parts = self.terms() if rates['terms'] else []
                for k, part in enumerate(parts if len(parts) > 1 else []):
                    ax.plot(1.0/fine if inverse else fine,
                            show(self.curve(fine, rates['pars'], part)),
                            lw=1.0, ls='--', alpha=0.9, color=f"C{(len(series) + k) % 10}",
                            label=part if len(part) <= 24 else part[:21] + "...")
        name = rates['name']
        ax.set_xlabel(f"1 / ({name})" if inverse else name)
        ax.set_ylabel(f"Rate, 1/{unit}" if rates['ymode'] == 0 else f"Time constant, {unit}")
        ax.set_xscale('log' if rates['logx'] else 'linear')
        ax.set_yscale('log' if rates['logy'] else 'linear')
        ax.grid(True, alpha=0.3)
        if drawn:
            ax.legend(fontsize='small', frameon=False)
        else:
            ax.text(0.5, 0.5, "No points yet: fit some datasets, and give each one a "
                              f"{name.split(',')[0].strip().lower()} in the Datasets list",
                    ha='center', va='center', transform=ax.transAxes, fontsize='small', color='0.5')
        theme.style_figure(self.figure, self.figure.axes)
        self.canvas.draw_idle()

        line = [f"{drawn} point{'s' if drawn != 1 else ''} from {len(series)} "
                f"series, {held} with a held T"]
        if dropped:
            line.append(f"{dropped} off the axis")
        line += warn
        if self.result is not None:
            r = self.result
            line.append(f"fit on {r['n']} points, {'weighted' if r['weighted'] else 'unweighted'}"
                        + (f", in {'1/T' if rates['fitmode'] == 0 else 'T'}" if flip else "")
                        + ("" if fresh else ", typed over since"))
        self.lbl_stat.SetLabel("   ".join(line))
        if not rates['pars']:
            self.lbl_fit.SetLabel("")
        else:
            # the interval belongs to the fit, so it is only shown while the values still are it
            err = self.result['err'] if fresh else {}
            text = "\n".join(f"{n} = {v:.6g}"
                             + (f" ± {err[n]:.3g}" if n in err else
                                "  (held)" if rates['hold'].get(n) else "")
                             for n, v in rates['pars'].items())
            if resid.size:
                spread = float(np.sum((py - py.mean())**2))
                r2 = 1 - float(np.nansum(resid**2))/spread if spread > 0 else float('nan')
                text += (f"\nRMSD = {float(np.sqrt(np.nanmean(resid**2))):.4g}   R² = {r2:.5f}"
                         f"   on {resid.size} point{'s' if resid.size != 1 else ''} of "
                         f"{'1/T' if rates['fitmode'] == 0 else 'T'}")
            if fresh:
                text += f"\n95 % intervals, {self.result['nu']} degrees of freedom"
            elif self.result is not None:
                text += "\nthe values have been typed over, so the fit's intervals are gone"
            self.lbl_fit.SetLabel(text)
        self.Layout()


# --------------------------------------------------------------------------- #
# Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    """The window, the datasets, and the analysis itself."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Data = []
        self.current = -1
        # shared by every dataset: comparing distributions means nothing unless they were made alike
        self.npts = 80            # grid points
        self.penalty = 0          # 0 curvature (CONTIN), 1 amplitude (Tikhonov)
        self.amode = 3            # 0 manual, 1 L-curve, 2 GCV, 3 F-test (CONTIN's own)
        self.prob = 0.5           # the F-test's probability to reject
        self.autoalpha = True     # bracket the alpha sweep from each trace rather than fix it
        self.alpha_lo = -5.0      # log10 of the swept range when it is not bracketed
        self.alpha_hi = 5.0
        self.rate_window = None   # RateParameterFrame while it is open; its settings live here
        self.rates = {'name': "Parameter", 'series': 0, 'ymode': 0, 'xmode': 0, 'logx': False,
                      'logy': False, 'bars': True, 'weight': True, 'draw': True, 'terms': False,
                      'fitmode': 0, 'target': 0,
                      'equation': "a*exp(-b/x) + c", 'pars': {}, 'hold': {}}
        self.SetTitle(f"EPRelax {VERSION}")
        self.SetSize((1280, 820))
        self.theme_name = "light"

        # ---- Top icon toolbar
        self.top_toolbar = wx.Panel(self)
        toolbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.tb_change_folder = theme.IconButton(self.top_toolbar, theme.draw_folder_icon, tooltip="Change folder...")
        self.tb_up = theme.IconButton(self.top_toolbar, theme.draw_up_icon, tooltip="Up")
        self.tb_load_all = theme.IconButton(self.top_toolbar, theme.draw_loadall_icon,
                                            tooltip="Load all\nLoad every file the browser lists")
        self.tb_save = theme.IconButton(
            self.top_toolbar, lambda gc, w, h, color: theme.draw_text_icon("CSV")(gc, w, h, color),
            tooltip="Export the current dataset...\nTrace, fit, residual and distribution as one CSV table")
        self.tb_session = theme.IconButton(self.top_toolbar, theme.draw_save_icon,
                                           tooltip="Save session...\nEvery dataset, every setting and the "
                                                   "view in one XML file,\nwhich opens again from the browser "
                                                   "without the Bruker files")
        self.tb_rates = theme.IconButton(self.top_toolbar, theme.draw_saturation_icon,
                                         tooltip="Rates against a parameter...\nThe discrete fits' relaxation "
                                                 "rates against a number you give each dataset,\nwith their 95 % "
                                                 "intervals and a fit of your own equation")
        self.tb_grid = theme.IconToggleButton(self.top_toolbar, theme.draw_grid_icon, tooltip="Axes grid", value=True)
        self.btn_theme_icon = theme.IconButton(self.top_toolbar, theme.draw_theme_icon,
                                               tooltip="Toggle dark mode (Ctrl+D)")
        # plot-row toggles, each icon the row's own line shape
        self.tb_rows = {
            'data': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: 2*np.exp(-2.5*(t + 1)) - 1), tooltip="Trace and fits", value=True),
            'resid': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: 0.45*np.sin(9*t)), tooltip="Residual", value=True),
            'dist': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: 2*np.exp(-14*t*t) - 1), tooltip="Distribution of relaxation times", value=True),
            'lcurve': theme.IconToggleButton(
                self.top_toolbar, lambda gc, w, h, color, on: theme.draw_text_icon("L")(gc, w, h, color),
                tooltip="L-curve of the current dataset\nResidual norm against ||L a|| over the alpha sweep")}
        for w in (self.tb_change_folder, self.tb_up, self.tb_load_all, None,
                  self.tb_save, self.tb_session, None, self.tb_rates, None,
                  *self.tb_rows.values(), 'stretch', None, self.tb_grid, self.btn_theme_icon):
            if w is None:
                toolbar_sizer.Add(wx.StaticLine(self.top_toolbar, style=wx.LI_VERTICAL), 0,
                                  wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
            elif isinstance(w, str):
                toolbar_sizer.AddStretchSpacer(1)
            else:
                toolbar_sizer.Add(w, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        self.top_toolbar.SetSizer(toolbar_sizer)
        self.tb_save.Bind(wx.EVT_BUTTON, self.on_save)
        self.tb_session.Bind(wx.EVT_BUTTON, lambda e: self.dumpXML())
        self.tb_rates.Bind(wx.EVT_BUTTON, self.on_rate_window)
        self.tb_grid.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_rows)
        for btn in self.tb_rows.values():
            btn.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_rows)
        self.btn_theme_icon.Bind(wx.EVT_BUTTON, self.on_toggle_theme)
        theme_id = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_toggle_theme, id=theme_id)
        self.SetAcceleratorTable(wx.AcceleratorTable([(wx.ACCEL_CMD, ord('D'), theme_id)]))

        self.statusbar = theme.StatusStrip(self, widths=[-1, 320])
        self.statusbar.SetStatusText("Ready", 0)

        # ---- file browser | (plot | cards), as in pyCWEPR
        self.splitter_main = wx.SplitterWindow(self, style=wx.SP_3D)
        self.splitter_main.SetMinimumPaneSize(20)
        self.file_browser = FileBrowserPanel(self.splitter_main, self)
        self.right_panel = wx.Panel(self.splitter_main)
        self.splitter_right = wx.SplitterWindow(self.right_panel, style=wx.SP_3D)
        self.plot_panel = PlotPanel(self.splitter_right, self)
        self.control_panel = ControlPanel(self.splitter_right, self)
        self.splitter_right.SetMinimumPaneSize(max(150, self.control_panel.floor))
        self.splitter_right.SplitVertically(self.plot_panel, self.control_panel)
        self.splitter_right.SetSashGravity(1)
        self.splitter_right.SetSashPosition(-330)
        right_sizer = wx.BoxSizer(wx.HORIZONTAL)
        right_sizer.Add(self.splitter_right, 1, wx.EXPAND)
        self.right_panel.SetSizer(right_sizer)
        self.left_sash = 220
        self.splitter_main.SplitVertically(self.file_browser, self.right_panel, sashPosition=self.left_sash)
        self.splitter_main.SetSashGravity(0)

        # the collapse toggle lives outside the splitter, so it stays visible once the pane is gone
        self.toggle_btn = theme.IconButton(self, theme.draw_chevron_icon("left"), tooltip="Hide file browser",
                                           size=wx.Size(18, -1), fill_key="accent")
        self.toggle_btn.Bind(wx.EVT_BUTTON, self.on_toggle_left)
        content_sizer = wx.BoxSizer(wx.HORIZONTAL)
        content_sizer.Add(self.toggle_btn, 0, wx.EXPAND)
        content_sizer.Add(self.splitter_main, 1, wx.EXPAND)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.top_toolbar, 0, wx.EXPAND)
        main_sizer.Add(content_sizer, 1, wx.EXPAND)
        main_sizer.Add(self.statusbar, 0, wx.EXPAND)
        self.SetSizer(main_sizer)

        self.resize_freeze = theme.ResizeFreezeGuard(self)
        self.Bind(wx.EVT_SIZE, self.resize_freeze.on_event)
        self.splitter_main.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self.resize_freeze.on_event)
        self.splitter_right.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self.resize_freeze.on_event)

        self.Centre()
        self.apply_theme()
        self.control_panel.update_controls()     # the shared settings, before anything is loaded
        self.plot_panel.update_graph()

    # -----------------------------------------------------------------------
    def on_toggle_theme(self, event):
        self.theme_name = "dark" if self.theme_name == "light" else "light"
        self.apply_theme()

    def on_toggle_left(self, event):
        """Hide or bring back the file browser, keeping where its sash was."""
        if self.splitter_main.IsSplit():
            self.left_sash = self.splitter_main.GetSashPosition()
            self.splitter_main.Unsplit(self.file_browser)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("right"))
            self.toggle_btn.SetToolTip("Show file browser")
        else:
            self.splitter_main.SplitVertically(self.file_browser, self.right_panel, self.left_sash)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("left"))
            self.toggle_btn.SetToolTip("Hide file browser")

    def on_rate_window(self, event):
        """Open the rates window, or raise it if it is already there."""
        if self.rate_window is None:
            self.rate_window = RateParameterFrame(self, self)
            self.rate_window.Show()
        else:
            self.rate_window.Raise()
        self.rate_window.update_graph()

    def on_toggle_rows(self, event):
        for name, btn in self.tb_rows.items():
            self.plot_panel.show[name] = btn.GetValue()
        self.plot_panel.showGrid = self.tb_grid.GetValue()
        self.plot_panel.update_graph()

    def apply_theme(self):
        theme.set_theme(self.theme_name)
        theme.apply_theme_to_window(self, self.theme_name)
        for ii, dd in enumerate(self.Data):     # the theme has just painted over the row colours
            self.control_panel.filelist.SetItemTextColour(ii, wx.Colour(to_hex(self.trace_color(dd))))
        theme.style_figure(self.plot_panel.figure, self.plot_panel.figure.axes)
        self.plot_panel.canvas.draw_idle()
        if self.rate_window is not None:
            self.rate_window.apply_theme()
        self.Refresh()

    def set_current(self, ii):
        if ii == self.current:        # clicking the row that is already current changes nothing
            return
        self.current = ii
        lst = self.control_panel.filelist
        if 0 <= ii < lst.GetItemCount() and not lst.IsSelected(ii):
            lst.Select(ii)                 # clicking a trace moves the list with it, and not back
        self.control_panel.update_controls()
        self.plot_panel.update_graph()

    def update_everything(self):
        """After datasets were added or removed: the list, the cards, the plot and the browser markers."""
        self.control_panel.update_list()
        self.control_panel.update_controls()
        self.plot_panel.update_graph()
        self.file_browser.refresh_file_list()

    def trace_color(self, dd):
        return dd.get('color') or f"C{dd['autocolor'] % 10}"

    def choose_color(self, ii, window):
        """Right-click on a trace, or on its row in the Datasets list: give that dataset a colour
        of its own, or hand it back to the automatic cycle."""
        if not 0 <= ii < len(self.Data):
            return
        menu = wx.Menu()
        pick, auto = menu.Append(wx.ID_ANY, "Colour..."), menu.Append(wx.ID_ANY, "Automatic")
        choice = window.GetPopupMenuSelectionFromUser(menu)
        if choice == auto.GetId():
            self.Data[ii]['color'] = ''
        elif choice == pick.GetId():
            data = wx.ColourData()
            data.SetColour(wx.Colour(to_hex(self.trace_color(self.Data[ii]))))
            with wx.ColourDialog(window, data) as dlg:
                if dlg.ShowModal() != wx.ID_OK:
                    return
                self.Data[ii]['color'] = dlg.GetColourData().GetColour().GetAsString(wx.C2S_HTML_SYNTAX)
        else:
            return
        # that row's colour, not the whole list: a rebuild would lose where the list is scrolled to
        self.control_panel.filelist.SetItemTextColour(ii, wx.Colour(to_hex(self.trace_color(self.Data[ii]))))
        self.plot_panel.update_graph()

    def harmonize_units(self):
        """Put every dataset on one time unit, chosen so the longest time of them all stays under
        1000. Overlaid traces share an axis, so they have to share a unit: left to pick their own,
        a 2 ms recovery would be drawn as 2 next to a 20 microsecond decay drawn as 20, and the
        longest trace would look the shortest. A dataset whose file names a unit off the ladder
        keeps its own numbers and takes no part. Returns the unit if anything was rescaled."""
        ladder = dict(TIME_UNITS)
        known = [dd for dd in self.Data if dd['tunit'] in ladder and dd['t'].size]
        if not known:
            return None
        seconds = max(float(np.max(np.abs(dd['t'])))*ladder[dd['tunit']] for dd in known)
        name, factor = TIME_UNITS[-1]
        for name, factor in TIME_UNITS:
            if seconds < 1000*factor:
                break
        moved = False
        for dd in known:
            ratio = ladder[dd['tunit']]/factor
            if ratio == 1.0:
                continue
            dd['t'] = dd['t']*ratio
            dd['win'] = [v*ratio for v in dd['win']]
            dd['Tmin'], dd['Tmax'], dd['swin'] = dd['Tmin']*ratio, dd['Tmax']*ratio, dd['swin']*ratio
            if dd['disc'] is not None:
                dd['disc']['T'], dd['disc']['dT'] = dd['disc']['T']*ratio, dd['disc']['dT']*ratio
            dd['tunit'], dd['sig'], dd['res'] = name, None, None
            moved = True
        return name if moved else None

    # -----------------------------------------------------------------------
    # Data handling
    # -----------------------------------------------------------------------
    def loadData(self, path, report=True, refresh=True):
        """Read one Bruker dataset and add it. Returns the exception on failure, None on success."""
        if path.lower().endswith('.xml'):
            return self.load_session(path, report=report)
        try:
            ax, data, dsc = Mybr.brukerread(path, return_ax=True, return_dsc=True)
            t = np.asarray(ax['x'], dtype=float).ravel()
            data = np.squeeze(np.asarray(data))
            if data.ndim > 2 or data.shape[0] != t.size:
                raise ValueError(f"not a relaxation trace (data shape {np.shape(data)}, {t.size} time points)")
            if t.size < 8:
                raise ValueError(f"only {t.size} time points")
            label = str(ax.get('xlabel', ''))
            unit = label.split(',')[-1].strip() if ',' in label else label.split()[-1] if label else ''
            # the file's own unit, spelled the way the ladder spells it; harmonize_units then puts
            # this dataset and every other one on whichever unit suits the longest trace loaded
            unit = {'us': 'µs', 'usec': 'µs', 'µsec': 'µs', 'nsec': 'ns', 'msec': 'ms',
                    'sec': 's', 'psec': 'ps'}.get(unit, unit)
            span = float(t[-1] - t[0])
            title = str(ax.get('title') or os.path.basename(path)).strip() or os.path.basename(path)
            # a trace whose magnitude grows towards the end is a recovery, not a decay
            probe = np.abs(np.real(data if data.ndim == 1 else data[:, 0]).astype(float))
            edge = max(probe.size//10, 1)
            recovery = bool(probe[-edge:].mean() > probe[:edge].mean())
            used = {d['autocolor'] for d in self.Data}
            dd = {'t': t, 'data': data, 'dsc': dsc if isinstance(dsc, dict) else {},
                  'fullpath': path, 'title': title, 'tunit': unit or 'a.u.',
                  'show': True, 'autocolor': next(k for k in range(len(used) + 1) if k not in used),
                  'color': '',                   # chosen by hand; empty means the automatic cycle
                  'param': float('nan'),         # the number this dataset is plotted against

                  'comp': 0,                     # 0 real as stored, 1 imaginary, 2 magnitude
                  'slice': 0, 'sum': False,      # which slice of a 2D set, or all of them added
                  'win': [float(t[0]), float(t[-1])],
                  'recovery': recovery, 'offset': True,
                  'smooth': False, 'swin': 0.0,  # boxcar width, in the dataset's own time unit
                  'autogrid': True, 'Tmin': span/100, 'Tmax': 3*span,
                  'logalpha': -2.0,              # the grid density, penalty and alpha mode are shared
                  'resample': False, 'perdec': 16.0, 'rweight': 'count',
                  'ndisc': 2, 'fixT': False, 'disc': None, 'res': None, 'sig': None, 'sweep': None}
            while dd['title'] in [d['title'] for d in self.Data]:
                dd['title'] += "'"
            self.Data.append(dd)
            self.current = len(self.Data) - 1
            rescaled = self.harmonize_units()
            if refresh:
                self.update_everything()
            self.statusbar.SetStatusText(
                f"Loaded '{dd['title']}'" + (", read as a recovery" if recovery else "")
                + (f"; shown in {rescaled}" if rescaled else ""), 0)
            return None
        except Exception as fault:
            traceback.print_exc()
            if report:
                wx.MessageBox(f"Loading '{os.path.basename(path)}' failed:\n\n{type(fault).__name__}: {fault}",
                              "Cannot load", wx.OK | wx.ICON_ERROR)
            return fault

    def on_save(self, event):
        """Write the current dataset's trace, fit, residual and distribution as one CSV table."""
        if self.current < 0 or self.Data[self.current]['res'] is None:
            wx.MessageBox("Nothing has been decomposed yet.", "Export", wx.OK | wx.ICON_INFORMATION)
            return
        dd = self.Data[self.current]
        r = dd['res']
        with wx.FileDialog(self, "Export the current dataset", defaultFile=f"{dd['title']}_ILT.csv",
                           wildcard="CSV files (*.csv)|*.csv",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            path = dlg.GetPath()
        rows = max(r['t'].size, r['T'].size)
        pad = lambda v: list(v) + ['']*(rows - len(v))
        cols = [pad(r['t']), pad(r['y']), pad(r['fit']), pad(r['y'] - r['fit']), pad(r['T']), pad(r['a'])]
        names = [f"t ({dd['tunit']})", "data", "fit", "residual", f"T ({dd['tunit']})", "amplitude"]
        if dd['disc'] is not None:
            cols.append(pad(self.discrete_curve(dd, r['t'])))
            names.append("discrete fit")
        try:
            # the byte-order mark so a spreadsheet reads the µ of a microsecond axis as written
            with open(path, 'w', encoding='utf-8-sig') as f:
                f.write(f"# EPRelax {VERSION} export of: {dd['title']}\n")
                f.write(f"# source: {dd['fullpath']}\n")
                f.write(f"# channel: {('real', 'imaginary', 'magnitude')[dd['comp']]}, "
                        f"window: {min(dd['win']):.6g} to {max(dd['win']):.6g} {dd['tunit']}, "
                        f"recovery: {dd['recovery']}, offset fitted: {dd['offset']}\n")
                boxcar = f"{dd['swin']:.6g} {dd['tunit']}" if dd['smooth'] else "off"
                f.write(f"# boxcar: {boxcar}\n")
                f.write(f"# penalty: {('curvature', 'amplitude')[self.penalty]}, "
                        f"alpha: {r['alpha']:.6g} ({('manual', 'L-curve', 'GCV', 'F-test')[self.amode]}), "
                        f"effective d.o.f.: {r['df']:.3f}\n")
                f.write(f"# offset: {r['c']:.6g}, total amplitude: {r['a'].sum():.6g}, "
                        f"RMS residual: {r['rms']:.6g}\n")
                d = dd['disc']
                if d is not None:
                    f.write(f"# discrete fit: offset={d['c']:.6g}+-{d['dc']:.3g}, RMS={d['rms']:.6g}, "
                            f"95 % intervals from the trust region, {d['dof']} d.o.f.\n")
                    f.write(f"# component, T ({dd['tunit']}), +-95%, a, +-95%, iLT low, iLT high, "
                            f"note\n")
                    edge = lambda v: f"{v:.6g}" if np.isfinite(v) else ""
                    for i, Td in enumerate(d['T']):
                        held = ("~" if d['Thold'] else "") + f"{d['dT'][i]:.3g}"
                        f.write(f"# {i + 1}, {Td:.6g}, {held}, {d['a'][i]:.6g}, {d['da'][i]:.3g}, "
                                f"{edge(Td/d['wide'][i])}, {edge(Td*d['wide'][i])}, "
                                f"{d['note'][i] or '-'}\n")
                f.write(",".join(names) + "\n")
                for k in range(rows):
                    f.write(",".join(f"{c[k]:.8g}" if c[k] != '' else '' for c in cols) + "\n")
        except OSError as fault:
            wx.MessageBox(f"Could not write '{path}':\n{fault}", "Export", wx.OK | wx.ICON_ERROR)
            return
        self.statusbar.SetStatusText(f"Exported to {path}", 0)

    def save_discrete(self):
        """Every component of every fitted dataset in one CSV table. Relaxation times go in seconds
        throughout, whatever unit each trace reads in, so the column sorts and plots across samples;
        the unit the dataset was shown in is named in a column of its own."""
        fitted = [dd for dd in self.Data if dd['disc'] is not None]
        if not fitted:
            wx.MessageBox("No dataset has been fitted yet.", "Save discrete fits",
                          wx.OK | wx.ICON_INFORMATION)
            return
        with wx.FileDialog(self, "Save the discrete fits", defaultFile="EPRelax_discrete.csv",
                           defaultDir=self.file_browser.path, wildcard="CSV files (*.csv)|*.csv",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            path = dlg.GetPath()
        # no commas inside a column name, or the header splits into more fields than there are columns
        names = ("dataset", "components", "component", "T (s)", "+-95% (s)", "a", "+-95%",
                 "iLT low (s)", "iLT high (s)", "note",
                 "offset", "offset +-95%", "RMS residual", "fit d.o.f.", "alpha", "shown in", "source")
        try:
            with open(path, 'w', encoding='utf-8-sig') as f:
                f.write(f"# EPRelax {VERSION} discrete fits of {len(fitted)} datasets\n")
                f.write(f"# penalty: {('curvature', 'amplitude')[self.penalty]}, alpha from: "
                        f"{('manual', 'L-curve', 'GCV', 'F-test')[self.amode]}, grid points: {self.npts}\n")
                f.write("# the 95 % intervals on T and a are the non-linear fit's own trust region; "
                        "'iLT low' and 'iLT high' are the ends of the distribution peak the component "
                        "sits in, which is its spread of relaxation times rather than an error on T\n")
                f.write(",".join(names) + "\n")
                for dd in fitted:
                    d = dd['disc']
                    sec = dict(TIME_UNITS).get(dd['tunit'], 1.0)   # unknown unit: the file's own numbers
                    quoted = lambda text: '"' + str(text).replace('"', "'") + '"'
                    edge = lambda v: f"{v:.8g}" if np.isfinite(v) else ""   # no peak, no peak width
                    tail = [f"{d['c']:.8g}", f"{d['dc']:.6g}", f"{d['rms']:.6g}", str(d['dof']),
                            f"{dd['res']['alpha']:.6g}" if dd['res'] else '',
                            dd['tunit'], quoted(dd['fullpath'])]
                    for i, Td in enumerate(d['T']):
                        # the combined table stays numeric; the note column says T was held
                        held = f"{d['dT'][i]*sec:.6g}"
                        f.write(",".join([quoted(dd['title']), str(len(d['T'])), str(i + 1),
                                          f"{Td*sec:.8g}", held,
                                          f"{d['a'][i]:.8g}", f"{d['da'][i]:.6g}",
                                          edge(Td*sec/d['wide'][i]), edge(Td*sec*d['wide'][i]),
                                          quoted(d['note'][i])] + tail) + "\n")
        except OSError as fault:
            wx.MessageBox(f"Could not write '{path}':\n{fault}", "Save discrete fits",
                          wx.OK | wx.ICON_ERROR)
            return
        self.statusbar.SetStatusText(f"{sum(len(dd['disc']['T']) for dd in fitted)} components from "
                                     f"{len(fitted)} datasets written to {path}", 0)
        self.file_browser.refresh_file_list()

    def dumpXML(self):
        """The whole session as one self-contained XML: every trace as a base64 numpy array, every
        setting, and what the plot is showing, so it reopens without the original Bruker files."""
        with wx.FileDialog(self, "Save Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultFile="EPRelax_session.xml", defaultDir=self.file_browser.path,
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            pathname = dlg.GetPath()
        pp = self.plot_panel
        session = {'program': 'EPRelax', 'version': VERSION,
                   # 'res', 'sig' and 'sweep' are caches, and 'res' holds the whole kernel
                   'Data': [{k: v for k, v in dd.items() if k not in ('res', 'sig', 'sweep')} for dd in self.Data],
                   'shared': {'npts': self.npts, 'penalty': self.penalty, 'amode': self.amode,
                              'prob': self.prob, 'autoalpha': self.autoalpha,
                              'alpha_lo': self.alpha_lo, 'alpha_hi': self.alpha_hi},
                   # the rates window keeps its settings here, so they are saved whether it is open or not
                   'rates': dict(self.rates),
                   'view': {'show': dict(pp.show), 'showGrid': pp.showGrid, 'norm': pp.norm,
                            'xmode': pp.xmode, 'sep': pp.sep, 'stack': pp.stack,
                            'stackmode': pp.stackmode, 'logt': pp.logt, 'snorm': pp.snorm,
                            'showfit': pp.showfit, 'showdisc': pp.showdisc, 'legend': pp.legend,
                            'current': self.current}}

        def create_element(parent, key, value):
            element = ET.SubElement(parent, 'item')
            element.set('key', str(key))
            if isinstance(value, np.ndarray):
                array_element = ET.SubElement(element, 'numpy_array')
                array_element.set('dtype', str(value.dtype))
                array_element.set('shape', str(value.shape))
                ET.SubElement(array_element, 'data').text = base64.b64encode(value.tobytes()).decode('utf-8')
            elif isinstance(value, dict):
                dict_element = ET.SubElement(element, 'dict')
                for k, v in value.items():
                    create_element(dict_element, k, v)
            elif isinstance(value, (list, tuple)):
                list_element = ET.SubElement(element, 'list')
                list_element.set('type', type(value).__name__)
                for item in value:
                    create_element(list_element, 'item', item)
            else:
                primitive_element = ET.SubElement(element, 'primitive')
                primitive_element.set('type', type(value).__name__)
                primitive_element.text = str(value)
        try:
            root = ET.Element('dictionary')
            for key, value in session.items():
                create_element(root, key, value)
            pretty = minidom.parseString(ET.tostring(root, encoding='unicode')).toprettyxml(indent="  ")
            with open(pathname, 'w', encoding='utf-8') as f:
                f.write(pretty)
        except (OSError, TypeError, ValueError) as fault:
            traceback.print_exc()
            wx.MessageBox(f"Could not write '{pathname}':\n{type(fault).__name__}: {fault}",
                          "Save Session", wx.OK | wx.ICON_ERROR)
            return
        self.statusbar.SetStatusText(f"Session saved to {pathname}", 0)
        self.file_browser.refresh_file_list(select_name=os.path.basename(pathname))

    def load_mixed_dict_xml(self, filename):
        """Dictionary back from the session XML encoding, with original types."""
        def parse_element(element):
            key = element.get('key')
            numpy_array = element.find('numpy_array')
            dict_element = element.find('dict')
            list_element = element.find('list')
            primitive_element = element.find('primitive')
            if numpy_array is not None:
                dtype = numpy_array.get('dtype')
                shape = ast.literal_eval(numpy_array.get('shape'))
                try:
                    array = np.frombuffer(base64.b64decode(numpy_array.find('data').text),
                                          dtype=dtype).reshape(shape).copy()
                except Exception as fault:
                    print(f"[EPRelax] Could not reconstruct array '{key}' (dtype={dtype}, shape={shape}) "
                          f"from '{filename}': {type(fault).__name__}: {fault}")
                    array = None
                return key, array
            if dict_element is not None:
                result = {}
                for item in dict_element.findall('item'):
                    k, v = parse_element(item)
                    if k is not None:
                        result[k] = v
                return key, result
            if list_element is not None:
                items = [parse_element(item)[1] for item in list_element.findall('item')]
                return key, tuple(items) if list_element.get('type') == 'tuple' else items
            if primitive_element is not None:
                value_type = primitive_element.get('type')
                value_text = primitive_element.text
                if value_type in ('int', 'int32', 'int64'):
                    return key, int(value_text)
                if value_type in ('float', 'float32', 'float64'):
                    return key, float(value_text)
                if value_type == 'NoneType':
                    return key, None
                if value_type in ('bool', 'bool_'):
                    return key, value_text.lower() == 'true'
                return key, value_text if value_text is not None else ''
            return key, None

        result = {}
        for item in ET.parse(filename).getroot().findall('item'):
            key, value = parse_element(item)
            if key is not None:
                result[key] = value
        return result

    def load_session(self, path, report=True):
        """Replace everything with a saved session. Returns the exception on failure, None on success."""
        try:
            session = self.load_mixed_dict_xml(path)
            data = session.get('Data')
            if not isinstance(data, list) or not all(isinstance(dd, dict) and 'data' in dd for dd in data):
                raise ValueError("no datasets in it; is this an EPRelax session?")
        except Exception as fault:
            traceback.print_exc()
            if report:
                wx.MessageBox(f"Could not read the session '{os.path.basename(path)}':\n\n"
                              f"{type(fault).__name__}: {fault}", "Open Session", wx.OK | wx.ICON_ERROR)
            return fault
        if self.Data and wx.MessageBox(f"'{os.path.basename(path)}' holds {len(data)} datasets.\n\n"
                                       f"Opening it drops the {len(self.Data)} loaded here. Continue?",
                                       "Open Session", wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
            return None
        for dd in data:                      # the caches are rebuilt on the first redraw
            dd['res'], dd['sig'], dd['sweep'] = None, None, None
            dd.setdefault('disc', None)
            dd.setdefault('fixT', False)
            dd.setdefault('resample', False)
            dd.setdefault('perdec', 16.0)
            dd.setdefault('rweight', 'count')
            dd.setdefault('param', float('nan'))
            dd['color'] = dd.get('color') or ''
        self.Data = data
        shared = session.get('shared') or {}
        self.npts = int(shared.get('npts', self.npts))
        self.penalty = int(shared.get('penalty', self.penalty))
        self.amode = int(shared.get('amode', self.amode))
        self.prob = float(shared.get('prob', self.prob))
        self.autoalpha = bool(shared.get('autoalpha', self.autoalpha))
        self.alpha_lo = float(shared.get('alpha_lo', self.alpha_lo))
        self.alpha_hi = float(shared.get('alpha_hi', self.alpha_hi))
        saved = session.get('rates') or {}
        for key in ('name', 'equation'):
            self.rates[key] = str(saved.get(key, self.rates[key]))
        for key in ('series', 'ymode', 'xmode', 'fitmode', 'target'):
            self.rates[key] = int(saved.get(key, self.rates[key]))
        for key in ('logx', 'logy', 'bars', 'weight', 'draw', 'terms'):
            self.rates[key] = bool(saved.get(key, self.rates[key]))
        for key in ('pars', 'hold'):
            self.rates[key] = dict(saved[key]) if isinstance(saved.get(key), dict) else {}
        view = session.get('view') or {}
        pp = self.plot_panel
        pp.show.update({k: bool(v) for k, v in (view.get('show') or {}).items() if k in pp.show})
        pp.showGrid = bool(view.get('showGrid', pp.showGrid))
        pp.norm = str(view.get('norm', pp.norm))
        pp.xmode = str(view.get('xmode', pp.xmode))
        pp.sep = float(view.get('sep', pp.sep))
        pp.stack = bool(view.get('stack', pp.stack))
        pp.stackmode = str(view.get('stackmode', pp.stackmode))
        pp.logt = bool(view.get('logt', pp.logt))
        pp.snorm = str(view.get('snorm', pp.snorm))
        pp.showfit = bool(view.get('showfit', pp.showfit))
        pp.showdisc = bool(view.get('showdisc', pp.showdisc))
        pp.legend = bool(view.get('legend', pp.legend))
        self.current = min(int(view.get('current', 0)), len(self.Data) - 1)
        cp = self.control_panel
        cp.seg_xmode.SetSelection(1 if pp.xmode == 'rate' else 0)
        cp.seg_norm.SetSelection({'area': 1, 'fraction': 2}.get(pp.norm, 0))
        cp.spin_sep.SetValue(pp.sep)
        cp.chk_stack.SetValue(pp.stack)
        cp.seg_stackmode.SetSelection(1 if pp.stackmode == 'bands' else 0)
        cp.chk_logt.SetValue(pp.logt)
        cp.seg_snorm.SetSelection(1 if pp.snorm == 'span' else 0)
        cp.chk_showfit.SetValue(pp.showfit)
        cp.chk_showdisc.SetValue(pp.showdisc)
        cp.chk_legend.SetValue(pp.legend)
        if self.rate_window is not None:
            self.rate_window.pull()
        self.harmonize_units()        # sessions written before the units were shared may be mixed
        self.tb_grid.SetValue(pp.showGrid)
        for name, btn in self.tb_rows.items():
            btn.SetValue(pp.show[name])
        self.update_everything()
        self.statusbar.SetStatusText(f"Session '{os.path.basename(path)}' opened: "
                                     f"{len(self.Data)} datasets", 0)
        return None

    # -----------------------------------------------------------------------
    # Analysis
    # -----------------------------------------------------------------------
    def iLT_decompose(self, A, L, y, alpha, n, w=None):
        """One non-negative regularized solve of the augmented system, each point weighted by w
        where resampling gave the points unequal precision. Returns amplitudes, offset, fit,
        residual sum of squares, the penalty seminorm, and the degrees of freedom the penalty
        actually leaves (the trace of the hat matrix restricted to the amplitudes NNLS kept)."""
        Aw, yw = (A, y) if w is None else (A*w[:, None], y*w)
        x = nnls(np.vstack([Aw, alpha*L]), np.concatenate([yw, np.zeros(L.shape[0])]))[0]
        fit = A @ x
        res = (y - fit) if w is None else (y - fit)*w
        rss = float(np.sum(res*res))
        eta = float(np.linalg.norm(L @ x))
        active = x > 0
        df = float(active.sum())
        if active.any():
            As, Ls = Aw[:, active], L[:, active]
            M = As.T @ As
            try:
                df = float(np.trace(np.linalg.solve(M + alpha*alpha*(Ls.T @ Ls), M)))
            except np.linalg.LinAlgError as fault:
                print(f"[EPRelax] the effective d.o.f. at alpha={alpha:g} fell back to the active-set "
                      f"count: {fault}")
        c = float(x[n] - x[n + 1]) if x.size > n else 0.0
        return x[:n], c, fit, rss, eta, df

    def analyse(self, dd):
        """Prepare the trace, choose alpha, decompose, and leave everything the plots need in dd['res'].
        Skipped when nothing that matters has changed since the last pass."""
        # the kernel part first: what the sweep over alpha depends on. Which alpha is then wanted --
        # the knob, the mode, the F-test probability -- costs one solve, not another thirty-three
        ksig = (dd['comp'], dd['slice'], dd['sum'], tuple(dd['win']), dd['recovery'], dd['offset'],
                dd['smooth'], dd['swin'], dd['autogrid'], self.npts, self.penalty,
                dd['resample'], dd['perdec'], dd['rweight'],
                self.autoalpha, self.alpha_lo, self.alpha_hi,
                None if dd['autogrid'] else (dd['Tmin'], dd['Tmax']))
        sig = ksig + (self.amode, dd['logalpha'], self.prob)
        # drawing the L-curve is not a change to the trace, so it is no part of the signature: it
        # only asks for the sweep, which a dataset decomposed without it has not got yet
        if dd['sig'] == sig and dd['res'] is not None and not (self.plot_panel.show['lcurve']
                                                               and dd['res']['rho'] is None):
            return
        dd['sig'] = sig

        data = dd['data']
        col = data if data.ndim == 1 else (data.sum(axis=1) if dd['sum']
                                           else data[:, min(dd['slice'], data.shape[1] - 1)])
        y = {0: np.real, 1: np.imag, 2: np.abs}[dd['comp']](col).astype(float)
        t = dd['t']
        keep = (t >= min(dd['win'])) & (t <= max(dd['win']))
        t, y = t[keep], y[keep]
        if t.size < 6 or not np.all(np.isfinite(y)):
            print(f"[EPRelax] '{dd['title']}' has {t.size} usable points in the window; nothing to decompose")
            dd['res'] = None
            return
        if dd['smooth'] and dd['swin'] > 0:
            # a boxcar one modulation period wide cancels the modulation and only rescales an exponential
            width = int(round(dd['swin']/max(abs(float(np.mean(np.diff(t)))), 1e-30)))
            y = uniform_filter1d(y, min(max(width, 3), t.size), mode='nearest')

        # log-spaced bins across the window, each point the mean of what fell in its bin. A bin
        # holding nothing is dropped, so the early part is never upsampled, and a bin holding m
        # points has averaged its noise down by sqrt(m) -- which is the weight it then carries
        weight = None
        if dd['resample'] and dd['perdec'] > 0:
            # zero and any negative time have no logarithm and cannot be binned, so they are kept
            # as they are: a recovery measured from t = 0 keeps its first points, one per point
            early = t <= 0
            tp, yp = t[~early], y[~early]
            if tp.size < 2:
                print(f"[EPRelax] '{dd['title']}' has no positive times in the window, so there is "
                      f"nothing to resample")
            else:
                decades = float(np.log10(tp[-1]/tp[0]))
                edges = np.logspace(np.log10(tp[0]), np.log10(tp[-1]),
                                    max(int(round(dd['perdec']*decades)), 2) + 1)
                idx = np.clip(np.searchsorted(edges, tp, side='right') - 1, 0, edges.size - 2)
                m = np.bincount(idx, minlength=edges.size - 1)
                keep = m > 0
                tb = np.bincount(idx, weights=tp, minlength=edges.size - 1)[keep]/m[keep]
                yb = np.bincount(idx, weights=yp, minlength=edges.size - 1)[keep]/m[keep]
                if tb.size + early.sum() < 6:
                    print(f"[EPRelax] resampling '{dd['title']}' at {dd['perdec']:g} points per "
                          f"decade leaves only {tb.size} points; the raw trace is used instead")
                else:
                    counts = np.concatenate([np.ones(int(early.sum())), m[keep].astype(float)])
                    t, y = np.concatenate([t[early], tb]), np.concatenate([y[early], yb])
                    weight = np.sqrt(counts) if dd['rweight'] == 'count' else None

        scale = max(float(np.max(np.abs(y))), 1e-300)     # alpha means the same whatever the signal units
        ys = y/scale
        span = max(float(t[-1] - t[0]), 1e-30)
        if dd['autogrid']:
            # the shortest time the trace can argue about is its first point, or its sampling
            # step where that is coarser. A hundredth of the span would do for an evenly sampled
            # trace, but on a logarithmic time list the span is all last point and the whole fast
            # half of the decay would fall off the end of the grid
            step = float(np.min(np.diff(t))) if t.size > 1 else span
            dd['Tmin'], dd['Tmax'] = max(float(t[0]), step, span*1e-6), 3*span
        lo, hi = sorted((max(dd['Tmin'], 1e-30), max(dd['Tmax'], 1e-29)))
        T = np.logspace(np.log10(lo), np.log10(hi), self.npts)
        n = T.size

        K = np.exp(-t[:, None]/T[None, :])
        if dd['recovery']:
            K = -K
        A = np.hstack([K, np.ones((t.size, 1)), -np.ones((t.size, 1))]) if dd['offset'] else K
        L = np.zeros((n, A.shape[1]))
        j = np.arange(n)
        if self.penalty == 0:            # CONTIN: the curvature along the log-T grid, with the
            L[j, j] = -2.0                # distribution taken as zero outside it, so that the two end
            L[j[:-1], j[:-1] + 1] = 1.0   # points cost what the middle does instead of 2.4 times less
            L[j[1:], j[1:] - 1] = 1.0
        else:                             # Tikhonov: the amplitudes themselves
            L[j, j] = 1.0

        alpha = 10.0**dd['logalpha']
        rho = eta = None
        if self.amode or self.plot_panel.show['lcurve']:
            cached = dd.get('sweep')
            if cached is not None and cached[0] == ksig:
                rss_s, eta_s, df_s, alphas = cached[1:]
            else:
                lo, hi = sorted((10.0**self.alpha_lo, 10.0**self.alpha_hi))
                hi = max(hi, lo*10)
                if self.autoalpha:
                    # one probe per decade: climb until the penalty has flattened the distribution
                    # to a single effective component, and start where the solve stops changing.
                    # The corner moved seven decades between test traces, so a fixed ladder puts
                    # the whole informative part of the curve outside the window for some of them
                    probe = np.logspace(-12, 12, 25)
                    rss_p, df_p, top = [], [], probe[-1]
                    for al in probe:
                        out = self.iLT_decompose(A, L, ys, al, n, weight)
                        rss_p.append(out[3])
                        df_p.append(out[5])
                        if df_p[-1] <= 1.01:      # one effective component left: nothing beyond here
                            top = al
                            break
                    rss_p = np.array(rss_p)
                    plateau = np.flatnonzero(rss_p <= 1.001*rss_p[0])
                    lo = float(probe[max(int(plateau[-1]) - 1, 0)]) if plateau.size else float(probe[0])
                    hi = max(float(top), lo*10)
                alphas = np.logspace(np.log10(lo), np.log10(hi), SWEEP)
                sweep = np.array([self.iLT_decompose(A, L, ys, al, n, weight)[3:] for al in alphas])
                rss_s, eta_s, df_s = sweep[:, 0], np.maximum(sweep[:, 1], 1e-300), sweep[:, 2]
                dd['sweep'] = (ksig, rss_s, eta_s, df_s, alphas)
            rho, eta = np.sqrt(np.maximum(rss_s, 1e-300)), eta_s
            m = float(t.size)
            if self.amode == 1:
                # Hansen's corner: the maximum curvature of the log-log L-curve, over the stretch
                # where it actually moves -- below some alpha the solve stops changing at all
                p, q = np.log(rho), np.log(eta)
                move = np.concatenate([[True], (np.diff(p)**2 + np.diff(q)**2) > 1e-12])
                if move.sum() < 5:
                    print(f"[EPRelax] the L-curve of '{dd['title']}' barely moves over the sweep, so it "
                          f"has no corner to find; alpha stays at the knob")
                else:
                    p, q = p[move], q[move]
                    d1p, d1q = np.gradient(p), np.gradient(q)
                    d2p, d2q = np.gradient(d1p), np.gradient(d1q)
                    kappa = (d1p*d2q - d1q*d2p)/np.power(d1p**2 + d1q**2, 1.5)
                    alpha = float(alphas[move][int(np.nanargmax(kappa[1:-1])) + 1])
            elif self.amode == 2:
                alpha = float(alphas[int(np.argmin(m*rss_s/np.maximum(m - df_s, 1e-6)**2))])
            elif self.amode == 3:
                # named rather than left as the last case: the sweep above also runs in Manual when
                # the L-curve row is drawn, and a bare else would then overrule the knob
                # CONTIN's F-test: the largest alpha whose fit is still no worse than the unregularized
                # one. Regularizing never buys back less than one degree of freedom, so dnu floors at 1
                rss0, df0 = rss_s[0], df_s[0]
                dnu = np.maximum(df0 - df_s, 1.0)
                fstat = ((rss_s - rss0)/dnu)/max(rss0/max(m - df0, 1e-6), 1e-300)
                crit = fdist.ppf(self.prob, dnu, max(m - df0, 1.0))
                passed = np.flatnonzero((fstat <= crit) & (df_s > 0))
                alpha = float(alphas[passed[-1]]) if passed.size else float(alphas[0])

        a, c, fit, rss, eta_at, df = self.iLT_decompose(A, L, ys, alpha, n, weight)
        if not np.any(a > 0):
            print(f"[EPRelax] '{dd['title']}' came out with every amplitude at zero, so only the offset "
                  f"is left: the trace runs the other way than the kernel. Try the Recovery switch, "
                  f"or another channel.")
        dd['res'] = {'t': t, 'y': y, 'T': T, 'a': a*scale, 'c': c*scale, 'fit': fit*scale,
                     # the scatter of the points that are drawn, not the weighted criterion
                     'rms': float(np.sqrt(np.mean((ys - fit)**2)))*scale, 'df': df, 'alpha': alpha,
                     'rho': rho, 'eta': eta,
                     # where this solve itself sits on that curve, which is not a swept point
                     'point': (float(np.sqrt(max(rss, 1e-300))), eta_at),
                     # kept so the discrete fit can put a clean exponential through the same decomposition
                     'A': A, 'L': L, 'scale': scale, 'w': weight}

    def discrete_curve(self, dd, t):
        """The discrete fit drawn on whatever time axis is current. The parameters are what the fit
        found; the curve is not stored, because the window can be dragged after a fit and a curve
        sampled on the old one is both the wrong length and in the wrong place."""
        d = dd['disc']
        sign = -1.0 if dd['recovery'] else 1.0
        return d['c'] + sign*sum(a*np.exp(-t/T) for a, T in zip(d['a'], d['T']))

    def peakwidth(self, T, a, Td):
        """The distribution peak that Td falls in, as (sigma, centre, note). The centre is the
        amplitude-weighted mean of log10 T over the peak, read back as a time -- the geometric
        mean, which is what a logarithmic grid makes of an average, and the same quantity the
        width below is a spread about. The width itself is a Gaussian in log10 T.
        The peak region runs out to the minima either side. A Td that is on no peak at all -- in a
        valley, or on the flank of one -- has no width to report, and sigma comes back NaN with the
        note saying so. Otherwise an answer always comes back, and the note says when it should not
        be trusted: a peak clipped by the end of the grid, too few points under it, or a fit that
        would not converge."""
        x = np.log10(T)
        j = int(np.argmin(np.abs(x - np.log10(Td))))
        if a[j] <= 0:
            return float('nan'), float(Td), "no peak here"
        # a fitted time constant can land anywhere -- on the flank of a peak, or in the valley
        # between two. There is no peak there to have a width, so none is reported rather than
        # measuring whatever happens to surround the point
        if (j > 0 and a[j] < a[j - 1]) or (j < T.size - 1 and a[j] < a[j + 1]):
            return float('nan'), float(Td), "not on a peak"
        lo = hi = j                       # out to the minima either side, stopping at the first zero
        while lo > 0 and 0 < a[lo - 1] <= a[lo]:
            lo -= 1
        while hi < T.size - 1 and 0 < a[hi + 1] <= a[hi]:
            hi += 1
        notes = []
        if lo == 0 or hi == T.size - 1:
            notes.append("clipped by the grid")
        if hi - lo + 1 < 5:
            notes.append(f"{hi - lo + 1} points")
        # the fit also takes two zero points either side: a narrow peak is only a couple of grid
        # points wide, and three parameters need somewhere to pin the shoulders down
        seg = slice(max(lo - 2, 0), min(hi + 3, T.size))
        # the width is bounded below by half a grid step: a peak narrower than the grid is simply
        # unresolved, and without the bound the optimizer chases sigma to zero and never converges
        step = float(x[1] - x[0])
        # the fit runs in grid steps away from the peak and in units of the peak's own height, so
        # all three parameters are of order one whatever the signal units and however fine the grid.
        # Left in the raw units the three are of wildly different size, the trust region is set by
        # the largest of them, and the optimizer crawls along a bound until it runs out of tries
        height = max(float(np.max(np.abs(a[seg]))), 1e-300)
        u, v = (x[seg] - x[j])/step, a[seg]/height
        guess = [max(float(a[lo:hi + 1].max()), 1e-300)/height, 0.0,
                 max((x[hi] - x[lo])/4, step)/step]
        limits = ([0.0, float(u[0]), 0.5], [np.inf, float(u[-1]), float(x[-1] - x[0])/step])
        try:
            with warnings.catch_warnings():       # the note already says when the fit is thin
                warnings.simplefilter("ignore", OptimizeWarning)
                popt = curve_fit(lambda uu, A, u0, w: A*np.exp(-0.5*((uu - u0)/w)**2),
                                 u, v, p0=np.clip(guess, *limits), bounds=limits,
                                 maxfev=20000)[0]
            sigma = abs(float(popt[2]))*step
        except (RuntimeError, ValueError, TypeError) as fault:
            print(f"[EPRelax] the Gaussian on the peak at {Td:.4g} did not converge ({fault}); "
                  f"its width is the starting guess instead")
            sigma, notes = guess[2]*step, notes + ["no convergence"]
        weight = a[lo:hi + 1]
        centre = 10.0**float((x[lo:hi + 1]*weight).sum()/max(float(weight.sum()), 1e-300))
        return sigma, centre, "; ".join(notes)

    def discrete(self, dd, report=True):
        """Fit dd['ndisc'] exponentials, started at the largest maxima of the regularized
        distribution, and put two kinds of margin on each. The first is the 95 % interval of the
        non-linear fit's own trust region: how tightly the trace pins that component down. The
        second is the 95 % width of the distribution peak the component sits in, from a Gaussian
        in log10 T, kept as the factor the ends sit at either side of T -- which is the spread of
        relaxation times, not an error on T."""
        r = dd['res']
        if r is None:
            if report:
                wx.MessageBox("Nothing has been decomposed yet.", "Discrete fit", wx.OK | wx.ICON_INFORMATION)
            return False
        a, T, t, y = r['a'], r['T'], r['t'], r['y']
        n = dd['ndisc']
        # tallest peak first: that order is what decides which peak becomes which component
        peaks = [j for j in range(1, T.size - 1) if a[j] > a[j - 1] and a[j] >= a[j + 1]]
        peaks.sort(key=lambda j: -a[j])
        fixed = dd['fixT']
        if fixed and not peaks:
            if report:
                wx.MessageBox("The distribution has no peak to fix a time constant to.",
                              "Discrete fit", wx.OK | wx.ICON_INFORMATION)
            return False
        if fixed:
            n = min(n, len(peaks))        # a peak each: there is nothing to hold the extras at
            Tfix = np.array([self.peakwidth(T, a, float(T[j]))[1] for j in peaks[:n]])
        T0 = [float(T[j]) for j in peaks[:n]]
        T0 += list(np.logspace(np.log10(T[0]), np.log10(T[-1]), n + 2)[1:-1])[len(T0):]
        s = -1.0 if dd['recovery'] else 1.0
        free = dd['offset']
        amp = max(float(a.sum()), float(np.max(np.abs(y))))/n
        if fixed:
            model = lambda tt, *p: ((p[-1] if free else 0.0)
                                    + s*sum(p[i]*np.exp(-tt/Tfix[i]) for i in range(n)))
            p0 = [amp]*n + ([float(r['c'])] if free else [])
            bounds = ([0.0]*n + ([-np.inf] if free else []),
                      [np.inf]*n + ([np.inf] if free else []))
        else:
            model = lambda tt, *p: ((p[-1] if free else 0.0)
                                    + s*sum(p[i]*np.exp(-tt/p[n + i]) for i in range(n)))
            p0 = [amp]*n + sorted(T0) + ([float(r['c'])] if free else [])
            # a decade of slack either side of the grid: outside it the trace cannot argue about T
            bounds = ([0.0]*n + [float(T[0])*0.1]*n + ([-np.inf] if free else []),
                      [np.inf]*n + [float(T[-1])*10]*n + ([np.inf] if free else []))
        try:
            # the same weights the decomposition used, so both analyses answer the same question of
            # resampled data; sigma is relative, so curve_fit still scales the covariance itself
            popt, pcov = curve_fit(model, t, y, p0=p0, bounds=bounds, maxfev=40000,
                                   sigma=None if r['w'] is None else 1.0/r['w'])
        except (RuntimeError, ValueError) as fault:
            dd['disc'] = None
            print(f"[EPRelax] the {n}-exponential fit of '{dd['title']}' did not converge: "
                  f"{type(fault).__name__}: {fault}")
            if report:
                wx.MessageBox(f"The {n}-exponential fit did not converge:\n\n"
                              f"{type(fault).__name__}: {fault}",
                              "Discrete fit", wx.OK | wx.ICON_WARNING)
            return False
        fit = model(t, *popt)
        # curve_fit already scales the covariance by the residual variance, so this is the
        # linearised 95 % interval around the optimum, Student t with the fit's own d.o.f.
        dof = max(t.size - popt.size, 1)
        err = np.sqrt(np.abs(np.diag(pcov)))*float(tdist.ppf(0.975, dof))
        if fixed:
            # T was not a parameter, so the fit gives it no interval of its own. This is the one it
            # would have had if it had been free and this value were the optimum: the curvature of
            # the residual surface along T at the held point, with the amplitudes free to move with
            # it, so the correlation between a component's time constant and its amplitude is in
            # there. It says how tightly the trace pins T down around that value -- not whether the
            # value is the right one, which is the question holding it begged in the first place.
            Tout = Tfix
            amps = np.asarray(popt[:n])
            columns = [s*np.exp(-t/Tfix[i]) for i in range(n)]
            columns += [s*amps[i]*(t/Tfix[i]**2)*np.exp(-t/Tfix[i]) for i in range(n)]
            if free:
                columns.append(np.ones_like(t))
            jac = np.column_stack(columns)
            if r['w'] is not None:
                jac = jac*r['w'][:, None]
            asif = max(t.size - jac.shape[1], 1)
            try:
                res = (y - fit) if r['w'] is None else (y - fit)*r['w']
                cov = np.linalg.inv(jac.T @ jac)*float(np.sum(res*res))/asif
                dTout = np.sqrt(np.abs(np.diag(cov)[n:2*n]))*float(tdist.ppf(0.975, asif))
            except np.linalg.LinAlgError as fault:
                print(f"[EPRelax] the held time constants of '{dd['title']}' have no interval to "
                      f"estimate -- the components are not separable there: {fault}")
                dTout = np.full(n, np.inf)
        else:
            Tout = np.asarray(popt[n:2*n])
            dTout = np.asarray(err[n:2*n])
        order = np.argsort(Tout)

        # a held component is measured on the peak it was held at, not at the centre of mass that
        # came out of it: an asymmetric peak puts that centre a point or two off its own maximum
        lookup = np.array([float(T[j]) for j in peaks[:n]]) if fixed else Tout
        wide, note = [], []
        for Td, at in zip(Tout[order], lookup[order]):
            sigma, _, flag = self.peakwidth(T, a, at)
            wide.append(10.0**(TWO_SIDED_95*sigma))
            note.append("; ".join(x for x in ("T held" if fixed else "", flag) if x))
        dd['disc'] = {'T': Tout[order], 'a': np.asarray(popt[:n])[order],
                      'dT': dTout[order], 'da': np.asarray(err[:n])[order],
                      'c': float(popt[-1]) if free else 0.0, 'dc': float(err[-1]) if free else 0.0,
                      'rms': float(np.sqrt(np.mean((y - fit)**2))), 'dof': dof, 'Thold': fixed,
                      'sig': dd['sig'],        # the trace it was fitted on, to notice it moving
                      'wide': np.array(wide), 'note': note}
        if report:
            self.statusbar.SetStatusText(
                f"{n}-exponential fit of '{dd['title']}' converged"
                + (f", time constants held at {n} of the {len(peaks)} peaks" if fixed else ""), 0)
        return True


# --------------------------------------------------------------------------- #
# Application entry point
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame(None)
    frame.Show()
    app.MainLoop()
