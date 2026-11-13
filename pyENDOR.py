import os
import ast
import copy
import time
import base64
import traceback
import xml.etree.ElementTree as ET
from xml.dom import minidom

import wx
import wx.lib.mixins.listctrl as listmix
import wx.lib.scrolledpanel as scrolled
import wx.propgrid as wxpg

import numpy as np
from scipy.interpolate import CubicSpline

from matplotlib.backend_tools import Cursors
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

from theme import theme
from theme.theme import PillButton, PillCheckBox, PillTabBar, SegmentedPill
from theme import classPropGridPanel as MypgPanel
from readers import brukerread as Mybr
from spin_system.SysPar import sysPar
from spin_system.endor_sim import expENDOR, optENDOR, ENDORsim

# pip install wxPython numpy matplotlib scipy
VERSION = 0.1
G_FACTOR = 71.4477      # g = G_FACTOR * f[GHz] / B[mT]   (h/muB = 0.0714477 mT/MHz)
HERE = os.path.dirname(os.path.abspath(__file__))   # icons are found from here, never from the CWD
SEQUENCES = ['ideal', 'Mims', 'Davies']
# default colours of the per-nucleus simulation components, kept apart from the C0..C9 of the spectra
COMPONENT_COLORS = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']


class DictPopup(wx.Dialog):
    """Key/value view of a file header (Bruker DSC, session summary)."""
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


class AxisRangeDialog(wx.Dialog):
    """Min / max of one plot axis, or autoscale. After OK, `result` is (lo, hi), or None for autoscale."""
    def __init__(self, parent, title, lo, hi, auto, which):
        super().__init__(parent, title=title)
        self.result = None
        self.txt_lo = wx.TextCtrl(self, value=f"{lo:.6g}")
        self.txt_hi = wx.TextCtrl(self, value=f"{hi:.6g}")
        self.chk_auto = wx.CheckBox(self, label="Autoscale")
        self.chk_auto.SetValue(auto)
        self.chk_auto.Bind(wx.EVT_CHECKBOX, lambda e: [t.Enable(not self.chk_auto.GetValue())
                                                       for t in (self.txt_lo, self.txt_hi)])
        for t in (self.txt_lo, self.txt_hi):
            t.Enable(not auto)
        # laid out like the axis: x puts Min and Max side by side, y stacks Max above Min
        pairs = (("Min:", self.txt_lo), ("Max:", self.txt_hi)) if which == 'x' else \
                (("Max:", self.txt_hi), ("Min:", self.txt_lo))
        grid = wx.FlexGridSizer(cols=4 if which == 'x' else 2, vgap=6, hgap=8)
        for label, t in pairs:
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(t, 1, wx.EXPAND)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.chk_auto, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        sizer.Add(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL), 0, wx.EXPAND | wx.ALL, 6)
        self.SetSizerAndFit(sizer)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        theme.apply_theme_to_window(self, theme.get_theme())
        # open at the mouse, kept inside the screen it is on
        pos = wx.GetMousePosition()
        area = wx.Display(max(wx.Display.GetFromPoint(pos), 0)).GetClientArea()
        w, h = self.GetSize()
        self.Move(min(max(pos.x, area.x), area.Right - w), min(max(pos.y, area.y), area.Bottom - h))

    def on_ok(self, event):
        if not self.chk_auto.GetValue():
            try:
                lo, hi = sorted((float(self.txt_lo.GetValue()), float(self.txt_hi.GetValue())))
            except ValueError:
                wx.MessageBox("Min and Max must be numbers.", "Axis range", wx.OK | wx.ICON_WARNING)
                return
            if lo == hi:
                wx.MessageBox("Min and Max must differ.", "Axis range", wx.OK | wx.ICON_WARNING)
                return
            self.result = (lo, hi)
        self.EndModal(wx.ID_OK)


# --------------------------------------------------------------------------- #
# File Browser Panel
# --------------------------------------------------------------------------- #
class FileBrowserPanel(wx.Panel, listmix.ColumnSorterMixin):

    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.parent = mainWindow
        self.ini_file = os.path.join(wx.StandardPaths.Get().GetUserConfigDir(), "pyENDOR.ini")
        self.path = os.getcwd()
        self.filters = {"ENDOR": ('.dsc', '.xml'), "DSC": ('.dsc',), "XML": ('.xml',), "All": None}
        self.current_filter = "ENDOR"

        vbox = wx.BoxSizer(wx.VERTICAL)
        # Editable path; Enter or focus loss navigates to what was typed
        self.lbl_path = wx.TextCtrl(self, value="", style=wx.TE_PROCESS_ENTER | wx.TE_RIGHT | wx.BORDER_NONE)
        self.lbl_path.Bind(wx.EVT_TEXT_ENTER, self.on_path_entered)
        self.lbl_path.Bind(wx.EVT_KILL_FOCUS, self.on_path_entered)
        vbox.Add(self.lbl_path, 0, wx.EXPAND)

        hbox_filter = wx.BoxSizer(wx.HORIZONTAL)
        hbox_filter.Add(wx.StaticText(self, label="Filter:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.choice_filter = wx.ComboBox(self, choices=list(self.filters), value=self.current_filter,
                                         style=wx.CB_READONLY | wx.BORDER_NONE)
        self.choice_filter.Bind(wx.EVT_COMBOBOX, self.on_filter_changed)
        hbox_filter.Add(self.choice_filter, 1)
        vbox.Add(hbox_filter, 0, wx.EXPAND | wx.ALL, 5)

        self.list = wx.ListCtrl(self, style=wx.LC_REPORT | wx.LC_VRULES | wx.LC_HRULES | wx.BORDER_NONE)
        self.list.InsertColumn(0, "Name", width=300)
        self.list.InsertColumn(1, "Title", width=300)
        vbox.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)

        img_list = wx.ImageList(16, 16)
        icon = lambda name: wx.Bitmap(os.path.join(HERE, 'icons', name), wx.BITMAP_TYPE_ICO)
        self.folder_img_id = img_list.Add(icon('Folder.ico'))
        self.xml_img_id = img_list.Add(icon('Script.ico'))
        self.file_img_id = img_list.Add(icon('New document.ico'))
        self.data_img_id = img_list.Add(icon('List.ico'))
        self.loaded_img_id = img_list.Add(icon('Green bookmark.ico'))
        self.loadedBG_img_id = img_list.Add(icon('Equipment.ico'))
        self.sort_up_img_id = img_list.Add(wx.ArtProvider.GetBitmap(wx.ART_GO_UP, wx.ART_OTHER, (16, 16)))
        self.sort_down_img_id = img_list.Add(wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN, wx.ART_OTHER, (16, 16)))
        self.list.AssignImageList(img_list, wx.IMAGE_LIST_SMALL)

        btns = wx.BoxSizer(wx.VERTICAL)
        self.btn_load = PillButton(self, "Load", color_key="pill_load", draw_badge=theme.draw_squiggle_badge)
        self.btn_load.SetToolTip("Load the selected DSC file(s), or an XML session.\n"
                                 "Ctrl- or Shift-click selects several files.")
        self.btn_load.Bind(wx.EVT_BUTTON, self.on_load_clicked)
        btns.Add(self.btn_load, 0, wx.EXPAND | wx.ALL, 2)
        self.btn_loadbg = PillButton(self, "+ BG", color_key="pill_load", draw_badge=theme.draw_flatline_badge)
        self.btn_loadbg.SetToolTip("Load the selected file(s) as background spectra.")
        self.btn_loadbg.Bind(wx.EVT_BUTTON, self.on_loadBG_clicked)
        btns.Add(self.btn_loadbg, 0, wx.EXPAND | wx.ALL, 2)
        self.btn_header = PillButton(self, "Header", color_key="pill_calm", draw_badge=theme.draw_textlines_badge)
        self.btn_header.SetToolTip("Show the parameters stored in the selected file without loading it.")
        self.btn_header.Bind(wx.EVT_BUTTON, self.on_header_clicked)
        btns.Add(self.btn_header, 0, wx.EXPAND | wx.ALL, 2)
        vbox.Add(btns, 0, wx.EXPAND)
        self.SetSizer(vbox)

        # ColumnSorterMixin: two sortable columns, keys in self.itemDataMap
        listmix.ColumnSorterMixin.__init__(self, 2)
        self.items = []          # (name, full_path, is_dir); row data id = index here
        self.itemDataMap = {}

        self.parent.tb_change_folder.Bind(wx.EVT_BUTTON, self.on_change_folder)
        self.parent.tb_up.Bind(wx.EVT_BUTTON, self.on_up_clicked)
        self.parent.tb_save_session.Bind(wx.EVT_BUTTON, lambda e: self.parent.dumpXML())
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_item_activated)
        self.Bind(wx.EVT_WINDOW_DESTROY, lambda e: self.save_path())

        try:
            with open(self.ini_file, "r", encoding="utf-8") as f:
                stored = f.readline().strip()
            if os.path.isdir(stored):
                self.path = stored
        except OSError as fault:
            print(f"[pyENDOR] Could not read the last folder from {self.ini_file}: {fault}")
        self.lbl_path.SetValue(self.path)
        self.lbl_path.SetInsertionPointEnd()
        self.refresh_file_list()

    def save_path(self):
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.path + "\n")
        except OSError as e:
            print(f"[pyENDOR] Could not save last folder to {self.ini_file}: {e}")

    def refresh_file_list(self, select_name=None):
        """Rebuild the list; keeps `select_name` (a name or a list of names) selected if still there."""
        try:
            raw = os.listdir(self.path)
        except OSError as e:
            self.parent.statusbar.SetStatusText(f"Could not list folder: {e}", 0)
            raw = []
        exts = self.filters[self.current_filter]
        dirs = sorted([n for n in raw if os.path.isdir(os.path.join(self.path, n))], key=str.lower)
        files = sorted([n for n in raw if not os.path.isdir(os.path.join(self.path, n))
                        and (exts is None or n.lower().endswith(exts))], key=str.lower)
        # a Bruker .DSC/.DTA pair is one dataset
        same = lambda p: os.path.splitext(p)[0].lower() + ('.dsc' if p.lower().endswith(('.dsc', '.dta'))
                                                           else os.path.splitext(p)[1].lower())
        loaded = {same(dd['fullpath']) for dd in self.parent.Data}
        loadedBG = {same(dd['fullpath']) for dd in self.parent.BGData}

        self.list.DeleteAllItems()
        self.items = []
        self.itemDataMap = {}
        for name in dirs + files:
            full = os.path.join(self.path, name)
            is_dir = name in dirs
            title = ''
            low = name.lower()
            if is_dir:
                image_id = self.folder_img_id
            elif low.endswith('.xml'):
                image_id = self.xml_img_id
            elif low.endswith('.dsc'):
                image_id = self.data_img_id
                try:
                    with open(full, 'r', encoding='latin-1') as f:
                        for line in f:
                            if line.startswith('TITL'):
                                title = line.split("'")[1] if "'" in line else line[4:].strip()
                                break
                except (OSError, IndexError) as e:
                    print(f"[pyENDOR] Could not read title from '{full}': {e}")
            else:
                image_id = self.file_img_id
            if same(full) in loaded:
                image_id = self.loaded_img_id
            if same(full) in loadedBG:
                image_id = self.loadedBG_img_id
            idx = self.list.InsertItem(self.list.GetItemCount(), name, image_id)
            self.list.SetItem(idx, 1, title)
            self.items.append((name, full, is_dir))
            self.list.SetItemData(idx, len(self.items) - 1)
            prefix = "0_" if is_dir else "1_"    # folders stay on top whatever the sort
            self.itemDataMap[len(self.items) - 1] = (prefix + low, prefix + (title or name).lower())
        self.list.SetColumnWidth(0, -1)
        if self.GetSortState()[0] >= 0:      # re-apply the user's column sort, if any
            self.SortListItems()
        names = {select_name} if isinstance(select_name, str) else set(select_name or ())
        focused = False
        for ridx in range(self.list.GetItemCount()):
            if self.items[self.list.GetItemData(ridx)][0] in names:
                self.list.Select(ridx)
                if not focused:
                    self.list.Focus(ridx)
                    self.list.EnsureVisible(ridx)
                    focused = True

    # ColumnSorterMixin contract
    def GetListCtrl(self):
        return self.list

    def GetSortImages(self):
        return self.sort_down_img_id, self.sort_up_img_id

    def selected_items(self):
        """[(name, full_path, is_dir), ...] of the selected rows in list order; [] with a warning if none."""
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
                self.parent.loadData(full, BG=False)

    def on_loadBG_clicked(self, event):
        for name, full, is_dir in self.selected_items():
            if not is_dir:
                self.parent.loadData(full, BG=True)

    def on_header_clicked(self, event):
        items = self.selected_items()
        if not items or items[0][2]:
            return
        name, full, _ = items[0]
        ext = os.path.splitext(name)[1].lower()
        try:
            if ext in ('.dsc', '.dta'):
                info = Mybr.parse_dsc(os.path.splitext(full)[0] + '.dsc')
            elif ext == '.xml':
                result = self.parent.load_mixed_dict_xml(full)
                info = {}
                for key, val in result.items():
                    if isinstance(val, list):
                        for ii, dd in enumerate(val):
                            info[f"{key} {ii + 1}"] = dd.get('title', '?') if isinstance(dd, dict) else type(dd).__name__
                    else:
                        info[key] = val if isinstance(val, (str, int, float)) else type(val).__name__
            else:
                wx.MessageBox(f"No header reader for '{ext}' files.", "Header", wx.OK | wx.ICON_INFORMATION)
                return
        except Exception as e:
            traceback.print_exc()
            wx.MessageBox(f"Could not read the header of '{name}':\n{type(e).__name__}: {e}",
                          "Header", wx.OK | wx.ICON_ERROR)
            return
        if not info:
            wx.MessageBox(f"'{name}' holds no header entries.", "Header", wx.OK | wx.ICON_INFORMATION)
            return
        with DictPopup(self, info, title=name) as dlg:
            dlg.ShowModal()

    def on_filter_changed(self, event):
        self.current_filter = self.choice_filter.GetStringSelection()
        self.refresh_file_list()

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
            self.parent.loadData(full, BG=False)

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
# Matplotlib Canvas Panel
# --------------------------------------------------------------------------- #
class MatplotlibPanel(wx.Panel):
    """One plot for everything: every shown spectrum overlaid, each with its own simulation on top.
    The current spectrum is drawn thicker and carries the raw/background/baseline extras."""
    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.parent = mainWindow
        self.showRawBG = False
        self.showBaselines = True
        self.showSim = True
        self.showParts = True     # one trace per simulated nucleus, next to the summed simulation
        self.showGrid = False
        self.xnuc = ''            # isotope whose Larmor frequency is the zero of the x axis; '' = plain frequency axis
        self.norm = 'none'        # 'none' | 'max'
        self.simfit = 'max'       # simulation amplitude: 'max' = match the maxima, 'lsq' = least squares
        self.offset = 0.0         # stacking step, fraction of the tallest trace
        self.keepx = False        # the x-zoom of the last draw survives the next one
        self.drag = None          # line / region drag in progress (see on_press)
        self.artists = {}         # region index -> grey patch of the current spectrum's baseline region
        self.height = 1.0         # tallest trace of the last draw: unit of the stack step and of every Y shift
        self.lineoff = {}         # (axes, spectrum index) -> vertical offset of that line
        self.xfix = {}            # x-axis mode (the xnuc it was set in) -> fixed (lo, hi); absent = autoscale
        self.yfix = None          # fixed (lo, hi) of the y axis, None = autoscale
        self.zoom_once = False    # Ctrl-Z zoom mode: leave it after the next zoom

        self.figure = Figure(figsize=(2, 2), dpi=100, layout='constrained')
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('figure_leave_event', lambda e: self.parent.statusbar.SetStatusText("", 1))
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        mainWindow.Bind(wx.EVT_CHAR_HOOK, self.on_key)

        self.chk_RawBG = PillCheckBox(self, "Raw & BG")
        self.chk_RawBG.SetValue(self.showRawBG)
        self.chk_Baselines = PillCheckBox(self, "Baselines")
        self.chk_Baselines.SetValue(self.showBaselines)
        self.chk_Sim = PillCheckBox(self, "Simulation")
        self.chk_Sim.SetValue(self.showSim)
        self.chk_Parts = PillCheckBox(self, "Components")
        self.chk_Parts.SetValue(self.showParts)
        sizerH = wx.WrapSizer(wx.HORIZONTAL)
        for chk in (self.chk_RawBG, self.chk_Baselines, self.chk_Sim, self.chk_Parts):
            chk.Bind(wx.EVT_CHECKBOX, self.on_check)
            sizerH.Add(chk, 0, wx.EXPAND)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerH, 0, wx.EXPAND)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()
        self.SetSizer(sizer)

        self.redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self.redraw_debouncer.on_size)

    def on_key(self, event):
        """Ctrl-Z (anywhere but a text box, which keeps its undo): zoom-to-rectangle for one zoom only."""
        if (event.GetKeyCode() == ord('Z') and event.GetModifiers() == wx.MOD_CONTROL
                and not isinstance(wx.Window.FindFocus(), wx.TextCtrl)):
            if str(self.toolbar.mode) != 'zoom rect':
                self.toolbar.zoom()
            self.zoom_once = True
            return
        event.Skip()

    def larmor(self, dd):
        """Larmor frequency [MHz] of the nucleus the x axis is centred on, at this spectrum's own field
        (the Field of 'Exp Over.' when that is set); 0 on the plain frequency axis. Plot x = RF - larmor."""
        mf = self.parent
        if not self.xnuc:
            return 0.0
        field = mf.tabulated_panel.Exp_param.parameters['Field']
        return mf.Sys.isotopes(self.xnuc)[2][0]*(field if field >= 0 else dd['field'])

    def hit(self, event):
        """Baseline region of the current spectrum under the mouse, as (j, part): part 'lo'/'hi' = an edge
        within 5 px, 'move' = inside. None if there is none."""
        mf = self.parent
        if event.inaxes is None or mf.current < 0 or not self.showBaselines:
            return None
        sh = self.larmor(mf.Data[mf.current])
        px = lambda f: event.inaxes.transData.transform((f - sh, 0))[0]
        regs = mf.Data[mf.current]['blregions']
        for j, (lo, hi) in enumerate(regs):
            for part, f in (('lo', lo), ('hi', hi)):
                if abs(px(f) - event.x) <= 5:
                    return (j, part)
        for j, (lo, hi) in enumerate(regs):
            if min(px(lo), px(hi)) <= event.x <= max(px(lo), px(hi)):
                return (j, 'move')
        return None

    def on_mouse_move(self, event):
        mf = self.parent
        if (self.drag is not None and self.drag['key'] in ('line', 'scale')
                and (self.drag['ax'], self.drag['ii']) not in self.lineoff):
            print(f"[pyENDOR] The plot was rebuilt during a {self.drag['key']} drag of spectrum "
                  f"{self.drag['ii'] + 1}: its axes and offsets are gone, so the drag was dropped")
            self.drag = None
        if self.drag is not None and self.drag['key'] == 'line':
            # Vertical line drag: the spectrum and its simulation move by the same number of tallest-trace units
            y = self.drag['ax'].transData.inverted().transform((event.x, event.y))[1]
            d = (y - self.drag['y0'])/self.height
            mf.Data[self.drag['ii']]['yshift'] = float(self.drag['orig'] + d)
            for ln, base in self.drag['lines']:
                ln.set_ydata(base + d*self.height)
            self.canvas.draw_idle()
        elif self.drag is not None and self.drag['key'] == 'scale':
            # Scale drag: the grabbed point follows the mouse, the line's zero (its offset) stays put.
            # Grabbed near zero there is nothing to follow: then up = larger, one plot height = x e^3
            ax0, ii = self.drag['ax'], self.drag['ii']
            y = ax0.transData.inverted().transform((event.x, event.y))[1]
            off, r0 = self.lineoff[(ax0, ii)], self.drag['y0'] - self.lineoff[(ax0, ii)]
            if abs(r0) > 0.02*self.height:
                factor = max((y - off)/r0, 0.01)
            else:
                factor = np.exp(3*(y - self.drag['y0'])/self.height)
            mf.Data[ii]['scale'] = float(self.drag['orig']*factor)
            for ln, base in self.drag['lines']:
                ln.set_ydata(off + (base - off)*factor)
            mf.statusbar.SetStatusText(f"Scale = {mf.Data[ii]['scale']:.4g}", 0)
            self.canvas.draw_idle()
        elif self.drag is not None and event.inaxes is not None and event.xdata is not None:
            # Live region drag: move the stored edge(s), redraw only the dragged patch
            j, part = self.drag['j'], self.drag['part']
            sh = self.larmor(mf.Data[mf.current])
            f = float(event.xdata) + sh
            reg = mf.Data[mf.current]['blregions'][j]
            if part == 'move':
                d = f - self.drag['start']
                reg[:] = [self.drag['orig'][0] + d, self.drag['orig'][1] + d]
            else:
                reg[0 if part == 'lo' else 1] = f
            ax = self.artists[j].axes
            self.artists[j].remove()
            self.artists[j] = ax.axvspan(reg[0] - sh, reg[1] - sh, color='0.5', alpha=0.15, lw=0)
            self.canvas.draw_idle()
        elif not str(self.toolbar.mode) and event.inaxes is not None:
            h = self.hit(event)
            if (h is None or h[1] == 'move') and any(str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]
                                                     for ln in event.inaxes.get_lines()):
                self.canvas.set_cursor(Cursors.RESIZE_VERTICAL)     # a spectrum line: drag up / down
            else:
                self.canvas.set_cursor(Cursors.POINTER if h is None else
                                       Cursors.MOVE if h[1] == 'move' else Cursors.RESIZE_HORIZONTAL)
        if event.inaxes is None or event.xdata is None:
            mf.statusbar.SetStatusText("", 1)
            return
        if self.xnuc:       # the shift differs between spectra: the plain frequency given is the current spectrum's
            text = f"RF − ν({self.xnuc}) = {event.xdata:.4f}" + (
                f",  RF = {event.xdata + self.larmor(mf.Data[mf.current]):.4f}" if mf.current >= 0 else "")
        else:
            text = f"RF = {event.xdata:.4f}"
        mf.statusbar.SetStatusText(f"{text} MHz,  y = {event.ydata:.4g}", 1)

    def on_press(self, event):
        """Double click on a spectrum line: make it the current spectrum; elsewhere: Home view.
        Left drag on a baseline region of the current spectrum: inside = move, edge = resize.
        Drag on a spectrum line: move it up / down (Y shift); shift-drag on it: scale it.
        Shift-drag on empty space adds a baseline region."""
        mf = self.parent
        if event.button != 1 or event.inaxes is None or event.xdata is None or mf.current < 0:
            return
        if event.dblclick:
            self.drag = None             # the first click of the pair may have started one
            hits = [ln for ln in event.inaxes.get_lines()
                    if str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]]
            # where lines overlap, prefer one that is not current yet (the current one is drawn on top)
            hits = [h for h in hits if h.get_gid() != f'spec{mf.current}'] or hits
            if hits:
                ii = int(max(hits, key=lambda h: h.get_zorder()).get_gid()[4:])
                if ii != mf.current:
                    mf.set_current(ii)
            else:
                self.toolbar.home()
            return
        if str(self.toolbar.mode):       # pan / zoom owns plain drags, but not the double click above
            return
        dd = mf.Data[mf.current]
        f = float(event.xdata) + self.larmor(dd)
        h = self.hit(event)
        if h is None or h[1] == 'move':
            # a spectrum line (topmost) under the cursor beats the inside of a region
            hits = [ln for ln in event.inaxes.get_lines()
                    if str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]]
            if hits:
                ii = int(max(hits, key=lambda h: h.get_zorder()).get_gid()[4:])
                scale = wx.GetKeyState(wx.WXK_SHIFT)
                self.drag = {'key': 'scale' if scale else 'line', 'ii': ii, 'ax': event.inaxes, 'y0': event.ydata,
                             'orig': mf.Data[ii]['scale' if scale else 'yshift'],
                             'lines': [(ln, np.asarray(ln.get_ydata(), dtype=float).copy())
                                       for ln in event.inaxes.get_lines()
                                       if ln.get_gid() in (f'spec{ii}', f'sim{ii}', f'part{ii}')]}
                return
        new = h is None
        if new:
            if not self.showBaselines or not wx.GetKeyState(wx.WXK_SHIFT):
                return
            dd['blregions'].append([f, f])
            h = (len(dd['blregions']) - 1, 'hi')
            self.update_graph()          # creates the patch the drag then reshapes
        j, part = h
        self.drag = {'key': 'region', 'j': j, 'part': part, 'start': f, 'new': new, 'orig': list(dd['blregions'][j])}

    def on_release(self, event):
        """Right click: on an axis (ticks / label) its range dialog; on a spectrum line the colour menu;
        on a baseline region deletes it. Left release ends a drag."""
        mf = self.parent
        if self.zoom_once and str(self.toolbar.mode) != 'zoom rect':
            self.zoom_once = False       # zoom was switched off on the plot toolbar meanwhile
        if self.zoom_once and event.inaxes is not None:
            # the toolbar finishes this zoom after us (connected later), so leave zoom mode once it has
            self.zoom_once = False
            wx.CallAfter(lambda: str(self.toolbar.mode) == 'zoom rect' and self.toolbar.zoom())
            return
        if event.button == 3 and event.inaxes is None and not str(self.toolbar.mode) and self.figure.axes:
            ax = self.figure.axes[0]
            renderer = self.canvas.get_renderer()
            for which, axis in (('x', ax.xaxis), ('y', ax.yaxis)):
                bb = axis.get_tightbbox(renderer)
                if bb is None or not bb.contains(event.x, event.y):
                    continue
                lo, hi = sorted(ax.get_xlim() if which == 'x' else ax.get_ylim())
                label = ax.get_xlabel() if which == 'x' else ax.get_ylabel()
                auto = self.xnuc not in self.xfix if which == 'x' else self.yfix is None
                with AxisRangeDialog(self, f"{which.upper()} axis: {label}", lo, hi, auto, which) as dlg:
                    if dlg.ShowModal() != wx.ID_OK:
                        return
                    rng = dlg.result
                if which == 'y':
                    self.yfix = rng
                elif rng is None:
                    self.xfix.pop(self.xnuc, None)
                else:
                    self.xfix[self.xnuc] = rng        # kept separately for every choice of the x-axis zero
                self.update_graph()
                return
            return
        if event.button == 3 and event.inaxes is not None and not str(self.toolbar.mode) and mf.current >= 0:
            region = self.hit(event)
            hits = [ln for ln in event.inaxes.get_lines()
                    if str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]]
            if not hits:
                if region is not None:
                    del mf.Data[mf.current]['blregions'][region[0]]
                    mf.tabulated_panel.update_controls()
                    self.update_graph()
                return
            # A spectrum line under the cursor (topmost wins): its colour, plus delete if inside a region
            ln = max(hits, key=lambda h: h.get_zorder())
            dd = mf.Data[int(ln.get_gid()[4:])]
            menu = wx.Menu()
            item_pick = menu.Append(wx.ID_ANY, f"Colour of '{dd['title']}'…")
            item_auto = menu.Append(wx.ID_ANY, "Default colour")
            item_auto.Enable(bool(dd['color']))
            id_del = None
            if region is not None:
                menu.AppendSeparator()
                id_del = menu.Append(wx.ID_ANY, "Delete baseline region").GetId()
            id_pick, id_auto = item_pick.GetId(), item_auto.GetId()    # the items die with the menu
            choice = self.canvas.GetPopupMenuSelectionFromUser(menu)
            menu.Destroy()
            if choice == id_del:
                del mf.Data[mf.current]['blregions'][region[0]]
                mf.tabulated_panel.update_controls()
            elif choice == id_pick:
                cdata = wx.ColourData()
                cdata.SetColour(wx.Colour(to_hex(ln.get_color())))
                with wx.ColourDialog(self, cdata) as dlg:
                    if dlg.ShowModal() != wx.ID_OK:
                        return
                    dd['color'] = dlg.GetColourData().GetColour().GetAsString(wx.C2S_HTML_SYNTAX)
            elif choice == id_auto:
                dd['color'] = ''
            else:
                return
            self.update_graph()
            return
        if event.button != 1 or self.drag is None:
            return
        if self.drag['key'] in ('line', 'scale'):
            key = 'yshift' if self.drag['key'] == 'line' else 'scale'
            ii, moved = self.drag['ii'], mf.Data[self.drag['ii']][key] != self.drag['orig']
            self.drag = None
            if moved:            # a plain click leaves everything as it was
                mf.tabulated_panel.filelist.set_value(ii, key, mf.Data[ii][key])
                self.update_graph()
            return
        j, new = self.drag['j'], self.drag['new']
        self.drag = None
        dd = mf.Data[mf.current]
        regs = dd['blregions']
        regs[j].sort()
        if new and regs[j][1] - regs[j][0] < 0.002*np.ptp(dd['ax']['x']):    # shift-click, not a drag
            del regs[j]
        mf.tabulated_panel.update_controls()
        self.update_graph()

    def on_check(self, event):
        self.showRawBG = self.chk_RawBG.GetValue()
        self.showBaselines = self.chk_Baselines.GetValue()
        self.showSim = self.chk_Sim.GetValue()
        self.showParts = self.chk_Parts.GetValue()
        self.update_graph()

    def update_graph(self):
        mf = self.parent
        keep = self.keepx and len(self.figure.axes) == 1
        xlim = self.figure.axes[0].get_xlim() if keep else None
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        self.keepx = True
        self.artists = {}
        self.lineoff = {}

        procs = [mf.process(dd) for dd in mf.Data]
        warnings = sorted({w for p in procs for w in p['warn']})
        cur = mf.current
        markcol = theme.theme_colors()["ax_fg"]      # simulations and fitted baselines: the axis colour of the theme

        # display factor of every shown spectrum, then the tallest scaled trace
        facs = {}
        for ii, (dd, p) in enumerate(zip(mf.Data, procs)):
            if not dd['show']:
                continue
            facs[ii] = dd['scale']
            if self.norm == 'max' and np.max(np.abs(p['y'])) > 0:
                facs[ii] /= np.max(np.abs(p['y']))
        self.height = max([np.ptp(procs[ii]['y']*fac) for ii, fac in facs.items()], default=0.0) or 1.0

        nucnames = [val['Nucs'][0] for key, val in mf.tabulated_panel.Sys_param.parameters.items() if 'nuc' in key]
        labelled = set()                 # legend entries shared by all spectra appear once
        for kk, (ii, fac) in enumerate(facs.items()):
            dd, p = mf.Data[ii], procs[ii]
            col = mf.trace_color(dd)
            off = (kk*self.offset + dd['yshift'])*self.height
            self.lineoff[(ax, ii)] = off
            xp = p['x'] - self.larmor(dd)      # the spectrum, its simulation and everything else of it move together
            ax.plot(xp, p['y']*fac + off, color=col, lw=1.6 if ii == cur else 0.9, gid=f'spec{ii}',
                    label=f"{ii + 1}: {dd['title']}", zorder=3 if ii == cur else 2)
            if self.showSim and dd['simdata'] is not None:
                # components of the nuclei still in the spin system, each times its scale; their sum is the simulation
                parts = [(k, mf.components[k], np.interp(p['x'], dd['simax']['x'], part, left=0.0, right=0.0))
                         for k, part in enumerate(dd['simparts'] or []) if part is not None and k < len(mf.components)]
                sim = sum(comp['scale']*part for k, comp, part in parts) if parts else \
                    np.interp(p['x'], dd['simax']['x'], dd['simdata'], left=0.0, right=0.0)
                amp = 0.0
                if self.simfit == 'lsq' and np.sum(sim**2) > 0:
                    amp = np.sum(p['y']*sim)/np.sum(sim**2)
                elif self.simfit == 'max' and np.max(np.abs(sim)) > 0:
                    # signed: inverted spectra stay inverted
                    amp = p['y'][np.argmax(np.abs(p['y']))]/sim[np.argmax(np.abs(sim))]
                amp = amp*dd['simscale']
                ax.plot(xp, amp*sim*fac + off, color=markcol, lw=1.0, gid=f'sim{ii}',
                        label='_nolegend_' if 'sim' in labelled else 'simulation', zorder=5)
                labelled.add('sim')
                for k, comp, part in (parts if self.showParts else []):
                    name = f"{nucnames[k] if k < len(nucnames) else '?'} ({k + 1})"
                    # dashed and on top: where one nucleus carries the whole simulation, both lines stay visible
                    ax.plot(xp, amp*comp['scale']*part*fac + off, color=comp['color'], lw=1.0, ls='--', gid=f'part{ii}',
                            label='_nolegend_' if k in labelled else name, zorder=6)
                    labelled.add(k)
            if ii != cur:
                continue
            if self.showRawBG:
                ax.plot(xp, p['raw']*fac + off, color='0.6', lw=0.8, label='raw')
                if dd['bglist']:
                    ax.plot(xp, p['bg']*fac + off, color='0.4', lw=0.8, ls='--', label='background')
            if self.showBaselines and dd['blorder'] >= 0:
                ax.plot(xp, p['y0']*fac + off, color=col, lw=0.6, alpha=0.6)
                ax.plot(xp, p['bl']*fac + off, color=markcol, lw=0.8, ls=':', label='baseline')
            if self.showBaselines:
                for j, (lo, hi) in enumerate(dd['blregions']):
                    self.artists[j] = ax.axvspan(lo - self.larmor(dd), hi - self.larmor(dd), color='0.5', alpha=0.15, lw=0)
        ax.grid(self.showGrid)
        ax.set_xlabel(f"RF − ν({self.xnuc}), MHz" if self.xnuc else "RF, MHz")
        ax.set_ylabel("ENDOR intensity")
        if facs:
            ax.legend(fontsize='small', frameon=False, loc='best', labelcolor=markcol)
        else:
            ax.text(0.5, 0.5, "Load a spectrum from the file browser", ha='center', va='center',
                    transform=ax.transAxes)
        theme.style_figure(self.figure, [ax])

        # Home = the full autoscaled view. A fixed x range (axis right-click) or a kept x-zoom limits x, and
        # y autoscales to what is visible -- unless the y range is fixed
        self.toolbar.update()
        self.toolbar.push_current()
        xfix = self.xfix.get(self.xnuc)
        if xfix is not None:
            xlim = xfix
        if xfix is not None or keep:
            ax.set_xlim(xlim)
            xlo, xhi = sorted(xlim)
            ys = [np.asarray(ln.get_ydata())[(np.asarray(ln.get_xdata()) >= xlo) & (np.asarray(ln.get_xdata()) <= xhi)]
                  for ln in ax.get_lines()]
            ys = np.concatenate(ys) if ys else np.array([])
            if ys.size and np.ptp(ys) > 0:
                ax.set_ylim(ys.min() - 0.05*np.ptp(ys), ys.max() + 0.05*np.ptp(ys))
        if self.yfix is not None:
            ax.set_ylim(self.yfix)
        self.canvas.draw_idle()

        mf.tabulated_panel.filelist.refresh_headers()     # a colour picked on the plot shows on its swatch too
        mf.statusbar.SetStatusText("; ".join(warnings) if warnings else "Ready", 0)
        if mf.epr_window is not None:        # its markers are these spectra; its curve follows the same parameters
            mf.epr_window.refresh()


# --------------------------------------------------------------------------- #
# Data tab spectrum list: one collapsible entry per spectrum, accordion style
# --------------------------------------------------------------------------- #
class SpectrumHeader(wx.Panel):
    """Title bar of one spectrum entry: [show] [colour] #: title [x], painted by hand.
    Clicks are reported to the owning SpectrumList as (ii, part, event)."""

    PARTS_TIPS = {'show': "Show / hide this spectrum",
                  'color': "Line colour: click to pick",
                  'delete': "Remove this spectrum"}

    def __init__(self, parent, owner, ii):
        super().__init__(parent, style=wx.BORDER_NONE)
        self.owner, self.ii = owner, ii
        self.hover = None
        self.SetMinSize((-1, self.FromDIP(24)))
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Refresh(), e.Skip()))
        self.Bind(wx.EVT_MOTION, self.on_motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, lambda e: self.set_hover(None))
        self.Bind(wx.EVT_LEFT_UP, lambda e: self.owner.on_header(self.ii, self.part_at(e.GetPosition()), e))

    def rects(self):
        w, h = self.GetClientSize()
        pad, box = self.FromDIP(6), self.FromDIP(15)
        y = (h - box)//2
        show = wx.Rect(pad, y, box, box)
        color = wx.Rect(show.GetRight() + self.FromDIP(6), y, self.FromDIP(22), box)
        delete = wx.Rect(w - pad - box, y, box, box)
        title = wx.Rect(color.GetRight() + self.FromDIP(8), 0, max(0, delete.GetLeft() - color.GetRight() - self.FromDIP(14)), h)
        return {'show': show, 'color': color, 'delete': delete, 'title': title}

    def part_at(self, pos):
        for part, r in self.rects().items():
            if part != 'title' and r.Inflate(2, 2).Contains(pos):
                return part
        return 'title'

    def set_hover(self, part):
        if part != self.hover:
            self.hover = part
            self.SetToolTip(self.PARTS_TIPS.get(part, self.owner.title(self.ii)))
            self.Refresh()

    def on_motion(self, event):
        self.set_hover(self.part_at(event.GetPosition()))
        event.Skip()

    def on_paint(self, event):
        t = theme.theme_colors()
        dd = self.owner.mf.Data[self.ii]
        w, h = self.GetClientSize()
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        gc = wx.GraphicsContext.Create(dc)
        if gc is None or w <= 0 or h <= 0:
            return
        accent, bg = wx.Colour(t["accent"]), wx.Colour(t["win_bg"])
        current = self.ii == self.owner.mf.current
        fill = accent if current else PillCheckBox._blend(bg, accent, 0.55 if self.hover == 'title' else 0.45)
        white = wx.Colour("#FFFFFF")
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.SetBrush(wx.Brush(fill))
        gc.DrawRoundedRectangle(0, 0, w, h, self.FromDIP(5))
        r = self.rects()

        # show checkbox: white box, accent tick when on
        b = r['show']
        if dd['show']:
            gc.SetBrush(wx.Brush(white))
            gc.DrawRoundedRectangle(b.x, b.y, b.width, b.height, self.FromDIP(4))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(fill).Width(max(1.4, b.width*0.14))
                                   .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
            tick = gc.CreatePath()
            tick.MoveToPoint(b.x + b.width*0.24, b.y + b.height*0.52)
            tick.AddLineToPoint(b.x + b.width*0.42, b.y + b.height*0.72)
            tick.AddLineToPoint(b.x + b.width*0.78, b.y + b.height*0.28)
            gc.StrokePath(tick)
        else:
            gc.SetBrush(wx.TRANSPARENT_BRUSH)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(white).Width(1.4)))
            gc.DrawRoundedRectangle(b.x + 1, b.y + 1, b.width - 2, b.height - 2, self.FromDIP(4))

        # colour swatch
        c = r['color']
        gc.SetBrush(wx.Brush(wx.Colour(self.owner.color(self.ii))))
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(white).Width(2.0 if self.hover == 'color' else 1.2)))
        gc.DrawRoundedRectangle(c.x + 1, c.y + 1, c.width - 2, c.height - 2, self.FromDIP(3))

        # delete cross, with a soft disc behind it on hover
        d = r['delete']
        if self.hover == 'delete':
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(t["pill_warn"])))
            gc.DrawEllipse(d.x - 1, d.y - 1, d.width + 2, d.height + 2)
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(white).Width(1.6).Cap(wx.CAP_ROUND)))
        m = d.width*0.3
        gc.StrokeLine(d.x + m, d.y + m, d.GetRight() + 1 - m, d.GetBottom() + 1 - m)
        gc.StrokeLine(d.x + m, d.GetBottom() + 1 - m, d.GetRight() + 1 - m, d.y + m)

        # title, shortened with an ellipsis to fit
        tr = r['title']
        font = self.GetFont()
        if current:
            font = font.Bold()
        gc.SetFont(font, white)
        text = self.owner.title(self.ii)
        if gc.GetTextExtent(text)[0] > tr.width:
            while text and gc.GetTextExtent(text + "…")[0] > tr.width:
                text = text[:-1]
            text += "…"
        th = gc.GetTextExtent(text)[1]
        gc.DrawText(text, tr.x, (h - th)/2)


class SpectrumList(scrolled.ScrolledPanel):
    """Accordion of spectra: a SpectrumHeader per spectrum, and below it a small PropertyGrid
    (properties named "<key> <index>") that only the open entry shows."""

    def __init__(self, parent, tabulated):
        super().__init__(parent, style=wx.BORDER_NONE)
        self.tp = tabulated
        self.mf = tabulated.parent
        self.headers, self.grids = [], []
        self.open = None                # index of the expanded entry, None = all collapsed
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)
        self.SetupScrolling(scroll_x=False, rate_y=10)

    def title(self, ii):
        return f"{ii + 1}: {self.mf.Data[ii]['title']}"

    def color(self, ii):
        return to_hex(self.mf.trace_color(self.mf.Data[ii]))

    def rebuild(self):
        self.Freeze()
        try:
            self.sizer.Clear(delete_windows=True)
            self.headers, self.grids = [], []
            for ii, dd in enumerate(self.mf.Data):
                header = SpectrumHeader(self, self, ii)
                self.sizer.Add(header, 0, wx.EXPAND | wx.TOP, 2)
                grid = wxpg.PropertyGrid(self, style=wxpg.PG_DEFAULT_STYLE | wxpg.PG_SPLITTER_AUTO_CENTER
                                         | wxpg.PG_HIDE_MARGIN | wxpg.PG_TOOLTIPS | wxpg.PG_NO_INTERNAL_BORDER
                                         | wx.BORDER_NONE)
                if not grid.GetEditorByName("SpinFloat"):
                    grid.RegisterEditor(MypgPanel.SpinCtrlDoubleEditor, "SpinFloat")
                prop = grid.Append(wxpg.StringProperty("f.name", f"fname {ii}", value=dd['fname']))
                grid.SetPropertyReadOnly(prop)
                prop = grid.Append(wxpg.FloatProperty("MW Freq [GHz]", f"freq {ii}", value=dd['freq']))
                prop.SetHelpString("Microwave frequency; with B0 it decides which orientations are excited.")
                prop = grid.Append(wxpg.FloatProperty("B0 [mT]", f"field {ii}", value=dd['field']))
                prop.SetHelpString("Static magnetic field the spectrum was recorded at.")
                seq = dd['seq'] if dd['seq'] in SEQUENCES else 'ideal'
                prop = grid.Append(wxpg.EnumProperty("Sequence", f"seq {ii}", labels=SEQUENCES,
                                                     value=SEQUENCES.index(seq)))
                prop.SetHelpString("ENDOR sequence of this spectrum, used by the simulation\n"
                                   "unless Sequence in 'Exp Over.' overrides it.")
                prop = grid.Append(wxpg.FloatProperty("tau [ns]", f"tau {ii}", value=dd['tau']))
                prop.SetHelpString("Mims ENDOR: delay between the first two microwave pulses.")
                prop = grid.Append(wxpg.FloatProperty("t.inv [ns]", f"tinv {ii}", value=dd['tinv']))
                prop.SetHelpString("Davies ENDOR: length of the inversion pulse.")
                for key, label, helpstr in (
                        ('scale', "Scale", "Display multiplier (applied on top of the Compare normalisation).\n"
                                           "Or shift-drag the spectrum's line on the plot."),
                        ('yshift', "Y shift", "Vertical display offset, in units of the tallest trace.\n"
                                              "Or drag the spectrum's line up or down on the plot."),
                        ('simscale', "Sim scale", "Extra multiplier of this spectrum's simulation, on top of the\n"
                                                  "amplitude chosen by 'Simulation scaling' in the Options tab.")):
                    prop = grid.Append(wxpg.FloatProperty(label, f"{key} {ii}", value=dd[key]))
                    grid.SetPropertyEditor(prop, "SpinFloat")
                    prop.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.1)
                    prop.SetHelpString(helpstr)
                grid.SetMinSize((-1, grid.GetRowHeight()*grid.GetRoot().GetChildCount() + 2))
                grid.Bind(wxpg.EVT_PG_CHANGED, self.tp.on_filetree)
                grid.Bind(wxpg.EVT_PG_SELECTED, self.tp.on_filetree_selected)
                grid.Bind(wx.EVT_MOUSEWHEEL, self.on_wheel)    # the grid fits its rows; scroll the list instead
                theme.theme_propgrid(grid, caption_key="accent")
                self.sizer.Add(grid, 0, wx.EXPAND | wx.LEFT, self.FromDIP(6))
                self.headers.append(header)
                self.grids.append(grid)
            self.open = self.mf.current if 0 <= self.mf.current < len(self.grids) else None
            for ii, grid in enumerate(self.grids):
                grid.Show(ii == self.open)
            self.SetupScrolling(scroll_x=False, rate_y=10, scrollToTop=False)
        finally:
            self.Thaw()          # never leave the list frozen: it would stop repainting for good

    def expand(self, ii):
        """Open entry ii (None: none) and close the others."""
        self.open = ii
        self.Freeze()
        try:
            for jj, grid in enumerate(self.grids):
                grid.Show(jj == ii)
            self.SetupScrolling(scroll_x=False, rate_y=10, scrollToTop=False)
        finally:
            self.Thaw()              # never leave the list frozen: it would stop repainting for good
        if ii is not None:
            self.ScrollChildIntoView(self.headers[ii])

    def set_current(self, ii):
        if ii != self.open:
            self.expand(ii)
        self.refresh_headers()

    def refresh_headers(self):
        for header in self.headers:
            header.Refresh()

    def set_value(self, ii, key, value):
        if 0 <= ii < len(self.grids):
            self.grids[ii].SetPropertyValue(f"{key} {ii}", value)

    def apply_theme(self):
        for grid in self.grids:
            theme.theme_propgrid(grid, caption_key="accent")
        self.refresh_headers()

    def on_wheel(self, event):
        lines = -event.GetWheelRotation()/event.GetWheelDelta()*event.GetLinesPerAction()
        x, y = self.GetViewStart()
        self.Scroll(x, y + int(round(lines)))

    def on_header(self, ii, part, event):
        mf = self.mf
        dd = mf.Data[ii]
        if part == 'show':
            dd['show'] = not dd['show']
            self.headers[ii].Refresh()
            mf.matplotlib_panel.update_graph()
        elif part == 'color':
            cdata = wx.ColourData()
            cdata.SetColour(wx.Colour(self.color(ii)))
            with wx.ColourDialog(self, cdata) as dlg:
                if dlg.ShowModal() != wx.ID_OK:
                    return
                dd['color'] = dlg.GetColourData().GetColour().GetAsString(wx.C2S_HTML_SYNTAX)
            mf.matplotlib_panel.update_graph()
        elif part == 'delete':
            wx.CallAfter(self.tp.delete_spectra, [ii])    # not from inside the header being destroyed
        elif part == 'title':
            if ii == mf.current:
                self.expand(None if self.open == ii else ii)
            else:
                mf.set_current(ii)


# --------------------------------------------------------------------------- #
# Tabulated panel: Data | Background | Baseline | SHam | Options | Exp Over.
# --------------------------------------------------------------------------- #
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow):
        super().__init__(parent, size=(280, 400))
        self.parent = mainWindow

        # One-shot simulation button above the tabs, as in pyHYSCORE
        actions = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_doSIM = PillButton(self, "SIM", color_key="pill_accent_2", draw_badge=theme.draw_play_badge,
                                    size=wx.Size(self.FromDIP(90), self.FromDIP(30)))
        self.btn_doSIM.SetToolTip("Simulate every loaded spectrum")
        self.btn_doSIM.Bind(wx.EVT_BUTTON, lambda e: self.parent.runSim())
        actions.Add(self.btn_doSIM, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(4))

        # The Background and Baseline tabs edit the spectrum picked here
        row_cur = wx.BoxSizer(wx.HORIZONTAL)
        row_cur.Add(wx.StaticText(self, label="Spectrum:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        self.choice_current = wx.Choice(self)
        self.choice_current.Bind(wx.EVT_CHOICE, lambda e: self.parent.set_current(self.choice_current.GetSelection()))
        row_cur.Add(self.choice_current, 1, wx.EXPAND)

        self.nb_book = wx.Simplebook(self)

        # ---------------- Data tab
        self.data_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.data_panel, "Data")
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_delete = PillButton(self.data_panel, "Data", color_key="pill_warn", draw_badge=theme.draw_cross_badge)
        self.btn_delete.SetToolTip("Delete loaded spectra.")
        self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete)
        sizerBtns.Add(self.btn_delete, 1, wx.ALL, 3)
        self.btn_header = PillButton(self.data_panel, "Header", color_key="pill_load", draw_badge=theme.draw_textlines_badge)
        self.btn_header.SetToolTip("Show the file header (parameters) of the current spectrum.")
        self.btn_header.Bind(wx.EVT_BUTTON, lambda e: self.show_header(self.parent.Data, self.parent.current))
        sizerBtns.Add(self.btn_header, 1, wx.ALL, 3)
        self.btn_moveUp = wx.Button(self.data_panel, label="▲", size=(24, 22), style=wx.BORDER_NONE)
        self.btn_moveDown = wx.Button(self.data_panel, label="▼", size=(24, 22), style=wx.BORDER_NONE)
        self.btn_moveUp.SetToolTip("Move the current spectrum up the list")
        self.btn_moveDown.SetToolTip("Move the current spectrum down the list")
        self.btn_moveUp.Bind(wx.EVT_BUTTON, lambda e: self.on_move_data(-1))
        self.btn_moveDown.Bind(wx.EVT_BUTTON, lambda e: self.on_move_data(+1))
        sizerBtns.Add(self.btn_moveUp, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)
        sizerBtns.Add(self.btn_moveDown, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 2)
        sizer.Add(sizerBtns, 0, wx.EXPAND | wx.ALL, 3)
        row_show = wx.BoxSizer(wx.HORIZONTAL)
        for label, show in (("Show all", True), ("Hide all", False)):
            btn = PillButton(self.data_panel, label, color_key="pill_calm")
            btn.SetToolTip(f"{label.split()[0]} every spectrum on the plot.")
            btn.Bind(wx.EVT_BUTTON, lambda e, s=show: self.on_show_all(s))
            row_show.Add(btn, 1, wx.ALL, 3)
        sizer.Add(row_show, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 3)

        self.filelist = SpectrumList(self.data_panel, self)
        sizer.Add(self.filelist, 1, wx.EXPAND)

        grid = wx.FlexGridSizer(cols=2, vgap=4, hgap=6)
        grid.AddGrowableCol(1)
        grid.Add(wx.StaticText(self.data_panel, label="Compare:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.seg_norm = SegmentedPill(self.data_panel, ["As is", "Max"])
        self.seg_norm.SetToolTip("As is: data x Scale\nMax: largest absolute value of the corrected spectrum = 1")
        self.seg_norm.Bind(wx.EVT_RADIOBOX, self.on_compare)
        grid.Add(self.seg_norm, 1, wx.EXPAND)
        grid.Add(wx.StaticText(self.data_panel, label="Stack offset:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.spin_offset = wx.SpinCtrlDouble(self.data_panel, min=0, max=5, inc=0.05, initial=0)
        self.spin_offset.SetDigits(2)
        self.spin_offset.SetToolTip("Vertical step between spectra, as a fraction of the tallest trace.\n"
                                    "Right-click to set the arrow step")
        self.spin_offset.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_compare)
        theme.bind_spin_enter(self.spin_offset, self.on_compare)
        for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
            self.spin_offset.Bind(evt, lambda e: self.on_step_menu(self.spin_offset, (0.01, 0.05, 0.1, 0.5, 1), 2))
        grid.Add(self.spin_offset, 1, wx.EXPAND)
        grid.Add(wx.StaticText(self.data_panel, label="X axis zero:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.choice_xnuc = wx.Choice(self.data_panel, choices=["off"] + list(self.parent.Sys.isotopes().keys()))
        self.choice_xnuc.SetSelection(0)
        self.choice_xnuc.SetToolTip("Put the Larmor frequency of this nucleus at 0: every spectrum is shifted by it,\n"
                                    "computed at that spectrum's own B0 (or the Field of 'Exp Over.' when set).\n"
                                    "Simulations and their components move with their spectrum, whatever nucleus they are of.")
        self.choice_xnuc.Bind(wx.EVT_CHOICE, self.on_xnuc)
        grid.Add(self.choice_xnuc, 1, wx.EXPAND)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 6)
        self.data_panel.SetSizer(sizer)

        # ---------------- Background tab
        self.bg_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.bg_panel, "Background")
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(self.bg_panel, label="Loaded backgrounds (file browser '+ BG'):"), 0, wx.ALL, 4)
        self.list_bg = wx.ListBox(self.bg_panel, style=wx.LB_SINGLE | wx.BORDER_NONE)
        sizer.Add(self.list_bg, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 4)
        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_delBG = PillButton(self.bg_panel, "BG", color_key="pill_warn", draw_badge=theme.draw_cross_badge)
        self.btn_delBG.SetToolTip("Unload the selected background (also removes it from every spectrum).")
        self.btn_delBG.Bind(wx.EVT_BUTTON, self.on_delete_bg)
        sizerBtns.Add(self.btn_delBG, 1, wx.ALL, 3)
        self.btn_headerBG = PillButton(self.bg_panel, "Header", color_key="pill_load", draw_badge=theme.draw_textlines_badge)
        self.btn_headerBG.SetToolTip("Show the file header of the selected background.")
        self.btn_headerBG.Bind(wx.EVT_BUTTON, lambda e: self.show_header(self.parent.BGData, self.list_bg.GetSelection()))
        sizerBtns.Add(self.btn_headerBG, 1, wx.ALL, 3)
        sizer.Add(sizerBtns, 0, wx.EXPAND)

        sizer.Add(wx.StaticText(self.bg_panel, label="Subtracted from the current spectrum:"), 0, wx.ALL, 4)
        self.bgtree = wxpg.PropertyGrid(self.bg_panel, style=wxpg.PG_DEFAULT_STYLE | wxpg.PG_SPLITTER_AUTO_CENTER
                                        | wxpg.PG_HIDE_MARGIN | wxpg.PG_TOOLTIPS | wxpg.PG_NO_INTERNAL_BORDER
                                        | wx.BORDER_NONE)
        if not self.bgtree.GetEditorByName("SpinFloat"):
            self.bgtree.RegisterEditor(MypgPanel.SpinCtrlDoubleEditor, "SpinFloat")
        self.bgtree.Bind(wxpg.EVT_PG_CHANGED, self.on_bgtree)
        sizer.Add(self.bgtree, 1, wx.EXPAND)
        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_addBG = PillButton(self.bg_panel, "Add", color_key="pill_load", draw_badge=theme.draw_plus_badge)
        self.btn_addBG.SetToolTip("Subtract another background (the one selected in the list above).")
        self.btn_addBG.Bind(wx.EVT_BUTTON, self.on_add_bg)
        sizerBtns.Add(self.btn_addBG, 1, wx.ALL, 3)
        self.btn_removeBG = PillButton(self.bg_panel, "Remove", color_key="pill_warn", draw_badge=theme.draw_minus_badge)
        self.btn_removeBG.SetToolTip("Stop subtracting the selected (or last) background.")
        self.btn_removeBG.Bind(wx.EVT_BUTTON, self.on_remove_bg)
        sizerBtns.Add(self.btn_removeBG, 1, wx.ALL, 3)
        sizer.Add(sizerBtns, 0, wx.EXPAND)
        self.btn_applyBG = PillButton(self.bg_panel, "Apply to all spectra", color_key="pill_calm")
        self.btn_applyBG.SetToolTip("Copy this background list (with scales and shifts) to every spectrum.")
        self.btn_applyBG.Bind(wx.EVT_BUTTON, lambda e: self.on_apply_all(('bglist',)))
        sizer.Add(self.btn_applyBG, 0, wx.EXPAND | wx.ALL, 3)
        txt = wx.StaticText(self.bg_panel, label=(
            "The background's frequency axis is moved by its shift,\n"
            "then cubic-spline interpolated onto the spectrum's axis;\n"
            "scale x background is subtracted."))
        sizer.Add(txt, 0, wx.ALL, 4)
        self.bg_panel.SetSizer(sizer)

        # ---------------- Baseline tab
        self.bl_panel = wx.ScrolledWindow(self.nb_book)     # scrolls once there are many fit ranges
        self.bl_panel.SetScrollRate(0, 10)
        self.nb_book.AddPage(self.bl_panel, "Baseline")
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.region_card = theme.Card(self.bl_panel, "Polynomial baseline (after the background)")
        box = self.region_card.body
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(self.region_card, label="Polynomial order:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.spin_order = wx.SpinCtrl(self.region_card, min=-1, max=12, initial=0)
        self.spin_order.SetToolTip("−1 = no baseline subtraction; 0 = offset; 1 = line; ...")
        self.spin_order.Bind(wx.EVT_SPINCTRL, lambda e: self.on_order())
        theme.bind_spin_enter(self.spin_order, self.on_order)
        row.Add(self.spin_order, 1, wx.EXPAND)
        box.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        # one 'from … to … ×' row per fit range, rebuilt by update_controls when their number changes
        self.region_grid = wx.FlexGridSizer(cols=5, vgap=3, hgap=4)
        self.region_grid.AddGrowableCol(1)
        self.region_grid.AddGrowableCol(3)
        self.region_spins = []        # (from spin, to spin) per range
        box.Add(self.region_grid, 0, wx.EXPAND | wx.ALL, 3)
        row = wx.BoxSizer(wx.HORIZONTAL)
        for text, action, tip in (("+ Add range", 'add', "Add a fit range after the last one"),
                                  ("Clear all", 'clear', "Delete all fit ranges")):
            btn = wx.Button(self.region_card, label=text, style=wx.BORDER_NONE)
            btn.SetToolTip(tip)
            btn.Bind(wx.EVT_BUTTON, lambda e, a=action: self.on_region(a))
            row.Add(btn, 0, wx.ALIGN_CENTER_VERTICAL)
            row.AddStretchSpacer()
        box.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        sizer.Add(self.region_card, 0, wx.EXPAND | wx.ALL, 5)
        hint = wx.StaticText(self.bl_panel, label=(
            "The polynomial is fitted to the points inside the\n"
            "ranges and subtracted from the whole spectrum.\n\n"
            "On the plot (Zoom/Pan off, 'Baselines' checked):\n"
            "• drag a grey region to move it, or its edge to resize\n"
            "• Shift-drag on empty space to add a region\n"
            "• right-click a region to delete it"))
        sizer.Add(hint, 0, wx.EXPAND | wx.ALL, 6)
        self.btn_applyBL = PillButton(self.bl_panel, "Apply to all spectra", color_key="pill_calm")
        self.btn_applyBL.SetToolTip("Copy the polynomial order and the fit ranges to every spectrum.")
        self.btn_applyBL.Bind(wx.EVT_BUTTON, lambda e: self.on_apply_all(('blorder', 'blregions')))
        sizer.Add(self.btn_applyBL, 0, wx.EXPAND | wx.ALL, 4)
        self.bl_panel.SetSizer(sizer)

        # ---------------- SHam tab
        self.Sys_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Sys_panel, "SHam")
        sizer = wx.BoxSizer(wx.VERTICAL)
        btnszr = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_addNuc = PillButton(self.Sys_panel, "Nuc", color_key="accent", draw_badge=theme.draw_plus_badge)
        self.btn_addNuc.SetToolTip("Add a nucleus to the spin system.")
        self.btn_addNuc.Bind(wx.EVT_BUTTON, self.on_add_Nuc)
        btnszr.Add(self.btn_addNuc, 0, wx.EXPAND)
        self.btn_delNuc = PillButton(self.Sys_panel, "Nuc", color_key="pill_warn", draw_badge=theme.draw_minus_badge)
        self.btn_delNuc.SetToolTip("Delete a nucleus from the spin system.")
        self.btn_delNuc.Bind(wx.EVT_BUTTON, self.on_delete_nuc)
        btnszr.Add(self.btn_delNuc, 0, wx.EXPAND)
        sizer.Add(btnszr, 0, wx.EXPAND)
        self.Sys_param = MypgPanel.PropGridPanel(self.Sys_panel, Prop_Dict={}, onChangeFunc=self.on_Sys, showModFunc=True)
        self.Sys_param.tooltips = self.parent.Sys.getToolTips()
        first = self.getNuc('1H')
        first['useFor'] = ['all', ['sim.only', 'ori.sel.only', 'all']]
        self.Sys_param.SetFromParClean({'spin(1)': {'S': 1/2, 'g': np.array([2, 2, 2], dtype=float)}, 'nuc(1)': first})
        self.parent.Sys.setFromCtrl(self.Sys_param.parameters)
        sizer.Add(self.Sys_param, 1, wx.EXPAND)

        # Simulation components: one row per nucleus of the spin system, rebuilt by update_components
        self.comp_card = theme.Card(self.Sys_panel, "Simulation components")
        self.comp_grid = wx.FlexGridSizer(cols=3, vgap=3, hgap=6)
        self.comp_grid.AddGrowableCol(1)
        self.comp_ctrls = []          # (scale spin, colour picker) per nucleus
        self.comp_card.body.Add(self.comp_grid, 0, wx.EXPAND | wx.ALL, 3)
        sizer.Add(self.comp_card, 0, wx.EXPAND | wx.ALL, 4)
        self.update_components()
        self.Sys_panel.SetSizer(sizer)

        # ---------------- Options tab
        self.Opt_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Opt_panel, "Options")
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.Opt_param = MypgPanel.PropGridPanel(self.Opt_panel, Prop_Dict={}, onChangeFunc=self.on_Opt)
        self.Opt_param.tooltips = self.parent.Opt.getToolTips()
        self.Opt_param.SetFromParClean(self.parent.Opt.getDefaultDict())
        sizer.Add(self.Opt_param, 1, wx.EXPAND)

        # Target grid points: nKnots means a different thing in every Grid method, this count is comparable
        card = theme.Card(self.Opt_panel, "Target grid points")
        self.spin_gridTarget = wx.SpinCtrl(card, min=1, max=10000000, initial=int(self.parent.Opt.nGridPoints() or 1),
                                           style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_gridTarget.SetToolTip("Requested number of orientations. nKnots is integer, so the grid "
                                        "snaps to the nearest achievable size -- the result is shown below.")
        card.body.Add(self.spin_gridTarget, 0, wx.EXPAND | wx.ALL, 2)
        self.txt_gridActual = wx.TextCtrl(card, style=wx.TE_READONLY | wx.TE_CENTRE | wx.BORDER_NONE)
        self.txt_gridActual.SetToolTip("Grid actually produced for the current Grid method and nKnots.")
        card.body.Add(self.txt_gridActual, 0, wx.EXPAND | wx.ALL, 2)
        self.spin_gridTarget.Bind(wx.EVT_SPINCTRL, self.on_gridTarget)
        self.spin_gridTarget.Bind(wx.EVT_TEXT_ENTER, self.on_gridTarget)
        self.updateGridInfo()
        self.btn_doOriSel = PillButton(card, "Update Ori.Sel. Grid", color_key="pill_load", draw_badge=theme.draw_play_badge)
        self.btn_doOriSel.SetToolTip("OriSelType 'precalculated': compute the orientation selection of every\n"
                                     "spectrum by full diagonalisation and keep it with the spectrum.")
        self.btn_doOriSel.Bind(wx.EVT_BUTTON, self.parent.on_update_orisel)
        card.body.Add(self.btn_doOriSel, 0, wx.EXPAND | wx.ALL, 2)
        sizer.Add(card, 0, wx.EXPAND | wx.ALL, 4)

        card = theme.Card(self.Opt_panel, "Simulation scaling")
        self.seg_simfit = SegmentedPill(card, ["Match maxima", "Least squares"])
        self.seg_simfit.SetToolTip("Amplitude of every simulation on its own corrected spectrum.\n"
                                   "Match maxima: the largest absolute value of the simulation is put on the\n"
                                   "   largest absolute value of the spectrum.\n"
                                   "Least squares: the amplitude that minimises the squared difference\n"
                                   "   over the whole spectrum.\n"
                                   "'Sim scale' of each spectrum (Data tab) multiplies on top of either.")
        self.seg_simfit.Bind(wx.EVT_RADIOBOX, self.on_simfit)
        card.body.Add(self.seg_simfit, 0, wx.EXPAND | wx.ALL, 3)
        sizer.Add(card, 0, wx.EXPAND | wx.ALL, 4)
        self.Opt_panel.SetSizer(sizer)

        # ---------------- Exp Over. tab
        self.Exp_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Exp_panel, "Exp Over.")
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.Exp_param = MypgPanel.PropGridPanel(self.Exp_panel, Prop_Dict={}, onChangeFunc=self.on_Exp)
        self.Exp_param.tooltips = self.parent.Exp.getToolTips()
        self.Exp_param.SetFromParClean(self.parent.Exp.getDefaultDict())
        sizer.Add(self.Exp_param, 1, wx.EXPAND)
        self.Exp_panel.SetSizer(sizer)

        self.nb = PillTabBar(self, self.nb_book, ["Data", "Background", "Baseline", "SHam", "Options", "Exp Over."],
                             tooltips=["Loaded spectra", "Background spectra and their subtraction",
                                       "Polynomial baseline of the current spectrum", "Spin Hamiltonian parameters",
                                       "Simulation options", "Global override of experimental parameters"])
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(actions, 0, wx.EXPAND | wx.TOP | wx.BOTTOM, self.FromDIP(6))
        sizer.Add(wx.StaticLine(self, style=wx.LI_HORIZONTAL), 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(4))
        sizer.Add(row_cur, 0, wx.EXPAND | wx.ALL, 4)
        sizer.Add(self.nb, 0, wx.EXPAND)
        sizer.Add(self.nb_book, 1, wx.EXPAND)
        self.SetSizer(sizer)

    # ------------------------------------------------------------------
    def update_filelist(self):
        """Rebuild the spectrum chooser, the Data list and the background list."""
        mf = self.parent
        self.choice_current.SetItems([f"{ii + 1}: {dd['title']}" for ii, dd in enumerate(mf.Data)])
        if mf.current >= 0:
            self.choice_current.SetSelection(mf.current)
        self.filelist.rebuild()
        self.list_bg.Set([bd['title'] for bd in mf.BGData])
        if mf.BGData:
            self.list_bg.SetSelection(0)
        self.update_controls()

    def update_controls(self):
        """Push the current spectrum's processing settings into the Background and Baseline tabs."""
        mf = self.parent
        self.bgtree.Clear()
        if mf.current < 0:
            self.region_grid.Clear(delete_windows=True)
            self.region_spins.clear()
            self.bl_panel.FitInside()
            return
        dd = mf.Data[mf.current]
        titles = [bd['title'] for bd in mf.BGData]
        for jj, b in enumerate(dd['bglist']):
            cat = self.bgtree.Append(wxpg.PropertyCategory(f"Background {jj + 1}", f"cat {jj}"))
            labels = titles if b['bgdata'] in titles else titles + [b['bgdata'] + " (missing)"]
            self.bgtree.AppendIn(cat, wxpg.EnumProperty("Background", f"bgdata {jj}", labels=labels,
                                                        value=labels.index(b['bgdata']) if b['bgdata'] in titles else len(titles)))
            for key, label, step, helpstr in (
                    ('scale', "Scale", 0.01, "Multiplier of this background before subtraction."),
                    ('shift', "RF shift [MHz]", 0.01, "Added to the background's frequency axis.")):
                prop = self.bgtree.AppendIn(cat, wxpg.FloatProperty(label, f"{key} {jj}", value=b[key]))
                self.bgtree.SetPropertyEditor(prop, "SpinFloat")
                prop.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, step)
                prop.SetHelpString(helpstr)
        self.bgtree.ExpandAll()

        self.spin_order.SetValue(dd['blorder'])
        regs = dd['blregions']
        if len(self.region_spins) != len(regs):
            self.region_grid.Clear(delete_windows=True)
            self.region_spins.clear()
            for j in range(len(regs)):
                pair = []
                for text in ("from", "to"):
                    self.region_grid.Add(wx.StaticText(self.region_card, label=text), 0, wx.ALIGN_CENTER_VERTICAL)
                    spin = wx.SpinCtrlDouble(self.region_card, min=-1e6, max=1e6, inc=0.1, initial=0)
                    spin.SetDigits(3)
                    spin.SetToolTip("Fit range in MHz (Enter to apply). Right-click to set the arrow step.")
                    spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, j=j: self.on_region('edit', j))
                    theme.bind_spin_enter(spin, lambda j=j: self.on_region('edit', j))
                    for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                        spin.Bind(evt, lambda e, s=spin: self.on_step_menu(s, (0.001, 0.01, 0.1, 1, 10, 100), 3))
                    self.region_grid.Add(spin, 1, wx.EXPAND)
                    pair.append(spin)
                btn = wx.Button(self.region_card, label="×", size=(24, 22), style=wx.BORDER_NONE)
                btn.SetToolTip("Delete this range")
                btn.Bind(wx.EVT_BUTTON, lambda e, j=j: self.on_region('delete', j))
                self.region_grid.Add(btn, 0, wx.ALIGN_CENTER_VERTICAL)
                self.region_spins.append(tuple(pair))
            theme.apply_theme_to_window(self.region_card, theme.get_theme())
            self.bl_panel.Layout()
            self.bl_panel.FitInside()
        for (s_lo, s_hi), (lo, hi) in zip(self.region_spins, regs):
            s_lo.SetValue(lo)
            s_hi.SetValue(hi)

    def show_header(self, datalist, idx):
        if idx is None or idx < 0 or idx >= len(datalist):
            wx.MessageBox("Nothing selected.", "Header", wx.OK | wx.ICON_INFORMATION)
            return
        if not datalist[idx]['dsc']:
            wx.MessageBox("This spectrum carries no header entries.", "Header", wx.OK | wx.ICON_INFORMATION)
            return
        with DictPopup(self, datalist[idx]['dsc'], title=datalist[idx]['title']) as dlg:
            dlg.ShowModal()

    # ---------------- Data tab handlers
    def on_filetree(self, event):
        key, ii = event.GetPropertyName().split()
        # an EnumProperty's value is the index into its choice list; the spectrum stores the sequence name
        self.parent.Data[int(ii)][key] = event.GetProperty().GetValueAsString() if key == 'seq' else event.GetValue()
        if key in ('freq', 'field', 'seq', 'tau', 'tinv') and self.parent.tb_autoSIM.GetValue():
            self.parent.runSim()       # redraws the figure itself at the end
        else:
            self.parent.matplotlib_panel.update_graph()

    def on_filetree_selected(self, event):
        prop = event.GetProperty()
        if prop is None:
            return
        ii = int(prop.GetName().split()[-1])
        if ii != self.parent.current:
            self.parent.set_current(ii)

    def on_delete(self, event):
        mf = self.parent
        with wx.MultiChoiceDialog(self, "Select the spectra to delete:", "Delete spectra",
                                  [dd['title'] for dd in mf.Data]) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            selected = dlg.GetSelections()
        self.delete_spectra(selected)

    def on_show_all(self, show):
        for dd in self.parent.Data:
            dd['show'] = show
        self.filelist.refresh_headers()
        self.parent.matplotlib_panel.update_graph()

    def delete_spectra(self, indices):
        mf = self.parent
        for ii in sorted(indices, reverse=True):
            mf.Data.pop(ii)
            if mf.current > ii:
                mf.current -= 1
        mf.current = min(mf.current, len(mf.Data) - 1)
        mf.update_everything()

    def on_move_data(self, direction):
        mf = self.parent
        ii, new = mf.current, mf.current + direction
        if ii < 0 or new < 0 or new >= len(mf.Data):
            return
        mf.Data[ii], mf.Data[new] = mf.Data[new], mf.Data[ii]
        mf.current = new
        mf.update_everything()

    def on_compare(self, event=None):
        mp = self.parent.matplotlib_panel
        mp.norm = ('none', 'max')[self.seg_norm.GetSelection()]
        offset = self.spin_offset.GetValue()
        if event is None and offset == mp.offset:   # Enter on an unchanged value
            return
        mp.offset = offset
        mp.update_graph()

    def on_xnuc(self, event):
        mp = self.parent.matplotlib_panel
        mp.xnuc = self.choice_xnuc.GetStringSelection() if self.choice_xnuc.GetSelection() > 0 else ''
        mp.keepx = False             # a zoom made on the other axis means nothing on this one
        mp.update_graph()

    def on_apply_all(self, keys):
        mf = self.parent
        if mf.current < 0:
            return
        for dd in mf.Data:
            for key in keys:
                dd[key] = copy.deepcopy(mf.Data[mf.current][key])
        mf.matplotlib_panel.update_graph()

    # ---------------- Background tab handlers
    def on_bgtree(self, event):
        mf = self.parent
        if mf.current < 0:
            return
        key, jj = event.GetPropertyName().split()
        b = mf.Data[mf.current]['bglist'][int(jj)]
        b[key] = event.GetProperty().GetValueAsString() if key == 'bgdata' else float(event.GetValue())
        mf.matplotlib_panel.update_graph()

    def on_add_bg(self, event):
        mf = self.parent
        if mf.current < 0:
            return
        if not mf.BGData:
            wx.MessageBox("Load a background first: select its file in the browser and press '+ BG'.",
                          "No background loaded", wx.OK | wx.ICON_INFORMATION)
            return
        sel = max(self.list_bg.GetSelection(), 0)
        mf.Data[mf.current]['bglist'].append({'bgdata': mf.BGData[sel]['title'], 'scale': 1.0, 'shift': 0.0})
        self.update_controls()
        mf.matplotlib_panel.update_graph()

    def on_remove_bg(self, event):
        mf = self.parent
        if mf.current < 0 or not mf.Data[mf.current]['bglist']:
            return
        prop = self.bgtree.GetSelection()
        jj = int(prop.GetName().split()[-1]) if prop is not None else -1
        mf.Data[mf.current]['bglist'].pop(jj)
        self.update_controls()
        mf.matplotlib_panel.update_graph()

    def on_delete_bg(self, event):
        mf = self.parent
        sel = self.list_bg.GetSelection()
        if sel < 0:
            return
        title = mf.BGData.pop(sel)['title']
        for dd in mf.Data:
            dd['bglist'] = [b for b in dd['bglist'] if b['bgdata'] != title]
        mf.update_everything()

    # ---------------- Baseline tab handlers
    def on_order(self):
        mf = self.parent
        if mf.current < 0 or self.spin_order.GetValue() == mf.Data[mf.current]['blorder']:
            return
        mf.Data[mf.current]['blorder'] = self.spin_order.GetValue()
        mf.matplotlib_panel.update_graph()

    def on_region(self, action, j=None):
        """Fit ranges of the baseline: 'edit' range j from its spins, 'delete' range j,
        'add' one after the last (5 % of the frequency span; mid-spectrum if there is none), 'clear' all."""
        mf = self.parent
        if mf.current < 0:
            return
        dd = mf.Data[mf.current]
        regs = dd['blregions']
        if action == 'edit':
            if j >= len(regs):              # the rows are about to be rebuilt
                return
            s_lo, s_hi = self.region_spins[j]
            new = sorted((s_lo.GetValue(), s_hi.GetValue()))
            if new == regs[j]:              # Enter or focus change on an unchanged value
                return
            regs[j] = new
            s_lo.SetValue(new[0])
            s_hi.SetValue(new[1])
        else:
            if action == 'delete':
                del regs[j]
            elif action == 'clear':
                regs.clear()
            else:
                x0, x1 = float(np.min(dd['ax']['x'])), float(np.max(dd['ax']['x']))
                w = 0.05*(x1 - x0)
                lo = max(hi for _, hi in regs) + 0.01*(x1 - x0) if regs else 0.5*(x0 + x1 - w)
                lo = min(lo, x1 - w)
                regs.append([lo, lo + w])
            wx.CallAfter(self.update_controls)    # rebuilds the rows, the clicked button among them
        mf.matplotlib_panel.update_graph()

    def on_step_menu(self, spin, steps=(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 1e4, 1e5, 1e6), digits=4):
        """Right-click on a spin box: choose its arrow step (as on the property-grid spins)."""
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
        v = spin.GetValue()
        spin.SetDigits(max(digits, int(np.ceil(-np.log10(ids[choice])))))
        spin.SetIncrement(ids[choice])
        spin.SetValue(v)

    # ---------------- SHam tab handlers
    def getNuc(self, nucname):
        """Property-grid block of one nucleus, with the isotope table as the choice list."""
        isotopes = list(self.parent.Sys.isotopes().keys())
        return {'Nucs': [nucname, isotopes],
                'nNucs': 1,
                'A': np.array([1, 1, 2], dtype=float),
                'Apa': np.array([0, 0, 0], dtype=float),
                'Q_K': 0.0,
                'Q_eta': 0.0,
                'Qpa': np.array([0, 0, 0], dtype=float),
                'lwEndor': 0.1,
                'useFor': ['sim.only', ['sim.only', 'ori.sel.only', 'all']]}

    def on_add_Nuc(self, event):
        isotopes = self.parent.Sys.isotopes()
        names = list(isotopes.keys())
        with wx.SingleChoiceDialog(self, "Select isotope to include:", "Choose Nucleus",
                                   [f'{nn:>8}\tI={isotopes[nn]["I"]:>4}' for nn in names]) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            selected = dlg.GetSelection()
        Params = self.Sys_param.parameters
        nNucs = sum(1 for key in Params if 'nuc' in key)
        Params[f'nuc({nNucs + 1})'] = self.getNuc(names[selected])
        self.parent.Sys.setFromCtrl(Params)
        self.Sys_param.SetFromParClean(Params)
        self.update_components()
        if self.parent.tb_autoSIM.GetValue():
            self.parent.runSim()

    def on_delete_nuc(self, event):
        Params = self.Sys_param.parameters
        keys = [key for key in Params if 'nuc' in key]
        if len(keys) < 2:
            wx.MessageBox("The spin system needs at least one nucleus.", "Delete nucleus", wx.OK | wx.ICON_INFORMATION)
            return
        with wx.SingleChoiceDialog(self, "Select the nucleus to delete:", "Delete nucleus",
                                   [f"{key}:  {Params[key]['Nucs'][0]}" for key in keys]) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            sel = dlg.GetSelection()
        del keys[sel]
        del self.parent.components[sel]        # scale and colour stay with the nuclei that remain
        # the remaining nuclei are renumbered 1, 2, ... in their old order
        new = {key: val for key, val in Params.items() if 'nuc' not in key}
        for kk, key in enumerate(keys):
            new[f'nuc({kk + 1})'] = Params[key]
        self.parent.Sys.setFromCtrl(new)
        self.Sys_param.SetFromParClean(new)
        self.update_components()
        for dd in self.parent.Data:            # the stored components are numbered for the old spin system
            dd['simdata'] = dd['simax'] = dd['simparts'] = None
        if self.parent.tb_autoSIM.GetValue():
            self.parent.runSim()
        else:
            self.parent.matplotlib_panel.update_graph()

    def update_components(self):
        """Rebuild the 'Simulation components' rows: a scale box and a colour button per nucleus of the spin system."""
        mf = self.parent
        names = [val['Nucs'][0] for key, val in self.Sys_param.parameters.items() if 'nuc' in key]
        while len(mf.components) < len(names):
            mf.components.append({'scale': 1.0, 'color': COMPONENT_COLORS[len(mf.components) % len(COMPONENT_COLORS)]})
        del mf.components[len(names):]
        self.comp_grid.Clear(delete_windows=True)
        self.comp_ctrls = []
        for text in ("nucleus", "scale", "colour"):
            self.comp_grid.Add(wx.StaticText(self.comp_card, label=text), 0, wx.ALIGN_CENTER_VERTICAL)
        for k, name in enumerate(names):
            self.comp_grid.Add(wx.StaticText(self.comp_card, label=f"{name} ({k + 1})"), 0, wx.ALIGN_CENTER_VERTICAL)
            spin = wx.SpinCtrlDouble(self.comp_card, min=-1e6, max=1e6, inc=0.1, initial=mf.components[k]['scale'])
            spin.SetDigits(3)
            spin.SetToolTip("Multiplier of this nucleus's component in the summed simulation, in every spectrum.\n"
                            "Right-click to set the arrow step")
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=k: self.on_component(k))
            theme.bind_spin_enter(spin, lambda k=k: self.on_component(k))
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                spin.Bind(evt, lambda e, s=spin: self.on_step_menu(s, (0.001, 0.01, 0.1, 1, 10), 3))
            self.comp_grid.Add(spin, 1, wx.EXPAND)
            pick = wx.ColourPickerCtrl(self.comp_card, colour=wx.Colour(mf.components[k]['color']))
            pick.SetToolTip("Colour of this nucleus's component trace")
            pick.Bind(wx.EVT_COLOURPICKER_CHANGED, lambda e, k=k: self.on_component(k))
            self.comp_grid.Add(pick, 0, wx.ALIGN_CENTER_VERTICAL)
            self.comp_ctrls.append((spin, pick))
        theme.apply_theme_to_window(self.comp_card, theme.get_theme())
        self.Sys_panel.Layout()

    def on_component(self, k):
        """Scale or colour of nucleus k's component changed: the stored components are re-summed, nothing is re-simulated."""
        spin, pick = self.comp_ctrls[k]
        new = {'scale': spin.GetValue(), 'color': pick.GetColour().GetAsString(wx.C2S_HTML_SYNTAX)}
        if new == self.parent.components[k]:       # Enter on an unchanged value
            return
        self.parent.components[k] = new
        self.parent.matplotlib_panel.update_graph()

    def on_Sys(self, mainname, name, val):
        self.parent.Sys.setFromCtrl(self.Sys_param.parameters)
        if name == 'Nucs':
            self.update_components()       # the isotope of a nucleus changed: its row is named after it
        if self.parent.tb_autoSIM.GetValue():
            self.parent.runSim()

    # ---------------- Options and Exp Over. tab handlers
    def on_Exp(self, mainname, name, val):
        self.parent.Exp.setFromCtrl(self.Exp_param.parameters)
        if name == 'Field' and self.parent.matplotlib_panel.xnuc:
            self.parent.matplotlib_panel.update_graph()      # the Larmor shift of the x axis follows the field

    def on_Opt(self, mainname, name, val):
        self.parent.Opt.setFromCtrl(self.Opt_param.parameters)
        if name == 'OriSelType':
            self.parent.on_update_orisel(None)
        if name in ('Grid', 'nKnots'):
            self.updateGridInfo()

    def on_simfit(self, event):
        mp = self.parent.matplotlib_panel
        mp.simfit = ('max', 'lsq')[self.seg_simfit.GetSelection()]
        mp.update_graph()

    def on_gridTarget(self, event):
        """Convert a requested grid size into nKnots and push it to the grid (trigger_run=False: no re-entry into on_Opt)."""
        if event is not None:
            event.Skip()
        Opt = self.parent.Opt
        k = Opt.nKnotsForGridPoints(self.spin_gridTarget.GetValue())
        if k is None:
            self.updateGridInfo()
            return
        prop = self.Opt_param.pg.GetPropertyByName('nKnots')
        if prop is not None and prop.GetValue() != k:
            prop.SetValue(k)
            self.Opt_param.OnValueChanged(None, prop=prop, trigger_run=False)
            Opt.setFromCtrl(self.Opt_param.parameters)
        self.updateGridInfo(keepTarget=True)

    def updateGridInfo(self, keepTarget=False):
        """Keep nKnots, Target grid points and the actual count consistent; keepTarget leaves the spin as typed."""
        Opt = self.parent.Opt
        n = Opt.nGridPoints()
        gname = Opt.Grid[0] if isinstance(Opt.Grid, (list, tuple)) else Opt.Grid
        txt = Opt.getToolTips('nKnots')
        if n is not None:
            txt = f"{txt}\n\nCurrently: '{gname}' with nKnots={Opt.nKnots}  ->  {n} grid points"
        self.Opt_param.tooltips['nKnots'] = txt
        prop = self.Opt_param.pg.GetPropertyByName('nKnots')
        if prop is not None:
            prop.SetHelpString(txt)
        self.txt_gridActual.SetValue('grid: n/a' if n is None else f'actual: {n} points  (nKnots = {Opt.nKnots})')
        if not keepTarget and n is not None and self.spin_gridTarget.GetValue() != n:
            self.spin_gridTarget.SetValue(int(n))     # nKnots edited directly, or a session loaded -> follow it


# --------------------------------------------------------------------------- #
# EPR spectrum window
# --------------------------------------------------------------------------- #
class EPRFrame(wx.Frame):
    """Separate window: the EPR spectrum as the orientation selection sees it -- at every field the sum over the
    orientation grid of the selection weights -- with the field position of every shown ENDOR spectrum marked on it.
    Its settings live in MainFrame.epr, so they survive closing the window and are saved in sessions."""

    def __init__(self, parent, mainWindow):
        super().__init__(parent, title="EPR spectrum and field positions", size=(980, 620))
        self.mf = mainWindow
        self.key = None                   # everything the stored spectrum was computed from; None = no spectrum
        self.fields = np.array([])        # field axis [mT] of the stored spectrum
        self.spectrum = np.array([])
        self.fault = ''                   # why there is no spectrum, written into the empty plot

        splitter = wx.SplitterWindow(self, style=wx.SP_3D | wx.SP_LIVE_UPDATE)
        splitter.SetMinimumPaneSize(140)
        left = scrolled.ScrolledPanel(splitter, style=wx.BORDER_NONE)
        right = wx.Panel(splitter)
        cards = wx.BoxSizer(wx.VERTICAL)

        card = theme.Card(left, "Microwave frequency")
        self.choice_source = wx.Choice(card)
        self.choice_source.SetToolTip("The spectrum whose microwave frequency the EPR spectrum is computed at.\n"
                                      "Every other spectrum is marked twice: dashed at the field it was recorded at,\n"
                                      "solid at that field times f(selected)/f(its own), where the same g-value sits\n"
                                      "at the selected frequency.")
        self.choice_source.Bind(wx.EVT_CHOICE, self.on_change)
        card.body.Add(self.choice_source, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        card = theme.Card(left, "X axis")
        self.seg_axis = SegmentedPill(card, ["Field", "g"])
        self.seg_axis.SetToolTip("Field in mT, or g = h·f/(μB·B) at the selected microwave frequency.")
        self.seg_axis.Bind(wx.EVT_RADIOBOX, self.on_change)
        card.body.Add(self.seg_axis, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        card = theme.Card(left, "Orientation grid")
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(wx.StaticText(card, label="nKnots:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
        self.spin_knots = wx.SpinCtrl(card, min=2, max=5000, initial=max(2, int(self.mf.epr['nKnots'])),
                                      style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_knots.SetToolTip("Grid resolution of this window only; Grid, OriSelType and everything else\n"
                                   "are those of the Options, SHam and 'Exp Over.' tabs.")
        self.spin_knots.Bind(wx.EVT_SPINCTRL, self.on_change)
        theme.bind_spin_enter(self.spin_knots, self.on_change)
        row.Add(self.spin_knots, 1, wx.EXPAND)
        card.body.Add(row, 0, wx.EXPAND | wx.ALL, 3)
        self.lbl_grid = wx.StaticText(card, label="")
        card.body.Add(self.lbl_grid, 0, wx.EXPAND | wx.ALL, 3)
        cards.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        note = wx.StaticText(left, label=(
            "The curve is the sum over the orientation\n"
            "grid of the weights the ENDOR simulation\n"
            "uses, at every field: its height at a marker\n"
            "is the share of orientations excited there."))
        cards.Add(note, 0, wx.EXPAND | wx.ALL, 8)
        left.SetSizer(cards)
        left.SetupScrolling(scroll_x=True, rate_y=10)    # a narrow pane scrolls instead of clipping a card

        self.figure = Figure(figsize=(2, 2), dpi=100, layout='constrained')
        self.canvas = FigureCanvasWxAgg(right, -1, self.figure)
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.canvas, 1, wx.EXPAND)
        box.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        right.SetSizer(box)
        self.toolbar.update()

        self.statusbar = theme.StatusStrip(self, widths=[-1, 300])
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('figure_leave_event', lambda e: self.statusbar.SetStatusText("", 1))

        splitter.SplitVertically(left, right, sashPosition=290)
        splitter.SetSashGravity(0)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(splitter, 1, wx.EXPAND)
        sizer.Add(self.statusbar, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self.redraw_debouncer.on_size)
        self.Bind(wx.EVT_CLOSE, lambda e: (setattr(self.mf, 'epr_window', None), self.Destroy()))
        self.apply_theme()

    def apply_theme(self):
        theme.apply_theme_to_window(self, theme.get_theme())
        theme.style_figure(self.figure, self.figure.axes)
        self.canvas.draw_idle()
        self.Refresh()

    def positions(self):
        """(microwave frequency [GHz], field [mT]) of every loaded spectrum, as its simulation takes them:
        the value of 'Exp Over.' when that is set, else the spectrum's own."""
        E = self.mf.tabulated_panel.Exp_param.parameters
        return [(E['mwFreq'] if E['mwFreq'] >= 0 else dd['freq'], E['Field'] if E['Field'] >= 0 else dd['field'])
                for dd in self.mf.Data]

    def refresh(self):
        """After every redraw of the main plot: the controls follow MainFrame.epr and the loaded spectra,
        the spectrum is recomputed if anything it depends on changed, and the markers are redrawn."""
        mf = self.mf
        labels = [f"{ii + 1}: {dd['title']}   {freq:.6g} GHz" for ii, (dd, (freq, field)) in enumerate(zip(mf.Data, self.positions()))]
        if labels != self.choice_source.GetStrings():
            self.choice_source.SetItems(labels)
        mf.epr['source'] = min(max(mf.epr['source'], 0), max(len(labels) - 1, 0))
        if labels:
            self.choice_source.SetSelection(mf.epr['source'])
        self.seg_axis.SetSelection(('field', 'g').index(mf.epr['xaxis']))
        if self.spin_knots.GetValue() != mf.epr['nKnots']:
            self.spin_knots.SetValue(int(mf.epr['nKnots']))
        self.compute()
        self.update_graph()

    def on_change(self, event=None):
        """Any control of this window changed."""
        mf = self.mf
        mf.epr['source'] = max(self.choice_source.GetSelection(), 0)
        mf.epr['xaxis'] = ('field', 'g')[self.seg_axis.GetSelection()]
        mf.epr['nKnots'] = self.spin_knots.GetValue()
        self.compute()
        self.update_graph()

    def compute(self):
        """Recompute the stored spectrum, unless everything it depends on is as it was when it was made."""
        mf = self.mf
        tp = mf.tabulated_panel
        pos = self.positions()
        self.fault = ''
        freq = pos[mf.epr['source']][0] if pos else 0.0
        width = tp.Exp_param.parameters['ExciteWidth']
        if not pos or not freq > 0 or not width > 0:
            self.fault = ("Load a spectrum first" if not pos else
                          "The selected spectrum has no microwave frequency: set 'MW Freq' in the Data tab" if not freq > 0 else
                          "ExciteWidth is 0: no orientation is selected, so there is no excitation profile to show")
            self.key = None
            return
        mf.Sys.setFromCtrl(tp.Sys_param.parameters)
        es = ENDORsim(Sys=mf.Sys, errorFunc=mf.RaiseError)      # with its own Exp and Opt: the main ones stay as they are
        es.verbose = False
        Opt = tp.Opt_param.parameters.copy()
        Opt['nKnots'] = int(mf.epr['nKnots'])
        es.Opt.setFromCtrl(Opt)
        es.Exp.set(mwFreq=freq, ExciteWidth=width, Field=1.0)
        kind = es.Opt.OriSelType[0]
        self.lbl_grid.SetLabel(f"{es.Opt.nGridPoints()} orientations\n'{es.grid_name()}' grid, '{kind}'")

        # field range: the powder pattern of g, widened by the hyperfine couplings that take part and by the
        # excitation width, and by any marker that lies near it; sampled at a sixth of the excitation width
        g = np.asarray(mf.Sys.g[0], dtype=float)
        mT = es.h_ov_beta_mTMHz/g.min()                  # mT per MHz, at the smallest g
        hfc = 0.0 if kind == 'g_eff' else sum(
            float(np.max(np.abs(A)))*I*max(1, int(n)) for A, I, n, use in zip(mf.Sys.A, mf.Sys.I, mf.Sys.nNucs, mf.Sys.useFor)
            if kind == 'g_eff+HFC' or use != 'sim.only')
        lo = G_FACTOR*freq/g.max() - (hfc + 3*width)*mT
        hi = G_FACTOR*freq/g.min() + (hfc + 3*width)*mT
        marks = [b for f, field in pos if f > 0 and field > 0 for b in (field, field*freq/f) if 0.5*lo < b < 2*hi]
        lo, hi = min([lo] + [b - width*mT for b in marks]), max([hi] + [b + width*mT for b in marks])
        n = int(np.clip((hi - lo)/(width*mT/6), 100, 1500))
        # compute_weights reads every nucleus for 'g_eff+HFC', which grid_signature does not cover
        key = (es.grid_signature(), lo, hi, n, tuple(np.asarray(A).tobytes() for A in mf.Sys.A) if kind == 'g_eff+HFC' else ())
        if key == self.key:
            return
        wait = wx.BusyCursor()
        try:
            self.fields = np.linspace(lo, hi, n)
            self.spectrum = es.computeEPR(self.fields)
            self.key = key
        except Exception as e:
            traceback.print_exc()
            print(f"[pyENDOR] The EPR spectrum could not be computed ({type(e).__name__}: {e}); the window shows none")
            self.fault = f"The EPR spectrum could not be computed:\n{type(e).__name__}: {e}"
            self.key = None
        finally:
            del wait

    def update_graph(self):
        mf = self.mf
        pos = self.positions()
        gmode = mf.epr['xaxis'] == 'g'
        markcol = theme.theme_colors()["ax_fg"]
        self.figure.clear()
        ax = self.figure.add_subplot(1, 1, 1)
        missing = []
        if self.key is not None and pos:
            fsel = pos[mf.epr['source']][0]
            toplot = lambda B: G_FACTOR*fsel/B if gmode else B        # g read at the selected frequency
            ax.plot(toplot(self.fields), self.spectrum/(np.max(self.spectrum) or 1.0), color=markcol, lw=1.3,
                    label=f"EPR at {fsel:.6g} GHz")
            for ii, (dd, (freq, field)) in enumerate(zip(mf.Data, pos)):
                if not dd['show']:
                    continue
                if not freq > 0 or not field > 0:
                    missing.append(str(ii + 1))
                    continue
                col, lw = mf.trace_color(dd), 1.8 if ii == mf.current else 1.0
                ax.axvline(toplot(field), color=col, lw=lw, ls='--')
                ax.axvline(toplot(field*fsel/freq), color=col, lw=lw, label=f"{ii + 1}: {dd['title']}")
            ax.plot([], [], color=markcol, lw=1.0, ls='--', label="dashed: field as recorded")
            ax.plot([], [], color=markcol, lw=1.0, label="solid: moved to the selected frequency")
            ax.set_xlim(toplot(self.fields[0]), toplot(self.fields[-1]))      # on the g axis this runs high -> low
            ax.set_ylim(0, 1.08)
            ax.legend(fontsize='small', frameon=False, loc='best', labelcolor=markcol)
        else:
            ax.text(0.5, 0.5, self.fault, ha='center', va='center', transform=ax.transAxes, color=markcol)
        ax.set_xlabel("g" if gmode else "B, mT")
        ax.set_ylabel("EPR: sum of the orientation weights")
        ax.grid(mf.matplotlib_panel.showGrid)
        theme.style_figure(self.figure, [ax])
        self.toolbar.update()
        self.canvas.draw_idle()
        self.statusbar.SetStatusText(f"no marker for spectrum {', '.join(missing)}: set its MW Freq and B0 in the Data tab"
                                     if missing else "Ready", 0)

    def on_mouse_move(self, event):
        """Cursor readout: field and g at the selected microwave frequency, whichever of them the axis shows."""
        pos = self.positions()
        if event.inaxes is None or event.xdata is None or not event.xdata > 0 or self.key is None or not pos:
            self.statusbar.SetStatusText("", 1)
            return
        other = G_FACTOR*pos[self.mf.epr['source']][0]/event.xdata
        B, g = (other, event.xdata) if self.mf.epr['xaxis'] == 'g' else (event.xdata, other)
        self.statusbar.SetStatusText(f"B = {B:.3f} mT,  g = {g:.5f},  y = {event.ydata:.3f}", 1)


# --------------------------------------------------------------------------- #
# Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Data = []
        self.BGData = []
        self.current = -1        # index in self.Data edited by the Background and Baseline tabs
        self.components = []     # per nucleus of the spin system: {'scale' in the summed simulation, 'color' of its trace}
        self.epr_window = None   # EPRFrame while it is open; refreshed by update_graph
        # its settings: index of the spectrum whose MW frequency is used, 'field' | 'g' axis, its own nKnots (0 = not opened yet)
        self.epr = {'source': 0, 'xaxis': 'field', 'nKnots': 0}
        self.Sys = sysPar()
        self.Exp = expENDOR()
        self.Opt = optENDOR()
        self.SetTitle(f"pyENDOR {VERSION}")
        self.SetSize((1280, 740))
        self.theme_name = "light"

        # ---- Top icon toolbar
        self.top_toolbar = wx.Panel(self)
        toolbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.tb_change_folder = theme.IconButton(self.top_toolbar, theme.draw_folder_icon, tooltip="Change folder…")
        self.tb_up = theme.IconButton(self.top_toolbar, theme.draw_up_icon, tooltip="Up")
        self.tb_save_session = theme.IconButton(self.top_toolbar, theme.draw_save_icon, tooltip="Save Session")
        self.tb_autoSIM = theme.IconToggleButton(self.top_toolbar, theme.draw_auto_toggle_icon("SIM"),
                                                 tooltip="auto SIM (re-simulate automatically on change)", value=True)
        self.tb_epr = theme.IconButton(self.top_toolbar, theme.draw_text_icon("EPR"),
                                       tooltip="EPR spectrum…\nThe EPR spectrum the orientation selection is computed from, "
                                               "with the field\nposition of every shown ENDOR spectrum marked on it")
        self.tb_grid = theme.IconToggleButton(self.top_toolbar, theme.draw_grid_icon, tooltip="Axes grid")
        self.btn_theme_icon = theme.IconButton(self.top_toolbar, theme.draw_theme_icon, tooltip="Toggle dark mode (Ctrl+D)")
        for w in (self.tb_change_folder, self.tb_up, None, self.tb_save_session, None, self.tb_autoSIM, None, self.tb_epr,
                  'stretch', None, self.tb_grid, self.btn_theme_icon):
            if w is None:
                toolbar_sizer.Add(wx.StaticLine(self.top_toolbar, style=wx.LI_VERTICAL), 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
            elif isinstance(w, str):
                toolbar_sizer.AddStretchSpacer(1)
            else:
                toolbar_sizer.Add(w, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        self.top_toolbar.SetSizer(toolbar_sizer)
        self.tb_grid.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_grid)
        self.tb_epr.Bind(wx.EVT_BUTTON, self.on_epr_window)
        self.btn_theme_icon.Bind(wx.EVT_BUTTON, self.on_toggle_theme)
        theme_id = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_toggle_theme, id=theme_id)
        self.SetAcceleratorTable(wx.AcceleratorTable([(wx.ACCEL_CMD, ord('D'), theme_id)]))

        # ---- Status strip: message | cursor readout | simulation progress gauge (hidden unless a run is going)
        self.statusbar = theme.StatusStrip(self, widths=[-1, 280])
        self.statusbar.SetStatusText("Ready", 0)
        self.sim_gauge = wx.Gauge(self.statusbar, range=100, style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
        self.statusbar.AddWidget(self.sim_gauge, width=170)
        self.sim_gauge.Hide()

        # ---- file browser | (plot | tabs), as in pyHYSCORE and pyCWEPR
        self.splitter_main = wx.SplitterWindow(self, style=wx.SP_3D)
        self.splitter_main.SetMinimumPaneSize(20)
        self.file_browser = FileBrowserPanel(self.splitter_main, self)
        self.right_panel = wx.Panel(self.splitter_main)
        self.splitter_right = wx.SplitterWindow(self.right_panel, style=wx.SP_3D)
        self.matplotlib_panel = MatplotlibPanel(self.splitter_right, self)
        self.tabulated_panel = TabulatedPanel(self.splitter_right, self)
        self.splitter_right.SplitVertically(self.matplotlib_panel, self.tabulated_panel)
        self.splitter_right.SetSashGravity(1)
        self.splitter_right.SetSashPosition(-430)    # room for six tabs
        self.right_sash = -430
        right_sizer = wx.BoxSizer(wx.HORIZONTAL)
        right_sizer.Add(self.splitter_right, 1, wx.EXPAND)
        self.right_panel.SetSizer(right_sizer)
        self.left_sash = 220
        self.splitter_main.SplitVertically(self.file_browser, self.right_panel, sashPosition=self.left_sash)
        self.splitter_main.SetSashGravity(0)

        # Collapse toggles live outside the splitters so they stay visible when their pane is hidden
        self.toggle_btn = theme.IconButton(self, theme.draw_chevron_icon("left"), tooltip="Hide file browser",
                                           size=wx.Size(18, -1), fill_key="accent")
        self.toggle_btn.Bind(wx.EVT_BUTTON, self.on_toggle_left)
        self.toggle_right_btn = theme.IconButton(self, theme.draw_chevron_icon("right"), tooltip="Hide tabulated panel",
                                                 size=wx.Size(18, -1), fill_key="accent")
        self.toggle_right_btn.Bind(wx.EVT_BUTTON, self.on_toggle_right)
        content_sizer = wx.BoxSizer(wx.HORIZONTAL)
        content_sizer.Add(self.toggle_btn, 0, wx.EXPAND)
        content_sizer.Add(self.splitter_main, 1, wx.EXPAND)
        content_sizer.Add(self.toggle_right_btn, 0, wx.EXPAND)

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

    # -----------------------------------------------------------------------
    def on_toggle_theme(self, event):
        self.theme_name = "dark" if self.theme_name == "light" else "light"
        self.apply_theme()

    def on_epr_window(self, event):
        """Open the EPR spectrum window, or raise it if it is already there."""
        if not self.epr['nKnots']:
            self.epr['nKnots'] = int(self.Opt.nKnots)      # first opening: the grid of the Options tab, then its own
        if self.epr_window is None:
            self.epr_window = EPRFrame(self, self)
            self.epr_window.Show()
        else:
            self.epr_window.Raise()
        self.epr_window.refresh()

    def on_toggle_left(self, event):
        if self.splitter_main.IsSplit():
            self.left_sash = self.splitter_main.GetSashPosition()
            self.splitter_main.Unsplit(self.file_browser)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("right"))
            self.toggle_btn.SetToolTip("Show file browser")
        else:
            self.splitter_main.SplitVertically(self.file_browser, self.right_panel, sashPosition=self.left_sash)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("left"))
            self.toggle_btn.SetToolTip("Hide file browser")

    def on_toggle_right(self, event):
        if self.splitter_right.IsSplit():
            self.right_sash = self.splitter_right.GetSashPosition()
            self.splitter_right.Unsplit(self.tabulated_panel)
            self.toggle_right_btn.SetDrawIcon(theme.draw_chevron_icon("left"))
            self.toggle_right_btn.SetToolTip("Show tabulated panel")
        else:
            self.splitter_right.SplitVertically(self.matplotlib_panel, self.tabulated_panel, sashPosition=self.right_sash)
            self.toggle_right_btn.SetDrawIcon(theme.draw_chevron_icon("right"))
            self.toggle_right_btn.SetToolTip("Hide tabulated panel")

    def on_toggle_grid(self, event):
        mp = self.matplotlib_panel
        mp.showGrid = self.tb_grid.GetValue()
        for ax in mp.figure.axes:
            ax.grid(mp.showGrid)
        mp.canvas.draw_idle()
        if self.epr_window is not None:
            self.epr_window.update_graph()

    def apply_theme(self):
        theme.set_theme(self.theme_name)
        theme.apply_theme_to_window(self, self.theme_name)
        tp = self.tabulated_panel
        tp.filelist.apply_theme()
        theme.theme_propgrid(tp.bgtree, caption_key="pill_load")
        for panel in (tp.Sys_param, tp.Opt_param, tp.Exp_param):
            theme.theme_propgrid(panel.pg, caption_key="accent")
        if self.epr_window is not None:
            self.epr_window.apply_theme()
        self.matplotlib_panel.update_graph()      # simulations and baselines are drawn in the theme's axis colour
        self.Refresh()

    def set_current(self, ii):
        self.current = ii
        self.tabulated_panel.choice_current.SetSelection(ii)
        self.tabulated_panel.filelist.set_current(ii)
        self.tabulated_panel.update_controls()
        self.matplotlib_panel.update_graph()

    def update_everything(self):
        """After spectra were added/removed/reordered: tabs, plot (zoom reset) and browser markers."""
        self.tabulated_panel.update_filelist()
        self.matplotlib_panel.keepx = False
        self.matplotlib_panel.update_graph()
        fb = self.file_browser
        selected, idx = [], fb.list.GetFirstSelected()
        while idx != -1:
            selected.append(fb.items[fb.list.GetItemData(idx)][0])
            idx = fb.list.GetNextSelected(idx)
        fb.refresh_file_list(select_name=selected)

    def RaiseError(self, message):
        """Error handler given to the simulation: show the message, then abort the run that hit it."""
        wx.MessageBox(message, "Error", wx.OK | wx.ICON_ERROR)
        raise Exception(message)

    # -----------------------------------------------------------------------
    # Data handling
    # -----------------------------------------------------------------------
    def loadData(self, path, BG=False):
        ext = os.path.splitext(path)[1].lower()
        try:
            if ext == '.xml':
                self.loadXML(path, BG)
            elif ext in ('.dsc', '.dta'):
                ax, data, dsc = Mybr.brukerread(path, return_ax=True, return_dsc=True)
                self.addDataset(self.newDataset(ax, data, dsc, path), BG)
            else:
                raise ValueError(f"unsupported file type '{ext}'")
        except Exception as e:
            traceback.print_exc()
            wx.MessageBox(f"Loading '{os.path.basename(path)}' failed:\n\n{type(e).__name__}: {e}",
                          "Cannot load", wx.OK | wx.ICON_ERROR)

    def newDataset(self, ax, data, dsc, path):
        """Dataset dict with default processing. Radio-frequency axis always ascending, in MHz; freq in GHz; field in mT."""
        dsc = dsc if isinstance(dsc, dict) else {}
        x = np.asarray(ax['x'], dtype=float).ravel()
        y = np.real(np.squeeze(np.asarray(data))).astype(float)
        if y.ndim != 1 or y.size != x.size:
            raise ValueError(f"not a 1D spectrum (data shape {np.shape(data)}, {x.size} frequency points)")
        units = str(ax.get('xlabel', '')).replace(',', ' ').split()
        x = x*(1e-3 if 'kHz' in units else 1e-6 if 'Hz' in units else 1e3 if 'GHz' in units else 1.0)
        if x[0] > x[-1]:
            x, y = x[::-1].copy(), y[::-1].copy()
        freq = float(ax.get('freq1', 0) or 0)
        if freq > 1e6:                       # Bruker MWFQ is in Hz
            freq *= 1e-9
        # static field: B0VL of the standard parameter layer, else the CenterField brukerread kept; both in tesla
        try:
            field = float(dsc.get('B0VL', 0) or 0)*1e3 or float(ax.get('cf', 0) or 0)*1e3
        except ValueError as fault:
            print(f"[pyENDOR] 'B0VL' in the header of '{path}' is not a number ({fault}); B0 stays 0, set it in the Data tab")
            field = 0.0
        text = (str(ax.get('title', '')) + ' ' + str(dsc.get('PlsSPELEXPSlct', ''))).lower()
        lo, hi, w = float(x[0]), float(x[-1]), float(x[-1] - x[0])
        return {'ax': {'x': x, 'xlabel': 'RF, MHz'},
                'data': y, 'dsc': dsc,
                'fname': os.path.basename(path), 'fullpath': path,
                'title': str(ax.get('title') or os.path.basename(path)),
                'freq': freq, 'field': field,
                # sequence named in the title or the PulseSPEL experiment, else none; its delays are typed in the Data tab
                'seq': 'Mims' if 'mims' in text else 'Davies' if 'davies' in text else 'ideal',
                'tau': 0.0,                       # Mims: delay between the first two pulses [ns]
                'tinv': 0.0,                      # Davies: length of the inversion pulse [ns]
                'show': True, 'scale': 1.0,
                'yshift': 0.0,                    # vertical display offset, in tallest-trace units
                'color': '',                      # '' = automatic: C<autocolor>
                'autocolor': -1,                  # default colour index, fixed when added (-1 = not yet)
                'bglist': [],                     # [{'bgdata': background title, 'scale', 'shift' [MHz]}]
                'blorder': 0, 'blregions': [[lo, lo + 0.1*w], [hi - 0.1*w, hi]],    # polynomial baseline
                'simdata': None, 'simax': None,   # last simulation and its radio-frequency axis
                'simparts': None,                 # its per-nucleus components, indexed as the nuclei (None = not simulated)
                'simscale': 1.0,                  # extra multiplier of the simulation
                'orisel': None}                   # precalculated orientation selection [phi, theta, weight]

    def addDataset(self, dd, BG=False):
        target = self.BGData if BG else self.Data
        # titles identify backgrounds, so keep them unique
        base, n = dd['title'], 2
        while dd['title'] in [d['title'] for d in target]:
            dd['title'] = f"{base} #{n}"
            n += 1
        target.append(dd)
        if not BG:
            self.assign_autocolor(dd)
            self.current = len(self.Data) - 1
        self.update_everything()

    def assign_autocolor(self, dd):
        """Give dd the lowest default colour index no other spectrum holds; it stays with the spectrum."""
        if dd.get('autocolor', -1) >= 0:
            return
        used = {d.get('autocolor', -1) for d in self.Data if d is not dd}
        dd['autocolor'] = next(k for k in range(len(used) + 1) if k not in used)

    def trace_color(self, dd):
        return dd['color'] or f"C{dd['autocolor'] % 10}"

    def polybaseline(self, x, y, regions, order):
        """Polynomial of `order` fitted to the points of y inside `regions` (list of [lo, hi] in x), evaluated on all of x."""
        m = np.zeros(x.shape, dtype=bool)
        for lo, hi in regions:
            m |= (x >= min(lo, hi)) & (x <= max(lo, hi))
        if order < 0 or m.sum() <= order:
            return np.zeros_like(y)
        xs = (x - x.mean())/np.ptp(x)        # centred/scaled: keeps high orders well conditioned
        return np.polynomial.polynomial.polyval(xs, np.polynomial.polynomial.polyfit(xs[m], y[m], order))

    def process(self, dd):
        """background -> polynomial baseline. x, raw, bg, y0 (raw - bg), bl and y (final spectrum) share the frequency axis."""
        x = dd['ax']['x']
        warn = []
        raw = dd['data']
        bg = np.zeros_like(x)
        for b in dd['bglist']:
            bd = next((d for d in self.BGData if d['title'] == b['bgdata']), None)
            if bd is None:
                warn.append(f"background '{b['bgdata']}' is not loaded")
                continue
            # cubic spline: the background may have its own range and point count
            bgi = CubicSpline(bd['ax']['x'] + b['shift'], bd['data'], extrapolate=False)(x)
            if np.isnan(bgi).any():
                warn.append(f"'{bd['title']}' does not cover the whole frequency range (set to 0 outside)")
            bg += b['scale']*np.nan_to_num(bgi)
        y0 = raw - bg
        bl = self.polybaseline(x, y0, dd['blregions'], dd['blorder'])
        return {'x': x, 'raw': raw, 'bg': bg, 'y0': y0, 'bl': bl, 'y': y0 - bl, 'warn': warn}

    # -----------------------------------------------------------------------
    # Simulation
    # -----------------------------------------------------------------------
    def on_update_orisel(self, event):
        """OriSelType 'precalculated': compute and keep every spectrum's orientation selection; any other type drops them."""
        if self.Opt.OriSelType[0] != 'precalculated':
            for dd in self.Data:
                dd['orisel'] = None
            return
        self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
        es = ENDORsim(Sys=self.Sys, errorFunc=self.RaiseError)
        wait = wx.BusyCursor()
        try:
            for dd in self.Data:
                self.actuateSimMethod(es, dd)
                dd['orisel'] = es.computeOrisel_eig(grid=es.grid_name(), nKnots=es.Opt.nKnots, epsilon=0.33,
                                                    returnMatrix=True, setActive=False).copy()
            self.statusbar.SetStatusText(f"Orientation selection updated for {len(self.Data)} spectra", 0)
        except Exception as e:
            print(f"[pyENDOR] Orientation selection failed: {type(e).__name__}: {e}")
            self.statusbar.SetStatusText(f"Orientation selection failed: {e}", 0)
        finally:
            del wait

    def actuateSimMethod(self, es, dd):
        """Resolve the 'Exp Over.' and Options tabs for one spectrum (-1 = take from the spectrum) into es."""
        tp = self.tabulated_panel
        Exp = tp.Exp_param.parameters.copy()
        x = dd['ax']['x']
        for key, value in (('mwFreq', dd['freq']), ('Field', dd['field']), ('tau', dd['tau']), ('Tinv', dd['tinv']),
                           ('MinFreq', float(x[0])), ('MaxFreq', float(x[-1])), ('nPoints', int(x.size))):
            if Exp[key] < 0:
                Exp[key] = value
        seq = Exp['Sequence'][0] if isinstance(Exp['Sequence'], (list, tuple)) else Exp['Sequence']
        Exp['Sequence'] = seq if seq in SEQUENCES else dd['seq']
        Opt = tp.Opt_param.parameters.copy()
        Opt['OriSelInp'] = dd['orisel']       # always this spectrum's own, never one left over from another
        self.Opt.setFromCtrl(Opt)
        self.Exp.setFromCtrl(Exp)
        es.Opt = self.Opt
        es.Exp = self.Exp

    def runSim(self):
        self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
        es = ENDORsim(Sys=self.Sys, errorFunc=self.RaiseError)
        es.verbose = False
        failed = []
        total = len(self.Data)
        if total:
            self.sim_gauge.SetRange(total)
            self.sim_gauge.SetValue(0)
            self.sim_gauge.Show()
            self.statusbar.Layout()
        start_time = time.time()
        try:
            for ii, dd in enumerate(self.Data):
                self.statusbar.SetStatusText(f"Simulating {ii + 1}/{total}...", 0)
                self.statusbar.Update()
                dd['simdata'] = dd['simax'] = dd['simparts'] = None
                try:
                    self.actuateSimMethod(es, dd)
                    # what a freshly loaded file usually still lacks: named in the status bar, not in a dialog per run
                    seq = es.sequence_name()
                    missing = ("B0 is not set" if not self.Exp.Field > 0 else
                               "Mims needs tau" if seq == 'Mims' and not self.Exp.tau > 0 else
                               "Davies needs t.inv" if seq == 'Davies' and not self.Exp.Tinv > 0 else "")
                    if missing:
                        failed.append(f"{dd['title']} ({missing})")
                    else:
                        es.reRun()
                        dd['simdata'] = es.Spectrum
                        dd['simparts'] = es.Parts
                        dd['simax'] = {'x': es.X, 'xlabel': 'RF, MHz'}
                except Exception as e:
                    # RaiseError already showed a dialog for the errors the simulation anticipated; this is the paper trail
                    print(f"[pyENDOR] Simulation failed for '{dd['title']}': {type(e).__name__}: {e}")
                    failed.append(f"{dd['title']} ({e})")
                self.sim_gauge.SetValue(ii + 1)
                self.sim_gauge.Update()
        finally:
            self.sim_gauge.Hide()            # however the series ends, the gauge must not stay at a stale value
            self.statusbar.Layout()
        self.matplotlib_panel.update_graph()
        if failed:
            self.statusbar.SetStatusText(f"No simulation for {len(failed)}/{total} spectra: " + "; ".join(failed), 0)
        elif total and self.statusbar.GetStatusText(0) == "Ready":      # processing warnings of the redraw stay visible
            self.statusbar.SetStatusText(f"Ready. It took {time.time() - start_time:.2f} s", 0)
        else:
            self.statusbar.SetStatusText("Nothing to simulate: load a spectrum first", 0)

    # -----------------------------------------------------------------------
    # Session XML (same encoding as pyHYSCORE / pyCWEPR / pyEPR)
    # -----------------------------------------------------------------------
    def dumpXML(self):
        if hasattr(wx, 'OSX_FILEDIALOG_ALWAYS_SHOW_TYPES'):
            wx.SystemOptions.SetOption(wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES, 1)
        with wx.FileDialog(self, "Save Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultDir=self.file_browser.path, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            pathname = dlg.GetPath()
        mp, tp = self.matplotlib_panel, self.tabulated_panel
        data_dict = {'program': 'pyENDOR', 'version': VERSION, 'Data': self.Data, 'BGData': self.BGData,
                     'Sys': tp.Sys_param.parameters, 'Exp': tp.Exp_param.parameters, 'Opt': tp.Opt_param.parameters,
                     'components': self.components, 'epr': self.epr,
                     'view': {'norm': mp.norm, 'simfit': mp.simfit, 'offset': float(mp.offset), 'xnuc': mp.xnuc,
                              'showRawBG': mp.showRawBG, 'showBaselines': mp.showBaselines, 'showSim': mp.showSim,
                              'showParts': mp.showParts, 'showGrid': mp.showGrid, 'xfix': mp.xfix, 'yfix': mp.yfix,
                              'autoSIM': self.tb_autoSIM.GetValue(), 'current': self.current}}

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
            for key, value in data_dict.items():
                create_element(root, key, value)
            pretty = minidom.parseString(ET.tostring(root, encoding='unicode')).toprettyxml(indent="  ")
            with open(pathname, 'w', encoding='utf-8') as f:
                f.write(pretty)
            self.statusbar.SetStatusText(f"Session saved to {pathname}", 0)
            self.file_browser.refresh_file_list(select_name=os.path.basename(pathname))
        except Exception as e:
            traceback.print_exc()
            wx.MessageBox(f"Could not save session to '{pathname}':\n{type(e).__name__}: {e}",
                          "Save Session Failed", wx.OK | wx.ICON_ERROR)

    def withDefaults(self, saved, defaults):
        """Saved Exp / Opt values laid over the current defaults, so older sessions gain parameters added since."""
        out = defaults.copy()
        for key, val in saved.items():
            dflt = defaults.get(key)
            if isinstance(dflt, list) and len(dflt) == 2 and isinstance(dflt[1], list):
                # [value, [choices]]: the saved selection, but always the current choice list
                sel = val[0] if (isinstance(val, list) and len(val) == 2) else val
                out[key] = [sel if sel in dflt[1] else dflt[0], list(dflt[1])]
            else:
                out[key] = val
        return out

    def loadXML(self, filename, BG=False):
        """A pyENDOR session replaces the workspace; loaded as background, or from another program,
        a session only contributes its 1D spectra."""
        result = self.load_mixed_dict_xml(filename)
        if result.get('program') != 'pyENDOR' or BG:
            n = 0
            for dd in result.get('Data') or []:
                if (isinstance(dd, dict) and isinstance(dd.get('ax'), dict) and isinstance(dd.get('data'), np.ndarray)
                        and np.squeeze(dd['data']).ndim == 1):
                    new = self.newDataset(dd['ax'], dd['data'], dd.get('dsc'), dd.get('fullpath') or filename)
                    new['title'] = str(dd.get('title') or new['title'])
                    self.addDataset(new, BG)
                    n += 1
            if n == 0:
                wx.MessageBox(f"'{os.path.basename(filename)}' contains no 1D spectra.", "Nothing loaded",
                              wx.OK | wx.ICON_INFORMATION)
            return
        if self.Data or self.BGData:
            if wx.MessageBox("Loading a session replaces all loaded spectra and settings. Continue?",
                             "Load session", wx.YES_NO | wx.ICON_QUESTION) != wx.YES:
                return
        restored = {'Data': [], 'BGData': []}
        for key in restored:
            for dd in result.get(key) or []:
                new = self.newDataset(dd['ax'], dd['data'], dd.get('dsc'), dd.get('fullpath', filename))
                # saved values over current defaults, so older sessions gain newer settings
                new.update({k: v for k, v in dd.items() if k in new and k not in ('ax', 'data')})
                restored[key].append(new)
        self.Data, self.BGData = restored['Data'], restored['BGData']
        for dd in self.Data:
            self.assign_autocolor(dd)
        mp, tp = self.matplotlib_panel, self.tabulated_panel
        self.tb_autoSIM.SetValue(False)       # filling the property grids must not start simulations
        if 'Sys' in result:
            tp.Sys_param.SetFromParClean(result['Sys'])
        if 'Exp' in result:
            tp.Exp_param.SetFromParClean(self.withDefaults(result['Exp'], self.Exp.getDefaultDict()))
        if 'Opt' in result:
            tp.Opt_param.SetFromParClean(self.withDefaults(result['Opt'], self.Opt.getDefaultDict()))
        self.Sys.setFromCtrl(tp.Sys_param.parameters)
        self.Exp.setFromCtrl(tp.Exp_param.parameters)
        self.Opt.setFromCtrl(tp.Opt_param.parameters)
        tp.updateGridInfo()
        # scales and colours of the components; update_components fits the list to the nuclei of the loaded spin system
        self.components = [{'scale': float(c.get('scale', 1.0)), 'color': str(c.get('color') or COMPONENT_COLORS[k % len(COMPONENT_COLORS)])}
                           for k, c in enumerate(result.get('components') or []) if isinstance(c, dict)]
        tp.update_components()
        saved = result.get('epr') or {}
        self.epr = {'source': int(saved.get('source', 0)), 'xaxis': 'g' if saved.get('xaxis') == 'g' else 'field',
                    'nKnots': int(saved.get('nKnots', 0))}
        view = result.get('view') or {}
        mp.norm = view.get('norm') if view.get('norm') in ('none', 'max') else 'none'
        mp.simfit = view.get('simfit') if view.get('simfit') in ('max', 'lsq') else 'max'
        mp.offset = float(view.get('offset', 0.0))
        mp.xnuc = view.get('xnuc') if view.get('xnuc') in self.Sys.isotopes() else ''
        tp.choice_xnuc.SetSelection(tp.choice_xnuc.FindString(mp.xnuc) if mp.xnuc else 0)
        for name, chk in (('showRawBG', mp.chk_RawBG), ('showBaselines', mp.chk_Baselines), ('showSim', mp.chk_Sim),
                          ('showParts', mp.chk_Parts), ('showGrid', self.tb_grid)):
            setattr(mp, name, bool(view.get(name, getattr(mp, name))))
            chk.SetValue(getattr(mp, name))
        # sessions saved before the x axis had a zero choice hold one plain (lo, hi)
        xfix = view.get('xfix') or {}
        mp.xfix = {k: tuple(v) for k, v in xfix.items()} if isinstance(xfix, dict) else {'': tuple(xfix)}
        mp.yfix = tuple(view['yfix']) if view.get('yfix') else None
        tp.seg_norm.SetSelection(('none', 'max').index(mp.norm))
        tp.seg_simfit.SetSelection(('max', 'lsq').index(mp.simfit))
        tp.spin_offset.SetValue(mp.offset)
        self.tb_autoSIM.SetValue(bool(view.get('autoSIM', True)))
        cur = view.get('current', 0)
        self.current = cur if isinstance(cur, int) and 0 <= cur < len(self.Data) else len(self.Data) - 1
        self.update_everything()
        self.statusbar.SetStatusText(f"Session {os.path.basename(filename)} loaded", 0)

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
                    array = np.frombuffer(base64.b64decode(numpy_array.find('data').text or ''), dtype=dtype).reshape(shape).copy()
                except Exception as e:
                    print(f"[pyENDOR] Could not reconstruct array '{key}' (dtype={dtype}, shape={shape}) "
                          f"from '{filename}': {type(e).__name__}: {e}")
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


# --------------------------------------------------------------------------- #
# Application entry point
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame(None)
    frame.Show()
    app.MainLoop()
