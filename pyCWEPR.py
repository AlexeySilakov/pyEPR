import os
import re
import ast
import copy
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
from scipy.signal import hilbert
from scipy.integrate import cumulative_trapezoid

from matplotlib.backend_tools import Cursors
from matplotlib.colors import to_hex
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

from theme import theme
from theme.theme import PillButton, PillCheckBox, PillTabBar, SegmentedPill
from theme import classPropGridPanel as MypgPanel
from readers import brukerread as Mybr
from readers.eprcsvread import csvread
from readers.eprxmlread import xmlread

# pip install wxPython numpy matplotlib scipy
VERSION = 0.1
G_FACTOR = 71.4477      # g = G_FACTOR * f[GHz] / B[mT]   (h/muB = 0.0714477 mT/MHz)
HERE = os.path.dirname(os.path.abspath(__file__))   # icons are found from here, never from the CWD


class DictPopup(wx.Dialog):
    """Key/value view of a file header (Magnettech CSV/XML recipe, Bruker DSC, session summary)."""
    def __init__(self, parent, data_dict, title="Header"):
        super().__init__(parent, title=title, size=(560, 640),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        self.list_ctrl = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.BORDER_NONE)
        self.list_ctrl.InsertColumn(0, "Key", width=220)
        self.list_ctrl.InsertColumn(1, "Value", width=300)
        for key, value in data_dict.items():
            # Magnettech CSV entries are [value, description], XML entries [value, unit]
            if isinstance(value, (list, tuple)) and len(value) == 2 and isinstance(value[1], str):
                value = f"{value[0]}   {value[1]}".strip()
            index = self.list_ctrl.InsertItem(self.list_ctrl.GetItemCount(), str(key))
            self.list_ctrl.SetItem(index, 1, str(value))
        vbox.Add(self.list_ctrl, 1, wx.EXPAND | wx.ALL, 8)
        vbox.Add(wx.Button(panel, wx.ID_OK, "OK"), 0, wx.ALIGN_CENTER | wx.BOTTOM, 8)
        panel.SetSizer(vbox)
        theme.apply_theme_to_window(self, theme.get_theme())


class AxisRangeDialog(wx.Dialog):
    """Min / max of one plot axis, or autoscale. After OK, `result` is (lo, hi), or None for autoscale.
    Laid out like the axis: which='x' puts Min and Max side by side, which='y' stacks Max above Min."""
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
        self.ini_file = os.path.join(wx.StandardPaths.Get().GetUserConfigDir(), "pyCWEPR.ini")
        self.path = os.getcwd()
        self.filters = {"EPR": ('.csv', '.dsc', '.xml'), "CSV": ('.csv',), "DSC": ('.dsc',),
                        "XML": ('.xml',), "All": None}
        self.current_filter = "EPR"

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
        self.btn_load.SetToolTip("Load the selected CSV / DSC / Magnettech XML file(s), or the spectra of an XML session.\n"
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
        except OSError:
            pass
        self.lbl_path.SetValue(self.path)
        self.lbl_path.SetInsertionPointEnd()
        self.refresh_file_list()

    def save_path(self):
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.path + "\n")
        except OSError as e:
            print(f"[pyCWEPR] Could not save last folder to {self.ini_file}: {e}")

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
        # a Bruker .DSC/.DTA pair is one dataset; a Magnettech .csv and .xml of the same name are two
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
                try:                             # Magnettech data XML: a spectrum, titled by <Measurement Name=...>
                    with open(full, 'r', encoding='utf-8-sig', errors='replace') as f:
                        head = f.read(4096)
                    if '<ESRXmlFile' in head:
                        image_id = self.data_img_id
                        found = re.search(r'<Measurement\s+Name="([^"]*)"', head)
                        title = found.group(1) if found else ''
                except OSError as e:
                    print(f"[pyCWEPR] Could not read title from '{full}': {e}")
            elif low.endswith(('.csv', '.dsc')):
                image_id = self.data_img_id
                try:
                    if low.endswith('.csv'):     # Magnettech: "Name;<title>"; pyCWEPR export: "# ... export of: <title>"
                        with open(full, 'r', encoding='utf-8-sig', errors='replace') as f:
                            first = f.readline()
                        if first.startswith('Name;'):
                            title = first.split(';', 1)[1].strip()
                        elif first.startswith('# pyCWEPR') and 'export of:' in first:
                            title = first.split('export of:', 1)[1].strip()
                    else:
                        with open(full, 'r', encoding='latin-1') as f:
                            for line in f:
                                if line.startswith('TITL'):
                                    title = line.split("'")[1] if "'" in line else line[4:].strip()
                                    break
                except (OSError, IndexError) as e:
                    print(f"[pyCWEPR] Could not read title from '{full}': {e}")
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
            if ext == '.csv':
                info = csvread(full)[2] or self.parent.readASCII(full)[2]
            elif ext in ('.txt', '.dat', '.asc'):
                info = self.parent.readASCII(full)[2]
            elif ext in ('.dsc', '.dta'):
                info = Mybr.parse_dsc(os.path.splitext(full)[0] + '.dsc')
            elif ext == '.xml' and next(ET.iterparse(full, events=('start',)))[1].tag == 'ESRXmlFile':
                spectra = xmlread(full)
                info = {'File.Measurements': len(spectra), **spectra[0][2]} if len(spectra) > 1 else spectra[0][2]
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
    """Rows (top -> bottom): derivative spectra, absorption (1st integral),
    running double integral. Every shown spectrum is overlaid in every row;
    the current one is drawn thicker and carries the raw/BG/baseline extras."""
    def __init__(self, parent, mainWindow):
        super().__init__(parent)
        self.parent = mainWindow
        self.showRawBG = False
        self.showBaselines = True
        self.showDeriv = True
        self.showAbs = True
        self.showDI = True
        self.showGrid = False
        self.norm = 'none'        # 'none' | 'max' | 'DI'
        self.xaxis = 'field'      # 'field' | 'g'
        self.offset = 0.0         # stacking step, fraction of the largest trace in a row
        self.layout = None        # (rows, xaxis) of the last draw: x-zoom survives redraws while it holds
        self.drag = None          # region / integration-limit drag in progress (see on_press)
        self.axrows = {}          # axes -> row name of the last draw
        self.artists = {}         # (key, j) -> patches / lines of the current spectrum's markup
        self.update_pending = False   # a coalesced redraw is queued (request_update)
        self.xfix = {}            # x-axis mode ('field' / 'g') -> fixed (lo, hi), shared by all rows
        self.yfix = {}            # row ('deriv' / 'abs' / 'dint') -> fixed (lo, hi)
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
        sizerH = wx.WrapSizer(wx.HORIZONTAL)      # the row toggles live on the main toolbar
        for chk in (self.chk_RawBG, self.chk_Baselines):
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

    def request_update(self):
        """update_graph for slider drags: ticks queued during a redraw only change values, one redraw follows.
        Renders at once (not draw_idle), so the ticks pile up behind a finished frame instead of between frames."""
        if not self.update_pending:
            self.update_pending = True
            wx.CallAfter(lambda: (setattr(self, 'update_pending', False), self.update_graph(), self.canvas.draw()))

    def on_key(self, event):
        """Ctrl-Z (anywhere but a text box, which keeps its undo): zoom-to-rectangle for one zoom only."""
        if (event.GetKeyCode() == ord('Z') and event.GetModifiers() == wx.MOD_CONTROL
                and not isinstance(wx.Window.FindFocus(), wx.TextCtrl)):
            if str(self.toolbar.mode) != 'zoom rect':
                self.toolbar.zoom()
            self.zoom_once = True
            return
        event.Skip()

    def xconv(self, x, f):
        """Field [mT] <-> plot x. g = G_FACTOR*f/B is its own inverse; field axis or no frequency: unchanged."""
        return G_FACTOR*f/x if self.xaxis == 'g' and f > 0 else x

    def hit(self, event):
        """Markup of the current spectrum under the mouse, as (key, j, part): key 'int' (integration
        limit j, part 'line') or 'bl1'/'bl2'/'bl3' (region j, part 'lo'/'hi' = an edge, 'move' = inside).
        None if there is nothing within 5 px."""
        mf = self.parent
        row = self.axrows.get(event.inaxes)
        if row is None or mf.current < 0:
            return None
        dd = mf.Data[mf.current]
        px = lambda b: event.inaxes.transData.transform((self.xconv(b, dd['freq']), 0))[0]
        for j, b in enumerate(dd['intrange']):
            if abs(px(b) - event.x) <= 5:
                return ('int', j, 'line')
        key = {'deriv': 'bl1', 'abs': 'bl2', 'dint': 'bl3'}.get(row)
        if key is None or not self.showBaselines:
            return None
        regs = dd[key + 'regions']
        for j, (lo, hi) in enumerate(regs):
            for part, b in (('lo', lo), ('hi', hi)):
                if abs(px(b) - event.x) <= 5:
                    return (key, j, part)
        for j, (lo, hi) in enumerate(regs):
            if min(px(lo), px(hi)) <= event.x <= max(px(lo), px(hi)):
                return (key, j, 'move')
        return None

    def on_mouse_move(self, event):
        mf = self.parent
        if self.drag is not None and self.drag['key'] == 'line':
            # Vertical line drag: mouse height read in the row where it started, so leaving that row is fine;
            # every row's line of this spectrum moves by the same number of tallest-trace units
            ax0 = self.drag['ax']
            y = ax0.transData.inverted().transform((event.x, event.y))[1]
            d = (y - self.drag['y0'])/self.rowheight[ax0]
            mf.Data[self.drag['ii']]['yshift'] = self.drag['orig'] + d
            for ln, base in self.drag['lines']:
                ln.set_ydata(base + d*self.rowheight[ln.axes])
            self.canvas.draw_idle()
        elif self.drag is not None and self.drag['key'] == 'scale':
            # Scale drag: the grabbed point follows the mouse, the line's zero (its offset) stays put.
            # Grabbed near zero there is nothing to follow: then up = larger, one row height = x e^3
            ax0, ii = self.drag['ax'], self.drag['ii']
            y = ax0.transData.inverted().transform((event.x, event.y))[1]
            off0, r0 = self.lineoff[(ax0, ii)], self.drag['y0'] - self.lineoff[(ax0, ii)]
            if abs(r0) > 0.02*self.rowheight[ax0]:
                factor = max((y - off0)/r0, 0.01)
            else:
                factor = np.exp(3*(y - self.drag['y0'])/self.rowheight[ax0])
            mf.Data[ii]['scale'] = self.drag['orig']*factor
            for ln, base in self.drag['lines']:
                off = self.lineoff[(ln.axes, ii)]
                ln.set_ydata(off + (base - off)*factor)
            mf.statusbar.SetStatusText(f"Scale = {mf.Data[ii]['scale']:.4g}", 0)
            self.canvas.draw_idle()
        elif self.drag is not None and event.inaxes is not None and event.xdata is not None and event.xdata > 0:
            # Live drag: move the stored edge(s), redraw only the dragged patch / limit lines
            dd = mf.Data[mf.current]
            key, j, part = self.drag['key'], self.drag['j'], self.drag['part']
            b = float(self.xconv(event.xdata, dd['freq']))
            if key == 'int':
                dd['intrange'][j] = b
                for ln in self.artists[(key, j)]:
                    ln.set_xdata([self.xconv(b, dd['freq'])]*2)
            else:
                reg = dd[key + 'regions'][j]
                if part == 'move':
                    d = b - self.drag['start']
                    reg[:] = [self.drag['orig'][0] + d, self.drag['orig'][1] + d]
                else:
                    reg[0 if part == 'lo' else 1] = b
                arts = self.artists[(key, j)]
                for k, patch in enumerate(arts):
                    ax = patch.axes
                    patch.remove()
                    arts[k] = ax.axvspan(self.xconv(reg[0], dd['freq']), self.xconv(reg[1], dd['freq']),
                                         color='0.5', alpha=0.15, lw=0)
            self.canvas.draw_idle()
        elif not str(self.toolbar.mode) and event.inaxes is not None:
            h = self.hit(event)
            if (h is None or h[2] == 'move') and any(str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]
                                                     for ln in event.inaxes.get_lines()):
                self.canvas.set_cursor(Cursors.RESIZE_VERTICAL)     # a spectrum line: drag up / down
            else:
                self.canvas.set_cursor(Cursors.POINTER if h is None else
                                       Cursors.MOVE if h[2] == 'move' else Cursors.RESIZE_HORIZONTAL)
        if event.inaxes is None or event.xdata is None:
            mf.statusbar.SetStatusText("", 1)
            return
        f = mf.Data[mf.current]['freq'] if mf.current >= 0 else 0.0
        x = event.xdata
        if self.xaxis == 'g':
            text = f"g = {x:.5f}" + (f",  B = {G_FACTOR*f/x:.3f} mT" if f > 0 and x > 0 else "")
        else:
            text = f"B = {x:.3f} mT" + (f",  g = {G_FACTOR*f/x:.5f}" if f > 0 and x > 0 else "")
        mf.statusbar.SetStatusText(f"{text},  y = {event.ydata:.4g}", 1)

    def on_press(self, event):
        """Double click on a spectrum line (any row): make it the current spectrum; elsewhere: Home view.
        Left drag on the current spectrum's markup: a region (inside = move, edge = resize) or an
        integration limit. Drag on a spectrum line: move it up / down (Y shift); shift-drag on it: scale it.
        Shift-drag on empty space in any row adds a region to that row's baseline."""
        mf = self.parent
        if (event.button != 1 or event.inaxes is None or event.xdata is None or str(self.toolbar.mode)
                or mf.current < 0):
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
        dd = mf.Data[mf.current]
        b = float(self.xconv(event.xdata, dd['freq']))
        h = self.hit(event)
        if h is None or h[2] == 'move':
            # a spectrum line (topmost) under the cursor -- it beats the inside of a region:
            # vertical drag, or with shift held a scale drag
            hits = [ln for ln in event.inaxes.get_lines()
                    if str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]]
            if hits:
                gid = max(hits, key=lambda h: h.get_zorder()).get_gid()
                ii = int(gid[4:])
                scale = wx.GetKeyState(wx.WXK_SHIFT)
                self.drag = {'key': 'scale' if scale else 'line', 'ii': ii, 'ax': event.inaxes, 'y0': event.ydata,
                             'orig': mf.Data[ii]['scale' if scale else 'yshift'],
                             'lines': [(ln, np.asarray(ln.get_ydata(), dtype=float).copy())
                                       for ax in self.figure.axes for ln in ax.get_lines() if ln.get_gid() == gid]}
                return
        new = h is None
        if new:
            key = {'deriv': 'bl1', 'abs': 'bl2', 'dint': 'bl3'}.get(self.axrows.get(event.inaxes))
            if key is None or not self.showBaselines or not wx.GetKeyState(wx.WXK_SHIFT):
                return
            dd[key + 'regions'].append([b, b])
            h = (key, len(dd[key + 'regions']) - 1, 'hi')
            self.update_graph()          # creates the patch the drag then reshapes
        key, j, part = h
        self.drag = {'key': key, 'j': j, 'part': part, 'start': b, 'new': new,
                     'orig': list(dd['intrange'] if key == 'int' else dd[key + 'regions'][j])}

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
        if event.button == 3 and event.inaxes is None and not str(self.toolbar.mode):
            renderer = self.canvas.get_renderer()
            for ax, row in self.axrows.items():
                for which, axis in (('x', ax.xaxis), ('y', ax.yaxis)):
                    bb = axis.get_tightbbox(renderer)
                    if bb is None or not bb.contains(event.x, event.y):
                        continue
                    lo, hi = sorted(ax.get_xlim() if which == 'x' else ax.get_ylim())
                    if which == 'x':
                        title, auto = f"X axis: {ax.get_xlabel()} (all rows)", self.xaxis not in self.xfix
                    else:
                        title, auto = f"Y axis: {ax.get_ylabel()}", row not in self.yfix
                    with AxisRangeDialog(self, title, lo, hi, auto, which) as dlg:
                        if dlg.ShowModal() != wx.ID_OK:
                            return
                        rng = dlg.result
                    target, key = (self.xfix, self.xaxis) if which == 'x' else (self.yfix, row)
                    if rng is None:
                        target.pop(key, None)
                    else:
                        target[key] = rng
                    self.update_graph()
                    return
            return
        if event.button == 3 and event.inaxes is not None and not str(self.toolbar.mode):
            h = self.hit(event)
            region = h if h is not None and h[0] != 'int' else None
            hits = [ln for ln in event.inaxes.get_lines()
                    if str(ln.get_gid()).startswith('spec') and ln.contains(event)[0]]
            if not hits:
                if region is not None:
                    del mf.Data[mf.current][region[0] + 'regions'][region[1]]
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
                del mf.Data[mf.current][region[0] + 'regions'][region[1]]
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
        key, j, new = self.drag['key'], self.drag['j'], self.drag['new']
        self.drag = None
        dd = mf.Data[mf.current]
        if key == 'int':
            dd['intrange'].sort()
        else:
            regs = dd[key + 'regions']
            regs[j].sort()
            if new and regs[j][1] - regs[j][0] < 0.002*np.ptp(dd['ax']['x']):    # shift-click, not a drag
                del regs[j]
        mf.tabulated_panel.update_controls()
        self.update_graph()

    def on_check(self, event):
        self.showRawBG = self.chk_RawBG.GetValue()
        self.showBaselines = self.chk_Baselines.GetValue()
        self.update_graph()

    def update_graph(self):
        mf = self.parent
        tp = mf.tabulated_panel
        rows = ((['deriv'] if self.showDeriv else []) + (['abs'] if self.showAbs else [])
                + (['dint'] if self.showDI else [])) or ['deriv']
        layout = (tuple(rows), self.xaxis)
        keep = layout == self.layout and len(self.figure.axes) == len(rows)
        xlim = self.figure.axes[0].get_xlim() if keep else None
        self.figure.clear()
        axes = [self.figure.add_subplot(len(rows), 1, 1)]
        for kk in range(1, len(rows)):
            axes.append(self.figure.add_subplot(len(rows), 1, kk + 1, sharex=axes[0]))
        self.layout = layout
        self.axrows = dict(zip(axes, rows))
        self.artists = {}

        # every spectrum (the DI readouts need them all); loaded data first, so math spectra reuse their sources
        mf.proc_cache = {}
        for dd in sorted(mf.Data, key=lambda d: d['math'] is not None):
            mf.proc_cache[id(dd)] = mf.process(dd)
        procs = [mf.proc_cache[id(dd)] for dd in mf.Data]
        mf.proc_cache = {}
        warnings = sorted({w for p in procs for w in p['warn']})
        gmode = self.xaxis == 'g'
        toplot = self.xconv
        if gmode and any(dd['show'] and dd['freq'] <= 0 for dd in mf.Data):
            warnings.append("spectra without MW frequency stay on the field axis -- set 'MW Freq' in the Data tab")
        cur = mf.current

        # Scaled traces per row, then the stacking step per row
        traces = {'deriv': [], 'abs': [], 'dint': []}
        facs = {}
        for ii, (dd, p) in enumerate(zip(mf.Data, procs)):
            if not dd['show']:
                continue
            fac = dd['scale']
            if self.norm == 'max' and np.ptp(p['y']) > 0:
                fac /= np.ptp(p['y'])
            elif self.norm == 'DI' and p['DI'] != 0:
                fac /= abs(p['DI'])
            facs[ii] = fac
            traces['deriv'].append((ii, toplot(p['x'], dd['freq']), p['y']*fac))
            traces['abs'].append((ii, toplot(p['xi'], dd['freq']), (p['abs'] - p['bl2'])*fac))
            traces['dint'].append((ii, toplot(p['xi'], dd['freq']), p['dint']*fac))

        self.rowheight = {}
        self.lineoff = {}                # (axes, spectrum index) -> vertical offset of that line, for shift-drag scaling
        for ax, row in zip(axes, rows):
            # tallest trace of the row: unit of the stack step and of each spectrum's own Y shift
            height = max([np.ptp(y) for _, _, y in traces[row]], default=0.0) or 1.0
            self.rowheight[ax] = height
            step = self.offset*height
            for kk, (ii, xp, y) in enumerate(traces[row]):
                col = mf.trace_color(mf.Data[ii])
                off = kk*step + mf.Data[ii]['yshift']*height
                self.lineoff[(ax, ii)] = off
                ax.plot(xp, y + off, color=col, lw=1.6 if ii == cur else 0.9, gid=f'spec{ii}',
                        label=f"{ii + 1}: {mf.Data[ii]['title']}", zorder=3 if ii == cur else 2)
                if ii != cur:
                    continue
                dd, p, fac = mf.Data[ii], procs[ii], facs[ii]
                if row == 'deriv' and self.showRawBG:
                    ax.plot(xp, p['raw']*fac + off, color='0.6', lw=0.8, label='raw')
                    if dd['bglist']:
                        ax.plot(xp, p['bg']*fac + off, color='0.4', lw=0.8, ls='--', label='background')
                if row == 'deriv' and self.showBaselines and dd['bl1order'] >= 0:
                    ax.plot(xp, p['y0']*fac + off, color=col, lw=0.6, alpha=0.6)
                    ax.plot(xp, p['bl1']*fac + off, color='k', lw=0.8, ls=':', label='deriv. baseline')
                if row == 'abs' and self.showBaselines and dd['bl2order'] >= 0:
                    ax.plot(toplot(p['xi'], dd['freq']), p['abs']*fac + off, color=col, lw=0.6, alpha=0.6)
                    ax.plot(toplot(p['xi'], dd['freq']), p['bl2']*fac + off, color='k', lw=0.8, ls=':',
                            label='abs. baseline')
                if row == 'dint' and self.showBaselines and dd['bl3order'] >= 0:
                    ax.plot(toplot(p['xi'], dd['freq']), p['dint0']*fac + off, color=col, lw=0.6, alpha=0.6)
                    ax.plot(toplot(p['xi'], dd['freq']), p['bl3']*fac + off, color='k', lw=0.8, ls=':',
                            label='DI baseline')
                if row == 'dint':
                    # the level DI is read at (last point of the corrected double integral), with its value
                    level = p['dint'][-1]*fac + off
                    ax.axhline(level, color=col, lw=0.8, ls='-.', alpha=0.8, gid='DIlevel')
                    ax.annotate(f"DI = {p['DI']:.5g}", (toplot(p['xi'][-1], dd['freq']), level),
                                xytext=(-4, 3), textcoords='offset points', ha='right', va='bottom',
                                color=col, fontsize='small')
                key = {'deriv': 'bl1', 'abs': 'bl2', 'dint': 'bl3'}.get(row)
                if self.showBaselines and key is not None:
                    for j, (lo, hi) in enumerate(dd[key + 'regions']):
                        self.artists[(key, j)] = [ax.axvspan(toplot(lo, dd['freq']), toplot(hi, dd['freq']),
                                                             color='0.5', alpha=0.15, lw=0)]
                for j, xr in enumerate(dd['intrange']):
                    self.artists.setdefault(('int', j), []).append(
                        ax.axvline(toplot(xr, dd['freq']), color='0.4', lw=0.8, ls='--', gid='int'))
            ax.grid(self.showGrid)
        for ax, row in zip(axes, rows):
            ax.set_ylabel({'deriv': "dχ''/dB", 'abs': "Absorption", 'dint': "Double integral"}[row])
        axes[-1].set_xlabel("g" if gmode else "B, mT")
        if traces['deriv']:
            axes[0].legend(fontsize='small', frameon=False, loc='best')
        else:
            axes[0].text(0.5, 0.5, "Load a spectrum from the file browser", ha='center', va='center',
                         transform=axes[0].transAxes)
        if gmode:
            axes[0].invert_xaxis()
        theme.style_figure(self.figure, axes)

        # Home = the full autoscaled view. A fixed x range (axis right-click) or a kept x-zoom limits x, and
        # y autoscales to what is visible -- except rows whose y range is fixed
        self.toolbar.update()
        self.toolbar.push_current()
        xfix = self.xfix.get(self.xaxis)
        if xfix is not None:
            xlim = (xfix[1], xfix[0]) if gmode else xfix     # g axis runs high -> low
        if xfix is not None or keep:
            axes[0].set_xlim(xlim)
            xlo, xhi = sorted(xlim)
            for ax in axes:
                ys = [np.asarray(ln.get_ydata())[(np.asarray(ln.get_xdata()) >= xlo) & (np.asarray(ln.get_xdata()) <= xhi)]
                      for ln in ax.get_lines() if ln.get_gid() not in ('int', 'DIlevel')]
                ys = np.concatenate(ys) if ys else np.array([])
                if ys.size and np.ptp(ys) > 0:
                    ax.set_ylim(ys.min() - 0.05*np.ptp(ys), ys.max() + 0.05*np.ptp(ys))
        for ax, row in zip(axes, rows):
            if row in self.yfix:
                ax.set_ylim(self.yfix[row])
        self.canvas.draw_idle()

        # Readouts
        for ii, p in enumerate(procs):
            tp.filelist.set_value(ii, 'DI', f"{p['DI']:.6g}")
        tp.filelist.refresh_headers()     # a colour picked on the plot shows on its swatch too
        tp.lbl_DI.SetLabel(f"DI = {procs[cur]['DI']:.6g}" if cur >= 0 else "DI = –")
        mf.statusbar.SetStatusText("; ".join(warnings) if warnings else "Ready", 0)


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
            prop.SetHelpString("Used for the g axis and for the background frequency ratio.")
            prop = grid.Append(wxpg.FloatProperty("Scale", f"scale {ii}", value=dd['scale']))
            grid.SetPropertyEditor(prop, "SpinFloat")
            prop.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.1)
            prop.SetHelpString("Display multiplier (applied on top of the Compare normalisation).\n"
                               "Or shift-drag the spectrum's line on the plot.")
            prop = grid.Append(wxpg.FloatProperty("Y shift", f"yshift {ii}", value=dd['yshift']))
            grid.SetPropertyEditor(prop, "SpinFloat")
            prop.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.1)
            prop.SetHelpString("Vertical display offset, in units of the tallest trace of each plot row\n"
                               "(display only). Or drag the spectrum's line up or down on the plot.")
            prop = grid.Append(wxpg.StringProperty("DI", f"DI {ii}", value=""))
            grid.SetPropertyReadOnly(prop)
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
        self.Thaw()

    def expand(self, ii):
        """Open entry ii (None: none) and close the others."""
        self.open = ii
        self.Freeze()
        for jj, grid in enumerate(self.grids):
            grid.Show(jj == ii)
        self.SetupScrolling(scroll_x=False, rate_y=10, scrollToTop=False)
        self.Thaw()
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
# Tabulated panel: Data | Background | Phase | Integrate
# --------------------------------------------------------------------------- #
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow):
        super().__init__(parent, size=(280, 400))
        self.parent = mainWindow

        # Every tab but Data edits the spectrum picked here
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
        self.seg_norm = SegmentedPill(self.data_panel, ["As is", "Max", "DI"])
        self.seg_norm.SetToolTip("As is: data x Scale\nMax: peak-to-peak of the derivative = 1\n"
                                 "DI: double integral = 1 (equal spin counts)")
        self.seg_norm.Bind(wx.EVT_RADIOBOX, self.on_compare)
        grid.Add(self.seg_norm, 1, wx.EXPAND)
        grid.Add(wx.StaticText(self.data_panel, label="Stack offset:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.spin_offset = wx.SpinCtrlDouble(self.data_panel, min=0, max=5, inc=0.05, initial=0)
        self.spin_offset.SetDigits(2)
        self.spin_offset.SetToolTip("Vertical step between spectra, as a fraction of the largest trace in each row.")
        self.spin_offset.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_compare)
        theme.bind_spin_enter(self.spin_offset, self.on_compare)
        grid.Add(self.spin_offset, 1, wx.EXPAND)
        sizer.Add(grid, 0, wx.EXPAND | wx.ALL, 6)
        self.btn_copyDI = PillButton(self.data_panel, "Copy DI table", color_key="pill_calm",
                                     draw_badge=theme.draw_textlines_badge)
        self.btn_copyDI.SetToolTip("Copy title, MW frequency and double integral of every spectrum\n"
                                   "to the clipboard (tab separated, pastes into Excel/Origin).")
        self.btn_copyDI.Bind(wx.EVT_BUTTON, self.on_copyDI)
        sizer.Add(self.btn_copyDI, 0, wx.EXPAND | wx.ALL, 4)
        self.btn_export = PillButton(self.data_panel, "Export ASCII", color_key="pill_calm",
                                     draw_badge=theme.draw_textlines_badge)
        self.btn_export.SetToolTip("Save traces of the current spectrum (derivative, absorption, double integral, ...)\n"
                                   "as a comma-separated .csv that pyCWEPR, Excel, Origin or pandas can read.")
        self.btn_export.Bind(wx.EVT_BUTTON, lambda e: self.parent.exportASCII())
        sizer.Add(self.btn_export, 0, wx.EXPAND | wx.ALL, 4)
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
        self.btn_applyBG.SetToolTip("Copy this background list (with scales, ratios, shifts) to every spectrum.")
        self.btn_applyBG.Bind(wx.EVT_BUTTON, lambda e: self.on_apply_all(('bglist',)))
        sizer.Add(self.btn_applyBG, 0, wx.EXPAND | wx.ALL, 3)
        txt = wx.StaticText(self.bg_panel, label=(
            "Background field axis is mapped as B·ratio + shift, then\n"
            "cubic-spline interpolated onto the spectrum's axis.\n"
            "Freq. ratio −1 = f(spectrum)/f(background), i.e. the\n"
            "same g-value lands on the same point."))
        sizer.Add(txt, 0, wx.ALL, 4)
        self.bg_panel.SetSizer(sizer)

        # ---------------- Phase tab
        self.phase_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.phase_panel, "Phase")
        sizer = wx.BoxSizer(wx.VERTICAL)
        # Two identical blocks: 'lockin' recombines the lock-in quadratures (first step, Magnettech XML only),
        # 'phase' mixes in the Hilbert-transform dispersion (after background and derivative baseline)
        self.slider_phase = {}
        self.spin_phase = {}
        self.lbl_lockin0 = self.btn_lockin0 = None
        for key, label, info in (
                ('lockin', "Lock-in (modulation) phase θ",
                 "Magnettech XML Quad-detect enabled Bruker DSC"),
                ('phase', "Absorption / dispersion phase φ",
                 "Mixes in the dispersion (Hilbert transform):\n"
                 "   y = s·cos φ − H[s]·sin φ\n")):
            parent = theme.Card(self.phase_panel, label)
            box = parent.body
            box.Add(wx.StaticText(parent, label=info), 0, wx.ALL, 3)
            slider = wx.Slider(parent, value=0, minValue=-1800, maxValue=1800)
            slider.SetToolTip("0.1° per tick (drag for a live update)")
            slider.Bind(wx.EVT_SLIDER, lambda e, k=key: self.on_phase_slider(k))
            box.Add(slider, 0, wx.EXPAND | wx.ALL, 3)
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(wx.StaticText(parent, label="[deg]:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
            spin = wx.SpinCtrlDouble(parent, min=-360, max=360, inc=0.5, initial=0)
            spin.SetDigits(2)
            spin.SetToolTip("Right-click to set the arrow step")
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=key: self.on_phase_spin(k))
            theme.bind_spin_enter(spin, lambda k=key: self.on_phase_spin(k))
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                spin.Bind(evt, lambda e, s=spin: self.on_step_menu(s, (0.01, 0.1, 0.5, 1, 5, 10, 45, 90), 2))
            row.Add(spin, 1, wx.EXPAND)
            if key == 'lockin':
                self.lbl_lockin0 = wx.StaticText(parent, label="file: –")
                row.Add(self.lbl_lockin0, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 6)
                self.btn_lockin0 = wx.Button(parent, label="↺", size=(24, 22), style=wx.BORDER_NONE)
                self.btn_lockin0.SetToolTip("Back to the phase stored in the file")
                self.btn_lockin0.Bind(wx.EVT_BUTTON, lambda e: self.on_phase_spin('lockin', reset=True))
                row.Add(self.btn_lockin0, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 3)
            box.Add(row, 0, wx.EXPAND | wx.ALL, 3)
            sizer.Add(parent, 0, wx.EXPAND | wx.ALL, 5)
            self.slider_phase[key] = slider
            self.spin_phase[key] = spin
        self.btn_applyPhase = PillButton(self.phase_panel, "Apply to all spectra", color_key="pill_calm")
        self.btn_applyPhase.SetToolTip("Copy both phases to every spectrum (θ matters only for XML spectra).")
        self.btn_applyPhase.Bind(wx.EVT_BUTTON, lambda e: self.on_apply_all(('lockin', 'phase')))
        sizer.Add(self.btn_applyPhase, 0, wx.EXPAND | wx.ALL, 6)
        self.phase_panel.SetSizer(sizer)

        # ---------------- Integrate tab
        self.int_panel = wx.ScrolledWindow(self.nb_book)     # scrolls once there are many fit ranges
        self.int_panel.SetScrollRate(0, 10)
        self.nb_book.AddPage(self.int_panel, "Integrate")
        sizer = wx.BoxSizer(wx.VERTICAL)
        card = theme.Card(self.int_panel, "Integration range")
        grid = wx.FlexGridSizer(cols=2, vgap=4, hgap=6)
        grid.AddGrowableCol(1)
        grid.Add(wx.StaticText(card, label="From [mT]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.spin_from = wx.SpinCtrlDouble(card, min=-1e6, max=1e6, inc=0.1, initial=0)
        grid.Add(self.spin_from, 1, wx.EXPAND)
        grid.Add(wx.StaticText(card, label="To [mT]:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.spin_to = wx.SpinCtrlDouble(card, min=-1e6, max=1e6, inc=0.1, initial=0)
        grid.Add(self.spin_to, 1, wx.EXPAND)
        for spin in (self.spin_from, self.spin_to):
            spin.SetDigits(3)
            spin.SetToolTip("Right-click to set the arrow step")
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_range)
            theme.bind_spin_enter(spin, self.on_range)
            for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                spin.Bind(evt, lambda e, s=spin: self.on_step_menu(s, (0.001, 0.01, 0.1, 1, 10, 100), 3))
        card.body.Add(grid, 0, wx.EXPAND | wx.ALL, 3)
        sizer.Add(card, 0, wx.EXPAND | wx.ALL, 5)

        # Three identical baseline blocks: 'bl1' on the derivative, 'bl2' on the absorption, 'bl3' on the double integral
        self.spin_order = {}
        self.region_rows = {}
        for key, label in (('bl1', "Derivative baseline (before 1st integration)"),
                           ('bl2', "Absorption baseline (before 2nd integration)"),
                           ('bl3', "Double-integral baseline (after 2nd integration)")):
            card = theme.Card(self.int_panel, label)
            box = card.body
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(wx.StaticText(card, label="Polynomial order:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
            spin = wx.SpinCtrl(card, min=-1, max=12, initial=0)
            spin.SetToolTip("−1 = no baseline subtraction; 0 = offset; 1 = line; ...")
            spin.Bind(wx.EVT_SPINCTRL, lambda e, k=key: self.on_order(k))
            theme.bind_spin_enter(spin, lambda k=key: self.on_order(k))
            row.Add(spin, 1, wx.EXPAND)
            box.Add(row, 0, wx.EXPAND | wx.ALL, 3)
            # one 'from … to … ×' row per fit range, rebuilt by update_controls when their number changes
            grid = wx.FlexGridSizer(cols=5, vgap=3, hgap=4)
            grid.AddGrowableCol(1)
            grid.AddGrowableCol(3)
            box.Add(grid, 0, wx.EXPAND | wx.ALL, 3)
            row = wx.BoxSizer(wx.HORIZONTAL)
            for text, action, tip in (("+ Add range", 'add', "Add a fit range after the last one"),
                                      ("Clear all", 'clear', "Delete all fit ranges of this baseline")):
                btn = wx.Button(card, label=text, style=wx.BORDER_NONE)
                btn.SetToolTip(tip)
                btn.Bind(wx.EVT_BUTTON, lambda e, k=key, a=action: self.on_region(k, a))
                row.Add(btn, 0, wx.ALIGN_CENTER_VERTICAL)
                row.AddStretchSpacer()
            box.Add(row, 0, wx.EXPAND | wx.ALL, 3)
            sizer.Add(card, 0, wx.EXPAND | wx.ALL, 5)
            self.spin_order[key] = spin
            self.region_rows[key] = (card, grid, [])      # [] -> (from spin, to spin) per range

        hint = wx.StaticText(self.int_panel, label=(
            "On the plot (Zoom/Pan off):\n"
            "• drag a grey region to move it, or its edge to resize\n"
            "• Shift-drag in any row to add one to that row's baseline\n"
            "• right-click a region to delete it\n"
            "• drag a dashed line to change the integration range"))
        sizer.Add(hint, 0, wx.EXPAND | wx.ALL, 6)

        self.lbl_DI = wx.StaticText(self.int_panel, label="DI = –")
        font = self.lbl_DI.GetFont()
        font.SetPointSize(font.GetPointSize() + 3)
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.lbl_DI.SetFont(font)
        self.lbl_DI.SetToolTip("Double integral of the current spectrum over the integration range\n"
                               "(data units × mT²), after both baseline subtractions.")
        sizer.Add(self.lbl_DI, 0, wx.EXPAND | wx.ALL, 8)
        self.btn_applyInt = PillButton(self.int_panel, "Apply to all spectra", color_key="pill_calm")
        self.btn_applyInt.SetToolTip("Copy range, polynomial orders and regions to every spectrum.")
        self.btn_applyInt.Bind(wx.EVT_BUTTON, lambda e: self.on_apply_all(
            ('intrange', 'bl1order', 'bl1regions', 'bl2order', 'bl2regions', 'bl3order', 'bl3regions')))
        sizer.Add(self.btn_applyInt, 0, wx.EXPAND | wx.ALL, 4)
        self.int_panel.SetSizer(sizer)

        # ---------------- Math tab
        self.math_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.math_panel, "Math")
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.txt_expr = wx.TextCtrl(self.math_panel, style=wx.TE_PROCESS_ENTER)
        self.txt_expr.SetHint("e.g. s1 - a*s2")
        self.txt_expr.SetToolTip("s1, s2 … = spectra by their number in the Data list.\n"
                                 "Any other name (a, k, c1 …) is a constant with a slider.\n"
                                 "Also: numbers, + - * / **, abs max min mean sqrt,\n"
                                 "diff (derivative d/dB, per mT).\n"
                                 "e.g.   s1 - a*s2      (s1 + s2)/2      k*(s3 - s4)\n"
                                 "Enter = Update.")
        self.txt_expr.Bind(wx.EVT_TEXT_ENTER, self.on_expr)
        sizer.Add(self.txt_expr, 0, wx.EXPAND | wx.ALL, 6)
        # quick-insert buttons: text goes in at the cursor (functions with the cursor between the brackets)
        for tokens, color_key in ((("s1", "s2", "s3", "s4"), "pill_load"),     # spectra
                                  (("a", "b", "c", "d"), "pill_calm"),         # constants
                                  (("abs", "diff", "(", ")"), "accent"),
                                  (("+", "-", "/", "*"), "pill_warn")):
            keys = wx.BoxSizer(wx.HORIZONTAL)         # one group of four per row, equal widths
            for token in tokens:
                btn = PillButton(self.math_panel, token, color_key=color_key)
                btn.Bind(wx.EVT_BUTTON, lambda e, t=token: self.on_math_insert(t))
                keys.Add(btn, 1, wx.EXPAND | wx.ALL, 1)
            sizer.Add(keys, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 5)
        row = wx.BoxSizer(wx.HORIZONTAL)
        for kind, label, tip in (('mean', "Average", "(s1 + s2 + …)/n"),
                                 ('weighted', "Weighted", "c1*s1 + c2*s2 + …, a slider per weight,")):
            btn = PillButton(self.math_panel, label, color_key="pill_calm")
            btn.SetToolTip(f"Fill in {tip} of all loaded spectra (math spectra excluded)")
            btn.Bind(wx.EVT_BUTTON, lambda e, k=kind: self.on_math_fill(k))
            row.Add(btn, 1, wx.EXPAND | wx.ALL, 1)
        sizer.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 5)
        # Box sizers, not a FlexGridSizer: a SegmentedPill's min width is -1, and a FlexGridSizer column
        # holding only such items counts as hidden -- the pills would never be placed (left at the top-left)
        self.seg_align = SegmentedPill(self.math_panel, ["Field", "g"])
        self.seg_align.SetToolTip("Field: points matched in mT.\n"
                                  "g: each source's field is first scaled by f(result)/f(source).")
        self.seg_stage = SegmentedPill(self.math_panel, ["Raw", "Processed"])
        self.seg_stage.SetToolTip("Raw: the loaded data (after the lock-in phase).\n"
                                  "Processed: each source's final derivative (after its background,\n"
                                  "derivative baseline and phase).")
        for label, seg in (("Align:", self.seg_align), ("Stage:", self.seg_stage)):
            seg.Bind(wx.EVT_RADIOBOX, self.on_math_mode)
            row = wx.BoxSizer(wx.HORIZONTAL)
            row.Add(wx.StaticText(self.math_panel, label=label, size=(45, -1)), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 6)
            row.Add(seg, 1, wx.EXPAND)
            sizer.Add(row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 6)
        row = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_math_new = PillButton(self.math_panel, "Create", color_key="pill_load", draw_badge=theme.draw_plus_badge)
        self.btn_math_new.SetToolTip("Add the expression as a new math spectrum. It stays live: its sources\n"
                                     "and constants can still be changed, and it follows them.")
        self.btn_math_new.Bind(wx.EVT_BUTTON, lambda e: self.on_expr(create=True))
        row.Add(self.btn_math_new, 1, wx.EXPAND | wx.ALL, 1)
        self.btn_math = PillButton(self.math_panel, "Update", color_key="pill_load")
        self.btn_math.SetToolTip("Apply the expression to the current math spectrum (or the last one in the\n"
                                 "Data list); creates one if there is none yet. Enter does the same.")
        self.btn_math.Bind(wx.EVT_BUTTON, lambda e: self.on_expr())
        row.Add(self.btn_math, 1, wx.EXPAND | wx.ALL, 1)
        sizer.Add(row, 0, wx.EXPAND | wx.ALL, 4)
        self.lbl_refs = wx.StaticText(self.math_panel, label="")
        sizer.Add(self.lbl_refs, 0, wx.EXPAND | wx.ALL, 6)
        # constants of the current math spectrum: name | slider | value | ± slider range (rebuilt per expression)
        self.const_panel = wx.Panel(self.math_panel)
        self.const_widgets = {}
        self.const_key = None
        sizer.Add(self.const_panel, 1, wx.EXPAND | wx.ALL, 4)
        self.math_panel.SetSizer(sizer)

        self.nb = PillTabBar(self, self.nb_book, ["Data", "Background", "Phase", "Integrate", "Math"])
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(row_cur, 0, wx.EXPAND | wx.ALL, 4)
        sizer.Add(self.nb, 0, wx.EXPAND)
        sizer.Add(self.nb_book, 1, wx.EXPAND)
        self.SetSizer(sizer)

    # ------------------------------------------------------------------
    def update_filelist(self):
        """Rebuild the spectrum chooser, the Data grid and the background list."""
        mf = self.parent
        titles = [f"{ii + 1}: {dd['title']}" for ii, dd in enumerate(mf.Data)]
        self.choice_current.SetItems(titles)
        if mf.current >= 0:
            self.choice_current.SetSelection(mf.current)

        self.filelist.rebuild()

        self.list_bg.Set([bd['title'] for bd in mf.BGData])
        if mf.BGData:
            self.list_bg.SetSelection(0)
        self.update_controls()

    def update_controls(self):
        """Push the current spectrum's processing settings into the Background/Phase/Integrate tabs."""
        mf = self.parent
        self.bgtree.Clear()
        if mf.current < 0:
            for key, (card, grid, spins) in self.region_rows.items():
                grid.Clear(delete_windows=True)
                spins.clear()
            self.int_panel.FitInside()
            self.lbl_DI.SetLabel("DI = –")
            self.const_panel.DestroyChildren()
            self.const_widgets, self.const_key = {}, None
            self.lbl_refs.SetLabel("")
            return
        dd = mf.Data[mf.current]
        titles = [bd['title'] for bd in mf.BGData]
        for jj, b in enumerate(dd['bglist']):
            ratio, bd = mf.bgratio(dd, b)
            cat = self.bgtree.Append(wxpg.PropertyCategory(f"BG {jj + 1}:  B·{ratio:.6f}", f"cat {jj}"))
            labels = titles if bd is not None else titles + [b['bgdata'] + " (missing)"]
            self.bgtree.AppendIn(cat, wxpg.EnumProperty("Background", f"bgdata {jj}", labels=labels,
                                                        value=labels.index(b['bgdata']) if bd is not None else len(titles)))
            for key, label, step, helpstr in (
                    ('scale', "Scale", 0.01, "Multiplier of this background before subtraction."),
                    ('fratio', "Freq. ratio", 1e-5, "Field axis of the background is multiplied by this.\n"
                                                    "−1 = f(spectrum)/f(background) from the file headers."),
                    ('shift', "Field shift [mT]", 0.01, "Added to the (scaled) background field axis.")):
                prop = self.bgtree.AppendIn(cat, wxpg.FloatProperty(label, f"{key} {jj}", value=b[key]))
                self.bgtree.SetPropertyEditor(prop, "SpinFloat")
                prop.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, step)
                prop.SetHelpString(helpstr)
        self.bgtree.ExpandAll()

        for key in self.spin_phase:
            self.slider_phase[key].SetValue(int(round(np.clip(dd[key], -180, 180)*10)))
            self.spin_phase[key].SetValue(dd[key])
        has_quad = dd['quad'] is not None
        self.slider_phase['lockin'].Enable(has_quad)
        self.spin_phase['lockin'].Enable(has_quad)
        self.btn_lockin0.Enable(has_quad)
        self.lbl_lockin0.SetLabel(f"file: {dd['lockin0']:g}°" if has_quad else "no sin/cos traces")
        self.spin_from.SetValue(min(dd['intrange']))
        self.spin_to.SetValue(max(dd['intrange']))
        for key, (card, grid, spins) in self.region_rows.items():
            self.spin_order[key].SetValue(dd[key + 'order'])
            regs = dd[key + 'regions']
            if len(spins) != len(regs):
                grid.Clear(delete_windows=True)
                spins.clear()
                for j in range(len(regs)):
                    pair = []
                    for text in ("from", "to"):
                        grid.Add(wx.StaticText(card, label=text), 0, wx.ALIGN_CENTER_VERTICAL)
                        spin = wx.SpinCtrlDouble(card, min=-1e6, max=1e6, inc=0.1, initial=0)
                        spin.SetDigits(3)
                        spin.SetToolTip("Fit range in mT (Enter to apply). Right-click to set the arrow step.")
                        spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=key, j=j: self.on_region(k, 'edit', j))
                        theme.bind_spin_enter(spin, lambda k=key, j=j: self.on_region(k, 'edit', j))
                        for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                            spin.Bind(evt, lambda e, s=spin: self.on_step_menu(s, (0.001, 0.01, 0.1, 1, 10, 100), 3))
                        grid.Add(spin, 1, wx.EXPAND)
                        pair.append(spin)
                    btn = wx.Button(card, label="×", size=(24, 22), style=wx.BORDER_NONE)
                    btn.SetToolTip("Delete this range")
                    btn.Bind(wx.EVT_BUTTON, lambda e, k=key, j=j: self.on_region(k, 'delete', j))
                    grid.Add(btn, 0, wx.ALIGN_CENTER_VERTICAL)
                    spins.append(tuple(pair))
                theme.apply_theme_to_window(card, theme.get_theme())
                self.int_panel.Layout()
                self.int_panel.FitInside()
            for (s_lo, s_hi), (lo, hi) in zip(spins, regs):
                s_lo.SetValue(lo)
                s_hi.SetValue(hi)

        # Math tab: expression (renumbered to the current Data order), sources, and a slider row per constant
        mm = dd['math']
        consts = []
        if mm is None:
            self.lbl_refs.SetLabel("")
        else:
            titles = [d['title'] for d in mf.Data]
            renumber = lambda m: (f"s{titles.index(mm['refs'][m.group()]) + 1}"
                                  if mm['refs'].get(m.group()) in titles else m.group())
            self.txt_expr.SetValue(re.sub(r'\bs\d+\b', renumber, mm['expr']))
            self.seg_align.SetSelection(('field', 'g').index(mm['align']))
            self.seg_stage.SetSelection(('raw', 'processed').index(mm['stage']))
            self.lbl_refs.SetLabel("\n".join(f"s{titles.index(t) + 1} = {t}" if t in titles else f"{s} = {t}  (not loaded)"
                                             for s, t in mm['refs'].items()))
            try:
                consts = mf.parse_expr(mm['expr'])[2]
            except (SyntaxError, ValueError):
                consts = []
        if (id(dd), tuple(consts)) != self.const_key:
            self.const_panel.DestroyChildren()
            self.const_widgets = {}
            grid = wx.FlexGridSizer(cols=4, vgap=2, hgap=4)
            grid.AddGrowableCol(1)
            if consts:
                for label in ("", "drag", "value", "± range"):
                    grid.Add(wx.StaticText(self.const_panel, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            for c in consts:
                grid.Add(wx.StaticText(self.const_panel, label=c), 0, wx.ALIGN_CENTER_VERTICAL)
                slider = wx.Slider(self.const_panel, value=500, minValue=0, maxValue=1000)
                slider.Bind(wx.EVT_SLIDER, lambda e, k=c: self.on_const(k, 'slider'))
                spin = wx.SpinCtrlDouble(self.const_panel, min=-1e9, max=1e9, inc=0.01, initial=0, size=(95, -1))
                spin.SetDigits(4)
                spin.SetToolTip("Right-click to set the arrow step")
                spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=c: self.on_const(k, 'value'))
                theme.bind_spin_enter(spin, lambda k=c: self.on_const(k, 'value'))
                # right-click step menu, as on the property-grid spin editors
                for evt in (wx.EVT_CONTEXT_MENU, wx.EVT_RIGHT_DOWN):
                    spin.Bind(evt, lambda e, s=spin, k=c: self.on_step_menu(s, key=k))
                rng = wx.SpinCtrlDouble(self.const_panel, min=1e-6, max=1e9, inc=0.5, initial=2, size=(75, -1))
                rng.SetDigits(3)
                rng.SetToolTip("The slider spans −range … +range")
                rng.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e, k=c: self.on_const(k, 'range'))
                theme.bind_spin_enter(rng, lambda k=c: self.on_const(k, 'range'))
                for w, flag in ((slider, wx.EXPAND), (spin, 0), (rng, 0)):
                    grid.Add(w, 0, flag | wx.ALIGN_CENTER_VERTICAL)
                self.const_widgets[c] = (slider, spin, rng)
            self.const_panel.SetSizer(grid)
            theme.apply_theme_to_window(self.const_panel, theme.get_theme())
            self.math_panel.Layout()
            self.const_key = (id(dd), tuple(consts))
        for c, (slider, spin, rng) in self.const_widgets.items():
            v, R = mm['params'].get(c, 1.0), mm['ranges'].get(c, 2.0)
            inc = mm.setdefault('steps', {}).get(c, 0.01)      # sessions saved before 'steps' existed
            rng.SetValue(R)
            spin.SetDigits(max(4, int(np.ceil(-np.log10(inc)))))   # enough decimals for the step to show
            spin.SetIncrement(inc)
            spin.SetValue(v)
            slider.SetValue(int(round(np.clip((v + R)/(2*R), 0, 1)*1000)))

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
        self.parent.Data[int(ii)][key] = event.GetValue()
        if key == 'freq':
            self.update_controls()    # the BG ratio labels follow the frequency
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
        mp.norm = ('none', 'max', 'DI')[self.seg_norm.GetSelection()]
        offset = self.spin_offset.GetValue()
        if event is None and offset == mp.offset:   # Enter on an unchanged value
            return
        mp.offset = offset
        mp.update_graph()

    def on_copyDI(self, event):
        mf = self.parent
        lines = ["#\ttitle\tMW freq [GHz]\tDI"]
        for ii, dd in enumerate(mf.Data):
            lines.append(f"{ii + 1}\t{dd['title']}\t{dd['freq']:.6f}\t{mf.process(dd)['DI']:.8g}")
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(wx.TextDataObject("\n".join(lines)))
            wx.TheClipboard.Close()
            mf.statusbar.SetStatusText(f"DI table of {len(mf.Data)} spectra copied to the clipboard", 0)

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
        dd = mf.Data[mf.current]
        b = dd['bglist'][int(jj)]
        b[key] = event.GetProperty().GetValueAsString() if key == 'bgdata' else float(event.GetValue())
        # no grid rebuild from inside its own change event -- only relabel the category
        self.bgtree.SetPropertyLabel(f"cat {jj}", f"BG {int(jj) + 1}:  B·{mf.bgratio(dd, b)[0]:.6f}")
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
        mf.Data[mf.current]['bglist'].append({'bgdata': mf.BGData[sel]['title'], 'scale': 1.0,
                                              'fratio': -1.0, 'shift': 0.0})
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

    # ---------------- Phase tab handlers
    def on_phase_slider(self, key):
        mf = self.parent
        if mf.current < 0:
            return
        mf.Data[mf.current][key] = self.slider_phase[key].GetValue()/10
        self.spin_phase[key].SetValue(mf.Data[mf.current][key])
        mf.matplotlib_panel.request_update()

    def on_phase_spin(self, key, reset=False):
        """key 'lockin' or 'phase'; reset puts the lock-in phase back to the file's value."""
        mf = self.parent
        if mf.current < 0:
            return
        dd = mf.Data[mf.current]
        if reset:
            self.spin_phase[key].SetValue(dd['lockin0'])
        if self.spin_phase[key].GetValue() == dd[key]:
            return
        dd[key] = self.spin_phase[key].GetValue()
        self.slider_phase[key].SetValue(int(round(np.clip(dd[key], -180, 180)*10)))
        mf.matplotlib_panel.update_graph()

    # ---------------- Integrate tab handlers
    def on_range(self, event=None):
        mf = self.parent
        if mf.current < 0:
            return
        new = [self.spin_from.GetValue(), self.spin_to.GetValue()]
        if new == mf.Data[mf.current]['intrange']:
            return
        mf.Data[mf.current]['intrange'] = new
        mf.matplotlib_panel.update_graph()

    def on_order(self, key):
        mf = self.parent
        if mf.current < 0 or self.spin_order[key].GetValue() == mf.Data[mf.current][key + 'order']:
            return
        mf.Data[mf.current][key + 'order'] = self.spin_order[key].GetValue()
        mf.matplotlib_panel.update_graph()

    def on_region(self, key, action, j=None):
        """Fit ranges of baseline key ('bl1'/'bl2'/'bl3'): 'edit' range j from its spins, 'delete' range j,
        'add' one after the last (5 % of the field span; mid-spectrum if there is none), 'clear' all."""
        mf = self.parent
        if mf.current < 0:
            return
        dd = mf.Data[mf.current]
        regs = dd[key + 'regions']
        if action == 'edit':
            if j >= len(regs):              # the rows are about to be rebuilt
                return
            s_lo, s_hi = self.region_rows[key][2][j]
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

    # ---------------- Math tab handlers
    def math_check(self, expr):
        """(constant names, refs {sN: source title}) of an expression typed with the current Data numbering.
        Raises SyntaxError / ValueError with a readable message."""
        mf = self.parent
        _, specs, consts = mf.parse_expr(expr)
        refs = {}
        for s in specs:
            n = int(s[1:])
            if not 1 <= n <= len(mf.Data):
                raise ValueError(f"there is no spectrum {n}")
            if mf.Data[n - 1]['math'] is not None:
                raise ValueError(f"{s} is a math spectrum -- use its sources instead")
            refs[s] = mf.Data[n - 1]['title']
        return consts, refs

    def on_expr(self, event=None, create=False):
        """Create (create=True) adds a new math spectrum. Update (button or Enter) changes the current math
        spectrum, else the last one in the Data list, and creates one only when there is none yet."""
        mf = self.parent
        expr = self.txt_expr.GetValue().strip()
        try:
            consts, refs = self.math_check(expr)
        except (SyntaxError, ValueError) as e:
            wx.MessageBox(f"Cannot use '{expr}':\n{e}", "Math", wx.OK | wx.ICON_WARNING)
            return
        targets = [ii for ii, d in enumerate(mf.Data) if d['math'] is not None]
        if create or not targets:
            # the first spectrum in the expression gives the field axis and the MW frequency
            first = next(d for d in mf.Data if d['title'] == next(iter(refs.values())))
            dd = mf.newDataset({'x': first['ax']['x'].copy(), 'xlabel': 'B [mT]', 'freq1': first['freq'],
                                'title': '= ' + expr}, np.zeros(first['ax']['x'].size), {'expression': expr}, '')
            dd['fname'] = 'math'
            dd['math'] = {'expr': expr, 'refs': refs, 'params': {c: 1.0 for c in consts},
                          'ranges': {c: 2.0 for c in consts}, 'steps': {c: 0.01 for c in consts},
                          'align': ('field', 'g')[self.seg_align.GetSelection()],
                          'stage': ('raw', 'processed')[self.seg_stage.GetSelection()]}
            mf.addDataset(dd)
            return
        target = mf.current if mf.current in targets else targets[-1]
        dd = mf.Data[target]
        mm = dd['math']
        mm['expr'], mm['refs'] = expr, refs
        mm['align'] = ('field', 'g')[self.seg_align.GetSelection()]
        mm['stage'] = ('raw', 'processed')[self.seg_stage.GetSelection()]
        for c in consts:
            mm['params'].setdefault(c, 1.0)
            mm['ranges'].setdefault(c, 2.0)
        dd['title'] = '= ' + expr
        mf.current = target             # its constant sliders show up
        self.update_filelist()
        mf.matplotlib_panel.update_graph()

    def on_math_mode(self, event=None):
        mf = self.parent
        if mf.current < 0 or mf.Data[mf.current]['math'] is None:
            return                 # used by the next Create / Update
        mm = mf.Data[mf.current]['math']
        mm['align'] = ('field', 'g')[self.seg_align.GetSelection()]
        mm['stage'] = ('raw', 'processed')[self.seg_stage.GetSelection()]
        mf.matplotlib_panel.update_graph()

    def on_math_insert(self, token):
        """Quick-insert button: operators padded with spaces, functions as 'f()' with the cursor inside."""
        text = {'abs': 'abs()', 'diff': 'diff()'}.get(token, f" {token} " if token in '+-/*' else token)
        self.txt_expr.WriteText(text)
        if text.endswith('()'):
            self.txt_expr.SetInsertionPoint(self.txt_expr.GetInsertionPoint() - 1)
        self.txt_expr.SetFocus()

    def on_math_fill(self, kind):
        """Average / weighted sum of every loaded spectrum (math spectra cannot be sources) into the expression box."""
        mf = self.parent
        ids = [ii + 1 for ii, d in enumerate(mf.Data) if d['math'] is None]
        if not ids:
            return
        expr = " + ".join(f"c{i}*s{i}" if kind == 'weighted' else f"s{i}" for i in ids)
        self.txt_expr.SetValue(f"({expr})/{len(ids)}" if kind == 'mean' else expr)

    def on_step_menu(self, spin, steps=(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 1e4, 1e5, 1e6),
                     digits=4, key=None):
        """Right-click on a spin box: choose its arrow step (as on the property-grid spins). Shows at least
        `digits` decimals, more if the step needs them. key: a math constant whose step the spectrum keeps."""
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
        mf = self.parent
        if key is not None and mf.current >= 0 and mf.Data[mf.current]['math'] is not None:
            mf.Data[mf.current]['math'].setdefault('steps', {})[key] = ids[choice]
        v = spin.GetValue()
        spin.SetDigits(max(digits, int(np.ceil(-np.log10(ids[choice])))))
        spin.SetIncrement(ids[choice])
        spin.SetValue(v)

    def on_const(self, c, source):
        """A constant's slider ('slider'), value box ('value') or slider half-range ('range') changed."""
        mf = self.parent
        if mf.current < 0 or mf.Data[mf.current]['math'] is None or c not in self.const_widgets:
            return
        mm = mf.Data[mf.current]['math']
        slider, spin, rng = self.const_widgets[c]
        R = max(rng.GetValue(), 1e-12)
        if source == 'slider':
            v = -R + 2*R*slider.GetValue()/1000
            spin.SetValue(v)
        elif source == 'value':
            v = spin.GetValue()
        else:                      # new range: re-place the slider, value unchanged
            mm['ranges'][c] = R
            v = mm['params'].get(c, 1.0)
        if source != 'slider':
            slider.SetValue(int(round(np.clip((v + R)/(2*R), 0, 1)*1000)))
        if v == mm['params'].get(c):
            return
        mm['params'][c] = v
        mf.matplotlib_panel.request_update()


# --------------------------------------------------------------------------- #
# Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Data = []
        self.BGData = []
        self.current = -1        # index in self.Data edited by the Background/Phase/Integrate tabs
        self.proc_cache = {}     # id(dataset) -> process() result, filled only during one update_graph pass
        self.SetTitle(f"pyCWEPR {VERSION}")
        self.SetSize((1200, 720))
        self.theme_name = "light"

        # ---- Top icon toolbar
        self.top_toolbar = wx.Panel(self)
        toolbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.tb_change_folder = theme.IconButton(self.top_toolbar, theme.draw_folder_icon, tooltip="Change folder…")
        self.tb_up = theme.IconButton(self.top_toolbar, theme.draw_up_icon, tooltip="Up")
        self.tb_save_session = theme.IconButton(self.top_toolbar, theme.draw_save_icon, tooltip="Save Session")
        self.tb_grid = theme.IconToggleButton(self.top_toolbar, theme.draw_grid_icon, tooltip="Axes grid")
        self.btn_theme_icon = theme.IconButton(self.top_toolbar, theme.draw_theme_icon, tooltip="Toggle dark mode (Ctrl+D)")
        # plot-row toggles, each icon the row's line shape: derivative, absorption peak, double-integral step
        self.tb_rows = {
            'showDeriv': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: -t*np.exp(-4*t*t)/0.2144), tooltip="Derivative row", value=True),
            'showAbs': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: 2*np.exp(-8*t*t) - 1), tooltip="Absorption row", value=True),
            'showDI': theme.IconToggleButton(self.top_toolbar, theme.draw_curve_icon(
                lambda t: np.tanh(3*t)), tooltip="Double-integral row", value=True)}
        # x axis: off = field (B), on = g; the glyph shows the axis in use
        self.tb_xaxis = theme.IconToggleButton(
            self.top_toolbar, lambda gc, w, h, color, on: theme.draw_text_icon("g" if on else "B")(gc, w, h, color),
            tooltip="X axis: field (B) or g\n(g uses each spectrum's own MW frequency)")
        for w in (self.tb_change_folder, self.tb_up, None, self.tb_save_session, None,
                  *self.tb_rows.values(), 'stretch', None, self.tb_xaxis, self.tb_grid, self.btn_theme_icon):
            if w is None:
                toolbar_sizer.Add(wx.StaticLine(self.top_toolbar, style=wx.LI_VERTICAL), 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
            elif isinstance(w, str):
                toolbar_sizer.AddStretchSpacer(1)
            else:
                toolbar_sizer.Add(w, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        self.top_toolbar.SetSizer(toolbar_sizer)
        self.tb_grid.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_grid)
        for btn in self.tb_rows.values():
            btn.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_rows)
        self.tb_xaxis.Bind(wx.EVT_TOGGLEBUTTON, lambda e: (
            setattr(self.matplotlib_panel, 'xaxis', 'g' if self.tb_xaxis.GetValue() else 'field'),
            self.matplotlib_panel.update_graph()))
        self.btn_theme_icon.Bind(wx.EVT_BUTTON, self.on_toggle_theme)
        theme_id = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_toggle_theme, id=theme_id)
        self.SetAcceleratorTable(wx.AcceleratorTable([(wx.ACCEL_CMD, ord('D'), theme_id)]))

        # ---- Status strip: message | cursor readout
        self.statusbar = theme.StatusStrip(self, widths=[-1, 320])
        self.statusbar.SetStatusText("Ready", 0)

        # ---- file browser | (plot | tabs), as in pyHYSCORE
        self.splitter_main = wx.SplitterWindow(self, style=wx.SP_3D)
        self.splitter_main.SetMinimumPaneSize(20)
        self.file_browser = FileBrowserPanel(self.splitter_main, self)
        self.right_panel = wx.Panel(self.splitter_main)
        self.splitter_right = wx.SplitterWindow(self.right_panel, style=wx.SP_3D)
        self.matplotlib_panel = MatplotlibPanel(self.splitter_right, self)
        self.tabulated_panel = TabulatedPanel(self.splitter_right, self)
        self.splitter_right.SplitVertically(self.matplotlib_panel, self.tabulated_panel)
        self.splitter_right.SetSashGravity(1)
        self.splitter_right.SetSashPosition(-360)    # room for five tabs and the constant sliders
        self.right_sash = -360
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
        self.matplotlib_panel.update_graph()

    # -----------------------------------------------------------------------
    def on_toggle_theme(self, event):
        self.theme_name = "dark" if self.theme_name == "light" else "light"
        self.apply_theme()

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

    def on_toggle_rows(self, event):
        """Derivative / Absorption / Double-integral row toggles; at least one row stays on."""
        if not any(btn.GetValue() for btn in self.tb_rows.values()):
            event.GetEventObject().SetValue(True)
        for name, btn in self.tb_rows.items():
            setattr(self.matplotlib_panel, name, btn.GetValue())
        self.matplotlib_panel.update_graph()

    def apply_theme(self):
        theme.set_theme(self.theme_name)
        theme.apply_theme_to_window(self, self.theme_name)
        self.tabulated_panel.filelist.apply_theme()
        theme.theme_propgrid(self.tabulated_panel.bgtree, caption_key="pill_load")
        theme.style_figure(self.matplotlib_panel.figure, self.matplotlib_panel.figure.axes)
        self.matplotlib_panel.canvas.draw_idle()
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
        self.matplotlib_panel.layout = None
        self.matplotlib_panel.update_graph()
        fb = self.file_browser
        selected, idx = [], fb.list.GetFirstSelected()
        while idx != -1:
            selected.append(fb.items[fb.list.GetItemData(idx)][0])
            idx = fb.list.GetNextSelected(idx)
        fb.refresh_file_list(select_name=selected)

    # -----------------------------------------------------------------------
    # Data handling
    # -----------------------------------------------------------------------
    def loadData(self, path, BG=False):
        ext = os.path.splitext(path)[1].lower()
        try:
            if ext == '.xml':
                self.loadXML(path, BG)
                return
            if ext == '.csv':
                ax, data, dsc = csvread(path)
                if ax is None:      # not a Magnettech result: plain columns (pyCWEPR export, other programs)
                    ax, data, dsc = self.readASCII(path)
            elif ext in ('.txt', '.dat', '.asc'):
                ax, data, dsc = self.readASCII(path)
            elif ext in ('.dsc', '.dta'):
                ax, data, dsc = Mybr.brukerread(path, return_ax=True, return_dsc=True)
            else:
                raise ValueError(f"unsupported file type '{ext}'")
            self.addDataset(self.newDataset(ax, data, dsc, path), BG)
        except Exception as e:
            traceback.print_exc()
            wx.MessageBox(f"Loading '{os.path.basename(path)}' failed:\n\n{type(e).__name__}: {e}",
                          "Cannot load", wx.OK | wx.ICON_ERROR)

    def readASCII(self, path):
        """(ax, data, dsc) from plain delimited text: a pyCWEPR export or any numeric CSV / TXT.
        '#' lines are metadata ('# key: value' -> dsc), an optional row of column names, then numbers
        separated by ',' ';' tab or spaces (with ';' or tab, a decimal comma is accepted).
        First column = field (mT, or G if its name says so); the spectrum is the 'derivative'
        column if there is one, otherwise the second column."""
        dsc, names, rows = {}, None, []
        with open(path, 'r', encoding='utf-8-sig', errors='replace') as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith('#'):
                    key, _, value = s[1:].partition(':')
                    if value:
                        dsc[key.strip()] = value.strip()
                    continue
                sep = ';' if ';' in s else '\t' if '\t' in s else ',' if ',' in s else None
                parts = [q.strip() for q in s.split(sep)]
                try:
                    rows.append([float(q.replace(',', '.')) if q else np.nan for q in parts])
                except ValueError:
                    if rows:
                        raise ValueError(f"non-numeric line after the data started: '{s[:60]}'")
                    names = parts           # column names (the last text row before the numbers)
        if not rows or min(len(r) for r in rows) < 2:
            raise ValueError("no table of at least two numeric columns")
        ncol = min(len(r) for r in rows)
        table = np.array([r[:ncol] for r in rows])
        names = names[:ncol] if names is not None and len(names) >= ncol else ['B [mT]'] + [f"column {k + 1}" for k in range(1, ncol)]
        col = next((k for k in range(1, ncol) if 'deriv' in names[k].lower()), 1)
        good = np.isfinite(table[:, 0]) & np.isfinite(table[:, col])
        try:
            freq = float(dsc.get('MW frequency [GHz]', 0))
        except ValueError:
            freq = 0.0
        ax = {'x': table[good, 0], 'xlabel': names[0], 'ylabel': names[col], 'freq1': freq,
              'title': dsc.get('title') or os.path.splitext(os.path.basename(path))[0]}
        return ax, table[good, col], dsc

    def newDataset(self, ax, data, dsc, path):
        """Dataset dict with default processing. Field axis always ascending, in mT; freq in GHz."""
        x = np.asarray(ax['x'], dtype=float).ravel()
        y = np.real(np.squeeze(np.asarray(data))).astype(float)
        if y.ndim != 1 or y.size != x.size:
            raise ValueError(f"not a 1D spectrum (data shape {np.shape(data)}, {x.size} field points)")
        quad = ax.get('quad')                # lock-in quadrature (Magnettech XML), same grid as data
        quad = np.asarray(quad, dtype=float).ravel() if quad is not None else None
        if quad is not None and quad.size != x.size:
            quad = None
        if 'G' in re.split(r'[\s,\[\]()]+', str(ax.get('xlabel', ''))):    # Bruker field in Gauss
            x = x/10
        if x[0] > x[-1]:
            x, y = x[::-1].copy(), y[::-1].copy()
            quad = quad[::-1].copy() if quad is not None else None
        freq = float(ax.get('freq1', 0) or 0)
        if freq > 1e6:                       # Bruker MWFQ is in Hz
            freq *= 1e-9
        lo, hi, w = float(x[0]), float(x[-1]), float(x[-1] - x[0])
        edges = [[lo, lo + 0.1*w], [hi - 0.1*w, hi]]
        return {'ax': {'x': x, 'xlabel': 'B, mT', 'ylabel': str(ax.get('ylabel', ''))},
                'data': y, 'dsc': dsc if isinstance(dsc, dict) else {},
                'fname': os.path.basename(path), 'fullpath': path,
                'title': str(ax.get('title') or os.path.basename(path)),
                'freq': freq, 'show': True, 'scale': 1.0,
                'yshift': 0.0,                    # vertical display offset, in tallest-trace units of each plot row
                'color': '',                      # '' = automatic: C<autocolor>
                'autocolor': -1,                  # default colour index, fixed when added (-1 = not yet)
                # [{'bgdata': BG title, 'scale', 'fratio' (-1 = f/f_BG), 'shift' [mT]}]
                'bglist': [], 'phase': 0.0,
                # lock-in phase [deg]: data*cos(lockin - lockin0) + quad*sin(lockin - lockin0); quad None = fixed
                'quad': quad, 'lockin0': float(ax.get('phase0', 0) or 0), 'lockin': float(ax.get('phase0', 0) or 0),
                # math spectrum (Math tab): {'expr', 'refs' {s1: source title}, 'params' {name: value},
                # 'ranges' {name: slider half-range}, 'steps' {name: value-box increment}, 'align' 'field'|'g', 'stage' 'raw'|'processed'}; None = loaded data
                'math': None,
                'bl1order': 0, 'bl1regions': copy.deepcopy(edges),    # derivative baseline
                'bl2order': 1, 'bl2regions': copy.deepcopy(edges),    # absorption baseline
                'bl3order': 0, 'bl3regions': copy.deepcopy(edges[:1]),    # double-integral baseline: left edge only
                'intrange': [lo, hi]}

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
        """Give dd the lowest default colour index no other spectrum holds; it stays with the spectrum,
        so removing or reordering others never recolours it."""
        if dd.get('autocolor', -1) >= 0:
            return
        used = {d.get('autocolor', -1) for d in self.Data if d is not dd}
        dd['autocolor'] = next(k for k in range(len(used) + 1) if k not in used)

    def trace_color(self, dd):
        return dd['color'] or f"C{dd['autocolor'] % 10}"

    def bgratio(self, dd, b):
        """(field-axis ratio, background dict or None) for one entry of dd['bglist']."""
        bd = next((d for d in self.BGData if d['title'] == b['bgdata']), None)
        if b['fratio'] > 0:
            return b['fratio'], bd
        if bd is not None and dd['freq'] > 0 and bd['freq'] > 0:
            return dd['freq']/bd['freq'], bd
        return 1.0, bd

    def polybaseline(self, x, y, regions, order):
        """Polynomial of `order` fitted to the points of y inside `regions` (list of [lo, hi] in x), evaluated on all of x."""
        m = np.zeros(x.shape, dtype=bool)
        for lo, hi in regions:
            m |= (x >= min(lo, hi)) & (x <= max(lo, hi))
        if order < 0 or m.sum() <= order:
            return np.zeros_like(y)
        xs = (x - x.mean())/np.ptp(x)        # centred/scaled: keeps high orders well conditioned
        return np.polynomial.polynomial.polyval(xs, np.polynomial.polynomial.polyfit(xs[m], y[m], order))

    def parse_expr(self, expr):
        """Checked math expression -> (code, spectrum names, constant names), names in order of appearance.
        Allowed: numbers, + - * / **, s1, s2 ... (spectra by Data-list number), calls of
        abs/max/min/mean/sqrt/diff, and any other name as an adjustable constant. Nothing else gets compiled."""
        tree = ast.parse(expr.strip().replace('−', '-').replace('×', '*').replace('·', '*'), mode='eval')
        funcs = ('abs', 'max', 'min', 'mean', 'sqrt', 'diff')
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if not (isinstance(node.func, ast.Name) and node.func.id in funcs) or node.keywords:
                    raise ValueError(f"only {', '.join(funcs)} can be called")
            elif isinstance(node, ast.Constant):
                if not isinstance(node.value, (int, float)) or isinstance(node.value, bool):
                    raise ValueError(f"'{node.value!r}' is not a number")
            elif not isinstance(node, (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Name, ast.Load,
                                       ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd)):
                raise ValueError(f"'{type(node).__name__}' is not allowed")
        names = sorted({(n.col_offset, n.id) for n in ast.walk(tree)
                        if isinstance(n, ast.Name) and n.id not in funcs})
        specs, consts = [], []
        for _, name in names:
            target = specs if re.fullmatch(r's\d+', name) else consts
            if name not in target:
                target.append(name)
        if not specs:
            raise ValueError("the expression uses no spectrum (s1, s2, ...)")
        return compile(tree, '<math>', 'eval'), specs, consts

    def process(self, dd):
        """[math] -> lock-in phase -> background -> derivative baseline -> phase -> 1st integral -> absorption baseline
        -> 2nd integral -> double-integral baseline. x, raw, bg, y0 (raw-bg), bl1, y (final derivative) are on the
        full field axis; xi, abs, bl2, dint0, bl3, dint (= dint0 - bl3) on the integration range only.
        DI = dint[-1] (data units x mT^2)."""
        x = dd['ax']['x']
        warn = []
        mm = dd['math']
        if mm is not None:
            # Math spectrum: re-evaluated on every pass, so its sources and constants act live.
            # Sources are found by title (robust to reordering); each is spline-interpolated onto x.
            try:
                code, specs, consts = self.parse_expr(mm['expr'])
                ns = {'abs': np.abs, 'max': np.max, 'min': np.min, 'mean': np.mean, 'sqrt': np.sqrt,
                      'diff': lambda v: np.gradient(np.broadcast_to(v, x.shape).astype(float), x)}   # d/dB [per mT]
                ns.update({c: mm['params'].get(c, 1.0) for c in consts})
                for s in specs:
                    src = next((d for d in self.Data if d['title'] == mm['refs'].get(s)), None)
                    if src is None:
                        raise ValueError(f"{s} ('{mm['refs'].get(s, '?')}') is not loaded")
                    if src['math'] is not None:
                        raise ValueError(f"{s} is itself a math spectrum")
                    ratio = dd['freq']/src['freq'] if mm['align'] == 'g' and dd['freq'] > 0 and src['freq'] > 0 else 1.0
                    ps = self.proc_cache.get(id(src)) or self.process(src)
                    v = CubicSpline(src['ax']['x']*ratio, ps['y' if mm['stage'] == 'processed' else 'raw'],
                                    extrapolate=False)(x)
                    if np.isnan(v).any():
                        warn.append(f"'{src['title']}' does not cover the whole field range of '{dd['title']}' (0 outside)")
                    ns[s] = np.nan_to_num(v)
                dd['data'] = np.broadcast_to(np.asarray(eval(code, {'__builtins__': {}}, ns), dtype=float), x.shape).copy()
            except Exception as e:      # keeps the last good result
                warn.append(f"'{dd['title']}': {e}")
        raw = dd['data']
        if dd['quad'] is not None:
            d = np.deg2rad(dd['lockin'] - dd['lockin0'])
            raw = raw*np.cos(d) + dd['quad']*np.sin(d)
        bg = np.zeros_like(x)
        for b in dd['bglist']:
            ratio, bd = self.bgratio(dd, b)
            if bd is None:
                warn.append(f"background '{b['bgdata']}' is not loaded")
                continue
            # cubic spline: the background may have its own range, point count and frequency
            bgi = CubicSpline(bd['ax']['x']*ratio + b['shift'], bd['data'], extrapolate=False)(x)
            if np.isnan(bgi).any():
                warn.append(f"'{bd['title']}' does not cover the whole field range (set to 0 outside)")
            bg += b['scale']*np.nan_to_num(bgi)
        y0 = raw - bg
        bl1 = self.polybaseline(x, y0, dd['bl1regions'], dd['bl1order'])
        # zero-padded to 2N so the FFT-based Hilbert transform does not wrap one end into the other
        y = (np.exp(1j*np.deg2rad(dd['phase']))*hilbert(y0 - bl1, N=2*x.size)[:x.size]).real
        lo, hi = min(dd['intrange']), max(dd['intrange'])
        m = (x >= lo) & (x <= hi)
        if m.sum() < 3:
            warn.append(f"integration range of '{dd['title']}' holds no data -- using all")
            m[:] = True
        xi = x[m]
        absn = cumulative_trapezoid(y[m], xi, initial=0)
        bl2 = self.polybaseline(xi, absn, dd['bl2regions'], dd['bl2order'])
        dint0 = cumulative_trapezoid(absn - bl2, xi, initial=0)
        bl3 = self.polybaseline(xi, dint0, dd['bl3regions'], dd['bl3order'])
        dint = dint0 - bl3
        return {'x': x, 'raw': raw, 'bg': bg, 'y0': y0, 'bl1': bl1, 'y': y, 'xi': xi, 'abs': absn, 'bl2': bl2,
                'dint0': dint0, 'bl3': bl3, 'dint': dint, 'DI': float(dint[-1]), 'warn': warn}

    # -----------------------------------------------------------------------
    # Session XML (same encoding as pyHYSCORE / pyEPR)
    # -----------------------------------------------------------------------
    def dumpXML(self):
        if hasattr(wx, 'OSX_FILEDIALOG_ALWAYS_SHOW_TYPES'):
            wx.SystemOptions.SetOption(wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES, 1)
        with wx.FileDialog(self, "Save Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultDir=self.file_browser.path, style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            pathname = dlg.GetPath()
        mp = self.matplotlib_panel
        data_dict = {'program': 'pyCWEPR', 'version': VERSION, 'Data': self.Data, 'BGData': self.BGData,
                     'view': {'norm': mp.norm, 'xaxis': mp.xaxis, 'offset': float(mp.offset),
                              'showRawBG': mp.showRawBG, 'showBaselines': mp.showBaselines,
                              'showDeriv': mp.showDeriv, 'showAbs': mp.showAbs, 'showDI': mp.showDI,
                              'showGrid': mp.showGrid,
                              'xfix': mp.xfix, 'yfix': mp.yfix, 'current': self.current}}

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

    def exportASCII(self):
        """Chosen traces of the current spectrum as comma-separated text: '# key: value' metadata lines,
        one row of column names, then B [mT] and the traces. Absorption and double integral exist only
        inside the integration range; their cells are empty outside it."""
        if self.current < 0:
            wx.MessageBox("Load a spectrum first.", "Export", wx.OK | wx.ICON_INFORMATION)
            return
        dd = self.Data[self.current]
        p = self.process(dd)
        inside = np.isin(p['x'], p['xi'])
        full = lambda v: np.where(inside, np.interp(p['x'], p['xi'], v), np.nan)    # range-only -> full axis
        traces = [("raw", p['raw']), ("background", p['bg']), ("derivative", p['y']),
                  ("absorption", full(p['abs'] - p['bl2'])), ("double integral", full(p['dint']))]
        with wx.MultiChoiceDialog(self, "Traces to save (field axis B [mT] is always the first column):",
                                  "Export ASCII", [name for name, _ in traces]) as dlg:
            dlg.SetSelections([2, 3, 4])
            if dlg.ShowModal() != wx.ID_OK or not dlg.GetSelections():
                return
            chosen = [traces[i] for i in dlg.GetSelections()]
        with wx.FileDialog(self, "Export ASCII", defaultDir=self.file_browser.path,
                           defaultFile=re.sub(r'[\\/:*?"<>|]+', '_', dd['title']) + '_pyCWEPR.csv',
                           wildcard="CSV (*.csv)|*.csv|All files (*.*)|*.*",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as dlg:
            if dlg.ShowModal() == wx.ID_CANCEL:
                return
            pathname = dlg.GetPath()
        meta = [f"pyCWEPR {VERSION} export of: {dd['title']}", f"title: {dd['title']}", f"source: {dd['fullpath']}",
                f"MW frequency [GHz]: {dd['freq']:.9g}"]
        if dd['quad'] is not None:
            meta.append(f"lock-in phase [deg]: {dd['lockin']:g} (file {dd['lockin0']:g})")
        meta += [f"background {jj + 1}: {b['bgdata']} x {b['scale']:g}, ratio {self.bgratio(dd, b)[0]:.6f}, "
                 f"shift {b['shift']:g} mT" for jj, b in enumerate(dd['bglist'])]
        meta += [f"phase [deg]: {dd['phase']:g}",
                 f"derivative baseline: order {dd['bl1order']}, regions [mT] "
                 + "; ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in dd['bl1regions']),
                 f"absorption baseline: order {dd['bl2order']}, regions [mT] "
                 + "; ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in dd['bl2regions']),
                 f"double-integral baseline: order {dd['bl3order']}, regions [mT] "
                 + "; ".join(f"{lo:.3f}-{hi:.3f}" for lo, hi in dd['bl3regions']),
                 f"integration range [mT]: {min(dd['intrange']):.3f}-{max(dd['intrange']):.3f}",
                 f"DI: {p['DI']:.10g}"]
        try:
            with open(pathname, 'w', encoding='utf-8', newline='') as f:
                f.write("".join(f"# {line}\n" for line in meta))
                f.write(",".join(["B [mT]"] + [name for name, _ in chosen]) + "\n")
                for row in np.column_stack([p['x']] + [v for _, v in chosen]):
                    f.write(",".join("" if np.isnan(v) else f"{v:.10g}" for v in row) + "\n")
            self.statusbar.SetStatusText(f"Exported {', '.join(n for n, _ in chosen)} to {pathname}", 0)
            self.file_browser.refresh_file_list(select_name=os.path.basename(pathname))
        except OSError as e:
            wx.MessageBox(f"Could not write '{pathname}':\n{e}", "Export failed", wx.OK | wx.ICON_ERROR)

    def loadXML(self, filename, BG=False):
        """A Magnettech data XML adds its measurement(s); a pyCWEPR session replaces the workspace;
        any other session XML (pyEPR / pyHYSCORE) only contributes its 1D spectra.
        Spectra are appended as data or as backgrounds."""
        if next(ET.iterparse(filename, events=('start',)))[1].tag == 'ESRXmlFile':
            for ax, data, dsc in xmlread(filename):
                self.addDataset(self.newDataset(ax, data, dsc, filename), BG)
            return
        result = self.load_mixed_dict_xml(filename)
        if result.get('program') == 'pyCWEPR' and not BG:
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
            for dd in self.Data:              # sessions saved before 'autocolor': by position, as they looked
                self.assign_autocolor(dd)
            view = result.get('view') or {}
            mp, tp = self.matplotlib_panel, self.tabulated_panel
            mp.norm = view.get('norm', 'none') if view.get('norm') in ('none', 'max', 'DI') else 'none'
            mp.xaxis = 'g' if view.get('xaxis') == 'g' else 'field'
            mp.offset = float(view.get('offset', 0.0))
            for name, chk in (('showRawBG', mp.chk_RawBG), ('showBaselines', mp.chk_Baselines),
                              *self.tb_rows.items()):
                setattr(mp, name, bool(view.get(name, getattr(mp, name))))
                chk.SetValue(getattr(mp, name))
            mp.showGrid = bool(view.get('showGrid', False))
            mp.xfix = {k: tuple(v) for k, v in (view.get('xfix') or {}).items()}
            mp.yfix = {k: tuple(v) for k, v in (view.get('yfix') or {}).items()}
            self.tb_grid.SetValue(mp.showGrid)
            tp.seg_norm.SetSelection(('none', 'max', 'DI').index(mp.norm))
            self.tb_xaxis.SetValue(mp.xaxis == 'g')
            tp.spin_offset.SetValue(mp.offset)
            cur = view.get('current', 0)
            self.current = cur if isinstance(cur, int) and 0 <= cur < len(self.Data) else len(self.Data) - 1
            self.update_everything()
            self.statusbar.SetStatusText(f"Session {os.path.basename(filename)} loaded", 0)
            return
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
                    array = np.frombuffer(base64.b64decode(numpy_array.find('data').text), dtype=dtype).reshape(shape).copy()
                except Exception as e:
                    print(f"[pyCWEPR] Could not reconstruct array '{key}' (dtype={dtype}, shape={shape}) "
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
