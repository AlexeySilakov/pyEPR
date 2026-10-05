import os
#import sys
import wx
import wx.lib.mixins.listctrl as listmix
# import wx.grid as gridlib
import wx.propgrid as wxpg

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit
from scipy import special as sp
#import re

import theme
from theme import PillButton, PillCheckBox, PillTabBar, ThemedRadioButton

# Matplotlib imports
import matplotlib
#matplotlib.use('WXAgg')                 # Force the WXAgg backend
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
import matplotlib.tri as mtri
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
# from matplotlib import cm
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

import brukerread as Mybr
import mathfunctions as mf
import classPropGridPanel as MypgPanel
from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim
from colormap_editor import (ColormapEditorDialog, stops_to_cmap,
                              colormap_to_stops)
import grid_search
import time

# import wx.lib.agw.customtreectrl as CT
VERSION=0.9

# pip install wxPython numpy matplotlib scipy
"""
TO DO:
    Recognize FFT File
    Make a list of background files loaded somewhere
    Allow to load orisel files (probably directly in the file list)
    Make SpinCtrlDoubleEditor smarter, so it does not show unnecessary digits 
        ... but the number of values that can be typed is unrestricted (something to do in the validator?)
        one suggestion is to change to wx.lib.agw.FloatSpin
    Make Delete Nucs functional
    Change all lists to the classPropGridPanel
        - figure out the EnumThing ... modify the logic of setting up the control from param classes
        Fix MultipleChoice to only allow one entry
        enum_choices = ["Option 1", "Option 2", "Option 3", "Option 4"]
        self.pg.Append(wxpg.EnumProperty("Color", labels=enum_choices, value=0))
    Add to plot the possibility of overlaying data and simdat.
    Limit zero filling based on data size
    Pretty up the PropertyGrid CTRLS
    Try to recognize frequency from title and see if there is an IF ON flag to add 24.5GHz
    
    Make a warning about earsing the progress when loading xml session
    Fix no-simulation when getting the no resonances
    
"""

class DictPopup(wx.Dialog):
    def __init__(self, parent, data_dict):
        super().__init__(parent, title="Dictionary Contents", size=(350, 300), style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # List control (two columns: Key, Value)
        self.list_ctrl = wx.ListCtrl(panel, style=wx.LC_REPORT | wx.BORDER_SUNKEN)
        self.list_ctrl.InsertColumn(0, "Key", width=120)
        self.list_ctrl.InsertColumn(1, "Value", width=200)

        # Populate the list
        for key, value in data_dict.items():
            index = self.list_ctrl.InsertItem(self.list_ctrl.GetItemCount(), str(key))
            self.list_ctrl.SetItem(index, 1, str(value))

        # OK button
        ok_btn = wx.Button(panel, wx.ID_OK, "OK")

        vbox.Add(self.list_ctrl, 1, wx.EXPAND | wx.ALL, 10)
        vbox.Add(ok_btn, 0, wx.ALIGN_CENTER | wx.BOTTOM, 10)

        panel.SetSizer(vbox)

# --------------------------------------------------------------------------- #
# File Browser Panel
# --------------------------------------------------------------------------- #
class FileBrowserPanel(wx.Panel, listmix.ListCtrlAutoWidthMixin,
                       listmix.ColumnSorterMixin):

    INI_FILENAME = os.path.join(os.path.expanduser("~"), ".filebrowser.ini")

    def __init__(self, parent, mainWindow, path=os.getcwd(), *args, **kw):
        super().__init__(parent, *args, **kw)
        self.parent = mainWindow
        
        self.ini_file = self.INI_FILENAME
        self.path = os.path.abspath(path)          # may be overwritten later
        self.filter_choices = ["All", "DSC", "XML"]
        self.current_filter = "DSC"

        
        mainbox = wx.BoxSizer(wx.HORIZONTAL)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        # 1️⃣ Current path (editable; Enter or focus-loss navigates to the typed path)
        hbox_path = wx.BoxSizer(wx.HORIZONTAL)

        self.lbl_path = wx.TextCtrl(self, value="",
                                    style=wx.TE_PROCESS_ENTER | wx.TE_RIGHT | wx.BORDER_NONE)
        self.lbl_path.Bind(wx.EVT_TEXT_ENTER, self.on_path_entered)
        self.lbl_path.Bind(wx.EVT_KILL_FOCUS, self.on_path_entered)

        hbox_path.Add(self.lbl_path, 1,
                      wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0)

        vbox.Add(hbox_path, 0, wx.EXPAND | wx.ALL, 0)

        # Change folder / Up now live as toolbar icons on MainFrame
        # (self.parent.tb_change_folder / tb_up), bound near the end of
        # this constructor once self.on_change_folder/on_up_clicked exist.

        # 3️⃣ Filter combobox
        hbox_filter = wx.BoxSizer(wx.HORIZONTAL)
        lbl = wx.StaticText(self, label="Filter:")
        self.choice_filter = wx.ComboBox(
            self, choices=self.filter_choices, value=self.current_filter,
            style=wx.CB_READONLY|wx.BORDER_NONE)
        hbox_filter.Add(lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        hbox_filter.Add(self.choice_filter, 1)
        vbox.Add(hbox_filter, 0, wx.EXPAND | wx.ALL, 5)

        # 4️⃣ One‑column list (must be created *before* loading a path!)
        self.list = wx.ListCtrl(
            self,
            style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_VRULES | wx.LC_HRULES|wx.BORDER_NONE)
        self.list.InsertColumn(0, "Name", width=400)
        self.list.InsertColumn(1, "Title", width=400)
        vbox.Add(self.list, 1,
                 wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)

        # Optional: folder icon
        img_list = wx.ImageList(16, 16)
        file_bitmap = wx.Bitmap('./icons/New document.ico', wx.BITMAP_TYPE_ICO)
        xml_bitmap = wx.Bitmap('./icons/Script.ico', wx.BITMAP_TYPE_ICO)
        folder_bitmap = wx.Bitmap('./icons/Folder.ico', wx.BITMAP_TYPE_ICO)
        dsc_bitmap = wx.Bitmap('./icons/List.ico', wx.BITMAP_TYPE_ICO)
        loaded_bitmap = wx.Bitmap('./icons/Green bookmark.ico', wx.BITMAP_TYPE_ICO)
        loadedBG_bitmap = wx.Bitmap('./icons/Equipment.ico', wx.BITMAP_TYPE_ICO)
        # folder_bitmap = wx.ArtProvider.GetBitmap(wx.ART_FOLDER,
        #                                         wx.ART_OTHER, (16, 16))
        # file_bitmap = wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE,
        #                                         wx.ART_OTHER, (16, 16))
        
        self.folder_img_id = img_list.Add(folder_bitmap)
        self.xml_img_id = img_list.Add(xml_bitmap)
        self.file_img_id = img_list.Add(file_bitmap)
        self.loaded_img_id = img_list.Add(loaded_bitmap)
        self.loadedBG_img_id = img_list.Add(loadedBG_bitmap)
        self.dsc_img_id = img_list.Add(dsc_bitmap)

        # Sort indicator arrows for the column headers. ColumnSorterMixin
        # will call SetImage(...) on the sorted column's header with one of
        # these IDs; if GetSortImages() returned (None, None) as it did
        # originally, newer wx raises TypeError. Providing real image IDs
        # both fixes that crash and gives the user a visible sort indicator.
        up_arrow   = wx.ArtProvider.GetBitmap(wx.ART_GO_UP,   wx.ART_OTHER, (16, 16))
        down_arrow = wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN, wx.ART_OTHER, (16, 16))
        self.sort_up_img_id   = img_list.Add(up_arrow)
        self.sort_down_img_id = img_list.Add(down_arrow)

        self.list.AssignImageList(img_list, wx.IMAGE_LIST_SMALL)

        # 5️⃣ Load button
        gbox_Load = wx.BoxSizer(wx.VERTICAL)
        self.btn_load   = PillButton(self, "Load", color_key="pill_load",
                                     draw_badge=theme.draw_squiggle_badge)
        self.btn_load.Bind(wx.EVT_BUTTON, self.on_load_clicked)
        gbox_Load.Add(self.btn_load, 0, wx.EXPAND | wx.ALL, 2)

        self.btn_loadbg = PillButton(self, "+ BG", color_key="pill_load",
                                     draw_badge=theme.draw_flatline_badge)
        self.btn_loadbg.Bind(wx.EVT_BUTTON, self.on_loadBG_clicked)
        gbox_Load.Add(self.btn_loadbg, 0, wx.EXPAND | wx.ALL, 2)

        # Save Session now lives as a toolbar icon on MainFrame
        # (self.parent.tb_save_session), bound below once
        # self.on_save_session exists.

        vbox.Add(gbox_Load, 0, wx.EXPAND | wx.ALL, 0)

        mainbox.Add(vbox, 1, wx.EXPAND | wx.ALL, 0)

        # The collapse toggle used to live here, but a button *inside* this
        # panel would disappear along with the rest of it once the panel is
        # actually hidden (see MainFrame.on_toggle_left) -- it now lives in
        # MainFrame itself, to the left of the splitter, so it stays visible.

        self.SetSizer(mainbox)

        # ------------------------------------------------------------------
        # 3️⃣ Sorting helpers / data storage
        # ------------------------------------------------------------------
        # The wx.ListCtrl has two columns (Name, Title), so pass 2 here.
        # Passing 1 caused IndexError in ColumnSorterMixin.__OnColClick when
        # the user clicked the "Title" header (col index 1 out of range).
        listmix.ColumnSorterMixin.__init__(self, 2)
        self.items = []          # (name, full_path, is_dir)

        # ------------------------------------------------------------------
        # 4️⃣ Event bindings (now that all widgets exist)
        # ------------------------------------------------------------------
        self.choice_filter.Bind(wx.EVT_COMBOBOX, self.on_filter_changed)

        # Change folder / Up / Save Session are toolbar icons on MainFrame
        # (built before FileBrowserPanel -- see MainFrame.__init__), bound
        # here now that this panel's handler methods exist.
        self.parent.tb_change_folder.Bind(wx.EVT_BUTTON, self.on_change_folder)
        self.parent.tb_up.Bind(wx.EVT_BUTTON, self.on_up_clicked)
        self.parent.tb_save_session.Bind(wx.EVT_BUTTON, self.on_save_session)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_item_activated)
        
        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy)   # persistence

        self._load_path_from_file()          # this may call set_path()

        # ------------------------------------------------------------------
        # 6️⃣ Show the (possibly overwritten) current path
        # ------------------------------------------------------------------
        self.lbl_path.SetValue(self.path)
        self.lbl_path.SetInsertionPointEnd()
        self.refresh_file_list()
        
    def on_save_session(self, event):
        self.parent.dumpXML()
        
    def on_load_clicked(self, event):
        idx = self.list.GetFirstSelected()
        if idx == -1:
            wx.MessageBox("Please select a file to load.", "No file selected",
                          wx.OK | wx.ICON_WARNING)
            return
        # After a header-click sort the visual row index no longer matches
        # self.items's insertion order; use the stable data ID attached to
        # the row instead.
        data_id = self.list.GetItemData(idx)
        name, full, is_dir = self.items[data_id]
        if is_dir:
            wx.MessageBox(f"{name} is a directory.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
        elif name.endswith('xml'):
            self.parent.loadXML(full)
        elif name.endswith('DSC'):
            self.parent.loadData(full, BG=False)
        else:
            wx.MessageBox(f"The file {name} does not appear to be to have the right extension.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)  
        # Keep the user on the same row after the list is rebuilt, so
        # successive loads don't jump the cursor back to the top.
        self.refresh_file_list(select_name=name)
            
    # ----------------------------------------------------------------------
    # Persistence helpers (plain‑text .ini file)
    # ----------------------------------------------------------------------
    def _load_path_from_file(self):
        """Read the first line of ~/.filebrowser.ini; if it is a valid
        directory, use it."""
        if not os.path.exists(self.ini_file):
            return
        try:
            with open(self.ini_file, "r", encoding="utf-8") as f:
                stored = f.readline().strip()
        except Exception as e:
            print(f"[pyHYSCORE] Could not read {self.ini_file}: {e}")
            return

        if stored and os.path.isdir(stored):
            # `set_path` will update lbl_path and list for us.
            self.set_path(stored)

    def _save_path_to_file(self):
        """Write the current path to the file."""
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.path + "\n")
        except Exception as e:
            print(f"[pyHYSCORE] Could not save last folder to {self.ini_file}: {e}")

    # ----------------------------------------------------------------------
    # Core listing / filtering logic
    # ----------------------------------------------------------------------
    def get_entries(self):
        try:
            raw = os.listdir(self.path)
        except OSError as e:
            print(f"[pyHYSCORE] Could not list folder '{self.path}': {e}")
            if hasattr(self.parent, 'statusbar'):
                self.parent.statusbar.SetStatusText(f"Could not list folder: {e}", 0)
            raw = []

        dirs, files = [], []
        for name in raw:
            full = os.path.join(self.path, name)
            if os.path.isdir(full):
                dirs.append((name, full, True))
            else:
                files.append((name, full, False))

        if self.current_filter != "All":
            ext = self.current_filter.lower()
            files = [(n, f, d) for n, f, d in files
                     if n.lower().endswith(f".{ext}")]

        dirs.sort(key=lambda x: x[0].lower())
        files.sort(key=lambda x: x[0].lower())
        return dirs + files

    def refresh_file_list(self, select_name=None):
        """
        Rebuild the file list. If `select_name` is given and that filename is
        still present after refresh, re-select it and scroll it into view; this
        lets callers preserve the user's position across a refresh.

        Falls back gracefully: if the name is gone (filter changed, file
        deleted, etc.) the list just renders unselected at the top.

        Also populates the itemDataMap used by ColumnSorterMixin so that
        clicking the "Name" or "Title" column headers actually sorts the
        list. Each row gets a stable data ID (its index in self.items)
        attached via SetItemData, which survives sort reordering so our
        lookup sites can map a visual row back to the underlying item.
        """
        # Remember where the top of the visible region is, so that if we
        # cannot restore a selection we at least keep the scroll position.
        try:
            top_item = self.list.GetTopItem()
        except Exception:
            top_item = -1

        self.list.DeleteAllItems()
        self.items.clear()
        # itemDataMap feeds ColumnSorterMixin: {data_id: (col0_key, col1_key)}
        self.itemDataMap = {}

        for name, full, is_dir in self.get_entries():
            idx = self.list.InsertItem(self.list.GetItemCount(), name)
            if is_dir:
                image_id=self.folder_img_id
            elif name.lower().endswith('.xml'):
                image_id=self.xml_img_id
            elif name.lower().endswith(".dsc"):
                title=''
                title = self.extract_TITL_from_dsc(full)
                image_id=self.dsc_img_id
                for dd in self.parent.Data:
                    if title in dd['title']:
                        image_id=self.loaded_img_id
                for dd in self.parent.BGData:
                    if title in dd['title']:
                        image_id=self.loadedBG_img_id        
                self.list.SetItem(idx, 1, title)

            else:
                image_id=self.file_img_id
                
            self.list.SetItemImage(idx, image_id)
            self.items.append((name, full, is_dir))

            # --- Sorting wiring -------------------------------------------------
            # data_id is the stable identity of this row. Its position in the
            # list may change as the user clicks column headers, but the data
            # value stays with the row, letting us recover (name, full, is_dir)
            # from self.items[data_id] no matter the current visual order.
            data_id = len(self.items) - 1
            self.list.SetItemData(idx, data_id)

            # Build the sort keys. Prefix "0_" for directories and "1_" for
            # files so directories always sort to the top regardless of which
            # column is sorted. Lowercase the rest so sort is case-insensitive.
            dir_prefix = "0_" if is_dir else "1_"
            name_key = dir_prefix + name.lower()
            if name.lower().endswith(".dsc") and not is_dir:
                # title was computed above for .dsc files
                title_key = dir_prefix + (title or "").lower()
            elif is_dir:
                # Make dirs sort by name in the Title column too, so they
                # still cluster at the top when user sorts by Title.
                title_key = dir_prefix + name.lower()
            else:
                title_key = dir_prefix + ""   # non-DSC files have no title
            self.itemDataMap[data_id] = (name_key, title_key)
            # --------------------------------------------------------------------

        self.list.SetColumnWidth(0, -1)

        # Re-apply whatever column sort the user had active before the refresh
        # (ColumnSorterMixin keeps that state in self._col / self._colSortFlag).
        # SortListItems() with no args re-uses the current settings.
        try:
            self.SortListItems()
        except Exception:
            pass

        # Restore selection / scroll position after rebuild.
        restored = False
        if select_name is not None:
            # After sorting, visual row order may differ from self.items order,
            # so walk visible rows and match by filename via the attached data
            # ID rather than by positional index.
            for ridx in range(self.list.GetItemCount()):
                data_id = self.list.GetItemData(ridx)
                if 0 <= data_id < len(self.items) and self.items[data_id][0] == select_name:
                    self.list.Select(ridx, on=1)
                    self.list.Focus(ridx)
                    self.list.EnsureVisible(ridx)
                    restored = True
                    break

        if not restored and top_item >= 0 and self.list.GetItemCount() > 0:
            # Best-effort fallback: scroll back to roughly where the user was.
            last = self.list.GetItemCount() - 1
            self.list.EnsureVisible(min(top_item, last))
        
    def extract_TITL_from_dsc(self, file_path):
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith("TITL"):
                        parts = line.split("'")  # Split by single quotes
                        if len(parts) > 1:
                            return parts[1]  # Return the string between the first two quotes
            return ''  # No matching line found
        except Exception as e:
            # Non-fatal: skip the title for this one file rather than aborting
            # the whole folder listing, but still surface why.
            print(f"[pyHYSCORE] Could not read title from '{file_path}': {e}")
            return ''
    # ----------------------------------------------------------------------
    # Sorting mixin helpers
    # ----------------------------------------------------------------------
    def GetListCtrl(self):  return self.list
    def GetSortImages(self):
        # Returned as (down-arrow-id, up-arrow-id) per ColumnSorterMixin's
        # contract: the down arrow marks a descending sort, the up arrow
        # marks an ascending sort. Returning (None, None) here would crash
        # SetImage in current wxPython.
        return self.sort_down_img_id, self.sort_up_img_id
    def GetColumnSorterData(self): return self.items

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------
    # Event handlers
    # ----------------------------------------------------------------------
    def on_filter_changed(self, event):
        self.current_filter = self.choice_filter.GetStringSelection()
        self.refresh_file_list()

    def on_loadBG_clicked(self, event):
        idx = self.list.GetFirstSelected()
        if idx == -1:
            wx.MessageBox("Please select a file to load.", "No file selected",
                          wx.OK | wx.ICON_WARNING)
            return
        data_id = self.list.GetItemData(idx)
        name, full, is_dir = self.items[data_id]
        if is_dir:
            wx.MessageBox(f"{name} is a directory.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
        else:
            self.parent.loadData(full, BG=True)
        self.refresh_file_list(select_name=name)
        
    def on_change_folder(self, event):
        dlg = wx.DirDialog(self, message="Select a folder",
                           defaultPath=self.path,
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.set_path(dlg.GetPath())
        dlg.Destroy()

    def on_up_clicked(self, event):
        parent = os.path.dirname(self.path)
        if parent and parent != self.path:
            self.set_path(parent)

    def on_item_activated(self, event):
        idx = event.GetIndex()
        data_id = self.list.GetItemData(idx)
        name, full, is_dir = self.items[data_id]
        if is_dir:
            self.set_path(full)

    def on_destroy(self, event):
        self._save_path_to_file()

    # ----------------------------------------------------------------------
    # Path changing helper (updates label, list & persistence)
    # ----------------------------------------------------------------------
    def set_path(self, new_path):
        new_path = os.path.abspath(new_path)
        if os.path.isdir(new_path):
            self.path = new_path
            self.lbl_path.SetValue(self.path)
            self.lbl_path.SetInsertionPointEnd()
            self.refresh_file_list()
            self._save_path_to_file()

    def on_path_entered(self, event):
        """Navigate to a path typed directly into the path TextCtrl."""
        typed = self.lbl_path.GetValue().strip()
        if typed and os.path.isdir(typed):
            self.set_path(typed)
        else:
            # Revert to current valid path on bad input
            self.lbl_path.SetValue(self.path)
            self.lbl_path.SetInsertionPointEnd()
        event.Skip()
# --------------------------------------------------------------------------- #
# Matplotlib Canvas Panel
# --------------------------------------------------------------------------- #
class MatplotlibPanel(wx.Panel):
    """
    A placeholder panel that hosts a matplotlib FigureCanvas.
    """
    def __init__(self, parent, mainWindow, *args, **kw):
        super().__init__(parent, *args, **kw)
        self.parent = mainWindow
        self.showFFT = True
        self.showSim = True
        self.showBG = False
        self.showSkyline = False
        self.showOriSel = False
        self.showDiagonalProj = False
        self.OverlaySim = False
        self.quadrant = 'all'   # 'all' | 'horizontal' | 'vertical'
        # Display-only widening of the f.min/f.max window: how many MHz of
        # context to draw around it. f.min/f.max keep their meaning as the
        # *analysis* window -- they are still the grid search's ROI and still
        # set the simulated frequency range -- and every number on this figure
        # (colour scale, trace normalization, noise floor) is still computed
        # from that window alone, so turning this up only ever adds
        # surroundings; it never moves a curve that was already on screen.
        # Session-local, like `quadrant` above and `noiseCompare` below.
        self.viewPad = 0.0
        self.showFitBox = True
        # How the projection panels reconcile a noiseless simulation with an
        # experiment that can never fall below its own noise floor. Set from
        # the right-click menu on any projection panel (on_projection_menu),
        # and deliberately session-local like `quadrant` above rather than a
        # saved setting.
        self.showNoiseFloor = True
        self.noiseCompare = 'none'   # 'none' (default) | 'lift' | 'subtract'
        # 2D map look; saved with the session, mirrored by toolbar toggles and the Opt "Contours" box
        self.mapMode = 'density'     # 'density' | 'contour'
        self.contourLevels = 10
        self.contourWidth = 0.8
        self.overlayColor = '#000000'
        self.showGrid = False
        self.projection_axes = []    # skyline/diagonal axes, for the hit test
        # One-shot: the unequal-step warning is raised from inside a redraw,
        # so without this it would reappear on every single one.
        self._warned_dxdy = False
        self.rebuildAxes = False  # set when the subplots must be recreated
                                  # even though their number hasn't changed
        self.figure = Figure(figsize=(2, 2), dpi=100)
        self.axes = []
        self.contours = []   # parallel to self.axes; holds QuadContourSet or None
        self.mapContours = []   # parallel to self.axes; colormap contours replacing the image in contour mode
        #self.axes.append(self.figure.add_subplot(111))
        #self.axes.set_title("Placeholder Figure")
        #self.axes.plot([0, 1, 2], [0, 1, 4], marker='o')

        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)

        # Track mouse motion over the figure so we can show data coordinates
        # in the main frame's status bar whenever the cursor is over an axis.
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('figure_leave_event', self.on_figure_leave)
        self.canvas.mpl_connect('button_press_event', self.on_projection_menu)
        
        self.chk_FFT = PillCheckBox(self, "FFT(y)")
        self.chk_FFT.SetValue(self.showFFT)

        self.chk_BG = PillCheckBox(self, "Background")
        self.chk_BG.SetValue(self.showBG)

        self.chk_Skyline = PillCheckBox(self, "Skyline")
        self.chk_Skyline.SetValue(self.showSkyline)

        self.chk_OriSel = PillCheckBox(self, "Orientat maps")
        self.chk_OriSel.SetValue(self.showOriSel)

        self.chk_DiagProj = PillCheckBox(self, "Diag.Projection")
        self.chk_DiagProj.SetValue(self.showDiagonalProj)

        self.chk_OverlaySim = PillCheckBox(self, "Overlay Sim")
        self.chk_OverlaySim.SetValue(self.OverlaySim)

        self.chk_FFT.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_BG.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_Skyline.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_OriSel.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_DiagProj.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_OverlaySim.Bind(wx.EVT_CHECKBOX, self.on_check)

        # View pad / Fit box / Quadrants controls live on the Opt tab (TabulatedPanel.add_optpanel)
        sizerH = wx.WrapSizer(wx.HORIZONTAL)
        sizerH.Add(self.chk_FFT, 0, wx.EXPAND)
        sizerH.Add(self.chk_BG, 0, wx.EXPAND)
        sizerH.Add(self.chk_Skyline, 0, wx.EXPAND)
        sizerH.Add(self.chk_OriSel, 0, wx.EXPAND)
        sizerH.Add(self.chk_DiagProj, 0, wx.EXPAND)
        sizerH.Add(self.chk_OverlaySim, 0, wx.EXPAND)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerH, 0, wx.EXPAND)
        
        sizer.Add(self.canvas, 1, wx.EXPAND)
        
        self.toolbar = theme.ThemedNavigationToolbar(self.canvas)
        self.toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        # update the axes menu on the toolbar
        self.toolbar.update()

        self.SetSizer(sizer)

        # Debounce the (expensive, full-figure) canvas redraw during a
        # live window resize or splitter sash drag -- see
        # theme.CanvasRedrawDebouncer for why this targets canvas.draw
        # specifically. Bound on this panel (not the canvas) since the
        # canvas's own size only actually changes as a *result* of this
        # panel's sizer being relaid out.
        self._redraw_debouncer = theme.CanvasRedrawDebouncer(self.canvas)
        self.Bind(wx.EVT_SIZE, self._redraw_debouncer.on_size)

    # ----- Status-bar coordinate readout -------------------------------------
    def on_mouse_move(self, event):
        """Show data-space x, y in the main frame's status bar while the
        cursor is hovering over a matplotlib axis."""
        sb = getattr(self.parent, 'statusbar', None)
        if sb is None:
            return
        ax = event.inaxes
        if ax is None or event.xdata is None or event.ydata is None:
            sb.SetStatusText("", 1)
            return
        # Use each axis's own formatter so units/precision match the ticks.
        try:
            xs = ax.format_xdata(event.xdata)
            ys = ax.format_ydata(event.ydata)
        except Exception:
            xs = f"{event.xdata:.4g}"
            ys = f"{event.ydata:.4g}"
        sb.SetStatusText(f"x = {xs},  y = {ys}", 1)

    def on_figure_leave(self, event):
        """Clear the coordinate readout when the cursor leaves the figure."""
        sb = getattr(self.parent, 'statusbar', None)
        if sb is not None:
            sb.SetStatusText("", 1)

    def on_projection_menu(self, event):
        if event.button != 3 or event.inaxes is None:
            return
        if str(getattr(self.toolbar, 'mode', '')):
            return


        menu = wx.Menu()
        modes = {}
        for label, mode in (("No noise correction", 'none'),
                            ("Lift simulation to noise floor", 'lift'),
                            ("Subtract noise from experiment", 'subtract')):
            item = menu.AppendRadioItem(wx.ID_ANY, label)
            item.Check(self.noiseCompare == mode)
            modes[item.GetId()] = mode
        floor_id = wx.ID_NONE
        if event.inaxes in self.projection_axes:
            menu.AppendSeparator()
            item = menu.AppendCheckItem(wx.ID_ANY, "Show noise floor")
            item.Check(self.showNoiseFloor)
            floor_id = item.GetId()

        # GetPopupMenuSelectionFromUser rather than PopupMenu + EVT_MENU: it
        # is modal and hands back the chosen id directly, so no handler has to
        # be bound (and then left behind) for menu items that only live for
        # the duration of one click.
        chosen = self.canvas.GetPopupMenuSelectionFromUser(
            menu, self.canvas.ScreenToClient(wx.GetMousePosition()))
        menu.Destroy()
        if chosen == wx.ID_NONE:
            return
        if chosen == floor_id:
                self.showNoiseFloor = not self.showNoiseFloor
        elif chosen in modes:
            self.noiseCompare = modes[chosen]
        self.update_graph()

    def on_quadrant(self, event):
        sel = self.get_quadrant()
        self.quadrant = ('all', 'horizontal', 'vertical')[sel]
        # The skyline panel shares its x axis with the data panel above only
        # while it projects onto that same axis, and which axis that is
        # depends on the quadrant -- so the subplots have to be rebuilt here,
        # not merely re-limited.
        self.rebuildAxes = True
        self.update_graph()

    def get_quadrant(self):
        tp = self.parent.tabulated_panel
        if tp.button_all4.GetValue():
            return 0
        if tp.button_horiz.GetValue():
            return 1
        return 2

    def on_viewpad(self, event=None):
        """New display pad in MHz. Bound to the arrows, to Enter (via
        theme.bind_spin_enter) and to losing focus, so it commits however the
        user leaves the field -- hence the no-op guard: the same value
        arriving twice must not cost a second full redraw."""
        pad = max(0.0, float(self.parent.tabulated_panel.spin_viewPad.GetValue()))
        if event is not None:
            event.Skip()
        if pad == self.viewPad:
            return
        self.viewPad = pad
        self.update_graph()

    def set_map_mode(self, mode):
        """Density image <-> colormap contours; axes rebuilt so no stale image lingers."""
        if mode == self.mapMode:
            return
        self.mapMode = mode
        self.parent.tb_mapmode.SetValue(mode == 'contour')
        self.parent.tb_mapmode.SetToolTip(
            "Contour map (click for density)" if mode == 'contour'
            else "Density map (click for contours)")
        self.parent.tabulated_panel.sync_contour_controls()
        self.rebuildAxes = True
        self.update_graph()

    def set_contour_style(self, levels=None, width=None, overlayColor=None):
        """Level count / line width / overlay colour from the Opt tab."""
        if levels is not None:
            self.contourLevels = int(levels)
        if width is not None:
            self.contourWidth = float(width)
        if overlayColor is not None:
            self.overlayColor = overlayColor
        self.update_graph()

    def set_show_grid(self, show):
        """Grid on/off touches only the existing axes -- no data redraw."""
        self.showGrid = bool(show)
        self.parent.tb_grid.SetValue(self.showGrid)
        for ax in self.axes:
            if ax is not None and ax.name != '3d':
                ax.grid(self.showGrid)
        self.canvas.draw_idle()

    def on_check(self, event):
        (self.showFFT, self.showBG, self.showSkyline, self.showOriSel, self.showDiagonalProj, self.OverlaySim, self.showFitBox) = (
        self.chk_FFT.GetValue(),
        self.chk_BG.GetValue(),
        self.chk_Skyline.GetValue(),
        self.chk_OriSel.GetValue(),
        self.chk_DiagProj.GetValue(),
        self.chk_OverlaySim.GetValue(),
        self.parent.tabulated_panel.chk_FitBox.GetValue(),
        )
        self.parent.tabulated_panel.sync_contour_controls()
        self.update_graph()
        
    def draw_fitbox(self, ax, corners, view):
        """Outline the region that is actually fitted """
        for p in list(ax.patches):
            if p.get_gid() == 'fitbox':
                p.remove()
        if not self.showFitBox or corners is None:
            return
        xmi, xma, ymi, yma = view
        outer = np.array([(xmi, ymi), (xma, ymi), (xma, yma), (xmi, yma)], dtype=float)
        inner = np.asarray(corners, dtype=float)
        shoelace = float(np.sum(inner[:, 0]*np.roll(inner[:, 1], -1)
                                - inner[:, 1]*np.roll(inner[:, 0], -1)))
        if shoelace > 0:          # same winding as `outer`, which is drawn
            inner = inner[::-1]   # counter-clockwise -- reverse to cut a hole

        verts = np.vstack([outer, outer[:1], inner, inner[:1]])
        codes = ([MplPath.MOVETO] + [MplPath.LINETO]*(len(outer)-1) + [MplPath.CLOSEPOLY] +
                 [MplPath.MOVETO] + [MplPath.LINETO]*(len(inner)-1) + [MplPath.CLOSEPOLY])
        wash = PathPatch(MplPath(verts, codes), facecolor='0.5', edgecolor='none',
                         alpha=0.18, zorder=2)
        wash.set_gid('fitbox')
        ax.add_patch(wash)

        box = plt.Polygon(inner, closed=True, fill=False, edgecolor='w',
                          linestyle='--', linewidth=0.8, zorder=3)
        box.set_gid('fitbox')
        ax.add_patch(box)

    def sync_diagproj(self):
        """Hide the Diag.Projection toggle unless a spectrum is on screen."""
        usable = any(dd['show'] and (dd['isfft'] or (self.showFFT and dd['fftactual']))
                     for dd in self.parent.Data)
        if not usable:
            self.showDiagonalProj = False
        if self.chk_DiagProj.IsShown() != usable:
            self.chk_DiagProj.Show(usable)
            self.Layout()

    def update_graph(self):
        self.sync_diagproj()
        shw = []
        for dd in self.parent.Data:
            shw.append(dd['show'])
        nData = shw.count(True)
        nrows = [self.showBG, self.showSim, self.showSkyline, self.showOriSel,self.showDiagonalProj].count(True)+1

        #print(nrows)
        makenew = False
        if len(self.axes)!=nData*nrows or self.rebuildAxes:
            for ax in self.axes:
                if ax is not None:
                    ax.clear()
                    ax.remove()
            self.axes = [None]*(nData*nrows)
            self.contours = [None]*(nData*nrows)
            self.mapContours = [None]*(nData*nrows)
            self.figure.clf()
            makenew = True
        self.rebuildAxes = False
        # Rebuilt from scratch each draw: the axes themselves are recreated
        # whenever the layout changes, so a stale list would make the
        # right-click menu miss (or, worse, fire on a 2D map).
        self.projection_axes = []
        axcnt= 0
        showLabels = True
        ### first row is always DATA
        for dd in self.parent.Data:
            if not dd['show']: continue
            rowcnt=0
            ########## ----- data ---------------------------------
            if self.showFFT & (not dd['isfft']) & dd['fftactual']:
                data = dd['fftdata']
                if isinstance(data, type(None)):
                    self.RaiseError("No FFT has been performed yet. Nothing to display")
                axEx = [np.min(dd['fftax']['x']), np.max(dd['fftax']['x']),
                        np.min(dd['fftax']['y']), np.max(dd['fftax']['y'])]
                freqDomain = True
                xfi,xfa, yfi, yfa = (dd['fmin'], dd['fmax'], dd['fmin'], dd['fmax'] )
                axlabel='Frequency, MHz'
            else:
                data = np.real(dd['data'][:,:,0])
                axEx = [np.min(dd['ax']['x']), np.max(dd['ax']['x']),
                        np.min(dd['ax']['y']), np.max(dd['ax']['y'])]
                freqDomain = dd['isfft']
                if not dd['isfft']:
                    xfi = np.min(dd['ax']['x'])
                    xfa = np.max(dd['ax']['x'])
                    yfi = np.min(dd['ax']['y'])
                    yfa = np.max(dd['ax']['y'])
                else:
                    xfi,xfa, yfi, yfa = (dd['fmin'], dd['fmax'], dd['fmin'], dd['fmax'] )
                axlabel=r'Time, $\mu$s'

            # -- Quadrant display mode ---------------------------------------
            cut_x = cut_y = False
            if self.quadrant == 'horizontal' and yfi < 0 < yfa:
                # Show full x range (±ν₁) but only positive y (ν₂ ≥ 0)
                yfi, cut_y = 0.0, True
            elif self.quadrant == 'vertical' and xfi < 0 < xfa:
                # Show only positive x (ν₁ ≥ 0) but full y range (±ν₂)
                xfi, cut_x = 0.0, True

            # -- Fit window vs view window -----------------------------------
            # xfi..xfa / yfi..yfa is the *analysis* window (f.min..f.max, as
            # cut by the quadrant)
            pad = self.viewPad
            xmi = xfi if cut_x else min(xfi, max(xfi - pad, axEx[0]))
            xma = max(xfa, min(xfa + pad, axEx[1]))
            ymi = yfi if cut_y else min(yfi, max(yfi - pad, axEx[2]))
            yma = max(yfa, min(yfa + pad, axEx[3]))

            # Index masks for the analysis window.
            x_axis = np.linspace(axEx[0], axEx[1], data.shape[1])
            y_axis = np.linspace(axEx[2], axEx[3], data.shape[0])
            x_fit = (x_axis >= xfi) & (x_axis <= xfa)
            y_fit = (y_axis >= yfi) & (y_axis <= yfa)

            # The tilted analysis box (optHYSCORE.diagonalBox). 
            box = None
            fitpoly = None
            activeBox = None
            simbox = None
            simx_axis = simy_axis = None
            if freqDomain:
                adspan = dd.get('antidiagspan', -1.0)
                box = self.parent.diagonalBox(y_axis, x_axis,
                                                  dd['fmin'], dd['fmax'], adspan)

                if type(dd['simdata']) != type(None) and type(dd['simax']) != type(None):
                    simx_axis = np.linspace(np.min(dd['simax']['x']),
                                            np.max(dd['simax']['x']),
                                            dd['simdata'].shape[1])
                    simy_axis = np.linspace(np.min(dd['simax']['y']),
                                            np.max(dd['simax']['y']),
                                            dd['simdata'].shape[0])
                    simbox = self.parent.diagonalBox(simy_axis, simx_axis,
                                                         dd['fmin'], dd['fmax'], adspan)
                if box is not None:
                    fitpoly = box['corners'][:, ::-1]
                    activeBox = data[box['idx0'], box['idx1']]

                dxa = x_axis[1]-x_axis[0] if x_axis.size > 1 else 0.0
                dya = y_axis[1]-y_axis[0] if y_axis.size > 1 else 0.0
                if not self._warned_dxdy and \
                        abs(dxa-dya) > 1e-6*max(abs(dxa), abs(dya), 1e-30):
                    self._warned_dxdy = True
                    wx.MessageBox(
                        f"This spectrum is sampled at {dxa:.6g} MHz along f1 but "
                        f"{dya:.6g} MHz along f2.\n\n"
                        " df1=df2 is expected",
                        wx.OK | wx.ICON_WARNING)

            # A negative z.max means "scale by the maximum of the fitted area" 
            if dd['zmax'] < 0:
                if activeBox is not None:
                    scale = np.max(activeBox)
                elif x_fit.any() and y_fit.any():
                    scale = np.max(data[np.ix_(y_fit, x_fit)])
            else:
                scale = np.max(data)*dd['zmax']
            
            if freqDomain:
                ma = scale
                mi = scale*dd['zmin']
            else:
                ma = np.max(np.max(data)); mi = np.min(np.min(data))
            if makenew:
                self.axes[axcnt]=self.figure.add_subplot(nrows, nData, axcnt+1)
            # levels shared by the colormap contours and the sim overlay (overlay in sim units, i.e. /ma)
            levels = np.linspace(mi, ma, self.contourLevels)
            if self.mapMode == 'contour':
                if self.mapContours[axcnt] is not None:
                    try:
                        self.mapContours[axcnt].remove()
                    except Exception:
                        pass
                # contour only the visible window (+1 cell): the full 2048^2 array takes seconds
                ix = np.flatnonzero((x_axis >= xmi) & (x_axis <= xma))
                iy = np.flatnonzero((y_axis >= ymi) & (y_axis <= yma))
                ix = slice(max(ix[0]-1, 0), ix[-1]+2) if ix.size else slice(None)
                iy = slice(max(iy[0]-1, 0), iy[-1]+2) if iy.size else slice(None)
                self.mapContours[axcnt] = self.axes[axcnt].contour(
                    x_axis[ix], y_axis[iy], data[iy, ix], levels, cmap=self.parent.current_cmap,
                    vmin=mi, vmax=ma, linewidths=self.contourWidth)
                self.axes[axcnt].set_aspect('equal')   # imshow does this itself; contour does not
            elif makenew:
                im = self.axes[axcnt].imshow(
                    data, extent=axEx, origin="lower",
                    cmap=self.parent.current_cmap, interpolation="nearest")
                im.set_clim(vmin=mi, vmax=ma)
            else:
                for ch in self.axes[axcnt].get_children():
                    if isinstance(ch, matplotlib.image.AxesImage):
                        ch.set_data(data)
                        ch.set_extent(axEx)
                        ch.set_clim(vmin=mi, vmax=ma)
                        ch.set_cmap(self.parent.current_cmap)
                        break
            self.axes[axcnt].grid(self.showGrid)
            self.axes[axcnt].set_xlim(xmi, xma)
            self.axes[axcnt].set_ylim(ymi, yma)
            self.draw_fitbox(self.axes[axcnt], fitpoly, (xmi, xma, ymi, yma))
            tf = dd['field']
            self.axes[axcnt].set_title(f'B$_0$={tf} mT', fontsize = 'small')
            if showLabels:
                # self.axes[axcnt].set_xlabel(axlabel, fontsize = 'small')
                self.axes[axcnt].set_ylabel(axlabel, fontsize = 'small')
            

            rowcnt+=1
            
            if (dd.get('noise', 0)>0): 
                ratio = (dd['noisesky']/dd['noise']) 
            else:
                ratio = 1.0
            if dd.get('noiselevel', -1.0) <0:
                noiselevel = dd.get('noise') or 0.0   # older sessions have no noise entry
            else:
                noiselevel = np.max(data)*dd['noiselevel']


            ########## ----- simulation ---------------------------------
            upshift = None
            if self.showSim or self.showSkyline or self.showDiagonalProj:
                if type(dd['simdata'])!=type(None):
                    simdata = np.array(dd['simdata'])
                    simaxEx = [np.min(dd['simax']['x']), np.max(dd['simax']['x']),
                            np.min(dd['simax']['y']), np.max(dd['simax']['y'])]

                    # "Lift simulation to noise floor" 
                    upshift = noiselevel/ma
                    if self.noiseCompare == 'lift' and freqDomain is not None:
                        simdata *= (1-upshift)/np.max(np.max(simdata))
                        simdata += upshift
                    else:
                        simdata /=np.max(np.max(simdata))
                    sima = 1.0
                    #if dd['zmax'] < 0 and simbox is not None:
                    #    sima = np.max(simdata[simbox['idx0'], simbox['idx1']])
                    simi = sima*dd['zmin']/abs(dd['zmax']) if dd['zmax'] else 0.0
                else:
                    simdata = np.zeros_like(data)
                    simi,sima = (0.0, 1.0)
                    simaxEx= axEx
                    if type(dd['simax'])!=type(None):
                        realsim='simulation failed'
                    else:
                        realsim='no simulation'
            if self.showSim:
                simax = axcnt+rowcnt*nData
                realsim = ''
                if makenew:
                    self.axes[simax]=self.figure.add_subplot(nrows, nData, simax+1)
                else:
                    for ch in self.axes[simax].get_children():
                        if type(ch)==matplotlib.text.Annotation:
                            ch.remove()
                if self.mapMode == 'contour':
                    if self.mapContours[simax] is not None:
                        try:
                            self.mapContours[simax].remove()
                        except Exception:
                            pass
                    simx = np.linspace(simaxEx[0], simaxEx[1], simdata.shape[1])
                    simy = np.linspace(simaxEx[2], simaxEx[3], simdata.shape[0])
                    self.mapContours[simax] = self.axes[simax].contour(
                        simx, simy, simdata, np.linspace(simi, sima, self.contourLevels),
                        cmap=self.parent.current_cmap, vmin=simi, vmax=sima,
                        linewidths=self.contourWidth)
                    self.axes[simax].set_aspect('equal')
                elif makenew:
                    im = self.axes[simax].imshow(simdata, extent=simaxEx,
                        origin="lower", cmap=self.parent.current_cmap,
                        interpolation="nearest")
                    im.set_clim(vmin=simi, vmax=sima)
                else:
                    for ch in self.axes[simax].get_children():
                        if type(ch)== matplotlib.image.AxesImage:
                            ch.set_data(simdata)
                            ch.set_extent(simaxEx)
                            ch.set_clim(vmin=simi, vmax=sima)
                            ch.set_cmap(self.parent.current_cmap)
                self.axes[simax].grid(self.showGrid)

                ### Contours cannot have their data replaced in-place;
                ### we must remove the old QuadContourSet and draw a new one.
                if len(self.contours) > axcnt and self.contours[axcnt] is not None:
                    try:
                        self.contours[axcnt].remove()
                    except Exception:
                        pass
                    self.contours[axcnt] = None
                if len(self.contours) > axcnt+rowcnt*nData and self.contours[axcnt+rowcnt*nData] is not None:
                    try:
                        self.contours[axcnt+rowcnt*nData].remove()
                    except Exception:
                        pass
                    self.contours[axcnt+rowcnt*nData] = None 
                    
                if self.OverlaySim:
                    self.contours[axcnt] = self.axes[axcnt].contour(
                        dd['simax']['x'], dd['simax']['y'],
                        simdata, levels/ma,
                        colors=[self.overlayColor], linewidths=self.contourWidth)
                    self.contours[axcnt+rowcnt*nData] = self.axes[axcnt+rowcnt*nData].contour(
                        dd['simax']['x'], dd['simax']['y'],
                        simdata, levels/ma,
                        colors=[self.overlayColor], linewidths=self.contourWidth)

                # same visible window as the data panel above
                self.axes[simax].set_xlim(xmi, xma)
                self.axes[simax].set_ylim(ymi, yma)
                self.draw_fitbox(self.axes[simax], fitpoly, (xmi, xma, ymi, yma))
                if len(realsim)>0:
                    self.axes[simax].annotate("no simulation", (.0, .0), xycoords='axes points', color='w')
                rowcnt+=1
            ########## ----- skyline (max projection along each axis) ----------
            if self.showSkyline:
                skyax = axcnt+rowcnt*nData
                show_x = self.quadrant in ('all', 'horizontal')
                show_y = self.quadrant == 'vertical'
                if makenew:
                    self.axes[skyax] = self.figure.add_subplot(
                        nrows, nData, skyax+1,
                        sharex=self.axes[axcnt] if show_x else None)
                else:
                    self.axes[skyax].cla()
                self.projection_axes.append(self.axes[skyax])

                if show_x:
                    plo, phi, flo, fhi = xmi, xma, xfi, xfa
                else:
                    plo, phi, flo, fhi = ymi, yma, yfi, yfa

                sky_data = data[np.ix_(y_fit, x_fit)]

                if show_x:
                    proj_axis = x_axis[x_fit]
                    sky = np.max(sky_data, axis=0)
                    npoints = sky_data.shape[1]
                else:
                    proj_axis = y_axis[y_fit] 
                    sky = np.max(sky_data, axis=1)
                    npoints = sky_data.shape[0]
                    
                if self.noiseCompare == 'subtract':
                    sky = np.maximum(sky - noiselevel, 0.0)
                    
                if np.max(sky) > 0:
                    norm_val = np.max(sky)
                else:
                    norm_val = 1.0

                self.axes[skyax].plot(proj_axis, sky/norm_val, color='b', label='exp')
                
                if self.showNoiseFloor and self.noiseCompare != 'subtract':
                    nslvl = np.sqrt(np.log(npoints))*upshift
                    self.axes[skyax].axhline(nslvl, color='0.45', linestyle=':',
                                             linewidth=1.0, label='noise floor')
                                             
                ssx_fit = (np.max(dd['simax']['x']) >= xfi) & (np.min(dd['simax']['x']) <= xfa)
                ssy_fit = (np.max(dd['simax']['y']) >= yfi) & (np.min(dd['simax']['y']) <= yfa)
                if show_x and (np.max(dd['simax']['x'])>= xfi or np.min(dd['simax']['x'])<= xfa):
                    sproj_axis = dd['simax']['x'][ssx_fit], 
                    ssky_data = np.max(simdata[np.ix_(ssy_fit, ssx_fit)], axis=0)
                elif show_y and (np.max(dd['simax']['y'])>= yfi or np.min(dd['simax']['y'])<= yfa):
                    sproj_axis = dd['simax']['y'][ssy_fit]
                    ssky_data = np.max(simdata[np.ix_(ssy_fit, ssx_fit)], axis=1)
                else:
                    ssky_data = 0.0
                if np.max(ssky_data)>0:
                    self.axes[skyax].plot(sproj_axis, ssky_data/np.max(ssky_data) , color='r',
                                              linestyle='--', label='sim')
                    
                    self.axes[skyax].set_xlim(plo, phi)
                # Grey the padding, the same wash the 2D maps get.
                if self.showFitBox:
                    for a, b in ((plo, flo), (fhi, phi)):
                        if b > a:
                            self.axes[skyax].axvspan(a, b, facecolor='0.5',
                                                     edgecolor='none', alpha=0.18)
                self.axes[skyax].set_ylabel('Amplitude (norm.)', fontsize='small')
                self.axes[skyax].legend(fontsize='x-small', loc='upper right')
                self.axes[skyax].grid(self.showGrid)
                rowcnt += 1

            ########## ----- show orientation selection -----------
          
            if self.showOriSel:
                oriax = axcnt+rowcnt*nData
                if makenew:  
                    self.axes[oriax]=self.figure.add_subplot(nrows, nData, oriax+1, projection='3d')
                else:
                    self.axes[oriax].cla()
                if type(dd['simax'])!=type(None):
                    phi =dd['simax']['orisel'][0, :]
                    theta =dd['simax']['orisel'][1, :]
                    ak =dd['simax']['orisel'][2, :]
                    # Convert spherical to Cartesian for plotting
                    x = np.sin(theta)*np.cos(phi)
                    y = np.sin(theta)*np.sin(phi)
                    z = np.cos(theta)
                    amax = np.max(ak)
                    norm = plt.Normalize(vmin=0, vmax=amax if amax > 0 else 1.0)

                    # Triangulation from the simulation grid (make_grid). It is
                    # built on the sphere, so it stays correct for the 'spiral'
                    # grid, which covers the full sphere -- there the old
                    # Triangulation(x, y) folded the two hemispheres together.
                    triangles = dd['simax'].get('oriseltri')
                    if triangles is None or len(triangles) == 0 or np.max(triangles) >= len(x):
                        # precalculated grids carry no triangulation of their own
                        triangles = mtri.Triangulation(x, y).triangles

                    triangle_ak = ak[triangles].max(axis=1)          # one value per triangle
                    facecolors = self.parent.current_cmap(norm(triangle_ak))

                    surf = self.axes[oriax].plot_trisurf(
                        x, y, z,
                        triangles=triangles,
                        linewidth=0, antialiased=True,
                        shade=False,  # we supply colour via facecolors
                        )
                        #
                    self.axes[oriax].set_aspect('equal')
                    self.axes[oriax].set_xlabel('X'); self.axes[oriax].set_ylabel('Y'); self.axes[oriax].set_zlabel('Z')
                    surf.set_facecolor(facecolors)
                    # per-triangle, like the faces: ak is per-knot and would be
                    # the wrong length here
                    surf.set_edgecolor(facecolors)
                rowcnt+=1
            if self.showDiagonalProj:
                diagprj = axcnt+rowcnt*nData
                if makenew:  
                    self.axes[diagprj]=self.figure.add_subplot(nrows, nData, diagprj+1)
                else:
                    self.axes[diagprj].cla()
                self.projection_axes.append(self.axes[diagprj])
                if box is not None:
                    ma = scale*abs(dd['zmax']) or 1.0
                    skyprj = np.max(activeBox, axis=1)/ma
                    nslvl = np.sqrt(np.log(box['nalong']))*upshift

                    xw = x_axis[box['sel1']]
                    xaX = np.linspace(xw[0], xw[-1], box['ndiag']) - (xw[0]+xw[-1])/2
                    if self.noiseCompare == 'subtract':
                        skyprj = np.maximum(skyprj - nslvl, 0.0)/(1-nslvl)
                    self.axes[diagprj].plot(xaX, skyprj, color='b')
                    if self.showNoiseFloor and self.noiseCompare != 'subtract':
                        self.axes[diagprj].plot(xaX, np.full(xaX.shape, nslvl), color='0.45',
                                                linestyle=':', linewidth=1.0, label='noise floor')

                    if simbox is not None:
                        simdata = dd['simdata']
                        sma = np.max(simdata) or 1.0
                        skyprj1 = np.max(simdata[simbox['idx0'], simbox['idx1']]/sma, axis=1)
                        sxw = simx_axis[simbox['sel1']]
                        sxaX = (np.linspace(sxw[0], sxw[-1], simbox['ndiag'])
                                - (sxw[0]+sxw[-1])/2)
                        self.axes[diagprj].plot(sxaX, skyprj1, color='r')
                self.axes[diagprj].grid(self.showGrid)
                rowcnt+=1

            axcnt+=1
        if makenew:
            self.figure.tight_layout(pad=1.2)
        self.canvas.draw()
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kwargs):
        super().__init__(parent,size=(200, 400), *args, **kwargs)
        self.parent = mainWindow
        
        self.autoFFT = False
        self.autoSim = False
        # [label, default, choices, type] -- floats carry a fifth entry, the
        # spin editor's initial step (right-clicking the spin box changes it).
        self.FFTparams = {
            "polynomial":["Polynom", 2, [], 'int'],
            "apodization":["Apodization", ['Hamming'], ['Hamming', 'Lor-Gau', 'Gaussian'], 'list'],
            "awidth":["A.Width",  1, [], 'float', 0.05],
            "aalpha":["Lor.Width",0.6, [], 'float', 0.05],
            "ashift":["A.Shift",  0.0, [], 'float', 0.05],
            "zfill": ["Zero fill" , 2, [], 'int'],
            "imag":  ["Use 'imaginary'",  False, [], 'bool'],
            }
        
        self.Color_BG_FILE_Title = wx.Colour(theme.theme_colors()["accent"])
        self.Color_BG_FILE_Inact = wx.Colour(theme.theme_colors()["tab_inactive"])
        self.Color_FG_FILE_Title = wx.Colour(theme.theme_colors()["propgrid_caption_fg"])
        self.Color_BG_FILE_Main = wx.Colour(theme.theme_colors()["ctrl_bg"])
        self.Color_FG_FILE_Main = wx.Colour(theme.theme_colors()["text"])

        self.Color_BG_FFT_Title = wx.Colour(theme.theme_colors()["pill_load"])
        self.Color_BG_FFT_Inact = wx.Colour(theme.theme_colors()["tab_inactive"])
        self.Color_FG_FFT_Title = wx.Colour(theme.theme_colors()["propgrid_caption_fg"])
        self.Color_BG_FFT_Main = wx.Colour(theme.theme_colors()["ctrl_bg"])
        self.Color_FG_FFT_Main = wx.Colour(theme.theme_colors()["text"])

        # PillTabBar (tab strip) driving a wx.Simplebook (page container) --
        # the book must exist before add_*panel() below, since those call
        # self.nb_book.AddPage(); the tab bar itself is built afterwards,
        # once the pages (and hence their labels/order) exist.
        self.nb_book = wx.Simplebook(self)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.add_datapanel()
        self.add_syspanel()
        self.add_exppanel()
        self.add_optpanel()

        self.nb = PillTabBar(self, self.nb_book, ["Data", "Sys", "Exp", "Opt"])
        sizer.Add(self.nb, 0, wx.EXPAND)
        sizer.Add(self.nb_book, 1, wx.EXPAND)

        # Do FFT / Do SIM / auto FFT / auto SIM now live as toolbar icons on
        # MainFrame (self.parent.tb_doFFT etc.) -- see MainFrame.__init__.
        # Bind them here now that self.on_doFFT/self.on_update_sim exist.
        self.parent.tb_doFFT.Bind(wx.EVT_BUTTON, self.on_doFFT)
        self.parent.tb_doSIM.Bind(wx.EVT_BUTTON, self.on_update_sim)
        self.parent.tb_gridsearch.Bind(wx.EVT_BUTTON, self.on_open_grid_search)

        self.SetSizer(sizer)

        # Store for later
        self.param_definitions = None
        self.param_controls   = {}

        self.ffttree.Bind(wxpg.EVT_PG_CHANGED, self.on_ffttree)
        #self.filetree.Bind(wxpg.EVT_PG_DOUBLE_CLICK, self.on_filetree) #EVT_G_SELECTED EVT_PG_RIGHT_CLICK
        #self.tree.Bind(CT.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit)
        #self.tree.Bind(CT.EVT_TREE_END_LABEL_EDIT, self.on_end_edit)        
        self.filetree.Bind(wxpg.EVT_PG_CHANGED, self.on_filetree)
        self.filetree.Bind(wxpg.EVT_PG_SELECTED, self.on_filetree_clicked)
        self.filetree.Bind(wxpg.EVT_PG_DOUBLE_CLICK, self.on_filetree_dclicked)
        #self.filetree.Bind(wxpg.EVT_PG_DOUBLE_CLICK, self.on_filetree) #EVT_G_SELECTED EVT_PG_RIGHT_CLICK
        #self.tree.Bind(CT.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit)
        #self.tree.Bind(CT.EVT_TREE_END_LABEL_EDIT, self.on_end_edit)
        
    def add_optpanel(self):
        # ------------- Opt
        
        self.Opt_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Opt_panel, "Opt") 
        Optsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Opt_param = MypgPanel.PropGridPanel(self.Opt_panel, Prop_Dict = {}, onChangeFunc=self.on_Opt)
        self.Opt_param.tooltips =self.parent.Opt.getToolTips()
        self.Opt_param.SetFromParClean(self.parent.Opt.getDefaultDict())

        Optsizer.Add(self.Opt_param, 1, wx.EXPAND)

        # ---- Target Grid Points ------------------------------------------
        # Off-grid helper for nKnots: nKnots means a different thing in every
        # Grid method, so this lets the resolution be set by the one figure
        # that is comparable. Kept out of the property grid because it is not
        # an Opt parameter -- nKnots remains the stored one, and the two stay
        # in sync in both directions via updateGridInfo().
        #tgtBox = wx.StaticBoxSizer(
        #    wx.StaticBox(self.Opt_panel, label="Target Grid Points"), wx.VERTICAL)
        
        Optsizer.Add(wx.StaticText(self.Opt_panel, label="\u25a0 Target Grid Points"), 0, wx.LEFT, 0)
        Optsizer.Add(wx.StaticLine(self.Opt_panel, wx.LI_HORIZONTAL), 0, wx.EXPAND)
        tgtBox = wx.BoxSizer(wx.VERTICAL)
        self.spin_gridTarget = wx.SpinCtrl(
            self.Opt_panel, min=1, max=10000000,
            initial=int(self.parent.Opt.nGridPoints() or 1),
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_gridTarget.SetToolTip(
            "Requested number of orientations. nKnots is integer, so the grid "
            "snaps to the nearest achievable size -- the result is shown below.")
        tgtBox.Add(self.spin_gridTarget, 0, wx.EXPAND | wx.ALL, 2)

        self.txt_gridActual = wx.TextCtrl(
            self.Opt_panel, style=wx.TE_READONLY | wx.TE_CENTRE|wx.BORDER_NONE )
        self.txt_gridActual.SetToolTip(
            "Grid actually produced for the current Grid method and nKnots.")
        tgtBox.Add(self.txt_gridActual, 0, wx.EXPAND | wx.ALL, 2)
        Optsizer.Add(tgtBox, 0, wx.EXPAND | wx.TOP, 4)

        self.spin_gridTarget.Bind(wx.EVT_SPINCTRL, self.on_gridTarget)
        self.spin_gridTarget.Bind(wx.EVT_TEXT_ENTER, self.on_gridTarget)
        self.updateGridInfo()

        self.btn_doOriSel   = PillButton(self.Opt_panel, "Update Ori.Sel. Grid", color_key="pill_load",
                                         draw_badge=theme.draw_play_badge)
        Optsizer.Add(self.btn_doOriSel, 0, wx.EXPAND)
        self.btn_doOriSel.Bind(wx.EVT_BUTTON, self.parent.on_update_orisel)

        # ---- Contours -----------------------------------------------------
        # mirrors MatplotlibPanel.contourLevels/contourWidth/overlayColor; toolbar toggle picks the mode
        mp = self.parent.matplotlib_panel
        #cntBox = wx.StaticBoxSizer(
        #    wx.StaticBox(self.Opt_panel, label="Contours"), wx.VERTICAL)
        Optsizer.Add(wx.StaticText(self.Opt_panel, label="\u25a0 Contours"), 0, wx.LEFT, 0)
        Optsizer.Add(wx.StaticLine(self.Opt_panel, wx.LI_HORIZONTAL), 0, wx.EXPAND)
        cntBox = wx.BoxSizer(wx.VERTICAL)
        cntGrid = wx.FlexGridSizer(3, 2, 2, 4)
        cntGrid.AddGrowableCol(1, 1)
        self.spin_contourLevels = wx.SpinCtrl(
            self.Opt_panel, min=2, max=100, initial=mp.contourLevels,
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_contourLevels.SetToolTip(
            "Number of linear levels between z.min and z.max (also used by the sim overlay).")
        self.spin_contourWidth = wx.SpinCtrlDouble(
            self.Opt_panel, min=0.1, max=5.0, inc=0.1, initial=mp.contourWidth,
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_contourWidth.SetDigits(1)
        self.spin_contourWidth.SetToolTip("Contour line width, points.")
        self.pick_overlayColor = wx.ColourPickerCtrl(
            self.Opt_panel, colour=wx.Colour(mp.overlayColor))
        self.pick_overlayColor.SetToolTip("Uniform colour of the 'Overlay Sim' contour lines.")
        theme.theme_control(self.spin_contourLevels)
        theme.theme_control(self.spin_contourWidth)
        for label, ctrl in (("Levels", self.spin_contourLevels),
                            ("Line width", self.spin_contourWidth),
                            ("Overlay colour", self.pick_overlayColor)):
            cntGrid.Add(wx.StaticText(self.Opt_panel, label=label), 0,
                        wx.ALIGN_CENTER_VERTICAL | wx.LEFT, 2)
            cntGrid.Add(ctrl, 1, wx.EXPAND)
        cntBox.Add(cntGrid, 0, wx.EXPAND | wx.ALL, 2)
        Optsizer.Add(cntBox, 0, wx.EXPAND | wx.TOP, 4)

        self.spin_contourLevels.Bind(wx.EVT_SPINCTRL, self.on_contour_style)
        self.spin_contourLevels.Bind(wx.EVT_TEXT_ENTER, self.on_contour_style)
        self.spin_contourWidth.Bind(wx.EVT_SPINCTRLDOUBLE, self.on_contour_style)
        theme.bind_spin_enter(self.spin_contourWidth, lambda: self.on_contour_style(None))
        self.pick_overlayColor.Bind(wx.EVT_COLOURPICKER_CHANGED, self.on_contour_style)
        self.sync_contour_controls()

        # ---- View pad / Fit box -------------------------------------------
        # display only, state on MatplotlibPanel (viewPad / showFitBox); handlers live there too
        
        #viewBox = wx.StaticBoxSizer(
        #    wx.StaticBox(self.Opt_panel, label="View pad, MHz"), wx.HORIZONTAL)
        Optsizer.Add(wx.StaticText(self.Opt_panel, label="\u25a0 View pad, MHz"), 0, wx.LEFT, 0)
        Optsizer.Add(wx.StaticLine(self.Opt_panel, wx.LI_HORIZONTAL), 0, wx.EXPAND)
        viewBox = wx.BoxSizer(wx.HORIZONTAL)
        self.spin_viewPad = wx.SpinCtrlDouble(
            self.Opt_panel, min=0.0, max=10000.0, inc=1.0, initial=mp.viewPad,
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER)
        self.spin_viewPad.SetDigits(2)
        self.spin_viewPad.SetToolTip(
            "Extra MHz drawn around the f.min..f.max analysis window.\n"
            "Display only: the fit region, the colour scale, the projection\n"
            "normalization and the noise floor all stay on f.min..f.max.\n"
            "Clamped to the measured spectrum, and never crosses a quadrant cut.")
        theme.theme_control(self.spin_viewPad)
        self.chk_FitBox = PillCheckBox(self.Opt_panel, "Fit box")
        self.chk_FitBox.SetValue(mp.showFitBox)
        viewBox.Add(self.spin_viewPad, 1, wx.ALL, 3)
        viewBox.Add(self.chk_FitBox, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 3)
        Optsizer.Add(viewBox, 0, wx.EXPAND | wx.TOP, 4)
        self.spin_viewPad.Bind(wx.EVT_SPINCTRLDOUBLE, mp.on_viewpad)
        self.spin_viewPad.Bind(wx.EVT_KILL_FOCUS, mp.on_viewpad)
        theme.bind_spin_enter(self.spin_viewPad, lambda: mp.on_viewpad(None))   # Enter is unreliable on SpinCtrlDouble
        self.chk_FitBox.Bind(wx.EVT_CHECKBOX, mp.on_check)

        # ---- Quadrants ----------------------------------------------------
        #quadBox = wx.StaticBoxSizer(
        #    wx.StaticBox(self.Opt_panel, label="Quadrants"), wx.VERTICAL)   # stacked: three across clip at the default pane width
        quadBox = wx.BoxSizer(wx.VERTICAL)
        Optsizer.Add(wx.StaticLine(self.Opt_panel, wx.LI_HORIZONTAL), 0, wx.EXPAND)
        Optsizer.Add(wx.StaticText(self.Opt_panel, label="\u25a0 Quadrants"), 0, wx.LEFT, 0)
        
        self.quadrant_group = []
        self.button_all4 = ThemedRadioButton(self.Opt_panel, "All 4", self.quadrant_group, True)
        self.button_horiz = ThemedRadioButton(self.Opt_panel, "±ν₁  (horiz)", self.quadrant_group)
        self.button_vert = ThemedRadioButton(self.Opt_panel, "±ν₂  (vert)", self.quadrant_group)
        quadBox.Add(self.button_all4, 0, wx.ALL, 3)
        quadBox.Add(self.button_horiz, 0, wx.ALL, 3)
        quadBox.Add(self.button_vert, 0, wx.ALL, 3)
        self.button_all4.BindRadio(mp.on_quadrant)
        self.button_horiz.BindRadio(mp.on_quadrant)
        self.button_vert.BindRadio(mp.on_quadrant)
        self.button_all4.SetValue(True)
        Optsizer.Add(quadBox, 0, wx.EXPAND | wx.TOP, 4)

        self.Opt_panel.SetSizer(Optsizer)

    def on_contour_style(self, event):
        if event is not None:
            event.Skip()
        self.parent.matplotlib_panel.set_contour_style(
            levels=self.spin_contourLevels.GetValue(),
            width=self.spin_contourWidth.GetValue(),
            overlayColor=self.pick_overlayColor.GetColour().GetAsString(wx.C2S_HTML_SYNTAX))

    def sync_contour_controls(self):
        """Grey out what has no effect: levels/width need contours or the overlay, colour needs the overlay."""
        mp = self.parent.matplotlib_panel
        self.spin_contourLevels.Enable(mp.mapMode == 'contour' or mp.OverlaySim)
        self.spin_contourWidth.Enable(mp.mapMode == 'contour' or mp.OverlaySim)
        self.pick_overlayColor.Enable(mp.OverlaySim)

    def add_setpanel(self):
        # ------------- Opt
        self.Set_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Set_panel, "Opt")
        Optsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Set_param = MypgPanel.PropGridPanel(self.Set_panel, Prop_Dict = {}, onChangeFunc=self.on_Set)
        #self.Set_param.tooltips =self.parent.Opt.getToolTips()
        #self.Set_param.SetFromParClean(self.parent.Opt.getDefaultDict())

        Optsizer.Add(self.Set_param, 1, wx.EXPAND)

        self.Set_param.SetSizer(Optsizer)
        
    def add_exppanel(self):
        # ------------- Exp
        self.Exp_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Exp_panel, "Exp") 
        
        Expsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Exp_param = MypgPanel.PropGridPanel(self.Exp_panel, Prop_Dict = {}, onChangeFunc=self.on_Exp)
        self.Exp_param.tooltips =self.parent.Exp.getToolTips()
        self.Exp_param.SetFromParClean(self.parent.Exp.getDefaultDict())

        Expsizer.Add(self.Exp_param, 1, wx.EXPAND)
        self.Exp_panel.SetSizer(Expsizer)
         
    def add_syspanel(self):
        # ------------- Sys 
        self.Sys_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.Sys_panel, "Sys") 
        
        Syssizer = wx.BoxSizer(wx.VERTICAL)
        
        btnszr = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_add = PillButton(self.Sys_panel, "Nuc", color_key="accent",
                                  draw_badge=theme.draw_plus_badge)
        self.btn_add.SetToolTip("Add a nucleus to the spin system.")
        self.btn_add.Bind(wx.EVT_BUTTON, self.on_add_Nuc)
        btnszr.Add(self.btn_add, 0, wx.EXPAND)
        self.btn_delete = PillButton(self.Sys_panel, "Nuc", color_key="pill_warn",
                                     draw_badge=theme.draw_minus_badge)
        self.btn_delete.SetToolTip("Delete the selected nucleus from the spin system.")
        self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete_nuc)
        btnszr.Add(self.btn_delete, 0, wx.EXPAND)

        Syssizer.Add(btnszr, 0, wx.EXPAND)
        
        self.Sys_param = MypgPanel.PropGridPanel(self.Sys_panel, Prop_Dict = {},
                                                  onChangeFunc=self.on_Sys, showModFunc=True)
        self.Sys_param.tooltips =self.parent.Sys.getToolTips()
        self.Sys_param.SetFromParClean(self.parent.Sys.getDefaultDict())

        Syssizer.Add(self.Sys_param, 1, wx.EXPAND)
      
        self.Sys_panel.SetSizer(Syssizer)
               
    def add_datapanel(self):
        # ----------Data panel
        self.data_panel = wx.Panel(self.nb_book)
        self.nb_book.AddPage(self.data_panel, "Data")
        
        sizer = wx.BoxSizer(wx.VERTICAL)

        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_delete = PillButton(self.data_panel, "Data", color_key="pill_warn",
                                     draw_badge=theme.draw_cross_badge)
        self.btn_delete.SetToolTip("Delete the selected data set.")
        self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete)
        sizerBtns.Add(self.btn_delete, 1, wx.ALIGN_CENTER | wx.ALL, 5)

        self.btn_dsc = PillButton(self.data_panel, "DSC", color_key="pill_load",
                                  draw_badge=theme.draw_textlines_badge)
        self.btn_dsc.SetToolTip("Show the DSC parameter file of the selected data set.")
        self.btn_dsc.Bind(wx.EVT_BUTTON, self.on_dsc)
        sizerBtns.Add(self.btn_dsc, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        sizer.Add(sizerBtns, 0, wx.ALIGN_LEFT | wx.ALL, 5)

        row_expand = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_expanded = PillCheckBox(self.data_panel, "Expand all")
        self.chk_expanded.SetValue(False)
        self.chk_expanded.Bind(wx.EVT_CHECKBOX, self.on_expandchk)
        row_expand.Add(self.chk_expanded, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 12)
        row_expand.AddStretchSpacer(1)

        # bmp_up   = wx.ArtProvider.GetBitmap(wx.ART_GO_UP,   wx.ART_BUTTON, (16, 16))
        # bmp_down = wx.ArtProvider.GetBitmap(wx.ART_GO_DOWN, wx.ART_BUTTON, (16, 16))
        # Tiny icon-only utility buttons -- kept as plain native wx.Button
        # (theme-retinted by theme.apply_theme_to_window) rather than pill
        # widgets, since the label-pill shape doesn't suit a 24px glyph button.
        self.btn_moveUp   = wx.Button(self.data_panel, label="▲", size=(24, 22), style=wx.BORDER_NONE)
        self.btn_moveDown = wx.Button(self.data_panel, label="▼", size=(24, 22), style=wx.BORDER_NONE)
        self.btn_moveUp.SetToolTip("Move selected dataset up (earlier column)")
        self.btn_moveDown.SetToolTip("Move selected dataset down (later column)")
        self.btn_moveUp.Bind(  wx.EVT_BUTTON, lambda e: self.on_move_data(-1))
        self.btn_moveDown.Bind(wx.EVT_BUTTON, lambda e: self.on_move_data(+1))
        row_expand.Add(self.btn_moveUp,   0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 2)
        row_expand.Add(self.btn_moveDown, 0, wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(row_expand, 0, wx.EXPAND | wx.ALL, 4)
        # file tree
        self.filetree = wxpg.PropertyGrid(
            self.data_panel,
            style=wxpg.PG_DEFAULT_STYLE | 
                  wxpg.PG_SPLITTER_AUTO_CENTER |
                  wxpg.PG_HIDE_MARGIN |
                  wxpg.PG_TOOLTIPS |
                  wxpg.PG_NO_INTERNAL_BORDER |wx.BORDER_NONE
        )
        
        self.filetree.SetCellBackgroundColour(self.Color_BG_FILE_Main)
        self.filetree.SetCellTextColour(self.Color_FG_FILE_Main)
        self.filetree.SetMarginColour(self.Color_BG_FILE_Main)
        #self.filetree.SetSelectionBackgroundColour(wx.Colour(200, 255, 255))
        
        self.filetree.SetCaptionBackgroundColour(self.Color_BG_FILE_Inact)
        self.filetree.SetCaptionTextColour(self.Color_FG_FILE_Title)
        self.filetree.SetEmptySpaceColour(self.Color_BG_FILE_Main)

        # "BG Scale" uses the app's spin-with-step-menu float editor. Editors
        # live in one registry shared by every grid, but this panel is built
        # before any PropGridPanel is, so it can't rely on one of those
        # having registered it already.
        if not self.filetree.GetEditorByName("SpinFloat"):
            self.filetree.RegisterEditor(MypgPanel.SpinCtrlDoubleEditor, "SpinFloat")
        # cat1 = self.filetree.Append(wxpg.PropertyCategory("Nothing loaded yet"))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("MW Freq", value=9.43))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("B0", value=350.0))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("tau", value=120.0))
        
        sizer.Add(self.filetree, 1, wx.EXPAND | wx.ALL, 0)
        
        lbl = wx.StaticText(self.data_panel, label=" FFT Parameters:")
        lbl.SetBackgroundColour(self.Color_BG_FFT_Title)
        lbl.SetForegroundColour(self.Color_FG_FFT_Title)
        lfont=lbl.GetFont(); lfont.SetWeight(wx.FONTWEIGHT_BOLD); lbl.SetFont(lfont)              # start with the current font

        sizer.Add(lbl, 0, wx.EXPAND, 0)
        # FFT tree
        self.ffttree = wxpg.PropertyGrid(
            self.data_panel,
            style=wxpg.PG_DEFAULT_STYLE | 
                  wxpg.PG_SPLITTER_AUTO_CENTER |wxpg.PG_HIDE_MARGIN |
                  wxpg.PG_NO_INTERNAL_BORDER|wx.BORDER_NONE
        )
        self.ffttree.SetCellBackgroundColour(self.Color_BG_FFT_Main)
        self.ffttree.SetCellTextColour(self.Color_FG_FFT_Main)
        self.ffttree.SetMarginColour(self.Color_BG_FFT_Title)
        #self.filetree.SetSelectionBackgroundColour(wx.Colour(200, 255, 255))
        
        self.ffttree.SetCaptionBackgroundColour(self.Color_BG_FFT_Title)
        self.ffttree.SetCaptionTextColour(self.Color_FG_FFT_Title)
        #cat2 = self.ffttree.Append(wxpg.PropertyCategory("FFT"))
        

        for kk in self.FFTparams.keys():
            if self.FFTparams[kk][3]=='int':
                self.ffttree.Append(wxpg.IntProperty(self.FFTparams[kk][0], kk, value=self.FFTparams[kk][1]))
            elif self.FFTparams[kk][3]=='float':
                # Same spin-with-step-menu editor as "BG Scale" above -- these
                # are all knobs one nudges and re-transforms, not values typed
                # once, and the editor commits on Enter as well as on the
                # arrows (plain FloatProperty needed Enter to be caught by the
                # grid, which is exactly what it kept swallowing).
                pro = self.ffttree.Append(
                    wxpg.FloatProperty(self.FFTparams[kk][0], kk, value=self.FFTparams[kk][1]))
                self.ffttree.SetPropertyEditor(pro, "SpinFloat")
                pro.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, self.FFTparams[kk][4])
            elif self.FFTparams[kk][3]=='bool':
                self.ffttree.Append(wxpg.BoolProperty(self.FFTparams[kk][0], kk, value=self.FFTparams[kk][1]))
            elif self.FFTparams[kk][3]=='list':
                self.ffttree.Append(wxpg.MultiChoiceProperty(self.FFTparams[kk][0], kk, 
                                                                     choices=self.FFTparams[kk][2],
                                                                     value = self.FFTparams[kk][1]))
        nprops = len(self.FFTparams.keys())
        row_h = self.ffttree.GetRowHeight()
        self.ffttree.SetMinSize(( -1, nprops*row_h+4 ))
        #self.ffttree.Layout()
        #height = self.ffttree.GetBestSize().GetHeight()
        #self.ffttree.SetSize((-1, 2*height))

        # self.ffttree.AppendIn(cat2, wxpg.IntProperty("Polynom", "polynomial", value=2))
        # self.ffttree.AppendIn(cat2, wxpg.MultiChoiceProperty("Apodization", "apodization", 
        #                                                      choices=['Hamming', 'Lor-Gau', 'Gaussian'], 
        #                                                      value=['Hamming']))
        # self.ffttree.AppendIn(cat2, wxpg.FloatProperty("A.Width", "awidth", value=0.6))
        # self.ffttree.AppendIn(cat2, wxpg.FloatProperty("Lor.Width", "aalpha", value=0.6))
        # self.ffttree.AppendIn(cat2, wxpg.FloatProperty("A.Shift", "ashift", value=0))
        # self.ffttree.AppendIn(cat2, wxpg.IntProperty("Zero fill", "zfill", value=1))
        # self.ffttree.AppendIn(cat2, wxpg.BoolProperty("Use 'imaginary'", "imag", value=False))
        

        
        sizer.Add(self.ffttree, 0, wx.EXPAND | wx.ALL, 5)
        
        self.data_panel.SetSizer(sizer)
    def on_add_Nuc(self, event):
        items_to_show = []
        isotopes = self.parent.Sys.isotopes()
        isotopes_list = list(self.parent.Sys.isotopes().keys())
        for dd in isotopes_list:
            items_to_show.append(f'{dd:>8}\tI={isotopes[dd]["I"]:>4}')
        
        dlg = wx.SingleChoiceDialog(
            self,
            "Select isotope to include:",
            "Choose Nucleus",
            items_to_show,
            style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
        )
    
        # Show the dialog modally â it will block until the user presses OK/Cancel
        if dlg.ShowModal() == wx.ID_OK:
            # dlg.GetSelections() returns indices of the checked items
            selected_idx = dlg.GetSelection()
        else:
            selected_idx = -1

        dlg.Destroy()
        if selected_idx>=0:
            
            Params = self.Sys_param.parameters
            # Count first, then append exactly one entry. This assignment used
            # to sit inside the loop, so with N nuclei it wrote nuc(2)..nuc(N+1)
            # and overwrote nuclei 2..N with the newly picked isotope.
            nNucs = sum(1 for key in Params if 'nuc' in key)
            Params[f'nuc({nNucs+1})']=self.getNuc(isotopes_list[selected_idx])

            self.parent.Sys.setFromCtrl(Params)
            self.Sys_param.SetFromParClean(Params)

            # for de in selected_indices:
            #     self.parent.Data.pop(de)
            # self.parent.update_everything()
    def getNuc(self, nucname):

        isotopes = list(self.parent.Sys.isotopes().keys())
        selected_idx = isotopes.index(nucname)
        nucdict = {'Nucs':[isotopes[selected_idx], isotopes],
                                 'nNucs':1,
                                 'A':np.array([0, 0, 0]),
                                 'Apa':np.array([0, 0, 0])}
        _,I,_=self.parent.Sys.isotopes(isotopes[selected_idx])
        if I[0]>0.5:
           nucdict['Q']=np.array([0, 0, 0])
           nucdict['Qpa']=np.array([0, 0, 0])
        nucdict['useFor'] = ['sim.only', ['sim.only', 'ori.sel.only', 'all']]
        return nucdict

    def on_expandchk(self, event):
        if self.chk_expanded.GetValue():
            self.filetree.ExpandAll()
        else:
            self.filetree.CollapseAll()
        
        #self.update_filelist()
    def on_delete_nuc(self, event):
        pass
    def on_open_grid_search(self, event):
        frame = grid_search.GridSearchFrame(self.parent)
        frame.Show()
    def on_update_sim(self, event):
        self.parent.runSim()
    def on_Sys(self, mainname, name, val):
        self.parent.Sys.setFromCtrl(self.Sys_param.parameters) # This is more to double check that input works
        if self.parent.tb_autoSIM.GetValue():
            self.parent.runSim()
        
    def on_Exp(self, mainname, name, val):
        self.parent.Exp.setFromCtrl(self.Exp_param.parameters)                            
    def on_Opt(self, mainname, name, val):
        self.parent.Opt.setFromCtrl(self.Opt_param.parameters)
        if name=='OriSelType':
            self.parent.on_update_orisel(None)
        if name in ('Grid', 'nKnots'):
            self.updateGridInfo()

    def on_gridTarget(self, event):
        """Convert a requested grid size into nKnots and push it to the grid.

        Only a ladder of totals is reachable, so the value asked for and the
        value obtained generally differ; the read-only box below reports what
        was actually built. Written with trigger_run=False so syncing nKnots
        does not re-enter on_Opt while this handler is still running.
        """
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
        """Keep nKnots, Target Grid Points and the actual count consistent.

        nKnots means a different thing in every grid method (89 / 761 / 968 at
        nKnots=20), so the knot count is the only comparable figure -- and the
        only one that predicts run time, which is one diagonalisation per knot.
        Called whenever Grid or nKnots changes and after a session is loaded.

        keepTarget leaves the spin control alone, so a request of 800 that
        snapped to 761 still reads 800 rather than silently rewriting itself.
        """
        Opt = self.parent.Opt
        n = Opt.nGridPoints()
        gname = Opt.Grid[0] if isinstance(Opt.Grid, (list, tuple)) else Opt.Grid

        txt = Opt.getToolTips('nKnots')
        if n is not None:
            txt = (f"{txt}\n\n"
                   f"Currently: '{gname}' with nKnots={Opt.nKnots}"
                   f"  ->  {n} grid points")
        self.Opt_param.tooltips['nKnots'] = txt
        prop = self.Opt_param.pg.GetPropertyByName('nKnots')
        if prop is not None:
            prop.SetHelpString(txt)

        if getattr(self, 'txt_gridActual', None) is not None:
            self.txt_gridActual.SetValue(
                'grid: n/a' if n is None
                else f'actual: {n} points  (nKnots = {Opt.nKnots})')
        # nKnots edited directly, or a session loaded -> follow it
        if not keepTarget and getattr(self, 'spin_gridTarget', None) is not None:
            if n is not None and self.spin_gridTarget.GetValue() != n:
                self.spin_gridTarget.SetValue(int(n))
    def on_Set(self, mainname, name, val):
        pass
    def on_dsc(self, event):
        items_to_show = []
        # selected_items = []
        for dd in self.parent.Data:
            items_to_show.append(dd['title'])
        
        dlg = wx.SingleChoiceDialog(
            self,
            "Select the data you want to display:",
            "Choose data",
            items_to_show,
            style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
        )
    
        # Show the dialog modally â it will block until the user presses OK/Cancel
        if dlg.ShowModal() == wx.ID_OK:
            # dlg.GetSelections() returns indices of the checked items
            selected_indices = dlg.GetSelection()
            #print(self.parent.Data[selected_indices]['dsc'])
            dlg.Destroy()
            
            app = wx.App(False)
            frame = wx.Frame(None)  # parent frame (not shown)
            dlg = DictPopup(frame, self.parent.Data[selected_indices]['dsc'])
            dlg.ShowModal()
            dlg.Destroy()
            frame.Destroy()
            app.MainLoop()
            
    def on_delete(self, event):
        items_to_show = []
        # selected_items = []
        for dd in self.parent.Data:
            items_to_show.append(dd['title'])
        
        dlg = wx.MultiChoiceDialog(
            self,
            "Select the data you want to delete:",
            "Choose data to delete",
            items_to_show,
            style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
        )
    
        # Show the dialog modally â it will block until the user presses OK/Cancel
        if dlg.ShowModal() == wx.ID_OK:
            # dlg.GetSelections() returns indices of the checked items
            selected_indices = dlg.GetSelections()
        else:
            selected_indices = []

        dlg.Destroy()
        
        if len(selected_indices)>0:
            for de in selected_indices:
                self.parent.Data.pop(de)
            self.parent.update_everything()

    def _get_selected_data_index(self):
        """
        Return the index in self.parent.Data of the file that is currently
        selected in filetree, or None if nothing file-related is selected.

        Works whether the user has clicked the file's category header or any
        of its child properties, because children are named "<key> <index>"
        (e.g. "show 0", "fname 2") by update_filelist().
        """
        prop = self.filetree.GetSelection()
        if prop is None:
            return None

        # Case 1: a child property is selected -> name ends with " <index>"
        try:
            name = prop.GetName()
            parts = name.split()
            if len(parts) >= 2 and parts[-1].lstrip('-').isdigit():
                idx = int(parts[-1])
                if 0 <= idx < len(self.parent.Data):
                    return idx
        except Exception:
            pass

        # Case 2: a category header is selected -> its displayed label is
        # "<1-based-index>: <title>" (see update_filelist). Parse the prefix.
        if isinstance(prop, wxpg.PropertyCategory):
            try:
                label = prop.GetLabel()
                prefix, _, _ = label.partition(':')
                prefix = prefix.strip()
                if prefix.isdigit():
                    idx = int(prefix) - 1
                    if 0 <= idx < len(self.parent.Data):
                        return idx
            except Exception:
                pass

        return None

    def on_move_data(self, direction):
        """Move the currently selected dataset by direction (-1 or +1)."""
        idx = self._get_selected_data_index()
        if idx is None:
            return
        new_idx = idx + direction
        if new_idx < 0 or new_idx >= len(self.parent.Data):
            return
        data_list = self.parent.Data
        data_list[idx], data_list[new_idx] = data_list[new_idx], data_list[idx]
        self.update_filelist()
        self.parent.matplotlib_panel.update_graph()
        try:
            new_prop = self.filetree.GetPropertyByName(f"show {new_idx}")
            if new_prop is not None:
                self.filetree.SelectProperty(new_prop)
                self.filetree.EnsureVisible(new_prop)
        except Exception:
            pass


    def on_doFFT(self, event):
        fftmethod = {}
        for kk in self.FFTparams.keys():
            prop = self.ffttree.GetProperty(kk)
            fftmethod[kk] = prop.GetValue()

        for ii, dd in enumerate(self.parent.Data):
            self.parent.Data[ii]['fftmethod'] = fftmethod
            self.parent.Data[ii]['fftactual'] = False # make sure to remember that things got changed.

        self.parent.update_FFT(fftmethod)
        
    def add_item_with_value(self, parent, name, value):
        item = self.filetree.AppendItem(parent, name)
        txt = wx.TextCtrl(self.filetree, value=value, style=wx.BORDER_NONE)
        self.filetree.SetItemWindow(item, txt)
        txt.Bind(wx.EVT_TEXT, lambda e, it=item: self.on_value_changed(it, e))
        return item    
    
    def on_value_changed(self, item, event):
        value = event.GetEventObject().GetValue()
        name = self.filetree.GetItemText(item)
        print(f"{name} updated to: {value}")
    def on_filetree_dclicked(self, event):
        if isinstance(event.GetProperty(), wxpg.PropertyCategory):
            self.filetree.CollapseAll()
            self.chk_expanded.SetValue(False)
        
    def on_filetree_clicked(self, event):
        if isinstance(event.GetProperty(), wxpg.PropertyCategory):
            cat_prop = event.GetProperty()

            if not self.chk_expanded.GetValue():
                self.filetree.CollapseAll()
                self.filetree.Expand(cat_prop)
                it = self.filetree.GetIterator(wxpg.PG_ITERATE_ALL)
                while not it.AtEnd():
                    prop = it.GetProperty()
                    if isinstance(prop, wxpg.PropertyCategory):
                        prop.GetCell(0).SetBgCol(wx.Colour(100, 150, 215))
                    it.Next()
                cat_prop.GetCell(0).SetBgCol(wx.Colour(64, 120, 215))

    def on_filetree(self, event):
        (key, strval) = event.GetPropertyName().split()
        val = int(strval)
        if key == 'bgdata':
            # An EnumProperty's value is the index into its choice list; the
            # dataset stores the background's title, which is what
            # update_FFT() matches against.
            self.parent.Data[val][key] = event.GetProperty().GetValueAsString()
        else:
            self.parent.Data[val][key] = event.GetValue()
        if self.parent.tb_autoFFT.GetValue():
            self.parent.update_FFT()   # redraws the figure itself at the end
        else:
            # f.min/f.max, z.min/z.max and Show change the figure without
            # needing a new transform, and with auto-FFT off nothing else
            # would redraw it -- which made the f.min/f.max spin arrows look
            # dead until something unrelated triggered a repaint.
            self.parent.matplotlib_panel.update_graph()

    def on_ffttree(self, event):
        fftmethod = {}
        for kk in self.FFTparams.keys():
            prop = self.ffttree.GetProperty(kk)
            fftmethod[kk] = prop.GetValue()

        for ii, dd in enumerate(self.parent.Data):
            self.parent.Data[ii]['fftmethod'] = fftmethod
            self.parent.Data[ii]['fftactual'] = False # make sure to remember that things got changed.

        if self.parent.tb_autoFFT.GetValue():
            self.parent.update_FFT(fftmethod)
    # ------------------------------------------------------------------
    # # Simulations tab (empty for now)
    # # ------------------------------------------------------------------
    # def _create_simulation_tab(self):
    #     self.sim_panel = wx.Panel(self.nb_book)
    #     self.nb_book.AddPage(self.sim_panel, "Simulations")

    #     placeholder = wx.StaticText(
    #         self.sim_panel,
    #         label="Simulation UI goes here.\nAdd buttons, charts, etc.",
    #     )
    #     placeholder.Wrap(250)

    #     sim_sizer = wx.BoxSizer(wx.VERTICAL)
    #     sim_sizer.AddStretchSpacer()
    #     sim_sizer.Add(placeholder, 0, wx.ALIGN_CENTER)
    #     sim_sizer.AddStretchSpacer()
    #     self.sim_panel.SetSizer(sim_sizer)
        
    def update_filelist(self):
        #print('tbd')
        self.filetree.Clear()
        cat = []
        bgchoices = ['none']
        for bg in self.parent.BGData:
            bgchoices.append(bg['title'])
            
        for ii, dd in enumerate(self.parent.Data):
            # Prefix the displayed title with the 1-based order number so the
            # user can see at a glance which column this dataset occupies in
            # the graph. The underlying dd['title'] is left unmodified.
            disp_title = f"{ii + 1}: {dd['title']}"
            cat.append(self.filetree.Append(wxpg.PropertyCategory(disp_title)))
            cat[ii].SetBackgroundColour(wx.Colour(255, 255, 255))
            cat[ii].SetTextColour(wx.Colour(0, 0, 180))
            self.filetree.AppendIn(cat[ii], wxpg.BoolProperty("Show", f"show {ii}",value=dd['show']))
            self.filetree.AppendIn(cat[ii], wxpg.StringProperty("f.name", f"fname {ii}",value=dd['fname']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("MW Freq [GHz]", f"freq {ii}", value=dd['freq']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("B0 [mT]", f"field {ii}", value=dd['field']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("tau [ns]", f"tau {ii}",value=dd['tau']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("z.max", f"zmax {ii}",value=dd['zmax']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("z.min", f"zmin {ii}",value=dd['zmin']))
            # f.min/f.max get the spin editor too: they are the analysis
            # window, walked in and out a step at a time while watching the
            # maps (with "View pad" to see what falls outside), so the arrows
            # and a reliable Enter matter more here than typing an exact number.
            fmax = self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("f.max [MHz]", f"fmax {ii}",value=dd['fmax']))
            self.filetree.SetPropertyEditor(fmax, "SpinFloat")
            fmax.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.1)
            fmin = self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("f.min [MHz]", f"fmin {ii}",value=dd['fmin']))
            self.filetree.SetPropertyEditor(fmin, "SpinFloat")
            fmin.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.1)
            # Width of the tilted analysis box, walked in and out like f.min
            # and f.max above and for the same reason, so it gets the same
            # editor. .get() rather than [] -- sessions saved before this
            # option existed have no key.
            adspan = self.filetree.AppendIn(
                cat[ii], wxpg.FloatProperty("adiag span [MHz]", f"antidiagspan {ii}",
                                            value=dd.get('antidiagspan', -1.0)))
            self.filetree.SetPropertyEditor(adspan, "SpinFloat")
            adspan.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.01)
            noiselvl = self.filetree.AppendIn(
                cat[ii], wxpg.FloatProperty("noise/max(y)", f"noiselevel {ii}",
                                            value=dd.get('noiselevel', -1.0)))
            self.filetree.SetPropertyEditor(noiselvl, "SpinFloat")
            noiselvl.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.01)
            noiselvl.SetAttribute(wxpg.PG_ATTR_MIN, -1.0)
            noiselvl.SetAttribute(wxpg.PG_ATTR_MAX, 1.0)
            noiselvl.SetHelpString(
                "Per-point noise level of this spectrum, as a fraction of the\n"
                "peak inside the analysis box: 0 to 1, the same peak z.max=-1\n"
                "scales by, so it is the floor height read straight off the\n"
                "Skyline and Diag.Projection traces. Converted to FFT(y)\n"
                "amplitude internally, so it reads the same on spectra of any\n"
                "absolute intensity, and a bright feature outside the box (which\n"
                "nothing here analyses) cannot move it.\n"
                " -1 = measure it, as before, off the far high-frequency sliver\n"
                "      where a HYSCORE spectrum carries no signal\n"
                "Any other value replaces that measurement -- use it when the\n"
                "sliver is not actually signal-free. The level a max projection\n"
                "of pure noise sits at (the floor drawn on those panels) is\n"
                "carried with it, keeping the max-to-mean ratio this spectrum\n"
                "measured.\n"
                "The grid search reads the same number for its chi-square, and,\n"
                "with 'Use exp. noise' on, fits the simulation on top of it.")
            adspan.SetHelpString(
                "How far off the f2=f1 diagonal the tilted analysis box reaches,\n"
                "in MHz of |f2-f1|. The box's length along the diagonal follows\n"
                "from it: the corners always stay on the f.min..f.max border, so\n"
                "a wider box is a shorter one.\n"
                " -1  = square box, half the window wide (corners at the middle\n"
                "       of each window edge)\n"
                " f.max-f.min = the widest the window holds -- the box collapses\n"
                "       to the bare anti-diagonal line, one sample thick\n"
                "Anything larger is taken at that maximum. This box is what the\n"
                "Diag.Projection maximizes over, what z.max scales by, and what\n"
                "the grid search fits.")
            # A single pulldown, like "Show" -- not a multi-choice. Only one
            # background can be subtracted anyway, and MultiChoiceProperty
            # handed back a *list*, which then never matched the plain title
            # strings this list is rebuilt from, so every rebuild silently
            # reset the selection to 'none'.
            bgname = dd['bgdata']
            if isinstance(bgname, (list, tuple)):   # sessions saved earlier
                bgname = bgname[0] if len(bgname) else 'none'
            if bgname not in bgchoices:
                bgname = 'none'
            self.parent.Data[ii]['bgdata'] = bgname
            self.filetree.AppendIn(cat[ii], wxpg.EnumProperty("BG Data", f"bgdata {ii}",
                                                              labels=bgchoices,
                                                              value=bgchoices.index(bgname)))
            bgscale = self.filetree.AppendIn(
                cat[ii], wxpg.FloatProperty("BG Scale", f"bgscale {ii}", value=dd['bgscale']))
            self.filetree.SetPropertyEditor(bgscale, "SpinFloat")
            bgscale.SetAttribute(wxpg.PG_ATTR_SPINCTRL_STEP, 0.01)
        if self.chk_expanded.GetValue():
            self.filetree.ExpandAll()
        else:
            self.filetree.CollapseAll()
            if cat:
                self.filetree.Expand(cat[-1])
            #self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("scale [ns]", f"tau{ii}",value=dd['tau']))

                
# --------------------------------------------------------------------------- #
# Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Data = []
        self.BGData = []
        self.SetTitle(f"pyHYSCORE {VERSION}")
        self.SetSize((1000, 600))
        self.Sys = sysPar()
        self.Exp = expPar()
        self.Opt = optHYSCORE()
        self.BackgroundColor = wx.Colour(255, 255, 255)
        self.ForegroundColor = wx.Colour(0, 0, 0)
        self.ControlBackgroundColor = wx.Colour(255, 255, 255)

        # ---- Colormap state (default: jet) ---------------------------------
        # Unrelated to the UI theme below: this colours the HYSCORE spectra
        # themselves, not the app's chrome.
        self.current_cmap_name  = "jet"
        self.current_cmap_stops = colormap_to_stops("jet", 9)   # list of (pos, hex)
        self.current_cmap       = stops_to_cmap(self.current_cmap_stops, "jet")

        # ---- UI theme (light/dark; pill buttons, checkboxes, tab bar) ------
        self._theme = "light"

        self.settings = self.get_default_settings() ####  replace with ini loader in the future

        # ---- Top icon toolbar ------------------------------------------------
        # Hand-drawn glyphs (no icon asset files needed) replacing what used
        # to be scattered pill buttons (file browser, Data tab) and the
        # Colormap/View menu bar. FileBrowserPanel and TabulatedPanel bind
        # their own handlers to tb_change_folder/tb_up/tb_save_session/
        # tb_doFFT/tb_doSIM once *they* exist (built after this toolbar);
        # tb_autoFFT/tb_autoSIM need no binding -- other code just reads
        # their .GetValue() on demand.
        self.top_toolbar = wx.Panel(self)
        toolbar_sizer = wx.BoxSizer(wx.HORIZONTAL)

        def add_sep():
            line = wx.StaticLine(self.top_toolbar, style=wx.LI_VERTICAL)
            toolbar_sizer.Add(line, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)

        self.tb_change_folder = theme.IconButton(
            self.top_toolbar, theme.draw_folder_icon, tooltip="Change folder…")
        self.tb_up = theme.IconButton(
            self.top_toolbar, theme.draw_up_icon, tooltip="Up")
        toolbar_sizer.Add(self.tb_change_folder, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.tb_up, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        add_sep()

        self.tb_save_session = theme.IconButton(
            self.top_toolbar, theme.draw_save_icon, tooltip="Save Session")
        toolbar_sizer.Add(self.tb_save_session, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        add_sep()

        self.tb_doFFT = theme.IconButton(
            self.top_toolbar, theme.draw_text_icon("FFT"), tooltip="Do FFT")
        self.tb_autoFFT = theme.IconToggleButton(
            self.top_toolbar, theme.draw_auto_toggle_icon("FFT"),
            tooltip="auto FFT (recompute automatically on change)", value=True)
        toolbar_sizer.Add(self.tb_doFFT, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.tb_autoFFT, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        add_sep()

        self.tb_doSIM = theme.IconButton(
            self.top_toolbar, theme.draw_text_icon("SIM"), tooltip="Do SIM")
        self.tb_autoSIM = theme.IconToggleButton(
            self.top_toolbar, theme.draw_auto_toggle_icon("SIM"),
            tooltip="auto SIM (re-simulate automatically on change)", value=True)
        toolbar_sizer.Add(self.tb_doSIM, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.tb_autoSIM, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        add_sep()

        self.tb_gridsearch = theme.IconButton(
            self.top_toolbar, theme.draw_gridsearch_icon, tooltip="Sys 2D Grid Search…")
        toolbar_sizer.Add(self.tb_gridsearch, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)

        toolbar_sizer.AddStretchSpacer(1)
        add_sep()

        # grid / density-contour toggles; state lives on MatplotlibPanel, handlers forward to it
        self.tb_grid = theme.IconToggleButton(
            self.top_toolbar, theme.draw_grid_icon, tooltip="Axes grid")
        self.tb_mapmode = theme.IconToggleButton(
            self.top_toolbar, theme.draw_map_mode_icon,
            tooltip="Density map (click for contours)")
        self.btn_cmap_icon = theme.IconButton(
            self.top_toolbar, theme.draw_colormap_icon, tooltip="Edit Colormap…")
        self.btn_theme_icon = theme.IconButton(
            self.top_toolbar, theme.draw_theme_icon, tooltip="Toggle dark mode (Ctrl+D)")
        toolbar_sizer.Add(self.tb_grid, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.tb_mapmode, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.btn_cmap_icon, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        toolbar_sizer.Add(self.btn_theme_icon, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 2)
        self.top_toolbar.SetSizer(toolbar_sizer)
        self.tb_grid.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_grid)
        self.tb_mapmode.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_mapmode)
        self.btn_cmap_icon.Bind(wx.EVT_BUTTON, self.on_edit_colormap)
        self.btn_theme_icon.Bind(wx.EVT_BUTTON, self.on_toggle_theme)

        # Ctrl+D still toggles the theme without needing a menu.
        theme_id = wx.NewIdRef()
        self.Bind(wx.EVT_MENU, self.on_toggle_theme, id=theme_id)
        self.SetAcceleratorTable(wx.AcceleratorTable([
            (wx.ACCEL_CTRL, ord('D'), theme_id),
        ]))

        # ---- Status strip --------------------------------------------------
        # Field 0: general status messages
        # Field 1: x, y coordinates when hovering over a matplotlib axis
        # Then the simulation progress gauge, hidden unless a run is going.
        # theme.StatusStrip rather than a native wx.StatusBar, whose field
        # text is drawn by MSW in the system colour and stayed black in dark
        # mode however the theme was set.
        self.statusbar = theme.StatusStrip(self, widths=[-1, 260])
        self.statusbar.SetStatusText("Ready", 0)
        self.statusbar.SetStatusText("", 1)

        self.sim_gauge = wx.Gauge(self.statusbar, range=100,
                                   style=wx.GA_HORIZONTAL | wx.GA_SMOOTH)
        self.statusbar.AddWidget(self.sim_gauge, width=170)
        self.sim_gauge.Hide()

        # Root splitter:  left | right
        self.splitter_main = wx.SplitterWindow(self, style=wx.SP_3D)
        self.splitter_main.SetMinimumPaneSize(20)     # keep left pane from getting too small

        # ---- Left pane ----------------------------------------------------
        self.file_browser = FileBrowserPanel(self.splitter_main, self)

        # ---- Right pane (will contain a *second* splitter) ---------------
        self.right_panel = wx.Panel(self.splitter_main)
        right_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # Second splitter (horizontal) inside the right panel
        self.splitter_right = wx.SplitterWindow(self.right_panel, style=wx.SP_3D)

        # Panels for the second splitter
        self.matplotlib_panel = MatplotlibPanel(self.splitter_right, self)
        self.tabulated_panel   = TabulatedPanel(self.splitter_right, self)

        # Split the right panel horizontally
        self.splitter_right.SplitVertically(
            self.matplotlib_panel, self.tabulated_panel,   # initial height of the top pane
        )
        self.splitter_right.SetSashGravity(1)  # proportion of space for the top pane
        self.splitter_right.SetSashPosition(-250)

        # Put the second splitter into the right panel's sizer
        right_sizer.Add(self.splitter_right, 1, wx.EXPAND)
        self.right_panel.SetSizer(right_sizer)

        # ---- Attach the two top level panes --------------------------------
        self._left_sash = 200      # remembered width to restore on expand
        self.splitter_main.SplitVertically(
            self.file_browser, self.right_panel,
            sashPosition=self._left_sash
        )
        self.splitter_main.SetSashGravity(0)  # proportion of space for left pane

        # ----------- Collapse toggle for the left panel ---------------------
        # Lives outside the splitter (not inside file_browser) so it stays
        # visible even once the panel it controls is fully hidden.
        self.toggle_btn = theme.IconButton(
            self, theme.draw_chevron_icon("left"), tooltip="Hide file browser",
            size=wx.Size(18, -1), fill_key="accent")
        self.toggle_btn.Bind(wx.EVT_BUTTON, self.on_toggle_left)

        # ----------- Collapse toggle for the right panel ---------------------
        # Lives outside the splitter (not inside right_panel) so it stays
        # visible even once the panel it controls is fully hidden.
        self.toggle_right_btn = theme.IconButton(
            self, theme.draw_chevron_icon("right"), tooltip="Hide tabulated panel",
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

        # Freeze the whole window for the duration of a live resize or
        # splitter-sash drag, instead of letting every intermediate tick
        # repaint the ~25+ owner-drawn pill/icon controls individually.
        # A sash drag alone doesn't fire EVT_SIZE on the frame (the frame's
        # own size is unchanged), so both splitters feed the same guard too.
        self._resize_freeze = theme.ResizeFreezeGuard(self)
        self.Bind(wx.EVT_SIZE, self._resize_freeze.on_event)
        self.splitter_main.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self._resize_freeze.on_event)
        self.splitter_right.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGING, self._resize_freeze.on_event)

        self.Centre()
        self.apply_theme()
    # -----------------------------------------------------------------------
    def on_toggle_theme(self, event):
        self._theme = "dark" if self._theme == "light" else "light"
        self.apply_theme()

    def on_toggle_left(self, event):
        """Actually hide/show the file-browser pane (Unsplit/SplitVertically)
        rather than just shrinking its sash to the minimum -- the toggle
        button itself stays visible either way (see __init__)."""
        if self.splitter_main.IsSplit():
            self._left_sash = self.splitter_main.GetSashPosition()
            self.splitter_main.Unsplit(self.file_browser)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("right"))
            self.toggle_btn.SetToolTip("Show file browser")
        else:
            self.splitter_main.SplitVertically(
                self.file_browser, self.right_panel, sashPosition=self._left_sash)
            self.toggle_btn.SetDrawIcon(theme.draw_chevron_icon("left"))
            self.toggle_btn.SetToolTip("Hide file browser")

    def on_toggle_right(self, event):
        if self.splitter_right.IsSplit():
            self._right_sash = self.splitter_right.GetSashPosition()
            self.splitter_right.Unsplit(self.tabulated_panel)
            self.toggle_right_btn.SetDrawIcon(theme.draw_chevron_icon("left"))
            self.toggle_right_btn.SetToolTip("Show tabulated panel")
        else:
            self.splitter_right.SplitVertically(
                self.matplotlib_panel, self.tabulated_panel, sashPosition=self._right_sash)
            self.toggle_right_btn.SetDrawIcon(theme.draw_chevron_icon("right"))
            self.toggle_right_btn.SetToolTip("Hide right panel")

    def apply_theme(self):
        """Push self._theme through every themed surface: native controls
        (via theme.apply_theme_to_window, which also Refresh()es every pill
        widget so it repaints in the new colours), the property grids (not
        reached by that generic walk -- they need a caption-colour choice),
        and the matplotlib figures (separate from the per-file data
        colormap, which is untouched here)."""
        theme.set_theme(self._theme)
        theme.apply_theme_to_window(self, self._theme)

        tp = self.tabulated_panel
        tp.Color_BG_FILE_Title = wx.Colour(theme.theme_colors()["accent"])
        tp.Color_BG_FILE_Inact = wx.Colour(theme.theme_colors()["tab_inactive"])
        tp.Color_FG_FILE_Title = wx.Colour(theme.theme_colors()["propgrid_caption_fg"])
        tp.Color_BG_FILE_Main = wx.Colour(theme.theme_colors()["ctrl_bg"])
        tp.Color_FG_FILE_Main = wx.Colour(theme.theme_colors()["text"])
        tp.Color_BG_FFT_Title = wx.Colour(theme.theme_colors()["pill_load"])
        tp.Color_BG_FFT_Inact = wx.Colour(theme.theme_colors()["tab_inactive"])
        tp.Color_FG_FFT_Title = wx.Colour(theme.theme_colors()["propgrid_caption_fg"])
        tp.Color_BG_FFT_Main = wx.Colour(theme.theme_colors()["ctrl_bg"])
        tp.Color_FG_FFT_Main = wx.Colour(theme.theme_colors()["text"])
        theme.theme_propgrid(tp.filetree, caption_key="accent")
        theme.theme_propgrid(tp.ffttree, caption_key="pill_load")
        theme.theme_propgrid(tp.Sys_param.pg, caption_key="accent")
        theme.theme_propgrid(tp.Exp_param.pg, caption_key="accent")
        theme.theme_propgrid(tp.Opt_param.pg, caption_key="accent")
        tp.update_filelist()

        theme.style_figure(self.matplotlib_panel.figure, self.matplotlib_panel.axes)
        self.matplotlib_panel.canvas.draw_idle()

        self.Refresh()

    # -----------------------------------------------------------------------
    def on_edit_colormap(self, event):
        """Open the ColormapEditorDialog; apply the chosen colormap if OK."""
        with ColormapEditorDialog(
                parent=self,
                preset_stops=self.current_cmap_stops,
                preset_name=self.current_cmap_name,
                title="Colormap Editor") as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                cmap  = dlg.get_colormap()
                name  = dlg.get_colormap_name()
                stops = dlg.get_stops()
                if cmap is not None:
                    self.current_cmap       = cmap
                    self.current_cmap_name  = name
                    self.current_cmap_stops = stops
                    self.matplotlib_panel.update_graph()

    def on_toggle_mapmode(self, event):
        self.matplotlib_panel.set_map_mode(
            'contour' if self.tb_mapmode.GetValue() else 'density')

    def on_toggle_grid(self, event):
        self.matplotlib_panel.set_show_grid(self.tb_grid.GetValue())

    def get_view_state(self):
        mp = self.matplotlib_panel
        return {'mapMode': mp.mapMode, 'contourLevels': int(mp.contourLevels),
                'contourWidth': float(mp.contourWidth), 'overlayColor': str(mp.overlayColor),
                'showGrid': bool(mp.showGrid)}

    def set_view_state(self, view):
        """Apply a saved 'view' dict; every key optional and validated, bad ones keep the current value."""
        mp = self.matplotlib_panel
        tp = self.tabulated_panel
        if not isinstance(view, dict):
            view = {}
        mode = view.get('mapMode', mp.mapMode)
        mp.mapMode = mode if mode in ('density', 'contour') else 'density'
        try:
            mp.contourLevels = int(np.clip(int(view.get('contourLevels', mp.contourLevels)), 2, 100))
        except (TypeError, ValueError):
            pass
        try:
            mp.contourWidth = float(np.clip(float(view.get('contourWidth', mp.contourWidth)), 0.1, 5.0))
        except (TypeError, ValueError):
            pass
        colour = wx.Colour(str(view.get('overlayColor', mp.overlayColor)))
        if colour.IsOk():
            mp.overlayColor = colour.GetAsString(wx.C2S_HTML_SYNTAX)
        mp.showGrid = bool(view.get('showGrid', mp.showGrid))

        self.tb_mapmode.SetValue(mp.mapMode == 'contour')
        self.tb_mapmode.SetToolTip("Contour map (click for density)" if mp.mapMode == 'contour'
                                   else "Density map (click for contours)")
        self.tb_grid.SetValue(mp.showGrid)
        tp.spin_contourLevels.SetValue(mp.contourLevels)
        tp.spin_contourWidth.SetValue(mp.contourWidth)
        tp.pick_overlayColor.SetColour(wx.Colour(mp.overlayColor))
        tp.sync_contour_controls()
        mp.rebuildAxes = True



    def RaiseError(self, message):
        # Show a modal error dialog with a custom message
        dlg = wx.MessageDialog(
            None,
            message,
            "Error",
            wx.OK | wx.ICON_MASK
        )
        dlg.ShowModal()              # the call blocks until the user presses OK
        dlg.Destroy()                # free dialog resources

        raise Exception(message)
    # -----------------------------------------------------------------------
    def RaiseMessage(self, message):
        # Show a modal error dialog with a custom message
        dlg = wx.MessageDialog(
            None,
            message,
            "Message",
            wx.OK | wx.ICON_MASK
        )
        dlg.ShowModal()              # the call blocks until the user presses OK
        dlg.Destroy()                # free dialog resources
    def get_default_settings(self):
        settings = {
            'mtl.showFFT':(True, 'show FFT'),
            'mtl.showSim':(True, 'show Sim'),
            'mtl.showBG': (False, 'show BG'),
            'mtl.showSkyline': (False, 'show Skyline'),
            'mtl.showOriSel': (False, 'show Ori.Sel.'),
            'mtl.showDiagonalProj': (False, 'show DiagonalProj'),
            'mtl.OverlaySim' : (False, 'Overlay Sim'),
            'mtl.title': (['field', ['field', 'gvalue', 'file', 'none']], 'Title'),
            'mtl.cmap': (['jet', ['jet', 'fall']], 'Gradient Color'),
            'mtl.showLabels':(True, 'show Labels'),
            'mtl.overlay_nLev':(7, '# of overlay levels'),
            'mtl.overlay_color':(['white', ['white', 'yellow', 'red', 'green', 'blue', 'black']], 'Overlay Color'),
            'tab.autoFFT':(True, 'auto FFT'),
            'tab.autoSim':(True, 'auto Sim'),
            'tab.Color_BG_FILE_Title':(wx.Colour(64, 120, 215), 'Color_BG_FILE_Title'),
            'tab.Color_BG_FILE_Inact':(wx.Colour(64, 120, 215), 'Color_BG_FILE_Inact'),
            'tab.Color_FG_FILE_Title':( wx.Colour(255, 255, 255), 'Color_FG_FILE_Title'),
            'tab.Color_BG_FILE_Main':(wx.Colour(255, 255, 255), 'Color_BG_FILE_Main'),
            'tab.Color_FG_FILE_Main':(wx.Colour(0, 0, 0), 'Color_FG_FILE_Main'),
            'tab.Color_BG_FFT_Title':(wx.Colour(215, 120, 64), 'Color_BG_FFT_Title'),
            'tab.Color_BG_FFT_Inact':(wx.Colour(100, 150, 215), 'Color_BG_FFT_Inact'),
            'tab.Color_FG_FFT_Title':(wx.Colour(255, 255, 255), 'Color_FG_FFT_Title'),
            'tab.Color_BG_FFT_Main':(wx.Colour(255, 255, 255), 'Color_BG_FFT_Main'),
            'tab.Color_FG_FFT_Main':(wx.Colour(0, 0, 0), 'Color_FG_FFT_Main'),
            'brw.current_filter':(["DSC",  ["All", "DSC", "XML"]], 'Default Filter'),
            'ma.ForegroundColor':(wx.Colour(0, 0, 0), 'Foreground Color'),
            'ma.ControlBackgroundColor':(wx.Colour(255, 255, 255), 'Control Background Color'),
            }
        return settings
    
    def set_default_settings(self):
        settings = self.get_default_settings()
        self.update_setting(settings)

    def update_setting(self, settings):
        for compkey in settings.keys():
            if 'mtl.' in compkey:
                key =compkey.split('mtl.')[1]
                if hasattr(self.matplotlib_panel, key):
                    setattr(self.matplotlib_panel, key, settings[compkey][0])
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for matplotlib_panel")
            if 'tab.' in compkey:
                key =compkey.split('tab.')[1]
                if hasattr(self.tabulated_panel, key):
                    setattr(self.tabulated_panel, key, settings[compkey][0])
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for tabulated_panel")
            if 'brw.' in compkey:
                key =compkey.split('brw.')[1]
                if hasattr(self.file_browser, key):
                    setattr(self.file_browser, key, settings[compkey][0])
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for file_browser")
            if 'man.' in compkey:
                key =compkey.split('man.')[1]
                if hasattr(self, key):
                    setattr(self, key, settings[compkey][0])
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for main window")
    def collect_current_settings(self):
        settings = self.get_default_settings() # just so we have all the parameters
    #
        for compkey in settings.keys():
            if 'mtl.' in compkey:
                key =compkey.split('mtl.')[1]
                if hasattr(self.matplotlib_panel, key):
                    settings[compkey][0] =getattr(self.matplotlib_panel, key)
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for matplotlib_panel")
            if 'tab.' in compkey:
                key =compkey.split('tab.')[1]
                if hasattr(self.tabulated_panel, key):
                    settings[compkey][0] =getattr(self.tabulated_panel, key)
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for tabulated_panel")
            if 'brw.' in compkey:
                key =compkey.split('brw.')[1]
                if hasattr(self.file_browser, key):
                    settings[compkey][0] =getattr(self.file_browser, key)
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for file_browser")
            if 'man.' in compkey:
                key =compkey.split('man.')[1]
                if hasattr(self, key):
                    settings[compkey][0] =getattr(self, key)
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for main window")
        self.settings = settings

    def update_settings_tab(self):
        self.collect_current_settings()
        ######### TBD #################

    def loadData(self, path, BG=False):
        name, ext = os.path.splitext(path) 
        _, fname = os.path.split(path) 
        if ext.lower() in ('.dsc', '.dta'):
            ax, data, dsc = Mybr.brukerread(path, return_ax=True, return_dsc=True)
            freq = ax['freq1']*1e-9
            tau=0
            field=0
            
            title = ax['title']
            for ss in title.split():
                if 'FTEzDelay1' in dsc.keys():
                    res = dsc['FTEzDelay1'].split('ns')
                    try:
                        tau = int(res[0])
                    except Exception:
                        tau = 0
                elif ('tau' in ss)|('ns' in ss):
                    try:
                        flt = filter(str.isnumeric,ss)
                        tau=int("".join(flt))
                    except Exception:
                        tau = 0
                    continue
                if 'G' in ss:
                    try:
                        flt = filter(str.isnumeric,ss)
                        field=float("".join(flt))/10
                    except Exception:
                        field = 0
                    continue
                if 'mT' in ss:
                    try:
                        flt = filter(str.isnumeric,ss)
                        field=float("".join(flt))
                    except Exception:
                        field = 0
                    continue
                if ss.startswith('9'):
                    ss[0]
            # try to identify the units
            if ('ns' in ax['xlabel']):
                ax['x'] *= 1e-3
                ax['y'] *= 1e-3
                ax['xlabel']='Time, us'
                ax['ylabel']='Time, us'
                
            isfft = False
            
            if 'fft' in fname.lower():
                isfft = True
            if 'fft' in title.lower(): 
                isfft = True
            if 'DSCR' in dsc.keys():
                if dsc['DSRC'] == 'MAN' and dsc['PROCESS'] == 'prAbs':
                    isfft = True
                            
            if isfft:
                fmin = ax['x'][0]
                fmax = ax['x'][1]
            else:
                dx = np.abs(ax['x'][1] - ax['x'][0])
                fmin = -1/2/dx
                fmax = 1/2/dx
                    
            tData = {'ax':ax, 'data':data, 'dsc':dsc,
                     'bgdata':'none', 'bgscale':1.0,    
                    'freq':freq, 'tau':tau, 'field':field,
                    'fname':fname, 'fullpath':path, 'show':True, 
                    'zmax':1.0, 'zmin':0.0, 'title':ax['title'],
                    'fmax':fmax, 'fmin':fmin,
                    # How far off the f2=f1 diagonal the tilted analysis box
                    # reaches, in MHz of |f2-f1|; -1 = the square box, half the
                    # f.min..f.max window wide. See optHYSCORE.diagonalBox.
                    'antidiagspan':-1.0,
                    'isfft':False,
                    'fftdata':None, 'fftax':None, 'fftmethod':None, 'fftactual':False,
                    # measured by update_FFT once there is a spectrum to
                    # measure: per-point sigma, and the level a max projection
                    # of pure noise sits at
                    'noise':0.0, 'noisesky':0.0,
                    # -1 = use the two measured above; 0..1 replaces the
                    # per-point level, as a fraction of the peak inside the
                    # analysis box, and carries the skyline maximum with it.
                    # See MainFrame.dataset_noise.
                    'noiselevel':-1.0, # Noize level from max
                    'simdata':None, 'simax':None, 'simmethod':None, 'simactual':False,
                    'orisel':None,  'region':None
                    }
                    
            if BG:
                self.BGData.append(tData)
            else:
                self.Data.append(tData)
        self.update_everything()        
        
    def update_everything(self):
        self.tabulated_panel.update_filelist()
        self.matplotlib_panel.update_graph()
            

    def update_FFT(self, fftmethod=None):
        # Build a fallback fftmethod from the UI controls so we always have
        # valid parameters even when called without an explicit argument.
        ui_fftmethod = {}
        for kk in self.tabulated_panel.FFTparams.keys():
            prop = self.tabulated_panel.ffttree.GetProperty(kk)
            ui_fftmethod[kk] = prop.GetValue()

        for ii, dd in enumerate(self.Data):
            # Per-dataset fftmethod: prefer the explicit argument, then the
            # value stored on the dataset, then fall back to the UI controls.
            if fftmethod:
                method = fftmethod
            elif dd['fftmethod'] is not None:
                method = dd['fftmethod'].copy()
            else:
                method = ui_fftmethod.copy()

            if len(dd['data'].shape)>2: ### need to fix at the source. Just for now let's deal with it.
                Data = np.array(dd['data'][:,:,0])
            else:
                Data = np.array(dd['data'])
                
            if dd['isfft']:
                # No transform to do, but the noise level still has to be
                # measured -- these datasets arrive as spectra already, so
                # their own data/axes are what the border is read off.
                nn = len(dd['ax']['x'])
                # the sliver is a slice of 2D Data at high frequency
                sliver = np.real(Data[0:int(nn/4), 0])
                self.Data[ii]['noisesky'] = np.max(sliver)
                self.Data[ii]['noise'] = np.mean(sliver)

                continue          # skip already-FFT datasets; do NOT return

            # dd['bgdata'] holds the chosen background's title, or 'none'.
            # Sessions saved while it was a MultiChoice property carry a
            # one-element list instead.
            bgname = dd['bgdata']
            if isinstance(bgname, (list, tuple)):
                bgname = bgname[0] if len(bgname) else 'none'
            if bgname!='none':
                BGscale= dd['bgscale']

                gotBG = False
                for bb in self.BGData:
                    if bgname==bb['title']:
                        if len(bb['data'].shape)>2: ### need to fix at the source. Just for now let's deal with it.
                            BG = np.array(bb['data'][:,:,0])
                        else:
                            BG = np.array(bb['data'])
                        gotBG = True
                        break
                if gotBG:
                    if Data.shape !=BG.shape:
                        # raise ValueError("😭 Background and Data have different dimensions.")
                        self.RaiseError("😭 Background and Data have different dimensions.")
                    Data -=BGscale*BG
                
            if not method['imag']:
                Data=np.real(Data)
                
            X=np.array(dd['ax']['x'])
            Y=np.array(dd['ax']['y'])
                                       
            for cc in range(Data.shape[0]):
                pp = np.polyfit(Y, Data[cc,:], method['polynomial'])
                Data[cc,:]-=np.polyval(pp, Y)
            for cc in range(Data.shape[1]):
                pp = np.polyfit(X, Data[:, cc], method['polynomial'])
                Data[:, cc]-=np.polyval(pp, X)
            
            if method['apodization'][0]!='none':
                tX, tY = np.meshgrid(X,Y)
                awidth= method['awidth']
                ax = np.max(X)*awidth
                ay = np.max(Y)*awidth
                if method['apodization'][0]=='Hamming':
                    app = (27/50+23/50*np.cos(np.pi*tX/ax) )*(27/50+23/50*np.cos(np.pi*tY/ay) )
                elif method['apodization'][0]=='Gaussian':
                    # sig = ax/2.0/2.354820045 # 2.354820045 = 2sqrt(2ln2)
                    sig2 = (ax/2.354820045)**2
                    app = np.exp(-( tX**2 + tY**2)/sig2 )
                elif method['apodization'][0]=='Lor-Gau':
                    alpha = method['aalpha']
                    shift = method['ashift']
                    sig2 = (ax/2.354820045)**2 # 2.354820045 = 2sqrt(2ln2)
                    # app = ( np.exp(-0.5/sig**2* ( (tX-shift)**2 + (tY-shift)**2) )*
                    #         np.exp(tX/ax/alpha + tY/ay/alpha) 
                    #        )
                    app = ( np.exp(-( (tX-shift)**2 + (tY-shift)**2)/sig2 
                                   + tX/ax/alpha + tY/ay/alpha)
                           )
                Data*=app
                
            if method['zfill']>0:
                ndim = np.ceil(np.log2(Data.shape[0]))+method['zfill']
                if ndim>13:
                    # raise AttributeError("Zero Filling is set too high.")
                    self.RaiseError("😭 Zero Filling is set too high.")
                padn = int(2**ndim-Data.shape[0])
                Data=np.pad(Data, (0,padn))
            
            fftData = np.abs(np.fft.fftshift(np.fft.fft2(Data)))
            nyqX = 1/2/(X[2]-X[1])
            nyqY = 1/2/(Y[2]-Y[1])
            
            # endpoint=False, because the transform's top bin sits one step
            # *below* +Nyquist and never at it: the bins are spaced
            # 2*nyq/N, not 2*nyq/(N-1). Including the endpoint stretched the
            # whole axis by N/(N-1), reading every frequency high by that
            # factor -- 0.39% at N=256, 0.05% at N=2048. Since the factor
            # depends on N, changing the zero fill visibly moved the peaks,
            # and the axis also disagreed with the simulation's own
            # (HYSCOREsim builds X/Y with endpoint=False already).
            # This is exactly np.fft.fftshift(np.fft.fftfreq(N, dt)).
            fftX = np.linspace(-nyqX, nyqX, fftData.shape[0], endpoint=False)
            fftY = np.linspace(-nyqY, nyqY, fftData.shape[1], endpoint=False)
            self.Data[ii]['fftdata']=fftData
            self.Data[ii]['fftax']={'x':fftX, 'y':fftY, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz'}
            self.Data[ii]['fftactual'] = True
            # Measured on the spectrum this transform just produced, so the
            # level always matches the apodization/zero-fill/background
            # settings actually in force rather than some earlier pass's.
            # Stored on the dataset, which means it rides along into the
            # session XML (save_mixed_dict_xml serializes the Data dicts
            # wholesale) and is available to the grid search without a rescan.
            
            nn = fftData.shape[0]
            self.Data[ii]['noisesky'] = np.max(fftData[0:int(nn/4), 0])
            self.Data[ii]['noise']    = np.mean(fftData[0:int(nn/4), 0])

        self.matplotlib_panel.update_graph()
    def on_update_orisel(self, event):
        if self.Opt.OriSelType[0]=='precalculated':

            self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
            hs = HYSCOREsim(Sys=self.Sys, errorFunc=self.RaiseError)
            if len(self.Data)>0:
                for ii,dd in enumerate(self.Data):
                    self.actuateSimMethod(hs, self.Data[ii])
                    # honour the Grid choice on the Opt tab; this used to be
                    # hardwired to 'fibonacci' and ignored the selection
                    matr = hs.computeOrisel_eig(grid=hs.grid_name(), nKnots=hs.Opt.nKnots, epsilon=0.33, returnMatrix=True, setActive=False)
                    self.Data[ii]['orisel'] = matr.copy()
        else:
            for ii,dd in enumerate(self.Data):
                self.Data[ii]['orisel'] = None

        # if self.tb_autoSIM.GetValue():
        #     self.runSim()

    def actuateSimMethod(self, hs, data):
        self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
        hs.Sys=self.Sys
        methoddic = {'Sys': self.tabulated_panel.Sys_param.parameters.copy(), 
                     'Exp': self.tabulated_panel.Exp_param.parameters.copy(),
                     'Opt': self.tabulated_panel.Opt_param.parameters.copy(),}
        if methoddic['Exp']['tau']<0:
            methoddic['Exp']['tau'] = data['tau']
        if methoddic['Exp']['Field']<0:
            methoddic['Exp']['Field'] = data['field']
        if methoddic['Exp']['mwFreq']<0:
            methoddic['Exp']['mwFreq'] = data['freq']
        if methoddic['Exp']['MaxFreq']<0:
            methoddic['Exp']['MaxFreq'] = data['fmax']
        if methoddic['Exp']['nPoints']<0:
            methoddic['Exp']['nPoints'] = len(data['fftax']['x'])
        
        if type(data['orisel'])!=type(None):
            methoddic['Opt']['OriSelInp'] = data['orisel']
        
        data['simmethod']=methoddic
        data['simactual']=False
    
        self.Opt.setFromCtrl(data['simmethod']['Opt'])
        self.Exp.setFromCtrl(data['simmethod']['Exp'])
        #print(self.Exp.Field)
        hs.Opt = self.Opt
        hs.Exp = self.Exp

    # ------------------------------------------------------------------
    # Simulation progress gauge (lives in status bar field 2)
    # ------------------------------------------------------------------
    def _sim_progress_begin(self, total):
        """Show an empty gauge sized to the number of spectra about to run."""
        if total <= 0:
            return
        self.sim_gauge.SetRange(total)
        self.sim_gauge.SetValue(0)
        self.sim_gauge.Show()
        self.statusbar.Layout()
        self.statusbar.SetStatusText(f"Simulating 0/{total}...", 0)
        self.statusbar.Update()

    def _sim_progress_step(self, done, total):
        """One spectrum finished. Repaint the gauge directly rather than
        yielding to the event loop: the simulation runs on the main thread,
        so a Yield here would let the user re-enter Run (or close the
        window) half way through the series."""
        self.sim_gauge.SetValue(min(done, total))
        self.statusbar.SetStatusText(f"Simulating {done}/{total}...", 0)
        self.sim_gauge.Update()
        self.statusbar.Update()

    def _sim_progress_end(self):
        self.sim_gauge.SetValue(0)
        self.sim_gauge.Hide()
        self.statusbar.Layout()

    def runSim(self):
        self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
        
        hs = HYSCOREsim(Sys=self.Sys, errorFunc=self.RaiseError)
        hs.preCompute() ## get housekeeping stuff out of the way to speed up computations a bit
        hs.verbose = False
        failed_titles = []
        total = len(self.Data)
        self._sim_progress_begin(total)
        
        end_time = 0.0
        try:
            start_time = time.time()
            self._runSim_series(hs, failed_titles, total)
            end_time = time.time()-start_time
        finally:
            # However the series ends -- cleanly, or on an exception that got
            # past the per-dataset handler -- the gauge must not be left
            # sitting on screen at a stale value.
            self._sim_progress_end()
        
        if failed_titles:
            self.statusbar.SetStatusText(
                f"Simulation failed for {len(failed_titles)}/{len(self.Data)} dataset(s)", 0)
        else:
            self.statusbar.SetStatusText(f"Ready. It took {end_time:.2f} s", 0)

        self.matplotlib_panel.update_graph()

    def _runSim_series(self, hs, failed_titles, total):
        """The per-dataset simulation loop, split out of runSim so the
        progress gauge can be cleaned up in a finally block around it."""
        for ii,dd in enumerate(self.Data):
            self.actuateSimMethod(hs, self.Data[ii])

            try:
                hs.reRun()
                
                self.Data[ii]['simdata']=hs.Spectrum
                self.Data[ii]['simax']={'x':hs.X, 'y':hs.Y, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz',
                                        'orisel': np.array([hs.phi, hs.theta, hs.ak]),
                                        # knot triangulation from make_grid; None
                                        # for a precalculated grid, in which case
                                        # the plot falls back to its own
                                        'oriseltri': None if hs.tri is None else hs.tri.copy()}
                self.Data[ii]['simactual']=True
            except Exception as e:
                # errorFunc (RaiseError) already showed a dialog for errors it
                # anticipated (e.g. "no resonances"); this catch-all also has
                # to handle genuinely unexpected exceptions, which otherwise
                # would vanish here with zero feedback. Always print so there
                # is at least a paper trail on the command line.
                print(f"[pyHYSCORE] Simulation failed for '{dd.get('title', dd.get('fname', '?'))}': "
                      f"{type(e).__name__}: {e}")
                failed_titles.append(dd.get('title', dd.get('fname', f'dataset {ii}')))
                self.Data[ii]['simdata']=None
                self.Data[ii]['simax']={'x':hs.X, 'y':hs.Y, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz',
                                        'orisel': np.array([[0], [0], [0]]),
                                        'oriseltri': None}
                self.Data[ii]['simactual']=True

            # Step whether the spectrum succeeded or failed -- the bar
            # reports how far through the series we are, not how many
            # worked, and it must still reach the end if some fail.
            self._sim_progress_step(ii + 1, total)
    def calcSSR(self, Exp, Sim, Mask, offset):
        #Data and Sim must of the same XY
        # Lives in mathfunctions so a grid-search worker process can score a
        # simulation without importing this module (and with it wx).
        return mf.calcSSR(Exp, Sim, Mask, offset)
    def diagonalBox(self, axis0, axis1, fmin, fmax, antidiagspan=-1.0):
        """Samples of a 2D spectrum lying inside the 45-degree tilted rectangle
        inscribed in the [fmin, fmax] x [fmin, fmax] window
        """
        axis0, axis1 = np.asarray(axis0), np.asarray(axis1)
        sel0 = (axis0 >= fmin) & (axis0 <= fmax)
        sel1 = (axis1 >= fmin) & (axis1 <= fmax)
        if not sel0.any() or not sel1.any():
            return None
        off0, n0 = int(np.flatnonzero(sel0)[0]), int(sel0.sum())
        off1, n1 = int(np.flatnonzero(sel1)[0]), int(sel1.sum())

        # Samples available on diagonal d = i-j, for i/j inside the window.
        # Non-increasing in |d|, so the whole span |d| <= D is covered as soon
        # as its two outermost diagonals are.
        count = lambda d: np.minimum(n0-1, n1-1+d) - np.maximum(0, d) + 1

        # f2-f1 is d steps, so the requested span in MHz is a half-range in d
        # directly. Clamped to the widest the window holds, which is d = nwin-1
        # (the corner-to-corner anti-diagonal) and leaves L == 1 below.
        nwin = min(n0, n1)
        step = abs(axis0[1]-axis0[0]) if axis0.size > 1 else 0.0
        span = float(antidiagspan) if antidiagspan is not None else -1.0
        if span >= 0 and step > 0:
            D = int(round(span/step))
        else:
            D = (nwin-1)//2          # -1: the square box, half the window wide
        D = max(0, min(D, nwin-1))

        # The corners sit on the window's border, so the diagonal extent is the
        # longest run every diagonal in |d| <= D can still hold -- decided by
        # the two outermost ones, since count only falls with |d|.
        L = int(min(count(-D), count(D)))
        if L < 1:
            return None

        d = np.arange(-D, D+1)
        # Centred run: its midpoint then lands on the window's anti-diagonal
        # for every d, which is what makes the selection a rectangle rather
        # than a ragged strip.
        i = (np.maximum(0, d) + (count(d)-L)//2)[:, None] + np.arange(L)[None, :]
        j = i - d[:, None]

        mask = np.zeros((n0, n1), dtype=bool)
        mask[i, j] = True

        # Outline the rectangle the selected samples fill, rather than four of
        # the samples themselves. The run ends alternate between two
        # neighbouring anti-diagonals (the parity above), so a line drawn
        # through any particular sample lands a whole step to one side of half
        # the boundary; taking the midpoint of the alternation instead centres
        # it, leaving every selected sample within 0.354 of a sample of the
        # line. Walked as a closed ring, so the four are in polygon order.
        s = i + j
        s_lo = 0.5*(s[:, 0].min() + s[:, 0].max())
        s_hi = 0.5*(s[:, -1].min() + s[:, -1].max())
        step0 = (axis0[1]-axis0[0]) if axis0.size > 1 else 0.0
        step1 = (axis1[1]-axis1[0]) if axis1.size > 1 else 0.0
        corners = np.array([(axis0[off0] + 0.5*(sv+dv)*step0,
                             axis1[off1] + 0.5*(sv-dv)*step1)
                            for sv, dv in ((s_lo, -D), (s_hi, -D),
                                           (s_hi, D), (s_lo, D))])
        return {'sel0': sel0, 'sel1': sel1, 'mask': mask,
                'idx0': off0+i, 'idx1': off1+j, 'corners': corners,
                'ndiag': int(d.size),      # points in the projection
                'nalong': L,               # samples behind each of them
                'nwindow': min(n0, n1)}    # longest diagonal the window holds
    
    def dataset_region(self, dd, fit_noise):
        """Precalculates dataset's experimental spectrum restricted to its
        [fmin, fmax] x [fmin, fmax] region
        """
        if dd['isfft']:
            raw = np.asarray(dd['data'])
            expdata = np.real(raw[:, :, 0]) if raw.ndim > 2 else np.real(raw)
            expX = np.asarray(dd['ax']['x'])
            expY = np.asarray(dd['ax']['y'])
        else:
            if not dd['fftactual'] or dd['fftdata'] is None:
                dd['region']=None
                return None
            expdata = dd['fftdata']
            expX = np.asarray(dd['fftax']['x'])
            expY = np.asarray(dd['fftax']['y'])

        fmin, fmax = dd['fmin'], dd['fmax']
        # Array axis 0 tracks the 'x' coordinate array and axis 1 tracks 'y'
        box = self.diagonalBox(expX, expY, fmin, fmax,
                                        dd.get('antidiagspan', -1.0))
        if box is None:
            return None
        xsel, ysel, mask = box['sel0'], box['sel1'], box['mask']
        n_fit = int(mask.sum())
        if not n_fit:
            return None
        # Define region:
        exp_region = expdata[np.ix_(xsel, ysel)]
        # 'null' is this dataset's no-simulation sum of squares:
        null = float(np.sum(exp_region[mask] ** 2))
        # 'sigma' is this spectrum's per-point noise level
        # Eitehr measured on its far high-frequency corner
        # or set by "Noise">0 (dd["noiselevel"] behind the scene) 
        if dd.get('noiselevel', -1.0) <0:
            sigma = dd['noise']
        else:
            sigma = np.max(expdata)*dd['noiselevel']

        offset = float(sigma) if (fit_noise and sigma) else 0.0

        
        region =  {'exp': exp_region,
                'mask': mask,
                'n': n_fit,
                'corners': box['corners'],
                'x': expX[xsel],
                'y': expY[ysel],
                'null': null,
                'null_rms': float(np.sqrt(null / n_fit)),
                'sigma': float(sigma) if sigma else None,
                'offset': offset,
                }
        dd['region'] = region.copy()
        return region
    def adaptSim(self, Data_x, Data_y, Sim, SimX, SimY):
        # See calcSSR above -- kept in mathfunctions for the worker's sake.
        return mf.adaptSim(Data_x, Data_y, Sim, SimX, SimY)
    
    def calcRMSD(self, norm=True):
        ### TBD; get all RMSD calculations to main pyHYSCORE from grid_search
        rmsd = 0
        for dd in Data:
            if dd['simdata'] is None:
                continue
            if dd['fftdata'] is None:
                if dd['isfft'] is False:
                    continue
                else:
                    data = np.real(dd['data'][:,:,0])
                    axEx = [np.min(dd['ax']['x']), np.max(dd['ax']['x']),
                            np.min(dd['ax']['y']), np.max(dd['ax']['y'])]
            else:
                data = dd['fftdata'][:,:,0]
                axEx = [np.min(dd['fftax']['x']), np.max(dd['fftax']['x']),
                        np.min(dd['fftax']['y']), np.max(dd['fftax']['y'])]

            xfi,xfa, yfi, yfa = (dd['fmin'], dd['fmax'], dd['fmin'], dd['fmax'] )

            adspan = dd.get('antidiagspan', -1.0)
            x_axis = np.linspace(axEx[0], axEx[1], data.shape[1])
            y_axis = np.linspace(axEx[2], axEx[3], data.shape[0])
            x_fit = (x_axis >= xfi) & (x_axis <= xfa)
            y_fit = (y_axis >= yfi) & (y_axis <= yfa)

            box = self.parent.diagonalBox(y_axis, x_axis,
                                              dd['fmin'], dd['fmax'], adspan)
            activeBox = data[box['idx0'], box['idx1']]

            simx_axis = np.linspace(np.min(dd['simax']['x']),
                                    np.max(dd['simax']['x']),
                                    dd['simdata'].shape[1])
            simy_axis = np.linspace(np.min(dd['simax']['y']),
                                    np.max(dd['simax']['y']),
                                    dd['simdata'].shape[0])
            simbox = self.parent.diagonalBox(simy_axis, simx_axis,
                                                 dd['fmin'], dd['fmax'], adspan)

            if dd['zmax'] < 0:
                if activeBox is not None:
                    scale = np.max(activeBox)
                elif x_fit.any() and y_fit.any():
                    scale = np.max(data[np.ix_(y_fit, x_fit)])

            interp = RegularGridInterpolator((hs.X, hs.Y), hs.Spectrum, bounds_error=False, fill_value=0.0)
    def _skewed_gaussian(self, x, amp, x0, sigma, shape):
        # Function to be used for Credible Interval estimations. 
        # x0 - the highest point position (not the middle of skewed gaussian)
        # if x0 is beyond the margins, we cannot extract shape. 
        # AlSi: I need to think about this on:
        # also, I would like to lock x0 to the bounds. 
        # othereiwse that will have a large error margin too
        if x0<np.min(x):
            x0 = np.min(x)
            shape = 0.0
        if x0>np.max(x):
            x0 = np.max(x)
            shape = 0.0
            
        d = shape/np.sqrt(1+shape**2)
        mz=d*np.sqrt(2/np.pi)
        sz = np.sqrt(1-mz**2)
        g1 = (4-np.pi)/2*mz**3/((1-mz**2)**(3/2))
        if shape==0:
            zm = mz-g1*sz/2
        else:
            zm = mz-g1*sz/2-np.sign(shape)/2*np.exp(-2*np.pi/np.abs(shape))
        z = (x - (x0-sigma*zm)) / (sigma)
        return amp * np.exp(-z**2 / 2) * (1 + sp.erf(shape * z / np.sqrt(2)))
    
    def calcCredible1D(self, param, rmsd, phys_bounds=None, credible_interval=0.95):
        amp0 = np.max(rmsd)
        x00 = param[np.argmax(rmsd)]
        sigma0 = np.sum(param*rmsd)/np.sum(rmsd)
        shape0 = 0.0
        if phys_bounds is not None:
            bounds = ([0.1*amp0,max(param.min()-sigma0*5, phys_bounds[0]), param[ 1]-param[0], -3],
                       [10*amp0, min(param.max()+sigma0*5, phys_bounds[1]), param[-1]-param[0],  3])
        else:
            bounds = ([0.1*amp0,param.min()-sigma0*5, param[ 1]-param[0], -3],
                      [10*amp0, param.max()+sigma0*5, param[-1]-param[0],  3])  
            
        prop, _ = curve_fit(self._skewed_gaussian, param, rmsd, p0=[amp0, x00, sigma0, shape0],
            bounds=bounds,
            maxfev=10000)
        amp, x0, sigma, shape = prop
        if phys_bounds is not None:
            x = np.linspace(max(x0-sigma*10, phys_bounds[0]),
                            min(x0+sigma*10, phys_bounds[1]), 
                            1000);
        else:
            x = np.linspace(x0-sigma*10,
                            x0+sigma*10, 
                            1000);
        yfit = self._skewed_gaussian(x, amp, x0, sigma, shape)
        sorted_y = np.sort(yfit)[::-1]
        #wi = np.argsort(yfit)[::-1]
        res = np.cumsum(sorted_y)
        res/=np.max(res) ## we assume we sample the complete space
        ii = np.argmin(np.abs(res-credible_interval))
        
        threshold_y = sorted_y[ii]
        midI = np.argmin(np.abs(x-x0))
        pt1 = np.argmin(np.abs(yfit[0:midI]-threshold_y))
        pt2 = np.argmin(np.abs(yfit[midI:-1]-threshold_y))+midI
        
        return x[pt1], x[pt2], x[midI], threshold_y
        
    def dumpXML(self):
        # otherwise ask the user what new file to open
        """
        On Mac macOS in the open file dialog the filter choice box is not shown by default. 
        Instead all given wildcards are applied at the same time: 
            To enforce the display of the filter choice set the corresponding 
            wx.SystemOptions before calling the file open dialog:
        """
        if hasattr(wx, 'OSX_FILEDIALOG_ALWAYS_SHOW_TYPES'):
            wx.SystemOptions.SetOption(wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES, 1)
        with wx.FileDialog(self, "Save Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultDir=self.file_browser.path, 
                           style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT) as fileDialog:
        
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
        
            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            try:
                fftmethod = {}
                for kk in self.tabulated_panel.FFTparams.keys():
                    prop = self.tabulated_panel.ffttree.GetProperty(kk)
                    fftmethod[kk] = prop.GetValue()
                    
                data_dict = {'Data': self.Data, 'BGData':self.BGData, 
                             'fftmethod': fftmethod, 
                             'Sys': self.tabulated_panel.Sys_param.parameters, 
                             'Exp': self.tabulated_panel.Exp_param.parameters, 
                             'Opt': self.tabulated_panel.Opt_param.parameters,
                             'cmap_name': self.current_cmap_name,
                             'cmap_stops': self.current_cmap_stops,
                             'view': self.get_view_state(),
                             }
                self.save_mixed_dict_xml(data_dict, pathname)
            except Exception as e:
                import traceback
                traceback.print_exc()
                wx.MessageBox(f"Could not save session to '{pathname}':\n{type(e).__name__}: {e}",
                              "Save Session Failed", wx.OK | wx.ICON_ERROR)
                
    def save_mixed_dict_xml(self, data_dict, filename):
        import xml.etree.ElementTree as ET
        from xml.dom import minidom
        import base64
        """
        Save a dictionary with mixed data types to an XML file.
        - NumPy arrays are encoded as base64 strings with metadata
        - Other data types (str, int, float, bool) are saved as text
        - Nested dictionaries are handled recursively
        """
        def create_element(parent, key, value):
            """Create XML element for a value"""
            element = ET.SubElement(parent, 'item')
            element.set('key', str(key))
    
            if isinstance(value, np.ndarray):
                # Handle numpy arrays
                array_element = ET.SubElement(element, 'numpy_array')
                array_element.set('dtype', str(value.dtype))
                array_element.set('shape', str(value.shape))
    
                # Encode array data as base64
                data_b64 = base64.b64encode(value.tobytes()).decode('utf-8')
                data_element = ET.SubElement(array_element, 'data')
                data_element.text = data_b64
    
            elif isinstance(value, dict):
                # Handle nested dictionaries
                dict_element = ET.SubElement(element, 'dict')
                for k, v in value.items():
                    create_element(dict_element, k, v)
    
            elif isinstance(value, (list, tuple)):
                # Handle lists and tuples
                list_element = ET.SubElement(element, 'list')
                list_element.set('type', type(value).__name__)
                for item in value:
                    create_element(list_element, 'item', item)
    
            else:
                # Handle primitive types
                primitive_element = ET.SubElement(element, 'primitive')
                primitive_element.set('type', type(value).__name__)
                primitive_element.text = str(value)
    
        # Create root element
        root = ET.Element('dictionary')
    
        # Process each item in the dictionary
        for key, value in data_dict.items():
            create_element(root, key, value)
    
        # Create pretty-printed XML
        rough_string = ET.tostring(root, encoding='unicode')
        reparsed = minidom.parseString(rough_string)
        pretty_xml = reparsed.toprettyxml(indent="  ")
    
        # Write to file
        with open(filename, 'w') as f:
            f.write(pretty_xml)
        
    def withDefaults(self, saved, defaults):
        """Saved session values laid over the current defaults.

        A session written before a parameter existed simply has no key for it,
        and SetFromParClean builds the property grid from exactly the dict it
        is handed -- so without this merge an older session silently drops
        every parameter added since (Opt.Grid was the first casualty) and the
        control vanishes from the tab.

        Enum entries are stored as [value, [choices]]. The saved selection is
        kept, but the choice list always comes from the defaults, so options
        added later show up on old sessions too; a selection that is no longer
        offered falls back to the default.
        """
        out = defaults.copy()
        for key, val in saved.items():
            dflt = defaults.get(key)
            isEnum = (isinstance(dflt, list) and len(dflt) == 2
                      and isinstance(dflt[1], list))
            if isEnum:
                sel = val[0] if (isinstance(val, list) and len(val) == 2) else val
                out[key] = [sel if sel in dflt[1] else dflt[0], list(dflt[1])]
            else:
                out[key] = val
        return out

    def loadXML(self, filename):
        curAutoSim_val = self.tb_autoSIM.GetValue()
        self.tb_autoSIM.SetValue(False)
        success = False
        try:
            result = self.load_mixed_dict_xml(filename)
            if 'Data' in result.keys():
                self.Data = result['Data'].copy()
            if 'BGData' in result.keys():
                self.BGData = result['BGData'].copy()
            if 'Sys' in result.keys():
                self.tabulated_panel.Sys_param.SetFromParClean(result['Sys'].copy())
            if 'Exp' in result.keys():
                self.tabulated_panel.Exp_param.SetFromParClean(result['Exp'].copy())
            if 'Opt' in result.keys():
                # merged, so a session saved before Opt.Grid existed still
                # gets the Grid dropdown (at its default) instead of losing it
                self.tabulated_panel.Opt_param.SetFromParClean(
                    self.withDefaults(result['Opt'], self.Opt.getDefaultDict()))
                self.Opt.setFromCtrl(self.tabulated_panel.Opt_param.parameters)
                self.tabulated_panel.updateGridInfo()
            if 'fftmethod' in result.keys():
                fftmethod = result['fftmethod']
                for kk in fftmethod.keys():
                    prop=self.tabulated_panel.ffttree.GetPropertyByName(kk)
                    prop.SetValue(fftmethod[kk])

            # ---- Restore colormap: saved stops, else by name, else jet -------
            name = str(result.get('cmap_name', 'jet') or 'jet')
            stops = None
            try:
                stops = [(float(p), str(c)) for p, c in result['cmap_stops']]
                if len(stops) < 2:
                    raise ValueError("fewer than two stops")
                self.current_cmap = stops_to_cmap(stops, name)
            except Exception as e:
                if 'cmap_stops' in result:
                    print(f"[pyHYSCORE] Saved colormap stops in '{filename}' unusable "
                          f"({type(e).__name__}: {e}); trying by name '{name}'.")
                try:
                    stops = colormap_to_stops(name, 9)
                    self.current_cmap = stops_to_cmap(stops, name)
                except Exception:
                    print(f"[pyHYSCORE] Unknown colormap '{name}'; using jet.")
                    name = "jet"
                    stops = colormap_to_stops("jet", 9)
                    self.current_cmap = stops_to_cmap(stops, "jet")
            self.current_cmap_name  = name
            self.current_cmap_stops = stops
            self.set_view_state(result.get('view', {}))

            self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters) # This is more to double check that input works
            self.Exp.setFromCtrl(self.tabulated_panel.Exp_param.parameters)                            
            self.Opt.setFromCtrl(self.tabulated_panel.Opt_param.parameters)  
            self.tabulated_panel.update_filelist()    
            self.matplotlib_panel.update_graph() 
            success = True
        except Exception as e:
            import traceback
            traceback.print_exc()
            wx.MessageBox(
                f"Loading '{filename}' failed:\n\n{type(e).__name__}: {e}\n\n"
                f"(full traceback printed to the console)",
                "Cannot load", wx.OK | wx.ICON_ERROR)
            
        
        if success:
            if self.Data and type(self.Data[0]['simdata'])!=type(None):
                dlg = wx.MessageDialog(self,
                        "The loaded session contains a simulation resut. \t Would you like to refresh it by running simulaiton again?",
                        "Actualize simulation", wx.YES_NO | wx.ICON_QUESTION)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result == wx.ID_YES:
                    self.runSim()
        self.tb_autoSIM.SetValue(curAutoSim_val)
            
    def load_mixed_dict_xml(self, filename):
        import xml.etree.ElementTree as ET
        import base64
        """
        Load a dictionary from XML file, reconstructing original structure and types.
        """
    
        def parse_element(element):
            """Parse XML element and reconstruct value"""
            # Get the key from the element
            key = element.get('key')
    
            # Check what type of element this is
            numpy_array = element.find('numpy_array')
            dict_element = element.find('dict')
            list_element = element.find('list')
            primitive_element = element.find('primitive')
    
            if numpy_array is not None:
                # Reconstruct numpy array
                dtype = numpy_array.get('dtype')
                shape = eval(numpy_array.get('shape'))  # Safe for simple shapes
                data_b64 = numpy_array.find('data').text
                data = base64.b64decode(data_b64)
                try:
                    array = np.array(np.frombuffer(data, dtype=dtype)).reshape(shape)
                except Exception as e:
                    print(f"[pyHYSCORE] Could not reconstruct array '{key}' "
                          f"(dtype={dtype}, shape={shape}) from '{filename}': "
                          f"{type(e).__name__}: {e}")
                    array = None
                return key, array
    
            elif dict_element is not None:
                # Reconstruct dictionary
                result = {}
                for item in dict_element.findall('item'):
                    if item.get('key') == "orisel":
                        pass 
                    k, v = parse_element(item)
                    if k is not None:
                        result[k] = v
                return key, result
    
            elif list_element is not None:
                # Reconstruct list/tuple
                items = []
                for item in list_element.findall('item'):
                    _, v = parse_element(item)
                    items.append(v)
                list_type = list_element.get('type')
                if list_type == 'tuple':
                    return key, tuple(items)
                else:
                    return key, items
    
            elif primitive_element is not None:
                # Reconstruct primitive type
                value_type = primitive_element.get('type')
                value_text = primitive_element.text
    
                if value_type == 'int':
                    return key, int(value_text)
                elif value_type == 'float':
                    return key, float(value_text)
                elif value_type == 'NoneType':
                    return key, None
                elif value_type == 'float64':
                    return key, float(value_text)
                elif value_type == 'bool':
                    return key, value_text.lower() == 'true'
                else:
                    return key, value_text
    
            return key, None
    
        # Parse XML file
        tree = ET.parse(filename)
        root = tree.getroot()
    
        # Reconstruct dictionary
        result = {}
        for item in root.findall('item'):
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
    frame.Show(True)

    app.MainLoop()