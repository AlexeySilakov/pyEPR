import os
import sys
import wx
import wx.lib.mixins.listctrl as listmix
import wx.grid as gridlib
import wx.propgrid as wxpg

import numpy as np
import re
import wx.lib.agw.flatnotebook as fnb

# Matplotlib imports
import matplotlib
#matplotlib.use('WXAgg')                 # Force the WXAgg backend
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg 
import matplotlib.tri as mtri
from matplotlib import cm
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

import brukerread as Mybr
import classPropGridPanel as MypgPanel
from SysPar import sysPar, expPar
from hyscore_sim import optHYSCORE, HYSCOREsim
# import wx.lib.agw.customtreectrl as CT
VERSION=0.1

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
    Add to plot the possibility of overlaying data and simdat.
    Add skyline projections
    Limit zero filling based on data size
    Preaty up the PropertyGrid CTRLS
    Try to recognize frequency from title and see if there is an IF ON flag to add 24.5GHz
    Fix MultipleChoice to only allow one entry
        enum_choices = ["Option 1", "Option 2", "Option 3", "Option 4"]
        self.pg.Append(wxpg.EnumProperty("Color", labels=enum_choices, value=0))
    Fix Lor-Gau
    Clicking in file in l
    oaded list opens that one but closes all others ... kind of like an accordeon deal

    
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
        
        # 1️⃣ Current path label
        hbox_path = wx.BoxSizer(wx.HORIZONTAL)
        
        self.lbl_path = wx.StaticText(self, label="")

        hbox_path.Add(self.lbl_path, 1,
                      wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 0)

        vbox.Add(hbox_path, 0, wx.EXPAND | wx.ALL, 0)

        # 2️⃣ Navigation buttons
        hbox_nav = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_change = wx.Button(self, label="Change folder…")
        self.btn_up = wx.Button(self, label="Up")
        hbox_nav.Add(self.btn_change, 0, wx.RIGHT, 0)
        hbox_nav.Add(self.btn_up, 0)
        vbox.Add(hbox_nav, 0, wx.ALIGN_LEFT | wx.ALL, 0)

        # 3️⃣ Filter combobox
        hbox_filter = wx.BoxSizer(wx.HORIZONTAL)
        lbl = wx.StaticText(self, label="Filter:")
        self.choice_filter = wx.ComboBox(
            self, choices=self.filter_choices, value=self.current_filter,
            style=wx.CB_READONLY)
        hbox_filter.Add(lbl, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        hbox_filter.Add(self.choice_filter, 1)
        vbox.Add(hbox_filter, 0, wx.EXPAND | wx.ALL, 5)

        # 4️⃣ One‑column list (must be created *before* loading a path!)
        self.list = wx.ListCtrl(
            self,
            style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_VRULES | wx.LC_HRULES)
        self.list.InsertColumn(0, "Name", width=400)
        self.list.InsertColumn(1, "Title", width=400)
        vbox.Add(self.list, 1,
                 wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 5)

        # Optional: folder icon
        img_list = wx.ImageList(16, 16)
        file_bitmap = wx.Bitmap('New document.ico', wx.BITMAP_TYPE_ICO)
        xml_bitmap = wx.Bitmap('Script.ico', wx.BITMAP_TYPE_ICO)
        folder_bitmap = wx.Bitmap('Folder.ico', wx.BITMAP_TYPE_ICO)
        dsc_bitmap = wx.Bitmap('List.ico', wx.BITMAP_TYPE_ICO)
        loaded_bitmap = wx.Bitmap('Green bookmark.ico', wx.BITMAP_TYPE_ICO)
        loadedBG_bitmap = wx.Bitmap('Equipment.ico', wx.BITMAP_TYPE_ICO)
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
        
        self.list.AssignImageList(img_list, wx.IMAGE_LIST_SMALL)

        # 5️⃣ Load button
        hbox_Load = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_load   = wx.Button(self, label="+ Data")
        self.btn_load.Bind(wx.EVT_BUTTON, self.on_load_clicked)
        hbox_Load.Add(self.btn_load, 0, wx.ALIGN_CENTER | wx.ALL, 2)
        
        self.btn_loadbg = wx.Button(self, label="+ BG")
        self.btn_loadbg.Bind(wx.EVT_BUTTON, self.on_loadBG_clicked)
        hbox_Load.Add(self.btn_loadbg, 0, wx.ALIGN_CENTER | wx.ALL, 2)
        
        self.btn_sessionLoad = wx.Button(self, label="+ Session")
        self.btn_sessionLoad.Bind(wx.EVT_BUTTON, self.on_load_session)
        hbox_Load.Add(self.btn_sessionLoad, 0, wx.ALIGN_CENTER | wx.ALL, 2)
        
        self.btn_sessionSave = wx.Button(self, label="Save Session")
        self.btn_sessionSave.Bind(wx.EVT_BUTTON, self.on_save_session)
        hbox_Load.Add(self.btn_sessionSave, 0, wx.ALIGN_CENTER | wx.ALL, 2)
        
        vbox.Add(hbox_Load, 0, wx.ALIGN_CENTER | wx.ALL, 0)
        
        mainbox.Add(vbox, 1, wx.EXPAND | wx.ALL, 0)
        
        self.toggle_btn = wx.Button(self, label="<", size=(15,-1))
        
        mainbox.Add(self.toggle_btn, 0, wx.EXPAND, 0)
        
        
        self.SetSizer(mainbox)

        # ------------------------------------------------------------------
        # 3️⃣ Sorting helpers / data storage
        # ------------------------------------------------------------------
        listmix.ColumnSorterMixin.__init__(self, 1)
        self.items = []          # (name, full_path, is_dir)

        # ------------------------------------------------------------------
        # 4️⃣ Event bindings (now that all widgets exist)
        # ------------------------------------------------------------------
        self.choice_filter.Bind(wx.EVT_COMBOBOX, self.on_filter_changed)
        
        
        
        self.btn_change.Bind(wx.EVT_BUTTON, self.on_change_folder)
        self.btn_up.Bind(wx.EVT_BUTTON, self.on_up_clicked)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_item_activated)
        
        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy)   # persistence

        self._load_path_from_file()          # this may call set_path()

        # ------------------------------------------------------------------
        # 6️⃣ Show the (possibly overwritten) current path
        # ------------------------------------------------------------------
        self.lbl_path.SetLabel(self.path)
        self.refresh_file_list()
    def on_save_session(self, event):
        self.parent.dumpXML()
    def on_load_session(self, event):
        idx = self.list.GetFirstSelected()
        if idx == -1:
            wx.MessageBox("Please select a file to load.", "No file selected",
                          wx.OK | wx.ICON_WARNING)
            return
        name, full, is_dir = self.items[idx]
        if is_dir:
            wx.MessageBox(f"{name} is a directory.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
        elif name.endswith('xml'):
            self.parent.loadXML(full)
  
        else:
            wx.MessageBox(f"The file {name} does not appear to be to have the right extension.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)  
            
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
        except Exception:
            return

        if stored and os.path.isdir(stored):
            # `set_path` will update lbl_path and list for us.
            self.set_path(stored)

    def _save_path_to_file(self):
        """Write the current path to the file."""
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.path + "\n")
        except Exception:
            pass

    # ----------------------------------------------------------------------
    # Core listing / filtering logic
    # ----------------------------------------------------------------------
    def get_entries(self):
        try:
            raw = os.listdir(self.path)
        except OSError:
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

    def refresh_file_list(self):
        self.list.DeleteAllItems()
        self.items.clear()

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

        self.list.SetColumnWidth(0, -1)
    def extract_TITL_from_dsc(self, file_path):
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith("TITL"):
                        parts = line.split("'")  # Split by single quotes
                        if len(parts) > 1:
                            return parts[1]  # Return the string between the first two quotes
            return ''  # No matching line found
        except FileNotFoundError:
            print(f"Error: File not found at {file_path}")
            return ''
    # ----------------------------------------------------------------------
    # Sorting mixin helpers
    # ----------------------------------------------------------------------
    def GetListCtrl(self):  return self.list
    def GetSortImages(self): return None, None
    def GetColumnSorterData(self): return self.items

    # ----------------------------------------------------------------------
    
    # ----------------------------------------------------------------------
    # Event handlers
    # ----------------------------------------------------------------------
    def on_filter_changed(self, event):
        self.current_filter = self.choice_filter.GetStringSelection()
        self.refresh_file_list()

    def on_load_clicked(self, event):
        idx = self.list.GetFirstSelected()
        if idx == -1:
            wx.MessageBox("Please select a file to load.", "No file selected",
                          wx.OK | wx.ICON_WARNING)
            return
        name, full, is_dir = self.items[idx]
        if is_dir:
            wx.MessageBox(f"{name} is a directory.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
        else:
            self.parent.loadData(full, BG=False)
        self.refresh_file_list()
    def on_loadBG_clicked(self, event):
        idx = self.list.GetFirstSelected()
        if idx == -1:
            wx.MessageBox("Please select a file to load.", "No file selected",
                          wx.OK | wx.ICON_WARNING)
            return
        name, full, is_dir = self.items[idx]
        if is_dir:
            wx.MessageBox(f"{name} is a directory.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
        else:
            self.parent.loadData(full, BG=True)
        self.refresh_file_list()
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
        name, full, is_dir = self.items[idx]
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
            self.lbl_path.SetLabel(self.path)
            self.refresh_file_list()
            self._save_path_to_file()
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
        self.figure = Figure(figsize=(2, 2), dpi=100)
        self.axes = []
        
        #self.axes.append(self.figure.add_subplot(111))
        #self.axes.set_title("Placeholder Figure")
        #self.axes.plot([0, 1, 2], [0, 1, 4], marker='o')

        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        


        self.chk_FFT = wx.CheckBox(self, label="FFT(y)")
        self.chk_FFT.SetValue(self.showFFT)
        
        self.chk_BG = wx.CheckBox(self, label="Background")
        self.chk_BG.SetValue(self.showBG)
        
        self.chk_Skyline = wx.CheckBox(self, label="Skyline")
        self.chk_Skyline.SetValue(self.showSkyline)
        
        self.chk_OriSel = wx.CheckBox(self, label="Orientat maps")
        self.chk_OriSel.SetValue(self.showOriSel)
        
        self.chk_DiagProj = wx.CheckBox(self, label="Diag.Projection")
        self.chk_DiagProj.SetValue(self.showDiagonalProj)
        
        self.chk_FFT.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_BG.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_Skyline.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_OriSel.Bind(wx.EVT_CHECKBOX, self.on_check)
        self.chk_DiagProj.Bind(wx.EVT_CHECKBOX, self.on_check)
        
        sizerH = wx.BoxSizer(wx.HORIZONTAL)
        sizerH.Add(self.chk_FFT, 0, wx.EXPAND)
        sizerH.Add(self.chk_BG, 0, wx.EXPAND)
        sizerH.Add(self.chk_Skyline, 0, wx.EXPAND)
        sizerH.Add(self.chk_OriSel, 0, wx.EXPAND)
        sizerH.Add(self.chk_DiagProj, 0, wx.EXPAND)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerH, 0, wx.EXPAND)
        
        sizer.Add(self.canvas, 1, wx.EXPAND)
        
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version.
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        # update the axes menu on the toolbar
        self.toolbar.update()
        
        self.SetSizer(sizer)
    def on_check(self, event):
        (self.showFFT, self.showBG, self.showSkyline, self.showOriSel,self.showDiagonalProj) = (
        self.chk_FFT.GetValue(),
        self.chk_BG.GetValue(),
        self.chk_Skyline.GetValue(),
        self.chk_OriSel.GetValue(),
        self.chk_DiagProj.GetValue(),
        )
        self.update_graph()
    def update_graph(self):
        shw = []
        for dd in self.parent.Data:
            shw.append(dd['show'])
        nData = shw.count(True)
        # self.showBG, self.showSkyline, let's not .. for now
        nrows = [self.showSim, self.showOriSel,self.showDiagonalProj].count(True)+1
        print(nrows)
        makenew = False
        if len(self.axes)!=nData*nrows:
            for ax in self.axes:
                if ax is not None:
                    ax.clear()
                    ax.remove()
            self.axes = [None]*(nData*nrows)
            self.figure.clf()
            makenew = True
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
                    raise AttributeError("No FFT has been performed yet. Nothing to display")    
                axEx = [np.min(dd['fftax']['x']), np.max(dd['fftax']['x']), 
                        np.min(dd['fftax']['y']), np.max(dd['fftax']['y'])]
                
                ma = np.max(np.max(data))*dd['zmax']
                mi = ma*dd['zmin'] 
                xmi,xma, ymi, yma = (dd['fmin'], dd['fmax'], dd['fmin'], dd['fmax'] )
                axlabel='Time, $\mu$s'
            else:
                data = np.real(dd['data'][:,:,0])
                axEx = [np.min(dd['ax']['x']), np.max(dd['ax']['x']), 
                        np.min(dd['ax']['y']), np.max(dd['ax']['y'])]
                if not dd['isfft']:
                    ma = np.max(np.max(data)); mi = np.min(np.min(data))
                    xmi = np.min(dd['ax']['x'])
                    xma = np.max(dd['ax']['x'])
                    ymi = np.min(dd['ax']['y'])
                    yma = np.max(dd['ax']['y'])  
                else:
                    ma = np.max(np.max(data))*dd['zmax']
                    mi = ma*dd['zmin'] 
                    xmi,xma, ymi, yma = (dd['fmin'], dd['fmax'], dd['fmin'], dd['fmax'] )
                axlabel='Frequency, MHz'
            if makenew:                                             
                self.axes[axcnt]=self.figure.add_subplot(nrows, nData, axcnt+1)
                
                im = self.axes[axcnt].imshow(
                    data, extent=axEx, origin="lower", cmap="jet", interpolation="nearest")
                im.set_clim(vmin=mi, vmax=ma)
            else:
                for ch in self.axes[axcnt].get_children():
                    if isinstance(ch, matplotlib.image.AxesImage):
                        ch.set_data(data)
                        ch.set_extent(axEx)
                        ch.set_clim(vmin=mi, vmax=ma)
                        break
            self.axes[axcnt].set_xlim(xmi, xma)
            self.axes[axcnt].set_ylim(ymi, yma)
            self.axes[axcnt].set_title(f'B$_0$={dd['field']} mT', fontsize = 'small')
            if showLabels:
                # self.axes[axcnt].set_xlabel(axlabel, fontsize = 'small')
                self.axes[axcnt].set_ylabel(axlabel, fontsize = 'small')
            rowcnt+=1
            ########## ----- simulation ---------------------------------
            if self.showSim:
                simax = axcnt+rowcnt*nData
                realsim = ''
                if type(dd['simdata'])!=type(None):
                    simdata = dd['simdata']
                    sima = np.max(np.max(simdata))
                    simi = ma*dd['zmin']/dd['zmax']
                    simaxEx = [np.min(dd['simax']['x']), np.max(dd['simax']['x']), 
                            np.min(dd['simax']['y']), np.max(dd['simax']['y'])]
                    
                else:
                    simdata = np.zeros_like(data)
                    simi,sima = (0.0, 1.0)
                    simaxEx= axEx
                    if type(dd['simax'])!=type(None):
                        realsim='simulation failed'
                    else:
                        realsim='no simulation'
                    
                if makenew:                                             
                    self.axes[simax]=self.figure.add_subplot(nrows, nData, simax+1)
                    im = self.axes[simax].imshow(
                        simdata, extent=simaxEx, origin="lower", cmap="jet", interpolation="nearest")
                    im.set_clim(vmin=simi, vmax=sima)
                else:
                    for ch in self.axes[simax].get_children():
                        if type(ch)==matplotlib.text.Annotation:
                            ch.remove()
                        if type(ch)== matplotlib.image.AxesImage:
                            ch.set_data(simdata)
                            ch.set_extent(simaxEx)
                            ch.set_clim(vmin=simi, vmax=sima)
                            
                        
                # axes lim set to data specs
                self.axes[simax].set_xlim(xmi, xma)
                self.axes[simax].set_ylim(ymi, yma)
                if len(realsim)>0:
                    self.axes[simax].annotate("no simulation", (.0, .0), xycoords='axes points', color='w')
                rowcnt+=1
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
                    norm = plt.Normalize(vmin=0, vmax=np.max(ak))
                    tri = mtri.Triangulation(x, y)    
                    triangle_ak = ak[tri.triangles].max(axis=1)          # one value per triangle
                    facecolors = cm.jet(norm(triangle_ak))
                    
                    surf = self.axes[oriax].plot_trisurf(
                        x, y, z,
                        triangles=tri.triangles,
                        linewidth=0, antialiased=True,
                        shade=False,  # we supply colour via facecolors
                        )
                        #
                    self.axes[oriax].set_aspect('equal') 
                    self.axes[oriax].set_xlabel('X'); self.axes[oriax].set_ylabel('Y'); self.axes[oriax].set_zlabel('Z')
                    surf.set_facecolor(facecolors)
                    surf.set_edgecolor(cm.jet(norm(ak)))
                rowcnt+=1
            if self.showDiagonalProj:
                diagprj = axcnt+rowcnt*nData
                if makenew:  
                    self.axes[diagprj]=self.figure.add_subplot(nrows, nData, diagprj+1)   
                else:
                    self.axes[diagprj].cla()    
                npts = data.shape[0]
                dx = (axEx[1]-axEx[0])/(npts-1)
                dy = (axEx[3]-axEx[2])/(npts-1)
                idxmi = int(np.floor((xmi-axEx[0])/dx))
                idxma = int(np.floor((xma-axEx[0])/dx))
                idymi = int(np.floor((ymi-axEx[2])/dy))
                idyma = int(np.floor((yma-axEx[2])/dy))
                ma = np.max(np.max(data))*dd['zmax']

                rotdata = self.rotate45(data[idxmi:idxma, idymi:idyma]/ma)
                skyprj = np.max(rotdata, axis=0)
                xaX = np.linspace(idxmi*dx+axEx[0], idxma*dx+axEx[0], rotdata.shape[0]) - (xma+xmi)/2
                self.axes[diagprj].plot(xaX, skyprj, color='b')

                if type(dd['simdata'])!=type(None):
                    simdata = dd['simdata']
                    sima = np.max(np.max(simdata))
                    simi = ma*dd['zmin']/dd['zmax']
                    simaxEx = [np.min(dd['simax']['x']), np.max(dd['simax']['x']), 
                            np.min(dd['simax']['y']), np.max(dd['simax']['y'])]
                    npts = simdata.shape[0]
                    dx = (simaxEx[1]-simaxEx[0])/(npts-1)
                    dy = (simaxEx[3]-simaxEx[2])/(npts-1)
                    idxmi = int(np.floor((xmi-simaxEx[0])/dx))
                    idxma = int(np.floor((xma-simaxEx[0])/dx))
                    idymi = int(np.floor((ymi-simaxEx[2])/dy))
                    idyma = int(np.floor((yma-simaxEx[2])/dy))
                    ma = np.max(np.max(simdata))

                    rotdata = self.rotate45(simdata[idxmi:idxma, idymi:idyma]/ma)
                    skyprj = np.max(rotdata, axis=0)
                    xaX = np.linspace(idxmi*dx+simaxEx[0], idxma*dx+simaxEx[0], rotdata.shape[0]) - (xma+xmi)/2
                    self.axes[diagprj].plot(xaX, skyprj, color='r')


                rowcnt+=1
                
            axcnt+=1
        self.figure.tight_layout()
        self.canvas.draw()
    def rotate45(self, A):
        N = A.shape[0]
        out_size = 2 * N - 1
        out = np.zeros((out_size, out_size), dtype=A.dtype)
        ii, jj = np.indices((N, N))
        x = ii + jj
        y = ii - jj + (N - 1)
        out[x, y] = A
        return out
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kwargs):
        super().__init__(parent,size=(200, 400), *args, **kwargs)
        self.parent = mainWindow
        
        self.autoFFT = False
        self.autoSim = False
        self.FFTparams = {
            "polynomial":["Polynom", 2, [], 'int'],
            "apodization":["Apodization", ['Hamming'], ['Hamming', 'Lor-Gau', 'Gaussian'], 'list'],
            "awidth":["A.Width",  1, [], 'float'],
            "aalpha":["Lor.Width",0.6, [], 'float' ],
            "ashift":["A.Shift",  0.0, [], 'float'],
            "zfill": ["Zero fill" , 2, [], 'int'],
            "imag":  ["Use 'imaginary'",  False, [], 'bool'],
            }
        
        
        # self.nb = wx.Notebook(self)
        self.nb = fnb.FlatNotebook(self, agwStyle=fnb.FNB_NO_X_BUTTON|fnb.FNB_NO_NAV_BUTTONS|
                                  fnb.FNB_SMART_TABS)
        #self.nb.SetActiveTabColour(wx.Colour("black"))
        self.nb.SetActiveTabTextColour(wx.Colour("black"))
        self.nb.SetBackgroundColour(wx.Colour("white"))
        # bold_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        # self.nb.SetFont(bold_font)

        # Main horizontal sizer – grid on the left, params on the right
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.add_datapanel()
        self.add_syspanel()
        self.add_exppanel()
        self.add_optpanel()
        
        # font = wx.Font(20, wx.FONTFAMILY_TELETYPE,
        #        wx.FONTSTYLE_NORMAL,
        #        wx.FONTWEIGHT_BOLD)
        # self.nb.SetTabAreaColour(wx.Colour("white"))
        
        # self.nb.SetFont(font)

        
        hbox_FFT = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_doFFT   = wx.Button(self, label="Do FFT")
        self.chk_autoFFT = wx.CheckBox(self, label="auto FFT")
        self.chk_autoFFT.SetValue(True)
        hbox_FFT.Add(self.btn_doFFT, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        hbox_FFT.Add(self.chk_autoFFT, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        sizer.Add(hbox_FFT, 0, wx.EXPAND | wx.ALL, 5)
        
        hbox_SIM = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_doSIM   = wx.Button(self, label="Do SIM")
        self.chk_autoSIM = wx.CheckBox(self, label="auto SIM")
        self.chk_autoSIM.SetValue(True)
        hbox_SIM.Add(self.btn_doSIM, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        hbox_SIM.Add(self.chk_autoSIM, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        sizer.Add(hbox_SIM, 0, wx.EXPAND | wx.ALL, 5)
        self.btn_doSIM.Bind(wx.EVT_BUTTON, self.on_update_sim)
        
        self.SetSizer(sizer)

        # Store for later
        self.param_definitions = None
        self.param_controls   = {}
       
        
        self.btn_doFFT.Bind(wx.EVT_BUTTON, self.on_doFFT)
        
        self.ffttree.Bind(wxpg.EVT_PG_CHANGED, self.on_ffttree)
        #self.filetree.Bind(wxpg.EVT_PG_DOUBLE_CLICK, self.on_filetree) #EVT_G_SELECTED EVT_PG_RIGHT_CLICK
        #self.tree.Bind(CT.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit)
        #self.tree.Bind(CT.EVT_TREE_END_LABEL_EDIT, self.on_end_edit)        
        self.filetree.Bind(wxpg.EVT_PG_CHANGED, self.on_filetree)
        #self.filetree.Bind(wxpg.EVT_PG_DOUBLE_CLICK, self.on_filetree) #EVT_G_SELECTED EVT_PG_RIGHT_CLICK
        #self.tree.Bind(CT.EVT_TREE_BEGIN_LABEL_EDIT, self.on_begin_edit)
        #self.tree.Bind(CT.EVT_TREE_END_LABEL_EDIT, self.on_end_edit)
    def add_optpanel(self):
        # ------------- Opt
        
        self.Opt_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.Opt_panel, "Opt") 
        Optsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Opt_param = MypgPanel.PropGridPanel(self.Opt_panel, Prop_Dict = {}, onChangeFunc=self.on_Opt)
        self.Opt_param.SetFromParClean(self.parent.Opt.getDefaultDict())
        Optsizer.Add(self.Opt_param, 1, wx.EXPAND)
        self.Opt_panel.SetSizer(Optsizer)
        
    def add_exppanel(self):
        # ------------- Exp
        self.Exp_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.Exp_panel, "Exp") 
        
        Expsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Exp_param = MypgPanel.PropGridPanel(self.Exp_panel, Prop_Dict = {}, onChangeFunc=self.on_Exp)
        self.Exp_param.SetFromParClean(self.parent.Exp.getDefaultDict())
        Expsizer.Add(self.Exp_param, 1, wx.EXPAND)
        self.Exp_panel.SetSizer(Expsizer)
         
    def add_syspanel(self):
        # ------------- Sys 
        self.Sys_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.Sys_panel, "Sys") 
        
        Syssizer = wx.BoxSizer(wx.VERTICAL)
        
        btnszr = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_add = wx.Button(self.Sys_panel, label="Add Nuc")
        self.btn_add.Bind(wx.EVT_BUTTON, self.on_add_Nuc)
        btnszr.Add(self.btn_add, 0, wx.EXPAND)
        self.btn_delete = wx.Button(self.Sys_panel, label="Delete Nuc")
        self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete_nuc)
        btnszr.Add(self.btn_delete, 0, wx.EXPAND)
       
        Syssizer.Add(btnszr, 0, wx.EXPAND)
        
        self.Sys_param = MypgPanel.PropGridPanel(self.Sys_panel, Prop_Dict = {}, onChangeFunc=self.on_Sys)
        self.Sys_param.SetFromParClean(self.parent.Sys.getDefaultDict())
        
        Syssizer.Add(self.Sys_param, 1, wx.EXPAND)
      
        self.Sys_panel.SetSizer(Syssizer)
               
    def add_datapanel(self):
        # ----------Data panel
        self.data_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.data_panel, "Data")
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_delete = wx.Button(self.data_panel, label="Delete Data")
        self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete)
        sizerBtns.Add(self.btn_delete, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        
        self.btn_dsc = wx.Button(self.data_panel, label="Show DSC")
        self.btn_dsc.Bind(wx.EVT_BUTTON, self.on_dsc)
        sizerBtns.Add(self.btn_dsc, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        sizer.Add(sizerBtns, 0, wx.ALIGN_LEFT | wx.ALL, 5)
        
        # file tree
        self.filetree = wxpg.PropertyGrid(
            self.data_panel,
            style=wxpg.PG_DEFAULT_STYLE | 
                  wxpg.PG_SPLITTER_AUTO_CENTER |
                  wxpg.PG_HIDE_MARGIN |
                  wxpg.PG_NO_INTERNAL_BORDER
        )
        self.filetree.SetCellBackgroundColour(wx.Colour(255, 255, 255))
        self.filetree.SetCellTextColour(wx.Colour(0, 0, 0))
        self.filetree.SetMarginColour(wx.Colour(64, 120, 215))
        #self.filetree.SetSelectionBackgroundColour(wx.Colour(200, 255, 255))
        
        self.filetree.SetCaptionBackgroundColour(wx.Colour(64, 120, 215))
        self.filetree.SetCaptionTextColour(wx.Colour(255, 255, 255))
        cat1 = self.filetree.Append(wxpg.PropertyCategory("Nothing loaded yet"))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("MW Freq", value=9.43))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("B0", value=350.0))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("tau", value=120.0))
        
        sizer.Add(self.filetree, 1, wx.EXPAND | wx.ALL, 0)
        
        lbl = wx.StaticText(self.data_panel, label="FFT Parameters:")
        sizer.Add(lbl, 0, wx.EXPAND, 0)
        # FFT tree
        self.ffttree = wxpg.PropertyGrid(
            self.data_panel,
            style=wxpg.PG_DEFAULT_STYLE | 
                  wxpg.PG_SPLITTER_AUTO_CENTER |wxpg.PG_HIDE_MARGIN |
                  wxpg.PG_NO_INTERNAL_BORDER
        )
        self.ffttree.SetCellBackgroundColour(wx.Colour(255, 255, 255))
        self.ffttree.SetCellTextColour(wx.Colour(0, 0, 0))
        self.ffttree.SetMarginColour(wx.Colour(215, 120, 64))
        #self.filetree.SetSelectionBackgroundColour(wx.Colour(200, 255, 255))
        
        self.ffttree.SetCaptionBackgroundColour(wx.Colour(215, 120, 64))
        self.ffttree.SetCaptionTextColour(wx.Colour(255, 255, 255))
        #cat2 = self.ffttree.Append(wxpg.PropertyCategory("FFT"))
        

        for kk in self.FFTparams.keys():
            if self.FFTparams[kk][3]=='int':
                self.ffttree.Append(wxpg.IntProperty(self.FFTparams[kk][0], kk, value=self.FFTparams[kk][1]))
            elif self.FFTparams[kk][3]=='float':
                self.ffttree.Append(wxpg.FloatProperty(self.FFTparams[kk][0], kk, value=self.FFTparams[kk][1])) 
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
        isotopes = list(self.parent.Sys.isotopes().keys())
        for dd in isotopes:
            items_to_show.append(dd)
        
        dlg = wx.SingleChoiceDialog(
            self,
            "Select the parameters you want to show:",
            "Choose parameters",
            items_to_show,
            style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
        )
    
        # Show the dialog modally â it will block until the user presses OK/Cancel
        if dlg.ShowModal() == wx.ID_OK:
            # dlg.GetSelections() returns indices of the checked items
            selected_indice = dlg.GetSelection()
        else:
            selected_indice = -1

        dlg.Destroy()
        if selected_indice>=0:
            
            Params = self.Sys_param.parameters
            nNucs = 0
            for key in Params:
                if 'nuc' in key:
                    nNucs+=1
            Params[f'nuc({nNucs+1})']={'Nucs':isotopes[selected_indice],
                                     'nNucs':1,
                                     'A':np.array([0, 0, 0]),
                                     'Apa':np.array([0, 0, 0])}
            _,I,_=self.parent.Sys.isotopes(isotopes[selected_indice])
            if I[0]>0.5:
                Params[f'nuc({nNucs+1})']['Q']=np.array([0, 0, 0])
                Params[f'nuc({nNucs+1})']['Qpa']=np.array([0, 0, 0])
            
            self.parent.Sys.setFromCtrl(Params)
            self.parent.Sys_param.SetFromParClean(Params)
            
            # for de in selected_indices:
            #     self.parent.Data.pop(de)
            # self.parent.update_everything()
            
    def on_delete_nuc(self, event):
        pass
    def on_update_sim(self, event):
        self.parent.runSim()
    def on_Sys(self, mainname, name, val):
        self.parent.Sys.setFromCtrl(self.Sys_param.parameters) # This is more to double check that input works
        if self.chk_autoSIM.GetValue():
            self.parent.runSim()
        
    def on_Exp(self, mainname, name, val):
        self.parent.Exp.setFromCtrl(self.Exp_param.parameters)                            
    def on_Opt(self, mainname, name, val):
        self.parent.Opt.setFromCtrl(self.Opt_param.parameters)                            
    def on_dsc(self, event):
        items_to_show = []
        selected_items = []
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
            print(self.parent.Data[selected_indices]['dsc'])
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
        selected_items = []
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
    def on_filetree(self, event):
        (key, strval) = event.GetPropertyName().split()
        val = int(strval)
        self.parent.Data[val][key]=event.GetValue()
        if self.chk_autoFFT.GetValue():
            self.parent.update_FFT()
            
    def on_ffttree(self, event):
        fftmethod = {}
        for kk in self.FFTparams.keys():
            prop = self.ffttree.GetProperty(kk)
            fftmethod[kk] = prop.GetValue()

        for ii, dd in enumerate(self.parent.Data):
            self.parent.Data[ii]['fftmethod'] = fftmethod
            self.parent.Data[ii]['fftactual'] = False # make sure to remember that things got changed.

        if self.chk_autoFFT.GetValue():
            self.parent.update_FFT(fftmethod)
    # ------------------------------------------------------------------
    # Simulations tab (empty for now)
    # ------------------------------------------------------------------
    def _create_simulation_tab(self):
        self.sim_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.sim_panel, "Simulations")

        placeholder = wx.StaticText(
            self.sim_panel,
            label="Simulation UI goes here.\nAdd buttons, charts, etc.",
        )
        placeholder.Wrap(250)

        sim_sizer = wx.BoxSizer(wx.VERTICAL)
        sim_sizer.AddStretchSpacer()
        sim_sizer.Add(placeholder, 0, wx.ALIGN_CENTER)
        sim_sizer.AddStretchSpacer()
        self.sim_panel.SetSizer(sim_sizer)
    def update_filelist(self):
        print('tbd')
        self.filetree.Clear()
        cat = []
        bgchoices = ['none']
        for bg in self.parent.BGData:
            bgchoices.append(bg['title'])
            
        for ii, dd in enumerate(self.parent.Data):
            cat.append(self.filetree.Append(wxpg.PropertyCategory(dd['title'])))
            cat[ii].SetBackgroundColour(wx.Colour(255, 255, 255))
            cat[ii].SetTextColour(wx.Colour(0, 0, 180))
            self.filetree.AppendIn(cat[ii], wxpg.BoolProperty("Show", f"show {ii}",value=dd['show']))
            self.filetree.AppendIn(cat[ii], wxpg.StringProperty("f.name", f"fname {ii}",value=dd['fname']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("MW Freq [GHz]", f"freq {ii}", value=dd['freq']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("B0 [mT]", f"field {ii}", value=dd['field']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("tau [ns]", f"tau {ii}",value=dd['tau']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("z.max", f"zmax {ii}",value=dd['zmax']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("z.min", f"zmin {ii}",value=dd['zmin']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("f.max", f"fmax {ii}",value=dd['fmax']))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("f.min", f"fmin {ii}",value=dd['fmin']))
            if dd['bgdata'] in bgchoices:
                val = [dd['bgdata']]
            else:
                ### this is to make sure everything lines up
                val = ['none']
                self.parent.Data[ii]['bgdata']=val
                
            self.filetree.AppendIn(cat[ii], wxpg.MultiChoiceProperty("BG Data", f"bgdata {ii}",
                                                                     choices=bgchoices,
                                                                     value=val))
            self.filetree.AppendIn(cat[ii], wxpg.FloatProperty("BG Scale", f"bgscale {ii}",value=dd['bgscale']))
            
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
        # 1ï¸â£ Root splitter:  left | right
        self.splitter_main = wx.SplitterWindow(self, style=wx.SP_3D)
        self.splitter_main.SetMinimumPaneSize(20)     # keep left pane from getting too small

        # ---- Left pane ----------------------------------------------------
        self.file_browser = FileBrowserPanel(self.splitter_main, self)

        # ---- Right pane (will contain a *second* splitter) ---------------
        right_panel = wx.Panel(self.splitter_main)
        right_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # 2ï¸â£ Second splitter (horizontal) inside the right panel
        self.splitter_right = wx.SplitterWindow(right_panel, style=wx.SP_3D)

        # Panels for the second splitter
        self.matplotlib_panel = MatplotlibPanel(self.splitter_right, self)
        self.tabulated_panel   = TabulatedPanel(self.splitter_right, self)

        # Split the right panel horizontally
        self.splitter_right.SplitVertically(
            self.matplotlib_panel, self.tabulated_panel,
            sashPosition=400   # initial height of the top pane
        )
        self.splitter_right.SetSashGravity(1)  # proportion of space for the top pane
        
        # Put the second splitter into the right panel's sizer
        right_sizer.Add(self.splitter_right, 1, wx.EXPAND)
        right_panel.SetSizer(right_sizer)

        # ---- Attach the two topâlevel panes --------------------------------
        self.splitter_main.SplitVertically(
            self.file_browser, right_panel,
            sashPosition=200      # initial width of left pane
        )
        self.splitter_main.SetSashGravity(0)  # proportion of space for left pane

        # ----------- Enable collapsing the left panel ----------------------
        # Store the last sash position so we can restore it when expanding
        self._prev_left_sash = 200
        self.file_browser.toggle_btn.Bind(wx.EVT_BUTTON, self.on_toggle_left)
        # self.splitter_right.SetSashSize(10) 
        # self.splitter_main.SetSashSize(10) 
        self.Centre()
        # self.Fit()  
    # -----------------------------------------------------------------------
    def on_toggle_left(self, event):
        """Called when the leftâpanel collapse/expand button is pressed."""
        minsize = self.splitter_main.GetMinimumPaneSize()
        if self.splitter_main.GetSashPosition() > minsize:
            # Collapse: store current position and move sash to zero
            self._prev_left_sash = self.splitter_main.GetSashPosition()
            self.splitter_main.SetSashPosition(minsize)
            self.file_browser.toggle_btn.SetLabel(">")
        else:
            # Expand: restore previous position
            self.splitter_main.SetSashPosition(self._prev_left_sash)
            self.file_browser.toggle_btn.SetLabel("<")

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
                    'isfft':False,
                    'fftdata':None, 'fftax':None, 'fftmethod':None, 'fftactual':False,
                    'simdata':None, 'simax':None, 'simmethod':None, 'simactual':False,
                    'orisel':None,
                    
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

        for ii, dd in enumerate(self.Data):
            if not fftmethod:
                fftmethod=self.Data[ii]['fftmethod'].copy()
            
            if len(dd['data'].shape)>2: ### need to fix at the source. Just for now let's deal with it.
                Data = np.array(dd['data'][:,:,0])
            else:
                Data = np.array(dd['data'])
                
            if dd['isfft']:
                return
            
            # dd['bgdata'] is comming from a MultiChoice property, so it's value is a list
            if dd['bgdata'][0]!='none':
                BGscale= dd['bgscale']
                
                gotBG = False
                for bb in self.BGData:
                    if dd['bgdata'][0]==bb['title']:
                        if len(bb['data'].shape)>2: ### need to fix at the source. Just for now let's deal with it.
                            BG = np.array(bb['data'][:,:,0])
                        else:
                            BG = np.array(bb['data'])
                        gotBG = True
                        break
                if gotBG:
                    if Data.shape !=BG.shape:
                        raise ValueError("😭 Background and Data have different dimensions.")
                    Data -=BGscale*BG
                
            if not fftmethod['imag']:
                Data=np.real(Data)
                
            X=np.array(dd['ax']['x'])
            Y=np.array(dd['ax']['y'])
                                       
            for cc in range(Data.shape[0]):
                pp = np.polyfit(Y, Data[cc,:], fftmethod['polynomial'])
                Data[cc,:]-=np.polyval(pp, Y)
            for cc in range(Data.shape[1]):
                pp = np.polyfit(X, Data[:, cc], fftmethod['polynomial'])
                Data[:, cc]-=np.polyval(pp, X)
            
            if fftmethod['apodization'][0]!='none':
                tX, tY = np.meshgrid(X,Y)
                awidth= fftmethod['awidth']
                ax = np.max(X)*awidth
                ay = np.max(Y)*awidth
                if fftmethod['apodization'][0]=='Hamming':
                    app = (27/50+23/50*np.cos(np.pi*tX/ax) )*(27/50+23/50*np.cos(np.pi*tY/ay) )
                elif fftmethod['apodization']=='Gaussian':
                    sig = (ax+ay)/2/2.354820045 # 2.354820045 = 2sqrt(2ln2)
                    app = np.exp(-0.5/sig**2 *( tX**2 + tY**2) )
                elif fftmethod['apodization']=='Lor-Gau':
                    alpha = fftmethod['aalpha']
                    shift = fftmethod['ashift']
                    sig = (ax+ay)/2/2.354820045 # 2.354820045 = 2sqrt(2ln2)
                    app = ( np.exp(-0.5/sig**2* ( (tX-shift)**2 + (tY-shift)**2) )*
                            np.exp(tX/ax/alpha + tY/ay/alpha) 
                           )
                Data*=app
                
            if fftmethod['zfill']>0:
                ndim = np.ceil(np.log2(Data.shape[0]))+fftmethod['zfill']
                
                padn = int(2**ndim-Data.shape[0])
                Data=np.pad(Data, (0,padn))
            
            fftData = np.abs(np.fft.fftshift(np.fft.fft2(Data)))
            nyqX = 1/2/(X[2]-X[1])
            nyqY = 1/2/(Y[2]-Y[1])
            
            fftX = np.linspace(-nyqX, nyqX, fftData.shape[0])
            fftY = np.linspace(-nyqY, nyqY, fftData.shape[1])
            self.Data[ii]['fftdata']=fftData
            self.Data[ii]['fftax']={'x':fftX, 'y':fftY, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz'}
            self.Data[ii]['fftactual'] = True
        self.matplotlib_panel.update_graph()
    def runSim(self):

        self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters)
        
        hs = HYSCOREsim(Sys=self.Sys)
        hs.preCompute() ## get housekeeping stuff out of the way to speed up computations a bit
        
        for ii,dd in enumerate(self.Data):
            methoddic = {'Sys': self.tabulated_panel.Sys_param.parameters.copy(), 
                         'Exp': self.tabulated_panel.Exp_param.parameters.copy(),
                         'Opt': self.tabulated_panel.Opt_param.parameters.copy(),}
            if methoddic['Exp']['tau']<0:
                methoddic['Exp']['tau'] = self.Data[ii]['tau']
            if methoddic['Exp']['Field']<0:
                methoddic['Exp']['Field'] = self.Data[ii]['field']
            if methoddic['Exp']['mwFreq']<0:
                methoddic['Exp']['mwFreq'] = self.Data[ii]['freq']
            if methoddic['Exp']['MaxFreq']<0:
                methoddic['Exp']['MaxFreq'] = self.Data[ii]['fmax']
            if methoddic['Exp']['nPoints']<0:
                methoddic['Exp']['nPoints'] = len(self.Data[ii]['fftax']['x'])
            
            if type(self.Data[ii]['orisel'])!=type(None):
                methoddic['Opt']['OriSelInp'] = self.Data[ii]['orisel']
            
            self.Data[ii]['simmethod']=methoddic
            self.Data[ii]['simactual']=False
        
            self.Opt.setFromCtrl(self.Data[ii]['simmethod']['Opt'])
            self.Exp.setFromCtrl(self.Data[ii]['simmethod']['Exp'])
            print(self.Exp.Field)
            hs.Opt = self.Opt
            hs.Exp = self.Exp
            
            try:
                hs.reRun()
                
                self.Data[ii]['simdata']=hs.Spectrum
                self.Data[ii]['simax']={'x':hs.X, 'y':hs.Y, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz', 
                                        'orisel': np.array([hs.phi, hs.theta, hs.ak])}
                self.Data[ii]['simactual']=True
            except Exception:
                self.Data[ii]['simdata']=None
                self.Data[ii]['simax']={'x':hs.X, 'y':hs.Y, 'xlabel':'Frequency, MHz', 'ylabel':'Frequency, MHz', 
                                        'orisel': np.array([hs.phi, hs.theta, hs.ak])}
                self.Data[ii]['simactual']=True
                
            
        self.matplotlib_panel.update_graph()  
            
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
                             }
                self.save_mixed_dict_xml(data_dict, pathname)
            except IOError:
                wx.LogError("Cannot open file '%s'." % newfile)
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
    def loadXML(self, filename):
        
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
                self.tabulated_panel.Opt_param.SetFromParClean(result['Opt'].copy())
            if 'fftmethod' in result.keys():
                fftmethod = result['fftmethod']
                for kk in fftmethod.keys():
                    prop=self.tabulated_panel.ffttree.GetPropertyByName(kk)
                    prop.SetValue(fftmethod[kk])
    
            self.Sys.setFromCtrl(self.tabulated_panel.Sys_param.parameters) # This is more to double check that input works
            self.Exp.setFromCtrl(self.tabulated_panel.Exp_param.parameters)                            
            self.Opt.setFromCtrl(self.tabulated_panel.Opt_param.parameters)  
            self.tabulated_panel.update_filelist()    
            self.matplotlib_panel.update_graph() 
            
            if type(self.Data[0]['simdata'])!=type(None):
                dlg = wx.MessageDialog(self,
                        "The loaded session contains a simulation resut. \t Would you like to refresh it by running simulaiton again?",
                        "Actualize simulation", wx.YES_NO | wx.ICON_QUESTION)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result == wx.ID_YES:
                    self.runSim()
                
        except Exception:
            wx.MessageBox(f"Loading {filename} failed.", "Cannot load",
                          wx.OK | wx.ICON_INFORMATION)
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
                except Exception:
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
