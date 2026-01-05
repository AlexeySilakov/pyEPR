import os
import wx
import wx.lib.mixins.listctrl as listmix
import wx.grid as gridlib
import wx.propgrid as wxpg
import numpy as np

import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg 

import brukerread as Mybr
import classPropGridPanel as MypgPanel
from SysPar import sysPar, expPar
from EPR_sim import optEPR, EPRsim

import wx.lib.agw.flatnotebook as fnb

import time


VERSION=0.2

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
# Matplotlib Canvas Panel
# --------------------------------------------------------------------------- #
class MatplotlibPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kw):
        super().__init__(parent, *args, **kw)
        self.parent = mainWindow
        self.showSim = True
        self.showBG = False
        self.showDiff = False
        self.showComp = False

        self.figure = Figure(figsize=(2, 2), dpi=100)
        self.axes = []
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        
        self.btn_LoaddFile = wx.Button(self, label="Load File")
        self.btn_LoaddFile.Bind(wx.EVT_BUTTON, self.parent.on_loadfile)
        
        
        self.chk_BG = wx.CheckBox(self, label="Background")
        self.chk_BG.SetValue(self.showBG)
        self.chk_BG.Bind(wx.EVT_CHECKBOX, self.on_check)
        
        self.chk_SIM = wx.CheckBox(self, label="Simulation")
        self.chk_SIM.SetValue(self.showSim)
        self.chk_SIM.Bind(wx.EVT_CHECKBOX, self.on_check)

        self.chk_DIFF = wx.CheckBox(self, label="Difference")
        self.chk_DIFF.SetValue(self.showDiff)
        self.chk_DIFF.Bind(wx.EVT_CHECKBOX, self.on_check)

        self.chk_COMP = wx.CheckBox(self, label="Components")
        self.chk_COMP.SetValue(self.showComp)
        self.chk_COMP.Bind(wx.EVT_CHECKBOX, self.on_check)

        sizerH = wx.BoxSizer(wx.HORIZONTAL)
        sizerH.Add(self.btn_LoaddFile, 0, wx.EXPAND)
        sizerH.AddSpacer(20)
        sizerH.Add(self.chk_BG, 0, wx.EXPAND)
        sizerH.Add(self.chk_SIM, 0, wx.EXPAND)
        sizerH.Add(self.chk_DIFF, 0, wx.EXPAND)
        sizerH.Add(self.chk_COMP, 0, wx.EXPAND)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerH, 0, wx.EXPAND)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)

        sizer.Add(self.toolbar, 0, wx.LEFT|wx.EXPAND)

        
        self.SetSizer(sizer)
        self.toolbar.Realize()
        self.toolbar.update()
    def on_check(self, event): 
        pass
    def update_graph(self):
        shw = []
        if len(self.parent.Data)==0:
            return

        for dd in self.parent.Data:
            shw.append(dd['show'])

        nData = shw.count(True)
        makenew = False
        if len(self.axes)!=nData:
            for ax in self.axes:
                if ax is not None:
                    ax.clear()
                    ax.remove()
            self.axes = [None]*(nData)
            self.figure.clf()
            makenew = True
        axcnt= 0
        showLabels = True
        
        if self.parent.ScaleModeChoice == "Fit Peak-to-Peak":
            scalefnc = lambda data: np.real(np.max(data) - np.min(data))
        elif self.parent.ScaleModeChoice == "Fit Max Signal":
            scalefnc = lambda data: np.real(np.max(data))
        elif self.parent.ScaleModeChoice == "Fit Min Signal":
            scalefnc = lambda data: np.abs(np.min(data))
        elif self.parent.ScaleModeChoice == "Fit Double Integral":
            scalefnc = lambda data: np.abs(np.sum(np.cumsum(data)))
        else:
            scalefnc =  lambda data: 1.0
        #["Fit Peak-to-Peak", "Fit Max Signal", "Fit Min Signal", "Fit Double Integral", "none"]
        
        ### first row is always DATA
        for dd in self.parent.Data:
            if not dd['show']: continue
            fname =np.real(dd['fname'])
            tdata = []
            txdata = []
            tydata = []
            tcolor = []
            tlable = []
            drawData = False
            drawSim = False
            drawDiff = False
            datasc = 1.0
            if dd is not None:
                tdata.append(np.real(dd['data']))
                txdata.append(np.real(dd['ax']['x']))
                tlable.append(f'Data: {fname}')
                #tydata.append(np.real(dd['ax']['y']))
                tcolor.append(dd['colour'])
                drawData = True
                datasc = scalefnc(np.real(dd['data']))
                
            ### not sure if that is necessary, but just in case we have actual simualtion without data, I would separate these two
            if dd['simactual'] and self.showSim:
                for ii in range(len(dd['simdata'])):
                    relsc = dd['simscale'][ii]
                    ssc = scalefnc(np.real(dd['simdata'][ii]))
                    tdata.append(np.real(dd['simdata'][ii])*datasc/ssc*relsc)
                    txdata.append(np.real(dd['simax'][ii]['x']))
                    #tydata.append(np.real(dd['ax']['y']))
                    tcolor.append(dd['simcolor'][ii]) ### need a way to define color
                    tlable.append(f'Sim {ii}')
                    drawSim = True
            #### add difference
          
            if makenew:                                             
                self.axes[axcnt]=self.figure.add_subplot(nData, 1, axcnt+1)
                
            ## clean what is not necessary
            for line in self.axes[axcnt].lines:
                if line.get_label() not in tlable:
                    line.remove()
                    
            for ii in range(len(tdata)):
                # drawData = True
                # data = np.real(dd['data'])
                # fname =np.real(dd['fname'])
                # xax = np.real(dd['ax']['x'])
                # wxc = dd['colour']
                r, g, b, a = tcolor[ii].Red(), tcolor[ii].Green(), tcolor[ii].Blue(), tcolor[ii].Alpha()
                color = (r / 255.0, g / 255.0, b / 255.0, a / 255.0)
                hasData = False
                for line in self.axes[axcnt].lines:
                    if tlable[ii] == line.get_label():
                        hasData = True
                        line.set_xdata(txdata[ii])
                        line.set_ydata(tdata[ii])
                        line.set_label(tlable[ii])
                        line.set_color(color)
                        
                if not hasData:
                    data = self.axes[axcnt].plot(txdata[ii], tdata[ii],
                                     label=tlable[ii],
                                     color=color)
                

            axcnt+=1
        self.figure.tight_layout()
        self.canvas.draw()
# --------------------------------------------------------------------------- #
# Tabulated Options Panel
# --------------------------------------------------------------------------- #   
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kwargs):
        super().__init__(parent,size=(200, 400), *args, **kwargs)
        self.parent = mainWindow
        self.Sys_panel = []
        self.Sys_param = []
        self.Color_BG_FILE_Title = wx.Colour(64, 120, 215)
        self.Color_BG_FILE_Inact = wx.Colour(100, 150, 215)
        self.Color_FG_FILE_Title = wx.Colour(255, 255, 255)
        self.Color_BG_FILE_Main = wx.Colour(255, 255, 255)
        self.Color_FG_FILE_Main = wx.Colour(0, 0, 0)
        
        # self.nb = wx.Notebook(self)
        self.nb = fnb.FlatNotebook(self, agwStyle=fnb.FNB_NO_X_BUTTON|fnb.FNB_NO_NAV_BUTTONS|
                                  fnb.FNB_SMART_TABS | fnb.FNB_NODRAG)
        self.nb.SetActiveTabTextColour(wx.Colour("black"))
        self.nb.SetBackgroundColour(wx.Colour("white"))
        #self.nb.Bind(wx.EVT_RIGHT_DOWN, self.OnTabRightClick)
        self.nbMenu = wx.Menu()

        rename_id = wx.ID_ANY
        close_id  = wx.ID_ANY
        self.TAB_ADD_NAME = "Add a system"
        self.TAB_REN_NAME = "Rename system"
        self.TAB_DEL_NAME = "Delete system"
        add_item = self.nbMenu.Append(rename_id, self.TAB_ADD_NAME)
        #rename_item = self.nbMenu.Append(rename_id, self.TAB_REN_NAME)
        del_item  = self.nbMenu.Append(close_id,  self.TAB_DEL_NAME)
        self.nbMenu.Bind(wx.EVT_MENU, self.on_tabcontext)
        #self.nbMenu.Bind(wx.EVT_MENU, self.on_tabcontext)
        # close_item.Bind(wx.EVT_MENU, self.OnCloseTab)
        
        self.nb.SetRightClickMenu(self.nbMenu)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.nb, 1, wx.EXPAND)
        
        sizer_h1 = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_play = wx.Button(self, label="▶", size=(40, 40))  #🢱
        self.btn_play.Bind(wx.EVT_BUTTON, self.parent.on_play)
        sizer_h1.Add(self.btn_play, 0, wx.ALL)
        self.chk_autoSIM = wx.CheckBox(self, label="Auto Sim")
        self.chk_autoSIM.SetValue(False)
        sizer_h1.Add(self.chk_autoSIM, 0, wx.ALL|wx.ALIGN_CENTER) 
        sizer.Add(sizer_h1, 0, wx.EXPAND)
        
        self.cmb_Scale = wx.ComboBox(self, choices=self.parent.ScaleChoices, style=wx.CB_DROPDOWN)
        idx = self.parent.ScaleChoices.index(self.parent.ScaleModeChoice)
        self.cmb_Scale.SetSelection(idx)  # Set default selection
        sizer.Add(self.cmb_Scale, 0, wx.EXPAND)
        self.cmb_Scale.Bind(wx.EVT_COMBOBOX, self.on_Scale_change)
        
        self.lbl_Status = wx.StaticText(self, label="> ")
        sizer.Add(self.lbl_Status, 0, wx.EXPAND)
        self.add_datapanel()
        self.add_syspanel()
        self.add_optpanel()
        self.SetSizer(sizer)
    # def OnTabRightClick(self, event):
    #     print('click')
    #     click_pos = event.GetPosition() event.GetSelection()
    #     page, flags = self.nb.HitTest(click_pos)
    #     if page == wx.NOT_FOUND:
    #         # Click was on the empty area of the notebook (e.g. on the
    #         # background between tabs).  We simply ignore it.
    #         print(flags)
    #         event.Skip()
    #         return
    #     menu = wx.Menu()

    #     rename_id = wx.ID_ANY
    #     close_id  = wx.ID_ANY

    #     rename_item = menu.Append(rename_id, "Rename Tab")
    #     close_item  = menu.Append(close_id,  "Close Tab")

    #     # Bind the menu commands – we capture the *page* value with a lambda
    #     # so the handler knows which tab the user clicked on.
    #     self.Bind(wx.EVT_MENU,
    #               lambda evt, pg=page: self.OnRenameTab(evt, pg),
    #               rename_item)
    #     self.Bind(wx.EVT_MENU,
    #               lambda evt, pg=page: self.OnCloseTab(evt, pg),
    #               close_item)
    def printStatus(self, message):
        self.lbl_Status.SetLabel(f'> {message}')
        
    def on_tabcontext(self, event):
        page = self.nb.GetSelection() #self.nb.GetCurrentPage()
        if self.nb.GetPageText(page) in ['Data', 'Exp', 'Opt']:
            print('cannot change Data or Opt tabs')
            return
        mm = wx.GetMousePosition()
        print(mm)
        idx = event.GetId()
        menu = self.nbMenu.FindItemById(idx).GetItemLabel()
        """Prompt for a new name and set it on the selected tab."""

        
        if menu == self.TAB_REN_NAME:
            
            #idx = getattr(self.nb, "_rightClickedTab", None)
            dlg = wx.TextEntryDialog(self,
                                     "Enter a new title for the tab:",
                                     "Rename Tab")
            if dlg.ShowModal() == wx.ID_OK:
                new_title = dlg.GetValue()
                self.nb.SetPageText(page, new_title)
            dlg.Destroy()
        elif menu == self.TAB_ADD_NAME:
            self.add_syspanel()
        elif menu == self.TAB_DEL_NAME:
            print('TBD')
            ### figure out how to find the index corresponence
            #self.nb.GetPageText(page)
            # use del_syspanel
            
    def del_syspanel(self, idx):
        if idx in range(len(self.Sys[idx])):
            del(self.Sys[idx])
            del(self.Sys_panel[idx])
            del(self.Sys_param[idx])
        
    def add_optpanel(self):
        # ------------- Opt
        
        self.Opt_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.Opt_panel, "Opt") 
        Optsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.Opt_param = MypgPanel.PropGridPanel(self.Opt_panel, Prop_Dict = {}, onChangeFunc=self.on_Opt)
        self.Opt_param.SetFromParClean(self.parent.Opt.getDefaultDict())
        Optsizer.Add(self.Opt_param, 1, wx.EXPAND)
        self.Opt_panel.SetSizer(Optsizer) 
    def add_syspanel(self):
        if not hasattr(self, 'btn_col'):
            self.btn_col = []
        if not hasattr(self, 'btn_add'):
            self.btn_add = []
        if not hasattr(self, 'btn_nucdel'):
            self.btn_nucdel = []  
        if not hasattr(self, 'spn_sca'):
            self.spn_sca = []             
        # ------------- Sys 
        self.parent.Sys.append(sysPar())
        self.Sys_panel.append(wx.Panel(self.nb))
        
        cnt = len(self.Sys_panel)-1
        self.nb.AddPage(self.Sys_panel[cnt], f"Sys_{cnt+1}") 
        
        Syssizer = wx.BoxSizer(wx.VERTICAL)
        Syssizer.AddSpacer(10)
        #----------------------------------------------------------------------
        lbl = wx.StaticText(self.Sys_panel[cnt], label="Color ")
        self.btn_col.append(wx.Button(self.Sys_panel[cnt], label=""))
        self.btn_col[cnt].SetBackgroundColour(wx.Colour('Red'))
        self.btn_col[cnt].Bind(wx.EVT_BUTTON, self.on_SysColor)
        btnszr0 = wx.BoxSizer(wx.HORIZONTAL)
        
        btnszr0.Add(lbl, 0, wx.ALIGN_CENTER_VERTICAL)
        btnszr0.Add(self.btn_col[cnt], 1, wx.EXPAND)
        
        Syssizer.Add(btnszr0, 0, wx.EXPAND)
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        lbl = wx.StaticText(self.Sys_panel[cnt], label="Rel.Scale ")
        self.spn_sca.append(wx.SpinCtrlDouble(self.Sys_panel[cnt],
                                       style=wx.TE_PROCESS_ENTER|wx.TE_CENTER|wx.SP_ARROW_KEYS))
        self.spn_sca[cnt].SetRange(0, 100)
        self.spn_sca[cnt].SetDigits(3)
        self.spn_sca[cnt].SetIncrement(0.01)
        self.spn_sca[cnt].SetValue(1.0)
        self.spn_sca[cnt].Bind(wx.EVT_SPINCTRLDOUBLE, self.on_SysScale)
        btnszr0 = wx.BoxSizer(wx.HORIZONTAL)
        
        btnszr0.Add(lbl, 0, wx.ALIGN_CENTER_VERTICAL)
        btnszr0.Add(self.spn_sca[cnt], 1, wx.EXPAND)
        
        Syssizer.Add(btnszr0, 0, wx.EXPAND)
        #----------------------------------------------------------------------
        
        Syssizer.AddSpacer(10)
        
        btnszr = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_add.append(wx.Button(self.Sys_panel[cnt], label="+ Nuc"))
        #self.btn_add.Bind(wx.EVT_BUTTON, self.on_add_Nuc)
        btnszr.Add(self.btn_add[cnt], 1, wx.EXPAND)
        self.btn_nucdel.append(wx.Button(self.Sys_panel[cnt], label="- Nuc"))
        #self.btn_nucdel.Bind(wx.EVT_BUTTON, self.on_delete_nuc)
        btnszr.Add(self.btn_nucdel[cnt], 1, wx.EXPAND)
       
        Syssizer.Add(btnszr, 0, wx.EXPAND)
        
        self.Sys_param.append(MypgPanel.PropGridPanel(self.Sys_panel[cnt], Prop_Dict = {}, onChangeFunc=self.on_Sys))
        self.Sys_param[cnt].SetFromParClean(self.parent.Sys[cnt].getDefaultDictEPR())
        
        
        Syssizer.Add(self.Sys_param[cnt], 1, wx.EXPAND)
      
        self.Sys_panel[cnt].SetSizer(Syssizer)  
        
    def add_datapanel(self):
        # ----------Data panel
        self.data_panel = wx.Panel(self.nb)
        self.nb.AddPage(self.data_panel, "Data")
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        sizerBtns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_delete = wx.Button(self.data_panel, label="Delete Data")
        #self.btn_delete.Bind(wx.EVT_BUTTON, self.on_delete)
        sizerBtns.Add(self.btn_delete, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        
        self.btn_dsc = wx.Button(self.data_panel, label="Show DSC")
        #self.btn_dsc.Bind(wx.EVT_BUTTON, self.on_dsc)
        sizerBtns.Add(self.btn_dsc, 1, wx.ALIGN_CENTER | wx.ALL, 5)
        sizer.Add(sizerBtns, 0, wx.ALIGN_LEFT | wx.ALL, 5)

        self.chk_expanded = wx.CheckBox(self.data_panel, label="Expand all")
        self.chk_expanded.SetValue(False)
        #self.chk_expanded.Bind(wx.EVT_CHECKBOX, self.on_expandchk)
        sizer.Add(self.chk_expanded, 0, wx.ALIGN_LEFT | wx.ALL, 0)
        # file tree

        self.data_tree = MypgPanel.PropGridPanel(self.data_panel, Prop_Dict={}, onChangeFunc=self.on_Data)
        self.data_tree.SetFromParClean(self.parent.currentExp)
        # self.filetree = wxpg.PropertyGrid(
        #     self.data_panel,
        #     style=wxpg.PG_DEFAULT_STYLE | 
        #           wxpg.PG_SPLITTER_AUTO_CENTER |
        #           wxpg.PG_HIDE_MARGIN |
        #           wxpg.PG_TOOLTIPS |
        #           wxpg.PG_NO_INTERNAL_BORDER
        # )
        
        # self.filetree.SetCellBackgroundColour(self.Color_BG_FILE_Main)
        # self.filetree.SetCellTextColour(self.Color_FG_FILE_Main)
        # self.filetree.SetMarginColour(self.Color_BG_FILE_Main)
        # #self.filetree.SetSelectionBackgroundColour(wx.Colour(200, 255, 255))
        
        # self.filetree.SetCaptionBackgroundColour(self.Color_BG_FILE_Inact)
        # self.filetree.SetCaptionTextColour(self.Color_FG_FILE_Title)
        # cat1 = self.filetree.Append(wxpg.PropertyCategory("Nothing loaded yet"))


        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("MW Freq", value=9.43))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("B0", value=350.0))
        # self.filetree.AppendIn(cat1, wxpg.FloatProperty("tau", value=120.0))
        
        sizer.Add(self.data_tree, 1, wx.EXPAND | wx.ALL, 0)
        
        self.data_panel.SetSizer(sizer)        
    def on_Sys(self, panel, mainname, name, val):
        ### find which Systab is activated
        
        idx = self.get_SysID(panel.GetParent())
        self.parent.Sys[idx].setFromCtrl(panel.parameters) # This is more to double check that input works
        if self.chk_autoSIM.GetValue():
            self.parent.on_play(None)
        
    def on_Data(self, panel, mainname, name, val):
        ### not sure what to do, yet
        #print(name)
        self.parent.update_data_fromtree()
        if name=='Colour':
            self.parent.plot_panel.update_graph()
        # pass
        #self.parent.Exp.setFromCtrl(self.Exp_param.parameters)

    def on_Opt(self, panel, mainname, name, val):
        self.parent.Opt.setFromCtrl(self.Opt_param.parameters)    
        
    def on_Scale_change(self, event):
        self.parent.ScaleModeChoice = self.cmb_Scale.GetValue()
        self.parent.plot_panel.update_graph()
    def get_SysID(self, panel_handle):
        return self.Sys_panel.index(panel_handle)
    def set_tabColor(self, page, color):
        bmp = wx.Bitmap(16, 16)
        dc = wx.MemoryDC(bmp)
        dc.SetBackground(wx.Brush(color))
        dc.Clear()
        dc.SelectObject(wx.NullBitmap)
        self.imglst = wx.ImageList(16, 16, mask=True)
        self.imglst.Add(bmp)
        self.nb.SetImageList(self.imglst)
        self.nb.SetPageImage(page, 0)
    def on_SysColor(self, event):
        
        apr = event.EventObject.GetParent()
        ID = self.get_SysID(apr)
        current_color = self.btn_col[ID].GetBackgroundColour()  
        color_dialog = wx.ColourDialog(self)
               
        # Configure dialog options
        color_data = color_dialog.GetColourData()
        color_data.SetChooseFull(True)  # Allow full color selection
        color_data.SetColour(current_color)  # Set initial color

                # Add custom colors (you can pre-populate these)
        custom_colors = [
            wx.Colour(255, 0, 0),    # Red
            wx.Colour(0, 255, 0),    # Green
            wx.Colour(0, 0, 255),    # Blue
            wx.Colour(255, 255, 0),  # Yellow
            wx.Colour(255, 0, 255),  # Magenta
            wx.Colour(0, 255, 255),  # Cyan
            wx.Colour(128, 128, 128), # Gray
            wx.Colour(0, 0, 0),      # Black
        ]
        
        for ii, color in enumerate(custom_colors):
            if ii < 16:  # Max 16 custom colors
                color_data.SetCustomColour(ii, color)
                
        # Show dialog
        if color_dialog.ShowModal() == wx.ID_OK:
            # Get the selected color
            selected_color = color_dialog.GetColourData().GetColour()
            self.btn_col[ID].SetBackgroundColour(selected_color)#SetBackgroundColour
            self.Sys_panel[ID].SetBackgroundColour(selected_color)
            #self.update_color_display(selected_color)
            #self.status_text.SetLabel(f"Basic picker: {selected_color.GetAsString(wx.C2S_HTML_SYNTAX)}")
            pan = event.EventObject.GetParent()
            page = pan.GetParent().GetSelection()
            self.set_tabColor(page, selected_color)
            if len(self.parent.Data[0]['simdata'])>=(ID-1):
                self.parent.Data[0]['simcolor'][ID] = selected_color
                self.parent.plot_panel.update_graph()
                
        color_dialog.Destroy()
    def on_SysScale(self, event):
        apr = event.EventObject.GetParent()
        ID = self.get_SysID(apr)
        
        val = float(self.spn_sca[ID].GetValue())
        if len(self.parent.Data[0]['simdata'])>=(ID-1):
            self.parent.Data[0]['simscale'][ID] = val
            self.parent.plot_panel.update_graph()

# --------------------------------------------------------------------------- #
# Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.Data = []
        self.SetTitle(f"pyEPR {VERSION}")
        self.SetSize((1000, 600))
        self.Sys = []
        self.Exp = expPar()
        self.Opt = optEPR()
        self.currentExp = {}
        self.ScaleChoices = ["Fit Peak-to-Peak", "Fit Max Signal", "Fit Min Signal", "Fit Double Integral", "none"]
        self.ScaleModeChoice = "Fit Peak-to-Peak"
        path, fname = os.path.split(os.path.abspath(__file__))
        self.currentPath = path
        self.ini_file = os.path.join(path, "pyEPR.ini")

        self.BackgroundColor = wx.Colour(255, 255, 255)
        self.ForegroundColor = wx.Colour(0, 0, 0)
        self.ControlBackgroundColor = wx.Colour(255, 255, 255)
        
        right_sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.splitter_right = wx.SplitterWindow(self, style=wx.SP_3D)
        
        
        self.plot_panel = MatplotlibPanel(self.splitter_right, self)
        self.param_panel   = TabulatedPanel(self.splitter_right, self)
        
        # Split the right panel horizontally
        self.splitter_right.SplitVertically(
            self.plot_panel, self.param_panel,
            sashPosition=400   # initial height of the top pane
        )
        self.splitter_right.SetSashGravity(1)  # proportion of space for the top pane
        
        right_sizer.Add(self.splitter_right, 1, wx.EXPAND)
        
        self.SetSizer(right_sizer)

        self.Centre()
        self.Bind(wx.EVT_WINDOW_DESTROY, self.on_destroy)   # persistence

        self._load_path_from_file() ### make it more general
        self.update_datatree()
    # ----------------------------------------------------------------------
    # (plain‑text .ini file)
    # ----------------------------------------------------------------------
    def _load_path_from_file(self):
        """Read the first line of ~/pyEPR.ini; if it is a valid
        directory, use it."""
        if not os.path.exists(self.ini_file):
            path, fname = os.path.split(os.path.abspath(__file__))
            self.currentPath = path
            return
        try:
            with open(self.ini_file, "r", encoding="utf-8") as f:
                stored = f.readline().strip()
        except Exception:
            return
       
        if stored and os.path.isdir(stored):
            # `set_path` will update lbl_path and list for us.
            self.currentPath = stored
        else:
            path, fname = os.path.split(os.path.abspath(__file__))
            self.currentPath = path
    def _save_path_to_file(self):
        """Write the current path to the file."""
        try:
            with open(self.ini_file, "w", encoding="utf-8") as f:
                f.write(self.currentPath + "\n")
        except Exception:
            pass
    def on_destroy(self, event):
        self._save_path_to_file()


    # -----------------------------------------------------------------------
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

        raise AttributeError(message)
    def on_loadfile(self, evnt):
        file_path = ''

        with wx.FileDialog(self,
            message="Select a file to load",
            wildcard="Text files (*.csv)|*.csv|All files (*.*)|*.*",
            defaultDir=self.currentPath,
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as dlg:
    
            if dlg.ShowModal() == wx.ID_OK:
                file_path = dlg.GetPath()
            else:
                return
                # try:
                #     with open(file_path, "r", encoding="utf-8") as f:
                #         content = f.read()
                #     self.text.SetValue(content)
                # except Exception as exc:
                #     wx.MessageBox(f"Failed to read file:\n{exc}",
                #                   "Error",
                #                   wx.OK | wx.ICON_ERROR)
        _, ext = os.path.splitext(file_path)
        path, fname = os.path.split(file_path)
        self.currentPath =path
        ext = ext.lower()
        data = None
        ####  add an attempt to read it as a standard ascii file
        ###########################################################################
        if ext == '.csv':
            from eprcsvread import csvread
            ax, data, dsc = csvread(file_path, return_ax=True, return_dsc=True)
            freq = ax['freq1']
        ###########################################################################
        elif ext in ('.dsc', '.dta'):
            import brukerread as Mybr
            ax, data, dsc = Mybr.brukerread(file_path, return_ax=True, return_dsc=True)
            freq = ax['freq1']*1e-9
        ###########################################################################
        if type(data)==type(None):
            self.RaiseError(f'File {file_path} does not contain expected data structure')
        else:
            tData ={'ax':ax, 'data':data, 'dsc':dsc, 'fname':fname,
                     'mwFreq':freq, 'BMin': float(np.min(ax['x'])), 'BMax': float(np.max(ax['x'])),
                     'nPoints':max(ax['x'].shape), 'colour':wx.Colour(wx.BLUE),
                     'fullpath':path, 'show':True, 'scale':1.0, 'temp':-0.1,
                     'simdata':None, 'simax':None, 'simmethod':None, 'simactual':False,
                     }
            if len(self.Data)==0:
                self.Data.append(tData)
            else:
                # just so we just replace data without affecting anything else
                tData['colour']=self.Data[0]['colour']
                tData['simdata']=self.Data[0]['simdata']
                tData['simax']=self.Data[0]['simax']
                tData['simmethod']=self.Data[0]['simmethod']
                tData['simactual']=self.Data[0]['simactual']
                tData['show']=self.Data[0]['show']
                self.Data[0] = tData

        # print(file_path)
        #
        #

        self.update_datatree()
        self.plot_panel.update_graph()
    def update_datatree(self):
        if len(self.Data)==0:
            self.currentExp = {'no data loaded':
                               {'Colour':wx.Colour(wx.BLUE),'Freq [GHz]':9.500,
                                'Bmin [mT]':float(300.0),'Bmax [mT]':float(400.0), 'nPoints':int(1000),
                                'Temp. [K]':-1.0}
                           }
        else:
            td = self.Data[0]

            self.currentExp = {td['fname']:
                               {'Colour':td['colour'],'Freq [GHz]':td['mwFreq'],
                                'Bmin [mT]':float(td['BMin']),'Bmax [mT]':float(td['BMax']), 'nPoints':td['nPoints'],
                                'Temp. [K]':td['temp']}
                           }
        self.param_panel.data_tree.SetFromParClean(self.currentExp)
    def update_data_fromtree(self):
        self.currentExp = self.param_panel.data_tree.parameters
        cnt = 0
        td = {}
        for kk in self.currentExp.keys():
            dd = self.currentExp[kk]
            if len(self.Data)<(cnt+1):
                self.Data.append({'ax':None, 'data':None, 'dsc':None, 'fname':'no data loaded',
                         'mwFreq':freq, 'BMin': float(np.min(ax['x'])), 'BMax': float(np.max(ax['x'])),
                         'nPoints':max(ax['x'].shape), 'colour':wx.Colour(wx.BLUE),
                         'fullpath':path, 'show':True, 'scale':1.0, 'temp':-0.1,
                         'simdata':[], 'simax':[], 'simmethod':[], 'simactual':False, 
                         'simcolor':[], 'simscale': [],
                         }
                         )
            self.Data[cnt]['fname']=kk
            self.Data[cnt]['colour']=dd['Colour']
            self.Data[cnt]['mwFreq']=dd['Freq [GHz]']
            self.Data[cnt]['BMin']=dd['Bmin [mT]']
            self.Data[cnt]['BMax']=dd['Bmax [mT]']
            self.Data[cnt]['nPoints']=dd['nPoints']
            self.Data[cnt]['temp']=dd['Temp. [K]']

            cnt+=1

    def on_play(self, evnt):
        for ii in range(len(self.Sys)):
            self.Sys[ii].setFromCtrl(self.param_panel.Sys_param[ii].parameters)
        self.update_data_fromtree()
        self.runSim()
        self.plot_panel.update_graph()
    def runSim(self):
        self.Data[0]['simax'] = []
        self.Data[0]['simmethod'] = []
        self.Data[0]['simdata'] = []
        self.Data[0]['simcolor'] = []
        self.Data[0]['simscale'] = []
        self.Data[0]['simactual']=False
        try:
            #### TBD: do multiple simulations 
            es = EPRsim()
            
            for ii in range(len(self.Sys)):
                es.Sys = self.Sys[ii]
                # es.Sys.set( S = [1/2],
                #         g = [np.array([1.981, 1.979, 1.944])], 
                #         Nucs = ['51V'],
                #         A = [np.array([519, 185, 192])],
                #         lw = [1],
                #         )
        
                es.Exp.set(mwFreq=self.Data[0]['mwFreq'],
                           BMin=self.Data[0]['BMin'],
                           BMax=self.Data[0]['BMax'],
                           nPoints=self.Data[0]['nPoints'],
                           Harmonic=1)
                timest = time.time()
                es.run()
                timeed = time.time()-timest
                if timeed<0.1:
                    sc = 1000
                    un = 'ms'
                else:
                    sc = 1
                    un = 's'
                self.param_panel.printStatus(f'time: {timeed*sc:.3} {un}')
                # plt.figure(2)
                # plt.clf()
                # plt.plot(es.rawX, es.rawY,'r')
                # plt.plot(es.X, es.Y,'b')
                # plt.show()
                if 'ax' in self.Data[0].keys():
                    self.Data[0]['simax'].append(self.Data[0]['ax'].copy())
                else:
                    self.Data[0]['simax'].append({'x':es.X.copy(), 'xlable':'Magnetic Field, mT'})
                    # self.Data[0]['simax']['x']=es.X.copy()
                    # self.Data[0]['simax']['xlable']='Magnetic Field, mT'
        
                simmethod = {'Sys':self.param_panel.Sys_param[ii].parameters.copy(),
                             'Exp':{'mwFreq':self.Data[0]['mwFreq'],
                                        'BMin':self.Data[0]['BMin'],
                                        'BMax':self.Data[0]['BMax'],
                                        'nPoints':self.Data[0]['nPoints'],
                                        'Harmonic':1}
                             }
                self.Data[0]['simmethod'].append(simmethod.copy())
                self.Data[0]['simactual']=True
                self.Data[0]['simdata'].append(es.Y.copy())
                self.Data[0]['simcolor'].append(self.param_panel.btn_col[ii].GetBackgroundColour())
                self.Data[0]['simscale'].append(float(self.param_panel.spn_sca[ii].GetValue()))
        except Exception:
            return

if __name__ == "__main__":
    wx.SystemOptions.SetOption("msw.dpiAware", "1")
    app = wx.App(False)

    frame = MainFrame(None)
    frame.Show(True)

    app.MainLoop()        