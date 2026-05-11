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
import wx.lib.agw.aui as aui

import time
import brukerread as Mybr

# ############### for debugging purposes
# wx.Log.SetActiveTarget(wx.LogStderr())
# wx.Log.SetLogLevel(3)  # verbose
# import faulthandler
# faulthandler.enable()

# GIT commands :
#git add *   # stage all files
#git commit -m "message"  
#git status

VERSION=0.6

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
        
# Matplotlib Canvas Panel
class MatplotlibPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kw):
        super().__init__(parent, *args, **kw)
        self.parent = mainWindow
        self.showSim = True
        self.showBG = False
        self.showDiff = False
        self.showComp = True
        self.measure_mode = False
        self.stored_xy = None
        self.figure = Figure(figsize=(2, 2), dpi=100)
        self.axes = []
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        
        self.canvas.mpl_connect("pick_event", self.on_pick)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        #self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect("button_press_event", self.on_mouse_click)
        
        # self.canvas.mpl_connect("axes_enter_event", self.on_axes)
        self._dragging_line = None
        self._press_ydata = None
        self._orig_ydata = None
        self._dy = 0.0
        self._line_label = None
        self.shifts = {}
        
        menu_btnsize = 32
        
        self.bmp_off = self.parent.BMP_OFF
        self.bmp_on  = self.parent.BMP_ON
        
        self.toolbar = wx.ToolBar(self, style=wx.TB_HORIZONTAL | wx.TB_TEXT | wx.TB_DEFAULT_STYLE | wx.TB_FLAT)
        self.btn_LoaddFile = self.toolbar.AddTool(wx.ID_ANY, "Data", self.parent.BMP_LOADEPR, shortHelp ="Load Data")
        self.Bind(wx.EVT_TOOL, self.parent.on_loadfile, self.btn_LoaddFile)
        
        self.btn_LoaddFile = self.toolbar.AddTool(wx.ID_ANY, "Session", self.parent.BMP_LOADXML, shortHelp ="Load Session from xml")
        self.Bind(wx.EVT_TOOL, self.parent.on_loadSession, self.btn_LoaddFile)
        
        self.btn_LoaddFile = self.toolbar.AddTool(wx.ID_ANY, "Session", self.parent.BMP_SAVEXML, shortHelp ="Save Session as xml")
        self.Bind(wx.EVT_TOOL, self.parent.on_saveSession, self.btn_LoaddFile)
        
        self.sep1 = self.toolbar.AddSeparator()
        
        self.btn_Zoom = self.toolbar.AddTool(wx.ID_ANY, "Zoom", self.parent.BMP_ZOOM, shortHelp ="Zoom on data")
        self.Bind(wx.EVT_TOOL, self.on_zoom, self.btn_Zoom)
        
        self.btn_Reset = self.toolbar.AddTool(wx.ID_ANY, "Reset", self.parent.BMP_HOME, shortHelp ="Reset plot scale")
        self.Bind(wx.EVT_TOOL, self.on_reset, self.btn_Reset)
        
        self.sep2 = self.toolbar.AddSeparator()
        
        self.chk_SIM = self.toolbar.AddCheckTool(wx.ID_ANY, "Sim.", self.bmp_off, shortHelp ="Show complete simulation")
        self.toolbar.ToggleTool(self.chk_SIM.GetId(), self.showSim)
        self.update_check(self.chk_SIM, self.showSim)
        self.Bind(wx.EVT_TOOL, self.on_check, self.chk_SIM)


        self.chk_DIFF = self.toolbar.AddCheckTool(wx.ID_ANY, "Diff.", self.bmp_off, shortHelp ="Show Data-Sim")
        self.toolbar.ToggleTool(self.chk_DIFF.GetId(), self.showDiff)
        self.update_check(self.chk_DIFF, self.showDiff)
        self.Bind(wx.EVT_TOOL, self.on_check, self.chk_DIFF)

        self.chk_COMP = self.toolbar.AddCheckTool(wx.ID_ANY, "Comp.", self.bmp_off, shortHelp ="Show all individual sim. components")
        self.toolbar.ToggleTool(self.chk_COMP.GetId(), self.showComp)
        self.update_check(self.chk_COMP, self.showComp)
        self.Bind(wx.EVT_TOOL, self.on_check, self.chk_COMP)

        self.toolbar.Realize()
        
        self.Bind(wx.EVT_CHAR_HOOK, self.on_key)
        
        txtsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.txt_1 = wx.StaticText(self, wx.ID_ANY, '', style = wx.BORDER_NONE) #BORDER_SUNKEN
        self.txt_2 = wx.StaticText(self, wx.ID_ANY, '', style = wx.BORDER_NONE)
        self.txt_3 = wx.StaticText(self, wx.ID_ANY, '', style = wx.BORDER_NONE)
        txtsizer.Add(self.txt_1, 1, wx.EXPAND)
        txtsizer.Add(self.txt_2, 1, wx.EXPAND)
        txtsizer.Add(self.txt_3, 1, wx.EXPAND)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        # sizer.Add(sizerH, 0, wx.EXPAND)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(txtsizer, 0, wx.EXPAND)
        
        self.MLtool = NavigationToolbar2WxAgg(self.canvas)
        self.MLtool.Hide()
        


        self.SetSizer(sizer)
        #self.toolbar.Realize()
        #self.toolbar.update()
    def update_check(self, ctrl, value):
        if value:
            self.toolbar.SetToolNormalBitmap(ctrl.GetId(), self.bmp_on)
        else:
            self.toolbar.SetToolNormalBitmap(ctrl.GetId(), self.bmp_off)
    # ------------------------------------------------------------------
    # Toolbar actions
    # ------------------------------------------------------------------
    def on_zoom(self, event):
        #self._disable_measure()
        self.MLtool.zoom()

    def on_pan(self, event):
        #self._disable_measure()
        self.MLtool.pan()

    def on_reset(self, event):
        self.MLtool.home()

    def on_measure_toggle(self, event):
        self._cancel_navigation_modes()
        self.measure_mode = not self.measure_mode
        self.clear_measurement()
        self.SetStatusText("Measure mode ON" if self.measure_mode else "")

    # ------------------------------------------------------------------
    # Keyboard
    # ------------------------------------------------------------------
    def on_key(self, event):
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self._cancel_navigation_modes()
            #self._disable_measure()
        else:
            event.Skip()

    def _cancel_navigation_modes(self):
        if self.MLtool.mode:
            if "zoom" in self.MLtool.mode.lower():
                self.MLtool.zoom()
            elif "pan" in self.MLtool.mode.lower():
                self.MLtool.pan()

    def _disable_measure(self):
        self.measure_mode = False
        self.clear_measurement()

    # ------------------------------------------------------------------
    # Mouse
    # ------------------------------------------------------------------
#    def on_mouse_move(self, event):


    def on_mouse_click(self, event):
        if event.button==2:
            if event.inaxes:
                self.txt_3.SetLabel(f"stored: x={event.xdata:.3f}, y={event.ydata:.3f}")
                self.stored_xy = [event.xdata, event.ydata]
        if event.dblclick:
            if event.button>1:
                self.stored_xy = None
                self.txt_3.SetLabel(" ")
            else:
                self.MLtool.home()
            
        if not self.measure_mode or not event.inaxes or event.button != 1:
            return

        self.measure_points.append((event.xdata, event.ydata))

        if len(self.measure_points) == 2:
            self._draw_measurement()
            self.measure_points.clear()

    # ------------------------------------------------------------------
    # Measurement logic
    # ------------------------------------------------------------------
    def _draw_measurement(self):
        (x1, y1), (x2, y2) = self.measure_points

        dx = x2 - x1
        dy = y2 - y1
        dist = np.hypot(dx, dy)

        # Draw line
        line, = self.ax.plot([x1, x2], [y1, y2], "r--", lw=2)

        # Annotation
        label = (
            f"ΔX = {dx:.3f}\n"
            f"ΔY = {dy:.3f}\n"
            f"D = {dist:.3f}"
        )

        text = self.ax.text(
            (x1 + x2) / 2,
            (y1 + y2) / 2,
            label,
            color="red",
            bbox=dict(facecolor="white", alpha=0.7)
        )

        self.measure_artists.extend([line, text])
        self.canvas.draw()

    def clear_measurement(self):
        for artist in self.measure_artists:
            artist.remove()
        self.measure_artists.clear()
        self.canvas.draw_idle()            
    def on_pick(self, event):
        # Called when a line is clicked.
        if not isinstance(event.artist, type(self.axes[0].lines[0])):
            return
        if event.mouseevent.button==1:
            self._dragging_line = event.artist
            self._press_ydata = event.mouseevent.ydata
            self._orig_ydata = self._dragging_line.get_ydata().copy()
            self._line_label =self._dragging_line.get_label()
        
    def on_motion(self, event):

        # Called when mouse is moved.
        if self._dragging_line is not None:
            if event.ydata is None:
                return

            self._dy = event.ydata - self._press_ydata
            new_y = self._orig_ydata + self._dy

            self._dragging_line.set_ydata(new_y)
            self.canvas.draw_idle()
        else:
            if event.inaxes:
                if self.stored_xy is not None:
                    self.txt_1.SetLabel(f"  x: {event.xdata:.3f}, y: {event.ydata:.3f}")
                    self.txt_2.SetLabel(f"Δx: {event.xdata-self.stored_xy[0]:.3f}, Δy: {event.ydata-self.stored_xy[1]:.3f}")
                else:
                    self.txt_1.SetLabel(f"x: {event.xdata:.3f}, y: {event.ydata:.3f}")
                    self.txt_2.SetLabel("")
            else:
                self.txt_1.SetLabel("x: ---, y: ---")
                self.txt_2.SetLabel("")

    def on_release(self, event):
        # Called when mouse button is released. Also called on zoom, so be careful
        if self._dragging_line is not None:
            self.shifts[self._line_label]+=self._dy
            self._dragging_line = None
            self._press_ydata = None
            self._orig_ydata = None
            self._line_label = None
            self._dy = 0.0
        
    def on_check(self, event): 
        #self.chk_BG.SetValue(self.showBG)
        self.showSim = self.toolbar.GetToolState(self.chk_SIM.GetId())
        self.update_check(self.chk_SIM, self.showSim)
        
        self.showDiff = self.toolbar.GetToolState(self.chk_DIFF.GetId())
        self.update_check(self.chk_DIFF, self.showDiff)
        
        self.showComp = self.toolbar.GetToolState(self.chk_COMP.GetId())
        self.update_check(self.chk_COMP, self.showComp)
        
        #self.showSim = self.chk_SIM.GetValue()
        # self.showDiff = self.chk_DIFF.GetValue()
        # self.showComp = self.chk_COMP.GetValue()
        
        self.update_graph()
        
    def update_graph(self):
        ###### Tought about current algorithm: 
        ###### So far code below works quick, but I wonder if one can make it faster 
        ###### if we don't redraw what has not changed
        
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
        
        #["Fit Peak-to-Peak", "Fit Max Signal", "Fit Min Signal", "Fit Double Integral", "none"]
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
        
        ### first row is always DATA
        for dd in self.parent.Data:
            if not dd['show']: continue
            fname =np.real(dd['fname'])
            tdata = [] # list of 1D data 
            txdata = [] # list of X axes 
            tydata = [] # list of Y axes
            tcolor = [] # list of colors as wx.Colour 
            tlable = [] # list of labels, acts as data identifiers, so those should be unique
            drawData = False
            drawSim = False
            drawDiff = False
            datasc = 1.0 # used as a standard scale to scale simulations to
            if dd is not None:
                tdata.append(np.real(dd['data']))
                txdata.append(np.real(dd['ax']['x']))
                tlable.append(f'Data: {fname}')
                #tydata.append(np.real(dd['ax']['y']))
                tcolor.append(dd['colour'])
                drawData = True
                datasc = scalefnc(np.real(dd['data']))
                
            ### not sure if that is necessary, but just in case we have actual simualtion without data, I would separate these two
            if dd['simactual']:
                cnt = len(tdata)
                ttdata = []
                ttxdata = []
                ttydata = [] 
                ttcolor = []
                ttlable = []
                
                
                for ii in range(len(dd['simdata'])):
                    relsc = dd['simscale'][ii]
                    ssc = scalefnc(np.real(dd['simdata'][ii]))
                    ttdata.append(np.real(dd['simdata'][ii])*datasc/ssc*relsc)
                    ttxdata.append(np.real(dd['simax'][ii]['x']))
                    #tydata.append(np.real(dd['ax']['y']))
                    ttcolor.append(dd['simcolor'][ii]) 
                    ttlable.append(f'Sim {ii}')

                    if self.showSim or self.showDiff: 
                        if ii == 0:
                            tsum = ttdata[ii].copy()
                        else:
                            tsum += ttdata[ii]

                if self.showSim:
                    tdata.append(tsum)
                    txdata.append(np.real(dd['simax'][0]['x']))
                    if (len(ttdata)>1):
                        tcolor.append(self.parent.MPL_SUM_COLOR) 
                        tlable.append('Simulation')
                    else:
                        tcolor.append(ttcolor[0])
                        tlable.append(ttlable[0])
                    

                if self.showDiff and (dd is not None) and (len(dd['simdata'])>0):
                    if dd['data'].shape==tsum.shape:
                        tt = dd['data']-tsum
                        tdata.append(tt)
                        txdata.append(np.real(dd['simax'][0]['x']))
                        tcolor.append(self.parent.MPL_DIFF_COLOR) 
                        tlable.append('difference')
                        
                if self.showComp and (len(dd['simdata'])>1):
                    tdata+=ttdata ### not sure if this is the most pythonian way of adding lists together
                    txdata+=ttxdata
                    tcolor+=ttcolor
                    tlable+=ttlable
                    
            if makenew:                                             
                self.axes[axcnt]=self.figure.add_subplot(nData, 1, axcnt+1)
                
            for ii in range(len(tdata)):
                try:
                    r, g, b, a = tcolor[ii].Red(), tcolor[ii].Green(), tcolor[ii].Blue(), tcolor[ii].Alpha()
                except:
                    print(tcolor)
                    r, g, b, a = 0, 0, 0, 255
                color = (r / 255.0, g / 255.0, b / 255.0, a / 255.0)
                hasData = False
                for line in self.axes[axcnt].lines:
                    if tlable[ii] == line.get_label():
                        hasData = True
                        dy = self.shifts[tlable[ii]]
                        line.set_xdata(txdata[ii])
                        line.set_ydata(tdata[ii]+dy)
                        line.set_label(tlable[ii])
                        line.set_color(color)
                ##  adding new data if it is not there
                if not hasData:
                    if not tlable[ii] in self.shifts:
                        self.shifts[tlable[ii]] = 0.0
                    dy = self.shifts[tlable[ii]]
                    data = self.axes[axcnt].plot(txdata[ii], tdata[ii]+dy,
                                     label=tlable[ii],
                                     color=color,
                                     picker=5)

            ## clean what is not supposed to be there anymore
            for line in self.axes[axcnt].lines:
                if line.get_label() not in tlable:
                    line.remove()  

            axcnt+=1
        self.figure.tight_layout()
        self.canvas.draw()
 
# Tabulated Options Panel
class TabulatedPanel(wx.Panel):
    def __init__(self, parent, mainWindow, *args, **kwargs):
        super().__init__(parent,size=(300, 400), *args, **kwargs)
        self.parent = mainWindow
        self.Sys_panel = []
        self.Sys_param = []
        self.Color_BG_FILE_Title = wx.Colour(64, 120, 215)
        self.Color_BG_FILE_Inact = wx.Colour(100, 150, 215)
        self.Color_FG_FILE_Title = wx.Colour(255, 255, 255)
        self.Color_BG_FILE_Main = wx.Colour(255, 255, 255)
        self.Color_FG_FILE_Main = wx.Colour(0, 0, 0)
        self.PAGE_IMG_SIZE = 10
        self.PageImgList = wx.ImageList(self.PAGE_IMG_SIZE, self.PAGE_IMG_SIZE, mask=False)
        
        # self.nb = wx.Notebook(self)
        self.nb = fnb.FlatNotebook(self, agwStyle=fnb.FNB_NO_X_BUTTON|fnb.FNB_NO_NAV_BUTTONS|
                                  fnb.FNB_SMART_TABS) #  | fnb.FNB_NODRAG
        # self.nb = aui.AuiNotebook(
        #     self,
        #     style=(
        #         aui.AUI_NB_DEFAULT_STYLE |
        #         aui.AUI_NB_TAB_MOVE |      # <-- enable drag reordering
        #         aui.AUI_NB_SCROLL_BUTTONS
        #     )
        # )
        self.nb.SetActiveTabTextColour(wx.Colour("black"))
        self.nb.SetBackgroundColour(wx.Colour("white"))
        self.nb.SetImageList(self.PageImgList)
        #self.nb.Bind(wx.EVT_RIGHT_DOWN, self.OnTabRightClick)
        self.nbMenu = wx.Menu()

        rename_id = wx.ID_ANY
        close_id  = wx.ID_ANY
        self.TAB_ADD_NAME = "Add a system"
        self.TAB_DUP_NAME = "Duplicate system"
        self.TAB_REN_NAME = "Rename system"
        self.TAB_DEL_NAME = "Delete system"
        add_item = self.nbMenu.Append(rename_id, self.TAB_ADD_NAME)
        add_item = self.nbMenu.Append(rename_id, self.TAB_DUP_NAME)
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
        self.add_optpanel()
        self.add_syspanel()
        
        self.SetSizer(sizer)
    def printStatus(self, message):
        self.lbl_Status.SetLabel(f'> {message}')
        
    def on_tabcontext(self, event):
        page = self.nb.GetSelection() #self.nb.GetCurrentPage()
        if self.nb.GetPageText(page) in ['Data', 'Exp', 'Opt']:
            print('cannot change Data or Opt tabs')
            return
        #mm = wx.GetMousePosition()
        #print(mm)
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
        elif menu == self.TAB_DUP_NAME:
            # self.parameters
            sysidx = self.Sys_panel.index(self.nb.GetPage(page))
            params = self.Sys_param[sysidx].parameters.copy()
            cnt = len(self.Sys_param)
            
            self.add_syspanel()
            self.Sys_param[cnt].SetFromParClean(params)
            pass
        elif menu == self.TAB_DEL_NAME:
            self.del_syspanel(page)
            # for ii in range(self.nb.GetPageCount()):
            #     self.nb.GetPageText(ii)
            # print('TBD')
            ### figure out how to find the index corresponence
            #self.nb.GetPageText(page)
            # use del_syspanel
            
    def del_syspanel(self, page):
        sysidx = self.Sys_panel.index(self.nb.GetPage(page))
        if sysidx in range(len(self.parent.Sys)):
            ### need a more standartized ways to delete all children
            del(self.parent.Sys[sysidx])
            del(self.Sys_panel[sysidx])
            del(self.Sys_param[sysidx])
            del(self.btn_col[sysidx])
            del(self.btn_add[sysidx])
            del(self.btn_addPar[sysidx])
            del(self.btn_nucdel[sysidx]) 
            del(self.spn_sca[sysidx])
        self.nb.DeletePage(page)

    def add_optpanel(self):
        bmp = wx.Bitmap(self.PAGE_IMG_SIZE, self.PAGE_IMG_SIZE)
        self.PageImgList.Add(bmp)
        
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
        if not hasattr(self, 'btn_addPar'):
            self.btn_addPar = []            
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
        self.btn_addPar.append(wx.Button(self.Sys_panel[cnt], label="Par", size=(45, -1)))
        self.btn_addPar[cnt].SetBitmap(self.parent.BMP_PLUS)
        self.btn_addPar[cnt].Bind(wx.EVT_BUTTON, self.on_addPar)
        btnszr.Add(self.btn_addPar[cnt], 0, wx.ALL)
        
        btnszr.AddSpacer(10)
        self.btn_add.append(wx.Button(self.Sys_panel[cnt], label="Nuc", size=(45, -1)))
        self.btn_add[cnt].Bind(wx.EVT_BUTTON, self.on_add_Nuc)
        self.btn_add[cnt].SetBitmap(self.parent.BMP_PLUS)
        btnszr.Add(self.btn_add[cnt], 0, wx.ALL)
        self.btn_nucdel.append(wx.Button(self.Sys_panel[cnt], label="Nuc", size=(45, -1)))
        self.btn_nucdel[cnt].SetBitmap(self.parent.BMP_MINUS)
        self.btn_nucdel[cnt].Bind(wx.EVT_BUTTON, self.on_delete_nuc)
        btnszr.Add(self.btn_nucdel[cnt], 0, wx.ALL)
       
        Syssizer.Add(btnszr, 0, wx.EXPAND)
        
        self.Sys_param.append(MypgPanel.PropGridPanel(self.Sys_panel[cnt], Prop_Dict = {},
                                                      onChangeFunc=self.on_Sys, 
                                                      showModFunc = True))
        self.Sys_param[cnt].SetFromParClean(self.parent.Sys[cnt].getDefaultDictEPR())

        self.Sys_param[cnt].pg.Bind(wx.EVT_RIGHT_UP, self.on_sys_rightup)
        Syssizer.Add(self.Sys_param[cnt], 1, wx.EXPAND)
      
        self.Sys_panel[cnt].SetSizer(Syssizer)  
        
        bmp = wx.Bitmap(self.PAGE_IMG_SIZE, self.PAGE_IMG_SIZE)
        self.PageImgList.Add(bmp)
        page = self.nb.GetPageCount()-1
        self.set_tabColor(page, self.btn_col[cnt].GetBackgroundColour())
    def add_datapanel(self):
        bmp = wx.Bitmap(self.PAGE_IMG_SIZE, self.PAGE_IMG_SIZE)
        self.PageImgList.Add(bmp)
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
    def on_sys_rightup(self, event):
        pt = event.GetPosition()
        ### event object is property grid. We have that on a panel,
        # that panel is on anotehr panel,
        # which is on a splitwindow
        # which then is on a panel that is the palen this thing knows about
        apr = event.EventObject.GetParent().GetParent().GetParent().GetParent()
        ID = self.get_SysID(apr)
        HitTest = self.Sys_param[ID].pg.HitTest(pt)
        prop = HitTest.GetProperty()
        propname = prop.GetName()
        root_propname = prop.GetParent().GetName()
        if root_propname=='<Root>':
            return
        self.Sys_param[ID].pg.SelectProperty(propname)
        print(propname)
        splt = propname.split('(')
        tmpMenu = wx.Menu()
        tmpMenu.Append(wx.ID_ANY, f'Delete {root_propname}:{splt[0]}')
        tmpMenu.Append(wx.ID_ANY, 'Add Property')
        tmpMenu.Bind(wx.EVT_MENU, lambda e: self.on_sys_mod(e, ID, root_propname, propname))
        pos = wx.GetMousePosition() #apr.GetPosition()
        apr.PopupMenu(tmpMenu, apr.ScreenToClient(pos))
    def on_sys_mod(self, event, ID, root_propname, propname):
        aa = event.GetEventObject()
        ii = event.GetId()
        it = aa.FindItemById(ii)
        if 'del' in it.ItemLabel.lower():
            splt = propname.split('(')
            dlg = wx.MessageDialog(
                self,
                f"Are you sure about deleting {root_propname}:{splt[0]}? Data will be lost",
                "Sure about that?",
                style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE)
            if dlg.ShowModal() == wx.ID_OK:
                Params = self.Sys_param[ID].parameters
                del(Params[root_propname][splt[0]])
                self.parent.Sys[ID].setFromCtrl(Params)
                self.Sys_param[ID].SetFromParClean(Params)

        elif 'add' in it.ItemLabe.lower():
            self.on_addPar(None, ID=ID)

    def on_Opt(self, panel, mainname, name, val):
        self.parent.Opt.setFromCtrl(self.Opt_param.parameters)    
        
    def on_Scale_change(self, event):
        self.parent.ScaleModeChoice = self.cmb_Scale.GetValue()
        self.parent.plot_panel.update_graph()
    def on_delete_nuc(self, event):
        apr = event.EventObject.GetParent()
        ID = self.get_SysID(apr)
        Params = self.Sys_param[ID].parameters
        keys = list(Params.keys())
        nNucs = 0
        items = []
        for key in keys:
            if 'nuc' in key:
                nNucs+=1
        if nNucs == 0:
            return

        if nNucs>1:
            for ii in range(nNucs):
                items.append(f'nuc({ii+1})')

            dlg = wx.SingleChoiceDialog(
                self,
                "Select nucleus to delete \n All data will be lost",
                "Choose Nucleus",
                items,
                style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
            )
        else:
            dlg = wx.MessageDialog(
                self,
                "Are you sure about deleting the nucleus? \n All data will be lost",
                "Sure about that?",
                style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE)

        del_key = ''
        if dlg.ShowModal() == wx.ID_OK:
            if nNucs ==1:
                del_key = 'nuc(1)'
                selected_idx = 1
            else:
                selected_idx = dlg.GetSelection()
                del_key = items[selected_idx]

        dlg.Destroy()

        if del_key!='':
            del(Params[del_key])
            if nNucs>0:
                cnt = 0
                ###  let's resort nuc numbers
                keys = list(Params.keys())
                for kk in keys:
                    if 'nuc' in kk:
                        Params[f'nuc({cnt+1})'] = Params.pop(kk)
                        cnt+=1
            self.parent.Sys[ID].setFromCtrl(Params)
            self.Sys_param[ID].SetFromParClean(Params)
    def on_add_Nuc(self, event):
        apr = event.EventObject.GetParent()
        ID = self.get_SysID(apr)

        items_to_show = []
        isotopes = self.parent.Sys[ID].isotopes()
        isotopes_list = list(self.parent.Sys[ID].isotopes().keys())
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
            
            Params = self.Sys_param[ID].parameters
            nNucs = 0
            keys = list(Params.keys())
            for key in keys:
                if 'nuc' in key:
                    nNucs+=1
            Params[f'nuc({nNucs+1})']=self.getNuc(isotopes_list[selected_idx])
            
            self.parent.Sys[ID].setFromCtrl(Params)
            self.Sys_param[ID].SetFromParClean(Params)
    def getNuc(self, nucname):
        isotopes = list(self.parent.Sys[0].isotopes().keys())
        selected_idx = isotopes.index(nucname)
        nucdict = {'Nucs':[isotopes[selected_idx], isotopes],
                                 'nNucs':1,
                                 'A':np.array([0, 0, 0]),
                                 'Apa':np.array([0, 0, 0])}
        _,I,_=self.parent.Sys[0].isotopes(isotopes[selected_idx])
        if I[0]>0.5:
           nucdict['Q']=np.array([0, 0, 0])
           nucdict['Qpa']=np.array([0, 0, 0])
        #nucdict['useFor'] = ['sim.only', ['sim.only', 'ori.sel.only', 'all']]
        return nucdict
    def on_addPar(self, event, ID=None):
        if ID is None:
            apr = event.EventObject.GetParent()
            ID = self.get_SysID(apr)

        Params = self.Sys_param[ID].parameters
        keys = list(Params.keys())
        nNucs = 0
        nSpins = 0
        for key in keys:
            if 'nuc' in key:
                nNucs+=1
            if 'spin' in key:
                nSpins+=1

        list_params = self.parent.Sys[ID].getAll()

        items_to_show = []

        # order spins first , then nuclei
        for ll in list_params.keys():
            if list_params[ll][0]=='S':
                for ii in range(nSpins):
                    items_to_show.append(f'spin({ii+1}):'+ll)
        for ll in list_params.keys():
            if list_params[ll][0]=='N':
                for ii in range(nNucs):
                    items_to_show.append(f'nuc({ii+1}):'+ll)

        dlg = wx.SingleChoiceDialog(
            self,
            "Select the parameters you want to show:",
            "Choose parameters",
            items_to_show,
            style=wx.OK|wx.CANCEL|wx.DEFAULT_DIALOG_STYLE,
        )
        if dlg.ShowModal() == wx.ID_OK:
            # dlg.GetSelections() returns indices of the checked items
            selected_idx = dlg.GetSelection()
        else:
            selected_idx = -1
        dlg.Destroy()

        if selected_idx>=0:
            splt1 = items_to_show[selected_idx].split(':')
            key_root = splt1[0]
            key_item = splt1[1]
            splt2 =list_params[key_item].split('_')
            if splt2[1]=='float':
                value = float(0.0)
            elif splt2[1]=='int':
                value = int(0)
            elif splt2[1]=='string':
                value = ''
            elif splt2[1]=='array3':
                value = np.zeros(3)

            Params[key_root][key_item] = value

            self.parent.Sys[ID].setFromCtrl(Params)
            self.Sys_param[ID].SetFromParClean(Params)

    def get_SysID(self, panel_handle):
        return self.Sys_panel.index(panel_handle)
    def set_tabColor(self, page, color):
        #bmp = wx.Bitmap(self.PAGE_IMG_SIZE, self.PAGE_IMG_SIZE)
        bmp = self.PageImgList.GetBitmap(page)
        dc = wx.MemoryDC(bmp)
        dc.SetBackground(wx.Brush(color))
        dc.Clear()
        dc.SelectObject(wx.NullBitmap)
        
        self.PageImgList.Replace(page, bmp)
        #self.nb.SetImageList(self.PageImgList)
        self.nb.SetPageImage(page, page)
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
            ### the dumb thing about python is that color we get is just a pointer (no one tells you that, of cause). 
            ### when we delete dialog and try do a couple updates with deleting objects, we also delete 
            ### the memory cell with color ... so we get "Windows fatal exception: access violation" 
            ### when we try access it again ... Idiots. I freaking hate python 
            ### ow, and you cannot do selected_color.copy(), like you do with arrays;   

            selected_color = color_dialog.GetColourData().GetColour()
            self.btn_col[ID].SetBackgroundColour(wx.Colour(selected_color.GetRGBA()) )#SetBackgroundColour
            #self.Sys_panel[ID].SetBackgroundColour(selected_color)
            #self.update_color_display(selected_color)
            #self.status_text.SetLabel(f"Basic picker: {selected_color.GetAsString(wx.C2S_HTML_SYNTAX)}")
            pan = event.EventObject.GetParent()
            page = pan.GetParent().GetSelection()
            self.set_tabColor(page, selected_color)
            if len(self.parent.Data[0]['simdata'])>=(ID-1):
                self.parent.Data[0]['simcolor'][ID] = wx.Colour(selected_color.GetRGBA())
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
        
        self.BMP_LOADEPR = wx.Image('./icons/load_data24_1.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_LOADXML = wx.Image('./icons/load_xml24_2.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_SAVEXML = wx.Image('./icons/save_xml24_1.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_ZOOM = wx.Image('./icons/zoom.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_HOME = wx.Image('./icons/home.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        
        self.BMP_OFF = wx.Image('./icons/off.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()# wx.ArtProvider.GetBitmap(wx.ART_CROSS_MARK, wx.ART_TOOLBAR)
        self.BMP_ON  = wx.Image('./icons/on.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()#wx.ArtProvider.GetBitmap(wx.ART_TICK_MARK, wx.ART_TOOLBAR)
        
        self.BMP_ADDPAR = wx.Image('./icons/plus_param_21_14.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_PLUS = wx.Image('./icons/plus_14_2.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        self.BMP_MINUS = wx.Image('./icons/minus_14_2.png', wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        ###  "Set in stone" properties are in capital. Ideally pull them from ini
        self.MAIN_BACK_COLOR = wx.Colour(255, 255, 255)
        self.MAIN_FOREGROUND_COLOR = wx.Colour(0, 0, 0)
        self.CTRL_BACK_COLOR = wx.Colour(255, 255, 255)
        self.MPL_SUM_COLOR = wx.Colour('red')
        self.MPL_DIFF_COLOR = wx.Colour('green')
        
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
            wildcard="All supported files (*.csv;*.dsc;*.dta)|*.csv;*.dsc;*.dta|Bruker DSC/DTA (*.dsc;*.dta)|*.dsc;*.dta|CSV files (*.csv)|*.csv|All files (*.*)|*.*",
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
                     'simdata':[], 'simax':[], 'simmethod':[], 'simactual':False,
                     'simcolor':[], 'simscale': [],
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
                tData['simcolor']=self.Data[0]['simcolor']
                tData['simscale']=self.Data[0]['simscale']
                tData['show']=self.Data[0]['show']
                self.Data[0] = tData

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

            
    def on_saveSession(self, event):
        if hasattr(wx, 'OSX_FILEDIALOG_ALWAYS_SHOW_TYPES'):
            wx.SystemOptions.SetOption(wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES, 1)
        with wx.FileDialog(self, "Save Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultDir=self.currentPath, 
                           style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT) as fileDialog:
        
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
        
            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            self.dumpXMLSave(pathname)
        
    def dumpXMLSave(self, pathname):
        # otherwise ask the user what new file to open
            try:
                data_dict = {'Data': self.Data,
                             'Opt': self.param_panel.Opt_param.parameters, 
                             'Sys':[],
                             }
                for ii in range(len(self.param_panel.Sys_param)):
                    data_dict['Sys'].append(self.param_panel.Sys_param[ii].parameters.copy()) #### yes, this is a pointer if not copied .. how can you tell, you cannot ... freaking python... I hate it soo much. 
                    col = self.param_panel.btn_col[ii].GetBackgroundColour()
                    data_dict['Sys'][ii]['SysColor'] = col #.GetRGBA()
                    data_dict['Sys'][ii]['SysScale'] = self.param_panel.spn_sca[ii].GetValue()
                    
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
            elif isinstance(value, wx.Colour):
                # Handle primitive types
                primitive_element = ET.SubElement(element, 'wxcolour')
                primitive_element.set('type', type(value).__name__)
                primitive_element.text = str(value.GetRGBA())
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
            
    def on_loadSession(self, event):
        """
        On Mac macOS in the open file dialog the filter choice box is not shown by default. 
        Instead all given wildcards are applied at the same time: 
            To enforce the display of the filter choice set the corresponding 
            wx.SystemOptions before calling the file open dialog:
        """
        if hasattr(wx, 'OSX_FILEDIALOG_ALWAYS_SHOW_TYPES'):
            wx.SystemOptions.SetOption(wx.OSX_FILEDIALOG_ALWAYS_SHOW_TYPES, 1)
        with wx.FileDialog(self, "Load Session", wildcard="Session xml (*.xml)|*.xml|All files (*.*)|*.*",
                           defaultDir=self.currentPath, 
                           style=wx.FD_OPEN) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  
            pathname = fileDialog.GetPath()
            print(pathname)
            self.loadXML(pathname)            
    def loadXML(self, filename):    
        try:
            result = self.load_mixed_dict_xml(filename)
            if 'Data' in result.keys():
                self.Data = result['Data'].copy()
            else:
                self.Data = []
            #if len(self.Data)>0:
                
            self.update_datatree()
                
            if 'Sys' in result.keys():
                if type(result['Sys'])==list:
                    for ii in range(len(result['Sys'])):
                        if ii>=len(self.param_panel.Sys_param):
                            self.param_panel.add_syspanel()
                        self.param_panel.spn_sca[ii].SetValue(result['Sys'][ii]['SysScale'])
                        col = result['Sys'][ii]['SysColor']
                        #col.SetRGBA(result['Sys'][ii]['SysColor'])
                        if type(result['Sys'][ii]['SysColor'])==wx.Colour:
                            self.param_panel.btn_col[ii].SetBackgroundColour(col)
                            page_count = self.param_panel.nb.GetPageCount()
                            for zz in range(page_count):
                                if self.param_panel.nb.GetPage(zz) is self.param_panel.Sys_panel[ii]:
                                    self.param_panel.set_tabColor(zz, col)
                                    break
                            
                        del(result['Sys'][ii]['SysColor'])
                        del(result['Sys'][ii]['SysScale'])
                        self.param_panel.Sys_param[ii].SetFromParClean(result['Sys'][ii].copy())
                        self.Sys[ii].setFromCtrl(self.param_panel.Sys_param[ii].parameters)
                        

                else:
                    self.RaiseError('Wrong Sys in the loaded XML. I expect a list')
                    
            if 'Opt' in result.keys():
                self.param_panel.Opt_param.SetFromParClean(result['Opt'].copy())
                self.Opt.setFromCtrl(self.param_panel.Opt_param.parameters)  
                
            self.plot_panel.update_graph() 
            
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
            wxcolour_element = element.find('wxcolour')
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
                
            elif wxcolour_element is not None:
                # Reconstruct primitive type
                value_type = wxcolour_element.get('type')
                value_text = wxcolour_element.text
    
                if value_type == 'Colour':
                    val = wx.Colour()
                    val.SetRGBA(int(value_text))
                    return key, val
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
if __name__ == "__main__":
    wx.SystemOptions.SetOption("msw.dpiAware", "1")
    app = wx.App(False)

    frame = MainFrame(None)
    frame.Show(True)

    app.MainLoop()