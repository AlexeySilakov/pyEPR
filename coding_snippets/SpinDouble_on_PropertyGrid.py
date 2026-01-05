import wx
import wx.propgrid as pg


# -------------------------------------------------------
#   Custom SpinCtrlDouble Editor for Float Properties
# -------------------------------------------------------
class SpinCtrlDoubleEditor(pg.PGEditor):
    def CreateControls(self, propGrid, property, pos, size):
        self.propGrid = propGrid
        self.property = property
        self.spin = wx.SpinCtrlDouble(propGrid, value=str(property.GetValue()))
        self.spin.SetRange(-1e9, 1e9)
        self.spin.SetIncrement(0.1)
        self.spin.SetDigits(3)
        self.spin.SetPosition(pos)
        self.spin.SetSize(size)
        self.spin.Bind(wx.EVT_CONTEXT_MENU, self.ShowStepMenu)
        self.spin.Bind(wx.EVT_RIGHT_DOWN, self.ShowStepMenu)

        # Following is not needed, but was suggested by ChatGPT. 
        # Just keep in mind if some weird errors occur        
        # Finish editing on Enter
        #spin.Bind(wx.EVT_TEXT_ENTER, self.GetValueFromControl())
        # Finish editing on focus loss
        #spin.Bind(wx.EVT_KILL_FOCUS, self.GetValueFromControl())

        return pg.PGWindowList(self.spin)
    def ShowStepMenu(self, evt):
        menu = wx.Menu()
        STEP_VALUES = [     1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
                            1, 10, 100, 1000, 1e4, 1e5, 1e6
                            ]
        for step in STEP_VALUES:
            item = menu.Append(-1, f"Set step to {step:g}")
            # closure: bind menu selection
            self.spin.Bind(wx.EVT_MENU, lambda e, s=step: self.spin.SetIncrement(s), item)
        screen_pos = evt.GetPosition()
        local_pos = self.spin.ScreenToClient(screen_pos)

        self.spin.PopupMenu(menu, local_pos)
        menu.Destroy()
        
    def UpdateControl(self, property, ctrl):
        ctrl.SetValue(str(property.GetValue()))
    def OnEvent(self, propGrid, property, ctrl, event):
        et = event.GetEventType()
        if et == wx.EVT_SPINCTRLDOUBLE.evtType[0] or et == wx.EVT_TEXT.evtType[0]:
            # return True
            try:
                val = float(ctrl.GetValue())
                
                property.SetValue(val)
                par = propGrid.GetParent()
                par.OnValueChanged(event, property)
            except ValueError:
                pass #return False
        if et == pg.EVT_PG_CHANGED:
            print('aaa')
        return False ###  False is a correct answer to all, somehow. 
    def GetValueFromControl(self, property, ctrl):
        try:
            return float(ctrl.GetValue())
        except ValueError:
            return property.GetValue()

# -------------------------------------------------------
#   Demo Frame
# -------------------------------------------------------
class DemoFrame(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent,size=(200, 400), *args, **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.pg = pg.PropertyGrid(
            self,
            style=pg.PG_SPLITTER_AUTO_CENTER | pg.PG_AUTO_SORT | wx.TAB_TRAVERSAL
        )

        # self.pg.SetExtraStyle(pg.PG_EX_)
        ## this is important to suppress native menue and allow spinctrl to get its own goind
        self.pg.Bind(wx.EVT_CONTEXT_MENU, lambda evt: None) #
        self.pg.Bind(pg.EVT_PG_CHANGED, self.OnValueChanged)
        # Register custom editor
        self.pg.RegisterEditor(SpinCtrlDoubleEditor, "SpinFloat")

        # Add properties
        self.pg.Append(pg.FloatProperty("Float A", value=1.234))
        self.pg.SetPropertyEditor("Float A", "SpinFloat")

        self.pg.Append(pg.FloatProperty("Float B", value=42.0))
        self.pg.SetPropertyEditor("Float B", "SpinFloat")

        # Normal ones to show contrast
        self.pg.Append(pg.StringProperty("String", value="Hello"))
        self.pg.Append(pg.IntProperty("Integer", value=5))

        sizer.Add(self.pg, 1, wx.EXPAND)
        self.SetSizer(sizer)
    def OnValueChanged(self, event, prop=None):
        if type(event)==wx._core.SpinDoubleEvent:
            print(f"{prop.GetName()} = {prop.GetValue()}")
        else:
            prop = event.GetProperty()
            print(f"{prop.GetName()} = {prop.GetValue()}")

# -------------------------------------------------------
#   Run Application
# -------------------------------------------------------
if __name__ == "__main__":
    app = wx.App(False)
    frame = wx.Frame(None)
    DemoFrame(frame)
    
    frame.Show()
    app.MainLoop()
