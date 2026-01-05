import wx
# TripleSpinCtrl: Composite control with 3 SpinCtrlDouble
class TripleSpinCtrl(wx.Panel):
    def __init__(self, parent, labels=None, defaults=None, min_val=0, max_val=100,
                 inc=1.0, digits=2):
        super().__init__(parent)

        if defaults is None:
            defaults = [0, 0, 0]

        self.ctrls = []
        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        for i in range(3):
            sub_panel = wx.Panel(self)
            vs = wx.BoxSizer(wx.VERTICAL)

            # Optional sub-labels (X, Y, Z)
            if labels:
                lab = wx.StaticText(sub_panel, label=str(labels[i]))
                vs.Add(lab, 0, wx.BOTTOM | wx.ALIGN_CENTER, 2)

            spin = wx.SpinCtrlDouble(sub_panel)
            spin.SetRange(min_val, max_val)
            spin.SetIncrement(inc)
            spin.SetDigits(digits)
            spin.SetValue(defaults[i])

            vs.Add(spin, 0, wx.EXPAND)
            sub_panel.SetSizer(vs)
            hsizer.Add(sub_panel, 0, wx.RIGHT, 6)
            self.ctrls.append(spin)
        self.SetSizer(hsizer)
    def GetValues(self):
        return [c.GetValue() for c in self.ctrls]
    def BindAll(self, handler):
        for spin in self.spins:
            spin.Bind(wx.EVT_SPINCTRLDOUBLE, handler)
            spin.Bind(wx.EVT_TEXT_ENTER, handler)
            spin.Bind(wx.EVT_KILL_FOCUS, handler)
# DynamicControlPanel: NOW A ScrolledWindow with scrollbar
class DynamicControlPanel(wx.ScrolledWindow):
    def __init__(self, parent, config_dict, label_width=140, units_width=60):
        super().__init__(parent, style=wx.VSCROLL)

        self.SetScrollRate(5, 5)   # Make scrolling smooth
        self.controls = {}

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        grid = wx.FlexGridSizer(rows=len(config_dict), cols=3, vgap=6, hgap=12)
        grid.AddGrowableCol(1, proportion=1)

        for label_text, params in config_dict.items():

            ctrl_type = params.get("type")
            units_text = params.get("units", "")

            label = wx.StaticText(self, label=label_text, size=(label_width, -1))

            # ------------------------------------------------
            # Instantiate appropriate control type
            # ------------------------------------------------
            if ctrl_type == "SpinCtrlDouble":
                ctrl = wx.SpinCtrlDouble(self)
                ctrl.SetRange(params.get("min", 0), params.get("max", 100))
                ctrl.SetIncrement(params.get("inc", 1))
                ctrl.SetDigits(params.get("digits", 2))
                ctrl.SetValue(params.get("default", 0.0))
                

            elif ctrl_type == "TextCtrl":
                ctrl = wx.TextCtrl(self, value=str(params.get("default", "")))

            elif ctrl_type == "ComboBox":
                ctrl = wx.ComboBox(
                    self,
                    choices=params.get("choices", []),
                    style=wx.CB_READONLY
                )
                ctrl.SetValue(str(params.get("default", "")))

            elif ctrl_type == "TripleSpin":
                ctrl = TripleSpinCtrl(
                    self,
                    labels=params.get("labels"),
                    defaults=params.get("defaults", [0, 0, 0]),
                    min_val=params.get("min", 0),
                    max_val=params.get("max", 100),
                    inc=params.get("inc", 1),
                    digits=params.get("digits", 2)
                )

            else:
                raise ValueError(f"Unknown control type: {ctrl_type}")

            units = wx.StaticText(self, label=units_text, size=(units_width, -1))

            # Row layout: Label | Control | Units
            grid.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(ctrl, 1, wx.EXPAND)
            grid.Add(units, 0, wx.ALIGN_CENTER_VERTICAL)

            self.controls[label_text] = ctrl

        main_sizer.Add(grid, 1, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(main_sizer)

        # Recalculate scrollable virtual size
        self.Layout()
        self.FitInside()   # Adjust scrolling area
    def OnChange(self, event):
        print(event)
        
    def BindAllControls(self):
        """Bind change events for all control types."""
        for name, ctrl in self.controls.items():
            if isinstance(ctrl, wx.TextCtrl):
                ctrl.Bind(wx.EVT_TEXT_ENTER, self.OnChange)
                ctrl.Bind(wx.EVT_KILL_FOCUS, self.OnChange)

            elif isinstance(ctrl, wx.SpinCtrlDouble):
                ctrl.Bind(wx.EVT_SPINCTRLDOUBLE, self.OnChange)
                ctrl.Bind(wx.EVT_TEXT_ENTER, self.OnChange)
                ctrl.Bind(wx.EVT_KILL_FOCUS, self.OnChange)

            elif isinstance(ctrl, wx.ComboBox):
                ctrl.Bind(wx.EVT_COMBOBOX, self.OnChange)

            elif isinstance(ctrl, TripleSpin):
                ctrl.BindAll(self.OnChange)
                
    def get_values(self):
        out = {}
        for label, ctrl in self.controls.items():

            if isinstance(ctrl, wx.TextCtrl):
                out[label] = ctrl.GetValue()

            elif isinstance(ctrl, wx.ComboBox):
                out[label] = ctrl.GetValue()

            elif isinstance(ctrl, wx.SpinCtrlDouble):
                out[label] = ctrl.GetValue()

            elif isinstance(ctrl, TripleSpinCtrl):
                out[label] = ctrl.GetValues()

        return out