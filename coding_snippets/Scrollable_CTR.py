import wx


# ------------------------------------------------------------
# TripleSpinCtrl: Composite control with 3 SpinCtrlDouble
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# DynamicControlPanel: NOW A ScrolledWindow with scrollbar
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# Demo Frame
# ------------------------------------------------------------
class MyFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Dynamic Panel Demo (Scrollable)", size=(520, 350))

        # Example configuration with many rows to show scrolling
        config = {
            "Position Vector": {
                "type": "TripleSpin",
                "labels": ["X", "Y", "Z"],
                "defaults": [1.0, 2.5, -3.0],
                "min": -100,
                "max": 100,
                "inc": 0.1,
                "digits": 2,
                "units": "mm"
            },
            "Temperature": {
                "type": "SpinCtrlDouble",
                "default": 25.0,
                "min": 0,
                "max": 200,
                "inc": 0.5,
                "digits": 1,
                "units": "°C"
            },
            "Sample Name": {
                "type": "TextCtrl",
                "default": "Experiment_1",
                "units": ""
            },
            "Mode": {
                "type": "ComboBox",
                "choices": ["Standard", "High Speed", "Eco"],
                "default": "High Speed",
                "units": ""
            }
        }

        # Add more test rows to demonstrate scrolling automatically:
        for i in range(20):
            config[f"Extra {i+1}"] = {
                "type": "SpinCtrlDouble",
                "default": float(i),
                "min": 0,
                "max": 100,
                "inc": 1,
                "units": "units"
            }

        panel = DynamicControlPanel(self, config)

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(panel, 1, wx.EXPAND)

        btn = wx.Button(self, label="Print Values")
        btn.Bind(wx.EVT_BUTTON, lambda evt: print(panel.get_values()))
        main_sizer.Add(btn, 0, wx.ALL | wx.CENTER, 10)

        self.SetSizer(main_sizer)
        self.Centre()
        self.Show()


# ------------------------------------------------------------
# Run the App
# ------------------------------------------------------------
if __name__ == "__main__":
    app = wx.App(False)
    MyFrame()
    app.MainLoop()
