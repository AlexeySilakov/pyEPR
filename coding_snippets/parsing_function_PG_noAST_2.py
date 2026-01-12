import wx
import wx.propgrid as pg
import re


class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(
            None,
            title="Dynamic Function PropertyGrid",
            size=(850, 580)
        )

        panel = wx.Panel(self)
        self.panel = panel

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        # ----------------------------
        # PropertyGrid A(1)...A(10)
        # ----------------------------
        self.pg = pg.PropertyGrid(
            panel,
            style=pg.PG_SPLITTER_AUTO_CENTER
        )

        for i in range(1, 11):
            self.pg.Append(pg.FloatProperty(f"A({i})", value=0.0))

        self.vbox.Add(self.pg, 1, wx.EXPAND | wx.ALL, 5)

        # ----------------------------
        # Function definition row
        # ----------------------------
        func_sizer = wx.BoxSizer(wx.HORIZONTAL)

        func_sizer.Add(
            wx.StaticText(panel, label="Function:"),
            0, wx.ALL | wx.CENTER, 5
        )

        self.func_ctrl = wx.TextCtrl(
            panel,
            value="A(1) = A(3) + 0.5*B",
            size=(430, -1)
        )

        func_sizer.Add(self.func_ctrl, 0, wx.ALL, 5)

        self.btn_define = wx.Button(panel, label="Define Function")
        self.btn_define.Bind(wx.EVT_BUTTON, self.on_define)

        func_sizer.Add(self.btn_define, 0, wx.ALL | wx.CENTER, 5)

        self.vbox.Add(func_sizer, 0, wx.ALL, 5)

        # ----------------------------
        # Variables container
        # ----------------------------
        self.vars_container = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.vars_container, 0, wx.ALL, 5)

        panel.SetSizer(self.vbox)

        # name -> dict(sizer, label, spin)
        self.variables = {}
        self.target_index = None

        self.pg.Bind(pg.EVT_PG_CHANGED, self.recompute)

    # --------------------------------
    # Define / redefine function
    # --------------------------------
    def on_define(self, event):
        expr = self.func_ctrl.GetValue()

        try:
            lhs, rhs = expr.split("=")
            lhs = lhs.strip()
            rhs = rhs.strip()

            if not lhs.startswith("A(") or not lhs.endswith(")"):
                return

            self.target_index = int(lhs[2:-1])

            names = set(re.findall(r"[A-Za-z_]\w*", rhs))

            for i in range(1, 11):
                names.discard(f"A{i}")

            current = set(self.variables.keys())

            to_add = names - current
            to_remove = current - names

            for name in sorted(to_add):
                self.add_variable(name)

            for name in sorted(to_remove):
                self.remove_variable(name)

            self.panel.Layout()
            self.recompute()

        except Exception:
            pass

    # --------------------------------
    # Add variable control (proper sizers)
    # --------------------------------
    def add_variable(self, name):
        row = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self.panel, label=f"{name}:")
        spin = wx.SpinCtrlDouble(
            self.panel,
            min=-1000,
            max=1000,
            inc=0.1,
            initial=0.0
        )

        spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.recompute)

        row.Add(label, 0, wx.ALL | wx.CENTER, 5)
        row.Add(spin, 0, wx.ALL, 5)

        self.vars_container.Add(row, 0, wx.EXPAND)

        self.variables[name] = {
            "sizer": row,
            "label": label,
            "spin": spin
        }

    # --------------------------------
    # Remove variable control cleanly
    # --------------------------------
    def remove_variable(self, name):
        info = self.variables.pop(name)

        sizer = info["sizer"]
        label = info["label"]
        spin = info["spin"]

        self.vars_container.Detach(sizer)

        label.Destroy()
        spin.Destroy()

        sizer.Clear()
        del sizer

    # --------------------------------
    # Recompute expression
    # --------------------------------
    def recompute(self, event=None):
        if self.target_index is None:
            return

        try:
            _, rhs = self.func_ctrl.GetValue().split("=")
            rhs = rhs.strip()

            env = {}

            for i in range(1, 11):
                env[f"A{i}"] = self.pg.GetPropertyValue(f"A({i})")
                rhs = rhs.replace(f"A({i})", f"A{i}")

            for name, info in self.variables.items():
                env[name] = info["spin"].GetValue()

            result = eval(rhs, {"__builtins__": {}}, env)

            self.pg.SetPropertyValue(
                f"A({self.target_index})",
                float(result)
            )

        except Exception:
            pass


class App(wx.App):
    def OnInit(self):
        frame = MainFrame()
        frame.Show()
        return True


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()
