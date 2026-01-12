# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:12:53 2026

@author: Alexey
"""

import wx
import wx.propgrid as pg


class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="PropertyGrid Function Evaluator",
                         size=(750, 500))

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # ----------------------------
        # PropertyGrid with A(1)..A(10)
        # ----------------------------
        self.pg = pg.PropertyGrid(
            panel,
            style=pg.PG_SPLITTER_AUTO_CENTER
        )

        for i in range(1, 11):
            self.pg.Append(
                pg.FloatProperty(f"A({i})", value=0.0)
            )

        vbox.Add(self.pg, 1, wx.EXPAND | wx.ALL, 5)

        # ----------------------------
        # Button
        # ----------------------------
        self.btn_define = wx.Button(panel, label="Define Function")
        self.btn_define.Bind(wx.EVT_BUTTON, self.on_define)
        vbox.Add(self.btn_define, 0, wx.ALL, 5)

        # Placeholder
        self.func_ctrl = None
        self.spin_B = None

        panel.SetSizer(vbox)

    # --------------------------------
    # Create function controls
    # --------------------------------
    def on_define(self, event):
        if self.func_ctrl:
            return

        panel = self.btn_define.GetParent()
        sizer = panel.GetSizer()

        hbox = wx.BoxSizer(wx.HORIZONTAL)

        hbox.Add(wx.StaticText(panel, label="Function:"),
                 0, wx.ALL | wx.CENTER, 5)

        self.func_ctrl = wx.TextCtrl(
            panel,
            value="A(1) = A(3) + 0.5 * B",
            size=(350, -1)
        )

        hbox.Add(self.func_ctrl, 0, wx.ALL, 5)

        hbox.Add(wx.StaticText(panel, label="B:"),
                 0, wx.ALL | wx.CENTER, 5)

        self.spin_B = wx.SpinCtrlDouble(
            panel,
            min=-100.0,
            max=100.0,
            inc=0.1,
            initial=0.0
        )

        hbox.Add(self.spin_B, 0, wx.ALL, 5)

        sizer.Add(hbox, 0, wx.EXPAND | wx.ALL, 5)
        panel.Layout()

        # Bind updates
        self.func_ctrl.Bind(wx.EVT_TEXT, self.recompute)
        self.spin_B.Bind(wx.EVT_SPINCTRLDOUBLE, self.recompute)
        self.pg.Bind(pg.EVT_PG_CHANGED, self.recompute)

        self.recompute()

    # --------------------------------
    # Core evaluation logic
    # --------------------------------
    def recompute(self, event=None):
        if not self.func_ctrl:
            return

        text = self.func_ctrl.GetValue()

        try:
            lhs, rhs = text.split("=")
            lhs = lhs.strip()
            rhs = rhs.strip()

            # Validate LHS
            if not lhs.startswith("A(") or not lhs.endswith(")"):
                return

            target = int(lhs[2:-1])
            if not (1 <= target <= 10):
                return

            # Build evaluation namespace
            env = {"B": self.spin_B.GetValue()}

            for i in range(1, 11):
                env[f"A{i}"] = self.pg.GetPropertyValue(f"A({i})")

            # Normalize syntax: A(3) -> A3
            for i in range(1, 11):
                rhs = rhs.replace(f"A({i})", f"A{i}")

            # Very restricted eval
            result = eval(
                rhs,
                {"__builtins__": {}},
                env
            )

            self.pg.SetPropertyValue(f"A({target})", float(result))

        except Exception:
            # Deliberately silent; easy to add StatusBar feedback
            pass


class App(wx.App):
    def OnInit(self):
        frame = MainFrame()
        frame.Show()
        return True


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()
