# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:38:14 2026

@author: Alexey

Lets the user parametrize property-grid values (e.g. Sys.A(1)) as Python
expressions of custom scalar parameters (e.g. r, a) and/or other properties
in the same grid (e.g. "spin(1).S"). Expressions are re-evaluated in
dependency order (with cycle detection) whenever a custom parameter changes
or "Recompute All" is pressed; the driven properties become read-only in
the main grid so it's never ambiguous which value is authoritative.
"""
import decimal
import math

import wx
import wx.lib.scrolledpanel
import wx.propgrid as wxpg

import sys_functions as sysfun
import theme
from theme import PillButton


class FunctionModPanel(wx.lib.scrolledpanel.ScrolledPanel):

    # Names available inside expressions in addition to custom parameters
    # and referenced properties. No other builtins are exposed to eval().
    _SAFE_FUNCS = sysfun.SAFE_FUNCS

    @property
    def _ERROR_COLOUR(self):
        return wx.Colour(theme.theme_colors()["error_bg"])

    @property
    def _DRIVEN_COLOUR(self):
        return wx.Colour(theme.theme_colors()["driven_bg"])

    def __init__(self, parent, PropertyGridPanel=None, *args, **kwargs):
        super().__init__(
            parent,
            *args, **kwargs)
        self.pgpanel = PropertyGridPanel
        self.variables = {}
        self.functions = []
        self.function_ctrls = []
        self._driven_labels = set()
        self._resim_pending = False
        self._last_touched_prop_name = None
        # internal parameters
        self.SPACE = 2
        self.F_FUNC_ID = 0
        self.F_SIZE_ID = 1
        self.F_COMP_ID = 2
        self.F_TEXT_ID = 3
        self.F_BDEL_ID = 4
        self.F_BDEF_ID = 5

        self.V_SIZE_ID = 0
        self.V_TEXT_ID = 1
        self.V_SPIN_ID = 2

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        topbar = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_addf = PillButton(self, "add function", color_key="pill_load")
        self.btn_addf.Bind(wx.EVT_BUTTON, self.on_addfunc)
        topbar.Add(self.btn_addf, 0, wx.ALL | wx.CENTER, self.SPACE)

        self.btn_recompute = PillButton(self, "Recompute All", color_key="accent")
        self.btn_recompute.SetToolTip(
            "Re-run every function. Use this after editing a property that a\n"
            "formula references, since that isn't picked up automatically.")
        self.btn_recompute.Bind(wx.EVT_BUTTON, lambda e: self.recompute_all())
        topbar.Add(self.btn_recompute, 0, wx.ALL | wx.CENTER, self.SPACE)
        self.vbox.Add(topbar, 0, wx.ALL | wx.CENTER, self.SPACE)

        # containers
        self.func_container = wx.BoxSizer(wx.VERTICAL)
        self.vars_container = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.func_container, 0, wx.EXPAND, self.SPACE)
        self.vbox.Add(self.vars_container, 0, wx.ALL, self.SPACE)

        self.SetSizer(self.vbox)
        # name -> dict(sizer, label, spin)

        self.target_index = None
        self.SetAutoLayout(1)
        self.SetupScrolling()

    def on_addfunc(self, event=None, function="59/r**3 + a"):
        # Function definition row
        pg_names, pg_labels, pg_parents = self.get_pgNames()
        if not pg_labels:
            wx.MessageBox("No numeric properties are available to target yet.",
                          "Add function", wx.OK | wx.ICON_WARNING)
            return

        func_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_delf = wx.Button(self, label="❌", size=(20, 20))
        btn_delf.Bind(wx.EVT_BUTTON, self.on_delfunc)
        func_sizer.Add(btn_delf, 0, wx.ALL | wx.CENTER, self.SPACE)
        func_choice = wx.Choice(self, choices=pg_labels, size=(140, -1))
        func_choice.SetStringSelection(pg_labels[0])
        func_choice.Bind(wx.EVT_CHOICE, lambda e: self.recompute_all())
        theme.theme_control(func_choice)
        func_sizer.Add(func_choice, 0, wx.ALL | wx.CENTER, self.SPACE)

        txt_f = wx.StaticText(self, label="=")
        func_sizer.Add(txt_f, 0, wx.ALL | wx.CENTER, self.SPACE)

        func_ctrl = wx.TextCtrl(
            self,
            value=function,
            style=wx.TE_PROCESS_ENTER,
        )
        func_ctrl.Bind(wx.EVT_TEXT_ENTER, lambda e: self.recompute_all())
        theme.theme_control(func_ctrl)
        func_sizer.Add(func_ctrl, 1, wx.EXPAND, self.SPACE)

        btn_define = wx.Button(self, label="¶", size=(20, 20))  # pilcrow
        btn_define.SetToolTip("Parse & recompute")
        btn_define.Bind(wx.EVT_BUTTON, self.on_define)

        func_sizer.Add(btn_define, 0, wx.ALL | wx.CENTER, self.SPACE)

        self.func_container.Add(func_sizer, 0, wx.EXPAND, self.SPACE)
        ctrls = [None] * 6
        ctrls[self.F_FUNC_ID] = func_ctrl
        ctrls[self.F_SIZE_ID] = func_sizer
        ctrls[self.F_COMP_ID] = func_choice
        ctrls[self.F_TEXT_ID] = txt_f
        ctrls[self.F_BDEL_ID] = btn_delf
        ctrls[self.F_BDEF_ID] = btn_define
        self.function_ctrls.append(ctrls)
        if event is not None:
            self.Layout()
            self.SetupScrolling()

    def on_delfunc(self, event):
        idx = None
        for i, ctlist in enumerate(self.function_ctrls):
            if event.GetEventObject() in ctlist:
                idx = i
                break
        if idx is None:
            return
        self.del_function(idx)
        self.Layout()
        self.SetupScrolling()
        self.recompute_all()

    def del_function(self, idx):
        ctlist = self.function_ctrls[idx]
        sizer = ctlist[self.F_SIZE_ID]
        self.func_container.Detach(sizer)
        ctlist[self.F_FUNC_ID].Destroy()
        ctlist[self.F_COMP_ID].Destroy()
        ctlist[self.F_TEXT_ID].Destroy()
        ctlist[self.F_BDEL_ID].Destroy()
        ctlist[self.F_BDEF_ID].Destroy()
        sizer.Clear()
        del self.function_ctrls[idx]

    def clear_all(self):
        for idx in range(len(self.function_ctrls) - 1, -1, -1):
            self.del_function(idx)
        self.Layout()
        self.recompute_all()
        self.SetupScrolling()

    def set_functions(self, func_list, clean=False):
        if clean:
            self.clear_all()

        for var, func in func_list:
            self.on_addfunc()

            cnt = len(self.function_ctrls) - 1
            choices = self.function_ctrls[cnt][self.F_COMP_ID].GetItems()
            if var in choices:
                self.function_ctrls[cnt][self.F_COMP_ID].SetStringSelection(var)
            else:
                print(f'[FunctionMod] target "{var}" is not in the current property '
                      f'list ({choices}); leaving row unresolved until it exists again.')
            self.function_ctrls[cnt][self.F_FUNC_ID].SetValue(func)
        self.recompute_all()
        self.SetupScrolling()

    def get_dict(self):
        dct = {}
        dct['func'] = []
        for ctrls in self.function_ctrls:
            var = ctrls[self.F_COMP_ID].GetStringSelection()
            fun = ctrls[self.F_FUNC_ID].GetValue()
            dct['func'].append([var, fun])
        dct['var'] = {}
        for name, info in self.variables.items():
            dct['var'][name] = info[self.V_SPIN_ID].GetValue()
        return dct

    def set_dict(self, dct, clean=True):
        if 'func' not in dct.keys():
            print('[FunctionMod] set_dict: no "func" key in dct - nothing to restore.')
            return
        if clean:
            self.clear_all()
        self.set_functions(dct['func'])
        if 'var' in dct.keys():
            for name, num in dct['var'].items():
                self.add_variable(name)
                self.set_variable_value(
                    self.variables[name][self.V_SPIN_ID], num)
        self.Layout()
        self.recompute_all()
        self.SetupScrolling()

    def on_define(self, event):
        self.recompute_all()

    def get_pgNames(self):
        """Return (names, qualified_labels, parent_labels) for every
        FloatProperty currently in the target grid. Qualified labels are
        "<category>.<leaf>" (e.g. "nuc(1).A(1)") for nested properties, or
        just the leaf label for top-level ones -- this is also the
        canonical form used to reference a property from *another*
        function's expression."""
        name = []
        label = []
        parent = []
        it = self.pgpanel.pg.GetIterator()
        prop = it.GetProperty()

        while prop:
            if isinstance(prop, wxpg.FloatProperty):
                par = prop.GetParent()
                if par is not None and par.GetName() != '<Root>':
                    parent.append(par.GetLabel())
                    label.append(par.GetLabel() + '.' + prop.GetLabel())
                else:
                    parent.append('')
                    label.append(prop.GetLabel())
                name.append(prop.GetName())

            it.Next()
            prop = it.GetProperty()
        return name, label, parent

    def _build_prop_index(self):
        """qualified label -> live wxpg.PGProperty object."""
        names, labels, _parents = self.get_pgNames()
        pg = self.pgpanel.pg
        index = {}
        for nm, lbl in zip(names, labels):
            prop = pg.GetPropertyByName(nm)
            if prop is not None:
                index[lbl] = prop
        return index

    def _mark_row_error(self, ctrls, msg):
        ctrl = ctrls[self.F_FUNC_ID]
        ctrl.SetBackgroundColour(self._ERROR_COLOUR)
        ctrl.SetToolTip(f"Error: {msg}")
        ctrl.Refresh()

    def _mark_row_ok(self, ctrls):
        ctrl = ctrls[self.F_FUNC_ID]
        # Not wx.NullColour: that resolves to the native system default
        # (usually white), which would silently defeat dark mode every
        # time a row recomputes successfully.
        ctrl.SetBackgroundColour(wx.Colour(theme.theme_colors()["ctrl_bg"]))
        ctrl.SetToolTip("")
        ctrl.Refresh()

    def _sync_readonly(self, driven_labels, prop_index):
        driven_labels = set(driven_labels)
        pg = self.pgpanel.pg
        for lbl in self._driven_labels - driven_labels:
            prop = prop_index.get(lbl)
            if prop is not None:
                pg.SetPropertyReadOnly(prop, False)
                # Explicit ctrl_bg, not wx.NullColour -- same "quietly
                # reverts to the native default and defeats dark mode"
                # trap as the function-row textboxes above.
                prop.SetBackgroundColour(wx.Colour(theme.theme_colors()["ctrl_bg"]))
        for lbl in driven_labels:
            prop = prop_index.get(lbl)
            if prop is not None:
                pg.SetPropertyReadOnly(prop, True)
                prop.SetBackgroundColour(self._DRIVEN_COLOUR)
        self._driven_labels = driven_labels
        pg.Refresh()

    def _schedule_resim(self, prop):
        self._last_touched_prop_name = prop.GetName() if prop is not None else None
        if not self._resim_pending:
            self._resim_pending = True
            wx.CallAfter(self._fire_resim)

    def _fire_resim(self):
        self._resim_pending = False
        if self._last_touched_prop_name is None:
            return
        prop = self.pgpanel.pg.GetPropertyByName(self._last_touched_prop_name)
        if prop is not None:
            self.pgpanel.OnValueChanged(None, prop=prop, trigger_run=True)

    # --------------------------------------------------------------
    # Main entry point: (re)parse every function row and push results
    # into the grid, in dependency order, with cycle/error reporting.
    # --------------------------------------------------------------
    def recompute_all(self, event=None):
        prop_index = self._build_prop_index()
        known_labels = list(prop_index.keys())

        rows = [(ctrls[self.F_COMP_ID].GetStringSelection(), ctrls[self.F_FUNC_ID].GetValue())
                for ctrls in self.function_ctrls]

        # Sync custom-parameter controls to the union of free names actually
        # used across all rows (including ones currently flagged as errors,
        # so a variable's control doesn't vanish just because e.g. its
        # sibling row has a duplicate-target problem).
        used_vars = set()
        for lbl, expr in rows:
            _work, _tok, freevars = sysfun.extract_refs_and_freevars(expr, known_labels, self._SAFE_FUNCS)
            used_vars |= freevars
        for nn in sorted(used_vars):
            if nn not in self.variables:
                self.add_variable(nn)
        for nn in list(self.variables.keys()):
            if nn not in used_vars:
                self.remove_variable(nn)

        variable_values = {name: info[self.V_SPIN_ID].GetValue()
                           for name, info in self.variables.items()}

        def lookup(label):
            prop = prop_index.get(label)
            if prop is None:
                raise sysfun.UnresolvedReference(label)
            return prop.GetValue()

        result = sysfun.evaluate_functions(rows, variable_values, lookup,
                                            known_labels=known_labels, safe_funcs=self._SAFE_FUNCS)

        changed_any = False
        last_prop = None
        driven_labels = set()
        for idx, ctrls in enumerate(self.function_ctrls):
            lhs, expr = rows[idx]
            if idx in result['errors']:
                msg = result['errors'][idx]
                self._mark_row_error(ctrls, msg)
                print(f'[FunctionMod] row {idx} ("{lhs} = {expr}"): {msg}')
                continue

            target_prop = prop_index.get(lhs)
            if target_prop is None:
                # Shouldn't happen (lhs came from the choice list built off
                # prop_index), but be defensive.
                continue
            self._mark_row_ok(ctrls)
            value = result['values'][idx]
            driven_labels.add(lhs)
            if target_prop.GetValue() != value:
                changed_any = True
            self.pgpanel.pg.SetPropertyValue(target_prop, value)
            self.pgpanel.OnValueChanged(None, prop=target_prop, trigger_run=False)
            last_prop = target_prop

        self._sync_readonly(driven_labels, prop_index)

        self.pgpanel.parameters['functions'] = self.get_dict()
        if changed_any and last_prop is not None:
            self._schedule_resim(last_prop)

    # --------------------------------------------------------------
    # Add variable control (proper sizers)
    # --------------------------------------------------------------
    def add_variable(self, name):
        if name in self.variables:
            return   # idempotent: don't create a second, orphaned widget row
        row = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, label=f"{name}:", size=(45, -1))
        spin = wx.SpinCtrlDouble(
            self,
            min=-1e12,
            max=1e12,
            inc=0.001,
            initial=0.0,
            style=wx.SP_ARROW_KEYS | wx.TE_PROCESS_ENTER
        )
        spin.SetDigits(6)
        spin.Bind(wx.EVT_CONTEXT_MENU, self.ShowStepMenu)
        spin.Bind(wx.EVT_RIGHT_DOWN, self.ShowStepMenu)
        # Widen the control to whatever precision was typed before reading
        # the value back out, on both ways of committing an edit. Enter goes
        # through theme.bind_spin_enter rather than wx.EVT_TEXT_ENTER, which
        # this control never emitted (it was built without TE_PROCESS_ENTER)
        # and which is unreliable on compound controls even with it.
        theme.bind_spin_enter(
            spin, lambda sp=spin: self.on_variable_edited(None, sp))
        spin.Bind(wx.EVT_KILL_FOCUS,
                  lambda e, sp=spin: self.on_variable_edited(e, sp))
        spin.Bind(wx.EVT_SPINCTRLDOUBLE, lambda e: self.recompute_all())
        theme.theme_control(spin)

        row.Add(label, 0, wx.ALL | wx.CENTER, self.SPACE)
        row.Add(spin, 0, wx.ALL, self.SPACE)

        self.vars_container.Add(row, 0, wx.EXPAND)

        ctrls = [None] * 3
        ctrls[self.V_SIZE_ID] = row
        ctrls[self.V_TEXT_ID] = label
        ctrls[self.V_SPIN_ID] = spin
        self.variables[name] = ctrls
        self.Layout()

    # --------------------------------------------------------------
    # Variable value <-> control, without losing precision
    # --------------------------------------------------------------
    def set_variable_value(self, spin, value):
        """Store `value` in a variable's spin control without losing
        precision.

        SetValue rounds to the control's `digits`, and that defaults to 1 --
        so a fresh control turns 2.8374512 into 2.8 and 0.000123456 into 0.0
        outright. Widen it first: enough decimals for the value, and at
        least enough for the step to have any effect (a 1e-6 nudge on a
        one-decimal control is rounded straight back off).

        The decimals a value needs come from repr(), the shortest text that
        round-trips back to the same float; formatting to a fixed width
        instead would expose the binary representation, since 123456.789 to
        12 decimals is '123456.789000000004'. Past 12 decimals that noise
        shows up regardless, so the count is capped there."""
        try:
            v = float(value)
        except (TypeError, ValueError):
            return

        spin.SetDigits(6)
        spin.SetValue(v)

    def on_variable_edited(self, event, spin):
        """Accept whatever precision was typed, reading the raw text before
        wx rounds it to the control's current digits. `event` is None when
        called from the Enter binding, which does its own skipping."""
        if event is not None:
            event.Skip()      # EVT_KILL_FOCUS must keep propagating
        try:
            self.set_variable_value(spin, float(spin.GetTextValue()))
        except (TypeError, ValueError):
            pass
        self.recompute_all()

    def set_variable_step(self, spin, step):
        """Change the nudge size, re-applying the value so the control keeps
        enough decimals for the new step to move it."""
        spin.SetIncrement(step)
        self.set_variable_value(spin, spin.GetValue())

    def ShowStepMenu(self, evt):
        menu = wx.Menu()
        if type(evt.GetEventObject()) == wx.TextCtrl:
            spin = evt.GetEventObject().Parent
        elif type(evt.GetEventObject()) == wx.SpinCtrlDouble:
            spin = evt.GetEventObject()
        else:
            print(type(evt.GetEventObject()))
            return
        STEP_VALUES = [
            1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1,
            1, 10, 100, 1000, 1e4, 1e5, 1e6
        ]
        for step in STEP_VALUES:
            menu_item = menu.AppendCheckItem(wx.ID_ANY, f"Set step to {step:g}")
            if step == spin.GetIncrement():
                menu_item.Check(True)  # Checks the item
            else:
                menu_item.Check(False)  # Checks the item
            spin.Bind(wx.EVT_MENU,
                      lambda e, s=step, sp=spin: self.set_variable_step(sp, s),
                      menu_item)
        screen_pos = evt.GetPosition()
        local = spin.ScreenToClient(screen_pos)
        spin.PopupMenu(menu, local)
        menu.Destroy()

    def remove_variable(self, name):
        info = self.variables.pop(name)

        sizer = info[self.V_SIZE_ID]
        label = info[self.V_TEXT_ID]
        spin = info[self.V_SPIN_ID]

        self.vars_container.Detach(sizer)

        label.Destroy()
        spin.Destroy()
        sizer.Clear()
        del sizer
        self.Layout()
        self.SetupScrolling()


if __name__ == "__main__":
    import classPropGridPanel as MypgPanel
    import numpy as np

    app = wx.App(False)
    frame = wx.Frame(None, size=(800, 500))
    panel = wx.Panel(frame)
    vbox = wx.BoxSizer(wx.VERTICAL)

    mydic = {'nuc(1)': {'A': np.array([1.0, 2.0, 3.0]),
                    'B': float(123)},
            'nuc(2)': {'A': np.array([1.0, 2.0, 3.0]),
                    'B': float(123)},
            'spin(1)': {'g': np.array([1.0, 2.0, 3.0]),
                    'lw': float(123),
                    'isused': False}
            }

    inpg = MypgPanel.PropGridPanel(panel)
    inpg.SetFromParClean(mydic)

    vbox.Add(inpg, 1, wx.EXPAND | wx.ALL, 5)

    funcmod = FunctionModPanel(panel, PropertyGridPanel=inpg, style=wx.SUNKEN_BORDER, name="panel1",)
    vbox.Add(funcmod, 1, wx.EXPAND | wx.ALL, 5)
    panel.SetSizer(vbox)
    frame.Show()
    app.MainLoop()
