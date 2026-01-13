# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:38:14 2026

@author: Alexey
"""
import wx
import wx.lib.scrolledpanel
# wx.propgrid as pg
import re


class FunctionModPanel(wx.lib.scrolledpanel.ScrolledPanel):
    def __init__(self, parent, PropertyGridPanel=None, *args, **kwargs):
        super().__init__(
            parent,
            *args, **kwargs)
        self.pgpanel = PropertyGridPanel
        self.variables = {}
        self.functions = []
        self.function_ctrls = []
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
        
        self.btn_addf = wx.Button(self, label="add function", size=(-1, -1))
        self.btn_addf.Bind(wx.EVT_BUTTON, self.on_addfunc)
        self.vbox.Add(self.btn_addf, 0, wx.ALL | wx.CENTER, self.SPACE)
        
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
        #self.pgpanel.Bind(pg.EVT_PG_CHANGED, self.recompute)
    def on_addfunc(self, event=None, function="Aiso + 0.5*Aanis"):
        # Function definition row
        pg_names, pg_labels, pg_parents = self.get_pgNames()
        
        func_sizer = wx.BoxSizer(wx.HORIZONTAL)
        btn_delf = wx.Button(self, label="❌", size=(20, 20))
        btn_delf.Bind(wx.EVT_BUTTON, self.on_delfunc)
        func_sizer.Add(btn_delf, 0, wx.ALL | wx.CENTER, self.SPACE)
        func_choice = wx.Choice(self, choices=pg_labels, size=(90, -1))
        func_choice.SetStringSelection(pg_labels[0])
        func_sizer.Add(func_choice, 0, wx.ALL | wx.CENTER, self.SPACE)
        
        
        txt_f =wx.StaticText(self, label="=") 
        func_sizer.Add(txt_f, 0, wx.ALL | wx.CENTER, self.SPACE)
        
        #func_ctrl = stc.StyledTextCtrl(self, -1)
        
        func_ctrl = wx.TextCtrl(
            self,
            value=function
        )
        
        func_sizer.Add(func_ctrl, 1, wx.EXPAND, self.SPACE)

        btn_define = wx.Button(self, label="¶", size=(20, 20)) #✔▼
        btn_define.Bind(wx.EVT_BUTTON, self.on_define)
        
        func_sizer.Add(btn_define, 0, wx.ALL | wx.CENTER, self.SPACE)

        self.func_container.Add(func_sizer, 0, wx.EXPAND, self.SPACE)
        ctrls = [None]*6
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
        count= 0
        for ctlist in self.function_ctrls:
            if event.GetEventObject() in ctlist:
                self.del_function(count)
                break
            count+=1  
        self.Layout()
        self.SetupScrolling()
        self.recompute()
        
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
        del sizer
        del self.function_ctrls[count]
    def clear_all(self):
        for idx in range(len(self.function_ctrls,0,-1)): ### need to go backwards so that idx is always valid
            self.del_function(idx)
        self.Layout()
        self.recompute()
        self.SetupScrolling()
        
    def set_functions(self, func_list, clean=False):
        if clean:
            self.clear_all() 
                
        for var, func in func_list:
            self.on_addfunc()

            cnt = len(self.function_ctrls)-1
            choices = self.function_ctrls[cnt][self.F_COMP_ID]
            if var in choices:
                self.function_ctrls[cnt][self.F_COMP_ID].SetStringSelection(var)
            else:
                print(f'{var} is not in choices={choices}')
            self.function_ctrls[cnt][self.F_FUNC_ID].SetValue(func)
        self.recompute()
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
            dct['var'][name]=info[self.V_SPIN_ID].GetValue()
        return dct
    def set_dict(self, dct, clean=True):
        if 'func' not in dct.keys():
            print('set_dict: no func in dct. Add some funk')
            return
        if clean:
            self.clear_all()
        self.set_functions(dct['func'])
        if 'var' in dct.keys():
            for name, num in self.dct['var'].items(): 
                self.add_variable(name)
                self.variables[name][V_SPIN_ID].SetValue(float(num))
        self.Layout()
        self.recompute()
        self.SetupScrolling()
    def on_define(self, event):
        varnames=[]
        for ctlist in self.function_ctrls:
            expr = ctlist[self.F_FUNC_ID].GetValue()
            
            try:
                rhs = expr.split("=")
                #lhs = lhs.strip()
                rhs = rhs[-1].strip()
                
                pg_names, pg_labels, pg_parents = self.get_pgNames()
                #if lhs not in pg_labels:
                #    print(f'Variable {lhs} is not a part of the list of property labels')
                #    return
                
                #self.target_index = int(pg_labels.index(lhs))
   
                for ss in pg_labels:
                    if ss in rhs:
                        rhs = re.sub(ss, ' ', rhs)
                varnames += set(re.findall(r"[A-Za-z_]\w*", rhs))
                
                
                for nn in varnames:
                    if nn in pg_names:
                        print('variables from the list are not allowed')
                        return
                    if (nn not in self.variables.keys()) and (nn not in pg_names): 
                        self.add_variable(nn)
    
                for nn in self.variables.keys():
                    if nn not in varnames:
                        self.remove_variable(nn)
    
                self.Layout()
                
                self.recompute()
                self.SetupScrolling()
            except Exception:
                pass
    def get_pgNames(self):
        name = []
        label = []
        parent = []
        it = self.pgpanel.pg.GetIterator()
        prop = it.GetProperty()
        
        while prop:
            if type(prop)==wx._propgrid.FloatProperty:
                par = prop.GetParent()
                if par.GetLabel!='<Root>':
                    parent.append(par.GetLabel())
                    label.append(par.GetLabel()+'.'+prop.GetLabel())
                else:
                    parent.append('')
                    label.append(prop.GetLabel())
                name.append(prop.GetName())
                
            it.Next()
            prop = it.GetProperty()
        return name, label, parent
    # --------------------------------
    # Add variable control (proper sizers)
    # --------------------------------
    def add_variable(self, name):
        row = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, label=f"{name}:", size=(45, -1))
        spin = wx.SpinCtrlDouble(
            self,
            min=-1e12,
            max=1e12,
            inc=0.1,
            initial=0.0
        )
        spin.Bind(wx.EVT_CONTEXT_MENU, self.ShowStepMenu)
        spin.Bind(wx.EVT_RIGHT_DOWN, self.ShowStepMenu)
        #spin.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
        #spin.Bind(wx.EVT_TEXT, self.recompute)
        spin.Bind(wx.EVT_TEXT_ENTER, self.recompute)
        spin.Bind(wx.EVT_SPINCTRLDOUBLE, self.recompute)

        row.Add(label, 0, wx.ALL | wx.CENTER, self.SPACE)
        row.Add(spin, 0, wx.ALL, self.SPACE)

        self.vars_container.Add(row, 0, wx.EXPAND)

        ctrls = [None]*3
        ctrls[self.V_SIZE_ID] = row
        ctrls[self.V_TEXT_ID] = label
        ctrls[self.V_SPIN_ID] = spin
        self.variables[name]=ctrls
        self.Layout()

    def ShowStepMenu(self, evt):
        menu = wx.Menu()
        if type(evt.GetEventObject())==wx.TextCtrl:
            spin = evt.GetEventObject().Parent
        elif type(evt.GetEventObject())==wx.SpinCtrlDouble:
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
            if step==spin.GetIncrement():
                menu_item.Check(True) # Checks the item
            else:
                menu_item.Check(False) # Checks the item
            #item = wx.MenuItem(menu, wx.ID_ANY, f"Set step to {step:g}")
            #menu.Append(item)
            #self.spin.Bind(wx.EVT_MENU, lambda e, s=step: self.spin.SetIncrement(s), item)
            spin.Bind(wx.EVT_MENU, lambda e, s=step: spin.SetIncrement(s), menu_item)
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
    # --------------------------------
    # Recompute expression
    # --------------------------------
    def recompute(self, event=None):
        for ctlist in self.function_ctrls:
            expr = ctlist[self.F_FUNC_ID].GetValue()
            try:
                rhs = expr.split("=")
                lhs = ctlist[self.F_COMP_ID].GetStringSelection()
                rhs = rhs[-1].strip()
                
                varnames = set(re.findall(r"[A-Za-z_]\w*", rhs))
                
                env = {}
                pg_names, pg_labels, pg_parents = self.get_pgNames()
                rhs = rhs.strip()
                
                pg_names, pg_labels, pg_parents = self.get_pgNames()
                if lhs not in pg_labels:
                    print(f'Variable {lhs} is not a part of the list of property labels')
                    return
                
                target_index = int(pg_labels.index(lhs))
                
                #for i in range(1, 11):
                #    env[f"A{i}"] = self.pgpanel.pg.GetPropertyValue(f"A({i})")
                #    rhs = rhs.replace(f"A({i})", f"A{i}")
                
                #corfunc = re.sub(r'\((.*?)\)', r'_\1_', rhs)
                
                for name, info in self.variables.items():
                    corname = re.sub(r'\((.*?)\)', r'_\1_', name)
                    env[corname] = info[self.V_SPIN_ID].GetValue()
                for nn in varnames:
                    if nn in pg_labels:
                        idx= pg_labels.index(nn)
                        env[corname] = self.pgpanel.pg.GetPropertyValue(pg_names[idx])
                
                result = eval(rhs, {"__builtins__": {}}, env)
    
                self.pgpanel.pg.SetPropertyValue(
                    pg_names[target_index],
                    float(result)
                )
                
                prop = self.pgpanel.pg.GetPropertyByName(pg_names[target_index])
                evt = wx.propgrid.PropertyGridEvent(
                    wx.propgrid.wxEVT_PG_CHANGED,
                    self.pgpanel.pg.GetId()
                )
                evt.SetEventObject(self.pgpanel.pg)
                evt.SetProperty(prop)

                self.pgpanel.pg.GetEventHandler().ProcessEvent(evt)
                
            except Exception:
                pass
        self.pgpane.parameters['functions'] = self.get_dict()
if __name__ == "__main__":
    import classPropGridPanel as MypgPanel
    import numpy as np
    
    app = wx.App(False)
    frame = wx.Frame(None, size=(800, 500))
    panel = wx.Panel(frame)
    vbox = wx.BoxSizer(wx.VERTICAL)
    
    mydic =  {'nuc(1)':{'A': np.array([1.0, 2.0, 3.0]),
                    'B': float(123)},
            'nuc(2)': {'A': np.array([1.0, 2.0, 3.0]),
                    'B': float(123)},
            'spin(1)': {'g': np.array([1.0, 2.0, 3.0]),
                    'lw': float(123),
                    'isused':False}
            }
    
    inpg = MypgPanel.PropGridPanel(panel)
    #prop=inpg.GetDefaultDictionary()
    inpg.SetFromParClean(mydic)

    vbox.Add(inpg, 1, wx.EXPAND | wx.ALL, 5)
        
    funcmod = FunctionModPanel(panel, PropertyGridPanel=inpg, style = wx.SUNKEN_BORDER, name="panel1",)
    vbox.Add(funcmod, 1, wx.EXPAND | wx.ALL, 5)
    panel.SetSizer(vbox)
    frame.Show()
    app.MainLoop()