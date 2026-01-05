# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 08:10:20 2025

@author: Alexey
"""

import wx

# --- Minimal isotope database (extend as needed) ---
ISOTOPES = {
    "H": [
        {"iso": "¹H", "abundance": 99.9885, "spin": "1/2", "g": 5.5856947},
        {"iso": "²H", "abundance": 0.0115,  "spin": "1",   "g": 0.8574382},
        {"iso": "³H", "abundance": 0.0,     "spin": "1/2", "g": 5.957924},
    ],
    "C": [
        {"iso": "¹²C", "abundance": 98.93, "spin": "0",   "g": 0.0},
        {"iso": "¹³C", "abundance": 1.07,  "spin": "1/2", "g": 1.404825},
    ],
    "N": [
        {"iso": "¹⁴N", "abundance": 99.636, "spin": "1",   "g": 0.403761},
        {"iso": "¹⁵N", "abundance": 0.364,  "spin": "1/2", "g": -0.566378},
    ],
}


# --- Periodic table layout (None = empty cell) ---
PERIODIC_TABLE = [
    ["H",  None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, "He"],
    ["Li", "Be", None, None, None, None, None, None, None, None, None, None, "B",  "C",  "N",  "O",  "F",  "Ne"],
    ["Na", "Mg", None, None, None, None, None, None, None, None, None, None, "Al", "Si", "P",  "S",  "Cl", "Ar"],
    ["K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe"],
    ["Cs", "Ba", "La*", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac*", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    # Spacer row (optional visual separation)
    [None] * 18,
    # Lanthanides
    [None, None, "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", None],
    # Actinides
    [None, None, "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", None],
]

ELEMENT_GROUP = {
    # Alkali metals
    "Li": "alkali", "Na": "alkali", "K": "alkali", "Rb": "alkali",
    "Cs": "alkali", "Fr": "alkali",

    # Alkaline earth metals
    "Be": "alkaline", "Mg": "alkaline", "Ca": "alkaline",
    "Sr": "alkaline", "Ba": "alkaline", "Ra": "alkaline",

    # Transition metals
    "Sc": "transition", "Ti": "transition", "V": "transition",
    "Cr": "transition", "Mn": "transition", "Fe": "transition",
    "Co": "transition", "Ni": "transition", "Cu": "transition",
    "Zn": "transition", "Y": "transition", "Zr": "transition",
    "Nb": "transition", "Mo": "transition", "Tc": "transition",
    "Ru": "transition", "Rh": "transition", "Pd": "transition",
    "Ag": "transition", "Cd": "transition",
    "Hf": "transition", "Ta": "transition", "W": "transition",
    "Re": "transition", "Os": "transition", "Ir": "transition",
    "Pt": "transition", "Au": "transition", "Hg": "transition",
    "Rf": "transition", "Db": "transition", "Sg": "transition",
    "Bh": "transition", "Hs": "transition", "Mt": "transition",
    "Ds": "transition", "Rg": "transition", "Cn": "transition",

    # Post-transition metals
    "Al": "post", "Ga": "post", "In": "post", "Sn": "post",
    "Tl": "post", "Pb": "post", "Bi": "post", "Po": "post",
    "Nh": "post", "Fl": "post", "Mc": "post", "Lv": "post",

    # Metalloids
    "B": "metalloid", "Si": "metalloid", "Ge": "metalloid",
    "As": "metalloid", "Sb": "metalloid", "Te": "metalloid",

    # Nonmetals
    "H": "nonmetal", "C": "nonmetal", "N": "nonmetal",
    "O": "nonmetal", "P": "nonmetal", "S": "nonmetal",
    "Se": "nonmetal",

    # Halogens
    "F": "halogen", "Cl": "halogen", "Br": "halogen",
    "I": "halogen", "At": "halogen", "Ts": "halogen",

    # Noble gases
    "He": "noble", "Ne": "noble", "Ar": "noble",
    "Kr": "noble", "Xe": "noble", "Rn": "noble", "Og": "noble",

    # Lanthanides
    "La": "lanthanide", "Ce": "lanthanide", "Pr": "lanthanide",
    "Nd": "lanthanide", "Pm": "lanthanide", "Sm": "lanthanide",
    "Eu": "lanthanide", "Gd": "lanthanide", "Tb": "lanthanide",
    "Dy": "lanthanide", "Ho": "lanthanide", "Er": "lanthanide",
    "Tm": "lanthanide", "Yb": "lanthanide", "Lu": "lanthanide",

    # Actinides
    "Ac": "actinide", "Th": "actinide", "Pa": "actinide",
    "U": "actinide", "Np": "actinide", "Pu": "actinide",
    "Am": "actinide", "Cm": "actinide", "Bk": "actinide",
    "Cf": "actinide", "Es": "actinide", "Fm": "actinide",
    "Md": "actinide", "No": "actinide", "Lr": "actinide",
}

GROUP_COLORS = {
    "alkali":      wx.Colour(255, 180, 180),
    "alkaline":    wx.Colour(255, 220, 180),
    "transition":  wx.Colour(255, 230, 160),
    "post":        wx.Colour(200, 230, 200),
    "metalloid":   wx.Colour(180, 220, 220),
    "nonmetal":    wx.Colour(200, 220, 255),
    "halogen":     wx.Colour(220, 200, 255),
    "noble":       wx.Colour(230, 200, 240),
    "lanthanide":  wx.Colour(255, 200, 200),
    "actinide":    wx.Colour(255, 180, 220),
}

# ----------------------------------------------------------------------

class PeriodicTableDialog(wx.Dialog):
    def __init__(self, parent):
        super().__init__(parent, title="Select Isotope",
                         size=(900, 600),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)

        self.selected_isotope = None
        self.current_element = None
        self.bntSize = (30, 30)
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # --- Periodic table grid ---
        grid = wx.GridSizer(rows=len(PERIODIC_TABLE), cols=18, hgap=0, vgap=0)

        self.element_buttons = {}
        font = wx.Font(
            10,
            wx.FONTFAMILY_SWISS,
            wx.FONTSTYLE_NORMAL,
            wx.FONTWEIGHT_BOLD
        )

        for row in PERIODIC_TABLE:
            for element in row:
                if element is None:
                    grid.Add(self.bntSize)  # spacer
                else:
                    btn = wx.Button(self, label=element, size=self.bntSize)

                    # Disable elements without isotope data
                    if element not in ISOTOPES or not ISOTOPES[element]:
                        btn.Enable(False)
                    else:
                        btn.SetBackgroundColour(wx.Colour(230, 245, 255))
                    btn.SetFont(font)
                    # Color by periodic table group
                    group = ELEMENT_GROUP.get(element)
                    if group in GROUP_COLORS:
                        btn.SetBackgroundColour(GROUP_COLORS[group])

                    btn.Bind(wx.EVT_BUTTON, self.on_element_clicked)
                    self.element_buttons[btn.GetId()] = element
                    grid.Add(btn, 0, wx.EXPAND)

        main_sizer.Add(grid, 0, wx.ALL | wx.CENTER, 0)

        # --- Isotope list ---
        isotope_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Isotopes")

        self.list_ctrl = wx.ListCtrl(
            self,
            style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.BORDER_SUNKEN
        )
        
        self.list_ctrl.InsertColumn(0, "Isotope", width=90)
        self.list_ctrl.InsertColumn(1, "Abundance (%)", wx.LIST_FORMAT_RIGHT, 120)
        self.list_ctrl.InsertColumn(2, "Spin (I)", wx.LIST_FORMAT_CENTER, 80)
        self.list_ctrl.InsertColumn(3, "gₙ", wx.LIST_FORMAT_RIGHT, 90)

        isotope_box.Add(self.list_ctrl, 1, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(isotope_box, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 10)

        # --- OK / Cancel buttons ---
        btn_sizer = wx.StdDialogButtonSizer()

        ok_btn = wx.Button(self, wx.ID_OK)
        cancel_btn = wx.Button(self, wx.ID_CANCEL)

        ok_btn.Bind(wx.EVT_BUTTON, self.on_ok)

        btn_sizer.AddButton(ok_btn)
        btn_sizer.AddButton(cancel_btn)
        btn_sizer.Realize()

        main_sizer.Add(btn_sizer, 0, wx.ALIGN_RIGHT | wx.ALL, 10)

        self.SetSizer(main_sizer)
        self.Centre()

    # --------------------------------------------------

    def on_element_clicked(self, event):
        element = self.element_buttons[event.GetId()]
        self.current_element = element
    
        self.list_ctrl.DeleteAllItems()
        self.selected_isotope = None
    
        for iso in ISOTOPES.get(element, []):
            idx = self.list_ctrl.InsertItem(
                self.list_ctrl.GetItemCount(),
                iso["iso"]
            )
            self.list_ctrl.SetItem(idx, 1, f'{iso["abundance"]:.4g}')
            self.list_ctrl.SetItem(idx, 2, iso["spin"])
            self.list_ctrl.SetItem(idx, 3, f'{iso["g"]:.6f}')

    def on_isotope_selected(self, event):
        row = event.GetIndex()
        self.selected_isotope = {
            "element": self.current_element,
            "isotope": self.list_ctrl.GetItemText(row, 0),
            "abundance": float(self.list_ctrl.GetItemText(row, 1)),
            "spin": self.list_ctrl.GetItemText(row, 2),
            "g": float(self.list_ctrl.GetItemText(row, 3)),
        }

    def on_ok(self, event):
        if self.selected_isotope is None:
            wx.MessageBox(
                "Please select an isotope.",
                "No selection",
                wx.ICON_WARNING | wx.OK,
                self,
            )
            return

        self.EndModal(wx.ID_OK)

# ----------------------------------------------------------------------

if __name__ == "__main__":
    app = wx.App(False)

    dlg = PeriodicTableDialog(None)
    if dlg.ShowModal() == wx.ID_OK:
        print(dlg.selected_isotope)

    dlg.Destroy()
    app.MainLoop()