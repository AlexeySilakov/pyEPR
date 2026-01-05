import wx
import wx.propgrid as wxpg
import numpy as np  # included per your requirement


class PropertyPanel(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.pg = wxpg.PropertyGrid(
            self,
            style=wxpg.PG_DEFAULT_STYLE |
                  wxpg.PG_SPLITTER_AUTO_CENTER |
                  wxpg.PG_TOOLBAR
        )

        # -------- Main category: wxTextCtrl --------
        cat1 = self.pg.Append(wxpg.PropertyCategory("wxTextCtrl"))

        self.pg.AppendIn(cat1, wxpg.StringProperty("name", value="text"))
        self.pg.AppendIn(cat1, wxpg.StringProperty("style", value=""))
        self.pg.AppendIn(cat1, wxpg.StringProperty("value", value=""))
        self.pg.AppendIn(cat1, wxpg.IntProperty("maxlength", value=0))

        # -------- Second category: wxWindow --------
        cat2 = self.pg.Append(wxpg.PropertyCategory("wxWindow"))

        self.pg.AppendIn(cat2, wxpg.StringProperty("id", value="wxID_ANY"))
        self.pg.AppendIn(cat2, wxpg.StringProperty("pos", value="-1; -1"))

        sizer.Add(self.pg, 1, wx.EXPAND)
        self.SetSizer(sizer)


class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Property Panel Example", size=(360, 350))

        notebook = wx.Notebook(self)

        # Properties tab
        prop_tab = PropertyPanel(notebook)
        notebook.AddPage(prop_tab, "Properties")

        # Events tab (placeholder)
        events_tab = wx.Panel(notebook)
        notebook.AddPage(events_tab, "Events")

        self.Show()


class App(wx.App):
    def OnInit(self):
        frame = MainFrame()
        frame.Show()
        return True


if __name__ == "__main__":
    app = App()
    app.MainLoop()
