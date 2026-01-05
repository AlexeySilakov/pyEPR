import wx
import wx.lib.agw.customtreectrl as CT

class MyFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Name + Editable Value Tree", size=(600, 400))

        panel = wx.Panel(self)

        self.tree = CT.CustomTreeCtrl(
            panel,
            agwStyle=(
                CT.TR_HAS_BUTTONS |
                CT.TR_HIDE_ROOT |
                CT.TR_DEFAULT_STYLE
            )
        )

        root = self.tree.AddRoot("ROOT")

        # Add several example items
        self.add_item_with_value(root, "Temperature", "298 K")
        self.add_item_with_value(root, "Pressure", "1 atm")
        self.add_item_with_value(root, "pH", "7.4")

        #self.tree.Expand(root)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tree, 1, wx.EXPAND | wx.ALL, 5)
        panel.SetSizer(sizer)

    def add_item_with_value(self, parent, name, value):
        """
        Creates:
            [Name]   [ TextCtrl(Value) ]
        where Name is the tree label, and Value is user-editable.
        """

        # Create label item
        item = self.tree.AppendItem(parent, name)

        # Create the editable value control
        txt = wx.TextCtrl(self.tree, value=value, style=wx.BORDER_NONE)

        # Attach TextCtrl to the tree item
        self.tree.SetItemWindow(item, txt)

        # Optional: handle updates
        txt.Bind(wx.EVT_TEXT, lambda e, it=item: self.on_value_changed(it, e))

        return item

    def on_value_changed(self, item, event):
        value = event.GetEventObject().GetValue()
        name = self.tree.GetItemText(item)
        print(f"{name} updated to: {value}")


class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame()
        frame.Show()
        return True

if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
