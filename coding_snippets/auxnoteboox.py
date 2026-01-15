import wx
import wx.lib.agw.aui as aui


class MyFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="AuiNotebook tab reordering", size=(600, 400))

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # --- AuiNotebook with tab moving enabled ---
        self.nb = aui.AuiNotebook(
            panel,
            style=(
                aui.AUI_NB_DEFAULT_STYLE |
                aui.AUI_NB_TAB_MOVE |      # <-- enable drag reordering
                aui.AUI_NB_SCROLL_BUTTONS
            )
        )

        # Add pages
        for i in range(5):
            page = wx.Panel(self.nb)
            s = wx.BoxSizer(wx.VERTICAL)
            s.Add(wx.StaticText(page, label=f"This is page {i}"), 0, wx.ALL, 10)
            page.SetSizer(s)

            self.nb.AddPage(page, f"Tab {i}")

        # Buttons for programmatic reorder
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        btn_left = wx.Button(panel, label="Move Left")
        btn_right = wx.Button(panel, label="Move Right")

        hbox.Add(btn_left, 0, wx.ALL, 5)
        hbox.Add(btn_right, 0, wx.ALL, 5)

        vbox.Add(self.nb, 1, wx.EXPAND | wx.ALL, 5)
        vbox.Add(hbox, 0, wx.ALIGN_CENTER)

        panel.SetSizer(vbox)

        # Bind events
        btn_left.Bind(wx.EVT_BUTTON, self.on_move_left)
        btn_right.Bind(wx.EVT_BUTTON, self.on_move_right)

        self.Centre()
        self.Show()

    # ---------------------------
    # Programmatic tab movement
    # ---------------------------
    def on_move_left(self, event):
        idx = self.nb.GetSelection()
        if idx > 0:
            self.move_page(idx, idx - 1)

    def on_move_right(self, event):
        idx = self.nb.GetSelection()
        if idx < self.nb.GetPageCount() - 1:
            self.move_page(idx, idx + 1)

    def move_page(self, old_idx, new_idx):
        """
        Move a page by removing and reinserting.
        This preserves the window and its contents.
        """
        page = self.nb.GetPage(old_idx)
        text = self.nb.GetPageText(old_idx)
        bmp = self.nb.GetPageBitmap(old_idx)

        self.nb.Freeze()
        self.nb.RemovePage(old_idx)
        self.nb.InsertPage(new_idx, page, text, select=True, bitmap=bmp)
        self.nb.Thaw()


class MyApp(wx.App):
    def OnInit(self):
        MyFrame()
        return True


if __name__ == "__main__":
    app = MyApp()
    app.MainLoop()
