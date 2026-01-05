# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 22:12:26 2025

@author: Alexey Silakov
"""

import wx

# List of common wx.ART_* IDs
art_ids = [
    wx.ART_ERROR, wx.ART_WARNING, wx.ART_QUESTION, wx.ART_INFORMATION,
    wx.ART_MISSING_IMAGE, wx.ART_COPY, wx.ART_CUT, wx.ART_PASTE,
    wx.ART_DELETE, wx.ART_NEW, wx.ART_FILE_OPEN, wx.ART_FILE_SAVE,
    wx.ART_FILE_SAVE_AS, wx.ART_PRINT, wx.ART_HELP, wx.ART_FIND,
    wx.ART_FIND_AND_REPLACE, wx.ART_GO_BACK, wx.ART_GO_FORWARD,
    wx.ART_GO_UP, wx.ART_GO_DOWN, wx.ART_GO_TO_PARENT, wx.ART_GO_HOME,
    wx.ART_EXECUTABLE_FILE, wx.ART_NORMAL_FILE,
    wx.ART_TIP, wx.ART_REPORT_VIEW, wx.ART_LIST_VIEW, wx.ART_NEW_DIR,
    wx.ART_HARDDISK, wx.ART_FLOPPY, wx.ART_CDROM, wx.ART_REMOVABLE,
    wx.ART_FOLDER, wx.ART_FOLDER_OPEN, 
]

class IconFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="wxPython ART Icons 16x16", size=(800, 600))

        panel = wx.Panel(self)
        sizer = wx.GridSizer(rows=0, cols=4, gap=(10, 10))

        for art_id in art_ids:
            bmp = wx.ArtProvider.GetBitmap(art_id, wx.ART_OTHER, (16, 16))
            if bmp.IsOk():
                bmp_ctrl = wx.StaticBitmap(panel, bitmap=bmp)
                label = wx.StaticText(panel, label=str(art_id))
                vbox = wx.BoxSizer(wx.VERTICAL)
                vbox.Add(bmp_ctrl, 0, wx.ALIGN_CENTER)
                vbox.Add(label, 0, wx.ALIGN_CENTER)
                sizer.Add(vbox, 0, wx.ALL, 5)

        panel.SetSizer(sizer)
        self.Layout()

if __name__ == "__main__":
    app = wx.App(False)
    frame = IconFrame()
    frame.Show()
    app.MainLoop()
