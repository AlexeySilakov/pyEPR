# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 10:14:18 2026

@author: Alexey Silakov
"""

import wx
import wx.stc as stc

class StyledTextFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="StyledTextCtrl Example", size=(800, 600))
        
        # Create the StyledTextCtrl
        self.stc = stc.StyledTextCtrl(self, style=wx.TE_MULTILINE | wx.TE_RICH2)
        
        # Set up syntax highlighting for Python
        self.setup_python_syntax()
        
        # Add some sample Python code
        self.set_sample_code()
        
        # Set up sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.stc, 1, wx.EXPAND)
        self.SetSizer(sizer)
        
        # Bind events
        self.Bind(wx.EVT_CLOSE, self.on_close)
        
    def setup_python_syntax(self):
        """Configure syntax highlighting for Python"""
        # Set lexer to Python
        self.stc.SetLexer(stc.STC_LEX_PYTHON)
        
        # Set colors for different elements
        self.stc.StyleSetForeground(stc.STC_P_DEFAULT, wx.Colour(0, 0, 0))      # Default text
        self.stc.StyleSetForeground(stc.STC_P_COMMENTLINE, wx.Colour(0, 128, 0)) # Comment
        self.stc.StyleSetForeground(stc.STC_P_NUMBER, wx.Colour(0, 0, 255))      # Numbers
        self.stc.StyleSetForeground(stc.STC_P_STRING, wx.Colour(255, 0, 0))      # Strings
        self.stc.StyleSetForeground(stc.STC_P_CHARACTER, wx.Colour(255, 0, 0))   # Characters
        self.stc.StyleSetForeground(stc.STC_P_WORD, wx.Colour(0, 0, 255))        # Keywords
        self.stc.StyleSetForeground(stc.STC_P_DEFNAME, wx.Colour(0, 0, 255))     # Function names
        self.stc.StyleSetForeground(stc.STC_P_CLASSNAME, wx.Colour(128, 0, 128)) # Class names
        self.stc.StyleSetForeground(stc.STC_P_OPERATOR, wx.Colour(0, 0, 0))      # Operators
        self.stc.StyleSetForeground(stc.STC_P_IDENTIFIER, wx.Colour(0, 0, 0))    # Identifiers
        
        # Set font
        font = wx.Font(10, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        self.stc.StyleSetFont(stc.STC_P_DEFAULT, font)
        self.stc.StyleSetFont(stc.STC_P_COMMENTLINE, font)
        self.stc.StyleSetFont(stc.STC_P_NUMBER, font)
        self.stc.StyleSetFont(stc.STC_P_STRING, font)
        self.stc.StyleSetFont(stc.STC_P_CHARACTER, font)
        self.stc.StyleSetFont(stc.STC_P_WORD, font)
        self.stc.StyleSetFont(stc.STC_P_DEFNAME, font)
        self.stc.StyleSetFont(stc.STC_P_CLASSNAME, font)
        self.stc.StyleSetFont(stc.STC_P_OPERATOR, font)
        self.stc.StyleSetFont(stc.STC_P_IDENTIFIER, font)
        
        # Enable line numbers
        self.stc.SetMarginWidth(1, 30)
        self.stc.SetMarginType(1, stc.STC_MARGIN_NUMBER)
        self.stc.SetMarginWidth(0, 0)  # Disable first margin
        
        # Enable folding
        self.stc.SetProperty("fold", "1")
        self.stc.SetProperty("fold.compact", "1")
        self.stc.SetProperty("fold.comment", "1")
        self.stc.SetMarginWidth(2, 15)
        self.stc.SetMarginType(2, stc.STC_MARGIN_SYMBOL)
        self.stc.SetMarginMask(2, stc.STC_MASK_FOLDERS)
        self.stc.SetMarginSensitive(2, True)
        self.stc.SetMarginBackground(2, wx.Colour(224, 224, 224))
        
        # Set folding markers
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDEROPEN, stc.STC_MARK_MINUS, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDER, stc.STC_MARK_PLUS, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDERSUB, stc.STC_MARK_EMPTY, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDEREND, stc.STC_MARK_EMPTY, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDEROPENMID, stc.STC_MARK_EMPTY, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDERMIDTAIL, stc.STC_MARK_EMPTY, "white", "black")
        self.stc.MarkerDefine(stc.STC_MARKNUM_FOLDERTAIL, stc.STC_MARK_EMPTY, "white", "black")
        
        # Enable auto-indentation
        #self.stc.SetAutoIndent(True)
        self.stc.SetTabWidth(4)
        self.stc.SetUseTabs(False)
        
    def set_sample_code(self):
        """Add sample Python code with syntax highlighting"""
        sample_code = '''# This is a sample Python program
def fibonacci(n):
    """Calculate Fibonacci sequence"""
    if n <= 1:
        return n
    else:
        return fibonacci(n-1) + fibonacci(n-2)

class Calculator:
    """A simple calculator class"""
    def __init__(self):
        self.result = 0
    
    def add(self, x, y):
        return x + y
    
    def multiply(self, x, y):
        return x * y

# Main execution
if __name__ == "__main__":
    calc = Calculator()
    print("Fibonacci of 10:", fibonacci(10))
    print("Addition:", calc.add(5, 3))
    print("Multiplication:", calc.multiply(4, 6))
'''
        self.stc.SetText(sample_code)
        self.stc.EmptyUndoBuffer()
        
    def on_close(self, event):
        self.Destroy()

class StyledTextApp(wx.App):
    def OnInit(self):
        frame = StyledTextFrame()
        frame.Show()
        return True

if __name__ == "__main__":
    app = StyledTextApp()
    app.MainLoop()