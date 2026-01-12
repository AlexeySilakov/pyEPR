
import wx
import wx.stc as stc

# --------------------------------------------------------------------------- #
#  Helper: build a keyword list for a very small C++ subset
# --------------------------------------------------------------------------- #
CPP_KEYWORDS = (
    "int", "float", "double", "char", "void", "if", "else", "while",
    "for", "return", "class", "struct", "enum", "public", "private",
    "protected", "const", "static", "virtual", "operator"
)

# --------------------------------------------------------------------------- #
#  Main Frame
# --------------------------------------------------------------------------- #
class MainFrame(wx.Frame):
    """Main window containing a StyledTextCtrl."""

    def __init__(self, parent, title):
        super().__init__(parent, title=title, size=(800, 600))

        # ----- create the styled text control --------------------------------
        self.stc = stc.StyledTextCtrl(self, -1)
        self._setup_stc()

        # ----- menu ----------------------------------------------------------
        self._create_menu()

        # ----- status bar ----------------------------------------------------
        self.CreateStatusBar()
        self.SetStatusText("Ready")

    # ----------------------------------------------------------------------- #
    #  Set up the StyledTextCtrl (lexer, styles, margin, etc.)
    # ----------------------------------------------------------------------- #
    def _setup_stc(self):
        """Configure lexer, styles and other features of the STC."""
        # Use the built‑in C++ lexer (for demonstration)
        self.stc.SetLexer(stc.STC_LEX_CPP)

        # Keyword list
        self.stc.SetKeyWords(0, " ".join(CPP_KEYWORDS))

        # --------------------------------------------------------------------
        # Styles – each style index has a name in the lexer docs.
        # We are overriding the default styles that the lexer provides.
        # --------------------------------------------------------------------
        # 0 – default
        self.stc.StyleSetSpec(stc.STC_STYLE_DEFAULT, "fore:#000000,back:#FFFFFF")
        self.stc.StyleClearAll()   # reset all styles to the default we just defined

        # 1 – keyword
        self.stc.StyleSetSpec(stc.STC_C_DEFAULT, "fore:#000000")
        self.stc.StyleSetSpec(stc.STC_C_COMMENT, "fore:#007F00,italic")
        self.stc.StyleSetSpec(stc.STC_C_COMMENTLINE, "fore:#007F00,italic")
        self.stc.StyleSetSpec(stc.STC_C_NUMBER, "fore:#007F7F")
        self.stc.StyleSetSpec(stc.STC_C_STRING, "fore:#7F007F")
        self.stc.StyleSetSpec(stc.STC_C_CHARACTER, "fore:#7F007F")
        self.stc.StyleSetSpec(stc.STC_C_UUID, "fore:#7F007F")
        self.stc.StyleSetSpec(stc.STC_C_PREPROCESSOR, "fore:#7F7F00")
        self.stc.StyleSetSpec(stc.STC_C_OPERATOR, "fore:#000000")
        self.stc.StyleSetSpec(stc.STC_C_IDENTIFIER, "fore:#000000")
        self.stc.StyleSetSpec(stc.STC_C_STRINGEOL, "fore:#FF0000,back:#FFFFFF,eolfilled")

        # --------------------------------------------------------------------
        # 2 – margin for line numbers
        # --------------------------------------------------------------------
        self.stc.SetMarginType(0, stc.STC_MARGIN_NUMBER)
        self.stc.SetMarginWidth(0, 40)          # wide enough for 3‑digit numbers

        # --------------------------------------------------------------------
        # 3 – set tab width to 4 spaces
        # --------------------------------------------------------------------
        self.stc.SetTabWidth(4)

        # --------------------------------------------------------------------
        # 4 – add a simple example text
        # --------------------------------------------------------------------
        sample_code = """
// A tiny demo of a C++ program
#include <iostream>

int main() {
    std::cout << "Hello, World!" << std::endl;
    int x = 42;           // the answer
    double y = 3.14159;   // Pi
    return 0;
}
"""
        self.stc.SetText(sample_code)

        # --------------------------------------------------------------------
        # 5 – set a margin for the indicator (e.g., breakpoints)
        # --------------------------------------------------------------------
        self.stc.SetMarginType(1, stc.STC_MARGIN_SYMBOL)
        self.stc.SetMarginMask(1, stc.STC_MASK_FOLDERS)
        self.stc.SetMarginWidth(1, 12)

        # --------------------------------------------------------------------
        # 6 – fold markers (optional)
        # --------------------------------------------------------------------
        self.stc.SetFoldLevel(0, stc.STC_FOLDLEVELBASE)
        self.stc.SetProperty("fold", "1")

    # ----------------------------------------------------------------------- #
    #  Create menu bar
    # ----------------------------------------------------------------------- #
    def _create_menu(self):
        menubar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(wx.ID_OPEN, "&Open\tCtrl+O", "Open a file")
        file_menu.AppendSeparator()
        file_menu.Append(wx.ID_EXIT, "E&xit\tCtrl+Q", "Quit the application")
        menubar.Append(file_menu, "&File")

        self.SetMenuBar(menubar)

        # Bind menu events
        self.Bind(wx.EVT_MENU, self.on_open, id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU, self.on_exit, id=wx.ID_EXIT)

    # ----------------------------------------------------------------------- #
    #  Event handlers
    # ----------------------------------------------------------------------- #
    def on_open(self, event):
        dlg = wx.FileDialog(
            self, "Open File", "", "", "C++ files (*.cpp;*.h)|*.cpp;*.h|All files (*.*)|*.*",
            wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            try:
                with open(path, "r", encoding="utf-8") as f:
                    code = f.read()
                self.stc.SetText(code)
                self.SetTitle(f"StyledTextCtrl Demo – {path}")
                self.SetStatusText(f"Opened {path}")
            except Exception as e:
                wx.MessageBox(f"Could not open file:\n{e}", "Error", wx.OK | wx.ICON_ERROR)
        dlg.Destroy()

    def on_exit(self, event):
        self.Close(True)


# --------------------------------------------------------------------------- #
#  Application entry point
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    app = wx.App(False)          # False → don't redirect stdout/stderr
    frame = MainFrame(None, title="wxPython StyledTextCtrl Demo")
    frame.Show()
    app.MainLoop()