import wx
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar


class PlotFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="wxPython + Matplotlib", size=(900, 650))

        self.measure_mode = False
        self.measure_points = []
        self.measure_artists = []

        self._build_menu()
        self._build_toolbar()
        self._build_ui()
        self._bind_events()

        self.Centre()
        self.Show()

    # ------------------------------------------------------------------
    # Menu
    # ------------------------------------------------------------------
    def _build_menu(self):
        menubar = wx.MenuBar()

        view_menu = wx.Menu()
        self.menu_measure = view_menu.Append(wx.ID_ANY, "Measure Distance\tM")

        menubar.Append(view_menu, "&Tools")
        self.SetMenuBar(menubar)

    # ------------------------------------------------------------------
    # Toolbar
    # ------------------------------------------------------------------
    def _build_toolbar(self):
        tb = self.CreateToolBar(style=wx.TB_HORIZONTAL | wx.TB_TEXT)

        self.tool_zoom = tb.AddCheckTool(
            wx.ID_ANY, "Zoom",
            wx.ArtProvider.GetBitmap(wx.ART_FIND, wx.ART_TOOLBAR)
        )
        self.tool_pan = tb.AddTool(
            wx.ID_ANY, "Pan",
            wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_TOOLBAR)
        )
        self.tool_home = tb.AddTool(
            wx.ID_ANY, "Home",
            wx.ArtProvider.GetBitmap(wx.ART_GO_HOME, wx.ART_TOOLBAR)
        )
        self.tool_measure = tb.AddTool(
            wx.ID_ANY, "Measure",
            wx.ArtProvider.GetBitmap(wx.ART_TIP, wx.ART_TOOLBAR)
        )

        tb.Realize()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        self.figure = Figure()
        self.ax = self.figure.add_subplot(111)

        self.canvas = FigureCanvas(panel, -1, self.figure)

        self.mpl_toolbar = NavigationToolbar(self.canvas)
        self.mpl_toolbar.Hide()

        self.coord_text = wx.StaticText(panel, label="x: ---, y: ---")

        vbox.Add(self.canvas, 1, wx.EXPAND)
        vbox.Add(self.coord_text, 0, wx.EXPAND | wx.ALL, 5)

        panel.SetSizer(vbox)

        self._plot_data()

    def _plot_data(self):
        x = np.linspace(0, 10, 500)
        y = np.sin(x)

        self.ax.plot(x, y)
        self.ax.set_title("Sine Wave")

        self.canvas.draw()

    # ------------------------------------------------------------------
    # Events
    # ------------------------------------------------------------------
    def _bind_events(self):
        self.Bind(wx.EVT_TOOL, self.on_zoom, self.tool_zoom)
        self.Bind(wx.EVT_TOOL, self.on_pan, self.tool_pan)
        self.Bind(wx.EVT_TOOL, self.on_home, self.tool_home)
        self.Bind(wx.EVT_TOOL, self.on_measure_toggle, self.tool_measure)

        self.Bind(wx.EVT_MENU, self.on_measure_toggle, self.menu_measure)
        self.Bind(wx.EVT_CHAR_HOOK, self.on_key)

        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect("button_press_event", self.on_mouse_click)

    # ------------------------------------------------------------------
    # Toolbar actions
    # ------------------------------------------------------------------
    def on_zoom(self, event):
        self._disable_measure()
        self.mpl_toolbar.zoom()

    def on_pan(self, event):
        self._disable_measure()
        self.mpl_toolbar.pan()

    def on_home(self, event):
        self.mpl_toolbar.home()

    def on_measure_toggle(self, event):
        self._cancel_navigation_modes()
        self.measure_mode = not self.measure_mode
        self.clear_measurement()
        self.SetStatusText("Measure mode ON" if self.measure_mode else "")

    # ------------------------------------------------------------------
    # Keyboard
    # ------------------------------------------------------------------
    def on_key(self, event):
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self._cancel_navigation_modes()
            self._disable_measure()
        else:
            event.Skip()

    def _cancel_navigation_modes(self):
        if self.mpl_toolbar.mode:
            if "zoom" in self.mpl_toolbar.mode.lower():
                self.mpl_toolbar.zoom()
            elif "pan" in self.mpl_toolbar.mode.lower():
                self.mpl_toolbar.pan()

    def _disable_measure(self):
        self.measure_mode = False
        self.clear_measurement()

    # ------------------------------------------------------------------
    # Mouse
    # ------------------------------------------------------------------
    def on_mouse_move(self, event):
        if event.inaxes:
            self.coord_text.SetLabel(
                f"x: {event.xdata:.3f}, y: {event.ydata:.3f}"
            )
        else:
            self.coord_text.SetLabel("x: ---, y: ---")

    def on_mouse_click(self, event):
        if not self.measure_mode or not event.inaxes or event.button != 1:
            return

        self.measure_points.append((event.xdata, event.ydata))

        if len(self.measure_points) == 2:
            self._draw_measurement()
            self.measure_points.clear()

    # ------------------------------------------------------------------
    # Measurement logic
    # ------------------------------------------------------------------
    def _draw_measurement(self):
        (x1, y1), (x2, y2) = self.measure_points

        dx = x2 - x1
        dy = y2 - y1
        dist = np.hypot(dx, dy)

        # Draw line
        line, = self.ax.plot([x1, x2], [y1, y2], "r--", lw=2)

        # Annotation
        label = (
            f"ΔX = {dx:.3f}\n"
            f"ΔY = {dy:.3f}\n"
            f"D = {dist:.3f}"
        )

        text = self.ax.text(
            (x1 + x2) / 2,
            (y1 + y2) / 2,
            label,
            color="red",
            bbox=dict(facecolor="white", alpha=0.7)
        )

        self.measure_artists.extend([line, text])
        self.canvas.draw()

    def clear_measurement(self):
        for artist in self.measure_artists:
            artist.remove()
        self.measure_artists.clear()
        self.canvas.draw_idle()


class PlotApp(wx.App):
    def OnInit(self):
        self.frame = PlotFrame()
        self.frame.CreateStatusBar()
        return True


if __name__ == "__main__":
    app = PlotApp(False)
    app.MainLoop()
