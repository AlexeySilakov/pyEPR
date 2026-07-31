"""
Shared light/dark theme for the pyEPR wx apps -- pill buttons, custom
checkboxes, a folder-style pill tab bar, and a segmented pill choice
control, plus helpers to retint native controls (property grids, list
ctrls, matplotlib figures) that stay native for behavioural reasons.

Palette and pill-widget design follow pyIR_peaksearch_v7.3.py, so the two
apps read as one product family; a handful of new colour keys cover
widget types pyIR doesn't have (property grids, driven/error rows in the
Sys "Function Modifiers" panel, list ctrls, sashes, gauges).

Usage: import this module, build widgets with PillButton/PillCheckBox/
PillTabBar/SegmentedPill instead of the native wx equivalents, and call
apply_theme_to_window(top_level_window, theme) whenever the theme changes
(plus theme_propgrid()/style_figure() for wxpg grids and matplotlib
figures, which apply_theme_to_window doesn't reach into on its own).
"""
import math
import os
import sys

import wx
import wx.propgrid as wxpg

THEMES = {
    "light": {
        "win_bg": "#F4F3F8",
        "ctrl_bg": "#FFFFFF",
        "text": "#4C4A61",
        "text_dim": "#6b6885",
        "border": "#d9d7e6",
        "fig_bg": "#F4F3F8",
        "ax_bg": "#fbfcfe",
        "ax_fg": "#22262d",
        "grid": "#d4d8de",
        "accent": "#535398",
        "pill_load": "#2980B9",
        "pill_calm": "#8f8cc4",
        "pill_warn": "#c0507a",
        "tab_inactive": "#DDE1E6",
        "tab_inactive_icon": "#5b5f68",
        "propgrid_caption_fg": "#FFFFFF",
        "driven_bg": "#E7E6F5",
        "error_bg": "#F8D7DA",
        "error_border": "#c0507a",
        "listctrl_selection": "#535398",
    },
    "dark": {
        "win_bg": "#201E2C",
        "ctrl_bg": "#464D68",
        "text": "#E3E8F7",
        "text_dim": "#AEB4CC",
        "border": "#54597a",
        "fig_bg": "#201E2C",
        "ax_bg": "#262938",
        "ax_fg": "#CCDDFF",
        "grid": "#50545F",
        "accent": "#7A78C9",
        "pill_load": "#4FA6E0",
        "pill_calm": "#66BBC1",
        "pill_warn": "#c97ba0",
        "tab_inactive": "#3a3a46",
        "tab_inactive_icon": "#9aa0b8",
        "propgrid_caption_fg": "#FFFFFF",
        "driven_bg": "#3A3560",
        "error_bg": "#4A2A2E",
        "error_border": "#c97ba0",
        "listctrl_selection": "#7A78C9",
    },
}

# The single live theme name every pill/tab/checkbox reads through --
# updated by set_theme(). A module global (rather than threading a
# callback through every widget) since there's exactly one active theme
# for the whole app at a time.
_CURRENT_THEME = "light"


def theme_colors():
    """The current theme's colour dict. Every pill button/tab/checkbox
    reads its colour from here at paint time, so calling set_theme() and
    refreshing is enough to reskin everything instantly."""
    return THEMES[_CURRENT_THEME]


def set_theme(name):
    global _CURRENT_THEME
    if name not in THEMES:
        raise ValueError(f"Unknown theme {name!r}; expected one of {list(THEMES)}")
    _CURRENT_THEME = name


def get_theme():
    return _CURRENT_THEME


def strip_native_visual_style(win):
    """On Windows, native edit/combo controls (TextCtrl, SpinCtrl,
    Choice, ComboBox...) draw their border/chrome via UxTheme, which
    ignores wx SetBackgroundColour/SetForegroundColour entirely -- that's
    the "white border in dark mode" everyone hits. Telling the control to
    use no visual-style theme falls it back to classic flat GDI rendering,
    which *does* honour our colours. Flat borders also suit this app's
    already-flat pill design better than the native 3D sunken look. No-op
    (and harmless) on non-Windows or if the call fails for any reason."""
    if sys.platform != "win32":
        return
    try:
        import ctypes
        hwnd = win.GetHandle()
        if hwnd:
            ctypes.windll.uxtheme.SetWindowTheme(hwnd, "", "")
    except Exception:
        pass


def _theme_listctrl_header(lc, t):
    """wx.ListCtrl's column-header row is a separate native subcontrol
    that SetBackgroundColour on the list itself doesn't reach -- that's
    the white "Name / Title" header bar staying light in dark mode.
    SetHeaderAttr is the (wx 4.1+) API for it; silently no-ops on older
    wx builds where it doesn't exist."""
    if not hasattr(lc, "SetHeaderAttr"):
        return
    try:
        attr = wx.ItemAttr()
        attr.SetBackgroundColour(wx.Colour(t["ctrl_bg"]))
        attr.SetTextColour(wx.Colour(t["text"]))
        lc.SetHeaderAttr(attr)
    except Exception:
        pass


def theme_control(w):
    """Theme a single freshly-created TextCtrl/SpinCtrl/SpinCtrlDouble/
    Choice/ComboBox in place (colours + native visual-style strip). Use
    this at the point of creation for controls built dynamically *after*
    the initial apply_theme_to_window() pass -- e.g. a new Function-
    Modifier row, a new Grid-Search range box, a property-grid's in-place
    cell editor -- since the recursive walk only re-touches widgets that
    already existed at the time it ran, not ones added afterward."""
    # Strip the native visual style *first* -- switching a control out of
    # themed rendering can reset its background to the classic-mode
    # default (light grey), which would otherwise clobber the colour set
    # applied afterward.
    strip_native_visual_style(w)
    t = theme_colors()
    try:
        w.SetBackgroundColour(wx.Colour(t["ctrl_bg"]))
        w.SetForegroundColour(wx.Colour(t["text"]))
    except Exception:
        pass
    # wx.SpinCtrl/SpinCtrlDouble are composites on MSW: the visible text
    # area is a *separate* child TextCtrl (GetBackgroundColour on the
    # outer control always reports the native default no matter what you
    # set -- it's the child that actually renders). Reach in and theme it
    # too, or the number field silently stays native-grey.
    if isinstance(w, (wx.SpinCtrl, wx.SpinCtrlDouble)):
        for child in w.GetChildren():
            if isinstance(child, wx.TextCtrl):
                strip_native_visual_style(child)
                try:
                    child.SetBackgroundColour(wx.Colour(t["ctrl_bg"]))
                    child.SetForegroundColour(wx.Colour(t["text"]))
                except Exception:
                    pass


# --------------------------------------------------------------------------- #
# Native-control retheming (things that stay native for behavioural reasons:
# property grids, list ctrls, splitters, gauges, status bars, plain
# buttons/checkboxes that were never converted to pill widgets).
# --------------------------------------------------------------------------- #
def apply_theme_to_window(win, theme):
    """Recursively recolour a wx window and its children. Pill widgets
    (PillButton/PillCheckBox/PillTabBar/SegmentedPill) already read
    theme_colors() live at paint time and just need a Refresh(), handled
    by the walk below along with every other widget type."""
    t = THEMES[theme]
    win_bg = wx.Colour(t["win_bg"])
    ctrl_bg = wx.Colour(t["ctrl_bg"])
    fg = wx.Colour(t["text"])
    border = wx.Colour(t["border"])

    def walk(w):
        if isinstance(w, wxpg.PropertyGrid):
            pass   # themed separately via theme_propgrid() -- needs caption-key choice
        elif isinstance(w, wx.ListCtrl):
            # Strip *before* colouring -- switching visual style off can
            # reset the control to its classic-mode default colour, which
            # would otherwise clobber whatever we set below.
            strip_native_visual_style(w)
            try:
                w.SetBackgroundColour(ctrl_bg)
                w.SetForegroundColour(fg)
            except Exception:
                pass
            _theme_listctrl_header(w, t)
        elif isinstance(w, wx.Gauge):
            try:
                w.SetBackgroundColour(win_bg)
                w.SetForegroundColour(wx.Colour(t["accent"]))
            except Exception:
                pass
        elif isinstance(w, (wx.TextCtrl, wx.SpinCtrl, wx.SpinCtrlDouble,
                            wx.Choice, wx.ComboBox)):
            strip_native_visual_style(w)
            w.SetBackgroundColour(ctrl_bg)
            try:
                w.SetForegroundColour(fg)
            except Exception:
                pass
        elif isinstance(w, (wx.Button, wx.ToggleButton, wx.RadioButton)):
            strip_native_visual_style(w)
            w.SetBackgroundColour(win_bg)
            try:
                w.SetForegroundColour(fg)
            except Exception:
                pass
        elif isinstance(w, wx.SplitterWindow):
            try:
                w.SetBackgroundColour(border)
            except Exception:
                pass
            w.SetBackgroundColour(win_bg)
        else:
            # panels, static text, static-box frames, checkboxes, radios...
            w.SetBackgroundColour(win_bg)
            try:
                w.SetForegroundColour(fg)
            except Exception:
                pass
        if isinstance(w, wx.StaticBox):
            try:
                w.SetForegroundColour(fg)
                w.SetBackgroundColour(win_bg)
            except Exception:
                pass
        for child in w.GetChildren():
            walk(child)
        try:
            w.Refresh()
        except Exception:
            pass
    walk(win)

    if isinstance(win, wx.Frame):
        sb = win.GetStatusBar()
        if sb is not None:
            sb.SetBackgroundColour(win_bg)
            sb.SetForegroundColour(fg)
            sb.Refresh()


def theme_propgrid(pg, theme=None, caption_key="accent"):
    """Colour a wx.propgrid.PropertyGrid to match the current theme.
    caption_key selects which palette colour fills category headers --
    pass a second key (e.g. "pill_load") for a second, visually distinct
    grid in the same panel (mirrors the old two-tone File-tree/FFT-params
    scheme, now driven by the theme instead of a hardcoded blue/orange)."""
    t = THEMES[theme or _CURRENT_THEME]
    ctrl_bg = wx.Colour(t["ctrl_bg"])
    text = wx.Colour(t["text"])
    win_bg = wx.Colour(t["win_bg"])
    caption_bg = wx.Colour(t[caption_key])
    caption_fg = wx.Colour(t["propgrid_caption_fg"])
    selection = wx.Colour(t["listctrl_selection"])

    # The control's own base background (blank area past the last row,
    # and anywhere cell colouring doesn't reach) -- SetCellBackgroundColour
    # alone leaves this area at the native/system colour (usually white),
    # which is why it stayed light even after everything else went dark.
    pg.SetBackgroundColour(ctrl_bg)
    pg.SetCellBackgroundColour(ctrl_bg)
    pg.SetCellTextColour(text)
    pg.SetMarginColour(win_bg)
    pg.SetCaptionBackgroundColour(caption_bg)
    pg.SetCaptionTextColour(caption_fg)
    pg.SetSelectionBackgroundColour(selection)
    pg.SetSelectionTextColour(wx.Colour("#FFFFFF"))
    pg.SetLineColour(wx.Colour(t["border"]))
    pg.SetEmptySpaceColour(wx.Colour(win_bg))
    pg.Refresh()


def style_figure(figure, axes_list, theme=None):
    """Recolour a matplotlib figure and its axes for the given theme.
    Unrelated to the *data* colormap (jet/fall/custom) used to render
    spectra -- this only touches figure/axes chrome (background, spines,
    ticks, labels)."""
    t = THEMES[theme or _CURRENT_THEME]
    figure.set_facecolor(t["fig_bg"])
    for ax in axes_list:
        if ax is None:
            continue
        ax.set_facecolor(t["ax_bg"])
        for spine in ax.spines.values():
            spine.set_color(t["ax_fg"])
        ax.tick_params(colors=t["ax_fg"], which="both")
        ax.xaxis.label.set_color(t["ax_fg"])
        ax.yaxis.label.set_color(t["ax_fg"])
        ax.title.set_color(t["ax_fg"])


def recolor_icon_image(path, hex_color):
    """Loads a flat, single-tone icon PNG and repaints it in hex_color,
    keeping the original alpha channel (the glyph shape) intact."""
    img = wx.Image(path, wx.BITMAP_TYPE_PNG)
    if not img.HasAlpha():
        img.InitAlpha()
    alpha = img.GetAlpha()
    w, h = img.GetWidth(), img.GetHeight()
    colour = wx.Colour(hex_color)
    rgb = bytes((colour.Red(), colour.Green(), colour.Blue())) * (w * h)
    img.SetData(rgb)
    img.SetAlpha(alpha)
    return img


def _badge_icon_path(name):
    if not name:
        return None
    base_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base_dir, "assets", "badges", f"{name}.png")


def _load_icon(name):
    path = _badge_icon_path(name)
    if not path:
        return None
    try:
        if os.path.isfile(path):
            return wx.Image(path, wx.BITMAP_TYPE_PNG)
    except Exception:
        pass
    return None


def _shade(color, factor):
    clamp = lambda v: max(0, min(255, int(v * factor)))
    return wx.Colour(clamp(color.Red()), clamp(color.Green()), clamp(color.Blue()))


def _refresh_on_real_size_change(win, event):
    """Shared EVT_SIZE handler for the owner-drawn pill/icon controls below.

    A resize/sash-drag burst fires many EVT_SIZE events per second, and
    several of these controls only ever *move* (their sizer position
    shifts) without their own pixel size actually changing -- e.g. a
    right-aligned toolbar icon when a stretch spacer grows/shrinks. The
    plain `self.Refresh()` these used to call unconditionally repainted
    every one of them on every tick regardless. Comparing against the
    last size actually painted at skips that dead work; controls whose
    size genuinely does change (a pill button stretched by wx.EXPAND)
    still repaint exactly as before.
    """
    new_size = win.GetSize()
    if getattr(win, '_last_paint_size', None) != new_size:
        win._last_paint_size = new_size
        win.Refresh()
    event.Skip()


def _cached_paint_bitmap(win, key):
    """For the owner-drawn controls below: `key` must be a tuple capturing
    every input that determines the control's appearance (size, theme,
    hover/down/checked state, label text, ...). Returns (bitmap, fresh):
    on a cache hit (key unchanged since the control's last paint) `fresh`
    is False and the caller should skip straight to blitting the returned
    bitmap; on a miss (or first paint) it's True and the caller must
    render into the (blank) bitmap before blitting it.

    This exists because these controls redraw via wx.GraphicsContext
    (rounded rects, badges, text layout), which costs real CPU every
    single paint even though most repaints are
    triggered by something that didn't actually change how the control
    looks (a sibling's Refresh(), a stretch-spacer shifting a fixed-size
    icon's position, ...). Re-running that vector drawing only when `key`
    actually changed turns the rest into a cheap bitmap blit.
    """
    w, h = key[0], key[1]
    if w <= 0 or h <= 0:
        return None, False
    cached = getattr(win, '_paint_cache', None)
    if cached is not None and cached[0] == key:
        return cached[1], False
    scale = win.GetDPIScaleFactor()
    bmp = wx.Bitmap(int(w * scale), int(h * scale))
    bmp.SetScaleFactor(scale)
    win._paint_cache = (key, bmp)
    return bmp, True


# --------------------------------------------------------------------------- #
# Pill-style buttons/tabs/checkboxes -- owner-drawn with wx.GraphicsContext
# since native wx.Button/wx.CheckBox/wx.Notebook can't do rounded/pill
# shapes. Every colour is looked up live from THEMES via theme_colors() by
# string key, so switching set_theme()+Refresh() reskins everything.
# --------------------------------------------------------------------------- #
class PillButton(wx.Panel):
    """Flat pill-shaped button: solid colour fill, optional white circular
    icon badge on the left, bold white label. Mimics enough of the
    wx.Button API (Bind(wx.EVT_BUTTON, ...), Enable/Disable, SetLabel) to
    be a drop-in replacement at call sites. icon_name=None renders a
    plain label pill with no badge."""

    def __init__(self, parent, label, color_key="accent", icon_name=None,
                 draw_badge=None, size=None):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._label = label
        self._color_key = color_key
        self._icon = _load_icon(icon_name)
        # draw_badge(gc, bx, by, diameter, color): vector alternative to a
        # PNG icon_name -- draws directly inside the white badge circle at
        # (bx, by, diameter), in `color` (the button's own fill colour, for
        # contrast against the white circle). Takes priority over icon_name.
        self._draw_badge = draw_badge
        self._hover = False
        self._down = False
        if size is None:
            size = wx.Size(-1, self.FromDIP(26))
        self.SetMinSize(size)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_ENTER_WINDOW, self._on_enter)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave)
        self.Bind(wx.EVT_LEFT_DOWN, self._on_down)
        self.Bind(wx.EVT_LEFT_UP, self._on_up)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def SetLabel(self, label):
        self._label = label
        self.Refresh()

    def GetLabel(self):
        return self._label

    def SetColorKey(self, color_key):
        self._color_key = color_key
        self.Refresh()

    def Enable(self, enable=True):
        changed = super().Enable(enable)
        self.Refresh()
        return changed

    def Disable(self):
        return self.Enable(False)

    def _on_enter(self, evt):
        self._hover = True
        self.Refresh()
        evt.Skip()

    def _on_leave(self, evt):
        self._hover = False
        self._down = False
        self.Refresh()
        evt.Skip()

    def _on_down(self, evt):
        if self.IsEnabled():
            self._down = True
            if not self.HasCapture():
                self.CaptureMouse()
            self.Refresh()
        evt.Skip()

    def _on_up(self, evt):
        if self.HasCapture():
            self.ReleaseMouse()
        was_down = self._down
        self._down = False
        self.Refresh()
        if was_down and self.IsEnabled() and \
                self.GetClientRect().Contains(evt.GetPosition()):
            self._fire_clicked()
        evt.Skip()

    def _fire_clicked(self):
        cmd = wx.CommandEvent(wx.wxEVT_BUTTON, self.GetId())
        cmd.SetEventObject(self)
        self.GetEventHandler().ProcessEvent(cmd)

    def _current_color(self):
        base = wx.Colour(theme_colors()[self._color_key])
        if not self.IsEnabled():
            return _shade(base, 0.6)
        if self._down:
            return _shade(base, 0.82)
        if self._hover:
            return _shade(base, 1.12)
        return base

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        color = self._current_color()
        key = (w, h, theme_colors()["win_bg"], color.GetRGBA(), self._label)
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if fresh:
            mdc = wx.MemoryDC(bmp)
            gc = wx.GraphicsContext.Create(mdc)

            if gc is None:
                mdc.SelectObject(wx.NullBitmap)
                return
            mdc.SetBackground(wx.Brush(wx.Colour(theme_colors()["win_bg"])))
            mdc.Clear()

            radius = h / 2.0
            path = gc.CreatePath()
            path.AddRoundedRectangle(0.5, 0.5, max(w - 1, 1), h - 1, radius)
            gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawPath(path)

            has_badge = self._icon is not None or self._draw_badge is not None
            badge_d = h - self.FromDIP(4)
            bx = by = (h - badge_d) / 2.0
            if has_badge:
                gc.SetBrush(gc.CreateBrush(wx.Brush(wx.Colour("#FFFFFF"))))
                gc.DrawEllipse(bx, by, badge_d, badge_d)
                if self._draw_badge is not None:
                    self._draw_badge(gc, bx, by, badge_d, color)
                else:
                    icon_sz = max(1, int(badge_d * 0.6))
                    img = self._icon.Scale(icon_sz, icon_sz, wx.IMAGE_QUALITY_HIGH)
                    ix = bx + (badge_d - icon_sz) / 2.0
                    iy = by + (badge_d - icon_sz) / 2.0
                    gc.DrawBitmap(wx.Bitmap(img), ix, iy, icon_sz, icon_sz)

            font = self.GetFont()
            bold = font.Bold()
            gc.SetFont(bold, wx.Colour("#FFFFFF"))
            left_pad = (bx + badge_d + self.FromDIP(8)) if has_badge else self.FromDIP(12)
            tw, th = gc.GetTextExtent(self._label)[:2]
            avail = w - left_pad - self.FromDIP(8)
            label = self._label
            if avail > 0 and tw > avail:
                while label and gc.GetTextExtent(label + "...")[0] > avail:
                    label = label[:-1]
                if label != self._label:
                    label += "..."
                    tw, th = gc.GetTextExtent(label)[:2]
            text_x = left_pad if has_badge else max(left_pad, (w - tw) / 2.0)
            gc.DrawText(label, text_x, (h - th) / 2.0)
            mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0, useMask=False)

    def DoGetBestSize(self):
        dc = wx.ClientDC(self)
        dc.SetFont(self.GetFont().Bold())
        tw, _th = dc.GetTextExtent(self._label)
        h = self.GetMinSize().height
        if h <= 0:
            h = self.FromDIP(26)
        if self._icon is not None or self._draw_badge is not None:
            badge_d = h - self.FromDIP(4)
            w = badge_d + self.FromDIP(8) + tw + self.FromDIP(12)
        else:
            w = tw + self.FromDIP(24)
        return wx.Size(int(w), int(h))


class PillToggleButton(PillButton):
    """Toggle variant of PillButton: swaps label/colour between a calm OFF
    state and a warning ON state. Mimics the wx.ToggleButton API
    (GetValue/SetValue, fires wx.EVT_TOGGLEBUTTON)."""

    def __init__(self, parent, label_off, label_on,
                 color_key_off="pill_calm", color_key_on="pill_warn", size=None):
        self._label_off, self._label_on = label_off, label_on
        self._color_key_off, self._color_key_on = color_key_off, color_key_on
        self._value = False
        super().__init__(parent, label_off, color_key=color_key_off, size=size)

    def GetValue(self):
        return self._value

    def SetValue(self, value):
        self._value = bool(value)
        self._label = self._label_on if self._value else self._label_off
        self._color_key = self._color_key_on if self._value else self._color_key_off
        self.Refresh()

    def _fire_clicked(self):
        self.SetValue(not self._value)
        cmd = wx.CommandEvent(wx.wxEVT_TOGGLEBUTTON, self.GetId())
        cmd.SetEventObject(self)
        cmd.SetInt(int(self._value))
        self.GetEventHandler().ProcessEvent(cmd)


class PillTabBar(wx.Panel):
    """Folder-style tab strip driving a wx.Simplebook: tabs rounded on the
    top two corners only, filled flat single-tone. Native wx.Notebook
    tabs can't be shaped like this, hence owner-drawing with
    wx.GraphicsContext instead."""

    def __init__(self, parent, book, labels):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._book = book
        self._labels = list(labels)
        self._active = 0
        self._hover = -1
        self._rects = []
        self.SetMinSize(wx.Size(-1, self.FromDIP(32)))
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_MOTION, self._on_motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave)
        self.Bind(wx.EVT_LEFT_DOWN, self._on_click)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def GetSelection(self):
        return self._active

    def SetSelection(self, index, notify=True):
        if index == self._active:
            return
        old = self._active
        self._active = index
        self._book.SetSelection(index)
        self.Refresh()
        if notify:
            self._fire_changed(index, old)

    def _fire_changed(self, index, old):
        evt = wx.BookCtrlEvent(wx.wxEVT_NOTEBOOK_PAGE_CHANGED, self.GetId(), index, old)
        evt.SetEventObject(self)
        self.GetEventHandler().ProcessEvent(evt)

    def _tab_at(self, pos):
        for i, r in enumerate(self._rects):
            if r.Contains(pos):
                return i
        return -1

    def _on_motion(self, evt):
        i = self._tab_at(evt.GetPosition())
        if i != self._hover:
            self._hover = i
            self.Refresh()
        evt.Skip()

    def _on_leave(self, evt):
        if self._hover != -1:
            self._hover = -1
            self.Refresh()
        evt.Skip()

    def _on_click(self, evt):
        i = self._tab_at(evt.GetPosition())
        if i != -1:
            self.SetSelection(i)
        evt.Skip()

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        tc = theme_colors()
        key = (w, h, tc["win_bg"], tc["tab_inactive"], tc["accent"],
              tc["tab_inactive_icon"], tuple(self._labels), self._active, self._hover)
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if not fresh:
            # self._rects (used for click hit-testing) only depends on
            # w/h/labels, all part of `key` -- still valid on a cache hit.
            dc = wx.AutoBufferedPaintDC(self)
            dc.DrawBitmap(bmp, 0, 0)
            return

        mdc = wx.MemoryDC(bmp)
        gc = wx.GraphicsContext.Create(mdc)
        if gc is None:
            mdc.SelectObject(wx.NullBitmap)
            return
        mdc.SetBackground(wx.Brush(wx.Colour(tc["win_bg"])))
        mdc.Clear()

        inactive_fill = wx.Colour(tc["tab_inactive"])
        radius = self.FromDIP(10)
        gap = self.FromDIP(3)
        pad_x = self.FromDIP(14)
        font = self.GetFont().Bold()

        x = 0.0
        self._rects = []
        for i, label in enumerate(self._labels):
            gc.SetFont(font, wx.Colour("#FFFFFF"))
            tw, th = gc.GetTextExtent(label)[:2]
            tab_w = pad_x * 2 + tw
            active = (i == self._active)
            color = wx.Colour(tc["accent"]) if active else inactive_fill
            if i == self._hover and not active:
                color = _shade(color, 1.10)

            path = gc.CreatePath()
            path.MoveToPoint(x, h)
            path.AddLineToPoint(x, radius)
            path.AddArcToPoint(x, 0, x + radius, 0, radius)
            path.AddLineToPoint(x + tab_w - radius, 0)
            path.AddArcToPoint(x + tab_w, 0, x + tab_w, radius, radius)
            path.AddLineToPoint(x + tab_w, h)
            path.CloseSubpath()
            gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawPath(path)

            text_color = (wx.Colour("#FFFFFF") if active
                         else wx.Colour(tc["tab_inactive_icon"]))
            gc.SetFont(font, text_color)
            gc.DrawText(label, x + (tab_w - tw) / 2.0, (h - th) / 2.0)

            self._rects.append(wx.Rect(int(x), 0, int(round(tab_w)), int(h)))
            x += tab_w + gap
        mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0)

    def DoGetBestSize(self):
        return wx.Size(-1, self.FromDIP(32))


class PillCheckBox(wx.Panel):
    """Flat rounded-square checkbox: outlined/empty when unchecked, solid
    theme-accent fill with a white checkmark when checked. Mimics the
    wx.CheckBox API (GetValue/SetValue, fires wx.EVT_CHECKBOX)."""

    def __init__(self, parent, label, color_key="accent"):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._label = label
        self._color_key = color_key
        self._value = False
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_LEFT_UP, self._on_click)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def GetValue(self):
        return self._value

    def SetValue(self, value):
        self._value = bool(value)
        self.Refresh()

    def SetLabel(self, label):
        self._label = label
        self.Refresh()

    def GetLabel(self):
        return self._label

    def Enable(self, enable=True):
        changed = super().Enable(enable)
        self.Refresh()
        return changed

    def Disable(self):
        return self.Enable(False)

    def _on_click(self, evt):
        if self.IsEnabled() and self.GetClientRect().Contains(evt.GetPosition()):
            self._value = not self._value
            self.Refresh()
            cmd = wx.CommandEvent(wx.wxEVT_CHECKBOX, self.GetId())
            cmd.SetEventObject(self)
            cmd.SetInt(int(self._value))
            self.GetEventHandler().ProcessEvent(cmd)
        evt.Skip()

    @staticmethod
    def _blend(c1, c2, t):
        return wx.Colour(int(c1.Red() + (c2.Red() - c1.Red()) * t),
                         int(c1.Green() + (c2.Green() - c1.Green()) * t),
                         int(c1.Blue() + (c2.Blue() - c1.Blue()) * t))

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        tc = theme_colors()
        key = (w, h, tc["win_bg"], tc["text"], tc[self._color_key],
              self._value, self.IsEnabled(), self._label)
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if not fresh:
            dc = wx.AutoBufferedPaintDC(self)
            dc.DrawBitmap(bmp, 0, 0)
            return

        mdc = wx.MemoryDC(bmp)
        gc = wx.GraphicsContext.Create(mdc)
        if gc is None:
            mdc.SelectObject(wx.NullBitmap)
            return
        bg = wx.Colour(tc["win_bg"])
        fg = wx.Colour(tc["text"])
        accent = wx.Colour(tc[self._color_key])
        mdc.SetBackground(wx.Brush(bg))
        mdc.Clear()

        box = self.FromDIP(15)
        bx, by = 0.0, (h - box) / 2.0
        border = self._blend(bg, fg, 0.45 if self.IsEnabled() else 0.22)
        radius = self.FromDIP(4)
        path = gc.CreatePath()
        path.AddRoundedRectangle(bx + 1, by + 1, box - 2, box - 2, radius)

        if self._value:
            color = accent if self.IsEnabled() else self._blend(bg, accent, 0.5)
            gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawPath(path)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#FFFFFF"))
                                   .Width(max(1.4, box * 0.14))
                                   .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
            tick = gc.CreatePath()
            tick.MoveToPoint(bx + box * 0.24, by + box * 0.52)
            tick.AddLineToPoint(bx + box * 0.42, by + box * 0.72)
            tick.AddLineToPoint(bx + box * 0.78, by + box * 0.28)
            gc.StrokePath(tick)
        else:
            gc.SetBrush(wx.TRANSPARENT_BRUSH)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(border).Width(1.4)))
            gc.DrawPath(path)

        font = self.GetFont()
        gc.SetFont(font, fg if self.IsEnabled() else self._blend(bg, fg, 0.4))
        tw, th = gc.GetTextExtent(self._label)[:2]
        gc.DrawText(self._label, bx + box + self.FromDIP(6), (h - th) / 2.0)
        mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0)

    def DoGetBestSize(self):
        dc = wx.ClientDC(self)
        dc.SetFont(self.GetFont())
        tw, th = dc.GetTextExtent(self._label)
        box = self.FromDIP(15)
        return wx.Size(int(box + self.FromDIP(6) + tw + self.FromDIP(2)),
                       int(max(box, th) + self.FromDIP(2)))


class SegmentedPill(wx.Panel):
    """Flat exclusive-choice control: a single pill-shaped bar divided
    into segments, the selected one filled with the theme accent colour.
    Drop-in-ish replacement for a single-row wx.RadioBox -- GetSelection()
    /SetSelection() and fires wx.EVT_RADIOBOX."""

    def __init__(self, parent, choices):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._choices = list(choices)
        self._selection = 0
        self._hover = -1
        self._rects = []
        self.SetMinSize(wx.Size(-1, self.FromDIP(26)))
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_MOTION, self._on_motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave)
        self.Bind(wx.EVT_LEFT_UP, self._on_click)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def GetSelection(self):
        return self._selection

    def SetSelection(self, index):
        if 0 <= index < len(self._choices):
            self._selection = index
            self.Refresh()

    def GetString(self, index):
        return self._choices[index]

    def _fire_changed(self):
        cmd = wx.CommandEvent(wx.wxEVT_COMMAND_RADIOBOX_SELECTED, self.GetId())
        cmd.SetEventObject(self)
        cmd.SetInt(self._selection)
        self.GetEventHandler().ProcessEvent(cmd)

    def _seg_at(self, pos):
        for i, r in enumerate(self._rects):
            if r.Contains(pos):
                return i
        return -1

    def _on_motion(self, evt):
        i = self._seg_at(evt.GetPosition())
        if i != self._hover:
            self._hover = i
            self.Refresh()
        evt.Skip()

    def _on_leave(self, evt):
        if self._hover != -1:
            self._hover = -1
            self.Refresh()
        evt.Skip()

    def _on_click(self, evt):
        i = self._seg_at(evt.GetPosition())
        if i != -1 and i != self._selection:
            self._selection = i
            self.Refresh()
            self._fire_changed()
        evt.Skip()

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        tc = theme_colors()
        key = (w, h, tc["win_bg"], tc["tab_inactive"], tc["accent"],
              tc["tab_inactive_icon"], tuple(self._choices), self._selection, self._hover)
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if not fresh:
            # self._rects (used for click/hover hit-testing) only depends
            # on w/h/choices, all part of `key` -- still valid on a hit.
            dc = wx.AutoBufferedPaintDC(self)
            dc.DrawBitmap(bmp, 0, 0)
            return

        mdc = wx.MemoryDC(bmp)
        gc = wx.GraphicsContext.Create(mdc)
        if gc is None:
            mdc.SelectObject(wx.NullBitmap)
            return
        mdc.SetBackground(wx.Brush(wx.Colour(tc["win_bg"])))
        mdc.Clear()

        inactive_fill = wx.Colour(tc["tab_inactive"])
        radius = h / 2.0
        font = self.GetFont().Bold()

        n = len(self._choices)
        gc.SetFont(font, wx.Colour("#FFFFFF"))
        pad = self.FromDIP(12)
        widths = [gc.GetTextExtent(label)[0] + pad * 2 for label in self._choices]
        total_w = sum(widths)

        path = gc.CreatePath()
        path.AddRoundedRectangle(0.5, 0.5, max(total_w - 1, 1), h - 1, radius)
        gc.SetBrush(gc.CreateBrush(wx.Brush(inactive_fill)))
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.DrawPath(path)

        x = 0.0
        self._rects = []
        for i, label in enumerate(self._choices):
            seg_w = widths[i]
            active = (i == self._selection)
            if active:
                color = wx.Colour(tc["accent"])
                seg_path = gc.CreatePath()
                if n == 1:
                    seg_path.AddRoundedRectangle(0.5, 0.5, max(seg_w - 1, 1), h - 1, radius)
                elif i == 0:
                    seg_path.MoveToPoint(x + seg_w, 0)
                    seg_path.AddLineToPoint(x + radius, 0)
                    seg_path.AddArcToPoint(x, 0, x, radius, radius)
                    seg_path.AddLineToPoint(x, h - radius)
                    seg_path.AddArcToPoint(x, h, x + radius, h, radius)
                    seg_path.AddLineToPoint(x + seg_w, h)
                    seg_path.CloseSubpath()
                elif i == n - 1:
                    seg_path.MoveToPoint(x, 0)
                    seg_path.AddLineToPoint(x + seg_w - radius, 0)
                    seg_path.AddArcToPoint(x + seg_w, 0, x + seg_w, radius, radius)
                    seg_path.AddLineToPoint(x + seg_w, h - radius)
                    seg_path.AddArcToPoint(x + seg_w, h, x + seg_w - radius, h, radius)
                    seg_path.AddLineToPoint(x, h)
                    seg_path.CloseSubpath()
                else:
                    seg_path.MoveToPoint(x, 0)
                    seg_path.AddLineToPoint(x + seg_w, 0)
                    seg_path.AddLineToPoint(x + seg_w, h)
                    seg_path.AddLineToPoint(x, h)
                    seg_path.CloseSubpath()
                gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.DrawPath(seg_path)
                text_color = wx.Colour("#FFFFFF")
            elif i == self._hover:
                text_color = wx.Colour(tc["accent"])
            else:
                text_color = wx.Colour(tc["tab_inactive_icon"])

            gc.SetFont(font, text_color)
            tw, th = gc.GetTextExtent(label)[:2]
            gc.DrawText(label, x + (seg_w - tw) / 2.0, (h - th) / 2.0)

            self._rects.append(wx.Rect(int(x), 0, int(round(seg_w)), int(h)))
            x += seg_w
        mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0)

    def DoGetBestSize(self):
        return wx.Size(-1, self.FromDIP(26))

class ThemedRadioButton(wx.Panel):
    def __init__(self, parent, label, group=None, value=False):
        super().__init__(parent)

        self.label = label
        self.checked = value
        self.group = group if group is not None else []
        self.group.append(self)

        self.callback = None

        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)

        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.Bind(wx.EVT_LEFT_DOWN, self.on_click)

        #Size to fit radio circle and text
        dc = wx.ClientDC(self)
        dc.SetFont(self.GetFont())
        text_w,text_h = dc.GetTextExtent(self.label)
        self.SetMinSize((text_w + 35, max(text_h + 6, 20)))

    def BindRadio(self, callback):
        self.callback = callback

    def on_click(self, event):
        # uncheck group
        for button in self.group:
            button.checked = False

        self.checked = True

        for button in self.group:
            button.Refresh()

        if self.callback:
            self.callback(self)

    def GetValue(self):
        return self.checked

    def SetValue(self, value):
        self.checked = value
        self.Refresh()

    def GetLabel(self):
        return self.label

    def on_paint(self, event):
        dc = wx.AutoBufferedPaintDC(self)

        bg = self.GetBackgroundColour()
        fg = self.GetForegroundColour()

        dc.SetBackground(wx.Brush(bg))
        dc.Clear()

        cy = self.GetSize().height // 2

        # outer circle
        dc.SetPen(wx.Pen(fg, 1))
        dc.SetBrush(wx.Brush(bg))
        dc.DrawCircle(8, cy, 7)

        # filled dot
        if self.checked:
            dc.SetBrush(wx.Brush(fg))
            dc.DrawCircle(8, cy, 4)

        # text
        dc.SetTextForeground(fg)
        dc.DrawText(self.label, 22, cy - 8)

class IconButton(wx.Panel):
    """Small square icon-only toolbar button: a hand-drawn vector glyph
    (no PNG assets needed) on a flat background, with a soft accent-tinted
    circle behind it on hover/press. Mimics wx.Button (Bind(wx.EVT_BUTTON,
    ...), Enable/Disable). draw_icon(gc, w, h, color) is called at paint
    time so the glyph can react live to theme/state changes (e.g. the
    sun/moon swap in draw_theme_icon)."""

    def __init__(self, parent, draw_icon, tooltip="", size=None, fill_key=None):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._draw_icon = draw_icon
        self._hover = False
        self._down = False
        # fill_key: optional palette key (e.g. "accent") for a solid,
        # always-on pill-style fill instead of the default transparent
        # background + soft hover tint -- used for the left-panel collapse
        # toggle so it reads as a clear, permanent accent-coloured control.
        self._fill_key = fill_key
        if size is None:
            size = wx.Size(self.FromDIP(32), self.FromDIP(28))
        self.SetMinSize(size)
        if tooltip:
            self.SetToolTip(tooltip)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_ENTER_WINDOW, self._on_enter)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave)
        self.Bind(wx.EVT_LEFT_DOWN, self._on_down)
        self.Bind(wx.EVT_LEFT_UP, self._on_up)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def SetDrawIcon(self, draw_icon):
        self._draw_icon = draw_icon
        self.Refresh()

    def Enable(self, enable=True):
        changed = super().Enable(enable)
        self.Refresh()
        return changed

    def Disable(self):
        return self.Enable(False)

    def _on_enter(self, evt):
        self._hover = True
        self.Refresh()
        evt.Skip()

    def _on_leave(self, evt):
        self._hover = False
        self._down = False
        self.Refresh()
        evt.Skip()

    def _on_down(self, evt):
        if self.IsEnabled():
            self._down = True
            self.Refresh()
        evt.Skip()

    def _on_up(self, evt):
        was_down = self._down
        self._down = False
        self.Refresh()
        if was_down and self.IsEnabled() and \
                self.GetClientRect().Contains(evt.GetPosition()):
            cmd = wx.CommandEvent(wx.wxEVT_BUTTON, self.GetId())
            cmd.SetEventObject(self)
            self.GetEventHandler().ProcessEvent(cmd)
        evt.Skip()

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        tc = theme_colors()
        # id(self._draw_icon) matters because SetDrawIcon() can swap the
        # glyph callback itself (e.g. the file-browser collapse chevron
        # flips between draw_chevron_icon("left")/("right")) without
        # touching size/theme/hover -- the key must change then too.
        key = (w, h, tc["win_bg"], self._fill_key, tc.get(self._fill_key),
              self.IsEnabled(), self._down, self._hover, tc["accent"], tc["text"],
              id(self._draw_icon))
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if not fresh:
            dc = wx.AutoBufferedPaintDC(self)
            dc.DrawBitmap(bmp, 0, 0)
            return

        mdc = wx.MemoryDC(bmp)
        gc = wx.GraphicsContext.Create(mdc)
        if gc is None:
            mdc.SelectObject(wx.NullBitmap)
            return
        mdc.SetBackground(wx.Brush(wx.Colour(tc["win_bg"])))
        mdc.Clear()

        if self._fill_key:
            base = wx.Colour(tc[self._fill_key])
            if not self.IsEnabled():
                fill = _shade(base, 0.6)
            elif self._down:
                fill = _shade(base, 0.82)
            elif self._hover:
                fill = _shade(base, 1.12)
            else:
                fill = base
            r = min(w, h) * 0.28
            gc.SetBrush(gc.CreateBrush(wx.Brush(fill)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRoundedRectangle(1, 1, w - 2, h - 2, r)
            color = wx.Colour("#FFFFFF")
        else:
            if self.IsEnabled() and (self._hover or self._down):
                accent = wx.Colour(tc["accent"])
                alpha = 180 if self._down else 120
                fill = wx.Colour(accent.Red(), accent.Green(), accent.Blue(), alpha)
                gc.SetBrush(gc.CreateBrush(wx.Brush(fill)))
                #alpha = 210 if self._down else 130
                #gc.SetBrush(gc.CreateBrush(wx.Brush(
                    #wx.Colour(accent.Red(), accent.Green(), accent.Blue(), alpha))))
                gc.SetPen(wx.TRANSPARENT_PEN)
                r = self.FromDIP(8)
                gc.DrawRoundedRectangle(1, 1, w - 2, h - 2, r)
                #d = min(w, h) - self.FromDIP(2)
                #gc.DrawEllipse((w - d) / 2.0, (h - d) / 2.0, d, d)
            color = wx.Colour(tc["text"]) if self.IsEnabled() else _shade(wx.Colour(tc["text"]), 1.6)

        if self._draw_icon:
            self._draw_icon(gc, w, h, color)
        mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0)

    def DoGetBestSize(self):
        return self.GetMinSize()


class IconToggleButton(wx.Panel):
    """Persistent on/off icon toggle for the toolbar (auto-FFT/auto-SIM):
    like IconButton, but keeps state and calls
    draw_icon(gc, w, h, color, is_on) so the glyph itself can change with
    the state (e.g. an infinity symbol when on, an arrow when off).
    Mimics wx.ToggleButton (GetValue/SetValue, fires wx.EVT_TOGGLEBUTTON)."""

    def __init__(self, parent, draw_icon, tooltip="", size=None, value=False):
        super().__init__(parent, style=wx.BORDER_NONE)
        self._draw_icon = draw_icon
        self._value = bool(value)
        self._hover = False
        self._down = False
        if size is None:
            size = wx.Size(self.FromDIP(36), self.FromDIP(36))
        self.SetMinSize(size)
        if tooltip:
            self.SetToolTip(tooltip)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, lambda e: None)
        self.Bind(wx.EVT_ENTER_WINDOW, self._on_enter)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._on_leave)
        self.Bind(wx.EVT_LEFT_DOWN, self._on_down)
        self.Bind(wx.EVT_LEFT_UP, self._on_up)
        self.Bind(wx.EVT_SIZE, lambda e, w=self: _refresh_on_real_size_change(w, e))

    def GetValue(self):
        return self._value

    def SetValue(self, value):
        self._value = bool(value)
        self.Refresh()

    def Enable(self, enable=True):
        changed = super().Enable(enable)
        self.Refresh()
        return changed

    def Disable(self):
        return self.Enable(False)

    def _on_enter(self, evt):
        self._hover = True
        self.Refresh()
        evt.Skip()

    def _on_leave(self, evt):
        self._hover = False
        self._down = False
        self.Refresh()
        evt.Skip()

    def _on_down(self, evt):
        if self.IsEnabled():
            self._down = True
            self.Refresh()
        evt.Skip()

    def _on_up(self, evt):
        was_down = self._down
        self._down = False
        self.Refresh()
        if was_down and self.IsEnabled() and \
                self.GetClientRect().Contains(evt.GetPosition()):
            self._value = not self._value
            self.Refresh()
            cmd = wx.CommandEvent(wx.wxEVT_TOGGLEBUTTON, self.GetId())
            cmd.SetEventObject(self)
            cmd.SetInt(int(self._value))
            self.GetEventHandler().ProcessEvent(cmd)
        evt.Skip()

    def _on_paint(self, _evt):
        w, h = self.GetClientSize()
        tc = theme_colors()
        key = (w, h, tc["win_bg"], tc["accent"], tc["tab_inactive"], tc["tab_inactive_icon"],
              self._value, self.IsEnabled(), self._down, self._hover, id(self._draw_icon))
        bmp, fresh = _cached_paint_bitmap(self, key)
        if bmp is None:
            return
        if not fresh:
            dc = wx.AutoBufferedPaintDC(self)
            dc.DrawBitmap(bmp, 0, 0)
            return

        mdc = wx.MemoryDC(bmp)
        gc = wx.GraphicsContext.Create(mdc)
        if gc is None:
            mdc.SelectObject(wx.NullBitmap)
            return
        mdc.SetBackground(wx.Brush(wx.Colour(tc["win_bg"])))
        mdc.Clear()

        accent = wx.Colour(tc["accent"])
        r = self.FromDIP(8)
        if self._value:
            alpha = 130 if self._down else (100 if self._hover else 72)
            fill = wx.Colour(accent.Red(), accent.Green(), accent.Blue(), alpha)
            gc.SetBrush(gc.CreateBrush(wx.Brush(fill)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRoundedRectangle(1, 1, w - 2, h - 2, r)
        elif self._hover or self._down:
            base = wx.Colour(tc["tab_inactive"])
            alpha = 160 if self._down else 100
            fill = wx.Colour(base.Red(), base.Green(), base.Blue(), alpha)
            gc.SetBrush(gc.CreateBrush(wx.Brush(fill)))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRoundedRectangle(1, 1, w - 2, h - 2, r)

        color = accent if self._value else wx.Colour(tc["tab_inactive_icon"])
        if not self.IsEnabled():
            color = _shade(color, 1.6)
        if self._draw_icon:
            self._draw_icon(gc, w, h, color, self._value)
        mdc.SelectObject(wx.NullBitmap)

        dc = wx.AutoBufferedPaintDC(self)
        dc.DrawBitmap(bmp, 0, 0)

    def DoGetBestSize(self):
        return self.GetMinSize()


def draw_colormap_icon(gc, w, h, color):
    """Small bordered swatch of colour bands -- reads as 'palette/colormap'
    without needing an icon asset file."""
    bw, bh = w * 0.56, h * 0.34
    bx, by = (w - bw) / 2.0, (h - bh) / 2.0
    bands = ["#e0524d", "#e8a33d", "#3fae5c", "#3a76b8"]
    seg = bw / len(bands)
    for i, hexcolor in enumerate(bands):
        gc.SetBrush(gc.CreateBrush(wx.Brush(wx.Colour(hexcolor))))
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.DrawRectangle(bx + i * seg, by, seg + 0.5, bh)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(1.2)))
    gc.DrawRectangle(bx, by, bw, bh)


def draw_theme_icon(gc, w, h, color):
    """Sun/moon glyph representing the theme you'll switch *to* -- moon
    while in light mode (click for dark), sun while in dark mode (click
    for light), matching the common dark-mode-toggle convention."""
    cx, cy = w / 2.0, h / 2.0
    r = min(w, h) * 0.19
    if get_theme() == "dark":
        gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.DrawEllipse(cx - r, cy - r, 2 * r, 2 * r)
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(1.5).Cap(wx.CAP_ROUND)))
        for i in range(8):
            ang = i * math.pi / 4.0
            x1, y1 = cx + math.cos(ang) * r * 1.5, cy + math.sin(ang) * r * 1.5
            x2, y2 = cx + math.cos(ang) * r * 2.15, cy + math.sin(ang) * r * 2.15
            gc.StrokeLine(x1, y1, x2, y2)
    else:
        bg = wx.Colour(theme_colors()["win_bg"])
        gc.SetBrush(gc.CreateBrush(wx.Brush(color)))
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.DrawEllipse(cx - r, cy - r, 2 * r, 2 * r)
        gc.SetBrush(gc.CreateBrush(wx.Brush(bg)))
        gc.DrawEllipse(cx - r + r * 0.65, cy - r - r * 0.3, 2 * r, 2 * r)


def draw_folder_icon(gc, w, h, color):
    """Plain folder outline (open a different data folder)."""
    pts = [(0.08, 0.25), (0.42, 0.25), (0.50, 0.33), (0.875, 0.33),
           (0.875, 0.83), (0.08, 0.83)]
    path = gc.CreatePath()
    path.MoveToPoint(pts[0][0] * w, pts[0][1] * h)
    for px, py in pts[1:]:
        path.AddLineToPoint(px * w, py * h)
    path.CloseSubpath()
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(1.5).Join(wx.JOIN_ROUND)))
    gc.StrokePath(path)


def draw_up_icon(gc, w, h, color):
    """Up chevron over a stem (go to the parent folder)."""
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(2.0)
                           .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
    gc.StrokeLine(0.5 * w, 0.79 * h, 0.5 * w, 0.21 * h)
    path = gc.CreatePath()
    path.MoveToPoint(0.21 * w, 0.5 * h)
    path.AddLineToPoint(0.5 * w, 0.21 * h)
    path.AddLineToPoint(0.79 * w, 0.5 * h)
    gc.StrokePath(path)


def draw_save_icon(gc, w, h, color):
    """Floppy-disk glyph (save the current session to XML)."""
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(1.5).Join(wx.JOIN_ROUND)))
    body = gc.CreatePath()
    body.MoveToPoint(0.21 * w, 0.125 * h)
    body.AddLineToPoint(0.71 * w, 0.125 * h)
    body.AddLineToPoint(0.875 * w, 0.29 * h)
    body.AddLineToPoint(0.875 * w, 0.875 * h)
    body.AddLineToPoint(0.21 * w, 0.875 * h)
    body.CloseSubpath()
    gc.StrokePath(body)
    gc.DrawRectangle(0.33 * w, 0.125 * h, 0.34 * w, 0.25 * h)
    gc.DrawRectangle(0.29 * w, 0.54 * h, 0.42 * w, 0.29 * h)


def draw_gridsearch_icon(gc, w, h, color):
    """3x3 heatmap-like grid with a magnifying glass over it (2D Sys
    parameter grid search / fitness map)."""
    gx, gy, gs = 0.10 * w, 0.10 * h, 0.16 * min(w, h)
    shades = [0.9, 0.55, 0.75, 0.4, 0.2, 0.6, 0.7, 0.45, 0.85]
    for i in range(3):
        for j in range(3):
            a = shades[i * 3 + j]
            gc.SetBrush(gc.CreateBrush(wx.Brush(
                wx.Colour(color.Red(), color.Green(), color.Blue(), int(255 * a)))))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRectangle(gx + i * gs, gy + j * gs, gs * 0.92, gs * 0.92)

    cx, cy, r = 0.62 * w, 0.62 * h, 0.22 * min(w, h)
    gc.SetBrush(wx.TRANSPARENT_BRUSH)
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(1.6).Cap(wx.CAP_ROUND)))
    gc.DrawEllipse(cx - r, cy - r, 2 * r, 2 * r)
    ang = math.pi / 4.0
    x1, y1 = cx + math.cos(ang) * r, cy + math.sin(ang) * r
    x2, y2 = cx + math.cos(ang) * r * 1.7, cy + math.sin(ang) * r * 1.7
    gc.StrokeLine(x1, y1, x2, y2)


def draw_text_icon(label):
    """Returns a draw_icon(gc, w, h, color) callable that just draws
    `label` (e.g. "FFT"/"SIM") bold and centred -- a plain text badge
    instead of a pictorial glyph."""
    def _draw(gc, w, h, color):
        font = wx.Font(max(8, int(h * 0.42)), wx.FONTFAMILY_DEFAULT,
                       wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        gc.SetFont(font, color)
        tw, th = gc.GetTextExtent(label)[:2]
        gc.DrawText(label, (w - tw) / 2.0, (h - th) / 2.0)
    return _draw


def draw_auto_toggle_icon(label):
    """Returns a draw_icon(gc, w, h, color, on) callable for an
    IconToggleButton: `label` ('FFT'/'SIM') drawn large, shifted down,
    overlapping a big infinity symbol when on or a forward arrow when
    off underneath -- the pair reads as "auto-repeat" vs "one-shot"."""
    def _draw(gc, w, h, color, on):
        label_font = wx.Font(max(8, int(h * 0.28)), wx.FONTFAMILY_DEFAULT,
                             wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)
        sym = "∞" if on else "→"
        sym_font = wx.Font(max(10, int(h * 0.36)), wx.FONTFAMILY_DEFAULT,
                           wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD)

        # Measure both pieces first and lay them out from the *actual* text
        # extents (font metrics vary by platform) rather than guessed
        # fractions of h -- that guessing is what let the symbol clip past
        # the bottom of the button before. The two overlap by about half
        # the label's height, and the whole overlapping group is centred
        # in the button with a slight downward nudge.
        gc.SetFont(label_font, color)
        tw, th = gc.GetTextExtent(label)[:2]
        gc.SetFont(sym_font, color)
        sw, sh = gc.GetTextExtent(sym)[:2]

        overlap = th * 0.5
        total_h = th + sh - overlap
        start_y = max(0.0, (h - total_h) / 2.0 + h * 0.06)

        gc.SetFont(label_font, color)
        gc.DrawText(label, (w - tw) / 2.0, start_y)
        gc.SetFont(sym_font, color)
        gc.DrawText(sym, (w - sw) / 2.0, start_y + th - overlap)
    return _draw


def draw_chevron_icon(direction):
    """Returns a draw_icon(gc, w, h, color) callable drawing a simple
    left/right chevron -- used for the file-browser collapse toggle."""
    def _draw(gc, w, h, color):
        cx, cy = w / 2.0, h / 2.0
        size = min(w, h) * 0.22
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(2.0)
                               .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
        path = gc.CreatePath()
        if direction == "left":
            path.MoveToPoint(cx + size * 0.5, cy - size)
            path.AddLineToPoint(cx - size * 0.5, cy)
            path.AddLineToPoint(cx + size * 0.5, cy + size)
        else:
            path.MoveToPoint(cx - size * 0.5, cy - size)
            path.AddLineToPoint(cx + size * 0.5, cy)
            path.AddLineToPoint(cx - size * 0.5, cy + size)
        gc.StrokePath(path)
    return _draw


def draw_squiggle_badge(gc, bx, by, d, color):
    """Wavy-line badge for the 'Load' pill button (a loaded waveform)."""
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(max(1.2, d * 0.09))
                           .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
    path = gc.CreatePath()
    n = 16
    amp = d * 0.16
    y0 = by + d * 0.5
    x0 = bx + d * 0.22
    xw = d * 0.56
    path.MoveToPoint(x0, y0)
    for i in range(1, n + 1):
        t = i / n
        x = x0 + xw * t
        y = y0 - amp * math.sin(t * 2 * math.pi * 1.5)
        path.AddLineToPoint(x, y)
    gc.StrokePath(path)


def draw_flatline_badge(gc, bx, by, d, color):
    """Flat horizontal-line badge for the '+ BG' (background) pill button."""
    y = by + d * 0.5
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(color).Width(max(1.4, d * 0.11))
                           .Cap(wx.CAP_ROUND)))
    gc.StrokeLine(bx + d * 0.22, y, bx + d * 0.78, y)


# --------------------------------------------------------------------------- #
# Resize/sash-drag performance helpers.
#
# Neither of these change what gets drawn -- only when. A live window
# resize or splitter sash drag fires a burst of many EVT_SIZE (or
# EVT_SPLITTER_SASH_POS_CHANGING) events per second; without debouncing,
# each one redraws the matplotlib canvas (a full Agg re-render of the 2D
# spectrum/contours/colorbar -- genuinely expensive) and every owner-drawn
# pill/icon control repaints too. Both helpers collapse a whole burst down
# to a single redraw/repaint once the burst settles.
# --------------------------------------------------------------------------- #
class CanvasRedrawDebouncer:
    """Suppresses a matplotlib FigureCanvasWxAgg's actual (expensive) Agg
    re-render during a burst of resize events, doing exactly one real
    redraw ~delay_ms after the last event in the burst instead of once per
    intermediate tick.

    matplotlib's own canvas._on_size() already resizes the Figure and asks
    to redraw on every tick (see backend_wx._FigureCanvasWxBase._on_size /
    draw_idle), but the actual rendering only happens lazily in
    _on_paint(), which calls canvas.draw() whenever canvas._isDrawn is
    False (true throughout an active resize/redraw burst). Swapping
    canvas.draw for a no-op for the duration of the burst makes that
    per-tick paint effectively free -- the last-rendered bitmap just stays
    on screen (stretched/clipped as the widget resizes, with any newly
    exposed margin blank until the real redraw) -- then restoring the real
    draw() and calling it once when the burst settles gives one correct,
    full-quality render instead of dozens.

    Usage: deb = CanvasRedrawDebouncer(canvas); bind deb.on_size to
    EVT_SIZE on the panel hosting the canvas (or the canvas itself).
    """

    def __init__(self, canvas, delay_ms=150):
        self.canvas = canvas
        self._delay_ms = delay_ms
        self._suspended = False
        self._timer = wx.Timer(canvas)
        canvas.Bind(wx.EVT_TIMER, self._on_settle, self._timer)
        canvas.Bind(wx.EVT_WINDOW_DESTROY, self._on_destroy)

    def on_size(self, event):
        if not self._suspended:
            self._suspended = True
            self.canvas.draw = lambda *a, **k: None
        self._timer.Start(self._delay_ms, wx.TIMER_ONE_SHOT)
        event.Skip()

    def _resume(self):
        if self._suspended:
            try:
                del self.canvas.draw
            except AttributeError:
                pass
            self._suspended = False

    def _on_settle(self, event):
        self._resume()
        try:
            self.canvas.draw()
        except RuntimeError:
            pass  # canvas was destroyed while the settle timer was pending

    def _on_destroy(self, event):
        self._timer.Stop()
        self._resume()
        event.Skip()


class ResizeFreezeGuard:
    """Freeze()s a top-level window once at the start of a burst of
    resize/sash-drag events and Thaw()s it once ~delay_ms after the last
    one, instead of letting every intermediate tick repaint the window's
    owner-drawn pill/icon controls individually. Complementary to
    CanvasRedrawDebouncer: this suppresses the *screen blit* of everything
    else in the window, while that one suppresses the expensive matplotlib
    *rendering* work specifically.

    Usage: guard = ResizeFreezeGuard(frame); bind guard.on_event to
    EVT_SIZE on the frame and to EVT_SPLITTER_SASH_POS_CHANGING on any
    splitters inside it (a sash drag alone doesn't resize the frame, so it
    needs its own trigger). Freeze()/Thaw() calls are balanced exactly
    once per burst; EVT_CLOSE is handled defensively so the window can
    never be left frozen if it's closed mid-drag.
    """

    def __init__(self, window, delay_ms=150):
        self.window = window
        self._delay_ms = delay_ms
        self._frozen = False
        self._timer = wx.Timer(window)
        window.Bind(wx.EVT_TIMER, self._on_settle, self._timer)
        window.Bind(wx.EVT_CLOSE, self._on_close)

    def on_event(self, event):
        if not self._frozen:
            self._frozen = True
            self.window.Freeze()
        self._timer.Start(self._delay_ms, wx.TIMER_ONE_SHOT)
        event.Skip()

    def _thaw(self):
        if self._frozen:
            self._frozen = False
            try:
                self.window.Thaw()
            except RuntimeError:
                pass  # window was destroyed

    def _on_settle(self, event):
        self._thaw()

    def _on_close(self, event):
        self._timer.Stop()
        self._thaw()
        event.Skip()
