"""Generates the small glyphs used inside the white circle badge of the
pill-style form buttons (load actions, spike-select toggle, auto-send
toggle). Same supersample-then-downsample approach as make_icon.py /
make_toolbar_icons.py, but each glyph is baked in a fixed accent colour
since it always sits on a plain white badge.
"""
import os

import numpy as np
from PIL import Image, ImageDraw

S = 512
OUT = 96

LOAD = (41, 128, 185)   # action blue, used for the three "Load ..." buttons
CALM = (39, 174, 96)    # toggle OFF -- calm green
WARN = (216, 89, 30)    # toggle ON -- warning orange
TEAL = (26, 148, 151)   # brand teal, used for the "Save baselined data" button

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "badges")
os.makedirs(OUT_DIR, exist_ok=True)


def new_mask():
    return Image.new("L", (S, S), 0)


def save_colored(mask, color, name):
    layer = Image.new("RGBA", (S, S), (*color, 255))
    layer.putalpha(mask)
    img = layer.resize((OUT, OUT), Image.LANCZOS)
    img.save(os.path.join(OUT_DIR, f"{name}.png"))
    print("saved", name, img.size)


def make_folder():
    """Open-folder glyph, for the Load buttons."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    d.rounded_rectangle([S * 0.12, S * 0.16, S * 0.46, S * 0.27],
                         radius=S * 0.02, fill=255)
    d.rounded_rectangle([S * 0.12, S * 0.22, S * 0.88, S * 0.34],
                         radius=S * 0.02, fill=255)
    d.polygon([(S * 0.06, S * 0.36), (S * 0.94, S * 0.36),
               (S * 0.84, S * 0.86), (S * 0.16, S * 0.86)], fill=255)
    return m


def make_target():
    """Crosshair / reticle glyph, for the spike-select toggle."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    cx, cy = S * 0.5, S * 0.5
    r_out = S * 0.26
    ring_w = S * 0.065
    d.ellipse([cx - r_out, cy - r_out, cx + r_out, cy + r_out], fill=255)
    r_in = r_out - ring_w
    d.ellipse([cx - r_in, cy - r_in, cx + r_in, cy + r_in], fill=0)
    r_dot = S * 0.06
    d.ellipse([cx - r_dot, cy - r_dot, cx + r_dot, cy + r_dot], fill=255)
    gap, tick_len, tick_w = S * 0.045, S * 0.11, S * 0.05
    for ang in (0, 90, 180, 270):
        a = np.radians(ang)
        x0, y0 = cx + (r_out + gap) * np.cos(a), cy + (r_out + gap) * np.sin(a)
        x1 = cx + (r_out + gap + tick_len) * np.cos(a)
        y1 = cy + (r_out + gap + tick_len) * np.sin(a)
        d.line([x0, y0, x1, y1], fill=255, width=int(tick_w))
    return m


def make_send():
    """Paper-plane / dart glyph, for the auto-send toggle."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    pts = [(S * 0.86, S * 0.50), (S * 0.16, S * 0.20),
           (S * 0.42, S * 0.50), (S * 0.16, S * 0.80)]
    d.polygon(pts, fill=255)
    return m


def make_floppy():
    """Floppy-disk glyph, for the Save baselined data button."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    margin = S * 0.15
    d.rounded_rectangle([margin, margin, S - margin, S - margin],
                         radius=S * 0.07, fill=255)
    d.rectangle([S * 0.36, margin, S * 0.72, S * 0.37], fill=0)
    d.rounded_rectangle([S * 0.30, S * 0.54, S * 0.70, S * 0.82],
                         radius=S * 0.03, fill=0)
    return m


save_colored(make_folder(), LOAD, "load")
save_colored(make_target(), CALM, "spikes_off")
save_colored(make_target(), WARN, "spikes_on")
save_colored(make_send(), CALM, "send_off")
save_colored(make_send(), WARN, "send_on")
save_colored(make_floppy(), TEAL, "save_teal")
