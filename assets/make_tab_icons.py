"""Generates the gear / mountain glyphs used in the folder-style tab
bar's icon badges. Same supersample-then-downsample approach as the
other assets/make_*_icons.py scripts. Each glyph is baked in two
tones: the active-tab brand teal, and a single neutral grey used for
the inactive tab on both light and dark backgrounds.
"""
import os

import numpy as np
from PIL import Image, ImageDraw

S = 512
OUT = 96

ACTIVE = (26, 148, 151)     # matches the toolbar icon teal
INACTIVE = (128, 134, 143)  # neutral grey, legible on both inactive fills

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


def make_gear():
    m = new_mask()
    d = ImageDraw.Draw(m)
    cx, cy = S * 0.5, S * 0.5
    r_rim = S * 0.27
    r_tooth_out = S * 0.365
    tooth_w = S * 0.11
    teeth = 8
    for i in range(teeth):
        ang = np.radians(i * 360 / teeth)
        cos_a, sin_a = np.cos(ang), np.sin(ang)
        hw = tooth_w / 2
        pts = []
        for lx, ly in ((-hw, r_rim * 0.55), (hw, r_rim * 0.55),
                       (hw, r_tooth_out), (-hw, r_tooth_out)):
            rx = lx * cos_a - ly * sin_a
            ry = lx * sin_a + ly * cos_a
            pts.append((cx + rx, cy + ry))
        d.polygon(pts, fill=255)
    d.ellipse([cx - r_rim, cy - r_rim, cx + r_rim, cy + r_rim], fill=255)
    r_hole = S * 0.12
    d.ellipse([cx - r_hole, cy - r_hole, cx + r_hole, cy + r_hole], fill=0)
    return m


def make_mountain():
    m = new_mask()
    d = ImageDraw.Draw(m)
    pts = [(S * 0.06, S * 0.82), (S * 0.34, S * 0.38), (S * 0.50, S * 0.60),
           (S * 0.66, S * 0.26), (S * 0.94, S * 0.82)]
    d.polygon(pts, fill=255)
    return m


gear_mask = make_gear()
mountain_mask = make_mountain()

save_colored(gear_mask, ACTIVE, "tab_gear_active")
save_colored(gear_mask, INACTIVE, "tab_gear_inactive")
save_colored(mountain_mask, ACTIVE, "tab_mountain_active")
save_colored(mountain_mask, INACTIVE, "tab_mountain_inactive")
