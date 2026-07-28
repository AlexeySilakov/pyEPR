"""Generates flat, one-tone toolbar icons for pyIR_plot (the FTIR
Spectrum Viewer companion app): loading a .pir (FTIR_peakfit session)
or .xir (this viewer's own arrangement) file, and saving an .xir
arrangement or exporting the figure as a PDF. Same technique as
make_toolbar_icons.py (mask -> recolour -> downsample), and reuses that
script's folder/floppy-disk base shapes so the two apps' icon sets read
as one family. Output goes into the same shared assets/toolbar/ folder.
"""
import os

from PIL import Image, ImageDraw, ImageFont

S = 512    # supersample canvas
OUT = 128  # saved master size
ACCENT = (75, 75, 150, 255)  # single tone, legible on light and dark bg

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "toolbar")
os.makedirs(OUT_DIR, exist_ok=True)


def new_mask():
    return Image.new("L", (S, S), 0)


def finalize(mask, name):
    color = Image.new("RGBA", (S, S), ACCENT)
    color.putalpha(mask)
    img = color.resize((OUT, OUT), Image.LANCZOS)
    img.save(os.path.join(OUT_DIR, f"{name}.png"))
    print("saved", name, img.size)
    return img


def folder_shape():
    """Same open-folder silhouette as make_toolbar_icons.make_load()."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    d.rounded_rectangle([S * 0.12, S * 0.16, S * 0.46, S * 0.27],
                         radius=S * 0.02, fill=255)
    d.rounded_rectangle([S * 0.12, S * 0.22, S * 0.88, S * 0.34],
                         radius=S * 0.02, fill=255)
    d.polygon([(S * 0.06, S * 0.36), (S * 0.94, S * 0.36),
               (S * 0.84, S * 0.86), (S * 0.16, S * 0.86)], fill=255)
    return m


def floppy_shape():
    """Same floppy-disk silhouette as make_toolbar_icons.make_save(),
    but WITHOUT cutting the label window -- callers cut their own mark
    into that area instead."""
    m = new_mask()
    d = ImageDraw.Draw(m)
    margin = S * 0.15
    d.rounded_rectangle([margin, margin, S - margin, S - margin],
                         radius=S * 0.07, fill=255)
    d.rectangle([S * 0.36, margin, S * 0.72, S * 0.37], fill=0)
    return m


def cut_mountain(mask, cx, cy, w, h):
    """A small twin-peak notch, echoing the Peak-fitting tab icon --
    marks a glyph as being about a .pir (FTIR_peakfit session) file."""
    d = ImageDraw.Draw(mask)
    pts = [(cx - w * 0.46, cy + h * 0.40), (cx - w * 0.14, cy - h * 0.30),
           (cx + w * 0.02, cy + h * 0.02), (cx + w * 0.24, cy - h * 0.40),
           (cx + w * 0.48, cy + h * 0.40)]
    d.polygon(pts, fill=0)


def cut_stack(mask, cx, cy, w, h, n=3):
    """A small stacked-lines notch, echoing this viewer's core
    "stack multiple spectra" concept -- marks a glyph as being about a
    .xir (this app's own arrangement) file."""
    d = ImageDraw.Draw(mask)
    bar_h = h * 0.16
    step = (h - bar_h) / (n - 1) if n > 1 else 0
    top = cy - h / 2
    for i in range(n):
        y = top + i * step
        d.rounded_rectangle([cx - w / 2, y, cx + w / 2, y + bar_h],
                             radius=bar_h * 0.3, fill=0)


def make_load_pir():
    m = folder_shape()
    cut_mountain(m, S * 0.5, S * 0.62, S * 0.30, S * 0.22)
    return m


def make_load_xir():
    m = folder_shape()
    cut_stack(m, S * 0.5, S * 0.62, S * 0.36, S * 0.24)
    return m


def make_save_xir():
    m = floppy_shape()
    d = ImageDraw.Draw(m)
    d.rounded_rectangle([S * 0.30, S * 0.54, S * 0.70, S * 0.82],
                         radius=S * 0.03, fill=255)
    cut_stack(m, S * 0.5, S * 0.68, S * 0.28, S * 0.20)
    return m


def make_save_pdf():
    """Same floppy-disk base as make_save_xir(), pairing the two "save"
    icons visually; the label window carries "PDF" text instead of the
    stacked-lines mark."""
    m = floppy_shape()
    d = ImageDraw.Draw(m)
    # label window widened so the bigger "PDF" text has room
    d.rounded_rectangle([S * 0.20, S * 0.52, S * 0.80, S * 0.84],
                         radius=S * 0.03, fill=255)
    font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", int(S * 0.24))
    text = "PDF"
    bbox = d.textbbox((0, 0), text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    tx = S * 0.5 - tw / 2 - bbox[0]
    ty = S * 0.68 - th / 2 - bbox[1]
    d.text((tx, ty), text, font=font, fill=0)
    return m


finalize(make_load_pir(), "load_pir")
finalize(make_load_xir(), "load_xir")
finalize(make_save_xir(), "save_xir")
finalize(make_save_pdf(), "save_pdf")
