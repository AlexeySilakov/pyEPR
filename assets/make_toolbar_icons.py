"""Generates flat, one-tone (single accent colour, transparent background)
toolbar icons for pyIR_peaksearch: undo, redo, save, load, theme toggle,
and the arPLS cutoff-scan tool. Filled-glyph style, rendered at high
resolution and downsampled for crisp edges (same approach as make_icon.py).
"""
import os

import numpy as np
from PIL import Image, ImageDraw

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


def stamped_stroke(draw, pts, width, n_samples=300):
    """Thick line as a union of stamped filled circles along a densely,
    arc-length-resampled path -- avoids notches at sharp turns."""
    pts = np.asarray(pts, dtype=np.float64)
    seglen = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seglen)])
    total = arc[-1]
    if total == 0:
        return
    targets = np.linspace(0.0, total, n_samples)
    rx = np.interp(targets, arc, pts[:, 0])
    ry = np.interp(targets, arc, pts[:, 1])
    r = width / 2.0
    for px, py in zip(rx, ry):
        draw.ellipse([px - r, py - r, px + r, py + r], fill=255)


def make_undo_arrow(mirror=False):
    m = new_mask()
    d = ImageDraw.Draw(m)
    cx, cy = S * 0.54, S * 0.60
    r = S * 0.26
    width = S * 0.10
    a0, a1 = -25.0, 235.0  # degrees, sweep from lower-right up over the top
    n = 240
    angles = np.radians(np.linspace(a0, a1, n))
    pts = [(cx + r * np.cos(a), cy + r * np.sin(a)) for a in angles]
    stamped_stroke(d, pts, width)

    # arrowhead at the a0 end, pointing tangentially "backwards"
    a_tip = np.radians(a0)
    tangent = np.array([-np.sin(a_tip), np.cos(a_tip)])
    normal = np.array([np.cos(a_tip), np.sin(a_tip)])
    base_pt = np.array([cx + r * np.cos(a_tip), cy + r * np.sin(a_tip)])
    head_len = S * 0.22
    head_w = S * 0.24
    tip = base_pt - tangent * head_len * 0.55
    left = base_pt + tangent * head_len * 0.55 + normal * head_w * 0.5
    right = base_pt + tangent * head_len * 0.55 - normal * head_w * 0.5
    d.polygon([tuple(tip), tuple(left), tuple(right)], fill=255)

    if mirror:
        m = m.transpose(Image.FLIP_LEFT_RIGHT)
    return m


def make_save():
    m = new_mask()
    d = ImageDraw.Draw(m)
    margin = S * 0.15
    d.rounded_rectangle([margin, margin, S - margin, S - margin],
                         radius=S * 0.07, fill=255)
    # top slot, open to the top edge
    d.rectangle([S * 0.36, margin, S * 0.72, S * 0.37], fill=0)
    # label window, lower half
    d.rounded_rectangle([S * 0.30, S * 0.54, S * 0.70, S * 0.82],
                         radius=S * 0.03, fill=0)
    return m


def make_load():
    m = new_mask()
    d = ImageDraw.Draw(m)
    # folder back (tab + body top edge)
    d.rounded_rectangle([S * 0.12, S * 0.16, S * 0.46, S * 0.27],
                         radius=S * 0.02, fill=255)
    d.rounded_rectangle([S * 0.12, S * 0.22, S * 0.88, S * 0.34],
                         radius=S * 0.02, fill=255)
    # front flap, trapezoid widening toward the bottom
    d.polygon([(S * 0.06, S * 0.36), (S * 0.94, S * 0.36),
               (S * 0.84, S * 0.86), (S * 0.16, S * 0.86)], fill=255)
    return m


def make_theme():
    m = new_mask()
    d = ImageDraw.Draw(m)
    cx, cy = S * 0.5, S * 0.5
    r = S * 0.22
    # right half: sun disc + rays
    d.pieslice([cx - r, cy - r, cx + r, cy + r], -90, 90, fill=255)
    for deg in (-60, -30, 0, 30, 60):
        a = np.radians(deg)
        x0, y0 = cx + r * 1.15 * np.cos(a), cy + r * 1.15 * np.sin(a)
        x1, y1 = cx + r * 1.62 * np.cos(a), cy + r * 1.62 * np.sin(a)
        d.line([x0, y0, x1, y1], fill=255, width=int(S * 0.065))
    # left half: crescent moon (circle minus an offset circle)
    d.pieslice([cx - r, cy - r, cx + r, cy + r], 90, 270, fill=255)
    dx = r * 0.38
    d.ellipse([cx + dx - r, cy - r, cx + dx + r, cy + r], fill=0)
    return m


def make_scan():
    m = new_mask()
    d = ImageDraw.Draw(m)
    lcx, lcy, lr = S * 0.42, S * 0.42, S * 0.27
    ring_w = S * 0.08
    d.ellipse([lcx - lr, lcy - lr, lcx + lr, lcy + lr], fill=255)
    d.ellipse([lcx - lr + ring_w, lcy - lr + ring_w,
               lcx + lr - ring_w, lcy + lr - ring_w], fill=0)
    # handle
    hx0, hy0 = lcx + lr * 0.80, lcy + lr * 0.80
    hx1, hy1 = S * 0.90, S * 0.90
    stamped_stroke(d, [(hx0, hy0), (hx1, hy1)], S * 0.095, n_samples=30)
    # small wave inside the lens, showing a scan sweep
    inner_r = lr - ring_w
    xs = np.linspace(lcx - inner_r * 0.70, lcx + inner_r * 0.70, 60)
    t = (xs - xs[0]) / (xs[-1] - xs[0])
    ys = lcy + np.sin(t * 2 * np.pi * 1.4) * inner_r * 0.36
    pts = list(zip(xs.tolist(), ys.tolist()))
    stamped_stroke(d, pts, S * 0.055)
    return m


finalize(make_undo_arrow(mirror=False), "undo")
finalize(make_undo_arrow(mirror=True), "redo")
finalize(make_save(), "save")
finalize(make_load(), "load")
finalize(make_theme(), "theme")
finalize(make_scan(), "scan")
