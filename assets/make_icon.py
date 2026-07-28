"""Generates the pyIR_peaksearch app icon: a spectrum curve with peaks on
a gradient tile. Renders at high resolution with PIL/numpy and downsamples
for a crisp multi-size .ico."""
import numpy as np
from PIL import Image, ImageDraw, ImageFilter

S = 1024  # master canvas size


def make_gradient_bg(size, c_tl, c_br):
    ys, xs = np.mgrid[0:size, 0:size].astype(np.float64)
    t = (xs + ys) / (2 * (size - 1))
    c_tl = np.array(c_tl, dtype=np.float64)
    c_br = np.array(c_br, dtype=np.float64)
    grad = c_tl[None, None, :] + (c_br - c_tl)[None, None, :] * t[..., None]
    return grad.astype(np.uint8)


def rounded_mask(size, radius):
    m = Image.new("L", (size, size), 0)
    d = ImageDraw.Draw(m)
    d.rounded_rectangle([0, 0, size - 1, size - 1], radius=radius, fill=255)
    return m


def stamped_stroke_mask(size, pts, width, n_samples):
    """Thick line as a union of stamped filled circles along a densely,
    arc-length-resampled path -- avoids the self-intersection notches a
    single offset-ribbon polygon gets at sharp peak tips."""
    pts = np.asarray(pts, dtype=np.float64)
    seglen = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seglen)])
    total = arc[-1]
    targets = np.linspace(0.0, total, n_samples)
    rx = np.interp(targets, arc, pts[:, 0])
    ry = np.interp(targets, arc, pts[:, 1])

    m = Image.new("L", (size, size), 0)
    d = ImageDraw.Draw(m)
    r = width / 2.0
    for px, py in zip(rx, ry):
        d.ellipse([px - r, py - r, px + r, py + r], fill=255)
    return m


# ---- background: rounded-square tile, deep navy -> teal gradient ----
bg_rgb = make_gradient_bg(S, (13, 42, 74), (12, 130, 132))
bg = Image.fromarray(bg_rgb, "RGB").convert("RGBA")
mask = rounded_mask(S, radius=int(S * 0.22))
tile = Image.new("RGBA", (S, S), (0, 0, 0, 0))
tile.paste(bg, (0, 0), mask)

# subtle darkening toward the bottom, clipped strictly to the rounded tile
shade_l = Image.new("L", (S, S), 0)
sd = ImageDraw.Draw(shade_l)
sd.ellipse([-S * 0.3, S * 0.55, S * 1.3, S * 1.6], fill=90)
shade_l = shade_l.filter(ImageFilter.GaussianBlur(S * 0.08))
shade_l = Image.composite(shade_l, Image.new("L", (S, S), 0), mask)
shade_rgba = Image.new("RGBA", (S, S), (0, 0, 0, 255))
shade_rgba.putalpha(shade_l)
tile = Image.alpha_composite(tile, shade_rgba)

# ---- spectrum curve: baseline + 3 Gaussian-ish peaks ----
n = 400
x = np.linspace(0.06, 0.94, n) * S


def gauss(xc, amp, fwhm, xs):
    sig = fwhm / 2.3548
    return amp * np.exp(-0.5 * ((xs - xc) / sig) ** 2)


xs_norm = np.linspace(0.06, 0.94, n)
baseline_y = 0.72
y = np.zeros(n)
y += gauss(0.22, 0.17, 0.11, xs_norm)
y += gauss(0.50, 0.44, 0.095, xs_norm)
y += gauss(0.76, 0.24, 0.08, xs_norm)

curve_y = (baseline_y - y) * S
pts = list(zip(x.tolist(), curve_y.tolist()))

# filled area under the curve (soft, low-opacity)
fill_layer = Image.new("RGBA", (S, S), (0, 0, 0, 0))
fd = ImageDraw.Draw(fill_layer)
poly = pts + [(x[-1], baseline_y * S), (x[0], baseline_y * S)]
fd.polygon(poly, fill=(255, 255, 255, 40))
tile = Image.alpha_composite(tile, fill_layer)

# baseline
draw = ImageDraw.Draw(tile)
draw.line([(x[0], baseline_y * S), (x[-1], baseline_y * S)],
          fill=(255, 255, 255, 65), width=max(2, S // 256))

# the curve itself: stamped-circle union, immune to self-intersection
stroke_w = S * 0.026
stroke_mask = stamped_stroke_mask(S, pts, stroke_w, n_samples=1600)
white_layer = Image.new("RGBA", (S, S), (255, 255, 255, 255))
white_layer.putalpha(stroke_mask)
tile = Image.alpha_composite(tile, white_layer)

# final: re-clip to the rounded tile so nothing can ever bleed past corners
out = Image.new("RGBA", (S, S), (0, 0, 0, 0))
out.paste(tile, (0, 0), mask)

out_path = "D:/Programming/Claude/assets/icon_master.png"
out.save(out_path)
print("saved", out_path, out.size)

# ---- multi-size ICO ----
sizes = [16, 24, 32, 48, 64, 128, 256]
imgs = [out.resize((s, s), Image.LANCZOS) for s in sizes]
ico_path = "D:/Programming/Claude/assets/icon.ico"
imgs[-1].save(ico_path, format="ICO", sizes=[(s, s) for s in sizes])
print("saved", ico_path)

png_path = "D:/Programming/Claude/assets/icon_256.png"
imgs[sizes.index(256)].save(png_path)
print("saved", png_path)
