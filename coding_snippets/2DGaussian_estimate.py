"""
2D Gaussian parameters from a log-probability grid, without curve_fit.

Method A (parabola / Laplace): fit a 6-term quadratic to log p near the mode.
                               Hessian H -> Sigma = -inv(H), mu = -inv(H) @ b
Method B (weighted moments):   first and second moments of exp(log p).

The surface here is deliberately NOT a pure Gaussian: it is a narrow core plus
a broad, slightly offset component. That is what makes A and B disagree --
A measures the curvature AT the mode, B measures the moments OF the whole
distribution. Set MIX = 0.0 to recover an exact Gaussian, where both methods
return Sigma_core to within grid resolution.

Linear script, no function definitions.
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# synthetic log-probability surface (built in log throughout)
# ----------------------------------------------------------------------

x0 = 1.5
y0 = -0.8
sigmaX = 0.8
sigmaY = 1.2
xg = np.linspace(-2.5, 5.5, 300)
yg = np.linspace(-3.8, 2.2, 280)
X, Y = np.meshgrid(xg, yg, indexing="ij")
p = np.exp(-(X-x0)**2/sigmaX**2/2)*np.exp(-(Y-y0)**2/sigmaY**2/2)
logp = np.log(p)

# ----------------------------------------------------------------------
# METHOD A: quadratic fit to log p near the mode
# ----------------------------------------------------------------------

ii, jj = np.unravel_index(np.argmax(logp), logp.shape)
cx = np.polyfit(xg, logp[:, jj], 2)
cy = np.polyfit(yg, logp[ii, :], 2)

x0_fit = -cx[1] / (2 * cx[0])
sigmaX_fit = np.sqrt(-1 / (2 * cx[0]))
y0_fit = -cy[1] / (2 * cy[0])
sigmaY_fit = np.sqrt(-1 / (2 * cy[0]))


print(f'DET sigmax = {sigmaX_fit:.3f}  REAL sigmax = {sigmaX:.3f}')
print(f'DET sigmay = {sigmaY_fit:.3f}  REAL sigmay = {sigmaY:.3f}')
print(f'DET x0 = {x0_fit:.3f}  REAL sigmax = {x0:.3f}')
print(f'DET y0 = {y0_fit:.3f}  REAL sigmay = {y0:.3f}')

# center and rescale before building the design matrix, or the columns
# 1, u, u**2 become badly conditioned for grids far from the origin
x0, y0 = X.mean(), Y.mean()
sx, sy = X.std(), Y.std()
u = (X.ravel() - x0) / sx
v = (Y.ravel() - y0) / sy
z = logp.ravel()




A = np.column_stack([np.ones(u.shape), u, v, u**2, v**2, u*v])
w = np.exp(z)                               # weight residuals toward the mode
coef, *_ = np.linalg.lstsq(A * w[:, None], z * w, rcond=None)
a_, b1, b2, c11, c22, c12 = coef

H = np.array([[2*c11, c12], [c12, 2*c22]])
evals = np.linalg.eigvalsh(H)
if not np.all(evals < 0):
    raise RuntimeError(f"Hessian not negative definite (eigenvalues {evals}) -- "
                       "the fit region is not a single peak")

Hinv = np.linalg.inv(H)
S = np.diag([sx, sy])
Sigma_par = S @ (-Hinv) @ S
mu_par = np.array([x0, y0]) + S @ (-Hinv @ np.array([b1, b2]))

sxp, syp = np.sqrt(np.diag(Sigma_par))
rhop = Sigma_par[0, 1] / (sxp * syp)


# ----------------------------------------------------------------------
# METHOD B: weighted moments of exp(log p)
# ----------------------------------------------------------------------

tot = p.sum()
mu_mom = np.array([np.sum(p*X), np.sum(p*Y)]) / tot
ex, ey = X - mu_mom[0], Y - mu_mom[1]
Sigma_mom = np.array([[np.sum(p*ex*ex), np.sum(p*ex*ey)],
                      [np.sum(p*ex*ey), np.sum(p*ey*ey)]]) / tot

sxm, sym = np.sqrt(np.diag(Sigma_mom))
rhom = Sigma_mom[0, 1] / (sxm * sym)

# ----------------------------------------------------------------------
# report
# ----------------------------------------------------------------------
# print(f"grid: {X.shape[0]} x {X.shape[1]}, "
#       f"cells per sigma_x = {sigmaX:.1f}, "
#       f"cells per sigma_y = {sigmaX:.1f}")
# print(f"{'':<22}{'mu_x':>8}{'mu_y':>8}{'sigma_x':>9}{'sigma_y':>9}{'rho':>8}")
# print("-" * 64)
# print(f"{'core component':<22}{mu_core[0]:>8.3f}{mu_core[1]:>8.3f}"
#       f"{sx_core:>9.3f}{sy_core:>9.3f}{rho_core:>8.3f}")
# print(f"{'A: parabola/Laplace':<22}{mu_par[0]:>8.3f}{mu_par[1]:>8.3f}"
#       f"{sxp:>9.3f}{syp:>9.3f}{rhop:>8.3f}")
# print(f"{'B: weighted moments':<22}{mu_mom[0]:>8.3f}{mu_mom[1]:>8.3f}"
#       f"{sxm:>9.3f}{sym:>9.3f}{rhom:>8.3f}")

# ----------------------------------------------------------------------
# sweep the fit region: a stable Sigma means the surface really is Gaussian
# ----------------------------------------------------------------------
cuts = np.linspace(-0.8, -12.0, 30)
sweep = np.full((cuts.size, 3), np.nan)
for i in range(cuts.size):
    mi = np.isfinite(logp) & (logp > cuts[i])
    if mi.sum() < 40:
        continue
    x0i, y0i = X[mi].mean(), Y[mi].mean()
    sxi, syi = X[mi].std(), Y[mi].std()
    ui = (X[mi] - x0i) / sxi
    vi = (Y[mi] - y0i) / syi
    zi = logp[mi]
    Ai = np.column_stack([np.ones(ui.size), ui, vi, ui**2, vi**2, ui*vi])
    wi = np.exp(zi)
    ci, *_ = np.linalg.lstsq(Ai * wi[:, None], zi * wi, rcond=None)
    Hi = np.array([[2*ci[3], ci[5]], [ci[5], 2*ci[4]]])
    if not np.all(np.linalg.eigvalsh(Hi) < 0):
        continue
    Si = np.diag([sxi, syi])
    Ci = Si @ (-np.linalg.inv(Hi)) @ Si
    sweep[i, 0] = np.sqrt(Ci[0, 0])
    sweep[i, 1] = np.sqrt(Ci[1, 1])
    sweep[i, 2] = Ci[0, 1] / np.sqrt(Ci[0, 0] * Ci[1, 1])

# ----------------------------------------------------------------------
# quadratic model evaluated back on the full grid (for the residual panel)
# ----------------------------------------------------------------------
U = (X - x0) / sx
V = (Y - y0) / sy
model = a_ + b1*U + b2*V + c11*U**2 + c22*V**2 + c12*U*V
resid = logp - model

# ----------------------------------------------------------------------
# 1-sigma ellipses, drawn from the eigen-decomposition of each Sigma
# ----------------------------------------------------------------------
t = np.linspace(0.0, 2*np.pi, 200)
circle = np.vstack([np.cos(t), np.sin(t)])

# ev, evec = np.linalg.eigh(Sigma_core)
# ell_core = mu_core[:, None] + evec @ np.diag(np.sqrt(ev)) @ circle
# ev, evec = np.linalg.eigh(Sigma_par)
# ell_par = mu_par[:, None] + evec @ np.diag(np.sqrt(ev)) @ circle
# ev, evec = np.linalg.eigh(Sigma_mom)
# ell_mom = mu_mom[:, None] + evec @ np.diag(np.sqrt(ev)) @ circle

# ----------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------
# fig, ax = plt.subplots(2, 2, num ='FU SPIDER', figsize=(12.5, 9))

# # --- the surface, the fit region, and the recovered ellipses
# im = ax[0, 0].pcolormesh(X, Y, np.where(logp > -18, logp, np.nan), shading="auto", cmap="viridis")
# fig.colorbar(im, ax=ax[0, 0], label="log p (peak at 0)")
# # ax[0, 0].contour(X, Y, logp, levels=[CUT], colors="w", linewidths=1.2, linestyles="--")
# # ax[0, 0].plot(ell_core[0], ell_core[1], color="w", lw=3.5, label="core component")
# # ax[0, 0].plot(ell_par[0], ell_par[1], color="crimson", lw=1.8, ls="-", label="A: parabola")
# # ax[0, 0].plot(ell_mom[0], ell_mom[1], color="orange", lw=2.0, ls="--", label="B: moments")
# ax[0, 0].plot(*mu_par, "r+", ms=10)
# ax[0, 0].plot(*mu_mom, "+", color="orange", ms=10)
# ax[0, 0].set_title("log p, fit region (dashed), 1-sigma ellipses")
# ax[0, 0].set_xlabel("x")
# ax[0, 0].set_ylabel("y")
# ax[0, 0].legend(fontsize=8, loc="upper right")

# # --- residual of the quadratic fit: structure here means non-Gaussian
# lim = np.nanpercentile(np.abs(resid), 97)   # clipped, or the region edge saturates it
# im = ax[0, 1].pcolormesh(X, Y, resid, shading="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
# fig.colorbar(im, ax=ax[0, 1], label="log p - quadratic")
# ax[0, 1].set_title("fit residual inside the region\n(structure = surface is not Gaussian)")
# ax[0, 1].set_xlabel("x")
# ax[0, 1].set_ylabel("y")

# # --- stability against the fit region threshold
# ax[1, 0].plot(-cuts, sweep[:, 0], lw=2, label="sigma_x")
# ax[1, 0].plot(-cuts, sweep[:, 1], lw=2, label="sigma_y")
# # ax[1, 0].axhline(sx_core, color="0.4", ls=":", lw=1)
# # ax[1, 0].axhline(sy_core, color="0.4", ls=":", lw=1)
# ax[1, 0].axvline(-CUT, color="k", ls="--", lw=1)
# ax[1, 0].set_xlabel("fit region depth  -log p")
# ax[1, 0].set_ylabel("recovered sigma")
# ax[1, 0].set_title("drift with fit region = non-Gaussianity\n(dotted: core component)")
# ax[1, 0].legend(fontsize=8)
# axr = ax[1, 0].twinx()
# axr.plot(-cuts, sweep[:, 2], lw=1.5, color="green", ls="-.")
# axr.set_ylabel("rho (green)", color="green")

# # --- slice through the mode
# j = np.argmin(np.abs(yg - mu_par[1]))
# ax[1, 1].plot(xg, logp[:, j], color="0.5", lw=1.2, label="log p slice")
# ax[1, 1].plot(xg, model[:, j], color="crimson", lw=2, label="A: fitted parabola")
# ax[1, 1].plot(xg, (lg_c - np.nanmax(logp))[:, j], color="k", lw=1.2, ls=":",
#               label="core component only")
# ax[1, 1].axhline(CUT, color="k", ls="--", lw=0.8)
# ax[1, 1].set_ylim(CUT - 6, 1)
# ax[1, 1].set_xlabel("x")
# ax[1, 1].set_ylabel("log p")
# ax[1, 1].set_title(f"slice at y = {yg[j]:.2f}  (dashed: fit cut)")
# ax[1, 1].legend(fontsize=8)

# plt.tight_layout()
# plt.show()