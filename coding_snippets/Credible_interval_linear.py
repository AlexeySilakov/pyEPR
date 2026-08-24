import numpy as np
import matplotlib.pyplot as plt

# ---- 1. make noisy data from y = a*x + b*x^2 (no constant term) ----
np.random.seed(0)
a_true = 2.0
b_true = -0.5
sigma_true = 3.0
conf_int = 0.98
x = np.linspace(-5, 5, 100)
y = a_true*x + b_true*x**2 + np.random.normal(0, sigma_true, size=x.size)
y_true =a_true*x + b_true*x**2
# ---- 2. grid search over a and b ----
a_grid = np.linspace(0.5, 3.5, 200)
b_grid = np.linspace(-1.2, 0.2, 200)
A, B = np.meshgrid(a_grid, b_grid)                       # shape (nb, na)
model = A[:, :, None]*x + B[:, :, None]*x**2             # (nb, na, N)
SSR = np.sum((y - model)**2, axis=2)                     # residual sum of squares
RMSD = np.sqrt(SSR / x.size)

# ---- 3. estimate noise sigma FROM the residuals of the best grid fit ----
imin = np.unravel_index(np.argmin(SSR), SSR.shape)
N, p = x.size, 2
sigma_hat = np.sqrt(SSR[imin] / (N - p))                 # NOT sigma_true

# ---- 4. likelihood + uniform prior -> normalised posterior on the grid ----
L = np.exp(-(SSR-SSR[imin]) / (2*sigma_hat**2))        # minus min = numerical stability
post = L / L.sum()                                       # uniform prior (1), then normalise

# ---- 5. highest-posterior-density region ----
flat = np.sort(post.ravel())[::-1]
cutoff = flat[np.searchsorted(np.cumsum(flat), conf_int)]
mask = post >= cutoff

# ---- 5b. marginal (per-parameter) 95% equal-tailed intervals ----
p_a = post.sum(axis=0); p_a /= p_a.sum()                 # marginalise over b
p_b = post.sum(axis=1); p_b /= p_b.sum()                 # marginalise over a
a_lo, a_hi = np.interp([0.025, 0.975], np.cumsum(p_a), a_grid)
b_lo, b_hi = np.interp([0.025, 0.975], np.cumsum(p_b), b_grid)

# ---- plots ----

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), num='bayes_grid', clear=True)

ax[0].plot(x, y, 'o', label='noisy data')
ax[0].plot(x, y_true, ':', label='data')
ax[0].plot(x, a_grid[imin[1]]*x + b_grid[imin[0]]*x**2, 'r-', lw=2, label='grid best fit')
ax[0].set_title('data & fit'); ax[0].set_xlabel('x'); ax[0].set_ylabel('y'); ax[0].legend()

pc = ax[1].pcolormesh(a_grid, b_grid, post, shading='auto', cmap='viridis')
ax[1].contour(a_grid, b_grid, mask.astype(float), levels=[0.5], colors='white', linewidths=2)
# ax[1].plot(a_grid[imin[1]], b_grid[imin[0]], 'r*', ms=14)
ax[1].axhline(b_true)
ax[1].axvline(a_true)

ax[1].set_title(f'a: ±{(a_hi-a_lo)/2:.2f}   b: ±{(b_hi-b_lo)/2:.2f} σ: {sigma_hat:.2f}')
ax[1].set_xlabel('a'); ax[1].set_ylabel('b'); fig.colorbar(pc, ax=ax[1])

rc = ax[2].pcolormesh(a_grid, b_grid, RMSD, shading='auto', cmap='magma')
ax[2].plot(a_grid[imin[1]], b_grid[imin[0]], 'c*', ms=14)
ax[2].set_title('RMSD'); ax[2].set_xlabel('a'); ax[2].set_ylabel('b'); fig.colorbar(rc, ax=ax[2])

plt.tight_layout()
plt.savefig('bayes_grid.png', dpi=130)

print(f"true a, b     = {a_true}, {b_true}")
print(f"best-fit a, b = {a_grid[imin[1]]:.3f}, {b_grid[imin[0]]:.3f}")
print(f"sigma_true    = {sigma_true}")
print(f"sigma_hat     = {sigma_hat:.3f}   (estimated from residuals)")