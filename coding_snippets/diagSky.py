import numpy as np
import matplotlib.pyplot as plt

# --- 45-degree discrete rotation (vectorized) ---
def rotate45_no_interpolation(A):
    N = A.shape[0]
    out_size = 2 * N - 1
    out = np.zeros((out_size, out_size), dtype=A.dtype)

    i, j = np.indices((N, N))
    x = i + j
    y = i - j + (N - 1)

    out[x, y] = A
    return out

# --- Create a stretched 2D Gaussian ---
N = 101
x = np.linspace(-2, 2, N)
y = np.linspace(-1, 1, N)   # different scale → stretched shape
X, Y = np.meshgrid(x, y)

A = np.exp(-(X**2 / 0.4 + Y**2 / 0.05))  # elliptical Gaussian

# Rotate
A_rot = rotate45_no_interpolation(A)

# --- Plotting ---
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

ax[0].imshow(A, cmap='hot', origin='lower')
ax[0].set_title("Original Stretched Gaussian")

ax[1].imshow(A_rot, cmap='hot', origin='lower')
ax[1].set_title("Discrete 45° Rotation (no interpolation)")

plt.tight_layout()
plt.show()
