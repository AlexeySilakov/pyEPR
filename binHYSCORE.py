#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bin_hyscore_demo.py

Demonstrates the speed of the vectorised `bin_hyscore` routine.
"""

from __future__ import annotations

import numpy as np
import time
from typing import Union
import matplotlib.pyplot as plt

def bin_hyscore(
    spec: np.ndarray,
    omega_a: np.ndarray,
    omega_b: np.ndarray,
    amp: np.ndarray,
    max_freq: float,
) -> None:
    """
    Bin a list of peaks into a 2‑D complex spectrum.

    Parameters
    ----------
    spec : ndarray, shape (N,N), dtype=complex128
        The spectrum that will be modified in place.
    omega_a : ndarray, shape (M,), dtype=float64
        First frequency axis of the peaks.
    omega_b : ndarray, shape (M,), dtype=float64
        Second frequency axis of the peaks. Must have same length as `omega_a`.
    amp : ndarray, shape (M,) or (M,2), dtype=complex128 or float64
        Peak amplitudes.  If real‑valued it is treated as a complex array with zero imaginary part.
    max_freq : float
        The maximum frequency that maps to the edge of the spectrum.

    Notes
    -----
    * `spec` must be a square, complex 2‑D array.
    * The routine works in place – no new array is allocated for the result.
    * The algorithm is fully vectorised: all peak indices are computed at once,
      then we use NumPy’s advanced indexing to add the amplitudes.

    Raises
    ------
    ValueError
        If input shapes or types do not match expectations.
    """
    # ---------- sanity checks ----------
    if spec.ndim != 2:
        raise ValueError("spec must be a 2‑D array")
    n_points = spec.shape[0]
    if spec.shape[1] != n_points:
        raise ValueError("spec must be square")
    if not np.iscomplexobj(spec):
        raise ValueError("spec must be complex")

    # Ensure all inputs are NumPy arrays of the right dtype
    #omega_a = np.asarray(omega_a, dtype=np.float64)
    #omega_b = np.asarray(omega_b, dtype=np.float64)

    if omega_a.ndim != 1 or omega_b.ndim != 1:
        raise ValueError("omega_a and omega_b must be 1‑D arrays")
    if omega_a.shape[0] != omega_b.shape[0]:
        raise ValueError("omega_a and omega_b must have the same length")

    # amplitude can be real or complex
    amp = np.asarray(amp, dtype=np.complex128)
    if amp.ndim != 1:
        raise ValueError("amp must be a 1‑D array")
    if amp.shape[0] != omega_a.shape[0]:
        raise ValueError("amp must have the same length as omega_a/b")

    # ---------- compute mapping ----------
    n_points_d = float(n_points)
    dx = (n_points_d - 1.0) / (2.0 * max_freq)

    # Convert frequencies to integer indices
    idx1 = np.round((omega_a + max_freq) * dx).astype(np.int64)
    idx2 = np.round((omega_b + max_freq) * dx).astype(np.int64)

    # Keep only peaks that fall inside the array bounds
    valid = (idx1 >= 0) & (idx1 < n_points) & (idx2 >= 0) & (idx2 < n_points)
    if not np.any(valid):
        return  # nothing to add

    idx1 = idx1[valid]
    idx2 = idx2[valid]
    amp_valid = amp[valid]

    # Linear index in column‑major order (as MATLAB does)
    linear_idx = idx1 + idx2 * n_points
    # Symmetric counterpart
    sym_linear_idx = idx2 + idx1 * n_points

    # ---------- add contributions ----------
    np.add.at(spec.ravel(), linear_idx, amp_valid)
    np.add.at(spec.ravel(), sym_linear_idx, amp_valid)


# ------------------------------------------------------------------
if __name__ == "__main__":
    
        # -------------------- demo parameters --------------------
    N = 512                     # spectrum size (N x N)
    M = 200_000                 # number of peaks
    max_freq = 1.0              # arbitrary units
    
    # ------------------------------------------------------------
    print(f"Creating a {N}x{N} complex spectrum and {M:,} random peaks…")
    # Random seed for reproducibility
    rng = np.random.default_rng(seed=42)

    spec = np.zeros((N, N), dtype=np.complex128)
    # -------------------- generate five‑star peaks ---------------
    # Parameters for the star
    n_points_star = 5000          # number of points per arm
    radius        = max_freq * 0.8   # how far from centre

    # Angles that produce a five‑pointed star (72° apart)
    angles = np.linspace(0, 2*np.pi, 5, endpoint=False)

    omega_a_list = []
    omega_b_list = []

    for ang in angles:
        # Each arm is a straight line from the centre to the tip
        t = np.linspace(0, radius, n_points_star)
        omega_a_list.append(t * np.cos(ang))
        omega_b_list.append(t * np.sin(ang))

    # Concatenate all arms into single arrays
    omega_a = np.concatenate(omega_a_list)
    omega_b = np.concatenate(omega_b_list)

    # Random amplitudes (all positive for a clean density plot)
    amp = np.ones(omega_a.size) + 0j

    # -------------------- binning ---------------------------------
    for ii in range(10):
        start = time.perf_counter()
        bin_hyscore(spec, omega_a, omega_b, amp, max_freq)
        elapsed = time.perf_counter() - start
        print(f"BinHYSCORE finished in {elapsed*1000:.2f} ms")

    # -------------------- plotting ---------------------------------
    data_to_plot = np.abs(spec)

    plt.figure(figsize=(6, 5))
    plt.imshow(
        data_to_plot,
        extent=[-max_freq, max_freq, -max_freq, max_freq],
        origin="lower",
        cmap="viridis",
        interpolation="nearest",
    )
    plt.colorbar(label=r"|Spec|")
    plt.title("Five‑pointed star – Binned HYSCORE spectrum")
    plt.xlabel(r"$\omega_a$")
    plt.ylabel(r"$\omega_b$")
    plt.tight_layout()
    plt.show()
