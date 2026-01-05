#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Eigenfield calculation for an S=1/2 – I=5/2 system
    • Zeeman interaction (electron only)
    • Hyperfine interaction (isotropic, A = 100 MHz)
    • g‑tensor = diag(2.0036, 2.0036, 2.0036)
    • Microwave frequency: 9.5 GHz (X‑band)

Author: ChatGPT
"""

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Physical constants
# ------------------------------------------------------------------
mu_B = 9.274009994e-24          # J/T  (Bohr magneton)
hbar = 1.054571817e-34          # J·s
GHz   = 1.0e9
MHz   = 1.0e6

# ------------------------------------------------------------------
# Spin operators
# ------------------------------------------------------------------
def spin_matrices(S):
    """Return Sx, Sy, Sz for spin S (float)."""
    dim = int(2 * S + 1)
    m_vals = np.arange(S, -S - 1, -1)   # |m> basis: S, S-1, … , -S
    Sx = np.zeros((dim, dim), dtype=complex)
    Sy = np.zeros((dim, dim), dtype=complex)
    Sz = np.diag(m_vals)

    for i, mi in enumerate(m_vals):
        if i + 1 < dim:
            mj = m_vals[i + 1]
            coeff = np.sqrt(S * (S + 1) - mi * (mi - 1))
            Sx[i, i + 1] = coeff / 2.0
            Sy[i, i + 1] = -1j * coeff / 2.0
        if i - 1 >= 0:
            mj = m_vals[i - 1]
            coeff = np.sqrt(S * (S + 1) - mi * (mi + 1))
            Sx[i, i - 1] = coeff / 2.0
            Sy[i, i - 1] = 1j * coeff / 2.0

    return Sx, Sy, Sz

# Electron (S = ½) and nuclear (I = 5/2) operators
Sx_e, Sy_e, Sz_e = spin_matrices(0.5)
Sx_n, Sy_n, Sz_n = spin_matrices(5.0 / 2.0)

# ------------------------------------------------------------------
# Hamiltonian builder
# ------------------------------------------------------------------
def hamiltonian(B_vec, g_tensor, A_tensor):
    """
    Full Hamiltonian for a given magnetic‑field vector B.

    Parameters
    ----------
    B_vec : array_like, shape (3,)
        Magnetic‑field vector (T).
    g_tensor : array_like, shape (3,3)
        g‑tensor (dimensionless).
    A_tensor : array_like, shape (3,3)
        Hyperfine tensor (MHz).

    Returns
    -------
    H : ndarray, shape (d, d)
        Hamiltonian in SI units (J).
    """
    B_vec = np.asarray(B_vec, dtype=float)

    # -------- Zeeman (electron only) --------
    # B·g = Σ_i (B_i * g_i·) → 3‑vector
    B_g = B_vec @ g_tensor               # (3,)
    H_Z = mu_B * (
        B_g[0] * np.kron(Sx_e, np.eye(Sx_n.shape[0])) +
        B_g[1] * np.kron(Sy_e, np.eye(Sx_n.shape[0])) +
        B_g[2] * np.kron(Sz_e, np.eye(Sx_n.shape[0]))
    )

    # -------- Hyperfine (S·A·I) --------
    # Convert A from MHz → J:  A_J = 2πħ * A_MHz * 1e6
    A_J = 2.0 * np.pi * hbar * A_tensor * MHz
    H_HF = (
        A_J[0, 0] * np.kron(Sx_e, Sx_n) +
        A_J[0, 1] * np.kron(Sx_e, Sy_n) +
        A_J[0, 2] * np.kron(Sx_e, Sz_n) +

        A_J[1, 0] * np.kron(Sy_e, Sx_n) +
        A_J[1, 1] * np.kron(Sy_e, Sy_n) +
        A_J[1, 2] * np.kron(Sy_e, Sz_n) +

        A_J[2, 0] * np.kron(Sz_e, Sx_n) +
        A_J[2, 1] * np.kron(Sz_e, Sy_n) +
        A_J[2, 2] * np.kron(Sz_e, Sz_n)
    )

    return H_Z + H_HF

# ------------------------------------------------------------------
# Eigenfield routine
# ------------------------------------------------------------------
def resonance_field(transition, freq, g_tensor, A_tensor,
                    theta=0.0, phi=0.0, B_max=0.5, n_points=2000):
    """
    Find the B‑magitude that satisfies hf = ΔE for a chosen transition.

    Parameters
    ----------
    transition : tuple(int, int)
        Indices of the two levels (i < j).
    freq : float
        Microwave frequency in Hz.
    g_tensor, A_tensor : see above
    theta, phi : float
        Orientation of B (radians).
    B_max : float
        Upper bound for B (T).
    n_points : int
        Grid resolution for initial scan.

    Returns
    -------
    B_res : list of float
        Field values (T) where the transition matches the frequency.
    """
    # Cartesian unit vector for B
    n_hat = np.array([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)
    ])

    # Initial scan to locate sign changes
    B_vals = np.linspace(0.0, B_max, n_points)
    delta_E = np.empty_like(B_vals)

    for k, Bmag in enumerate(B_vals):
        H = hamiltonian(Bmag * n_hat, g_tensor, A_tensor)
        w, _ = eigh(H)                      # eigenvalues in J
        w_GHz = w / (hbar * GHz)            # convert to GHz
        delta_E[k] = w_GHz[transition[1]] - w_GHz[transition[0]]

    target = freq / GHz  # GHz

    # Find intervals where ΔE crosses the target
    sign = np.sign(delta_E - target)
    idx = np.where(np.diff(sign))[0]

    if len(idx) == 0:
        raise ValueError(f"No resonance found for transition {transition} within B_max={B_max} T.")

    B_res = []

    # Brent’s method inside each interval
    for i in idx:
        B1, B2 = B_vals[i], B_vals[i + 1]

        def f(B):
            H = hamiltonian(B * n_hat, g_tensor, A_tensor)
            w, _ = eigh(H)
            w_GHz = w / (hbar * GHz)
            return w_GHz[transition[1]] - w_GHz[transition[0]] - target

        B_root = brentq(f, B1, B2)
        B_res.append(B_root)

    return B_res

# ------------------------------------------------------------------
# Example parameters
# ------------------------------------------------------------------
g_tensor = np.diag([2.0036, 2.0036, 2.0036])     # isotropic g
A_tensor = np.diag([100.0, 100.0, 100.0])        # isotropic A in MHz
freq = 9.5 * GHz                                 # X‑band

# ------------------------------------------------------------------
# Define the six allowed Δm_S=±1 transitions (m_I = 5/2 … –5/2)
# Basis order: |mS> first ( +½, –½ ), |mI> last ( +5/2 … –5/2 )
dim_e, dim_n = 2, 6
transitions = []
for n_idx in range(dim_n):                       # same nuclear index
    i = 0 * dim_n + n_idx                        # mS = +½
    j = 1 * dim_n + n_idx                        # mS = –½
    transitions.append((i, j))

# ------------------------------------------------------------------
# Compute resonance fields
# ------------------------------------------------------------------
theta, phi = 0.0, 0.0        # B along z
B_res_list = []

for trans in transitions:
    B_res = resonance_field(
        trans, freq, g_tensor, A_tensor,
        theta=theta, phi=phi, B_max=0.3
    )
    B_res_list.append((trans, B_res))

# ------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------
print("Allowed transitions (index pair) → resonance fields (T):")
for (trans, B_res) in B_res_list:
    print(f"  {trans}: {', '.join(f'{b:.4f}' for b in B_res)}")

# ------------------------------------------------------------------
# Plot energy levels vs B
# ------------------------------------------------------------------
B_scan = np.linspace(0.0, 0.4, 400)
energies_scan = []

for B in B_scan:
    H = hamiltonian(B * n_hat, g_tensor, A_tensor)
    w, _ = eigh(H)
    energies_scan.append(w / (hbar * GHz))  # GHz

energies_scan = np.array(energies_scan)

plt.figure(figsize=(8, 5))
for lvl in range(energies_scan.shape[1]):
    plt.plot(B_scan, energies_scan[:, lvl], lw=1)

# Mark resonance fields
for (_, B_res) in B_res_list:
    plt.scatter(B_res, [freq / GHz] * len(B_res), color='red', zorder=10)

plt.xlabel("B (T)")
plt.ylabel("Energy (GHz)")
plt.title("Spin Hamiltonian levels for S=½, I=5/2 (isotropic g, A)")
plt.grid(True)
plt.tight_layout()
plt.show()
