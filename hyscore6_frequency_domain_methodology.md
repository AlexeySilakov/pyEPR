# Frequency-domain simulation of 6-pulse HYSCORE

How to simulate 6-pulse HYSCORE (hyperfine sublevel correlation) spectra directly in the frequency domain, and what has to change in an existing 4-pulse HYSCORE simulation to turn it into a 6-pulse one.

The method works for any nuclear spin I ≥ ½, including nuclear quadrupole interactions. Every expression below was checked against a full density-matrix simulation that propagates the complete electron–nuclear density matrix through the pulse sequence and selects each pathway by projection. The tests covered I = ½, 1 and 3/2 with random anisotropic hyperfine and quadrupole tensors. The agreement is about 10⁻¹⁵, which is floating-point rounding. See section 8.

---

## 1. Assumptions

- Electron spin S = ½, high-field approximation for the electron (only the electron spin component along the field enters the hyperfine term).
- Ideal pulses: infinitely short, acting only on the electron spin, with exact flip angles.
- No relaxation and no excitation-bandwidth effects.
- One unpaired electron. Any number of nuclei can be included through the product rule (section 6).

Non-ideal pulses, and in particular a non-ideal π mixing pulse that produces intra-manifold (diagonal) peaks, are **not** covered.

---

## 2. Pulse sequences

| Sequence | Pulses and delays |
|---|---|
| 4-pulse HYSCORE | (π/2) – τ – (π/2) – T₁ – (π) – T₂ – (π/2) – τ – echo |
| 6-pulse HYSCORE | (π/2)ᵧ – τ₁ – (π)ᵧ – τ₁ – (π/2)ₓ – T₁ – (π)ₓ – T₂ – (π/2)ₓ – τ₂ – (π)ₓ – τ₂ – echo |

In both sequences T₁ and T₂ are incremented and τ (or τ₁, τ₂) are fixed.

The 6-pulse sequence is the 4-pulse sequence with each fixed τ period replaced by a refocused block: a (π/2 – τ₁ – π – τ₁) block at the start and a (π/2 – τ₂ – π – τ₂) block at the end. It was introduced by Song and co-workers (1995) and analysed by Kasumaj and Stoll (2008).

---

## 3. Spin system and notation

Nuclear Hamiltonians in the two electron manifolds (angular frequency units):

- H_α = ω_I I_z + ½ (A_zx I_x + A_zy I_y + A_zz I_z) + H_Q   (m_S = +½)
- H_β = ω_I I_z − ½ (A_zx I_x + A_zy I_y + A_zz I_z) + H_Q   (m_S = −½)

Here (A_zx, A_zy, A_zz) is the row of the hyperfine tensor along the magnetic field, and H_Q = I·P·I is the nuclear quadrupole term (zero for I = ½). Both are evaluated for the current orientation.

Definitions:

- n = 2I + 1, the number of nuclear sublevels per electron manifold.
- Eigen-decompositions: H_α = V_α diag(ω^α_i) V_α†, H_β = V_β diag(ω^β_i) V_β†.
- Propagators: U_α(t) = exp(−i H_α t), U_β(t) = exp(−i H_β t).
- Nuclear transition frequencies: ω^α_jl = ω^α_j − ω^α_l, and likewise for β. For I ≥ 1 these include double-quantum and higher transitions.

---

## 4. Pathway formalism

With ideal pulses, the echo is a sum of contributions from electron coherence transfer pathways: sequences of electron coherence orders (+1, −1, or population in the α or β manifold) during each free-evolution period.

A pulse only moves the nuclear operator from one electron block to another, multiplied by a scalar electron factor. Free evolution in block (x, y) acts as X → U_x(t) X U_y(t)†. Each pathway therefore reduces to a product of nuclear propagators, multiplied by an electron factor that depends only on the pulse flip angles and phases.

### 4.1 Pathways that contribute

| Sequence | Pathway label | Coherence orders during the free periods | Electron factor (pulse phases as in section 2) |
|---|---|---|---|
| 4-pulse | αβ | (+, 0_α, 0_β, −) | +i/8 |
| 4-pulse | βα | (+, 0_β, 0_α, −) | +i/8 |
| 6-pulse | αβ+ | (+, −, 0_α, 0_β, +, −) | −1/8 |
| 6-pulse | αβ− | (−, +, 0_α, 0_β, +, −) | +1/8 |
| 6-pulse | βα+ | (+, −, 0_β, 0_α, +, −) | −1/8 |
| 6-pulse | βα− | (−, +, 0_β, 0_α, +, −) | +1/8 |

0_α and 0_β are populations or nuclear coherences inside the α and β manifold. The measured echo is the **plain sum** of all contributing pathways, each multiplied by its electron factor.

The absolute factors (1/8, i/8) depend on the pulse-phase and detection conventions. What matters physically is the **relative sign** between the + and − pathways of the 6-pulse sequence. With that sign, the T-independent parts cancel between the + and − pathways for a single nucleus, and the modulation adds. This matches the statement in Kasumaj and Stoll (2008).

### 4.2 General trace expression

Every pathway in both sequences fits one template. Let x be the manifold during T₁ and y the manifold during T₂, so (x, y) = (α, β) or (β, α):

$$E^{xy}(T_1,T_2) = \frac{1}{n}\,\mathrm{tr}\Big[\,B\;U_y(T_2)\,U_x(T_1)\;X\;U_x(-T_1)\,U_y(-T_2)\Big]$$

- **X** is the *preparation operator*: everything before T₁. It depends only on the first fixed delay(s).
- **B** is the *detection operator*: everything after T₂, moved to the front using the cyclic property of the trace. It depends only on the last fixed delay(s).

| Sequence and pathway | Preparation operator X | Detection operator B |
|---|---|---|
| 4-pulse, both pathways | U_α(τ) U_β(−τ) | U_α(−τ) U_β(τ) |
| 6-pulse, + pathways | U_β(τ₁) U_α(τ₁) U_β(−τ₁) U_α(−τ₁) | U_β(−τ₂) U_α(−τ₂) U_β(τ₂) U_α(τ₂) |
| 6-pulse, − pathways | U_α(τ₁) U_β(τ₁) U_α(−τ₁) U_β(−τ₁) | U_β(−τ₂) U_α(−τ₂) U_β(τ₂) U_α(τ₂) |

**Important:** going from αβ to βα swaps **only** the manifolds during T₁ and T₂. X and B stay the same. Exchanging α and β *everywhere* (including inside X and B) gives the right answer for I = ½ only. It fails for I ≥ 1: in our test it was off by about 0.03 in a signal of order 0.1.

---

## 5. Factorization into amplitudes and frequencies

This step is what makes the density-matrix approach as fast as a frequency-domain calculation.

Write U_x(T₁) in the eigenbasis of H_x and U_y(T₂) in the eigenbasis of H_y. Define

- W = V_x† V_y, the overlap matrix between the two eigenbases (equal to M = V_α†V_β for αβ, and to M† for βα),
- X′ = V_x† X V_x, the preparation operator in the T₁ eigenbasis,
- B′ = V_y† B V_y, the detection operator in the T₂ eigenbasis.

The trace then becomes a sum of complex exponentials:

$$E^{xy}(T_1,T_2) = \sum_{j,l,q,t} A_{jl,qt}\;e^{-i\,\omega^{x}_{jl}\,T_1}\;e^{-i\,\omega^{y}_{qt}\,T_2}$$

$$A_{jl,qt} = \frac{1}{n}\;B'_{tq}\;\big(W^\dagger\big)_{qj}\;X'_{jl}\;W_{lt}$$

For each pathway and orientation:

- there are n⁴ complex amplitudes A (16 for I = ½, 81 for I = 1, 256 for I = 3/2);
- each amplitude belongs to a frequency pair (ω^x_jl, ω^y_qt);
- **nothing depends on T₁ or T₂.** The fixed delays enter only through X and B, so the blind-spot behaviour is built into the amplitudes.

These are the same "amplitude at a frequency pair" lists that frequency-domain HYSCORE programs work with. The density-matrix route is simply a general way to obtain them, for any sequence that fits the template in section 4.2.

---

## 6. Algorithm

### 6.1 Single nucleus, one orientation

1. Build H_α and H_β for this orientation.
2. Diagonalize both to get eigenvalues ω^α_i, ω^β_i and eigenvectors V_α, V_β.
3. Build the propagators needed for X and B at the fixed delays: U(t) = V diag(e^{−iωt}) V†.
4. For each pathway:
   1. form X and B (section 4.2 table);
   2. pick (x, y) from the pathway;
   3. compute W, X′, B′ and the amplitudes A_{jl,qt};
   4. compute the frequency pairs (ω^x_jl, ω^y_qt);
   5. multiply the amplitudes by the pathway's electron factor.
5. Add every amplitude into a 2D histogram at its frequency pair (frequency in MHz = ω / 2π).

### 6.2 Powder average

Loop over orientations, using a grid over the sphere with integration weights. Repeat 6.1 for each one, adding into the **same** histogram with the orientation weight. The cost per orientation is 4 × n⁴ amplitudes for the 6-pulse sequence (2 × n⁴ for the 4-pulse), independent of how fine the spectral grid is.

### 6.3 Final spectrum

- Convolve the 2D histogram with a 2D line-shape function (for example Gaussian) for line broadening.
- Bin the real and imaginary parts of the amplitudes separately. For I ≥ 1 the amplitudes have imaginary parts (up to about 0.02 in our I = 1 and I = 3/2 tests). For a magnitude spectrum, combine them after binning: |Re + i·Im|.
- Each amplitude with frequency +ω has a partner at −ω, so the histogram fills both the (+,+) and (+,−) quadrants. Which half is called "positive" depends on the Fourier-transform sign convention. Keep it consistent with the experimental processing.

### 6.4 Time-domain alternative

If you need a time trace (for example for apodization identical to the experiment):

$$E = F_1\,A\,F_2^{\mathsf T},\qquad (F_1)_{a,(jl)} = e^{-i\omega^x_{jl}T_{1,a}},\quad (F_2)_{b,(qt)} = e^{-i\omega^y_{qt}T_{2,b}}$$

This is two matrix multiplications per pathway.

### 6.5 Measured timing (one orientation, I = 1, 256 × 256 grid)

| Method | Time |
|---|---|
| Full density-matrix simulation, point by point | about 87 s |
| Trace expression (section 4.2) evaluated point by point | about 38 s |
| Factored amplitudes → time domain (section 6.4) | about 15 ms |
| Factored amplitudes → frequency histogram (section 6.1) | about 2.4 ms |

These were measured on a cloud machine running Python with NumPy; the first two rows are extrapolated from subsets of the grid. The ratios matter more than the absolute numbers.

### 6.6 Reference implementation (Python)

```python
import numpy as np


def propagator(eigenvalues, eigenvectors, t):
    phases = np.exp(-1j * eigenvalues * t)
    return eigenvectors @ np.diag(phases) @ eigenvectors.conj().T


def pathway_amplitudes(H_alpha, H_beta, X, B, order):
    """Amplitudes (n^2 x n^2) and frequency lists for one pathway.

    order = "ab": alpha manifold during T1, beta during T2; "ba": the reverse.
    X, B: preparation and detection operators (section 4.2), in the laboratory nuclear basis.
    """
    n = H_alpha.shape[0]
    e_alpha, V_alpha = np.linalg.eigh(H_alpha)
    e_beta, V_beta = np.linalg.eigh(H_beta)
    if order == "ab":
        e_first = e_alpha
        V_first = V_alpha
        e_second = e_beta
        V_second = V_beta
    else:
        e_first = e_beta
        V_first = V_beta
        e_second = e_alpha
        V_second = V_alpha
    W = V_first.conj().T @ V_second
    X_eigen = V_first.conj().T @ X @ V_first
    B_eigen = V_second.conj().T @ B @ V_second
    A = np.einsum("tq,jq,jl,lt->jlqt", B_eigen, W.conj(), X_eigen, W) / n
    w1 = (e_first[:, None] - e_first[None, :]).ravel()
    w2 = (e_second[:, None] - e_second[None, :]).ravel()
    return A.reshape(n * n, n * n), w1, w2


def six_pulse_pathways(H_alpha, H_beta, tau1, tau2):
    """List of (amplitudes, w1, w2) for the four 6-pulse pathways, electron factors included."""
    e_alpha, V_alpha = np.linalg.eigh(H_alpha)
    e_beta, V_beta = np.linalg.eigh(H_beta)

    def Ua(t):
        return propagator(e_alpha, V_alpha, t)

    def Ub(t):
        return propagator(e_beta, V_beta, t)

    X_plus = Ub(tau1) @ Ua(tau1) @ Ub(-tau1) @ Ua(-tau1)
    X_minus = Ua(tau1) @ Ub(tau1) @ Ua(-tau1) @ Ub(-tau1)
    B = Ub(-tau2) @ Ua(-tau2) @ Ub(tau2) @ Ua(tau2)
    result = []
    for order in ("ab", "ba"):
        A_plus, w1, w2 = pathway_amplitudes(H_alpha, H_beta, X_plus, B, order)
        result.append((-A_plus / 8, w1, w2))
        A_minus, w1, w2 = pathway_amplitudes(H_alpha, H_beta, X_minus, B, order)
        result.append((A_minus / 8, w1, w2))
    return result


def add_to_histogram(pathways, edges1_MHz, edges2_MHz, spectrum_real, spectrum_imag, weight=1.0):
    for A, w1, w2 in pathways:
        nu1 = np.repeat(w1 / (2 * np.pi), w2.size)
        nu2 = np.tile(w2 / (2 * np.pi), w1.size)
        values = weight * A.ravel()
        real_part, _, _ = np.histogram2d(nu1, nu2, bins=[edges1_MHz, edges2_MHz], weights=values.real)
        imag_part, _, _ = np.histogram2d(nu1, nu2, bins=[edges1_MHz, edges2_MHz], weights=values.imag)
        spectrum_real += real_part
        spectrum_imag += imag_part
```

---

## 7. Several nuclei

With ideal pulses the product rule holds for nuclei of any spin, **per pathway**:

$$E_{6p} = \sum_{p\,\in\,\{\alpha\beta+,\ \alpha\beta-,\ \beta\alpha+,\ \beta\alpha-\}} f_p \prod_{q=1}^{N} E^{p}_{q}$$

Here f_p is the electron factor of pathway p and E^p_q is the single-nucleus nuclear trace for nucleus q. The product must be taken **inside each pathway, before** the pathways are summed.

In the frequency domain, a product in time is a combination of amplitude lists. For each pathway, take every combination of one entry from each nucleus' list, multiply the amplitudes and add the frequencies (in both dimensions). The number of terms is ∏ n_q⁴ per pathway. When that becomes large, convolve the per-nucleus 2D histograms of each pathway, for example with a fast Fourier transform, keeping the complex values and the zero-frequency bins. Then sum the pathways.

Expected physical behaviour (Kasumaj and Stoll 2008): the inter-nuclear combination peaks that appear in 4-pulse HYSCORE cancel between the + and − pathways in 6-pulse HYSCORE, and cross-suppression of weakly modulating nuclei is reduced.

*Status:* the per-pathway product rule follows from the ideal-pulse factorization. The frequency-domain combination step described here has **not** been implemented or tested in this work.

---

## 8. Modifying a 4-pulse HYSCORE simulation into a 6-pulse one

If you already have a frequency-domain 4-pulse simulation built on the amplitude factorization of section 5, most of it carries over unchanged.

### 8.1 What stays the same

- Spin Hamiltonians H_α and H_β, their diagonalization, and the orientation grid with its weights.
- The overlap matrix W = V_x†V_y and the amplitude formula A_{jl,qt} = (1/n) B′_{tq} W†_{qj} X′_{jl} W_{lt}.
- The frequency pairs (ω^x_jl, ω^y_qt): 6-pulse HYSCORE puts peaks at the **same positions** as 4-pulse HYSCORE.
- Histogram binning, line broadening and powder summation.
- The way the βα pathway is obtained from αβ: swap only the T₁/T₂ manifolds.

### 8.2 What changes

| Item | 4-pulse | 6-pulse |
|---|---|---|
| Fixed delays | one, τ | two, τ₁ and τ₂ (set equal unless you are tuning blind spots) |
| Number of pathways | 2 (αβ, βα) | 4 (αβ+, αβ−, βα+, βα−) |
| Preparation operator X | U_α(τ) U_β(−τ) | + pathways: U_β(τ₁) U_α(τ₁) U_β(−τ₁) U_α(−τ₁)<br>− pathways: U_α(τ₁) U_β(τ₁) U_α(−τ₁) U_β(−τ₁) |
| Detection operator B | U_α(−τ) U_β(τ) | U_β(−τ₂) U_α(−τ₂) U_β(τ₂) U_α(τ₂) (same for all four pathways) |
| Electron factors | equal sign for both pathways (+i/8 each, with the phases in section 2) | −1/8 for + pathways, +1/8 for − pathways (relative sign is essential) |
| Amplitudes per orientation | 2 × n⁴ | 4 × n⁴ |
| Blind-spot factor (I = ½) | sin(ω_α τ/2) sin(ω_β τ/2) | sin(ω_α τ₁/2) sin(ω_α τ₂/2) sin(ω_β τ₁/2) sin(ω_β τ₂/2) |
| Multi-nucleus product rule | 2 products, summed | 4 products, combined with the signs of the electron factors |

### 8.3 Step-by-step changes to the code

1. **Inputs:** replace the single τ with τ₁ and τ₂.
2. **Preparation operator:** replace the one-propagator-pair X with the two refocused versions X₊ and X₋. That is four propagators each, built from U_α and U_β at ±τ₁.
3. **Detection operator:** replace B with the four-propagator τ₂ block.
4. **Pathway loop:** loop over (order ∈ {αβ, βα}) × (sign ∈ {+, −}), giving four pathways instead of two. Use X₊ or X₋ according to the sign; B is shared.
5. **Electron factors:** multiply + pathways by −1/8 and − pathways by +1/8 (or by the factors that follow from your own pulse-phase convention, keeping the relative sign).
6. **Summation:** add all four pathways into the histogram. Do **not** subtract the − pathways a second time. Their sign is already in the electron factor.
7. **Multiple nuclei:** form the per-nucleus products within each of the four pathways before summing.

### 8.4 Pitfalls found while developing this

- **Exchanging α and β everywhere for the βα pathways** is only correct for I = ½. For I ≥ 1 only the T₁/T₂ manifolds swap.
- **Double-counting the pathway sign** (applying the "+ − + −" product-rule signs on top of electron factors that already contain them) cancels the modulation. For a single nucleus the result is a flat, T-independent signal.
- **Blind spots:** with τ = 0.136 µs a proton at 14.5 MHz sits on a blind spot (ν_I τ ≈ 2) and its peaks disappear. The 6-pulse blind-spot factor is a product of four sines, so blind spots are sharper. Compute spectra at several (τ₁, τ₂) values when you compare with experiment.
- **Imaginary parts** of the amplitudes are not negligible for I ≥ 1. Keep both parts until the final magnitude step.

---

## 9. Validation performed

| Check | Result |
|---|---|
| 2- and 4-pulse expressions in this notation against full simulation | match to about 10⁻¹⁵ |
| 6-pulse trace expression, all four pathways, I = ½, 1, 3/2, random hyperfine and quadrupole tensors | match to about 10⁻¹⁵, with fixed electron factors ∓1/8 |
| Factored amplitudes (section 5), 4-pulse and 6-pulse, all pathways, I = ½, 1, 3/2 | match to about 10⁻¹⁵ |
| I = ½ result against the closed-form expression of Kasumaj and Stoll (2008), Appendix C | consistent (same simulation reproduces both) |
| Frequency-domain combination for several nuclei (section 7) | **not tested** |

The full simulation used for these checks builds the complete (2 × n)-dimensional electron–nuclear density matrix, applies exact pulse rotations, projects onto each coherence order after every pulse, and propagates with the full Hamiltonian. It shares no code with the factored formula, apart from the Hamiltonian construction.

---

## 10. References

- B. Kasumaj, S. Stoll, "5- and 6-pulse electron spin echo envelope modulation (ESEEM) of multi-nuclear spin systems", *J. Magn. Reson.* **190** (2008) 233–247. Closed-form 6-pulse HYSCORE expression for I = ½ (Appendix C, Equations 50–54), product rule, phase cycles. https://depts.washington.edu/stollgrp/publication/2008_kasumaj_and_stoll/2008_kasumaj_and_stoll.pdf
- R. Song, Y. C. Zhong, C. J. Noble, J. R. Pilbrow, D. R. Hutton, "A new six-pulse two-dimensional electron spin echo envelope modulation (ESEEM) correlation spectroscopy", *Chem. Phys. Lett.* **237** (1995) 86–90. Original 6-pulse sequence.
- P. Höfer, A. Grupp, H. Nebenführ, M. Mehring, *Chem. Phys. Lett.* **132** (1986) 279–282. Original 4-pulse HYSCORE.
- A. Schweiger, G. Jeschke, *Principles of Pulse Electron Paramagnetic Resonance*, Oxford University Press (2001). Density-matrix treatment of ESEEM and HYSCORE, including arbitrary nuclear spin.

Companion scripts from this work: `plot_6pulse_comparison.py` (full simulation against the trace expression, with plots) and `fast_hyscore6.py` (factored amplitudes, time-domain and histogram evaluation, timing).
