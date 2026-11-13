# 4-pulse HYSCORE simulation: needed fix for nuclei with I ≥ 1

**Status:** not fixed yet. The 6-pulse branch (`HYSCOREsim.amplitudes_6p`) is correct for every I. The existing 4-pulse branch is correct for I = ½ only.

## 1. Where

`spin_system/hyscore_sim.py`, `HYSCOREsim.compute_transitions`, the 4-pulse block (the `else` path after the `if seq == '6-pulse'` branch):

```python
ampar = np.array(
         M[ii, ll]*
        cM[ll, jj]*
         M[jj, nn]*
        cM[nn, kk]*
         M[kk, mm]*
        cM[mm, ii]
         )
exppar = np.array(np.exp(-2j*np.pi*(wa[ii, jj] + wb[ll, mm])*self.Exp.tau*1e-3) +
                  np.exp(+2j*np.pi*(wa[kk, jj] + wb[nn, mm])*self.Exp.tau*1e-3))*self.ak[n_ang]
...
tAmp[:, 0] = ampar*exppar            # binned at (wa[ii,kk], wb[ll,nn])
tAmp[:, 1] = np.conj(ampar)*exppar   # binned at (wb[ll,nn], wa[ii,kk])
```

## 2. What is wrong

Write the two exponentials of `exppar` as

- `e1 = exp(-2πi (wa[i,j] + wb[l,m]) τ)`
- `e2 = exp(+2πi (wa[k,j] + wb[n,m]) τ)`

The correct pathway amplitudes, derived from the trace formula of `hyscore6_frequency_domain_methodology.md` §4.2–5 with X = U_α(τ)U_β(−τ) and B = U_α(−τ)U_β(τ), are:

| Pathway | Bin position | Correct amplitude | Existing code |
|---|---|---|---|
| αβ (column 0) | (wa[i,k], wb[l,n]) | `conj(ampar) * e1` | `ampar * (e1 + e2)` |
| βα (column 1) | (wb[l,n], wa[i,k]) | `ampar * e2` | `conj(ampar) * (e1 + e2)` |

So the existing code gives **each** pathway both phase terms, and pairs them with the wrong one of `ampar` / `conj(ampar)`. Term `ampar*e2` is the βα amplitude, but it is binned at the αβ position (and vice versa).

For I = ½ this happens not to matter: `bin_hyscore` adds every amplitude at both (a, b) and (b, a), and after that symmetrisation the two versions are identical to rounding. For I ≥ 1 they are not.

## 3. Evidence

All checks are scratch scripts (not in the repo). The reference is a full (2n)-dimensional density-matrix simulation: ideal pulses, coherence-order projection after every pulse, random anisotropic A and Q tensors. It shares no code with `hyscore_sim.py`.

| Check | I = ½ (¹H) | I = 1 (¹⁴N) | I = 3/2 (¹¹B) |
|---|---|---|---|
| Generic trace formula vs full density matrix, 4-pulse | 3e-16 | 5e-16 | 3e-15 |
| Existing code vs generic, time domain, one orientation | exact (scale 4) | **~40 % off** | – |
| Existing code vs generic, powder spectrum (125 orientations, `Spectrum`) | 8e-16 | **33 % off** (2 parameter sets) | – |
| Corrected expression (table in §2) vs generic, time domain | 1e-15 | 1e-15 | 5e-16 |

No choice of sign convention, α/β mapping or complex conjugation reconciles the existing code with the reference for I = 1. Peak **positions** are right; relative peak **intensities** are wrong. That matters for intensity fitting (grid search) on ¹⁴N, ²H, ¹¹B, ⁶³Cu and similar nuclei.

## 4. The fix (drop-in, keeps the n⁶ index arrays)

```python
ampar = (M[ii, ll]*cM[ll, jj]*M[jj, nn]*cM[nn, kk]*M[kk, mm]*cM[mm, ii])
tau = self.Exp.tau*1e-3
e1 = np.exp(-2j*np.pi*(wa[ii, jj] + wb[ll, mm])*tau)
e2 = np.exp(+2j*np.pi*(wa[kk, jj] + wb[nn, mm])*tau)
scale = 4.0/mult*self.ak[n_ang]            # see "Scale" below
tAmp[:, 0] = scale*np.conj(ampar)*e1       # alpha-beta, at (wa[ii,kk], wb[ll,nn])
tAmp[:, 1] = scale*ampar*e2                # beta-alpha, at (wb[ll,nn], wa[ii,kk])
if self.Exp.tDead:
    tAmp *= np.exp(-2j*np.pi*(wa[ii, kk] + wb[ll, nn])*self.Exp.tDead)[:, None]
```

`tW`, the concatenation and both `bin_hyscore` calls stay unchanged.

**Scale.** The generic formula equals `(1/n) × corrected expression`. The current I = ½ output equals 4 × generic. With `4/mult` the I = ½ spectra stay numerically identical to today's (to be confirmed after the change), and every I uses the same normalisation as `amplitudes_6p` (which also uses `4/n`). That keeps 4-pulse and 6-pulse datasets directly comparable.

**Alternative.** Route 4-pulse through the same generic code as `amplitudes_6p`, with X = U_α(τ)U_β(−τ), B = U_α(−τ)U_β(τ), a single pathway sign +1 and `4/n`. This needs n⁴ instead of n⁶ amplitude terms per nucleus and orientation. It would also make the `cii…cmm` arrays in `build_spin_operators` unnecessary.

## 5. How to verify after fixing

1. ¹H-only systems: `rawSpec` identical to the version before the fix (max relative difference about 1e-15).
2. ¹⁴N with Q, one orientation: `existing_4p_td`-style time-domain sum vs the generic route / full density matrix, max relative error about 1e-15.
3. Re-run any grid-search fits on I ≥ 1 nuclei; their χ² values will change.

## 6. Not addressed here

- `bin_hyscore` symmetrises every amplitude about the diagonal ((a, b) and (b, a)). This is applied identically to both sequences and was left as is.
- The product rule for several nuclei (`Opt.ProdRule`) is still unimplemented for both sequences. Nuclei are summed, so no combination peaks are produced.
