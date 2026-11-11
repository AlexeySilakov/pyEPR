# EPRelax

A wxPython program for regularized exponential analysis of Bruker relaxation
data: it loads a BES3T (`.DSC`/`.DTA`) or ESP (`.PAR`/`.SPC`) relaxation trace
and decomposes it into a distribution of relaxation times, the way CONTIN does for
dynamic light scattering, with a conventional multi-exponential fit alongside.
It uses the same layout and widgets as pyCWEPR and pyHYSCORE.

Start it with:

```
python EPRelax.py
```

Requirements: `pip install wxPython numpy matplotlib scipy`

| File | Contents |
|---|---|
| `EPRelax.py` | The program: `MainFrame` (datasets, analysis, input/output), `FileBrowserPanel`, `PlotPanel`, `ControlPanel`, `DictPopup` |
| `theme/` | Shared with pyCWEPR and pyHYSCORE: pill widgets, `Card`, `PillTabBar`, `IconButton`, spin-box Enter handling, themed navigation toolbar |
| `readers/brukerread.py` | Bruker BES3T/ESP reader (`brukerread`, `parse_dsc`) |
| `icons/` | File-list icons, found relative to the script |

This document is in two parts. Part I is the theory and the algorithms as they
are implemented. Part II is the interface.

---

# Part I — Theory and algorithms

## 1. The model

Every trace is described as a sum of decaying exponentials on a constant
baseline:

$$
S(t) \;=\; c \;+\; s\sum_{j=1}^{n} a_j\,e^{-t/T_j},
\qquad a_j \ge 0,\quad c \ \text{free in sign},\quad s=\pm 1 .
$$

The time constants $T_j$ are **not** fitted. They are a fixed, logarithmically
spaced grid of $n$ points, so the only unknowns are the amplitudes $a_j$ and the
offset $c$, and the problem is **linear**. The amplitude vector $\mathbf a$ over
that grid is the result: the distribution of relaxation times.

$s=-1$ (the **Recovery** switch) turns the same solver into inversion or
saturation recovery, where the exponentials come up to the offset from below.
The amplitudes stay non-negative in either case; only the sign in front of the
sum turns over.

Written as a matrix problem, with $\mathbf y$ the measured trace at the $m$ time
points $t_i$:

$$
A\,\mathbf x \simeq \mathbf y,
\qquad
A=\begin{pmatrix} K & +\mathbf 1 & -\mathbf 1\end{pmatrix},
\qquad
K_{ij}=s\,e^{-t_i/T_j},
$$

$$
\mathbf x=(a_1,\dots,a_n,\,u,\,v)^{\mathsf T},
\qquad c=u-v,
\qquad \mathbf x\ge 0 .
$$

In `MainFrame.analyse`, that is three lines:

```python
K = np.exp(-t[:, None]/T[None, :])
if dd['recovery']:
    K = -K
A = np.hstack([K, np.ones((t.size, 1)), -np.ones((t.size, 1))]) if dd['offset'] else K
```

Non-negativity is required of the amplitudes on physical grounds — a relaxation
component cannot have negative weight — and it is what makes the decomposition
usable at all, being the strongest piece of prior information available. The
offset must stay free in sign, so it enters as the difference of two
non-negative columns, $+\mathbf 1$ and $-\mathbf 1$, rather than as one signed
parameter. Without **Fit offset** those two columns are absent and $c=0$.

Recovering $\mathbf a$ from $\mathbf y$ is the numerical inverse Laplace
transform, and it is ill-posed: the columns of $K$ are nearly parallel for
neighbouring $T_j$, so wildly different distributions fit the same data to
within the noise. The problem is therefore regularized (§5).

## 2. From file to fitted trace

`MainFrame.analyse` builds the trace the decomposition sees, in this order.

### 2.1 Channel and slice

A complex dataset is reduced to one real trace by the **Channel** switch: real,
imaginary or magnitude. No phase correction is applied anywhere, so a mis-set
spectrometer phase stays visible rather than being corrected silently. The
magnitude channel is positive-biased in the noise, which shows up as spurious
long-$T$ amplitude in the tail.

A two-dimensional dataset is reduced to one column either by taking slice $k$
(**Slice**) or by summing all slices (**Sum**), which is only meaningful when
every slice decays the same way. Both choices, at the top of `MainFrame.analyse`:

```python
col = data if data.ndim == 1 else (data.sum(axis=1) if dd['sum']
                                   else data[:, min(dd['slice'], data.shape[1] - 1)])
y = {0: np.real, 1: np.imag, 2: np.abs}[dd['comp']](col).astype(float)
```

### 2.2 Time window

Points outside $[\,t_\mathrm{from},\,t_\mathrm{to}\,]$ are dropped. The window
cuts dead time, ringdown, the rising edge of an echo, and a noisy tail. Both
limits are the dotted lines on the trace row and can be dragged there.

In `MainFrame.analyse`:

```python
keep = (t >= min(dd['win'])) & (t <= max(dd['win']))
t, y = t[keep], y[keep]
if t.size < 6 or not np.all(np.isfinite(y)):
    print(f"[EPRelax] '{dd['title']}' has {t.size} usable points in the window; nothing to decompose")
    dd['res'] = None
    return
```

The decomposition refuses to run on fewer than 6 points in the window, or on a trace
with any non-finite value, and says so on the console.

### 2.3 ESEEM suppression (boxcar)

2-pulse ESEEM carries deep nuclear modulation that no sum of decaying
exponentials can represent. By default the trace is decomposed exactly as recorded
and the modulation is left visible in the residual row. Optionally a boxcar
average of width $w$ is applied first
(`scipy.ndimage.uniform_filter1d`, `mode='nearest'`), over as many points as the
mean sampling step $\overline{\Delta t}$ makes of it:

$$
N_\mathrm{box}=\operatorname{round}\!\left(\frac{w}{\overline{\Delta t}}\right)
\quad\text{points, clamped to } 3\le N_\mathrm{box}\le m .
$$

A boxcar one modulation period wide cancels the modulation. It only rescales an
exponential, and does not shift it:

$$
\frac1w\int_{t-w/2}^{\,t+w/2} e^{-t'/T}\,\mathrm{d}t'
\;=\; e^{-t/T}\,\frac{\sinh\!\left(w/2T\right)}{w/2T},
$$

so the relaxation times come out unchanged while the amplitudes are inflated by
$\sinh(w/2T)/(w/2T)$, which is $\ge 1$ and grows for components faster than the
boxcar. In `MainFrame.analyse`:

```python
if dd['smooth'] and dd['swin'] > 0:
    # a boxcar one modulation period wide cancels the modulation and only rescales an exponential
    width = int(round(dd['swin']/max(abs(float(np.mean(np.diff(t)))), 1e-30)))
    y = uniform_filter1d(y, min(max(width, 3), t.size), mode='nearest')
```

### 2.4 Logarithmic resampling

A linearly sampled trace spends most of its points on the tail, where little is
happening. **Log t** resampling averages the trace into bins evenly spaced in
$\log t$ across the window, after the boxcar. With $\rho$ points per decade, the
bin count and the bin edges are

$$
N_\mathrm{bins}=\max\left(\operatorname{round}\left[\rho\,
\log_{10}\frac{t_\mathrm{last}}{t_\mathrm{first}}\right],\,2\right),
\qquad
\tau_k=t_\mathrm{first}
\left(\frac{t_\mathrm{last}}{t_\mathrm{first}}\right)^{k/N_\mathrm{bins}},
\quad k=0\dots N_\mathrm{bins},
$$

and bin $k$ holds the $m_k$ points with $\tau_k\le t_i<\tau_{k+1}$ (the last bin
includes its upper edge), contributing one point

$$
\bar t_k=\frac1{m_k}\sum_{t_i\in\mathrm{bin}\,k} t_i,
\qquad
\bar y_k=\frac1{m_k}\sum_{t_i\in\mathrm{bin}\,k} y_i .
$$

`MainFrame.analyse` does the binning with `np.bincount`, which both counts the
points per bin and sums them:

```python
decades = float(np.log10(tp[-1]/tp[0]))
edges = np.logspace(np.log10(tp[0]), np.log10(tp[-1]),
                    max(int(round(dd['perdec']*decades)), 2) + 1)
idx = np.clip(np.searchsorted(edges, tp, side='right') - 1, 0, edges.size - 2)
m = np.bincount(idx, minlength=edges.size - 1)
keep = m > 0
tb = np.bincount(idx, weights=tp, minlength=edges.size - 1)[keep]/m[keep]
yb = np.bincount(idx, weights=yp, minlength=edges.size - 1)[keep]/m[keep]
```

A bin holding no points is **dropped**, not interpolated, so the early decay is
never upsampled. Times at or below zero have no logarithm and cannot be binned:
they are kept exactly as they are, one point each, ahead of the bins — a
recovery measured from $t=0$, or a trace with a pre-trigger baseline, keeps those
points. If the result would leave fewer than 6 points altogether, the raw trace
is used instead and the reason is printed.

The density is declared **per decade** rather than as a total count, so it is the
same for every dataset however long the trace ran, and two traces stay
comparable. That is also what lets the card's **Apply to All** button hand the
three settings straight to every loaded dataset, with nothing to rescale on the
way (`ControlPanel.on_resample_all`):

```python
for other in mf.Data:
    other['resample'], other['perdec'], other['rweight'] = (dd['resample'], dd['perdec'],
                                                            dd['rweight'])
```

The settings are part of the kernel signature (§6.6), so every dataset that
changed is decomposed again on the next redraw and nothing else is.

The non-positive times keep a weight of 1 each, and the binned points carry the
precision their count bought:

```python
counts = np.concatenate([np.ones(int(early.sum())), m[keep].astype(float)])
t, y = np.concatenate([t[early], tb]), np.concatenate([y[early], yb])
weight = np.sqrt(counts) if dd['rweight'] == 'count' else None
```

A bin holding $m_k$ points has averaged its noise down by $\sqrt{m_k}$, so each
binned point is given a weight $w_k$:

| Weights | $w_k$ | What it means |
|---|---|---|
| **By count** (default) | $\sqrt{m_k}$, and $1$ for an unbinned non-positive time | The weighted fit is statistically the same question as fitting the raw trace. Resampling then changes the picture and the speed, but not the answer. |
| **Equal** | none, $W=I$ | Every bin counts once, so each decade of time has a say in proportion to its bins rather than to how densely it was sampled. That matches the logarithmic $T$ grid and stops a long flat tail from dominating — but it is a different question, and the answer moves. |

The same weights are handed to the discrete fit (§8), so both analyses answer the
same question.

### 2.5 Amplitude scaling

Before the solve the trace is divided by its own peak,

$$
\varsigma=\max_i\left\lvert y_i\right\rvert,
\qquad
\mathbf y_s=\mathbf y/\varsigma,
$$

so $\alpha$ means the same thing whatever the signal units (`MainFrame.analyse`):

```python
scale = max(float(np.max(np.abs(y))), 1e-300)     # alpha means the same whatever the signal units
ys = y/scale
```

The amplitudes, the offset, the fit and the RMS residual are multiplied back by
$\varsigma$ afterwards.

## 3. The relaxation-time grid

$n$ points (**Points**, shared by every dataset) logarithmically spaced between
$T_\mathrm{min}$ and $T_\mathrm{max}$:

$$
T_j=T_\mathrm{min}
\left(\frac{T_\mathrm{max}}{T_\mathrm{min}}\right)^{\frac{j-1}{n-1}},
\qquad j=1\dots n .
$$

With **Auto** on, the range is read off the trace as it will be fitted — after
windowing and resampling — with $\mathrm{span}=t_\mathrm{last}-t_\mathrm{first}$:

$$
T_\mathrm{min}=\max\!\left(t_\mathrm{first},\;\min\Delta t,\;
10^{-6}\,\mathrm{span}\right),
\qquad
T_\mathrm{max}=3\,\mathrm{span} .
$$

Both, in `MainFrame.analyse`:

```python
step = float(np.min(np.diff(t))) if t.size > 1 else span
dd['Tmin'], dd['Tmax'] = max(float(t[0]), step, span*1e-6), 3*span
lo, hi = sorted((max(dd['Tmin'], 1e-30), max(dd['Tmax'], 1e-29)))
T = np.logspace(np.log10(lo), np.log10(hi), self.npts)
```

The shortest time a trace can argue about is its first point, or its sampling
step where that is coarser. A fixed fraction of the span would do for an evenly
sampled trace, but on a logarithmic time list the span is all last point, and the
whole fast half of the decay would fall off the end of the grid. The upper limit
is three times the span: a relaxation time far longer than the measurement is the
same curve as the fitted offset, and the two then simply trade amplitude.

With **Auto** off the typed limits are used, sorted and floored at $10^{-30}$.
With **Auto** on, the limits actually used are written back into the two boxes, so
they are there to start from when it is switched off.

The regularizer, not the grid, sets how much structure the answer can have, so
more grid points cost time and change little above about 60. The default is 80.

## 4. The penalty matrix

$L$ has $n$ rows and as many columns as $A$, with the offset columns left at
zero — the baseline is never penalized.

**Curvature** (CONTIN's choice, the default) is the second difference along the
log-$T$ grid,

$$
(L\mathbf a)_j=a_{j-1}-2a_j+a_{j+1},
\qquad a_0\equiv a_{n+1}\equiv 0,
\qquad
L=\begin{pmatrix}
-2 & 1 & & & \\
1 & -2 & 1 & & \\
 & \ddots & \ddots & \ddots & \\
 & & 1 & -2 & 1\\
 & & & 1 & -2
\end{pmatrix}.
$$

Both penalties are built in `MainFrame.analyse`, the offset columns never
written to and so left at zero:

```python
L = np.zeros((n, A.shape[1]))
j = np.arange(n)
if self.penalty == 0:            # CONTIN: the curvature along the log-T grid, with the
    L[j, j] = -2.0                # distribution taken as zero outside it, so that the two end
    L[j[:-1], j[:-1] + 1] = 1.0   # points cost what the middle does instead of 2.4 times less
    L[j[1:], j[1:] - 1] = 1.0
else:                             # Tikhonov: the amplitudes themselves
    L[j, j] = 1.0
```

This asks for the smoothest distribution still consistent with the data, and
leaves broad features at full height. The Dirichlet boundary — taking the
distribution as zero outside the grid, which is what puts the $-2$ in the first
and last rows — is deliberate: it prices the two end points like the middle. Left
open, as rows with only two entries, the ends would be about 2.4 times cheaper
and would collect a spike.

**Amplitude** (plain Tikhonov) is $L=I$ on the amplitude block, penalizing the
amplitudes themselves. Very stable, but it shrinks a broad distribution instead
of smoothing it.

The penalty is shared by every dataset, because comparing one distribution with
another only means something when both were made the same way.

## 5. The regularized non-negative solve

### 5.1 Augmented non-negative least squares

For a given $\alpha$, `MainFrame.iLT_decompose` solves

$$
\min_{\mathbf x\ge 0}\;
\left\lVert W\left(A\mathbf x-\mathbf y_s\right)\right\rVert^2
+\alpha^2\left\lVert L\mathbf x\right\rVert^2 ,
$$

where $W=\operatorname{diag}(\mathbf w)$ are the resampling weights, absent
unless **By count** is in use. This is written as one ordinary non-negative
least-squares problem by stacking the penalty under the weighted kernel,

$$
\min_{\mathbf x\ge 0}\;
\left\lVert
\begin{pmatrix} W A\\[2pt] \alpha L\end{pmatrix}\mathbf x
-\begin{pmatrix} W\mathbf y_s\\[2pt] \mathbf 0\end{pmatrix}
\right\rVert^2 ,
$$

and handed to `scipy.optimize.nnls`. Non-negativity is therefore imposed exactly,
not by clipping and not through a barrier. The whole solve is four lines of
`MainFrame.iLT_decompose`:

```python
Aw, yw = (A, y) if w is None else (A*w[:, None], y*w)
x = nnls(np.vstack([Aw, alpha*L]), np.concatenate([yw, np.zeros(L.shape[0])]))[0]
fit = A @ x
res = (y - fit) if w is None else (y - fit)*w
rss = float(np.sum(res*res))
eta = float(np.linalg.norm(L @ x))
c = float(x[n] - x[n + 1]) if x.size > n else 0.0
```

The reported quantities are

$$
\hat{\mathbf y}=A\mathbf x,
\qquad
\mathrm{rss}=\left\lVert W\left(\mathbf y_s-\hat{\mathbf y}\right)\right\rVert^2,
\qquad
\rho=\sqrt{\mathrm{rss}},
\qquad
\eta=\left\lVert L\mathbf x\right\rVert,
$$

$$
\mathrm{RMS}=\varsigma\,
\sqrt{\frac1m\sum_{i=1}^{m}\left(y_{s,i}-\hat y_i\right)^2} .
$$

The fit $\hat{\mathbf y}$ is unweighted, and $\mathrm{RMS}$ is deliberately the
unweighted scatter of the points that are drawn, not the weighted criterion that
was minimized.

### 5.2 Effective degrees of freedom

The number of parameters the data actually paid for is the trace of the hat
matrix, restricted to the amplitudes NNLS kept — the **active set**
$S=\{\,j: x_j>0\,\}$, which includes whichever offset column is active:

$$
\nu(\alpha)=\operatorname{tr}\left[
\left(A_S^{\mathsf T}W^2A_S+\alpha^2L_S^{\mathsf T}L_S\right)^{-1}
A_S^{\mathsf T}W^2A_S\right].
$$

The rest of `MainFrame.iLT_decompose`:

```python
active = x > 0
df = float(active.sum())
if active.any():
    As, Ls = Aw[:, active], L[:, active]
    M = As.T @ As
    try:
        df = float(np.trace(np.linalg.solve(M + alpha*alpha*(Ls.T @ Ls), M)))
    except np.linalg.LinAlgError as fault:
        print(f"[EPRelax] the effective d.o.f. at alpha={alpha:g} fell back to the active-set "
              f"count: {fault}")
```

It runs from $\lvert S\rvert$ at $\alpha\to 0$ down towards 1 at large $\alpha$,
where the penalty has flattened the distribution to a single effective component.
If the linear solve is singular, $\lvert S\rvert$ is used instead and the
fallback is printed to the console.

$\nu$ is what the GCV score and the F-test are built on, and what the automatic
$\alpha$ bracket (§6.1) climbs until it reaches 1.

## 6. Choosing α

Four rules, shared by every dataset so that all of them are regularized the same
way. Three of them need a sweep over $\alpha$.

### 6.1 The sweep and its bracket

33 logarithmically spaced values of $\alpha$ are solved, and $\mathrm{rss}$,
$\eta$ and $\nu$ are kept at each. The sweep runs when $\alpha$ comes from any
automatic rule, and also whenever the L-curve row is drawn — the L-curve *is*
that sweep.

The corner of the L-curve moves by **decades** between traces, with the noise,
the point count and how many components there are, so a fixed ladder can leave
the whole informative part of the curve outside the window. That shows up in the
interface as an $\alpha$ chosen at an end of the L-curve with the residual still
at its floor. With **Auto range** on, the default, the bracket is therefore found
from each trace, by probing one $\alpha$ per decade,
$\alpha^\mathrm{p}_k=10^{k}$ for $k=-12\dots 12$:

$$
\alpha_\mathrm{top}=\min\left\{\alpha^\mathrm{p}_k:\
\nu\!\left(\alpha^\mathrm{p}_k\right)\le 1.01\right\},
\qquad
\mathcal{F}=\left\{k:\
\mathrm{rss}\!\left(\alpha^\mathrm{p}_k\right)\le
1.001\,\mathrm{rss}\!\left(\alpha^\mathrm{p}_{-12}\right)\right\},
$$

$$
\alpha_\mathrm{lo}=10^{\max(\mathcal F)-1},
\qquad
\alpha_\mathrm{hi}=\max\left(\alpha_\mathrm{top},\,10\,\alpha_\mathrm{lo}\right),
\qquad
\alpha_k=\alpha_\mathrm{lo}
\left(\frac{\alpha_\mathrm{hi}}{\alpha_\mathrm{lo}}\right)^{k/32},
\quad k=0\dots 32 .
$$

In `MainFrame.analyse`:

```python
probe = np.logspace(-12, 12, 25)
rss_p, df_p, top = [], [], probe[-1]
for al in probe:
    out = self.iLT_decompose(A, L, ys, al, n, weight)
    rss_p.append(out[3])
    df_p.append(out[5])
    if df_p[-1] <= 1.01:      # one effective component left: nothing beyond here
        top = al
        break
rss_p = np.array(rss_p)
plateau = np.flatnonzero(rss_p <= 1.001*rss_p[0])
lo = float(probe[max(int(plateau[-1]) - 1, 0)]) if plateau.size else float(probe[0])
hi = max(float(top), lo*10)
```

and then the sweep itself, whose three columns are everything the rules below
need:

```python
alphas = np.logspace(np.log10(lo), np.log10(hi), SWEEP)
sweep = np.array([self.iLT_decompose(A, L, ys, al, n, weight)[3:] for al in alphas])
rss_s, eta_s, df_s = sweep[:, 0], np.maximum(sweep[:, 1], 1e-300), sweep[:, 2]
dd['sweep'] = (ksig, rss_s, eta_s, df_s, alphas)
```

The climb stops at $\alpha_\mathrm{top}$, where the distribution has collapsed:
one effective component is left, and nothing can happen beyond it. The bottom is
a decade below the last point of the plateau $\mathcal F$, where the solve stops
changing at all. The 33 swept points then go between. Bracketing costs about 20
extra solves per kernel change, which is cached (§6.6).

With **Auto range** off, $\alpha_\mathrm{lo}$ and $\alpha_\mathrm{hi}$ are the
fixed range typed in the interface, by default $10^{-5}\dots10^{5}$, the same for
every dataset, which is what makes two L-curves directly comparable.

### 6.2 Manual

$\alpha=10^{\,\mathrm{logalpha}}$, from the knob, which belongs to the current
dataset. No sweep is needed unless the L-curve row is being drawn. The one line
of `MainFrame.analyse` every other mode then overwrites:

```python
alpha = 10.0**dd['logalpha']
```

### 6.3 L-curve (Hansen's corner)

The L-curve is $\eta$ against $\rho$ on log-log axes; its corner is the point of
maximum curvature. With $p=\ln\rho$ and $q=\ln\eta$ over the swept points, and
derivatives from `np.gradient`,

$$
\kappa=\frac{p'q''-q'p''}{\left(p'^2+q'^2\right)^{3/2}},
\qquad
\alpha=\arg\max_{\alpha_k}\ \kappa ,
$$

In `MainFrame.analyse`, the `amode == 1` branch:

```python
p, q = np.log(rho), np.log(eta)
move = np.concatenate([[True], (np.diff(p)**2 + np.diff(q)**2) > 1e-12])
if move.sum() < 5:
    print(f"[EPRelax] the L-curve of '{dd['title']}' barely moves over the sweep, so it "
          f"has no corner to find; alpha stays at the knob")
else:
    p, q = p[move], q[move]
    d1p, d1q = np.gradient(p), np.gradient(q)
    d2p, d2q = np.gradient(d1p), np.gradient(d1q)
    kappa = (d1p*d2q - d1q*d2p)/np.power(d1p**2 + d1q**2, 1.5)
    alpha = float(alphas[move][int(np.nanargmax(kappa[1:-1])) + 1])
```

The two endpoints are excluded by the `[1:-1]`. Only the stretch of the curve
that actually moves is used: points with
$\left(\Delta p\right)^2+\left(\Delta q\right)^2\le10^{-12}$ are dropped, because
below some $\alpha$ the solve stops changing at all and the spacing there is
numerical noise, which would otherwise own the curvature. If fewer than 5 points
move, the curve has no corner to find: $\alpha$ stays at the knob and the reason
is printed.

### 6.4 Generalized cross-validation

$$
V(\alpha)=\frac{m\,\mathrm{rss}(\alpha)}{\left[m-\nu(\alpha)\right]^2}
\;=\;\frac{\mathrm{rss}(\alpha)/m}{\left[1-\nu(\alpha)/m\right]^2},
\qquad
\alpha=\arg\min_{\alpha_k}\ V ,
$$

the standard GCV score, with $m$ the number of time points. One line of
`MainFrame.analyse`:

```python
alpha = float(alphas[int(np.argmin(m*rss_s/np.maximum(m - df_s, 1e-6)**2))])
```

### 6.5 F-test (CONTIN's own rule)

The largest $\alpha$ whose fit is still statistically indistinguishable from the
unregularized one. With $\mathrm{rss}_0$ and $\nu_0$ taken at the smallest swept
$\alpha$,

$$
\Delta\nu=\max\left(\nu_0-\nu(\alpha),\,1\right),
\qquad
F(\alpha)=\frac{\left[\mathrm{rss}(\alpha)-\mathrm{rss}_0\right]/\Delta\nu}
{\mathrm{rss}_0/\left(m-\nu_0\right)},
$$

$$
\alpha=\max\left\{\alpha_k:\
F(\alpha_k)\le F^{-1}\!\left(P;\,\Delta\nu,\,m-\nu_0\right)
\ \text{ and }\ \nu(\alpha_k)>0\right\},
$$

where $F^{-1}(P;\cdot,\cdot)$ is the F quantile at probability $P$, which
`MainFrame.analyse` takes from `scipy.stats.f`:

```python
rss0, df0 = rss_s[0], df_s[0]
dnu = np.maximum(df0 - df_s, 1.0)
fstat = ((rss_s - rss0)/dnu)/max(rss0/max(m - df0, 1e-6), 1e-300)
crit = fdist.ppf(self.prob, dnu, max(m - df0, 1.0))
passed = np.flatnonzero((fstat <= crit) & (df_s > 0))
alpha = float(alphas[passed[-1]]) if passed.size else float(alphas[0])
```

$\Delta\nu$ is
floored at 1 because regularizing never buys back less than one degree of
freedom, and an unfloored $\Delta\nu$ near zero makes $F$ diverge. $P$ — the
**Probability** box, shared — defaults to CONTIN's 0.5: the chosen solution is
the smoothest one that a coin toss would not call a worse fit than the
unregularized answer. Raising $P$ allows more smoothing. If nothing passes, the
smallest swept $\alpha$ is used.

This is the default rule.

### 6.6 Caching

Two signatures guard the work, at the head of `MainFrame.analyse`:

```python
ksig = (dd['comp'], dd['slice'], dd['sum'], tuple(dd['win']), dd['recovery'], dd['offset'],
        dd['smooth'], dd['swin'], dd['autogrid'], self.npts, self.penalty,
        dd['resample'], dd['perdec'], dd['rweight'],
        self.autoalpha, self.alpha_lo, self.alpha_hi,
        None if dd['autogrid'] else (dd['Tmin'], dd['Tmax']))
sig = ksig + (self.amode, dd['logalpha'], self.prob)
if dd['sig'] == sig and dd['res'] is not None and not (self.plot_panel.show['lcurve']
                                                       and dd['res']['rho'] is None):
    return
```


* the **kernel signature** — channel, slice, sum, window, recovery, offset,
  boxcar, resampling, grid, point count, penalty and the sweep range. The sweep
  (33 solves, or about 53 with bracketing) is cached against it.
* the **full signature** — the kernel signature plus the $\alpha$ mode, the
  manual $\alpha$ and the F-test probability. If it is unchanged and a result
  exists, `analyse` returns at once.

So turning the manual $\alpha$ knob with the L-curve on screen costs **one**
solve, not a new sweep: the marker moves along the already computed curve.
Drawing the L-curve is not a change to the trace, so it is no part of the
signature; it only asks for the sweep, which a dataset decomposed without it has
not got yet, and that case is checked for separately.

The marker drawn on the L-curve is $(\rho,\eta)$ of the solve actually being
displayed, not the nearest swept point, so it moves continuously as the knob
turns.

## 7. What one decomposition leaves behind

The last act of `MainFrame.analyse` is to leave everything the plots and the
discrete fit need in one dictionary, with the amplitudes, the offset and the fit
scaled back into signal units:

```python
a, c, fit, rss, eta_at, df = self.iLT_decompose(A, L, ys, alpha, n, weight)
...
dd['res'] = {'t': t, 'y': y, 'T': T, 'a': a*scale, 'c': c*scale, 'fit': fit*scale,
             # the scatter of the points that are drawn, not the weighted criterion
             'rms': float(np.sqrt(np.mean((ys - fit)**2)))*scale, 'df': df, 'alpha': alpha,
             'rho': rho, 'eta': eta,
             # where this solve itself sits on that curve, which is not a swept point
             'point': (float(np.sqrt(max(rss, 1e-300))), eta_at),
             # kept so the discrete fit can put a clean exponential through the same decomposition
             'A': A, 'L': L, 'scale': scale, 'w': weight}
```

So: the prepared trace, the grid $T_j$, the amplitudes $a_j$ and the offset $c$,
the fit, the RMS residual, $\nu$, the chosen $\alpha$, the swept $\rho$ and
$\eta$ (or `None` when no sweep was needed), the $(\rho,\eta)$ point of this
particular solve, and the kernel $A$, the penalty $L$, the scale factor
$\varsigma$ and the weights $\mathbf w$ — the last four kept so the discrete fit
can put a clean exponential through the same decomposition.

If every amplitude comes out at zero, so that only the offset is left, the trace
runs the other way than the kernel; the console says so and suggests the
**Recovery** switch or another channel.

## 8. The discrete fit

A conventional multi-exponential fit sits alongside the decomposition and is started
from it, so both analyses are asked the same question.

### 8.1 Starting values

The peaks of the distribution are its interior local maxima,

$$
\mathcal P=\left\{\,j:\ a_j>a_{j-1}\ \text{ and }\ a_j\ge a_{j+1}\,\right\},
\qquad\text{sorted by descending } a_j .
$$

the same list `update_stats` reports. In `MainFrame.discrete`:

```python
peaks = [j for j in range(1, T.size - 1) if a[j] > a[j - 1] and a[j] >= a[j + 1]]
peaks.sort(key=lambda j: -a[j])
...
T0 = [float(T[j]) for j in peaks[:n]]
T0 += list(np.logspace(np.log10(T[0]), np.log10(T[-1]), n + 2)[1:-1])[len(T0):]
amp = max(float(a.sum()), float(np.max(np.abs(y))))/n
```

The $N$ tallest become the $N$ starting time constants, **tallest first**, so the
tallest peak becomes component 1. If more components are asked for than the
distribution has peaks, the remainder start from an evenly log-spaced ladder
across the grid, where they rarely stay put. Every amplitude starts at

$$
a^{(0)}=\frac1N\max\left(\sum_j a_j,\ \max_i\lvert y_i\rvert\right),
$$

and the offset starts at the decomposition's $c$.

### 8.2 The fit

`scipy.optimize.curve_fit` — trust-region reflective, because there are
bounds — on

$$
y(t)=c+s\sum_{i=1}^{N}a_i\,e^{-t/T_i},
\qquad
a_i\ge 0,
\qquad
0.1\,T_\mathrm{min}\le T_i\le 10\,T_\mathrm{max},
\qquad
c\in\mathbb R .
$$

The decade of slack either side of the grid is as far as the trace can argue
about $T$ at all. `sigma` $=1/\mathbf w$ is passed when the trace was resampled
by count, so the discrete fit is weighted exactly as the decomposition was. The free
branch of `MainFrame.discrete`:

```python
model = lambda tt, *p: ((p[-1] if free else 0.0)
                        + s*sum(p[i]*np.exp(-tt/p[n + i]) for i in range(n)))
p0 = [amp]*n + sorted(T0) + ([float(r['c'])] if free else [])
# a decade of slack either side of the grid: outside it the trace cannot argue about T
bounds = ([0.0]*n + [float(T[0])*0.1]*n + ([-np.inf] if free else []),
          [np.inf]*n + [float(T[-1])*10]*n + ([np.inf] if free else []))
popt, pcov = curve_fit(model, t, y, p0=p0, bounds=bounds, maxfev=40000,
                       sigma=None if r['w'] is None else 1.0/r['w'])
```

### 8.3 The 95 % intervals on T and a

`curve_fit` already scales the covariance $\mathbf C$ by the residual variance,
so it is the linearized covariance at the optimum and, with $p$ fitted
parameters,

$$
\nu_\mathrm{fit}=\max\left(m-p,\,1\right),
\qquad
\Delta_{95}=t_{0.975,\,\nu_\mathrm{fit}}\sqrt{\operatorname{diag}\mathbf C},
$$

with $t_{0.975,\nu}$ the two-sided Student-$t$ quantile, from `scipy.stats.t`
(`MainFrame.discrete`):

```python
dof = max(t.size - popt.size, 1)
err = np.sqrt(np.abs(np.diag(pcov)))*float(tdist.ppf(0.975, dof))
```

This is the fit's own trust region: how tightly the trace pins that component down. It goes wide when
two components are not really separable.

### 8.4 Held time constants

With **Hold T at the peaks** the time constants are not fitted. Each is held at
the amplitude-weighted geometric mean of its distribution peak (§9), only the
amplitudes and the offset are fitted, and the number of components is capped at
$\lvert\mathcal P\rvert$ — there is nothing to hold the extras at.

A held component's peak width is then measured **at the peak it was held at**,
not at the centre of mass that came out of it: an asymmetric peak puts that
centre a grid point or two off its own maximum, which would otherwise read as
"not on a peak" (§9).

```python
lookup = np.array([float(T[j]) for j in peaks[:n]]) if fixed else Tout
```

A held $T$ is not a fitted parameter, so the fit gives it no interval. The one
reported, marked `~`, is the interval it *would* have had if it had been free and
that value were the optimum. It comes from the Jacobian of the full model at the
held point,

$$
\frac{\partial y}{\partial a_i}=s\,e^{-t/T_i},
\qquad
\frac{\partial y}{\partial T_i}=s\,a_i\,\frac{t}{T_i^{2}}\,e^{-t/T_i},
\qquad
\frac{\partial y}{\partial c}=1,
$$

with the amplitudes free to move with $T$, so the correlation between a
component's time constant and its amplitude is in there. The columns are weighted
as the fit was, $J\leftarrow WJ$, and then

$$
\mathbf C=\left(J^{\mathsf T}J\right)^{-1}
\frac{\sum_i r_i^2}{\nu_\mathrm{asif}},
\qquad
\nu_\mathrm{asif}=\max\left(m-\operatorname{cols}J,\,1\right),
\qquad
\Delta_{95}(T_i)=t_{0.975,\,\nu_\mathrm{asif}}\sqrt{C_{ii}} .
$$

The held branch of `MainFrame.discrete` builds exactly that:

```python
Tout = Tfix
amps = np.asarray(popt[:n])
columns = [s*np.exp(-t/Tfix[i]) for i in range(n)]
columns += [s*amps[i]*(t/Tfix[i]**2)*np.exp(-t/Tfix[i]) for i in range(n)]
if free:
    columns.append(np.ones_like(t))
jac = np.column_stack(columns)
if r['w'] is not None:
    jac = jac*r['w'][:, None]
asif = max(t.size - jac.shape[1], 1)
try:
    res = (y - fit) if r['w'] is None else (y - fit)*r['w']
    cov = np.linalg.inv(jac.T @ jac)*float(np.sum(res*res))/asif
    dTout = np.sqrt(np.abs(np.diag(cov)[n:2*n]))*float(tdist.ppf(0.975, asif))
except np.linalg.LinAlgError as fault:
    print(f"[EPRelax] the held time constants of '{dd['title']}' have no interval to "
          f"estimate -- the components are not separable there: {fault}")
    dTout = np.full(n, np.inf)
```

It says how tightly the trace pins $T$ down *around that value* — not whether the
value is the right one, which is the question holding it begged in the first
place. If $J^{\mathsf T}J$ is singular the components are not separable there,
the interval is reported as infinite, and the console says so.

### 8.5 Ordering and staleness

Components are sorted by ascending $T$ for display, after the fit. The fit stores
the trace signature it was made on; when the window is dragged afterwards, the
drawn curve still follows the parameters but the numbers belong to the old
window, and the interface says so and asks for a refit. The curve is always
re-evaluated on the current time axis rather than stored, because a curve sampled
on the old window is both the wrong length and in the wrong place. That is all
`MainFrame.discrete_curve` is:

```python
d = dd['disc']
sign = -1.0 if dd['recovery'] else 1.0
return d['c'] + sign*sum(a*np.exp(-t/T) for a, T in zip(d['a'], d['T']))
```

## 9. The width of a distribution peak — `iLT low` and `iLT high`

These are **not** error bars on $T$. They are the width of the distribution peak
the component sits in: a regularized decomposition returns a hump even for a
perfectly monoexponential sample, because the penalty smooths, and this says how
wide that hump is.

`MainFrame.peakwidth` works in $x=\log_{10}T$:

1. Find the grid point nearest the component's $T$.
2. If that point is not a local maximum — the amplitude rises on either side of
   it — the component is on the flank of a peak or in the valley between two.
   There is no peak there to have a width, so none is reported: $\sigma$ comes
   back as NaN, the two cells are left empty and the note says `not on a peak`.
   This is not rare. A fitted $T$ is free to land anywhere, and over an $\alpha$
   sweep of the test traces it misses the maxima about half the time.
3. Walk out to the minima either side, stopping at the first zero amplitude.
   That stretch, $j\in[\,j_\mathrm{lo},j_\mathrm{hi}\,]$, is the peak.
4. Fit a Gaussian in $x$ to the amplitudes over the peak plus two further points
   either side — a narrow peak is only a couple of grid points wide, and three
   parameters need somewhere to pin the shoulders down:

   $$
   a(x)\simeq A\exp\left[-\frac12\left(\frac{x-x_0}{\sigma}\right)^2\right],
   \qquad
   \frac{\Delta x}{2}\le\sigma\le x_n-x_1 .
   $$

   where $\Delta x=x_{j+1}-x_j$ is the grid step in $\log_{10}T$. That lower
   bound of half a grid step matters: a peak narrower than the grid is
   simply unresolved, and without it the optimizer chases $\sigma\to0$ and never
   converges.

   The fit is **done in normalized coordinates** — the abscissa in grid steps
   away from the peak, the ordinate in units of the peak's own height,

   $$
   u=\frac{x-x_j}{\Delta x},
   \qquad
   v=\frac{a}{\max\lvert a\rvert},
   \qquad
   \sigma=\Delta x\,\sigma_u .
   $$

   so that all three parameters are of order one whatever the signal units and
   however fine the grid. In the raw units they differ by many orders of
   magnitude, the trust region is set by the largest of them, and the optimizer
   crawls along a bound until it runs out of evaluations — which is what the
   console message about a Gaussian that "did not converge" used to be.

5. The reported ends are the two-sided 95 % points of that Gaussian, as a
   **factor**:

   $$
   f=10^{\,1.959963985\,\sigma},
   \qquad
   T_\mathrm{low}=\frac{T}{f},
   \qquad
   T_\mathrm{high}=T\,f .
   $$

Steps 2 to 4 are the body of `MainFrame.peakwidth`:

```python
if (j > 0 and a[j] < a[j - 1]) or (j < T.size - 1 and a[j] < a[j + 1]):
    return float('nan'), float(Td), "not on a peak"
lo = hi = j                       # out to the minima either side, stopping at the first zero
while lo > 0 and 0 < a[lo - 1] <= a[lo]:
    lo -= 1
while hi < T.size - 1 and 0 < a[hi + 1] <= a[hi]:
    hi += 1
...
seg = slice(max(lo - 2, 0), min(hi + 3, T.size))
step = float(x[1] - x[0])
height = max(float(np.max(np.abs(a[seg]))), 1e-300)
u, v = (x[seg] - x[j])/step, a[seg]/height
guess = [max(float(a[lo:hi + 1].max()), 1e-300)/height, 0.0,
         max((x[hi] - x[lo])/4, step)/step]
limits = ([0.0, float(u[0]), 0.5], [np.inf, float(u[-1]), float(x[-1] - x[0])/step])
popt = curve_fit(lambda uu, A, u0, w: A*np.exp(-0.5*((uu - u0)/w)**2),
                 u, v, p0=np.clip(guess, *limits), bounds=limits,
                 maxfev=20000)[0]
sigma = abs(float(popt[2]))*step
```

and step 5 is one line of `MainFrame.discrete`, with
`TWO_SIDED_95 = 1.959963985` a module constant:

```python
wide.append(10.0**(TWO_SIDED_95*sigma))
```

They are not symmetric about $T$, because a width in $\log T$ is a factor rather
than a difference.

The **centre** `peakwidth` also returns — used as the held value in §8.4 — is the
amplitude-weighted mean of $\log_{10}T$ across the peak, read back as a time:

$$
\log_{10}T_\mathrm{centre}=
\frac{\displaystyle\sum_{j=j_\mathrm{lo}}^{j_\mathrm{hi}} a_j\,x_j}
     {\displaystyle\sum_{j=j_\mathrm{lo}}^{j_\mathrm{hi}} a_j}\, ,
$$

the geometric mean, which is what a logarithmic grid makes of an average, and the
same quantity the width is a spread about. The last two lines of
`MainFrame.peakwidth`:

```python
weight = a[lo:hi + 1]
centre = 10.0**float((x[lo:hi + 1]*weight).sum()/max(float(weight.sum()), 1e-300))
```

An answer always comes back; the `note` row says when it should not be trusted:

| Note | Meaning |
|---|---|
| `clipped by the grid` | The peak runs into an end of the $T$ grid, so its width is a lower bound |
| `k points` | Fewer than 5 grid points under the peak |
| `no convergence` | The Gaussian fit failed; the starting guess, $\left(x_{j_\mathrm{hi}}-x_{j_\mathrm{lo}}\right)/4$, is reported instead |
| `not on a peak` | The component sits on a flank or in a valley; the two cells are empty, since there is no peak there to measure |
| `no peak here` | The amplitude at that $T$ is zero |
| `T held` | Added when the time constant was held rather than fitted |

In the rate view the two ends swap, since the slow end of the peak is the small
rate.

## 10. Derived numbers and display normalizations

The statistics block under the Regularization card reports $\alpha$ and $\nu$ out
of $n$ grid points, the RMS residual both absolutely and as a percentage of
$\max\lvert y\rvert$, the offset $c$, the total amplitude $\sum_j a_j$, the four
tallest peaks, and the amplitude-weighted mean

$$
\langle T\rangle=\frac{\sum_j a_j T_j}{\sum_j a_j},
$$

which is arithmetic in $T$, unlike the per-peak centre of §9.

The distribution's vertical axis has three settings:

| Setting | Plotted quantity |
|---|---|
| **Raw** | $a_j$, the fitted amplitudes in the signal's own units, so the height is physical and the area under the curve is the total decaying signal |
| **Area** | $P_j=\dfrac{a_j}{\left(\sum_k a_k\right)\Delta\ln T}$ with $\Delta\ln T=\dfrac{\ln\left(T_n/T_1\right)}{n-1}$: a density in $\ln T$, so that $\sum_j P_j\,\Delta\ln T=1$, which is what compares shapes between samples of different concentration. The axis is labelled "Density per unit ln T" — the curve runs well above 1, being a density and not a share |
| **Fraction** | $a_j\big/\sum_k a_k$, so the points add up to 1 and each reads as its share of the decaying signal. Not a density, so its height depends on $n$ |

All three are worked out in `PlotPanel.update_graph`, before anything is drawn,
because the **Separation** below is measured in units of the tallest of them:

```python
total = max(float(r['a'].sum()), 1e-300)
if self.norm == 'area':
    amps[ii] = r['a']/(total*float(np.log(r['T'][-1]/r['T'][0]))/max(r['T'].size - 1, 1))
elif self.norm == 'fraction':
    amps[ii] = r['a']/total
else:
    amps[ii] = r['a']
```

A density in $\ln T$ is unchanged by $T\to1/T$, since $\ln(1/T)=-\ln T$, so the
**Rate** view simply mirrors the curve against $1/T$ and rescales nothing.

**Separation** fans the distributions out vertically: the $k$-th shown trace is
lifted by $k\,\sigma_\mathrm{sep}\max_{i,j}P_j^{(i)}$, a multiple of the tallest
amplitude of them all. With a gap set, the vertical axis no longer reads as
amplitude for anything but the first trace.

The experimental trace has its own two settings:

| Setting | Plotted quantity |
|---|---|
| **Raw** | the signal as recorded |
| **0 to 1** | $\tilde y=\dfrac{y-\min y}{\max y-\min y}$ over the dataset's own window, so traces of different intensity lie on top of each other and a decay and a recovery share one scale |

One map, built per dataset in `PlotPanel.update_graph` and applied to the trace,
both fits and the stacked components alike:

```python
lo = float(r['y'].min()) if self.snorm == 'span' else 0.0
gain = 1.0/max(float(r['y'].max()) - lo, 1e-300) if self.snorm == 'span' else 1.0
scale = lambda v: (v - lo)*gain
```

The residual is multiplied by the gain $1/(\max y-\min y)$ only, so it stays
centred on zero. A
partial inversion is no longer visible as such under **0 to 1**, since every
trace is stretched to full height.

**Components** draws the discrete components stacked under the trace, slowest at
the bottom, starting at the fitted offset. With the components ordered by
descending $T$,

$$
\ell_0(t)=c,
\qquad
\ell_k(t)=\ell_{k-1}(t)+s\,a_{(k)}\,e^{-t/T_{(k)}},
$$

so each curve holds its level until its own time constant comes and then peels
away: the gap between two curves is what the component between them is worth, and
the topmost curve, $\ell_N$, is the whole fit. The horizontal line under each
curve is the plateau $c+s\sum_{i\le k}a_{(i)}$ it starts from. In
`PlotPanel.update_graph`:

```python
sign = -1.0 if dd['recovery'] else 1.0
level = np.full(r['t'].shape, d['c'])
plateau = d['c']
for k, i in enumerate(np.argsort(d['T'])[::-1]):
    top = level + sign*d['a'][i]*np.exp(-r['t']/d['T'][i])
    plateau += sign*d['a'][i]
    ...
    level = top
```

## 11. Units

The time axis is read from the file's own $x$ label and spelled onto a ladder:
ps, ns, µs, ms, s. Overlaid traces share an axis, so they must share a unit:
`MainFrame.harmonize_units` picks the unit for which the longest time of all
loaded datasets stays under 1000,

```python
seconds = max(float(np.max(np.abs(dd['t'])))*ladder[dd['tunit']] for dd in known)
name, factor = TIME_UNITS[-1]
for name, factor in TIME_UNITS:
    if seconds < 1000*factor:
        break
```

and rescales every dataset's time axis, window, grid limits,
boxcar width and fitted time constants onto it, invalidating the cached
decompositions. A dataset whose file names a unit off the ladder keeps its own numbers
and takes no part.

The combined discrete-fit CSV writes every time in **seconds** regardless, so the
column sorts and plots across samples, and names the unit each dataset was shown
in in a column of its own.

## 12. What is assumed, and what is not done

* The signal is a sum of exponentials on a constant baseline. Stretched
  exponentials, spectral diffusion with a time-dependent rate, and oscillatory
  components are not in the model — they appear in the residual, or as spurious
  width in the distribution.
* No phase correction, and no correction for the noise bias of magnitude mode.
* The noise is treated as uncorrelated and of uniform variance across the window,
  apart from the $\sqrt{m_k}$ weights that resampling introduces. No noise level
  is read from the file or estimated separately; the F-test takes its variance
  estimate from the unregularized fit itself.
* The distribution has no confidence band. `iLT low` and `iLT high` describe the
  width of a peak, not the uncertainty of the amplitudes.
* Where $\alpha$ comes from, the penalty, the grid density, the sweep range and
  the F-test probability are shared by every dataset; the channel, window,
  boxcar, resampling, grid limits, manual $\alpha$ and discrete-fit settings
  belong to one.

---

# Part II — The interface

## 13. Window layout

```
┌ toolbar: [folder][up][load all] | [CSV][session] | [trace][resid][dist][L] …… [grid][theme] ┐
│◂│ file browser │        plot + navigation toolbar        │ datasets + four tabs │
└ status: [message]                                      [cursor readout] ┘
```

Three panes in two splitters, as in pyCWEPR. The chevron at the far left hides
and restores the file browser, remembering where its sash was. The right sash
cannot be dragged narrower than the Datasets card needs, so none of its buttons
is ever clipped; the tab pages scroll sideways instead, and their number boxes
give up width first. `Ctrl+D` toggles dark mode.

## 14. Top toolbar

| Button | What it does |
|---|---|
| Folder | Choose the browser's folder |
| Up | The parent folder |
| Load all | Load every data file the browser lists, skipping folders and sessions; asks first above 20 files, and lists the failures afterwards rather than stopping |
| CSV | Export the current dataset — trace, fit, residual, distribution and discrete curve — as one CSV table (§23) |
| Session | Save every dataset, every setting and the view as one XML file (§23) |
| Rates | Open the rates window: the fitted relaxation rates against a parameter of your own (§24) |
| Trace / Residual / Distribution / L | Show or hide each plot row; each icon is its row's own line shape |
| Grid | Axes grid on all rows |
| Theme | Light or dark (`Ctrl+D`) |

## 15. File browser

A folder listing in which Bruker pairs appear once, under the `.DSC` or `.PAR`
name, with the title Xepr wrote in a second column so the list reads like the
experiment log. An EPRelax session XML is recognized by looking inside it and is
listed as `EPRelax session`. Files already loaded carry a green bookmark.

* **Filter** — which extensions are listed: `EPR` (`.dsc`, `.par`, `.xml`),
  `Bruker`, `Session`, `All`.
* **Title** — a second filter on the stored title, on top of the extension
  filter. A plain word matches anywhere, so `nitroxide` finds `Nitroxide 20 K`;
  `* ? [ ]` say exactly where, as in `T1*` or `*20 K`. `*`, the default, keeps
  everything, including files whose header carries no title. The status bar says
  how many files the keyword hides.
* **Load** — load the selected files; `Ctrl`- or `Shift`-click selects several,
  and a single selected folder is opened instead. Double-clicking a row does the
  same.
* **Header** — show the parsed header of the selected file in a key/value window,
  without loading it.

Both columns sort, folders always staying on top. The path box can be typed in.
The last folder is remembered in `EPRelax.ini` in the user's config directory.

On loading, a trace is rejected if it is not one-or-two-dimensional with one point
per time, or if it has fewer than 8 points. Whether it is a **recovery** is
guessed from the trace itself — the mean magnitude of the last tenth against the
first tenth — and said in the status bar; the switch is in the Signal tab if the
guess is wrong. Duplicate titles get an apostrophe appended. Every dataset is
given the next free colour of the matplotlib cycle.

## 16. Datasets card

Above the tabs, so it is always visible.

* The list: **tick** a dataset to draw it, **click** it to make it the one the
  tabs edit. Each row is written in that dataset's trace colour.
* **Double-click a row**, or press Enter on it, to type that dataset's
  **parameter** into the second column — the temperature, concentration or field
  the trace was measured at, which §24 plots the rates against. Enter or a click
  away keeps it, Escape leaves the old value, anything that is not a number
  clears the cell, and an empty cell simply keeps that dataset out of the rates
  window. The column is headed by whatever the parameter has been named there.
* **Right-click a row** — give that dataset a colour of its own, or hand it back
  to the automatic cycle.
* **Show all** / **Show only** — tick every dataset, or only the selected one.
* **Remove** / **Remove all** — drop the selected dataset, or all of them; the
  second asks first.

The list is **rebuilt only when the datasets themselves change** — one loaded,
one removed, a session opened. A tick, a colour or a parameter is written into
its own row, because a rebuilt list ticks and selects its rows as it fills, and
every one of those reads as a click: with twelve datasets loaded, typing a
parameter used to cost thirteen redraws and 156 analyses, half a second, and the
list scrolled back to the top. It now costs none of them. The rebuild that is
still needed sets a flag while it fills, so its own events are ignored, and it
puts the list back where it was scrolled to:

```python
self.filling = False
if 0 < top < self.filelist.GetItemCount():
    self.filelist.EnsureVisible(self.filelist.GetItemCount() - 1)
    self.filelist.EnsureVisible(top)      # back to the row the list was showing from
```

Each dataset keeps its own channel, slice, window, recovery and offset switches,
boxcar, resampling, grid limits, manual $\alpha$, discrete-fit settings, colour
and discrete fit. The grid density, the penalty, where $\alpha$ comes from, the
sweep range and the F-test probability are shared.

## 17. Signal tab

**Signal** card

| Control | Function |
|---|---|
| **Channel** — Real / Imag. / Magn. | Which channel of a complex dataset is decomposed (§2.1) |
| **Slice**, **Sum** | Which slice of a 2D dataset, or all of them added. Greyed for 1D data |
| **From**, **To** | The time window (§2.2). Both are the dotted lines on the trace and can be dragged there. The label under them recalls the full recorded range |
| **Full** | Both limits back to the ends of the trace |
| **Recovery** | The exponentials rise to the offset instead of decaying, $s=-1$. Set from the trace at load time |
| **Fit offset** | Fit a free constant. On by default: a residual constant otherwise has nowhere to go but the longest-$T$ end of the distribution |

**ESEEM suppression** card

| Control | Function |
|---|---|
| **Suppress** | Run the boxcar before decomposing (§2.3). Off by default — an envelope is a model choice sitting between the file and the result |
| **Width** | Boxcar width $w$ in the dataset's own time unit. Set it to one modulation period, $1/\nu_I$ for the nucleus you see |

**Resampling** card

| Control | Function |
|---|---|
| **Log t** | Bin the trace in $\log t$ (§2.4). The trace row then draws the binned points, which are what is fitted |
| **Per decade** | The bin density $\rho$; a trace spanning two decades at 16 per decade comes out as about 32 points |
| **Weights** — By count / Equal | How much each binned point counts for (§2.4) |
| **Apply to All** | Hand all three settings — the switch, the density and the weighting — to every loaded dataset, so one comparison is not half binned and half raw. Switching **Log t** off here switches it off everywhere the same way |

## 18. Regularization tab

**Relaxation-time grid** card

| Control | Function |
|---|---|
| **Auto** | Take $T_\mathrm{min}$ and $T_\mathrm{max}$ from the trace (§3). The limits used are left in the boxes to start from |
| **T min**, **T max** | The ends of the grid, in the dataset's own time unit — always a time, even when the plot is showing rates |
| **Points** | The grid point count $n$, shared by every dataset |

**Regularization** card

| Control | Function |
|---|---|
| **Penalty** — Curvature / Amplitude | §4. Shared |
| **alpha from** — Manual / L-curve / GCV / F-test | §6. Shared; F-test is the default |
| **Auto range** | Bracket the $\alpha$ sweep from each trace instead of fixing it (§6.1). Shared. Only enabled when a sweep is needed at all — an automatic mode, or the L-curve row on screen |
| **Sweep from**, **Sweep to** | $\log_{10}\alpha_\mathrm{lo}$ and $\log_{10}\alpha_\mathrm{hi}$, used when **Auto range** is off. Shared, so two L-curves stay comparable |
| **log10 alpha** | The manual knob, per dataset. The automatic modes write their choice back into it |
| **Probability** | The F-test's $P$ (§6.5). Shared. Enabled in F-test mode only |

Under them, the statistics of the last decomposition (§10).

Every number box takes its arrow step from a **right-click** menu, as everywhere
else in the suite, and accepts `Enter`.

## 19. Discrete fit tab

| Control | Function |
|---|---|
| **Components** | How many exponentials, $N=1\dots6$ |
| **Hold T at the peaks** | Do not fit the time constants; hold each at a peak of the distribution and fit only the amplitudes and the offset (§8.4) |
| **Fit** | Run the fit on the current dataset, from the distribution as it stands |
| **Fit all** | Hand the component count and the hold switch to every dataset and fit them all, each from its own distribution, so the table compares like with like. Datasets that do not converge are listed afterwards rather than stopping the run |
| **Clear** | Drop the current dataset's discrete fit |
| **Save** | Write every component of every fitted dataset to one CSV table (§23) |

The table has **one column per component** and these rows:

| Row | Meaning |
|---|---|
| `T` (or `k`) | The fitted time constant, in the dataset's unit, or the rate $1/T$ when the plot is in rate mode |
| `± 95%` | The fit's own trust region (§8.3), converted as $\Delta T/T^2$ in rate mode. Prefixed `~` when $T$ was held (§8.4) |
| `a` | The fitted amplitude |
| `± 95%` | Its trust region |
| `iLT low`, `iLT high` | The ends of the distribution peak the component sits in (§9) — a width, not an error. Empty when the component sits on no peak |
| `note` | Why a width should not be trusted, and whether $T$ was held |

Under the table: the fitted offset, the RMS residual, the degrees of freedom the
intervals were taken with, a line when $T$ was held, and a warning when the trace
has been changed since the fit — the curve still follows the parameters, but the
numbers are the old window's.

## 20. Plot tab

| Control | Function |
|---|---|
| **Axis** — Time constant / Rate | What the distribution is drawn against (§10). The grid boxes stay in time either way |
| **Distribution** — Raw / Area / Fraction | The distribution's vertical axis (§10) |
| **Signal** — Raw / 0 to 1 | The experimental trace's vertical axis (§10), the same for every dataset |
| **Regularized fit** | Draw the regularized fit over every trace, and its residual in the residual row |
| **Discrete fit** | Draw the discrete fit, dashed, over every trace that has one, and its residual. Untick both and the residual row is left empty |
| **Legend** | The box of dataset names on the trace row. The L-curve keeps its own $\alpha$ box, which is a readout rather than a key to the colours |
| **Separation** | Vertical gap $\sigma_\mathrm{sep}$ between distribution traces, in units of the tallest (§10) |
| **Components** | The stacked discrete components under each fitted trace (§10) |
| **Drawn as** — Lines / Bands | Thin cumulative curves with a horizontal where each starts, or the gaps filled in |
| **Log t** | Logarithmic time axis on the trace and residual rows |

## 21. The plot

Four rows, drawn in this order and each one toggled from the toolbar, with height
ratios 2.2 / 1.0 / 2.0 / 1.6:

1. **Trace** — every shown dataset's points, with the regularized fit, the
   discrete fit and the stacked components as ticked, the time window as two
   dotted lines, and the legend of dataset names unless it is unticked.
2. **Residual** — the same residuals, sharing the trace's time axis. ESEEM
   modulation, a wrong baseline and a window that still holds the ringdown all
   show up here first.
3. **Distribution** — the amplitudes against $T$ (or $1/T$) on a logarithmic
   axis, filled to the baseline, with a dashed vertical at each discrete time
   constant.
4. **L-curve** — $\eta$ against $\rho$ on log-log axes over the swept $\alpha$,
   with an open marker at the solve being displayed and its $\alpha$ in the
   legend. This is a diagnostic of one decomposition, so it always draws **only the
   current dataset**.

Mouse:

| Action | Effect |
|---|---|
| **Left-click** a trace, residual or distribution | Make that dataset current. The nearest drawn point within 20 pixels wins, so clicking empty space selects nothing — with traces overlaid, picking the wrong one is worse than picking none |
| **Right-click** one | The colour menu for that dataset |
| **Drag** a dotted line on the trace row | Move that window limit; the decomposition follows live, and the boxes update on release |
| **Hover** | The status bar's right field reads out the cursor: $t$ and $y$ on the trace and residual, $T$ (or the rate) and $P$ on the distribution, and the residual norm and $\lVert L\mathbf a\rVert$ on the L-curve |

The matplotlib navigation toolbar sits under the canvas; while a zoom or pan tool
is active, clicking does not select. Redraws during a window resize are debounced.

## 22. Status bar

Two fields: messages on the left — what was loaded, what was fitted, what was
exported, how many files the title filter hides — and the cursor readout on the
right.

Faults that are swallowed rather than shown are printed to the console, prefixed
`[EPRelax]`: a trace with too few points, a resampling that would leave too few,
an L-curve with no corner, a Gaussian peak fit that would not converge, a
singular hat matrix or Jacobian, a decomposition that came out all zeros, and any
file whose title or header could not be read.

## 23. Exports and session files

**CSV export of one dataset** (toolbar). A header of comment lines records the
source file, the channel, the window, the recovery and offset switches, the
boxcar, the penalty, $\alpha$ and where it came from, the effective degrees of
freedom, the offset, the total amplitude and the RMS residual, and — when there
is one — the whole discrete fit including the intervals and the notes. Then the
columns: $t$, data, fit, residual, $T$, amplitude, and the discrete curve. The
file is written with a byte-order mark so a spreadsheet reads the µ of a
microsecond axis as written.

**Combined discrete-fit CSV** (the **Save** button in the Discrete fit tab). One
row per component of every fitted dataset: dataset title, how many components,
which one, $T$, its interval, $a$, its interval, `iLT low`, `iLT high`, the note,
the offset and its interval, the RMS residual, the fit's degrees of freedom,
$\alpha$, the unit the dataset was shown in, and the source path. All times in
seconds (§11).

**Session XML** (toolbar). Every dataset — the raw array included, base64-encoded
with its dtype and shape — every per-dataset setting, the shared settings, and the
whole view state: which rows are shown, the grid, both normalizations, the axis
mode, the separation, the stacking, the log-time switch, the two fit ticks, the
legend tick and which dataset was current. The rates window's own settings (§24)
go in a block of their own, so they are kept whether that window is open or not,
and each dataset's parameter rides along in its own entry. The cached decompositions and sweeps
are left out, as they are rebuilt on the first redraw.

Nothing a card can set is left out of the file: `tmp/check_session_complete.py`
walks the live objects for any plain-valued setting missing from the XML, then
sets all of them away from their defaults, saves, loads into a fresh window and
compares — 12 view, 7 shared and 25 per-dataset settings at the time of writing,
including that the cards themselves show what was loaded, since otherwise the
next click would write the old value straight back. The file reopens from the browser like a data
file, without needing the original Bruker files, and replaces the loaded session
after asking. Settings missing from an older session fall back to the current
defaults.

## 24. The rates window

A separate window, opened from the toolbar: the relaxation rates of the discrete
fits against a number given to each dataset — a temperature, a concentration, a
field, whatever the series was measured along. It is `RateParameterFrame`, and
its settings live on the main window as `MainFrame.rates`, so they survive the
window being closed and are saved with the session.

### 24.1 Where the points come from

Every dataset that has **both** a discrete fit and a parameter contributes one
point per component. Ticked or not: the Datasets ticks say what the trace plot
draws, which is a different question from this one. A dataset without a fit, or
with an empty parameter cell, is counted in the line under the plot rather than
dropped silently.

`points` keeps each component's fitted $T$ and its interval as they stand, and
`quantity` makes the rate of them where a rate is wanted — once for the drawing,
once for the fit, which need not agree (§24.4). The conversion is the Discrete
fit table's own, so a number here can be found there:

$$
\text{Rate:}\quad y_i=\frac1{T_i},\quad \delta y_i=\frac{\delta T_i}{T_i^{2}}
\qquad\qquad
\text{Time constant:}\quad y_i=T_i,\quad \delta y_i=\delta T_i
$$

where $\delta T_i$ is the component's 95 % interval (§8.3) — or, for a held time
constant, the as-if-free interval of §8.4, which is drawn like any other. In
`RateParameterFrame.quantity`:

```python
T, dT = row['T'], row['dT']
return (1.0/T, dT/T**2) if mode == 0 else (T, dT)
```

The time unit is the one the datasets share (§11); a dataset whose file named a
unit off the ladder keeps its own numbers and the line under the plot says so.

### 24.2 Series

A dataset holds up to six components, and which of them is "the same" component
as in the next dataset is a judgement, not a fact. Three answers are offered:

| Series | What one series is |
|---|---|
| **By T** | Component 1 is the fastest in every dataset, 2 the next — the order the Discrete fit table already shows. Right while the components keep their order along the series |
| **By ampl.** | Series 1 is the tallest component of each dataset, whatever its $T$. Follows the dominant process, but jumps where two amplitudes cross |
| **All in one** | Every component as one set of points, matched up with nothing — what to draw when the components are not comparable |

```python
order = np.argsort(d['a'])[::-1] if rates['series'] == 1 else range(len(d['T']))
```

### 24.3 Axes

| Control | Effect |
|---|---|
| **Y axis** — Rate / Time const. | $1/T$ or $T$, with the interval converted to match (§24.1). This is the drawing alone; what the fit works in is its own control (§24.4) |
| **X axis** — Parameter / 1 / parameter | The parameter $p$ as typed, or $1/p$ — an Arrhenius plot. The fit is unaffected (§24.4) |
| **Log x**, **Log y** | Logarithmic axes, either or both. A point that cannot be drawn on a log axis, or a zero parameter under $1/p$, is left out and counted in the line under the plot |
| **Error bars** | The 95 % intervals, drawn with caps. Untick for bare points |
| **Name** | What the parameter is called: it heads the parameter column in the Datasets list and labels this axis. The unit is part of the name — "Temperature, K" — since nothing is ever converted |

### 24.4 The fit

The equation is typed in, as `y` in terms of `x` and any other names you like:
`a*exp(-b/x) + c`, `a + b*x`, `a*x**n`. Every name that is not `x` and not a
function becomes a fitted parameter with a value box and a **Hold** tick. The
suite's own expression reader does the parsing, so `^` works as `**` and the
function names are the ones the Sys function modifiers use:

```python
work, _tokens, free = extract_refs_and_freevars(self.mf.rates['equation'], (), FIT_FUNCS)
return work, sorted(n for n in free if n != 'x')
```

`x` is **always the parameter itself**, never its inverse, whatever the X axis is
showing: flipping the axis re-draws the same fit rather than changing what was
fitted. The curve is evaluated over the parameter and then drawn against
whichever axis is up.

**Fit in** — Rate / Time const. — says what `y` is to the equation, and it is
**independent of the Y axis**. The two are different questions: a rate law is
written in $1/T$, while a plot often reads better in $T$, and neither should
force the other. The control sets

* the quantity the equation describes, so its parameters mean what they mean in
  that quantity;
* the residuals the fit minimizes, and the RMSD and $R^2$ reported — the
  statistics line names the quantity whenever it is not the one on the axis;
* which interval becomes the weight, $\delta(1/T)=\delta T/T^2$ or $\delta T$.

Where the two differ the drawn curve is turned over, $y\to1/y$, so it still lies
on the points, and the legend reads `fit in 1/T` or `fit in T` rather than plain
`fit`. The fit itself is untouched by the axis: fitting `a + b*x` to the rates
gives the same $a$ and $b$ whether the plot shows $1/T$ or $T$, which the harness
checks by fitting both ways round.

Changing **Fit in** drops the kept fit, since it was a fit of the other quantity.

A term drawn by **Term by term** (§24.6) is a term of the *fitted* quantity and
is turned over with the whole curve, so on a flipped axis the term lines no
longer add up to the curve — $1/(u+v)\ne1/u+1/v$. They still show what each term
contributes to the quantity that was fitted.

**Weight by error bars** (on by default) gives each point the weight its own
interval earns. The stored interval is $t_{0.975,\nu_i}$ standard errors wide,
with $\nu_i$ the degrees of freedom of the discrete fit it came from, so the
standard error itself is

$$
\sigma_i=\frac{\Delta_{95}(T_i)}{t_{0.975,\,\nu_i}},
$$

which is what `curve_fit` is given, with `absolute_sigma=True` so the parameter
intervals below mean what they say:

```python
sigma = dy/np.array([float(tdist.ppf(0.975, max(v, 1))) for v in dof])
```

A point whose interval is zero or not finite is given the median of the others,
and if no point has a usable one the fit falls back to unweighted — both say so
on the console. Unticking the box fits every point equally.

```python
popt, pcov = curve_fit(lambda xx, *vals: self.curve(xx, {**rates['pars'], **dict(zip(free, vals))}),
                       x, y, p0=[rates['pars'][n] for n in free], sigma=sigma,
                       absolute_sigma=sigma is not None, maxfev=40000)
```

The reported interval on each fitted parameter is the linearized one at the
optimum, as in §8.3,

$$
\Delta_{95}(\theta)=t_{0.975,\,n-p}\sqrt{\operatorname{diag}\mathbf C},
$$

with $n$ the points fitted and $p$ the free parameters, beside the RMSD and

$$
R^2=1-\frac{\sum_i r_i^2}{\sum_i\left(y_i-\bar y\right)^2}.
$$

**Fit to** chooses the points: all of them, or one series on its own. One fit is
kept at a time. A hand-typed equation has no automatic starting guess, so the
value boxes are where the fit starts from — within an order of magnitude is
usually enough, and a fit that does not converge says so rather than moving the
boxes.

### 24.5 The curve follows the boxes

**Draw the equation** (on by default) draws the equation over the points with
whatever is in the value boxes, fitted or typed. Turn a value and the curve moves
with it at once, which is how a hand-typed equation is brought near enough to fit
from — there is no need to fit first, and a curve is drawn before any fit has
been run.

Because the drawn curve is no longer necessarily the fitted one, the numbers
under it say which it is:

| | Legend | Parameters | RMSD and R² | Intervals |
|---|---|---|---|---|
| Straight after a fit | `fit` | the fitted values | of the fitted curve | shown, with the degrees of freedom |
| After typing over them | `equation` | what is in the boxes | recomputed for the drawn curve | gone, with a line saying why |
| Before any fit | `equation` | what is in the boxes | of the drawn curve | none to show |

RMSD and R² are always those of the curve on screen, over the points **Fit to**
picks out, so they move with the boxes and can be watched while a starting value
is dialled in. The 95 % intervals belong to the fit alone and are dropped the
moment a value is typed over, rather than being left beside numbers they no
longer describe.

**Clear** drops the fit — its intervals and its claim on the curve — while the
equation and the values stay, so the curve remains until the tick is cleared.
Both the switch and the parameter values are saved with the session.

### 24.6 Term by term

**Term by term** draws each term of the equation on its own as well, dashed and
in its own colour, with the same parameter values as the whole curve — so a sum
of two processes can be seen as the two it is made of, which is usually the
reason for writing it as a sum:

$$
y=\underbrace{a\,e^{-b/x}}_{\text{one term}}+\underbrace{c}_{\text{another}} .
$$

`RateParameterFrame.terms` breaks the equation at every `+` and `-` that is not
inside brackets. The sign goes with the term that follows it, so the terms always
add back up to the whole curve — which the harness checks on every example below,
at three values of $x$:

| Equation | Terms |
|---|---|
| `a*exp(-b/x) + c` | `a*exp(-b/x)`, `+ c` |
| `a - b*x + c` | `a`, `- b*x`, `+ c` |
| `-a + b` | `-a`, `+ b` |
| `a*(x + c) - d` | `a*(x + c)`, `- d` |
| `a*x**-2 + b` | `a*x**-2`, `+ b` |
| `a/(1 + b*x) - 3.5e-2*x` | `a/(1 + b*x)`, `- 3.5e-2*x` |
| `a*exp(-(x-b)**2/c) + d*x - e` | `a*exp(-(x-b)**2/c)`, `+ d*x`, `- e` |
| `a*x**n` | `a*x**n` — one term, so nothing is drawn |

Three kinds of sign are **not** a break: the one in front of the first term, the
one in an exponent like `x**-2`, and the one inside a number like `3.5e-2`.

```python
elif depth == 0 and ch in '+-' and k > start:
    before = work[:k].rstrip()
    exponent = before[-1:] in ('e', 'E') and (before[-2:-1].isdigit() or before[-2:-1] == '.')
    if before and before[-1] not in '+-*/^(,' and not exponent:
        out.append(work[start:k])
        start = k
```

The term colours carry on past the ones the data series have taken, so a term
line is never the colour of a set of points. An equation of a single term has
nothing to break apart and draws nothing extra. The switch is saved with the
session, and is greyed while **Draw the equation** is off.

### 24.7 Export

**Save** writes one CSV: the equation as a comment line — headed `fit` or
`equation` by the same rule as §24.5 — with its parameter values, their intervals
where the fit still owns them, and the RMSD and R² of the curve as it stands;
then a row per point with its series, dataset title, parameter, plotted quantity,
95 % interval, and whether its time constant was held.
