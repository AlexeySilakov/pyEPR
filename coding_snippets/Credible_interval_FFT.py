import numpy as np
import matplotlib.pyplot as plt
np.random.seed(1)
fs = 100.0        # sample rate (Hz)
N = 256           # number of samples
Ngrid = 200
amp = 1
noiseamp = 0.1
t  = np.arange(N)/fs
t0 = 0            # FID: packet is largest at t=0 and decays
f_true = 10        # the frequency
fwhm_true = 3    # linewidth
tau  = 2.3548/(2*np.pi*fwhm_true)
packet   = amp*np.cos(2*np.pi*f_true*(t-t0))*np.exp(-(t-t0)**2/(2*tau**2))
baseline = 0.3*t + 0.02*t**2
raw      = packet + baseline + np.random.normal(0, noiseamp, N)

#### Analysis
raw_bs = raw - np.polyval(np.polyfit(t, raw, 2), t)

# --- apodization: HALF-Hamming (max at t=0, correct for an FID) + 4x zero-fill ---
zf   = 4
Nz   = zf*N
wham = np.hamming(2*N)[N:]                 # descending half -> ~1 at t=0, 0.08 at end
sig_apod = np.zeros(Nz);
sig_apod[:N] = raw_bs*wham

freq = np.fft.rfftfreq(Nz, 1/fs)
spec = np.abs(np.fft.rfft(sig_apod))

# amplitude KNOWN: propagate the noiseless line through the SAME apodization pipeline

cl_apod  = np.zeros(Nz);
cl_apod[:N] = packet*wham
spec_known = np.abs(np.fft.rfft(cl_apod))
amp_known = spec_known.max()-spec_known[-1]
# --- noise floor estimated from the UPPER HALF of the spectrum (signal-free) ---
floor_est = spec[freq > freq.max()/2].mean()

f0_grid = np.linspace(f_true-fwhm_true/2, f_true+fwhm_true/2, Ngrid)   # candidate centre frequencies
fw_grid = np.linspace(0.1*fwhm_true, 2*fwhm_true, Ngrid)              # candidate linewidths (FWHM), widened for apodization
F0, FW = np.meshgrid(f0_grid, fw_grid)
s = FW/2.3548
model = amp_known*np.exp(-(freq - F0[:, :, None])**2/(2*s[:, :, None]**2)) + floor_est
SSR = np.sum((spec - model)**2, axis=2)

imin = np.unravel_index(np.argmin(SSR), SSR.shape)
Neff = freq.size
sigma_hat = np.sqrt(SSR[imin]/(Neff-2))
L = np.exp(-(SSR - SSR[imin])/(2*zf*sigma_hat**2))
post = L/L.sum()
flat = np.sort(post.ravel())[::-1]
cutoff = flat[np.searchsorted(np.cumsum(flat), 0.95)]
mask = post >= cutoff
p_f = post.sum(axis=0); p_f /= p_f.sum()
p_w = post.sum(axis=1); p_w /= p_w.sum()
f_lo, f_hi = np.interp([0.025, 0.975], np.cumsum(p_f), f0_grid)
w_lo, w_hi = np.interp([0.025, 0.975], np.cumsum(p_w), fw_grid)
f0_hat, fw_hat = f0_grid[imin[1]], fw_grid[imin[0]]
f_half, w_half = (f_hi-f_lo)/2, (w_hi-w_lo)/2 # +/- half interval


# ---- alternative: marginalize sigma analytically with Jeffreys prior p(sigma) ∝ 1/sigma ----
#   ∫_0^inf sigma^-(n+1) exp(-SSR/2sigma^2) dsigma = (1/2) Gamma(n/2) (2/SSR)^(n/2)
#   => sigma-marginalized posterior  p(theta) ∝ SSR^(-n_eff/2)   (Student-t marginal)
n_eff=freq.size/zf
logJ = -(n_eff/2.0)*np.log(SSR); logJ -= logJ.max()
postJ = np.exp(logJ); postJ /= postJ.sum()
flatJ = np.sort(postJ.ravel())[::-1]
cutJ = flatJ[np.searchsorted(np.cumsum(flatJ), 0.95)]
maskJ = postJ >= cutJ
pfJ = postJ.sum(0); pfJ /= pfJ.sum()
pwJ = postJ.sum(1); pwJ /= pwJ.sum()
fJ_lo, fJ_hi = np.interp([0.025, 0.975], np.cumsum(pfJ), f0_grid)
wJ_lo, wJ_hi = np.interp([0.025, 0.975], np.cumsum(pwJ), fw_grid)
fJ_half, wJ_half = (fJ_hi-fJ_lo)/2, (wJ_hi-wJ_lo)/2 # +/- half interval


# ============================================================================
#  plots
# ============================================================================
fig, ax = plt.subplots(2, 2, figsize=(16, 4.6), num='fft_bayes', clear=True)
ax[0,0].plot(t, raw, lw=4, label='raw')
ax[0,0].plot(t, np.polyval(np.polyfit(t, raw, 2), t), 'r', lw=2, label='fitted baseline')

ax[0,0].plot(t, raw_bs, ':k', lw=.7, label='baselined')
ax[0,0].plot(t, wham, ':b', lw=.7, label='apodization')

ax[0,0].set_title(f'A) time domain, freq={f_true:.3f} Hz FWHH={fwhm_true:.3f} Hz');
ax[0,0].set_xlabel('t (s)');
ax[0,0].legend(fontsize=8)

ax[1,0].plot(freq, spec_known,
           'k', lw=1, label='fft(noiseless)')

ax[1,0].plot(freq, spec, '.', ms=3, label='fft(exp)')
ax[1,0].plot(freq, model[imin[0], imin[1], :],
           'r', lw=2, label='best-fit Gaussian + floor')
ax[1,0].plot(freq, spec - model[imin[0], imin[1], :], ':', ms=1, label='difference')
ax[1,0].axhline(floor_est, color='green', ls=':', label=f'noise floor = {floor_est:.2f}')
ax[1,0].set_title(f'B) FFT σ_noise={sigma_hat/amp_known:.2f}');
ax[1,0].set_xlabel('frequency (Hz)')
ax[1,0].legend(fontsize=8)

pc = ax[0,1].pcolormesh(f0_grid, fw_grid, post, shading='auto', cmap='viridis')
ax[0,1].contour(f0_grid, fw_grid, mask.astype(float), [0.5], colors='white', linewidths=2)
ax[0,1].axvline(f_true, color='gray', ls=':', label=f'noise floor = {floor_est:.2f}')
ax[0,1].axhline(fwhm_true, color='gray', ls=':', label=f'noise floor = {floor_est:.2f}')
#ax[0,1].plot(f0_hat, fw_hat, 'r*', ms=7)
ax[0,1].set_title(f'C) set σ. freq = {f0_hat:.3f} ± {f_half:.3f} Hz    '
                f'FWHM = {fw_hat:.3f} ± {w_half:.3f} Hz')
ax[0,1].set_xlabel('frequency (Hz)');
ax[0,1].set_ylabel('linewidth FWHM (Hz)')
fig.colorbar(pc, ax=ax[0,1])

pc1 = ax[1,1].pcolormesh(f0_grid, fw_grid, postJ, shading='auto', cmap='viridis')
ax[1,1].contour(f0_grid, fw_grid, maskJ.astype(float), [0.5], colors='white', linewidths=2)

# ax[1,1].plot(f0_hat, fw_hat, 'r*', ms=7)

ax[1,1].set_title(f'D) Jeffreys prior. freq = {f0_hat:.3f} ± {fJ_half:.3f} Hz    '
                f'FWHM = {fw_hat:.3f} ± {wJ_half:.3f} Hz')
ax[1,1].axvline(f_true, color='gray', ls=':', label=f'noise floor = {floor_est:.2f}')
ax[1,1].axhline(fwhm_true, color='gray', ls=':', label=f'noise floor = {floor_est:.2f}')

ax[1,1].set_xlabel('frequency (Hz)');
ax[1,1].set_ylabel('linewidth FWHM (Hz)')
fig.colorbar(pc, ax=ax[1,1])


plt.tight_layout(); plt.savefig('fft_gaussian_bayes.png', dpi=130)
print(f'true  : f = {f_true} Hz , FWHM = {fwhm_true} Hz (pre-apodization)')
print(f'fit   : f = {f0_hat:.3f} +/- {f_half:.3f} Hz , FWHM = {fw_hat:.3f} +/- {w_half:.3f} Hz')
print(f'amp_known = {amp_known:.2f} , noise floor = {floor_est:.3f} , sigma_hat = {sigma_hat:.3f}')