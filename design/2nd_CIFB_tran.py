import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
plt.style.use('ggplot')
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.signal import welch, find_peaks
from scipy.ndimage import median_filter

# ===================
# Parameters
# ===================
fs = 3.072e6
fin = 1e3
amplitude = 0.05
N = 262144

# ===================
# Loop parameter
# ===================
#A = 0.375
#B = 1.0

# ===================
# Input signal
# ===================
n = np.arange(N)
t = n/fs

#u = amplitude * np.sin(2*np.pi*fin*t)
u = np.ones(len(n))*amplitude

# ===================
# Dither random gen
# ===================
rng = np.random.default_rng(1234)
#sigma = 0.01
#dither = rng.normal(
#    loc=0.0,
#    scale = sigma,
#    size = N
#)


# ===================
# Second-order 1-bit delta-sigma modulator
# ===================
def dsm_sim(sigma, A, B):
    dither = rng.normal(
        loc=0.0,
        scale = sigma,
        size = N
    )

    x1 = np.zeros(N + 1)
    x2 = np.zeros(N + 1)
    y = np.zeros(N)

    for i in range(N):
        # 1-bit quantizer
        y[i] = 1.0 if x2[i] + dither[i] >= 0 else -1.0

        # First integrator
        x1[i + 1] = x1[i] + u[i] - y[i]

        # Second-integrator
        x2[i + 1] = x2[i] + A*x1[i] - B*y[i]

    return x1, x2, y

# ===================
# PSD calculation
# ===================
def calculated_PSD(y):
    f, psd = welch(
        y,
        fs = fs,
        window = "hann",
        nperseg = 65536,
        noverlap = 32768,
        detrend = "constant",
        scaling = "density"
    )

    psd_db = 10 * np.log10(
        np.maximum(psd/0.5, 1e-30)
    )

    return f, psd_db

def calculate_inband_noise(y, bandwidth):
    y = np.asarray(y, dtype=float)
    Num = len(y)

    y_ac = y - np.mean(y)

    Y = np.fft.rfft(y_ac)
    f = np.fft.rfftfreq(N, d=1/fs)

    # One-sided power spectrum
    power = np.abs(Y)**2 / N**2

    if N % 2 == 0:
        power[1:-1] *= 2
    else:
        power[1:] *= 2

    band = (f > 0) & (f <= bandwidth)
    noise_power = np.sum(power[band])

    noise_dbfs = 10 * np.log10(
        max(noise_power, 1e-30) / 0.5
    )

    return noise_power, noise_dbfs

def calculate_tone_power(
    y,
    search_range = (11000, 13000),
    nperseg = 65536
):
    f, psd = welch(
        y,
        fs = fs,
        window = "hann",
        nperseg = nperseg,
        noverlap = nperseg/2,
        detrend = "constant",
        scaling = "density"
    )

    df = f[1] - f[0]
    search = (
        (f >= search_range[0]) &
        (f <= search_range[1]))

    idx = np.flatnonzero(search)
    peak_idx = idx[np.argmax(psd[idx])]
    peak_freq = f[peak_idx]

    noise_region = (
        (np.abs(f - peak_freq) >= 400) &
        (np.abs(f - peak_freq) <= 1500)
    )

    noise_floor = np.median(psd[noise_region])

    prominence_db = 10 * np.log10(
        max(psd[peak_idx] / noise_floor, 1e-30)
    )
    
    lo = max(1, peak_idx - 3)
    hi = min(len(psd), peak_idx + 4)

    peak_power = np.sum(psd[lo:hi]) * df

    background_power = (hi - lo) * noise_floor * df

    tone_power = max(
        peak_power - background_power,
        0.0
    )

    tone_dbfs = 10 * np.log10(
        max(tone_power, 1e-30) / 0.5
    )

    return {
        "frequency": peak_freq,
        "tone_power": tone_power,
        "tone_dbfs": tone_dbfs,
        "prominence_db": prominence_db,
        "bin_width": df
    }

# ===================
# DSM simulation
# ===================

def dsm_simulation(A, B):
    x1, x2, y = dsm_sim(0.00, A, B)
    p0, db0 = calculate_inband_noise(y, 20e3)
    result = calculate_tone_power(y)
    print("====================")
    print("A: {0:f}".format(A))
    print("B: {0:f}".format(B))
    #print("NTF: {0:f}".format(
    print("Inband noise without Dither: {0:f}".format(db0))
    print("Tone Frequency: {0:f}".format(result["frequency"]))
    print("Tone power: {0:f}".format(result["tone_power"]))
    print("Tone level: {0:f}".format(result["tone_dbfs"]))
    print("Peak prominence: {0:f}".format(result["prominence_db"]))
    dc_output = np.mean(y)
    print("DC (mean): {0:f}".format(dc_output))
    print()
    fig = plt.figure(figsize=(5,7))
    gs_master = GridSpec(nrows=3, ncols=1, height_ratios=[1,1,1])
    gs_1 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[0, 0])
    gs_2 = GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_master[1, 0])
    gs_3 = GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs_master[2, 0])
    ax1 = fig.add_subplot(gs_1[:,:])
    ax2 = fig.add_subplot(gs_2[:,:])
    ax3 = fig.add_subplot(gs_3[:,:])

    M = 200

    ax1.step(t[:M] * 1e6, y[:M], where="post")
    ax1.plot(t[:M] * 1e6, u[:M])
#ax1.set_xlim()
#ax1.set_ylim()
    ax1.grid(ls='--')
    ax1.set_xlabel('Time [us]')
    ax1.set_ylabel('Amplitude [V]')
    ax1.set_title("PDM output")

    freq, psd_db = calculated_PSD(y)

    ax2.plot(freq/1e3, psd_db)
#ax2.set_xlim()
#ax2.set_ylim()
    ax2.grid(ls='--')
    ax2.set_xlabel('Frequency [kHz]')
    ax2.set_ylabel('Amplitude [dBFS, peak-reference]')
    ax2.set_title("PDM Output Spectrum")

    M = 2000

    ax3.plot(t[:M] * 1e6, x1[:M])
    ax3.plot(t[:M] * 1e6, x2[:M])
    ax3.set_xlim()
    ax3.set_ylim()
    ax3.grid(ls='--')
    ax3.set_xlabel('Time')
    ax3.set_title("Integrator internal states")

    fig.tight_layout()
#plt.savefig('.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    dsm_simulation(0.3, 0.9)
    dsm_simulation(0.35, 1)
    dsm_simulation(0.4, 1)
    dsm_simulation(0.45, 1)
    dsm_simulation(0.5, 1)
