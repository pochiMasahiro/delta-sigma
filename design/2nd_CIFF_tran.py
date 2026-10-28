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
fs = 12.288e6
#fin = 996.09375
amplitude = 0.5

N_FFT = 2**18
N_transient = 8192
N = N_FFT + N_transient

k = 21
fin = k*fs/N_FFT

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

u = amplitude * np.sin(2*np.pi*fin*t)
#u = np.ones(len(n))*amplitude

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
def dsm_sim(sigma, c1, c2, a1, a2, b1, b2, e1, e2):
    dither = rng.normal(
        loc=0.0,
        scale = sigma,
        size = N
    )

    x1 = np.zeros(N + 1)
    x2 = np.zeros(N + 1)
    v = np.zeros(N)
    y = np.zeros(N)

    for i in range(N):

        v[i] = c1*x1[i] + c2*x2[i]

        # 1-bit quantizer
        y[i] = 1.0 if v[i] + dither[i] >= 0 else -1.0

        # First integrator
        x1t = a1*x1[i] + b1*(u[i] - y[i])

        # Second-integrator
        x2t = a2*x2[i] + b2*x1[i]

        x1[i + 1] = x1[i] + e1 * (x1t - x1[i])
        x2[i + 1] = x2[i] + e2 * (x2t - x2[i])

        if (
            not np.isfinite(x1[i+1])
            or not np.isfinite(x2[i+1])
            or abs(x1[i+1]) > 1e6
            or abs(x2[i+1]) > 1e6
        ):
            raise RuntimeError(
                f"Excessive internal state at sample {i}"
            )
    return x1, x2, y, v

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
    f = np.fft.rfftfreq(Num, d=1/fs)

    # One-sided power spectrum
    power = np.abs(Y)**2 / Num**2

    if Num % 2 == 0:
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


def find_period(y, max_period=4096):
    y = np.asarray(y)

    for period in range(1, max_period+1):
        if np.array_equal(
            y[:-period],
            y[period:]
        ):
            return period
    return None

def calculate_sndr(y , bandwidth, tone_bin):
    y = np.asarray(y, dtype=float)
    Num = len(y)

    y_ac = y - np.mean(y)
    Y = np.fft.rfft(y_ac)
    f = np.fft.rfftfreq(Num, d=1/fs)

    power = np.abs(Y)**2 / Num**2

    if Num % 2 == 0:
        power[1:-1] *= 2
    else:
        power[1:] *= 2

    band = (f > 0) & (f <= bandwidth)

    if not band[tone_bin]:
        raise ValueError("Fundamental is outside bandwidth")

    signal_power = power[tone_bin]
    noise_mask = band.copy()
    noise_mask[tone_bin] = False

    noise_distortion_power = np.sum(power[noise_mask])

    sndr = 10 * np.log10(
        signal_power / noise_distortion_power
    )

    signal_dbfs = 10 * np.log10(signal_power / 0.5)
    noise_dbfs = 10 * np.log10(
        noise_distortion_power / 0.5
    )

    enob = (sndr - 1.76) / 6.02

    return{
        "sndr": sndr,
        "enob": enob,
        "signal_dbfs": signal_dbfs,
        "noise_dbfs": noise_dbfs,
        "frequency": f,
        "power": power
    }

# ===================
# DSM simulation
# ===================

def dsm_simulation(c1, c2, a1, a2, b1, b2, e1, e2):
    x1, x2, y, v = dsm_sim(0.00, c1, c2, a1, a2, b1, b2, e1, e2)
    #p0, db0 = calculate_inband_noise(y, 20e3)
    #y = y[N_transient:]
    #result = calculate_tone_power(y)
#    print("====================")
#    print("A: {0:f}".format(A))
#    print("B: {0:f}".format(B))
#    #print("NTF: {0:f}".format(
#    print("Inband noise without Dither: {0:f}".format(db0))
#    print("Tone Frequency: {0:f}".format(result["frequency"]))
#    print("Tone power: {0:f}".format(result["tone_power"]))
#    print("Tone level: {0:f}".format(result["tone_dbfs"]))
#    print("Peak prominence: {0:f}".format(result["prominence_db"]))
#    dc_output = np.mean(y)
#    print("DC (mean): {0:f}".format(dc_output))

    y_fft = y[N_transient:N_transient+N_FFT]
    result = calculate_sndr(
        y_fft,
        bandwidth = 20e3,
        tone_bin = k
    )

    print("===================================")
    print(f"C1: {c1:.3f}, C2: {c2:.3f}")
    print(f"a1: {a1:.3f}, b1: {b1:.3f}, e1: {e1:.3f}")
    print(f"a2: {a2:.3f}, b2: {b2:.3f}, e2: {e2:.3f}")
    print(f"Signal : {result["signal_dbfs"]:.3f} dBFS")
    print(f"Noise+D: {result["noise_dbfs"]:.3f} dBFS")
    print(f"SNDR   : {result["sndr"]:.3f} dB")
    print(f"ENOB   : {result["enob"]:.3f} bit")

    print()
#    fig = plt.figure(figsize=(5,7))
#    gs_master = GridSpec(nrows=3, ncols=1, height_ratios=[1,1,1])
#    gs_1 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[0, 0])
#    gs_2 = GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_master[1, 0])
#    gs_3 = GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs_master[2, 0])
#    ax1 = fig.add_subplot(gs_1[:,:])
#    ax2 = fig.add_subplot(gs_2[:,:])
#    ax3 = fig.add_subplot(gs_3[:,:])
#
#    M = 200
#
#    ax1.step(t[:M] * 1e6, y[:M], where="post")
#    ax1.plot(t[:M] * 1e6, u[:M])
#    ax1.grid(ls='--')
#    ax1.set_xlabel('Time [us]')
#    ax1.set_ylabel('Amplitude [V]')
#    ax1.set_title("PDM output")
#
#    freq, psd_db = calculated_PSD(y_fft)
#
#    ax2.plot(freq/1e3, psd_db)
#    ax2.grid(ls='--')
#    ax2.set_xlabel('Frequency [kHz]')
#    ax2.set_ylabel('PSD [dBFS/Hz]')
#    ax2.set_title("PDM Output Spectrum")
#
#    M = 2000
#
#    ax3.plot(t[:M] * 1e6, x1[:M])
#    ax3.plot(t[:M] * 1e6, x2[:M])
#    ax3.set_xlim()
#    ax3.set_ylim()
#    ax3.grid(ls='--')
#    ax3.set_xlabel('Time')
#    ax3.set_title("Integrator internal states")
#
#    fig.tight_layout()

    return result

if __name__ == "__main__":
    c1_list = [0.8, 1.0, 1.2, 1.4, 1.6]
    c2_list = [0.2, 0.4, 0.6, 0.8]
    a1_list = [0.999999, 0.99999, 0.9999, 0.999, 0.99];
    a2_list = [0.999999, 0.99999, 0.9999, 0.999, 0.99];

    sndr = np.zeros((len(a2_list), len(a1_list)))
    enob = np.zeros((len(a2_list), len(a1_list)))
    XX, YY = np.meshgrid(a1_list, a2_list)

    for x, a1 in enumerate(a1_list):
        for y, a2 in enumerate(a2_list):
            # def dsm_sim(sigma, c1, c2, a1, a2, b1, b2, e1, e2):
            result = dsm_simulation(1.0, 0.4, a1, a2, 1, 1, 1, 1)
            sndr[y, x] = result["sndr"]
            enob[y, x] = result["enob"]

    idx = np.argmax(sndr)
    coords = np.unravel_index(idx, sndr.shape)
    print(f"Max SNDR: {sndr[coords]:.3f} dB")


    fig = plt.figure(figsize=(10, 5))
    gs_master = GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])
    gs_1 = GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs_master[0, 0])
    gs_2 = GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs_master[0, 1])
    #gs_3 = GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_master[0, 2])
    ax1 = fig.add_subplot(gs_1[:,:])
    ax2 = fig.add_subplot(gs_2[:,:])
    #ax3 = fig.add_subplot(gs_3[:,:])

    conf = ax1.contourf(XX, YY, sndr)
    ax1.set_xlim()
    ax1.set_ylim()
    #ax1.grid(ls='--')
    ax1.set_xlabel('C1')
    ax1.set_ylabel('C2')
    ax1.set_title("SNDR")
    #ax1.set_aspect("equal")

    conf = ax2.contourf(XX, YY, enob)
    ax2.set_xlim()
    ax2.set_ylim()
    #ax2.grid(ls='--')
    ax2.set_xlabel('C1')
    ax2.set_ylabel('C2')
    ax2.set_title("ENOB")


    #fig.colorbar(conf, cax=ax2)

    fig.tight_layout()
    #plt.savefig('.png', format='png', dpi=300, bbox_inches='tight')
    plt.show()
