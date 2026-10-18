import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Parameters
# ============================================================
fs = 3.072e6       # Sampling frequency [Hz]
f_band = 20e3      # Signal bandwidth [Hz]

# Frequency axis: 0 ~ Nyquist
N = 200000
f = np.linspace(0, fs / 2, N)

# Digital angular frequency
w = 2 * np.pi * f / fs

# z = exp(jw)
z = np.exp(1j * w)

# ============================================================
# STF / NTF
# ============================================================
STF = 1 / (z**2 - z + 1)

NTF = (z - 1)**2 / (z**2 - z + 1)

# Magnitude [dB]
STF_dB = 20 * np.log10(np.abs(STF))
NTF_dB = 20 * np.log10(np.abs(NTF))

# ============================================================
# Plot 1: Full frequency range
# ============================================================
plt.figure(figsize=(10, 5))

plt.plot(f / 1e3, STF_dB, label="STF")
plt.plot(f / 1e3, NTF_dB, label="NTF")

plt.axvline(f_band / 1e3, linestyle="--", label="20 kHz")
plt.axhline(0, linestyle="--")

plt.xlabel("Frequency [kHz]")
plt.ylabel("Magnitude [dB]")
plt.title("STF / NTF Frequency Response")
plt.xlim(0, fs / 2 / 1e3)
plt.ylim(-100, 20)
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# ============================================================
# Plot 2: Signal bandwidth
# ============================================================
idx = f <= f_band

plt.figure(figsize=(10, 5))

plt.plot(f[idx] / 1e3, STF_dB[idx], label="STF")
plt.plot(f[idx] / 1e3, NTF_dB[idx], label="NTF")

plt.axhline(0, linestyle="--")

plt.xlabel("Frequency [kHz]")
plt.ylabel("Magnitude [dB]")
plt.title("STF / NTF in Signal Band (0-20 kHz)")
plt.xlim(0, f_band / 1e3)
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

# ============================================================
# Numerical values at DC and 20 kHz
# ============================================================
for freq in [0, 1e3, 10e3, 20e3]:
    ww = 2 * np.pi * freq / fs
    zz = np.exp(1j * ww)

    stf = 1 / (zz**2 - zz + 1)
    ntf = (zz - 1)**2 / (zz**2 - zz + 1)

    print(
        f"{freq/1e3:6.1f} kHz : "
        f"STF = {20*np.log10(abs(stf)):8.4f} dB, "
        f"NTF = {20*np.log10(abs(ntf)) if freq != 0 else -np.inf:8.4f} dB"
    )
