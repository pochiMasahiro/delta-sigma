import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
plt.style.use('ggplot')
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


# ===================
# Parameters
# ===================
fs = 3.072e6
fin = 1e3
amplitude = 0.5
N = 65536

# ===================
# Loop parameter
# ===================
A = 0.375
B = 1.0

# ===================
# Input signal
# ===================
n = np.arange(N)
t = n/fs

u = amplitude * np.sin(2*np.pi*fin*t)

# ===================
# Second-order 1-bit delta-sigma modulator
# ===================
x1 = np.zeros(N + 1)
x2 = np.zeros(N + 1)
y = np.zeros(N)

for i in range(N):
    # 1-bit quantizer
    y[i] = 1.0 if x2[i] >= 0 else -1.0

    # First integrator
    x1[i + 1] = x1[i] + u[i] - y[i]

    # Second-integrator
    x2[i + 1] = x2[i] + A*x1[i] - B*y[i]

# ===================  
# Plot PDM
# ===================

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

start = 4096
data = y[start:]

window = np.hanning(len(data))
Y = np.fft.rfft(data * window)

freq = np.fft.rfftfreq(len(data), 1/fs)

scale = np.sum(window) / 2
spectrum = np.abs(Y) / scale
spectrum[0] /= 2

spectrum_dB = 20 * np.log10(
    np.maximum(spectrum, 1e-15)
)

ax2.plot(freq/1e3, spectrum_dB)
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



