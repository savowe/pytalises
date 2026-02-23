#!/usr/bin/env python3
"""Generate a simple GIF animation of a Gaussian wavepacket for the README."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path

# Simulation parameters
N = 512
x = np.linspace(-15, 15, N)
dx = x[1] - x[0]
k = np.fft.fftfreq(N, dx) * 2 * np.pi

# Initial Gaussian wavepacket with momentum
sigma = 1.0
k0 = 3.0  # initial momentum
psi = np.exp(-x**2 / (2 * sigma**2)) * np.exp(1j * k0 * x)
psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)

# Time evolution parameters
dt = 0.02
n_frames = 120

# Kinetic propagator (split-step)
kinetic_prop = np.exp(-0.5j * k**2 * dt)

# Store frames
frames = []
psi_t = psi.copy()

for _ in range(n_frames):
    # Split-step propagation (free particle)
    psi_k = np.fft.fft(psi_t)
    psi_k *= kinetic_prop
    psi_t = np.fft.ifft(psi_k)
    frames.append(np.abs(psi_t)**2)

# Create animation
fig, ax = plt.subplots(figsize=(8, 3), dpi=100)
fig.patch.set_facecolor('#0d1117')  # GitHub dark mode
ax.set_facecolor('#0d1117')

line, = ax.plot([], [], color='#58a6ff', linewidth=2)
fill = ax.fill_between([], [], alpha=0.3, color='#58a6ff')

ax.set_xlim(-15, 15)
ax.set_ylim(0, 0.5)
ax.set_xlabel('Position', color='#c9d1d9', fontsize=12)
ax.set_ylabel('|ψ|²', color='#c9d1d9', fontsize=12)
ax.tick_params(colors='#c9d1d9')
for spine in ax.spines.values():
    spine.set_color('#30363d')

# Title
ax.set_title('Gaussian wavepacket propagation', color='#c9d1d9', fontsize=14, pad=10)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    global fill
    fill.remove()
    y = frames[i]
    line.set_data(x, y)
    fill = ax.fill_between(x, y, alpha=0.3, color='#58a6ff')
    return line, fill

anim = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=n_frames, interval=50, blit=False
)

# Save as GIF
output_path = Path(__file__).parent / 'wavepacket_animation.gif'
anim.save(output_path, writer='pillow', fps=20)
print(f"Saved animation to {output_path}")

plt.close()
