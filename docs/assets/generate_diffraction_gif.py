#!/usr/bin/env python3
"""Generate a diffraction grating GIF for the README."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numexpr as ne
from pathlib import Path

# Use legacy API for compatibility with notebook example
import pytalises as pt

# Create wavefunction - 2D Gaussian with momentum
psi = pt.legacy.Wavefunction(
    "exp(-((x-x0)/sigmax)**2)*exp(-((y-y0)/sigmay)**2)*exp(1j*ky*y)", 
    variables={'x0': 0, 'y0': -6, 'sigmax': 4, 'sigmay': 1, 'ky': 4},
    number_of_grid_points=(128, 320),
    spatial_ext=[(-10, 10), (-12, 12)],
)

# Periodic grating potential
v = "where(y<.2, 1, 0)*where(y>-.2, 1, 0)*where(cos(3*x)<0, 1, 0)*1000"
potential = ne.evaluate(v, local_dict=psi.default_var_dict)[:, :, 0]

# Create figure with dark theme
fig, ax = plt.subplots(figsize=(6, 6), dpi=80)
fig.patch.set_facecolor('#0d1117')
ax.set_facecolor('#0d1117')

# Initial plot
data = (np.abs(psi.amp**2) + potential * 0.0001).T
im = ax.imshow(
    data,
    origin='lower',
    extent=[-10, 10, -12, 12],
    vmax=np.max(np.abs(psi.amp**2)),
    cmap='inferno',
    aspect='equal',
)

ax.set_xlabel('x', color='#c9d1d9', fontsize=11)
ax.set_ylabel('y', color='#c9d1d9', fontsize=11)
ax.tick_params(colors='#c9d1d9')
for spine in ax.spines.values():
    spine.set_color('#30363d')

# Add subtle grating overlay
grating_vis = np.where(potential > 0, 0.3, 0)
ax.imshow(
    grating_vis.T,
    origin='lower',
    extent=[-10, 10, -12, 12],
    cmap='Greys',
    alpha=0.5,
    aspect='equal',
)

ax.set_title('Diffraction on periodic grating', color='#c9d1d9', fontsize=12, pad=10)

frames_data = []

# Pre-compute frames
n_frames = 200
for i in range(n_frames):
    psi.propagate(v, num_time_steps=6, delta_t=0.003, diag=True)
    frame = (np.abs(psi.amp**2) + potential * 0.00005).T.copy()
    frames_data.append(frame)
    if i % 30 == 0:
        print(f"Computing frame {i}/{n_frames}")

def animate(i):
    im.set_data(frames_data[i])
    return [im]

anim = animation.FuncAnimation(
    fig, animate, frames=n_frames, interval=40, blit=True
)

# Save as GIF
output_path = Path(__file__).parent / 'diffraction_grating.gif'
print("Saving animation...")
anim.save(output_path, writer='pillow', fps=25)
print(f"Saved animation to {output_path}")

plt.close()
