import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from matplotlib.animation import FuncAnimation

# Load the spherical coefficients data
file_path = '20250104_221829_simulation_data/displacement101.txt'

# Read data from the file
with open(file_path, 'r') as f:
    lines = f.readlines()

# Parse data into time steps and coefficients
data = {}
current_time = None
for line in lines:
    values = line.split()
    if len(values) == 1:  # Time step
        current_time = float(values[0])
        data[current_time] = []
    elif len(values) == 4:  # Coefficient data
        j, m, re, im = int(values[0]), int(values[1]), float(values[2]), float(values[3])
        data[current_time].append((j, m, re, im))

# INPUT PARAMETERS
l_max = 5  # Maximum degree of spherical harmonics
n_theta, n_phi = 100, 200  # Grid resolution
time_steps = len(data)  # Number of time steps

# Prepare coefficients array for animation
C = np.zeros((l_max + 1, 2 * l_max + 1, time_steps), dtype=complex)  # Complex array

# Process coefficients from data
times = sorted(data.keys())
for t_idx, time in enumerate(times):
    for j, m, re, im in data[time]:
        if j <= l_max and m >= 0:  # Only process m >= 0
            C[j, l_max + m, t_idx] = complex(re, im)  # Store complex coefficient

# Reconstruct coefficients for m < 0 using symmetry
for l in range(l_max + 1):
    for m in range(1, l + 1):  # Skip m = 0
        for t_idx in range(time_steps):
            # C[l, -m] = (-1)^m * conj(C[l, m])
            C[l, l_max - m, t_idx] = (-1) ** m * np.conj(C[l, l_max + m, t_idx])

# GRID CREATION
theta = np.linspace(0, np.pi, n_theta)
phi = np.linspace(0, 2 * np.pi, n_phi)
theta, phi = np.meshgrid(theta, phi)

# Longitude and latitude conversions
lon = phi - np.pi  # Longitude
lat = np.pi / 2 - theta  # Latitude

# Compute global min and max for consistent color scaling
field_min, field_max = float('inf'), float('-inf')

# Calculate the min and max across all timesteps
for t_idx in range(time_steps):
    field = np.zeros_like(theta, dtype=float)
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            Y_lm = sph_harm(m, l, phi, theta)
            field += (C[l, l_max + m, t_idx] * Y_lm).real  # Real part of harmonics

    # Update global min and max
    field_min = min(field_min, np.min(field))
    field_max = max(field_max, np.max(field))

print(f"Global Min: {field_min}, Global Max: {field_max}")

# Initialize plot for one layer with fixed color scale
fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'mollweide'})
plt.subplots_adjust(wspace=0.3)  # Add spacing

# Create color mesh for the single layer with fixed vmin and vmax
c_mesh = ax.pcolormesh(lon, lat, np.zeros_like(theta), cmap='inferno', 
                       shading='auto', vmin=field_min, vmax=field_max)
plt.colorbar(c_mesh, ax=ax, label='Intensity')
ax.set_title('Spherical Harmonics Distribution')
ax.grid(True)

# UPDATE FUNCTION FOR ANIMATION
def update(t):
    # Compute field for time step t
    field = np.zeros_like(theta, dtype=float)
    for l in range(l_max + 1):
        for m in range(-l, l + 1):  # Loop over full m range (-l to l)
            Y_lm = sph_harm(m, l, phi, theta)  # Compute spherical harmonic
            field += (C[l, l_max + m, t] * Y_lm).real  # Use the real part only

    # Update the plot
    c_mesh.set_array(field.ravel())
    ax.set_title(f'Time: {times[t]:.2f}s')
    return c_mesh,

# CREATE ANIMATION
ani = FuncAnimation(fig, update, frames=time_steps, interval=500, blit=False)

# DISPLAY THE ANIMATION
plt.show()

# Uncomment this line if you want to save the animation to a video file
ani.save('spherical_harmonics_animation.mp4', fps=40, writer='ffmpeg')
