"""
Simple example demonstrating coma grid creation with model profiles.

Creates a 100-layer coma using velocity, temperature, and density models.
Plots the resulting profiles from 1 km to 1e6 km with logarithmic spacing.
"""

import numpy as np
import matplotlib.pyplot as plt
from src import ComaGrid
# from src import make_coma_from_arrays
from src.profiles.velocity_model import *
from src.profiles.temperature_model import *
from src.profiles.density_model import *
from src.profiles.electron import *
from src import *


# Grid parameters
nx = 100
r_min = 1e3      # 1 km in meters
r_max = 1e9      # 1e6 km in meters

# Create logarithmic grid for cell centers
xc = np.logspace(np.log10(r_min), np.log10(r_max), nx)

# Create grid boundaries (nx+1 points)
xb = np.zeros(nx + 1)
xb[0] = r_min
xb[-1] = r_max
# Interpolate boundaries between cell centers
for i in range(1, nx):
    xb[i] = np.sqrt(xc[i-1] * xc[i])

# Compute model profiles on grid cell centers
print(f"Computing profiles for {nx} layers from {r_min/1e3:.1e} km to {r_max/1e3:.1e} km")

# Velocity profile: using tangent model
# v(r) = v0 * tanh(r / c)
v0 = 720.0          # terminal velocity in m/s
c_scale = 3.5e3      # scale parameter in meters
velocity = velocity_tangent(xc, v0=v0, c=c_scale)

# Temperature profile: inverse distance dependence
# T(r) = T0 * (a * r0 / r + b)
T0 = 10.0          # reference temperature in K
temperature = temperature_revesedist(xc, a=7, b=1, T0=T0, r0=2000)

# Density profile: isotropic radial outflow
# n(r) = Q / (4π r² v)
Q = 1e26            # production rate (molecules/s)
density = calc_density(Q, xc, velocity, r_h=2.5, model_type='isotropic', photo_dissociation=True)


# electron density and temperature profiles
nelec, telec = electron_Biver1997(Q, r_h=2.5, r=xc, v=velocity, Tkin=temperature, Temax=10000)

# Create ComaGrid with computed profiles
grid = ComaGrid(nx, nlevel=2, density=density, velocity=velocity, temperature=temperature, electron_density=nelec, electron_temperature=telec)

grid.x = xb
grid.xc = xc


print(f"\nComaGrid created:")
print(f"  Shape: {grid.props.shape}")
print(f"  Memory: {grid.nbytes / 1e3:.2f} kB")
print(f"  Contiguous: {grid.is_contiguous}")

# Create plots
fig, axes = plt.subplots(2,3, figsize=(12, 8), dpi=300)

axes = axes.flatten()

# Velocity plot
ax = axes[2]
ax.plot(xc / 1e3, grid.props[IVEL, :])
ax.set_xscale('log')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Velocity (m/s)')
ax.set_title('Velocity Profile')
ax.grid(True, alpha=0.3)

# Density plot
ax = axes[0]
ax.loglog(xc / 1e3, grid.props[IDEN, :], 'k-', label='total')
for i in range(grid.nlevel):
    ax.plot(xc / 1e3, grid.get_level_populations()[i, :], color=['r','b','g','y'][i], linestyle='--', label=f'L{i}')

for ix in range(grid.nx):
    grid.make_level_transition(ix, from_level=0, to_level=1, ratio=0.5)
for i in range(grid.nlevel):
    ax.plot(xc / 1e3, grid.get_level_populations()[i, :], color=['r','b','g','y'][i], linestyle='--', label=f"L{i}'")

ax.set_ylim(1e3,1e17)
ax.legend()
ax.set_xlabel('Distance (km)')
ax.set_ylabel('molecules (m$^{-3}$)')
ax.set_title('Density Profile (Isotropic)')
ax.grid(True, alpha=0.3)

# Temperature plot
ax = axes[1]
ax.plot(xc / 1e3, grid.props[ITEMP, :])
ax.set_xscale('log')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Temperature (K)')
ax.set_title('Temperature Profile')
ax.grid(True, alpha=0.3)



# Electron density plot
ax = axes[3]
ax.loglog(xc / 1e3, grid.props[IELE, :])
ax.set_xscale('log')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Electron Density (m$^{-3}$)')
ax.set_title('Electron Density Profile')
ax.grid(True, alpha=0.3)

# Electron temperature plot
ax = axes[4]
ax.loglog(xc / 1e3, grid.props[ITE, :])
ax.set_xscale('log')
ax.set_xlabel('Distance (km)')
ax.set_ylabel('Electron Temperature (K)')
ax.set_title('Electron Temperature Profile')
ax.grid(True, alpha=0.3)    



plt.tight_layout()
plt.savefig('coma_profiles.png', dpi=150)
print("\nPlot saved to: coma_profiles.png")
# plt.show()

