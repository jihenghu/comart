"""
Constants and subscripts for coma property arrays.
"""

# Basic property indices - these are aliases for accessing array elements
IDEN = 0          # Number density
IVEL = 1          # Velocity
ITEMP = 2         # Temperature
IELE= 3           # electron density 
ITE = 4           # Electron temperature

# Number of basic properties per grid point
NPROPS = 5   # Density, velocity, temperature, electron density, electron temperature

# Level population ratio indices (start at LEV0)
LEV0 = 5          # Level 0 population ratio index

NLEVEL = None

ORTHO_RATIO = 0.75  # Ortho-to-para ratio for water, can be adjusted based on observations

# molecular physics constants
# planck constant
_kB=1.38064852e-23  # Boltzmann constant (J/K)
_h=6.62607015e-34   # Planck constant (J*s)
_c=299792458         # Speed of light (m/s)

AU=1.495978707e11    # Astronomical Unit in meters

# level aliases
# for ortho water
J101=0
J110=1
J212=2
J221=3
J303=4
J312=5
J321=6

# for para water
# ...


