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


ORTHO_RATIO = 0.75  # Ortho-to-para ratio for water, can be adjusted based on observations

# molecular physics constants
# planck constant
BOLTZMANNC=1.38064852e-23  # Boltzmann constant (J/K)
PLANKC=6.62607015e-34   # Planck constant (J*s)
LIGHTSPEED=299792458         # Speed of light (m/s)

ELECTRONMASS=9.10938356e-31 # Electron mass (kg)
ELEMCHARGE=1.602176634e-19  # Elementary charge (C)
EPSILON0=8.854187817e-12  # Vacuum permittivity (F/m)




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

GFACTORS={
    J101: 9,
    J110: 9,
    J212: 15,
    J221: 15,
    J303: 21,
    J312: 21,
    J321: 21
}

# Predefined tables for common molecules/atoms
ORTHO_H2O_COEFFS_TABLE = {
    # ortho-H2O molecule transitions (Zakharov et al. 2007)
    # (nu_ij in Hz, A_ij in s^-1, G_ij in s^-1, σ_ij in m^2)
    (J110, J101): (556.936E9, 3.456E-3, 1.423E-5, 2.924E-18),       # Level 1 -> 0
    (J101, J110): (556.936E9,       0., 1.654E-5,        0.),       # Level 0 -> 1
}


# for para water
# ...


