"""
Example: Calculate Einstein coefficients for molecular transitions.

Demonstrates calculation of Einstein A and B coefficients for transitions
between vibrational-rotational states using the coefficient table.

Transitions:
- Level 1 (101): v1=1, v2=0, v3=1
- Level 0 (110): v1=1, v2=1, v3=0
"""

import numpy as np
from src.const import *
from src.error import ComartError
from src import ComaGrid
from src import EinsteinCoeffs
from src import Einstein_Ce_ij
from src import Einstein_Cm_ij


print("=" * 70)
print("Einstein Coefficients for Molecular Transitions")
print("=" * 70)

# Get Einstein coefficients for transition between levels 1 and 0
# EinsteinCoeffs(i, j) returns:
# nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji

try:
    from_level = J101  # Upper level (101)
    to_level = J110    # Lower level (110)
    r_h = 1.0       # Heliocentric distance in AU
    species = 'ortho-H2O'
    
    nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji = EinsteinCoeffs(
        from_level, to_level, r_h=r_h, species=species
    )
    
    print(f"\nTransition: Level {from_level} <-> Level {to_level}")
    print(f"Species: {species}")
    print(f"Heliocentric distance: {r_h} AU")
    
    # Calculate wavelength from frequency
    c = 299792458  # m/s
    wavelength = c / nu_ij
    wavenumber = 1 / wavelength * 100  # Convert to cm^-1
    
    print(f"\nTransition Properties:")
    print(f"  Frequency: {nu_ij:.3e} Hz")
    print(f"  Wavelength: {wavelength*1e6:.2f} μm")
    print(f"  Wavenumber: {wavenumber:.1f} cm^-1")
    
    print(f"\nEinstein A Coefficients:")
    print(f"  A_{from_level}->{to_level} (downward): {A_ij:.3e} s^-1")
    print(f"  A_{to_level}->{from_level} (upward):  {A_ji:.3e} s^-1")
    
    if A_ij > 0:
        print(f"  Lifetime τ = 1/A = {1/A_ij:.3e} s")
    
    print(f"\nEinstein B Coefficients:")
    print(f"  B_{from_level}->{to_level} (absorption): {B_ij:.3e} m^3 J^-1 s^-2")
    print(f"  B_{to_level}->{from_level} (stimulated): {B_ji:.3e} m^3 J^-1 s^-2")
    
    print(f"\nRadiative Pumping Rates (solar infrared):")
    print(f"  G_{from_level}->{to_level}: {G_ij:.3e} s^-1")
    print(f"  G_{to_level}->{from_level}: {G_ji:.3e} s^-1")
    
    print(f"\nCollisional Cross-sections:")
    print(f"  σ_{from_level}->{to_level}: {sigma_ij:.3e} m^2")
    print(f"  σ_{to_level}->{from_level}: {sigma_ji:.3e} m^2")
    
    print(f"\n" + "=" * 70)
    print("Comparison at different heliocentric distances:")
    print("=" * 70)
    
    for r_h_test in [1.0, 2.0, 5.0]:
        nu, A, A_rev, B, B_rev, G, G_rev, s, s_rev = EinsteinCoeffs(
            from_level, to_level, r_h=r_h_test, species=species
        )
        print(f"\nr_h = {r_h_test} AU:")
        print(f"  G_{from_level}->{to_level}: {G:.3e} s^-1")
        print(f"  G_{to_level}->{from_level}: {G_rev:.3e} s^-1")

except ComartError as e:
    print(f"\nError: {e}")
    print("\nNote: This example requires the Einstein coefficient table to be populated.")
    print("Please ensure ORTHO_H2O_COEFFS_TABLE is defined in the constants.")


# Test electron-molecule collision rates (Ce_ij)
print("\n" + "=" * 70)
print("Electron-Molecule Collision Rates (Ce_ij)")
print("=" * 70)

try:
    # Typical coma plasma conditions
    n_e = 1e10  # electron density (m^-3)
    T_e = 500   # electron temperature (K)
    
    from_level = J110
    to_level = J101
    
    Ce_ij = Einstein_Ce_ij(from_level, to_level, n_e, T_e)
    
    print(f"\nElectron-molecule collision rate:")
    print(f"  Transition: Level {from_level} -> Level {to_level}")
    print(f"  Electron density: {n_e:.2e} m^-3")
    print(f"  Electron temperature: {T_e} K")
    print(f"  Ce_ij (excitation rate coeff): {Ce_ij:.3e} m^-3 s^-1")
    
    # Test at different plasma conditions
    print(f"\n" + "-" * 70)
    print("Ce_ij variation with electron density and temperature:")
    print("-" * 70)
    
    for n_e_test in [1e9, 1e10, 1e11]:
        for T_e_test in [300, 500, 1000]:
            Ce = Einstein_Ce_ij(from_level, to_level, n_e_test, T_e_test)
            print(f"n_e={n_e_test:.0e} m^-3, T_e={T_e_test} K: Ce_ij = {Ce:.3e} m^-3 s^-1")

except ComartError as e:
    print(f"\nError calculating Ce_ij: {e}")


try:
    # Typical coma plasma conditions
    n_e = 1e10  # electron density (m^-3)
    T_e = 500   # electron temperature (K)
    
    to_level = J110
    from_level = J101
    
    Ce_ij = Einstein_Ce_ij(from_level, to_level, n_e, T_e)
    
    print(f"\nElectron-molecule collision rate:")
    print(f"  Transition: Level {from_level} -> Level {to_level}")
    print(f"  Electron density: {n_e:.2e} m^-3")
    print(f"  Electron temperature: {T_e} K")
    print(f"  Ce_ij (excitation rate coeff): {Ce_ij:.3e} m^-3 s^-1")
    
    # Test at different plasma conditions
    print(f"\n" + "-" * 70)
    print("Ce_ij variation with electron density and temperature:")
    print("-" * 70)
    
    for n_e_test in [1e9, 1e10, 1e11]:
        for T_e_test in [300, 500, 1000]:
            Ce = Einstein_Ce_ij(from_level, to_level, n_e_test, T_e_test)
            print(f"n_e={n_e_test:.0e} m^-3, T_e={T_e_test} K: Ce_ij = {Ce:.3e} m^-3 s^-1")

except ComartError as e:
    print(f"\nError calculating Ce_ij: {e}")

# Test molecular collision rates (Cm)
print("\n" + "=" * 70)
print("Molecular Collision Coefficients (Cm)")
print("=" * 70)

try:
    # Typical coma neutral atmosphere conditions
    sigma = 2.924e-18    # collision cross-section (m^2) - typical for H2O-H2O
    n_gas = 1e15     # neutral gas density (m^-3)
    T_gas = 100      # gas temperature (K)
    
    Cm = Einstein_Cm_ij(sigma, n_gas, T_gas)
    
    print(f"\nMolecular collision coefficient:")
    print(f"  Collision cross-section: {sigma:.2e} m^2")
    print(f"  Gas density: {n_gas:.2e} m^-3")
    print(f"  Gas temperature: {T_gas} K")
    print(f"  Cm (collision coefficient): {Cm:.3e} m^3 s^-1")
    
    # Test at different gas conditions
    print(f"\n" + "-" * 70)
    print("Cm variation with gas density and temperature:")
    print("-" * 70)
    
    for n_test in [1e14, 1e15, 1e16]:
        for T_test in [50, 100, 200]:
            Cm_val = Einstein_Cm_ij(sigma, n_test, T_test)
            print(f"n_gas={n_test:.0e} m^-3, T_gas={T_test} K: Cm = {Cm_val:.3e} m^3 s^-1")

except ComartError as e:
    print(f"\nError calculating Cm: {e}")
