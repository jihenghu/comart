"""
Einstein A coefficient tables for level-to-level transitions.

Stores and manages static Einstein A coefficients (A_ij) for transitions
from level i to level j. Can be precomputed from oscillator strengths or
loaded from spectroscopic databases.
"""

from src.error import ComartError
from src.const import *
import numpy as np

def EinsteinCoeffs(i, j, r_h=1.0, species='ortho-H2O'):
    """
    Compute the Einstein A coefficient for a transition from level i to level j.

    Parameters:
    - i: from_level index 
    - j: to_level index

    Returns:
    - nu_ij: transition frequency in Hz
    - A_ij: Einstein A coefficient in s^-1
    - A_ji: Einstein A coefficient for the reverse transition in s^-1
    - B_ij: Einstein B coefficient in m^3 J^-1 s^-2
    - B_ji: Einstein B coefficient for the reverse transition in m^3 J^-1 s^-2
    - G_ij: effective pumping rate due to the solar infrared radiation, in s^-1
    - G_ji: effective pumping rate for the reverse transition due to the solar infrared radiation, in s^-1
    - sigma_ij: molecular collisional cross-section for de-excitation, in m^2

    Raises: level indices out of range [0, NLEVEL-1], or species not recognized
    """
    
    if i == j:
        raise ComartError("Einstein A coefficient is not defined for transitions between the same level.")
    
    if any(np.array([i, j]) < 0):
        raise ComartError("Energy level indices must be non-negative integers.")
    
    if any(np.array([i, j]) >= len(GFACTORS)):
        raise ComartError(f"Energy level indices must be less than NLEVEL ({len(GFACTORS)}).")
    
    # find table according to species name

    if species == 'ortho-H2O':
        ATable = ORTHO_H2O_COEFFS_TABLE
    else:
        raise ComartError(f"Species '{species}' not recognized. Available species: 'ortho-H2O'...")
    
    if (i, j) in ATable:
        [nu_ij, A_ij, G_ij, sigma_ij] = ATable[(i, j)]
        [   _ , A_ji, G_ji, sigma_ji] = ATable[(j, i)]
    else:
        raise ComartError(f"Transition from level {i} to level {j} not found in the Einstein coefficients table for species '{species}'.")
        # Default values for undefined transitions - can be adjusted based on physical considerations

    G_ij = G_ij  * (r_h**(-2))  # Scale collisional cross-section with heliocentric distance (simplified assumption)
    G_ji = G_ji  * (r_h**(-2))  # Scale collisional cross-section with heliocentric distance (simplified assumption)

    # B12 and B21 can be computed from A_ij using Einstein relations

    # Fv = 2 * PLANKC * (nu_ij**3) / (LIGHTSPEED**2)  # Spectral energy density of the radiation field at frequency nu_ij (simplified as blackbody radiation)

    if i>j: # then, Aij is nonzero, downward transition, and we can compute both B_ij and B_ji
        # Downward transition
        B_ij = A_ij * LIGHTSPEED**2 / (2.0 * PLANKC * nu_ij**3)  # Einstein B coefficient for stimulated emission
        B_ji = B_ij * GFACTORS[i]  / GFACTORS[j]  # Einstein B coefficient for absorption

    else: # then, Aij is zero, upward transition, and we can only compute B_ij from B_ji using Einstein relations 
        # Upward transition
        B_ji = A_ji * LIGHTSPEED**2 / (2.0 * PLANKC * nu_ij**3)  # Einstein B coefficient for stimulated emission (from downward transition)
        B_ij = B_ji * GFACTORS[j] / GFACTORS[i]  # Einstein B coefficient for absorption (from downward transition)

    return nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji






