"""
Einstein A coefficient tables for level-to-level transitions.

Stores and manages static Einstein A coefficients (A_ij) for transitions
from level i to level j. Can be precomputed from oscillator strengths or
loaded from spectroscopic databases.
"""

from src.error import ComartError
from src.const import *

def ESERates(i, j, r_h=1.0, species='ortho-H2O'):
    """
    Compute the Einstein A coefficient for a transition from level i to level j.

    Parameters:
    - i: from_level index 
    - j: to_level index

    Returns:
    - nu_ij: transition frequency in Hz
    - A_ij: Einstein A coefficient in s^-1
    - G_ij: effective pumping rate due to the solar infrared radiation, in s^-1
    - sigma_ij: molecular collisional cross-section for de-excitation, in m^2

    Raises: level indices out of range [0, NLEVEL-1], or species not recognized
    """
    
    if i == j:
        raise ComartError("Einstein A coefficient is not defined for transitions between the same level.")
    
    if any(level < 0 for level in (i, j)):
        raise ComartError("Energy level indices must be non-negative integers.")
    
    if any(level >= NLEVEL for level in (i, j)):
        raise ComartError(f"Energy level indices must be less than NLEVEL ({NLEVEL}).")
    
    # find table according to species name

    if species == 'ortho-H2O':
        ATable = ORTHO_H2O_COEFFS_TABLE
    else:
        raise ComartError(f"Species '{species}' not recognized. Available species: 'ortho-H2O'...")
    
    if (i, j) in ATable:
        [nu_ij, A_ij, G_ij, sigma_ij] = ATable[(i, j)]
    else:
       [nu_ij, A_ij, G_ij, sigma_ij] = (0, 0, 0, 0)  # Default values for undefined transitions - can be adjusted based on physical considerations
    G_ij = sigma_ij  * r_h**(-2)  # Scale collisional cross-section with heliocentric distance (simplified assumption)

    return nu_ij, A_ij, G_ij, sigma_ij


# Predefined tables for common molecules/atoms
ORTHO_H2O_COEFFS_TABLE = {
    # ortho-H2O molecule transitions (Zakharov et al. 2007)
    # (nu_ij in Hz, A_ij in s^-1, G_ij in s^-1, σ_ij in m^2)
    (J110, J101): (556.936E9, 3.456E-3, 1.423E-5, 2.924E-18),       # Level 1 -> 0
    (J101, J110): (556.936E9,       0., 1.654E-5,        0.),       # Level 0 -> 1
}



