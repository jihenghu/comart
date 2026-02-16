"""
Molecular collision rates and coefficients.

Calculates collision rates and coefficients (Cm) for molecular collisions
in a gas. Used for computing collisional de-excitation rates and thermalization.

References: 
    Buffa, et al., 2000, H2O-H2O Collision Rate Coefficients, ApJS, 128 597

"""

import numpy as np
from src.error import ComartError



def thermal_velocity(temperature, mass1=18.01528, mass2=18.01528):
    """
    Calculate mean thermal velocity of gas molecules.
    
    v_th = sqrt(8 * k_B * T / (π * m))
    
    Args:
        temperature: Gas temperature (K)
        mass1: Mass of first molecule (amu, default: 18.01528 for H2O)
        mass2: Mass of second molecule (amu, default: 18.01528 for H2O)
    
    Returns:
        Mean thermal velocity (m/s)
    
    Raises:
        ComartError: If inputs are invalid
    """
    if temperature < 0:
        raise ComartError("Temperature must be non-negative")
    
    if mass1 <= 0 or mass2 <= 0:
        raise ComartError("Masses must be positive")

    # reduced mass for H2O-H2O collisions 
    mu = (mass1 * mass2) / (mass1 + mass2) * 1.66053906660e-27 
    
    k_B = 1.38064852e-23  # Boltzmann constant (J/K)
    
    if temperature == 0:
        return 0.0
    
    v_th = np.sqrt(8 * k_B * temperature / (np.pi * mu))
    return v_th

def Einstein_Cm_ij(sigma_ij, n_gas, temperature):
    """
    Calculate collision coefficient Cm for molecular collisions.
    
    Combines collision cross-section, gas density, and thermal velocity.
    This is the effective collision coefficient used in statistical equilibrium equations.
    
    Cm = σ * sqrt(8 * k_B * T / (π * m)) * n_gas
    
    Args:
        sigma_ij: Collision cross-section (m^2)
        n_gas: Gas number density (m^-3)
        temperature: Gas temperature (K)
    
    Returns:
        Collision coefficient Cm (m^3 s^-1)
    
    Raises:
        ComartError: If inputs are invalid
    """
    if sigma_ij < 0:
        raise ComartError("Collision cross-section must be non-negative")
    if n_gas < 0:
        raise ComartError("Gas density must be non-negative")
    if temperature < 0:
        raise ComartError("Temperature must be non-negative")
    
    v_th = thermal_velocity(temperature)
    Cm_ij = sigma_ij * v_th * n_gas
    return Cm_ij





