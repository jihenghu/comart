"""
Water gas density model for isotropic radial outflow in cometary comas.

Implements steady-state spherical expansion model for water molecules (H2O).
"""

import numpy as np
from src.error import ComartError

def _isotropic(q, r, v):
    """
    Calculate water (H2O) number density at distance r in isotropic radial outflow.
    
    Uses the continuity equation for steady-state expansion:
    n(r) = Q / (4π r² v)
    
    This assumes:
    - Spherical symmetry
    - Steady-state outflow
    - Constant expansion velocity
    - No recombination/destruction
    
    Parameters
    ----------
    q_h2o : float
        Total water production rate from nucleus [molecules/s]
    r : float
        Distance from nucleus center [m]
    v : float
        Radial expansion velocity [m/s]
    
    Returns
    -------
    float
        Number density of water molecules [m^-3]
    
    Raises
    ------
    ComartError
        If inputs are invalid (non-positive values, etc.)
    
    Notes
    -----
    The formula is derived from mass conservation in spherical coordinates:
    
    .. math::
        4 \\pi r^2 v n(r) = Q
        
    where Q is the constant outflow rate (production rate at equilibrium).
    """
    if np.any(q <= 0):
        raise ComartError(f"Production rate must be positive, got {q}")
    
    if np.any(r <= 0):
        raise ComartError(f"Distance must be positive, got {r} m")
    
    if np.any(v <= 0):
        raise ComartError(f"Expansion velocity must be positive, got {v} m/s")
    
    n = q / (4.0 * np.pi * r**2 * v)
    
    return n


def _jet(q, r, v, theta=0.0):
    """
    Calculate water density in a jet outflow model.
    
    This model assumes a collimated jet with a Gaussian angular distribution:
    n(r, θ) = (Q / (4π r² v)) * exp(-θ² / (2σ²))
    
    where σ controls the jet width and θ is the angle from the jet axis.
    """    
    
    pass


def _seasonal(q, r, v, season_factor=1.0):
    """
    Calculate water density in a seasonal outflow model.
    
    This model assumes that the production rate varies with season, which can
    be represented by a time-dependent factor:
    n(r, t) = (Q * season_factor(t) / (4π r² v))
    """
    pass

def calc_density(q, r, v, r_h=1.0, model_type='isotropic', photo_dissociation=True,  **kwargs):
    """
    General outgassing model selector.
    
    Parameters
    ----------
    model_type : str
        Type of outgassing model ('isotropic', 'jet', 'seasonal')
    q : float
        Total water production rate from nucleus [molecules/s]
    r_h : float
        Heliocentric distance [AU]
    r : float
        Distance from nucleus center [m]
    v : float
        Radial expansion velocity [m/s]
    photo_dissociation : bool
        Whether to include photodissociation effects in the coma modeling
    **kwargs
        Additional parameters specific to the chosen model (e.g., theta for jet, season_factor for seasonal)
    
    Returns
    -------
    float
        Number density of water molecules based on the selected outgassing model [m^-3]
    
    Raises
    ------
    ComartError
        If model_type is invalid or if inputs are invalid for the chosen model.
    
    Notes
    -----
    This function serves as a unified interface to select different outgassing models based on user input.
    It validates the model type and passes the appropriate parameters to the corresponding model function.
    """
    if photo_dissociation:
        beta = kwargs.get('beta', 1.042e-5)  # Default photodissociation rate at 1 AU  (Haser et al. 1957)
        dissoci_rate = np.exp(-beta / (r_h**2) * r / v)
    else:
        dissoci_rate = 1.0

    if model_type == 'isotropic':
        return _isotropic(q, r, v) * dissoci_rate
    
    elif model_type == 'jet':
        theta = kwargs.get('theta', 0.0)
        return _jet(q, r, v, theta) * dissoci_rate
    
    elif model_type == 'seasonal':
        season_factor = kwargs.get('season_factor', 1.0)
        return _seasonal(q, r, v, season_factor) * dissoci_rate
    
    else:
        raise ComartError(f"Invalid outgassing model '{model_type}'. Must be one of: 'isotropic', 'jet', 'seasonal'.")  