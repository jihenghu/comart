"""
Velocity profile models for comet coma.

Provides functions to compute velocity profiles as a function of distance
from the comet nucleus.
"""

import numpy as np
from src.error import ComartError

def velocity_tangent(r, v0=800, c=1e5):
    """
    Velocity profile for cometary coma expansion (hyperbolic tangent-like).
    
    References:
        S. Lee et al. (2015). A&A 583, A5.
        Spatial and diurnal variation of water outgassing on comet 67P/Churyumov-Gerasimenko observed from Rosetta/MIRO in August 2014
    
    The velocity profile is given by:
        v(r) = v0 * tanh(r / c)

    Parameters
    ----------
    r : float or ndarray
        Radial distance from nucleus center in meters.
        Can be scalar or array.
    v0 : float, optional
        Terminal (asymptotic) velocity at large distances in m/s.
        This is the maximum expansion velocity of the coma.
    c : float, optional
        Scale parameter controlling the velocity transition rate in meters.
        Smaller c means faster initial velocity increase.
        Must be positive.
    
    Returns
    -------
    float or ndarray
        Radial velocity at distance r in m/s. Same shape as input r.
    
    Raises
    ------
    ComartError
        If v0 < 0, c <= 0, or if r contains negative values
    
    """
    
    # Validate inputs
    if v0 < 0:
        raise ComartError(f"Terminal velocity v0 must be non-negative, got {v0}")
    if c <= 0:
        raise ComartError(f"Scale parameter c must be positive, got {c}")
    
    # Convert r to numpy array for uniform handling
    r_arr = np.atleast_1d(r)
    
    # Check for negative distances
    if np.any(r_arr < 0):
        raise ComartError(f"Distance r must be non-negative, got negative values")
    
    # Calculate velocity profile
    velocity = v0 * np.tanh(r_arr / c)
    
    # Return scalar if input was scalar
    if np.isscalar(r):
        return float(velocity[0])
    
    return velocity


def velocity_power_law(distance, v0, d0, alpha):
    pass


def velocity_thermal(temperature, mass):
    pass
