"""
Temperature profile models for comet coma.

Provides functions to compute temperature profiles as a function of distance
from the comet nucleus.
"""

import numpy as np
from src import const
from src.error import ComartError

def temperature_constant(r, T):
    """
    Constant temperature profile.
    
    Args:
        r: Array of distances from comet nucleus
        T: Constant temperature value
    
    Returns:
        Array of temperature values at given distances
    """
    return np.ones_like(r) * T



def temperature_revesedist(r, a=1, b=0, T0=150, r0=2e3):
    """
    Temperature profile for cometary coma (inverse distance dependence).
    
    The temperature profile is given by:
        T(r) = T0 * (a * r0 / r + b)
    
    Parameters
    ----------
    r : float or ndarray
        Radial distance from nucleus center in meters.
        Can be scalar or array. Must be positive.
    a : float, optional
        Scaling factor for the 1/r dependence of the temperature profile. Default is 1.
    b : float, optional
        Scaling factor for the terminal (asymptotic) temperature value
        as distance increases. Default is 0.
    T0 : float, optional
        Reference temperature constant in Kelvin. Default is 150 K.
        Must be positive.

        TODO: should be function of heliocentric distance, or set by user?

    r0 : float, optional
        Reference distance constant in meters. Default is 2 km (2e3 m).
        Must be positive.
    
    Returns
    -------
    float or ndarray
        Temperature at distance r in Kelvin. Same shape as input r.
    
    Raises
    ------
    ComartError
        If T0 <= 0, r0 <= 0, or if r contains non-positive values
    
    """
    
    # Validate inputs
    if T0 <= 0:
        raise ComartError(f"Reference temperature T0 must be positive, got {T0}")
    if r0 <= 0:
        raise ComartError(f"Reference distance r0 must be positive, got {r0}")
    
    # Convert r to numpy array for uniform handling
    r_arr = np.atleast_1d(r)
    
    # Check for non-positive distances
    if np.any(r_arr <= 0):
        raise ComartError(f"Distance r must be positive, got non-positive values")
    
    # Calculate temperature profile
    temperature = T0 * (a * r0 / r_arr + b)
    
    # Return scalar if input was scalar
    if np.isscalar(r):
        return float(temperature[0])
    
    return temperature
