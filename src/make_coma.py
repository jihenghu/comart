"""
Factory module for creating and initializing ComaGrid instances.

Provides functions to create ComaGrid objects from input arrays
of physical properties (density, velocity, temperature, level ratios).
"""

import numpy as np
from src import const
from src.coma import ComaGrid
from src.error import ComartError



def make_coma_from_arrays(xb, density, velocity, temperature, electron_density=None, electron_temperature=None):
    """
    Create and initialize a ComaGrid from input property arrays.
    
    Args:
        xb: 1D array of grid boundary values (length nx+1)
        density: 1D array of number density values (length nx)
        velocity: 1D array of velocity values (length nx)
        temperature: 1D array of temperature values (length nx)
        electron_density: 1D array of electron density values (length nx), optional
        electron_temperature: 1D array of electron temperature values (length nx), optional
    
    Returns:
        ComaGrid: Initialized grid with properties set from input arrays
    
    Raises:
        ComartError: If input array dimensions don't match nx
    """

    nx = len(density)
    # Create empty grid
    grid = ComaGrid(nx)
    
    # Set grid boundaries
    if len(xb) != nx + 1:
        raise ComartError(f"xb array size {len(xb)} does not match nx+1={nx+1}")
    grid.x[:] = np.asarray(xb, dtype=np.float64)
    grid.xc[:] = 0.5 * (grid.x[:-1] + grid.x[1:])
    
    # Validate input array shapes
    if len(velocity) != nx:
        raise ComartError(f"velocity array size {len(velocity)} does not match nx={nx}")
    if len(temperature) != nx:
        raise ComartError(f"temperature array size {len(temperature)} does not match nx={nx}")
    if electron_density is not None and len(electron_density) != nx:
        raise ComartError(f"electron_density array size {len(electron_density)} does not match nx={nx}")
    if electron_temperature is not None and len(electron_temperature) != nx:
        raise ComartError(f"electron_temperature array size {len(electron_temperature)} does not match nx={nx}")
    
    # Set properties
    grid.props[const.IDEN, :] = np.asarray(density, dtype=np.float64)
    grid.props[const.IVEL, :] = np.asarray(velocity, dtype=np.float64)
    grid.props[const.ITEMP, :] = np.asarray(temperature, dtype=np.float64)
    grid.init_level_populations()  # Initialize level populations based on density

    if electron_density is not None:
        grid.props[const.IELE, :] = np.asarray(electron_density, dtype=np.float64)
    if electron_temperature is not None:
        grid.props[const.TELE, :] = np.asarray(electron_temperature, dtype=np.float64)
 
    return grid
