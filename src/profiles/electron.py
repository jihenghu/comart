"""
Electron density and temperature profile models for comet coma.

Provides functions to compute electron density and electron temperature
profiles as a function of distance from the comet nucleus.
Electrons are produced by photoionization and impact ionization.
"""

import numpy as np
from src.error import ComartError

def electron_Biver1997(Q, r_h, r, v, Tkin, Temax=10000.):

    """Electron density profile model based on Biver et al. (1997) for comet 1P/Halley.
    
    Args:
        Q: Water production rate (molecules/s) at surface of nucleus
        r_h: Heliocentric distance (AU)
        r: Array of distances from comet nucleus (m)
        v: Gas expansion velocity (m/s)
        Tkin: Kinetic temperature of gas (K)
        Temax: Maximum electron temperature (K)
    
    Returns:
        Array of electron density values (m^-3)
        Array of electron temperature values (K)

    Raises:
        ComartError: If inputs are invalid
    """

    if r_h <= 0:
        raise ComartError("Heliocentric distance r_h must be positive")
    
    if Q < 0:
        raise ComartError("Surface water production rate Q must be non-negative")
    
    if np.any(r < 0):
        raise ComartError("Distance from comet nucleus r must be non-negative")
    
    if np.any(v <= 0):
        raise ComartError("Gas expansion velocity v must be positive")
    
    if np.any(Tkin < 0):
        raise ComartError("Kinetic temperature Tkin must be non-negative")
    
    if Temax <= 0:
        raise ComartError("Maximum electron temperature Temax must be positive")
    
    # Q29, water production rate scaled by heliocentric distance
    Q29 = 1e-29 * Q # Water production rate (molecules/s) scaled by heliocentric distance

    # contact surface distance (m)
    R_cs = 1.125E6 * (Q29**0.75)

    # recombination zone boundary distance (m)
    R_rec = 3.2E6 * (Q29**0.5)



    # telec, electron temperature profile (K)
    telec=np.zeros_like(r)

    index = r < R_cs
    telec[index] = Tkin[index]

    index = (r >= R_cs) & (r < 2*R_cs)
    telec[index] = Tkin[index] + (Temax - Tkin[index]) * (r[index] - R_cs) / R_cs

    index = r >= 2*R_cs
    telec[index] = Temax



    ###
    
    # K_ion, water photoionization rate (s^-1)
    K_ion = 4.1E-7 * (r_h**(-2))

    # K_rec, electron recombination rate with ions H3O+ (m^-3 s-1)
    K_rec = 7.0E-13 * ((300. / telec)**0.5)

    # n_e, electron density (m^-3)
    term1= Q* K_ion / (v* K_rec * (r_h**2))
    term2= telec / 300.
    term3= (R_rec / (r**2))
    term4= 1- np.exp(-r / R_rec)
    term5= 5.0E6*(r_h**(-2))

    nelec = (term1**0.5) * (term2**0.15) * term3 * term4 + term5

    return nelec, telec


