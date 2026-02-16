"""
Electron-molecule collision rates and coefficients.

Calculates collision rates and coefficients for electron-molecule collisions.
Electrons are much lighter than neutral molecules, leading to different collision
dynamics compared to neutral-neutral collisions.
"""

import numpy as np
from src.error import ComartError
import src.const as CONST
from src.ese.Einstein_coeffs import EinsteinCoeffs

from scipy.special import k0
# Physical constants



def calc_sigmae_ij(Aij, v_ij):
    """
    Calculate electron collision cross-section σ^e_ij from Einstein A coefficient.
    
    """

    if Aij <= 0:
        raise ComartError("Einstein A coefficient must be positive to calculate collision cross-section.")
    
    if v_ij <= 0:
        raise ComartError("Transition frequency must be positive to calculate collision cross-section.")

    determinator=16. * np.pi*np.pi * (CONST.PLANKC**2) * (v_ij**4) * CONST.EPSILON0

    Sigma_ij=CONST.ELECTRONMASS*(CONST.ELEMCHARGE**2)* (CONST.LIGHTSPEED**3) * Aij / determinator

    return Sigma_ij


def Einstein_Ce_ij(i,j, ne, te):
    """Calculate collisional excitation/de-excitation rate coefficient Ce_ij for electron collisions.
       param: 
            i: from energy level index
            j: to energy level index
            ne: electron number density (m^-3)
            te: electron temperature (K)
        reference:
            Zarkharov, A. M., 2007

    """
    if i==j:
        raise ComartError("Collision coefficient is not defined for transitions between the same level.")
    
    if any(np.array([i, j]) < 0):
        raise ComartError("Energy level indices must be non-negative integers.")
    
    if any(np.array([i, j]) >= len(CONST.GFACTORS)):
        raise ComartError(f"Energy level indices must be less than NLEVEL ({len(CONST.GFACTORS)}).")
    


    if i>j:    # de-excitation transition 
        A_ij=EinsteinCoeffs(i,j)[1] # A_ij
        v_ij=EinsteinCoeffs(i,j)[0] # nu_ij

        aij=CONST.PLANKC * v_ij/(2*CONST.BOLTZMANNC*te) # 

        # elec themal velocity
        ve= np.sqrt(8*CONST.BOLTZMANNC*te/(np.pi*CONST.ELECTRONMASS))
        
        # collisional cross-section for electron - water
        sigma_e_ij=calc_sigmae_ij(A_ij, v_ij)

        # collisional excitation rate coefficient for electron-water collisions
        Ce_ij= ne* ve * sigma_e_ij * 2. * aij * np.exp(aij) * k0(aij) 
        

    else:    # excitation transition
        A_ji=EinsteinCoeffs(j,i)[1] # A_ji
        v_ji=EinsteinCoeffs(j,i)[0] # nu_ji

        aji=CONST.PLANKC * v_ji/(2*CONST.BOLTZMANNC*te) # 
        ve= np.sqrt(8*CONST.BOLTZMANNC*te/(np.pi*CONST.ELECTRONMASS))

        sigma_e_ji=calc_sigmae_ij(A_ji, v_ji)

        Ce_ij= ne* ve * CONST.GFACTORS[j]/CONST.GFACTORS[i] * sigma_e_ji * 2. * aji * np.exp(-aji) * k0(aji)

    return Ce_ij