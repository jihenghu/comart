"""
Example: Test Einstein coefficients and collision rates on ComaGrid.

Demonstrates how to use Einstein coefficients, molecular collision rates (Cm),
and electron collision rates (Ce_ij) on a simulated comet coma grid.
Creates profiles and computes collision coefficients at different grid points.
"""

import numpy as np
import matplotlib.pyplot as plt
from src.const import *
from src.error import ComartError
from src.coma import ComaGrid
from src.profiles.velocity_model import velocity_tangent
from src.profiles.temperature_model import temperature_revesedist
from src.profiles.density_model import calc_density
from src.profiles.electron import electron_Biver1997

from src.ese.Einstein_coeffs import EinsteinCoeffs
from src.ese.molec_elec_collision import Einstein_Ce_ij
from src.ese.molecular_collision import Einstein_Cm_ij


print("=" * 80)
print("Einstein Coefficients + Collision Rates on ComaGrid")
print("=" * 80)

try:
    # Create a small grid for testing
    nx = 20
    r_min = 1e3      # 1 km
    r_max = 1e8      # 1E5 km
    
    # Grid boundaries and centers
    xc = np.logspace(np.log10(r_min), np.log10(r_max), nx)
    xb = np.zeros(nx + 1)
    xb[0] = r_min
    xb[-1] = r_max
    for i in range(1, nx):
        xb[i] = np.sqrt(xc[i-1] * xc[i])
    
    print(f"\nCreating ComaGrid with {nx} layers")
    print(f"  Distance range: {r_min/1e3:.1e} - {r_max/1e3:.1e} km")
    
    # Compute profiles
    print("\nComputing neutral and electron profiles...")
    v0 = 720.0
    c_scale = 3.5e3
    velocity = velocity_tangent(xc, v0=v0, c=c_scale)
    
    T0 = 10.0
    temperature = temperature_revesedist(xc, a=7, b=1, T0=T0, r0=2000)
    
    Q = 1e26
    density = calc_density(Q, xc, velocity, r_h=2.5, model_type='isotropic', photo_dissociation=True)
    
    # Electron profiles
    nelec, telec = electron_Biver1997(Q, r_h=2.5, r=xc, v=velocity, Tkin=temperature, Temax=10000)
    
    # Create grid with profiles
    grid = ComaGrid(nx, nlevel=2, density=density, velocity=velocity, 
                    temperature=temperature, electron_density=nelec, 
                    electron_temperature=telec)
    grid.x = xb
    grid.xc = xc
    
    print(f"\n✓ ComaGrid created successfully")
    print(f"  Neutral density range: {density.min():.2e} - {density.max():.2e} m⁻³")
    print(f"  Temperature range: {temperature.min():.1f} - {temperature.max():.1f} K")
    print(f"  Electron density range: {nelec.min():.2e} - {nelec.max():.2e} m⁻³")
    print(f"  Electron temperature range: {telec.min():.1f} - {telec.max():.1f} K")
    
    # Test coefficients at selected grid points
    print(f"\n" + "=" * 80)
    print("Collision Coefficients at Different Grid Points")
    print("=" * 80)
    print(f"{'ix':<4} {'r (km)':<12} {'ρ (m⁻³)':<12} {'T (K)':<8} {'n_e (m⁻³)':<12} {'Cm (m³/s)':<14} {'Ce_ij (m⁻³/s)':<14}")
    print("-" * 80)
    
    from_level = J110
    to_level = J101
    
    for ix in [0, 5, 10, 15, 19]:  # Sample points
        r_km = xc[ix] / 1e3
        rho = grid.props[IDEN, ix]
        T_neutral = grid.props[ITEMP, ix]
        n_e = grid.props[IELE, ix]
        T_e = grid.props[ITE, ix]
        
        # Get transition properties
        nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji = EinsteinCoeffs(
            from_level, to_level, r_h=2.5, species='ortho-H2O'
        )
        
        # Compute Cm (molecular collision coefficient)
        Cm = Einstein_Cm_ij(sigma_ij, rho, T_neutral)
        
        # Compute Ce_ij (electron collision rate)
        try:
            Ce = Einstein_Ce_ij(from_level, to_level, n_e, T_e)
        except ComartError as e:
            Ce = np.nan  # Skip if transition not defined
        
        print(f"{ix:<4} {r_km:<12.2e} {rho:<12.2e} {T_neutral:<8.1f} {n_e:<12.2e} {Cm:<14.3e} {Ce:<14.3e}")
    
    print("-" * 80)
    
    # Compute collision cooling/heating balance at each point
    print(f"\n" + "=" * 80)
    print("Collision Rate Analysis (Excitation vs De-excitation)")
    print("=" * 80)
    
    print(f"\nTransition: Level {from_level} <-> Level {to_level}")
    
    # Get transition frequency for energy calculation
    nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji = EinsteinCoeffs(
        from_level, to_level, r_h=2.5, species='ortho-H2O'
    )
    

    
    print(f"  Transition frequency: {nu_ij:.3e} Hz")
    print(f"  Collision cross-section (σ_ij): {sigma_ij:.3e} m²")
    
    print(f"\n{'ix':<4} {'r (km)':<12} {'Cm (neutral)':<14} {'Ce_ij (electron)':<16} {'Ratio (Ce/Cm)':<14}")
    print("-" * 80)
    
    for ix in [0, 5, 10, 15, 19]:
        r_km = xc[ix] / 1e3
        rho = grid.props[IDEN, ix]
        T_neutral = grid.props[ITEMP, ix]
        n_e = grid.props[IELE, ix]
        T_e = grid.props[ITE, ix]
        
        Cm = Einstein_Cm_ij(sigma_ij, rho, T_neutral)
        
        try:
            Ce = Einstein_Ce_ij(from_level, to_level, n_e, T_e)
            ratio = Ce / Cm if Cm > 0 else 0
        except ComartError:
            Ce = np.nan
            ratio = np.nan
        
        print(f"{ix:<4} {r_km:<12.2e} {Cm:<14.3e} {Ce:<16.3e} {ratio:<14.3e}")
    
    print("-" * 80)
    
    # Compute Cm, Ce_ij, and G_ij for all grid points and create plots
    print(f"\n" + "=" * 80)
    print("Computing Coefficients for Full Grid Profiles")
    print("=" * 80)
    
    # Arrays to store results
    Cm_array = np.zeros(nx)
    Ce_array = np.zeros(nx)

    
    # Get transition properties once (constant across grid)
    nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji = EinsteinCoeffs(
        from_level, to_level, r_h=2.5, species='ortho-H2O')
    
    # Compute coefficients at all grid points
    for ix in range(nx):
        rho = grid.props[IDEN, ix]
        T_neutral = grid.props[ITEMP, ix]
        n_e = grid.props[IELE, ix]
        T_e = grid.props[ITE, ix]
        

        # Cm: molecular collision coefficient
        Cm_array[ix] = Einstein_Cm_ij(sigma_ij, rho, T_neutral)
        
        # Ce_ij: electron collision rate
        Ce_array[ix] = Einstein_Ce_ij(from_level, to_level, n_e, T_e)


    print(f"\n✓ Computed coefficients for all {nx} grid points")

    
    # Create 3-panel plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=150)
    
    r_km = xc / 1e3  # Convert to km for plotting
    
    # Panel 1: Molecular collision coefficient (Cm)
    ax.axhline(A_ij, color='gray', linestyle='-.', label='Spontaneous emission')
    ax.loglog(r_km, Cm_array, c='dodgerblue', linestyle='-', label='H2O-H2O collision')
    ax.loglog(r_km, Ce_array, 'r--', label=r'H2O-e$^-$ collision')

    ax.axhline(G_ij, color='g', linestyle=':', label='Solar pumping')

    jv=1E-8/xc**2

    ax.loglog(r_km, jv*B_ij, 'm-x', label='Stimulated emission')
    ax.set_ylabel('Coeffs (/s)')
    ax.set_title(f"J110 -> J101 (r_h = 2.5 AU)")
    ax.grid(True, alpha=0.3, which='both')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('coma_Einstein_coeffs.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Plot saved to: coma_Einstein_coeffs.png")
    
    print("\n✓ Grid-based collision coefficient test completed successfully!")
    
except Exception as e:
    print(f"\n✗ Error in grid-based test: {e}")
    import traceback
    traceback.print_exc()
