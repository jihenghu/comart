"""
Example: Monte Carlo solver for level population equilibrium on ComaGrid.

Demonstrates solution of statistical equilibrium using MC method for
a two-level system on a cometary coma grid.
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
from src.MonteCarlo.mc_solver import MCLevelSolver


print("=" * 80)
print("Monte Carlo Level Population Solver on ComaGrid")
print("=" * 80)

try:
    # Create coma grid
    nx = 20
    r_min = 1e3      # 1 km
    r_max = 1e8      # 100,000 km
    
    xc = np.logspace(np.log10(r_min), np.log10(r_max), nx)
    xb = np.zeros(nx + 1)
    xb[0] = r_min
    xb[-1] = r_max
    for i in range(1, nx):
        xb[i] = np.sqrt(xc[i-1] * xc[i])
    
    print(f"\nCreating ComaGrid with {nx} layers")
    print(f"  Distance range: {r_min/1e3:.1e} - {r_max/1e3:.1e} km")
    
    # Compute profiles
    print("Computing neutral and electron profiles...")
    v0 = 720.0
    c_scale = 3.5e3
    velocity = velocity_tangent(xc, v0=v0, c=c_scale)
    
    T0 = 10.0
    temperature = temperature_revesedist(xc, a=7, b=1, T0=T0, r0=2000)
    
    Q = 1e26
    density = calc_density(Q, xc, velocity, r_h=2.5, model_type='isotropic', photo_dissociation=True)
    
    # Electron profiles
    nelec, telec = electron_Biver1997(Q, r_h=2.5, r=xc, v=velocity, Tkin=temperature, Temax=10000)
    
    # Create grid
    grid = ComaGrid(nx, nlevel=2, density=density, velocity=velocity, 
                    temperature=temperature, electron_density=nelec, 
                    electron_temperature=telec)
    grid.x = xb
    grid.xc = xc
    
    print(f"✓ ComaGrid created")
    
    # Get Einstein coefficients and collision rates
    level1 = J101
    level2 = J110
    
    nu_ij, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, sigma_ij, sigma_ji = EinsteinCoeffs(
        level1, level2, r_h=2.5, species='ortho-H2O'
    )
    
    print(f"\nTransition: Level {level1} <-> Level {level2}")
    print(f"  Frequency: {nu_ij:.3e} Hz")
    print(f"  A_ij (upward): {A_ij:.3e} s^-1")
    print(f"  A_ji (downward): {A_ji:.3e} s^-1")
    print(f"  G_ij (radiative excitation): {G_ij:.3e} s^-1")
    print(f"  G_ji (radiative de-excitation): {G_ji:.3e} s^-1")
    
    # Test MC solver at selected grid points
    print(f"\n" + "=" * 80)
    print("Monte Carlo Solution at Different Grid Points")
    print("=" * 80)
    
    # Storage for results
    mc_results = {}
    
    # Test at a few grid points
    test_indices = [1, 10, 15, 18]  # Near, intermediate, far
    
    for ix in test_indices:
        print(f"\n{'=' * 80}")
        print(f"Grid point {ix}: r = {xc[ix]/1e3:.2e} km")
        print(f"{'=' * 80}")
        
        r_km = xc[ix] / 1e3
        rho = grid.props[IDEN, ix]
        T_neutral = grid.props[ITEMP, ix]
        n_e = grid.props[IELE, ix]
        T_e = grid.props[ITE, ix]
        
        # Compute collision coefficients
        Cm_ij = Einstein_Cm_ij(sigma_ij, rho, T_neutral)
        Cm_ji = Einstein_Cm_ij(sigma_ji, rho, T_neutral)
    
        Ce_ij = Einstein_Ce_ij(level1, level2, n_e, T_e)
        Ce_ji = Einstein_Ce_ij(level2, level1, n_e, T_e)

        
        print(f"\nEnvironmental conditions:")
        print(f"  Neutral density: {rho:.2e} m^-3")
        print(f"  Temperature: {T_neutral:.1f} K")
        print(f"  Electron density: {n_e:.2e} m^-3")
        print(f"  Electron temperature: {T_e:.1f} K")
        print(f"\nCollision coefficients:")
        print(f"  Cm_ij (collision coeff): {Cm_ij:.3e} m^3/s")
        print(f"  Cm_ji (collision coeff): {Cm_ji:.3e} m^3/s")
        print(f"  Ce_ij (electron collision): {Ce_ij:.3e} m^-3/s")
        print(f"  Ce_ji (electron collision): {Ce_ji:.3e} m^-3/s")

        # print Aij, Aji, Gij, Gji for reference
        print(f"\nRadiative coefficients:")
        print(f"  A_ij (downward): {A_ij:.3e} s^-1")
        print(f"  A_ji (upward): {A_ji:.3e} s^-1")
        print(f"  B_ij (absorption): {B_ij:.3e} m^3/J/s^2")
        print(f"  B_ji (stimulated emission): {B_ji:.3e} m^3/J/s^2")
        print(f"  G_ij (radiative excitation): {G_ij:.3e} s^-1")
        print(f"  G_ji (radiative de-excitation): {G_ji:.3e} s^-1")
        
        # Initial populations (unequal distribution)
        initial_pop = np.array([rho *0.2, rho *0.8])
        
        # Create and run MC solver
        solver = MCLevelSolver(nlevel=2, max_steps=100000, tol=1e-8, verbose=False)
        jv_rad = 1e-15  # Radiation field intensity
        pop_final, converged, nsteps = solver.solve(
            initial_populations=initial_pop,
            jv=jv_rad,
            A_ij=A_ij,
            A_ji=A_ji,
            B_ij=B_ij,
            B_ji=B_ji,
            G_ij=G_ij,
            G_ji=G_ji,
            Cm_ij=Cm_ij,
            Cm_ji=Cm_ji,
            Ce_ij=Ce_ij,
            Ce_ji=Ce_ji,
            time_scale=1.0
        )
        
        # Store results
        mc_results[ix] = {
            'populations': pop_final,
            'converged': converged,
            'nsteps': nsteps,
            'populations_history': solver.populations_history,
            'history': solver.get_convergence_history(),
            'ratio1': pop_final[0] / np.sum(pop_final),
            'ratio2': pop_final[1] / np.sum(pop_final)
        }
        
        # Print results
        print(f"\nMC Solution (converged after {nsteps} steps):")
        print(f"  Initial populations: n0={initial_pop[0]:.2e}, n1={initial_pop[1]:.2e}")
        print(f"  Final populations:   n0={pop_final[0]:.2e}, n1={pop_final[1]:.2e}")
        print(f"  Equilibrium ratio n1/n0: {pop_final[1]/pop_final[0]:.3e}")
        print(f"  Convergence: {'✓' if converged else '✗'}")
    
    # Create convergence plots for selected points
    print(f"\n" + "=" * 80)
    print("Creating convergence plots...")
    print("=" * 80)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
    axes = axes.flatten()
    for idx, ix in enumerate(test_indices):
            ax = axes[idx]
            history = mc_results[ix]['populations_history']
            r_km = xc[ix] / 1e3
            
            indx= history[:, 0] + history[:, 1] > 0
            ax.plot(history[indx, 0], 'b-', linewidth=1.5, label='Level 0', alpha=0.7)
            ax.plot(history[indx, 1], 'r-', linewidth=1.5, label='Level 1', alpha=0.7)

            ax.set_xlabel('MC Step', fontsize=10, fontweight='bold')
            ax.set_ylabel('Population (m⁻³)', fontsize=10, fontweight='bold')
            ax.set_title(f'r = {r_km:.2e} km (converged in {mc_results[ix]["nsteps"]} steps)', 
                        fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9)


    
    plt.tight_layout()
    plt.savefig('mc_convergence_curves.png', dpi=150, bbox_inches='tight')
    print(f"✓ Plot saved to: mc_convergence_curves.png")
    
    # Summary
    print(f"\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nMC solver results at {len(test_indices)} grid points:")
    print(f"{'Point':<6} {'r (km)':<12} {'n0 (m⁻³)':<14} {'n1 (m⁻³)':<14} {'n1/N':<12} {'n2/N':<12} {'Steps':<8}")
    print("-" * 80)
    
    for ix in test_indices:
        pop = mc_results[ix]['populations']
        ratio1 = mc_results[ix]['ratio1']
        ratio2 = mc_results[ix]['ratio2']
        nsteps = mc_results[ix]['nsteps']
        r_km = xc[ix] / 1e3
        print(f"{ix:<6} {r_km:<12.2e} {pop[0]:<14.2e} {pop[1]:<14.2e} {ratio1:<12.3f} {ratio2:<12.3f} {nsteps:<8}")
    
    print("-" * 80)
    print("\n✓ Monte Carlo solver test completed successfully!")
    
except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
