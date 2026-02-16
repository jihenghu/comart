"""
Monte Carlo solver for statistical equilibrium in multi-level system.

Solves for equilibrium level populations using Monte Carlo sampling of
transition probabilities. Handles two-level system with transitions driven by:

  - Radiative decay (spontaneous emission via A_ij)
  - Absorption and stimulated emission (via B_ij and B_ji)
  - Radiative pumping (solar-pumped excitation via G_ij)
  - Molecular collisions (neutral-molecule collisions via Cm_ij)
  - Electron collisions (electron-molecule collisions via Ce_ij)

"""

import numpy as np
from src.error import ComartError


class MCLevelSolver:
    """
    Monte Carlo solver for level population equilibrium.
    
    Advances a two-level (or multi-level) system through MC steps until
    convergence to equilibrium populations.
    """
    
    def __init__(self, nlevel=2, max_steps=10000, tol=1e-4, verbose=False):
        """
        Initialize MC solver for level populations.
        
        Parameters
        ----------
        nlevel : int
            Number of energy levels (default 2 for two-level system)
        max_steps : int
            Maximum number of MC steps to run (default 10000)
        tol : float
            Convergence tolerance for population change (default 1e-6, i.e. 0.0001%)
        verbose : bool
            Print progress information (default False)
        """
        self.nlevel = nlevel
        self.max_steps = max_steps
        self.tol = tol
        self.verbose = verbose
        
        # History tracking
        self.populations_history = np.zeros((self.max_steps+1, self.nlevel), dtype=np.float64)  # Shape: (nsteps, nlevel)
        self.transition_counts = np.zeros((self.nlevel, self.nlevel), dtype=np.int32)    # Count transitions: (i, j) -> count


    def solve(self, initial_populations, jv, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, 
              Cm_ij, Cm_ji , Ce_ij, Ce_ji, time_scale=1.0):
        """
        Solve for equilibrium level populations using Monte Carlo.
        
        Parameters
        ----------
        initial_populations : array_like
            Initial level populations (length nlevel), sum should equal total density
        jv : float
            radiation field intensity at transition frequency (W/m^2/Hz/sr)
        A_ij : float
            Einstein A coefficient for transition i->j (downward, s^-1)
        A_ji : float
            Einstein A coefficient for transition j->i (upward, s^-1)
        B_ij : float
            Einstein B coefficient for transition i->j (absorption, m^3/J/s^2)
        B_ji : float
            Einstein B coefficient for transition j->i (stimulated emission, m^3/J/s^2)
        G_ij : float
            Radiative pumping rate (excitation via solar radiation, s^-1)
        G_ji : float
            Radiative pumping rate (de-excitation, s^-1)
        Cm_ij : float
            Molecular collision coefficient (m^3/s)
        Cm_ji : float
            Molecular collision coefficient (m^3/s)
        Ce_ij : float
            Electron collision rate coefficient (m^-3/s)
        Ce_ji : float
            Electron collision rate coefficient (m^-3/s)
        time_scale : float
            Time scale in seconds for one MC step (default 1.0)
        
        Returns
        -------
        populations_final : ndarray
            Final equilibrium level populations (length nlevel)
        converged : bool
            Whether solution converged within max_steps
        nsteps : int
            Number of MC steps taken
        """
        
        # Validate inputs
        if len(initial_populations) != self.nlevel:
            raise ComartError(f"initial_populations length must match nlevel {self.nlevel}")
        
        total_density = np.sum(initial_populations)
        if total_density <= 0:
            raise ComartError("Total density must be positive")
        
        # Initialize populations
        populations = np.asarray(initial_populations).copy()
        self.populations_history[0] = initial_populations.copy()
        self.transition_counts = np.zeros((self.nlevel, self.nlevel), dtype=np.int32)
        
        # Transition rate matrix (for two-level system)
        # rates[i][j] = rate coefficient for i->j transition
        rates = self._build_rate_matrix(jv, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, Cm_ij, Cm_ji, Ce_ij, Ce_ji)
        
        if self.verbose:
            print(f"\nMC Level Population Solver")
            print(f"  Levels: {self.nlevel}")
            print(f"  Max steps: {self.max_steps}")
            print(f"  Convergence tolerance: {self.tol}")
            print(f"  Initial populations: {initial_populations}")
            print(f"  Total density: {total_density:.2e}")
            print(f"\nRate coefficients:")
            print(f"  A_ij (downward emission): {A_ij:.3e} s^-1")
            print(f"  A_ji (upward emission): {A_ji:.3e} s^-1")
            print(f"  G_ij (radiative excitation): {G_ij:.3e} s^-1")
            print(f"  G_ji (radiative de-excitation): {G_ji:.3e} s^-1")
            print(f"  Cm_ij (collision coefficient): {Cm_ij:.3e} m^3/s")
            print(f"  Ce_ij (electron collision): {Ce_ij:.3e} m^-3/s")
        
        # Run MC steps
        converged = False
        for step in range(self.max_steps):
            # Save old populations for convergence check
            populations_old = populations.copy()
            
            # Perform MC step: sample transitions
            populations = self._mc_step(populations, rates, time_scale, total_density)
            
            # Store history
            self.populations_history[step+1] = populations.copy()
            
            # Check convergence: if the maximum change in all level populations is less than tol, we consider it converged
            delta_pop = (populations[0] - np.mean(self.populations_history[max(0, step+1-10):step+1,0])) / total_density
            
            # print(f"Step {step:5d}: populations = {populations}, delta = {delta_pop:.3e}")

            if self.verbose and (step % 100 == 0 or step < 10):
                print(f"Step {step:5d}: populations = {populations}, delta = {delta_pop:.3e}")
            
            # Convergence criterion
            if delta_pop < self.tol:
                converged = True
                if self.verbose:
                    print(f"\n✓ Converged at step {step} with delta = {delta_pop:.3e}")
                break
        
        if not converged and self.verbose:
            print(f"\n⚠ Did not converge after {self.max_steps} steps (delta = {delta_pop:.3e})")
        
        return populations, converged, step + 1


    def _build_rate_matrix(self, jv, A_ij, A_ji, B_ij, B_ji, G_ij, G_ji, Cm_ij, Cm_ji, Ce_ij, Ce_ji):
        """
        Build transition rate matrix.
        
        For a two-level system (0, 1):
        - rates[0] = rate for 0->1 transition (excitation)
        - rates[1] = rate for 1->0 transition (de-excitation)
        
        Returns
        -------
        rates : ndarray
            use np.array to store rates as values
        """

        rates = np.zeros((2), dtype=np.float64)
        
        if self.nlevel == 2:
            # Excitation rates (0 -> 1)   
            rates[0] = B_ij * jv + Cm_ij + G_ij + Ce_ij  # B_ij * J_nu + G_ij + collision terms
            
            # De-excitation rates (1 -> 0)
            rates[1] = A_ji + B_ji * jv + G_ji + Cm_ji + Ce_ji  # Cm_ij is bidirectional in detail balance
        else:
            raise ComartError("Only 2-level system implemented")
        
        return rates


    def _mc_step(self, populations, rates, time_scale, total_density):
        """
        Perform one Monte Carlo step: sample transitions from current population.
        
        Parameters
        ----------
        populations : ndarray
            Current level populations (length nlevel)
        rates : dict
            Transition rate coefficients
        time_scale : float
            Time interval for one step (seconds)
        total_density : float
            Total particle density (conserved)
        
        Returns
        -------
        populations_new : ndarray
            Updated level populations after one MC step
        """
        populations_new = populations.copy()
        
        # Two-level system (0, 1)
        pop0, pop1 = populations_new[0], populations_new[1]
        
        rate_01 = rates[0]  # 0->1 excitation
        rate_10 = rates[1]  # 1->0 de-excitation
        
        # Total transition rates
        trans_01 = rate_01 * pop0 * time_scale  # Number of 0->1 transitions
        trans_10 = rate_10 * pop1 * time_scale  # Number of 1->0 transitions
        
        # For small probabilities, use exact; for large, saturate at available population，This prevents negative populations
        n_01 = min(int(np.random.poisson(trans_01)), int(pop0))
        n_10 = min(int(np.random.poisson(trans_10)), int(pop1))
        
        # Update populations
        populations_new[0] = pop0 - n_01 + n_10
        populations_new[1] = pop1 + n_01 - n_10
        
        # Ensure total density is conserved
        populations_new = populations_new / np.sum(populations_new) * total_density
        
        # Track transitions (array-based indexing)
        self.transition_counts[0, 1] += n_01
        self.transition_counts[1, 0] += n_10
        
        return populations_new


    def get_convergence_history(self):
        """
        Get population evolution history.
        
        Returns
        -------
        history : ndarray
            Population history shape (nsteps, nlevel)
        """
        return np.array(self.populations_history)


    def print_summary(self, populations_final, converged, nsteps):
        """
        Print summary of MC solution.
        
        Parameters
        ----------
        populations_final : ndarray
            Final equilibrium populations
        converged : bool
            Convergence status
        nsteps : int
            Number of steps taken
        """
        print("\n" + "=" * 70)
        print("Monte Carlo Solution Summary")
        print("=" * 70)
        
        print(f"\nConvergence: {'✓ Yes' if converged else '✗ No'}")
        print(f"Steps taken: {nsteps}")
        print(f"Initial populations: {self.populations_history[0]}")
        print(f"Final populations: {populations_final}")
        print(f"Total density (conserved): {np.sum(populations_final):.3e}")
        
        # Population ratios
        if populations_final[0] > 0:
            ratio = populations_final[1] / populations_final[0]
            print(f"\nPopulation ratio n1/n0: {ratio:.3e}")
        
        # Transition statistics (array-based)
        print(f"\nTransition counts:")
        for i in range(self.nlevel):
            for j in range(self.nlevel):
                if i != j and self.transition_counts[i, j] > 0:
                    print(f"  {i}->{j}: {int(self.transition_counts[i, j])} transitions")
