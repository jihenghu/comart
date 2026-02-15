"""
Comet atmosphere (coma) 2D property array module.

Uses NumPy arrays with C-contiguous memory layout for cache efficiency.
Properties are stored as the last dimension for memory locality.
"""

import numpy as np
from src.const import *
from src.error import ComartError


class ComaGrid:
    """
    1D property grid for comet coma simulation.
    
    Stores: number density, velocity, temperature, level 0 population ratio
    Dimension: nx (radial/distance points)
    
    Memory layout: Continuous C-contiguous (row-major) array
    Shape: (nprops, nx) for efficient Disort computation
    """
    
    def __init__(self, nx, nlevel=2, density=None, velocity=None, temperature=None, electron_density=None, electron_temperature=None):
        """
        Initialize the coma property grid.
        
        Args:
            nx: Number of grid points
            nlevel: Number of energy levels (for level population ratios)
            density: Optional initial density array
            velocity: Optional initial velocity array
            temperature: Optional initial temperature array
            electron_density: Optional initial electron density array
            electron_temperature: Optional initial electron temperature array
        """

        self.nx = nx
        self.nlevel = nlevel  # Number of energy levels (for level population ratios)
        self.nprops = NPROPS + nlevel  # Number of properties per grid point (density, velocity, temperature, electron density, electron temperature, level0, .. level n-1)
 
        NLEVEL = nlevel  # Set global constant for number of levels, used in other modules like Einstein_Aij
        
        # Create 1D property array: (nprops, nx)
        # Shape optimized for efficient extension to 3D and Disort computation
        self.props = np.zeros(
            (self.nprops, nx),
            dtype=np.float64,
            order='C'  # Explicitly use C-contiguous (row-major) order
        )
        
        # Grid coordinate arrays
        self.x = np.zeros(nx + 1, dtype=np.float64)   # Grid boundaries
        self.xc = np.zeros(nx, dtype=np.float64)      # Grid cell centers

        # Set initial properties if provided
        if density is not None:
            if len(density) != nx:
                raise ComartError(f"Density array size {len(density)} does not match nx={nx}")
            self.props[IDEN, :] = np.asarray(density, dtype=np.float64)
            self.init_level_populations()  # Initialize level populations based on density

        if velocity is not None:
            if len(velocity) != nx:
                raise ComartError(f"Velocity array size {len(velocity)} does not match nx={nx}")
            self.props[IVEL, :] = np.asarray(velocity, dtype=np.float64)

        if temperature is not None:
            if len(temperature) != nx:
                raise ComartError(f"Temperature array size {len(temperature)} does not match nx={nx}")
            self.props[ITEMP, :] = np.asarray(temperature, dtype=np.float64)
        
        if electron_density is not None:
            if len(electron_density) != nx:
                raise ComartError(f"Electron density array size {len(electron_density)} does not match nx={nx}")
            self.props[IELE, :] = np.asarray(electron_density, dtype=np.float64)

        if electron_temperature is not None:
            if len(electron_temperature) != nx:
                raise ComartError(f"Electron temperature array size {len(electron_temperature)} does not match nx={nx}")
            self.props[ITE, :] = np.asarray(electron_temperature, dtype=np.float64)



    def init_level_populations(self):
        """ Initialize level populations equally distributed """
    
        initial_pop = self.props[IDEN, :] / self.nlevel
        for i in range(self.nlevel-1):
            self.props[LEV0 + i, :] = initial_pop
        self.props[LEV0 + self.nlevel - 1, :] = self.props[IDEN, :] - np.sum(self.props[LEV0:LEV0 + self.nlevel - 1, :], axis=0)
        
        # Auto-check that level populations sum to 1
        self.check_all_levels()

    def make_level_transition(self, ix, from_level, to_level, ratio):
        """
        Perform a level transition at a given cell.
        
        Args:
            ix: Grid index
            from_level: Level index to transition from (0-based)
            to_level: Level index to transition to (0-based)
            ratio: Fraction of population to transition (0 to 1)
        
        Raises:
            ComartError: If levels are out of range or ratio is invalid
        """
        if from_level < 0 or from_level >= self.nlevel:
            raise ComartError(f"from_level {from_level} is out of range [0, {self.nlevel-1}]")
        if to_level < 0 or to_level >= self.nlevel:
            raise ComartError(f"to_level {to_level} is out of range [0, {self.nlevel-1}]")
        if ratio < 0 or ratio > 1:
            raise ComartError(f"ratio {ratio} must be in the range [0, 1]")
        
        from_idx = LEV0 + from_level
        to_idx = LEV0 + to_level
        
        # Calculate population transfer amount
        transfer_amount = self.props[from_idx, ix] * ratio
        
        # Update level populations
        self.props[from_idx, ix] -= transfer_amount
        self.props[to_idx, ix] += transfer_amount
        
        # Auto-check that level populations sum to 1 ??
        # self.check_level_sum(ix) 
        # # Optional: check after each transition, or check all levels at the end of a series of transitions to reduce overhead

    def get_density(self):
        """Get the number density array."""
        return self.props[IDEN, :]
    
    def get_velocity(self):
        """Get the velocity array."""
        return self.props[IVEL, :]
    
    def get_temperature(self):
        """Get the temperature array."""
        return self.props[ITEMP, :]
    
    def get_electron_density(self):
        """Get the electron density array."""
        return self.props[IELE, :]
    
    def get_electron_temperature(self):
        """Get the electron temperature array."""
        return self.props[ITE, :]

    def get_level_populations(self):
        """Get the level populations"""
        return self.props[LEV0:, :]
    
    def set_property(self, prop_idx, ix, value):
        """
        Set a property value at a grid point.
        
        Args:
            prop_idx: Property index (use IDEN, IVEL, etc.)
            ix: Grid index
            value: Property value
        """

        self.props[prop_idx, ix] = value

        # Auto-check level populations sum to 1 if setting a level property
        if prop_idx >= LEV0:
            self.check_level_sum(ix)


    def get_property(self, prop_idx, ix):
        """
        Get a property value at a grid point.
        
        Args:
            prop_idx: Property index (use IDEN, IVEL, etc.)
            ix: Grid index
        
        Returns:
            Property value
        """
        return self.props[prop_idx, ix]
    
    def get_cell(self, ix):
        """
        Get all properties at a grid point.
        
        Args:
            ix: Grid index
        
        Returns:
            Array of [density, velocity, temperature, level0_ratio, ...]
        """
        return self.props[:, ix]
    
    def check_level_sum(self, ix, tolerance=0.1):
        """
        Check that level populations sum to 1 at a given cell.
        
        Args:
            ix: Grid index
            tolerance: Tolerance for deviation from 1.0 (default 1e-10)
        
        Raises:
            ComartError: If level populations don't sum to 1 within tolerance
        """
        level_sum = np.sum(self.props[LEV0:, ix])
        if not np.isclose(level_sum, self.props[IDEN, ix], atol=tolerance):
            raise ComartError(
                f"Level populations at cell {ix} sum to {level_sum:.6e}, "
                f"expected {self.props[IDEN, ix]:.6e} (tolerance={tolerance})"
            )
    

    def check_all_levels(self, tolerance=0.1):
        """
        Check that level populations sum to 1 at all grid points.
        
        Args:
            tolerance: Tolerance for deviation from expected density (default 0.1)
        
        Raises:
            ComartError: If any cell's level populations don't sum to expected density
        """
        level_sums = np.sum(self.props[LEV0:, :], axis=0)
        bad_cells = np.where(~np.isclose(level_sums, self.props[IDEN, :], atol=tolerance))[0]
        
        if len(bad_cells) > 0:
            raise ComartError(
                f"Level populations don't sum to {self.props[IDEN, 0]:.6e} at {len(bad_cells)} cells: {bad_cells[:]}... "
                f"(values: {level_sums[bad_cells[:]]})"
            )
        else:
            print(f"All level populations sum to expected density within tolerance {tolerance}.")
    
    @property
    def is_contiguous(self):
        """Check if array is C-contiguous in memory."""
        return self.props.flags['C_CONTIGUOUS']
    
    @property
    def shape(self):
        """Return array shape (nprops, nx)."""
        return self.props.shape
    
    @property
    def nbytes(self):
        """Return memory usage in bytes."""
        return self.props.nbytes
    
