import numpy as np

class Scatter:
    def __init__(self, grid_size=(100,100,100)):
        self.grid_size = np.array(grid_size)

    def propagate(self, max_step=100, rand_seed=2):
        """Propagate a photon through the grid
        
        Inputs
            max_step : int
                Maximum number of path steps to make before giving up.
            rand_seed : int
                Seed to generate reproducible random numbers.
                Used for initial position and tau value.
                
        """
        self._max_step = max_step
        self._rand     = np.random.RandomState(rand_seed)
        self.pos_vec   = np.append(rand[0], 0.) \
            * self.grid_size # x,y position (z=0)
        self.prp_vec   = _pol_to_cart( # entry angle (theta < pi/2)
                np.insert(rand[1]*np.array([.5, 2.])*np.pi, 0, 1.)) 
        self.max_tau   = -np.log(1 - rand.uniform()) # integrated path length
        self.sum_tau   = 0.
        
        for nStep in range(self._max_steps):
            self._cell = np.floor(self.pos_vec) # x,y,z indices of cell
            if (0 > any(self._cell)) and any(self._cell > self.grid_size)
            self._beta_ext   = self.abs_grid[tuple(self._cell)]
            stop_len = None
            for p in range(2): # iterate over cell point (0,0,0) and (1,1,1)
                p0 = self._cell + p
                for n in range(3): # iterate over x,y,z normal vectors 
                    norm    = np.array([0.,0.,0.], dtype=float64)
                    norm[n] = 1.
                    test_len = np.inner((p0 - self.pos_vec), norm) \
                                / np.inner(self.prp_vec, norm)
                    # shortest positive intercept distance is exit point
                    if 0. < test_len < exit_len:
                        stop_len = test_len
            self.sum_tau += stop_len*beta_ext
            # check whether max_tau is reached along this path
            if self.sum_tau >= self.max_tau:
                stop_len = (self.max_tau - self.sum_tau)/beta_ext
                self.sum_tau = self.max_tau
            self.pos_vec += stop_len*self.prp_vec
            if self.sum_tau is self.max_tau:
                self._scatter()
    
    def _pol_to_cart(self, pol_vec):
        """Convert vector from polar to Cartesian coordinates"""
        x = pol_vec[0]*np.sin(pol_vec[1])*np.cos(pol_vec[2])
        y = pol_vec[0]*np.sin(pol_vec[1])*np.sin(pol_vec[2])
        z = pol_vec[0]*np.cos(pol_vec[1])
        return np.array([x, y, z])
        



