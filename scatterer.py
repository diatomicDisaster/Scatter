import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
start_time = time.time()


# Possibly needs to become an inherited class to
# allow for different phase functions. Need to think about this.
class Molecule:
    def __init__(self, name="no name", phase_SD=0.5):
        """A class to define molecules with different width phase functions
        
        Inputs
        name : str
            A name for the molecule, default 'no name'.
        phase_SD : float
            Standard deviation (radians) of the molecule's phase function.
        """
        self.name = name
        self.phase_SD = phase_SD

    def phase_func(self, theta):
        """Returns probability of scattering at angle theta."""
        return 2*np.exp(-theta/2*self.phase_SD**2) # x2 because theta > 0

    def inv_phase_func(self, rand_state=np.random):
        """Draws random theta, phi from Gaussian scattering phase function.
        Returns
            theta : float
                A random azimuthal angle theta > 0 from a normal 
                distribution centred on theta = 0.
            phi : float
                A random polar angle phi from a uniform
                distribution in the range [0, 2*pi].
        """
        theta = np.abs(rand_state.normal(loc=0.0, scale=self.phase_SD))
        phi = rand_state.uniform()*2*np.pi
        return theta, phi

    def z_abund_func(self, z):
        """Simple exponential abundance function.
        Input
            z : float
                Height above the surface, z=0.
        Returns
            A(z) : float
                The relative abundance at height z. """
        return 5e-2 #np.exp(-z) 

class Grid:
    def __init__(self, grid_size=(100,100,100), molecules=['arb']):
        """A grid instance can be defined with a given grid size and
        populated with a given set of molecules.

        Inputs
            grid_size : array-like shape=(3)
                The number of grid cells along the (x, y, z) axes.
            molecules : list
                A list of molecules from the molculeDict to populate
                the grid with.
        """
        self.grid_size = grid_size
        self.abs_grid  = np.zeros(self.grid_size)
        self.molecules = []
        for mol in molecules:
            self.molecules.append(moleculeDict.get(mol))
        for i in range(self.grid_size[2]):
            z = (grid_size[2]-i-1)/self.grid_size[2]
            for molecule in self.molecules:
                self.abs_grid[:, :, i] += molecule.z_abund_func(z)
        self.abs_grid = self.abs_grid/len(self.molecules)
        self.sca_grid = self.abs_grid
        self.sum_grid = self.abs_grid + self.sca_grid

class Scatter:
    def __init__(self, grid):
        """A scattering instance on a given grid.
        Inputs:
            grid : Grid object
                An instance of a grid object on which to perform the
                scattering.
        """
        self.grid = grid

    def propagate(self, max_step=999, rand_seed=None, nu=None):
        """Propagate a photon through the grid.
        
        Inputs
            max_step : int
                Maximum number of path steps to make before giving up.
            rand_seed : int
                Seed to generate reproducible random numbers.
                Used for initial position and propagation vector.
            nu : float
                Wavenumber of photon to propagate.
                
        """
        self._max_step   = max_step
        self._rand       = np.random.RandomState(rand_seed)
        self.weights     = [[],[]]
        
        # Random initial x, y position (z=0) and propagation vector
        rand = self._rand.uniform(size=(2,2))
        self.pos_vec     = np.full((max_step+1, 3), np.nan)
        self.pos_vec[0]  = np.append(rand[0], 0.) * self.grid.grid_size
        self.prp_vec     = self._pol_to_cart(
            np.insert(rand[1]*np.array([.25, 2.])*np.pi,  0, 1.))
        
        ## Wilf: you had -np.log(1 - rho/beta_scat) but I'm not sure this makes
        ##       sense. That function is only defined for beta_scat>rho, so generates
        ##       infinite path lengths in many cases.
        # Initialise random optical depth integrated path length for first step
        self.max_abs     = -np.log(1 - self._rand.uniform()) 
        self.sum_abs     = 0.
        
        self.left_grid = False
        for n_step in range(max_step):
            self._cell = np.floor(self.pos_vec[n_step]).astype(int) # x,y,z indices of cell
            if (any(0 > self._cell)) or any(self._cell >= self.grid.grid_size):
                self.left_grid = True
                self.pos_vec[n_step] = np.array([np.nan, np.nan, np.nan])
                break
            stop_len = float('inf')
            grad = self.prp_vec < 0
            for p in range(2): # iterate over cell point (0,0,0) and (1,1,1)
                p0 = self._cell + p*(1 - 2*p*grad)
                for n in range(3): # iterate over x,y,z normal vectors 
                    norm    = np.array([0.,0.,0.], dtype=np.float64)
                    norm[n] = 1.
                    test_len = np.inner((p0 - self.pos_vec[n_step]), norm) \
                                / np.inner(self.prp_vec, norm)
                    # shortest positive intercept distance is exit point
                    if 0. < test_len < stop_len:
                        stop_len = test_len
            beta_abs   = self.grid.abs_grid[tuple(self._cell)]
            self.sum_abs += stop_len*beta_abs
            # check whether max_abs is reached along this path
            if self.sum_abs >= self.max_abs:
                stop_len = (self.sum_abs - self.max_abs)/beta_abs
                self.sum_abs = self.max_abs
            self.pos_vec[n_step + 1] = self.pos_vec[n_step] + stop_len*self.prp_vec
            if self.sum_abs is self.max_abs:
                theta = self._scatter()
                self._detect(theta, self.pos_vec[n_step+1])
                # New random optical depth and zero integrated path length
                self.max_abs = -np.log(1 - self._rand.uniform())
                self.sum_abs = 0.
        ## Wilf: Currently we get values of I_p > 1, is this right?
        ##       Maybe there needs to be some form of normalisation somewhere,
        ##       e.g dividing by the number of scattering orders or total path length?
        # Sum of weight products over scatter orders
        if len(self.weights[0]) == 0 and len(self.weights[1]) ==0:
            self.intens_p = "no weights" # means photon left grid before scattered
        else:
            self.intens_p = np.dot(self.weights[0], self.weights[1])
    
    def _pol_to_cart(self, pol_vec):
        """Convert vector from polar (r, theta phi) to Cartesian (x, y, z)"""
        x = pol_vec[0]*np.sin(pol_vec[1])*np.cos(pol_vec[2])
        y = pol_vec[0]*np.sin(pol_vec[1])*np.sin(pol_vec[2])
        z = pol_vec[0]*np.cos(pol_vec[1])
        return np.array([x, y, z])

    def _cart_to_pol(self, cart_vec):
        """Convert vector from Cartesian (x, y, z) to Polar (r, theta phi)"""
        r = np.sqrt(np.sum(cart_vec**2))
        theta = np.arccos(cart_vec[2]/r)
        phi = np.arctan(cart_vec[1]/cart_vec[0])
        if np.isnan(phi):
            phi = 0.
        elif cart_vec[0] < 0:
            phi = phi + np.pi
        return np.array([r, theta, phi])

    def _detect(self, theta, pos_vec):
        """Calculate weights for scattering and absorption paths at current scatter order.
        Inputs
            theta : float
                Phase function angle 0 <= theta <= pi.
            pos_vec : float, array-like
                Position vector of scattering event.
        """
        # Sum of absorption and scattering along path to detector
        detect_path  = ((self._cell[2] + 1) - pos_vec[2]) \
            * self.grid.sum_grid[tuple(self._cell)] \
            + np.sum(self.grid.sum_grid[self._cell[0], 
                                        self._cell[1], 
                                        self._cell[2]+1:])
        ## Wilf: which molecule's phase function do we use below?
        P = self.grid.molecules[0].phase_func(theta) # phase function probability
        self.weights[1].append(np.exp(-detect_path)*P) # scattering weight
        self.weights[0].append(np.exp(-self.sum_abs)) # absorption weight

    def _scatter(self):
        """Scatter the photon from a molecule following some phase function."""
        rand_mol  = self._rand.randint(0, high=len(self.grid.molecules))
        theta_, phi_ = self.grid.molecules[rand_mol].inv_phase_func(rand_state=self._rand)
        new_prp_ = self._pol_to_cart(np.array([1., theta_, phi_]))
        r, theta, phi = self._cart_to_pol(self.prp_vec)
        Ry = np.array([[np.cos(theta), 0., np.sin(theta)],
                       [0., 1., 0.],
                       [-np.sin(theta), 0., np.cos(theta)]])
        Rz = np.array([[np.cos(phi), -np.sin(phi), 0.],
                       [np.sin(phi), np.cos(phi), 0.],
                       [0., 0., 1.]])
        new_prp = np.dot(Rz, np.dot(Ry, new_prp_))
        self.prp_vec = new_prp
        return theta_

# A dictionary of molecules
moleculeDict = {
    "arb": Molecule(name = "arbitrary") # arbitrary molecules
    }
# Define a 100 x 100 x 100 grid with default molecules
gridA = Grid(grid_size=(100,100,100))

## For plotting the grid of absorption values
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# gridInds = np.indices(gridA.grid_size)+.5
# p = ax.scatter(gridInds[0], gridInds[1], gridInds[2], c=gridA.abs_grid.flatten(), marker='o')
# fig.colorbar(p)
# ax.set_xlim(0, gridA.grid_size[0])
# ax.set_ylim(0, gridA.grid_size[1])
# ax.set_zlim(0, gridA.grid_size[2])
# plt.savefig('C(z).png')
# plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

# Invoke 100 instances of 'Scatter' on the grid and plot.
scatters = []
nl = 'No scattering'
for i in range(50):
    scatters.append(Scatter(gridA))
    scatters[i].propagate()
    p = scatters[i].intens_p
    points = scatters[i].pos_vec.T
    if p == "no weights":
        ax.plot(points[0], points[1], points[2], c='red', label=nl)
        nl = '_nolegend_'
    else:
        ax.plot(points[0], points[1], points[2], c='blue')
ax.set_xlim(0, gridA.grid_size[0])
ax.set_ylim(0, gridA.grid_size[1])
ax.set_zlim(0, gridA.grid_size[2])
ax.set_zlabel('Atmospheric Depth (arb. units)')
ax.invert_zaxis()
plt.legend()
plt.savefig('traces2.png')
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()

