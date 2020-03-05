import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
start_time = time.time()


# Possibly needs to become an inherited class to
# allow for different phase functions. Need to think about this.
class Molecule:
    def __init__(self, name="no name", phase_SD=np.pi/4):
        """A class to define molecules with different width phase functions
        
        Inputs
        name : str
            A name for the molecule, default 'no name'.
        phase_SD : float
            Standard deviation (radians) of the molecule's phase function.
        """
        self.name = name
        self.phase_SD = phase_SD

    def phase_func(self):
        """Gaussian scattering phase function.
        Returns
            theta : float
                A random azimuthal angle theta > 0 from a normal 
                distribution centred on theta = 0.
            phi : float
                A random polar angle phi from a uniform
                distribution in the range [0, 2*pi].
        """
        theta = np.abs(np.random.normal(loc=0.0, scale=self.phase_SD))
        phi = np.random.uniform()*2*np.pi
        return theta, phi

    def z_abund_func(self, z):
        """Simple exponential abundance function.
        Input
            z : float
                Height above the surface, z=0.
        Returns
            A(z) : float
                The relative abundance at height z. """
        return np.exp(-z)

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

class Scatter:
    def __init__(self, grid):
        """A scattering instance on a given grid.
        Inputs:
            grid : Grid object
                An instance of a grid object on which to perform the
                scattering.
        """
        self.grid = grid

    def propagate(self, max_step=10000, rand_seed=None): #changed max step from 999 to 10000
        """Propagate a photon through the grid.
        
        Inputs
            max_step : int
                Maximum number of path steps to make before giving up.
            rand_seed : int
                Seed to generate reproducible random numbers.
                Used for initial position and tau value.
                
        """
        self._max_step   = max_step
        self._rand       = np.random.RandomState(rand_seed)
        
        rand = self._rand.uniform(size=(2,2))
        self.pos_vec     = np.full((max_step+1, 3), np.nan)
        self.pos_vec[0]  = np.append(rand[0], 0.) * self.grid.grid_size
        self.prp_vec     = self._pol_to_cart(np.insert(rand[1]*np.array([.25, 2.])*np.pi,  0, 1.))
        self.max_tau     = -np.log(1 - self._rand.uniform()) # integrated path length. LARA: need to scale second term in log by 1/beta_sca, which is a function of the cell!! :)
        self.sum_tau     = 0.
        
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
            beta_sca   = self.grid.abs_grid[tuple(self._cell)]
            self.sum_tau += stop_len*beta_sca
            # check whether max_tau is reached along this path
            if self.sum_tau >= self.max_tau:
                stop_len = (self.sum_tau - self.max_tau)/beta_sca
                self.sum_tau = self.max_tau
            self.pos_vec[n_step + 1] = self.pos_vec[n_step] + stop_len*self.prp_vec
            if self.sum_tau is self.max_tau:
                self._scatter()
    
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
        
    def _scatter(self):
        """Scatter the photon from a molecule following some phase function"""
        rand_mol  = np.random.randint(0, high=len(self.grid.molecules))
        theta_, phi_ = self.grid.molecules[rand_mol].phase_func()
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
        self.sum_tau = 0.
        return self

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

# Perform 100 scatter instances on the grid and plot. Lara: you mean 100 photons not scattering orders??
scatters = []
for i in range(1000):
    scatters.append(Scatter(gridA))
    scatters[i].propagate()
    points = scatters[i].pos_vec.T
    ax.plot(points[0], points[1], points[2])
ax.set_xlim(0, gridA.grid_size[0])
ax.set_ylim(0, gridA.grid_size[1])
ax.set_zlim(0, gridA.grid_size[2])
plt.savefig('traces2.png')
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()

