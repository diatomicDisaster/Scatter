import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

    def propagate(self, max_step=999, rand_seed=None):
        ## Lara: could we possibly add the beta_sca parameter in the arguments here? It will be the sum of all the beta_sca_j for 
        # molecule types j in the cells that are propagated through. For now we can take this to be homogeneous just to get it running! 
        # So let's say beta_sca_j = 0.01 for 10 types of molecule --> beta_sca = 0.1
        beta_sca = 0.1
        """Propagate a photon through the grid.
        
        Inputs
            max_step : int
                Maximum number of path steps to make before giving up. 
            rand_seed : int
                Seed to generate reproducible random numbers.
                Used for initial position and tau value.
                
        """
        ##Lara: I wonder if we should let the paths run and only stop them when they reach the top of the z axis?? as opposed to max step.

        self._max_step   = max_step
        self._rand       = np.random.RandomState(rand_seed)
        
        rand = self._rand.uniform(size=(2,2))
        self.pos_vec     = np.full((max_step+1, 3), np.nan)
        self.pos_vec[0]  = np.append(rand[0], 0.) * self.grid.grid_size
        self.prp_vec     = self._pol_to_cart(np.insert(rand[1]*np.array([.25, 2.])*np.pi,  0, 1.))
        self.max_tau     = -np.log(1 - (self._rand.uniform())/beta_sca) # draw a random integrated optical path length [0, infty]
        self.sum_tau     = 0. # summation of traversed optical path length

        self.left_grid = False
        for n_step in range(max_step):
            self._cell = np.floor(self.pos_vec[n_step]).astype(int) # x,y,z indices of cell
            # Check we are still in the grid
            if (any(0 > self._cell)) or any(self._cell >= self.grid.grid_size):
                self.left_grid = True
                self.pos_vec[n_step] = np.array([np.nan, np.nan, np.nan])
                break
            stop_len = float('inf')
            grad = self.prp_vec < 0
            # Iterate over cell walls to find exit point
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
            #beta_sca   = self.grid.abs_grid[tuple(self._cell)]  # integrate beta along path length
            
            self.sum_tau += stop_len*beta_sca # add to summation
            if self.sum_tau >= self.max_tau: # check if max_tau is reached along this path
                stop_len = (self.sum_tau - self.max_tau)/beta_sca # shorten stop length
                self.sum_tau = self.max_tau
            self.pos_vec[n_step + 1] = self.pos_vec[n_step] + stop_len*self.prp_vec
            if self.sum_tau is self.max_tau: # scatter if scattering optical depth reached
                self._scatter()


            ## Lara: store absorption weight! :)
            # do i need to initialise both of these as empty arrays?? 
            self.w_abs_this_step = np.exp(-stop_len*beta_sca)
            # needs:
            # length_in_cell = individual cell lengths traveresed x the beta_sca param in each cell
            # BUT FOR NOW we can just use stop_len multiplied by beta_sca in cells traversed so far?? need to check that stop_len is free length travelled and not starting from zero! 
            # Then we cumulate previous path weights to obtain an aborption weight for each scattering order:
            self.w_abs += self.w_abs_this_step # lol is this even right ??

    
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
        """Scatter the photon from a molecule, along a new propagation vector.
        Scattering angle defined by the a random molecule's phase function."""
        rand_mol  = np.random.randint(0, high=len(self.grid.molecules)) # select random molecules
        theta_, phi_ = self.grid.molecules[rand_mol].phase_func() # calculate scattering angle
        new_prp_ = self._pol_to_cart(np.array([1., theta_, phi_])) # new propagation vector k'
        r, theta, phi = self._cart_to_pol(self.prp_vec)
        Ry = np.array([[np.cos(theta), 0., np.sin(theta)],
                       [0., 1., 0.],
                       [-np.sin(theta), 0., np.cos(theta)]])
        Rz = np.array([[np.cos(phi), -np.sin(phi), 0.],
                       [np.sin(phi), np.cos(phi), 0.],
                       [0., 0., 1.]])
        new_prp = np.dot(Rz, np.dot(Ry, new_prp_)) # transform k' from k=z to global Cartesian frame
        self.prp_vec = new_prp
        self.sum_tau = 0.

        #Lara: store scattering weight :) 
        ## This is the probability of no scatter+absorption along the path traversed from this 
        ## point up until the detector.
        P = phase_func(theta_, phi_)
        theta_d = 0 ## this is the angle of the incoming radiation beam to the detector. since we're assuming up-down/left-right propagation, for now set to 0.
        z_length =      # the cartesian length of the path straight to the detector _before_ we scatter. total distance to det - distance weve travelled in total so far!! 
        beta_abs_to_det = f(z_length)    # the absorption coefficients summed over the cells traversed now _to_ the detector. lets for now put chosen number. z_length would determine no. of cells. 
        beta_sca_to_det = g(z_length)   # the scattering coefficients summed over the cells traversed now _to_ the detector. lets for now put chosen number. z_length would determine no. of cells. 

        self.w_sca = P*np.exp(-z_length*(beta_abs_to_det + beta_sca_to_det))/np.cos(theta_d) 

        # Now want to, at each scattering order, store the value of w_abs and w_sca, and multiply together to obtain intensity contribution at this step.
        # When we reach the end of the grid, so the location of the detector, we will sum all of these values to get intensity contribution for all possible paths. 
        # Then divide by the number of photons - argument of "for i in range", call it N_p - to obtain averaged intensity for this wavelength. YAY

         I_order = self.w_sca*self.w_abs # this is for ONE photon at EACH scattering order. 
         I_photon = sum_total_orders{I_order}     # Maybe these three lines come out of iteration? quite sure this isnt implemented correctly lol 
         I = (1/N_p)*sum_photons{I_photon}      # Not sure if this is making sense ??

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

# Perform 100 scatter instances on the grid and plot
scatters = []
for i in range(100):
    scatters.append(Scatter(gridA))
    scatters[i].propagate()
    points = scatters[i].pos_vec.T
    ax.plot(points[0], points[1], points[2])
ax.set_xlim(0, gridA.grid_size[0])
ax.set_ylim(0, gridA.grid_size[1])
ax.set_zlim(0, gridA.grid_size[2])
plt.savefig('traces2.png')
plt.show()
