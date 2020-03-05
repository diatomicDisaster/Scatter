## At each scattering order, we want to store a weight (w_abs) which corresponds 
## to the likelihood that a photon makes it this far without being absorbed, 
## and a weight (w_sca) which corresponds to the likelihood that the photon
## no longer undergoes any scattering or absorption between this point and 
## the detector. 

## At the end of photon path tracing we sum the product of the    
## weights for every photon to calulcate *intensity* :D 
## This is all for one wavelength!! we will add a spectral weight when we start importance sampling.

import numpy as np

## probability of no absorption along the path traversed in this scattering order.
w_abs = np.exp(-sum_i(length_in_cell_i*beta_sca_i))# sum over incoming cells: sum_i(length_in_cell_i*beta_sca_i) where beta_sca_i come from cells i traversed. 


## probability of no scatter+absorption along the path traversed from this 
## point up until the detector.
theta_, phi_ = self.grid.molecules[rand_mol].phase_func()
P = phase_func(theta_, phi_)
theta_d = 0 ## this is the angle of the incoming radiation beam to the detector. since we're assuming up-down/left-right propagation, for now set to 0.

w_sca = P*np.exp(-sum_i(length_in_cell_i*(beta_abs_i + beta_sca_i)))/np.cos(theta_d) 
# the sum of the cell lengths corresponds to the length of path straight along z axis to detector. 
# where i indexes the cells traveresed now TO the detector. 