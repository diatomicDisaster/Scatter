from time import process_time as time
import numpy as np

def pol_to_cart(pol_vec):
    """Convert vector from polar (r, theta phi) to Cartesian (x, y, z)"""
    x = pol_vec[0]*np.sin(pol_vec[1])*np.cos(pol_vec[2])
    y = pol_vec[0]*np.sin(pol_vec[1])*np.sin(pol_vec[2])
    z = pol_vec[0]*np.cos(pol_vec[1])
    return np.array([x, y, z])

def cart_to_pol(cart_vec):
    """Convert vector from Cartesian (x, y, z) to Polar (r, theta phi)"""
    r = np.sqrt(np.sum(cart_vec**2))
    theta = np.arccos(cart_vec[2]/r)
    phi = np.arctan(cart_vec[1]/cart_vec[0])
    if np.isnan(phi):
        phi = 0.
    elif cart_vec[0] < 0.:
        phi = phi + np.pi
    return np.array([r, theta, phi])

vec1 = np.array([-0.38434178,  0.16539601,  0.90825412])
vec1 = vec1/np.sqrt(np.sum(vec1**2))
x2 = np.array([1., 0., 0.])
y2 = np.array([0., 1., 0.])
z2 = np.array([0., 0., 1.])
vec2 = np.array([0, 0, 1.])
r, theta, phi = cart_to_pol(vec1)

Ry = np.array([[np.cos(theta), 0., np.sin(theta)],
                [0., 1., 0.],
                [-np.sin(theta), 0., np.cos(theta)]])
Rz = np.array([[np.cos(phi), -np.sin(phi), 0.],
                [np.sin(phi), np.cos(phi), 0.],
                [0., 0., 1.]])

new_vec = np.dot(Rz, np.dot(Ry, vec2))
new_x2 = np.dot(Rz, np.dot(Ry, x2))
new_y2 = np.dot(Rz, np.dot(Ry, y2))
new_z2 = np.dot(Rz, np.dot(Ry, z2))
# new_vec = np.dot(Ry, vec2)
# new_x2 = np.dot(Ry, x2)
# new_y2 = np.dot(Ry, y2)
# new_z2 = np.dot(Ry, z2)

print(vec2, new_vec)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')
for i, vec in enumerate([vec1, vec2, new_x2, new_y2, new_z2, new_vec]):
    ax.plot([0,vec[0]], [0, vec[1]], [0, vec[2]], label='{0}'.format(i))
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
plt.legend()
plt.savefig('traces2.png')
plt.show()
