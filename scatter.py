import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')

Nx = 5
Ny = 5
Nz = 5
size = np.array([Nx, Ny, Nz])

Lx = 47.893
Ly = Lx
Lz = 1
nMols = 2

# Dummy routine for assigning abundance to cells
def returnNum(z, nTot):
    A1 = z*z
    return A1*nTot
cells = np.zeros((Nx, Ny, Nz), dtype=np.float64)
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            z = k * Lz/Nz
            cells[i][j][k] = returnNum(z, z*1e6)

# Random initial position and angle
initAng = np.random.uniform(size=2)*np.pi
initPos = np.random.uniform(size=2)
initCell = np.append(np.floor(np.multiply(initPos, size[:2])), 0)

# Path length
tau = 5

# Garbage that needs cleaning
dead = False
sumPath = 0.
r0Vec = np.array([2.5, 3.5, 0.])
dVec  = np.array([1., 1., 1.])
dVec  = dVec / np.linalg.norm(dVec)
points = []

while not dead:
    path = float('inf')
    thisCell = np.floor(r0Vec) # x, y, z indices of current cell
    nextCell = thisCell + 1 # next (x+1, y+1, z+1) cell
    #unitVec = rNow/np.linalg.norm(rNow) # unit vector for current direction
    # need two p0 points to specify 6 planes - (0, 0, 0) and (1, 1, 1)
    for p0Ind, p0Vec in enumerate([thisCell, nextCell]):
        # and three unit vectors to specify 6 planes
        for normInd, norm in enumerate(np.array([[1, 0, 0], 
                                                [0, 1, 0], 
                                                [0, 0, 1]])):
            # calculate path length until plane intercept for all 6 planes
            t = np.inner((p0Vec - r0Vec), norm)/np.inner(dVec, norm)
            # find the first plane we pass through - i.e when we exit the cell
            if ((abs(t) < path) and (t > 0.0)):
                sumPath += t
                path = t
                # index the cell wall by it's normal vector and p0
                wall = [normInd, p0Ind]
    # save points for plotting
    points.append(r0Vec)
    # update position
    r0Vec = r0Vec + path*dVec
    # check path length
    if sumPath > 5:
      dead = True

# Plotting to check results
xs = np.array([point[0] for point in points])
ys = np.array([point[1] for point in points])
zs = np.array([point[2] for point in points])
ax.scatter(xs, ys, zs)
X, Y = np.meshgrid(np.arange(0, Nx+1, 1), np.arange(0, Ny+1, 1))
for k in range(Nz+1):
    Z = np.full((Nx+1, Ny+1), k)
    ax.plot_wireframe(X, Y, Z, lw=0.1)

ax.set_xlim(0, Nx)
ax.set_ylim(0, Ny)
ax.set_zlim(0, Nz)
plt.show()
