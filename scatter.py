from numpy.random import uniform, normal
from numpy import multiply, array, pi, floor, inner, sin, cos, float64, append,\
    empty, zeros, fromfunction, flip, ones_like

#### WILF S'IL TE PLAIT, TA GUEULE :)
# Dev libraries
from numpy import indices
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import process_time

gridSize = array([100, 100, 100])
stopLen = 5.
maxSteps = 1000
nParts = 10000

class Cells:
    # Initialise cells with grid size (Nx, Ny, Nz)
    def __init__(self, gridSize):
        self.gS  = gridSize

    # Fill grid with integration constants C(z)
    def fill(self, fill_function, **kwargs):
        # If statements can be added below to select the form of C(z)
        if fill_function in ["Pressure", "pressure"]:
            self.fillFunc = self.pressure
        if fill_function in ["Constant", "constant"]:
            self.fillFunc = self.constant

        # Generate Numpy array from user defined function C(z)
        # !! Currently we have to flip here because grid is defined with z=0 at 
        # !!    top of atmosphere, and functional form is usually z=0 at bottom. 
        # !! We can change the definition of particles' reference frame so that 
        # !!    z=0 is bottom an particles propagate from top down.
        self.abs = flip(fromfunction(self.fillFunc, gridSize, **kwargs))

    # Functional forms of C(z) can be defined below
    # A function form C(z) based on atmospheric pressure with Earth defaults
    def pressure(self, i, j, k, 
                 p0=1.01325e5, T0=288.15, g=9.80665, L=6.5e-3, M=0.0289654):
        return p0*(1 - L*k/(self.gS[2]-1)/T0)**(g*M/(8.31447*L))
    # Constant C(z) = 1
    def constant(self, i, j, k):
        return ones_like(k)

# Cells initialised by instance of 'Cells' class
# The C(z) function is chosen to be the "pressure" case
grid = Cells(gridSize)
grid.fill("pressure", p0=1, T0=1, g=1, L=1, M=8.31447)
#grid.fill("constant")

# Plot C(z) function for checking
# fig = plt.figure(1)
# ax = fig.gca(projection='3d')
# gridInds = indices(gridSize)+0.5
# p = ax.scatter(gridInds[0], gridInds[1], gridInds[2], c=grid.abs.flatten(), marker='o')
# fig.colorbar(p)
# ax.set_xlim(0, gridSize[0])
# ax.set_ylim(0, gridSize[1])
# ax.set_zlim(0, gridSize[2])
# plt.savefig('C(z).png')
# plt.show()

# Prepare figure for particles plot
fig = plt.figure(2)
ax = fig.gca(projection='3d')

# Loop over incident photons
for particle in range(nParts):
    plotPoints = []
    # Random x, y position (z=0) and propagation angle (theta < pi/2)
    posVec  = multiply(array([*uniform(size=2), 0]), gridSize)
    randAng = multiply(uniform(size=2), array([pi/2, 2*pi]))
    propVec = array([sin(randAng[0])*cos(randAng[1]), 
                    sin(randAng[0])*sin(randAng[1]),
                    cos(randAng[0])])
    # Initialise path length, begin propagation
    sumPath = 0.
    inGrid = True
    nStep = 0
    while inGrid and nStep < maxSteps:
        nStep += 1
        plotPoints.append(posVec.copy())
        thisCell = floor(posVec) # Retrieve cell x, y, z indices
        topCell  = thisCell  +1 # Retrieve n + 1 cell
        # Find cell exit point
        cellPath = float('inf')
        for pointInd, wallPoint in enumerate([thisCell, topCell]):
            for normInd in range(3):
                normVec = zeros(3, dtype=float64)
                normVec[normInd] = 1.
                # Calculate distance to exit on nth cell wall 
                exitLen = inner((wallPoint - posVec), normVec) \
                    / inner(propVec, normVec)
                # Choose shortest exit length in forward direction
                if (exitLen > 0.) and (exitLen < cellPath):
                    cellPath = exitLen
                    exitWall = [pointInd, normInd]
        # If total path length exceeds scatter length: scatter
        optPath = cellPath
        if sumPath + cellPath > stopLen:
            cellPath -= sumPath + cellPath - stopLen
            randAng   = multiply(uniform(size=2), array([pi, 2*pi]))
            posVec   += cellPath * propVec
            propVec   = array([sin(randAng[0])*cos(randAng[1]), 
                               sin(randAng[0])*sin(randAng[1]),
                               cos(randAng[0])])
            sumPath = 0.
        else:
            sumPath += cellPath
            posVec  += cellPath*propVec
        # Stop if we're outside the grid
        if not all([0 <= posVec[i] <= gridSize[i] for i in range(3)]):
            inGrid = False

    # Plot each particles path to check
    points = array(plotPoints).T
    ax.plot(points[0], points[1], points[2])

# Stuff for plotting to check
ax.set_xlim(0, gridSize[0])
ax.set_ylim(0, gridSize[1])
ax.set_zlim(0, gridSize[2])
plt.savefig('traces.png')
plt.show()

