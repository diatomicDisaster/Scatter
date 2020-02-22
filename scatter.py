from numpy.random import uniform, normal
from numpy import multiply, array, pi, floor, inner, sin, cos, float64, append, empty, zeros

# Dev libs
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import process_time
fig = plt.figure()
ax = fig.gca(projection='3d')

gridSize = array([100, 100, 100])
stopLen = 50.
maxSteps = 1000

# Loop over incident photons
for particle in range(10):
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
        topCell  = thisCell +1 # Retrieve n + 1 cell
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
        if not all([0 <= thisCell[i] <= gridSize[i] for i in range(3)]):
            inGrid = False

    # Stuff for plotting to check
    points = array(plotPoints).T
    ax.plot(points[0], points[1], points[2])

# Stuff for plotting to check
ax.set_xlim(0, gridSize[0])
ax.set_ylim(0, gridSize[1])
ax.set_zlim(0, gridSize[2])
plt.savefig('traces.png')
plt.show()

