import numpy as np
        
Nx = 10
Ny = 10
Nz = 10
size = np.array([Nx, Ny, Nz])

Lx = 47.893
Ly = Lx
Lz = 1
nMols = 2

def returnNum(z, nTot):
    A1 = z*z
    A2 = 1-z*z
    return np.array([A1, A2])*nTot

cells = np.zeros((Nx, Ny, Nz, nMols), dtype=np.float64)
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz):
            z = k * Lz/Nz
            cells[i][j][k] = returnNum(z, z*1e6)
            
initAng = np.random.uniform(size=2)*np.pi
initPos = np.random.uniform(size=2)
initCell = np.append(np.floor(np.multiply(initPos, size[:2])), 0)

tau = 0.5

norms = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

dead = False
while not dead:
    rNow = rNext
    thisCell = np.floor(rNow) # x, y, z indices of current cell
    nextCell = rNow + 1 # next (x+1, y+1, z+1) cell
    unitVec = rNow/np.linalg.norm(rNow) # unit vector for current direction
    # need two p0 points to specify 6 planes - (0, 0, 0) and (1, 1, 1)
    for pInd, pZero in enumerate([thisCell, nextCell]):
        # and three unit vectors to specify 6 planes
        for normInd, norm in enumerate(np.array([[1, 0, 0], 
                                                [0, 1, 0], 
                                                [0, 0, 1]])):
            # calculate path length until plane intercept for all 6 planes
            t = np.inner((rNow - pZero), norm))/np.inner(unitVec, norm)
            # find the first plane we pass through - i.e when we exit the cell
            if t < path:
                path = t
                # index the cell wall by it's normal vector and p0
                wall = [normInd, pInd]

print(initPos, np.multiply(initPos, [Lx, Ly]), initCell)
