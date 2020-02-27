import numpy as np

class Scatter:
    def __init__(self, gridSize=[100,100,100]):
        self.gridSize = gridSize

    def propagate(self, nSteps=100, angRange=[.5, 2], randSeed=2):
        """Propagate a photon through the grid"""
        
        self.randSeed_  = randSeed
        self.angRange_ = np.array(angleRange*np.pi)
        self.nSteps_   = maxSteps
        self.tau_      = -np.log(1 - np.random.uniform(self.randSeed_))
        
        # Generate random initial x, y position and propagation angle
        random = np.random.RandomState(self.randSeed_)
        randPos, randAng = random.uniform(size=(2,2))
        
        self.posVec_  = np.multiply(np.append(randPos, 0.), self.gridSize)
        self.propVec_ = np.multiply(randAng, self.angRange_)

        for nStep in range(self.maxSteps):
            self.cell_ = np.floor(self.posVec_)
            exitLen = None
            # loop (2 p0 values) * (3 normal vectors) = 6 cell walls
            for p in range(2):
                p0 = self.cell_ + p
                for n in range(3):
                    norm    = np.zeros(3, dtype=float64)
                    norm[n] = 1.
                    testLen = np.inner((p0 - self.posVec_), norm) \
                                / np.inner(self.propVec_, norm)
                    # Update if exit length is closer
                    exitLen = (testLen if 0. < testLen < exitLen)
        optDepth = exitLen*self.absGrid[*self.cell_]
            
    def rand_vec(self, scale):
            



