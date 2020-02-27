import numpy as np

class Scatter:
    def __init__(self, 
                 gridSize=[100,100,100]
                 ):
        self.gridSize = gridSize

    def propagate(self, 
                  randSeed=2, 
                  angRange=[np.pi/2, 2*np.pi]
                  ):
        """Propagate a photon through the grid"""
        
        self.randSeed_   = randSeed
        self.angRange_ = np.array(angleRange)
        
        # Generate random initial x, y position and angle
        random = np.random.RandomState(self.randSeed)
        randPos, randAng = random.uniform(size=(2,2))
        randPos = np.append(rand, 0.)
        self.initPos_ = np.multiply(randPos, gridSize)
        self.initAng_ = np.multiply(randAng, self.angRange_)

        # Initialise position and propagation vectors
        posVec  = self.initPos_
        propVec = self.initAng_
        for nStep in range(maxSteps):
            # Retrieve indices of current cell and +1 cell
            p0 = np.floor(pos)
            p1 = p0 + 1
            cellPath 

