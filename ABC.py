import numpy as np

class ABC:
    """
    class for approximate bayesian computation
    """

    def __init__(self):
        self.thetas = []
        self.dists = [] #distances between distributions of simulation and given data

    def abc(self, simulation, distFunc, eps, sampleSize=100, numParam=4, prior=np.random.uniform):
        '''
        performs approximate bayesian computation
        :param simulation: function that takes theta/-s and outputs a distribution
        :param distFunc: computes distance between distributions of simulation and given data
        :param eps: threshold for acceptance region for distance between distributions
        :param sampleSize: number of runs of simulation
        :param numParam: number of parameters for simulation
        :param prior: prior distribution for sampling theta
        :return:
        '''
        for i in range(sampleSize):
            param = prior(size=numParam)
            simData = simulation(param)
            dist = distFunc(simData)
            #if dist <= eps:
            self.thetas.append(param)
            self.dists.append(dist)
            print(dist, param)
