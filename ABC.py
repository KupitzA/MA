import numpy as np
from matplotlib import pyplot as plt

class ABC:
    """
    class for approximate bayesian computation
    """

    def __init__(self):
        self.thetas = []
        self.dists = [] #distances between distributions of simulation and given data
        self.simData = [] #accepted simulated distributions

    def abc(self, simulation, distFunc, eps, sampleSize=10000, numParam=4, prior=np.random.uniform):
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
            if dist <= eps:
                self.thetas.append(param)
                self.dists.append(dist)
                self.simData.append(simData)
                print(dist, param)
        #compute mean theta and distribution if data accepted
        if len(self.simData) != 0:
            theta = []
            for i in range(numParam):
                theta.append(np.mean([t[i] for t in self.thetas]))
            print(theta)
            #compute mean distribution
            accumulated = self.meanDistri()
            lists = sorted(accumulated.items()) # sorted by key, return a list of tuples
            x, y = zip(*lists) # unpack a list of pairs into two tuples
            plt.plot(x, y)
            plt.xlabel('pattern value')
            plt.ylabel('distribution value')
            plt.title('pattern distribution of simulation')
            plt.show()

    def meanDistri(self):
        accumulated = dict()
        for d in self.simData:
            for k, v in d.items():
                accumulated[k] = accumulated.get(k, 0) + v
        accumulated = {x: float(y/(len(self.simData)*10000)) for x, y in accumulated.items()}
        return accumulated
