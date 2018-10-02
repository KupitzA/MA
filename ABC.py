import numpy as np
from matplotlib import pyplot as plt
from simulation import Simulation

class ABC:
    """
    class for approximate bayesian computation
    """

    def __init__(self, distributionData, computePatternDistribution):
        '''
        :param distributionData: pattern distribution file to compare with
        :param computePatternDistribution: function that takes theta/-s and outputs a distribution
        '''
        self.distributionData = distributionData
        self.computePatternDistribution = computePatternDistribution
        self.thetas = []
        self.distances = [] #distances between distributions of simulation and given data
        self.distribution = [] #accepted simulated distributions

    def abc(self, distFunc, eps, sampleSize=10000, numParam=4, prior=np.random.uniform):
        '''
        performs approximate bayesian computation
        :param distFunc: computes distance between distributions of simulation and given data
        :param eps: threshold for acceptance region for distance between distributions
        :param sampleSize: number of runs of simulation
        :param numParam: number of parameters for simulation
        :param prior: prior distribution for sampling theta
        :return:
        '''
        improvements = 1000 #number of improvements, where drawing of thetas is improved by mean of accepted thetas
        for i in range(improvements):
            for j in range(int(sampleSize/improvements)):
                if i == 0:
                    param = prior(size=numParam)
                else:
                    param = self.ownPrior(numParam, self.thetas)
                distributionSim = self.computePatternDistribution(param)
                dist = distFunc(self.distributionData, distributionSim)
                if dist < eps:
                    self.thetas.append(param)
                    self.distances.append(dist)
                    self.distribution.append(distributionSim)
                    print(dist, param)
            eps /= 2.0
        #compute mean theta and distribution if data accepted
        if len(self.distribution) != 0:
            theta = []
            for i in range(numParam):
                theta.append(np.mean([t[i] for t in self.thetas]))
            print(theta)
            #compute mean distribution
            #accumulated = self.meanDistri()
            #lists = sorted(accumulated.items()) # sorted by key, return a list of tuples
            #x, y = zip(*lists) # unpack a list of pairs into two tuples
            #plt.plot(x, y)
            #plt.xlabel('pattern value')
            #plt.ylabel('distribution value')
            #plt.title('pattern distribution of simulation')
            #plt.show()

    def meanDistri(self):
        accumulated = dict()
        for d in self.distribution:
            for k, v in d.items():
                accumulated[k] = accumulated.get(k, 0) + v
        accumulated = {x: float(y/(len(self.distribution)*10000)) for x, y in accumulated.items()}
        return accumulated

    def ownPrior(self, numParam, thetas):
        prior = []
        for i in range(numParam):
            if len(thetas) >= 10:
                tail = 10
            else:
                tail = len(thetas)-1
            m = np.mean(thetas[-tail:][i])
            sd = np.std(thetas[-tail:][i])
            prior.append(np.random.normal(m, sd))
            print(m)
        return prior

    def distanceFunction(self, distributionData, distributionSim):
        '''
        A distance function for an ABC-algorithm, computes the distance of two distributions as sum of differences
        between all counts of patterns
        :param distributionSim: pattern distribution of simulation
        :param distributionData: pattern distribution of data
        :return: distance between pattern distribution of simulation and pattern distribution of given data
        '''
        dist = 0
        for k, v in enumerate(distributionData):
            dist += abs(distributionSim[k]-v) if k in distributionSim else v
        return dist

    def dist(self, distributionData, ditributionSim):
        """
        distance function for ABC, distance of to distributions is sum over weighted distances of all relative distances
        between the pairwise patterns from each distribution
        :param patterns: pattern distribution of simulation
        :return: distance between pattern distribution of simulation and pattern distribution of KO-file
        """
        dist = 0
        for k, v in enumerate(distributionData):
            for key, value in ditributionSim.items():
                dist += self.w(k, key) * (v - value)**2
        return dist

    def w(self, keyData, keySim):
        """
        computes a weight relative to two patterns
        :param keyData: pattern from file
        :param keySim: pattern from simulation
        :return: weight w
        """
        w = 0.0
        L = 0.0
        while keyData != 0 or keySim != 0:
            L += 1.0
            modData = keyData % 4
            modSim = keySim % 4
            if modData+modSim == 3: #then the methylation state at both strands is complementary
                w += 2.0
            elif modData == modSim: #one position of methylation pattern different
                w += 1.0
            keyData = (keyData-modData) / 4
            keySim = (keySim-modSim) / 4
        if L != 0:
            w /= L
        return w


#DNT1KO:
sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT1KO.txt", [13, 14], True)
distriData = sim.computePatternDistribution([0.1, 0.8, 0.8, 0])

#DNMT3KO:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT3abKO.txt", [13, 14], False, True)

#WT:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatWTJ1C.txt", [13, 14])
#distriData = sim.computePatternDistribution([[0.1, 0.8, 0.8, 0], [0.5, 0.5, 0, 1]])
abc = ABC(distriData, sim.computePatternDistribution)

abc.abc(abc.dist, 100000.0)
#plt.plot(sim.distributionKO)
#plt.xlabel('pattern value')
#plt.ylabel('distribution value')
#plt.title('pattern distribution of data')
#plt.show()
#print(sim.dist(sim.computePatternDistribution([0.5,          0.5,          0,              1])))
