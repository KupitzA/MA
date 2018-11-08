import numpy as np
from matplotlib import pyplot as plt
from simulation import Simulation

class ABC:
    """
    class for approximate bayesian computation
    """

    def __init__(self, distributionData, computePatternDistribution, L):
        '''
        :param distributionData: pattern distribution file to compare with
        :param computePatternDistribution: function that takes theta/-s and outputs a distribution
        :param L: number of CpGs
        '''
        self.distributionData = distributionData
        self.computePatternDistribution = computePatternDistribution
        self.L = L
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
        k = 20 #number of elements to be accepted
        for i in range(sampleSize):
            #draw k-times from prior
            if len(self.thetas) <= k:
                param = prior(size=numParam)
            else:
                param = self.ownPrior(numParam, self.thetas, k)
            distributionSim = self.computePatternDistribution(param)
            dist = distFunc(self.distributionData, distributionSim)
            #store values in acceptance region
            if dist < eps:
                self.thetas.append(param)
                self.distances.append(dist)
                self.distribution.append(distributionSim)
                print(dist, param)
                #keep only k parameters with lowest distance
                if len(self.thetas) > k:
                    todel = self.distances.index(max(self.distances))
                    self.thetas.pop(todel)
                    self.distances.pop(todel)
                    self.distribution.pop(todel)
                    eps = np.mean(self.distances) #resize epsilon
        #compute mean theta and distribution if data accepted
        if len(self.distribution) != 0:
            theta = []
            for i in range(numParam):
                theta.append(np.mean([t[i] for t in self.thetas[-k:]]))
            print(theta)
            #compute mean distribution
            accumulated = self.meanDistri(3)
            lists = sorted(accumulated.items()) # sorted by key, return a list of tuples
            print(lists)
            x, y = zip(*lists) # unpack a list of pairs into two tuples
            plt.plot(x, y)
            plt.xlabel('pattern value')
            plt.ylabel('distribution value')
            plt.title('pattern distribution of simulation')
            plt.show()

    def meanDistri(self, kbest):
        """
        compute mean value for all parameters
        :return: mean value for all parameters
        """
        accumulated = dict()
        for d in self.distribution[-kbest:]:
            for k, v in d.items():
                accumulated[k] = accumulated.get(k, 0) + v
        accumulated = {x: float(y/100) for x, y in accumulated.items()}
        return accumulated

    def ownPrior(self, numParam, thetas, k):
        """
        computes a posterior distribution based on k prior draws by building a normal distribution around one of them
        :param numParam: number of parameters to draw
        :param thetas: list of accepted values
        :param k: number of previously drawn values on which posterior is based
        :return: list of drawn parameters
        """
        prior = []
        bestk = thetas[-k:]
        for i in range(numParam):
            Sum = sum(self.distances[-k:])
            invSum = sum([Sum-self.distances[-l] for l in range(k)])
            prob = 0.0 #probbility for each accepted combination of parameters to form new distribution to draw from
            q = np.random.random()
            #choose from k values with probability prob
            for j in range(k):
                prob += (Sum - self.distances[-j])/invSum
                if q <= prob:
                    break
            m = bestk[j][i]
            sd = np.mean([abs(bestk[j][i] - bestk[l][i]) for l in range(k)]) #choose mean deviation from this to other
            # values as standard deviation for posterior distribution
            #m = np.mean(weightedSum)
            #sd = np.std(weightedSum)
            p = np.random.normal(m, sd) #built new posterior
            #keep p in range 0<=p<=1
            if p < 0:
                p = 0
            elif p > 1:
                p = 1
            prior.append(p)
            print(m,sd)
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
        for k in range(4**self.L):
            keyData = k if k in distributionData else 0
            keySim = k if k in distributionSim else 0
            dist += abs(keyData - keySim)
        return dist

    def dist(self, distributionData, ditributionSim):
        """
        distance function for ABC, distance of to distributions is sum over weighted distances of all relative distances
        between the pairwise patterns from each distribution
        :param patterns: pattern distribution of simulation
        :return: distance between pattern distribution of simulation and pattern distribution of KO-file
        """
        dist = 0
        for keyData in range(4**self.L):
            valueData = distributionData[keyData] if keyData in distributionData else 0
            for keySim in range(4**self.L):
                valueSim = ditributionSim[keySim] if keySim in ditributionSim else 0
                dist += self.w(keyData, keySim) * (valueData - valueSim)**2
        return dist

    def w(self, keyData, keySim):
        """
        computes a weight relative to two patterns
        :param keyData: pattern from file
        :param keySim: pattern from simulation
        :return: weight w
        """
        w = 0.0
        while keyData != 0 or keySim != 0:
            modData = keyData % 4
            modSim = keySim % 4
            if modData+modSim == 3: #then the methylation state at both strands is complementary
                w += 2.0
            elif modData != modSim: #one position of methylation pattern different
                w += 1.0
            keyData = (keyData-modData) / 4
            keySim = (keySim-modSim) / 4
        if self.L != 0:
            w /= self.L
        return w


#DNT1KO:
sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT1KO.txt", [13, 14], True)
distriData = sim.computePatternDistribution([0,          1,          0,              1])

#DNMT3KO:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT3abKO.txt", [13, 14], False, True)
#distriData = sim.computePatternDistribution([0,          1,          1,              0])

#WT:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatWTJ1C.txt", [13, 14])
#distriData = sim.computePatternDistribution([[0.1, 0.8, 0.8, 0], [0.5, 0.5, 0, 1]])
abc = ABC(distriData, sim.computePatternDistribution, sim.L)

abc.abc(abc.dist, 20.0)
lists = sorted(distriData.items()) # sorted by key, return a list of tuples
print(lists)
x, y = zip(*lists) # unpack a list of pairs into two tuples
plt.plot(x, y)
plt.xlabel('pattern value')
plt.ylabel('distribution value')
plt.title('pattern distribution of data')
plt.show()
#print(sim.dist(sim.computePatternDistribution([0.5,          0.5,          0,              1])))
