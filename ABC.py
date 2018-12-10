import statistics
import numpy as np
import math
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
        self.weights = self.buildWeights()

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
            if len(self.thetas) < k:
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
            sds = []
            for i in range(numParam):
                parami = [t[i] for t in self.thetas[-k:]]
                theta.append(np.mean(parami))
                sds.append(statistics.stdev(parami))
            print(theta, sds)
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

    def meanDistri(self, kbest=1):
        """
        compute mean value for all parameters
        :return: mean value for all parameters
        """
        accumulated = dict()
        distances2 = self.distances
        for d in range(kbest):
            minimum = distances2.index(min(distances2))
            d = self.distribution[minimum]
            for k, v in d.items():
                accumulated[k] = accumulated.get(k, 0) + v
            distances2.pop(minimum)
        accumulated = {x: float(y/kbest) for x, y in accumulated.items()}
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
            keyData = distributionData[k] if k in distributionData else 0
            keySim = distributionSim[k] if k in distributionSim else 0
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
        #for keyData, valueData in distributionData.items():
            for keySim in range(4**self.L):
                valueSim = ditributionSim[keySim] if keySim in ditributionSim else 0
            #for keySim, valueSim in ditributionSim.items():
                dist += self.w(keyData, keySim) * (valueData - valueSim+0.00000001)**2
        return dist

    def w(self, keyData, keySim):
        """
        computes a weight relative to two patterns
        :param keyData: pattern from file
        :param keySim: pattern from simulation
        :return: weight w
        """
        w = 0.0
        #iterations = self.L
        #while keyData != 0 or keySim != 0:
        for i in range(self.L):
            #iterations -= 1
            modData = keyData % 4
            modSim = keySim % 4
            if modData+modSim == 3: #then the methylation state at both strands is complementary
                w += 3.0
            elif modData != modSim: #one position of methylation pattern different
                w += 2.0
            else:
                w += 1.0
            keyData = (keyData-modData) / 4
            keySim = (keySim-modSim) / 4
        #w += iterations
        return w

    def balance(self, pattern):
        mid = self.L/2
        r = 0
        l = 0
        for i in range(self.L):
            modP = pattern % 4
            for j in range(0, modP+2, 2):
                if i <= mid-1:
                    l += 1
                elif i >= mid:
                    r += 1
            pattern = pattern // 4
        b = (r-l)/abs(r-l) if r-l != 0 else r-l
        return b


    def mahalonisDist(self, distributionData, distributionSim):
        X = [distributionData[k] if k in distributionData else 0 for k in range(4**self.L)]
        Y = [distributionSim[k] if k in distributionSim else 0 for k in range(4**self.L)]
        #mhy = np.average(X, weights=range(0,4**self.L))
        #ny = np.average(Y, weights=range(0,4**self.L))
        #zipped = zip([x-mhy for x in X], [y-ny for y in Y])
        zipped = zip(X, Y)
        #cov = np.cov(list(zipped))
        #inv = self.invert(cov)
        cov = self.weights
        subtr = np.subtract(X, Y)
        d = np.inner(np.dot(subtr, cov), subtr)
        return math.sqrt(d)

    def invert(self, cov):
        try:
            inv = np.linalg.inv(cov)
        except np.linalg.LinAlgError:
            #cov[0][0] += 0.001
            #self.invert(cov)
            pass
        else:
            pass
        return cov

    def buildWeights(self):
        #w = [[1, 2, 2, 4], [2, 1, 3, 2], [2, 3, 1, 2], [4, 2, 2, 1]]
        cov = np.zeros((4**self.L, 4**self.L))
        for i in range(0, 4**self.L):
            balanceI = self.balance(i)
            for j in range(i, 4**self.L):
                balanceJ = self.balance(j)
                balance = abs((balanceI-balanceJ)/2)
                w = (self.w(i, j) + balance) / self.L
                cov[i][j] = w
                cov[j][i] = w
        return cov


#DNT1KO:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT1KO.txt", [13, 14], True)
#distriData = sim.computePatternDistribution([0.5, 0.5, 0, 1])

#DNMT3KO:
sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT3abKO.txt", [13, 14], False, True)
distriData = sim.computePatternDistribution([0.1, 0.8, 0.8, 0])

#WT:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatWTJ1C.txt", [13, 14])
#distriData = sim.computePatternDistribution([[0.1, 0.8, 0.8, 0], [0.5, 0.5, 0, 1]])

#distriData = {i:sim.distributionKO[i] for i in range(len(sim.distributionKO))}
abc = ABC(distriData, sim.computePatternDistribution, sim.L)

abc.abc(abc.mahalonisDist, 0.1)
lists = sorted(distriData.items()) # sorted by key, return a list of tuples
print(lists)
x, y = zip(*lists) # unpack a list of pairs into two tuples
plt.plot(x, y)
plt.xlabel('pattern value')
plt.ylabel('distribution value')
plt.title('pattern distribution of data')
plt.show()
#lists = sorted(distriSim.items()) # sorted by key, return a list of tuples
#print(lists)
#x, y = zip(*lists) # unpack a list of pairs into two tuples
#plt.plot(x, y)
#plt.xlabel('pattern value')
#plt.ylabel('distribution value')
#plt.title('pattern distribution of sim')
#plt.show()
#print(abc.distanceFunction(distriData, distriSim))
