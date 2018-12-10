import math
import numpy as np

def distanceFunction(distributionData, distributionSim, L=1):
        '''
        A distance function for an ABC-algorithm, computes the distance of two distributions as sum of differences
        between all counts of patterns
        :param distributionSim: pattern distribution of simulation
        :param distributionData: pattern distribution of data
        :return: distance between pattern distribution of simulation and pattern distribution of given data
        '''
        dist = 0
        for k in range(4**L):
            keyData = distributionData[k] if k in distributionData else 0
            keySim = distributionSim[k] if k in distributionSim else 0
            dist += abs(keyData - keySim)
        return dist

def dist(distributionData, ditributionSim, L=1):
    """
    distance function for ABC, distance of to distributions is sum over weighted distances of all relative distances
    between the pairwise patterns from each distribution
    :param patterns: pattern distribution of simulation
    :return: distance between pattern distribution of simulation and pattern distribution of KO-file
    """
    dist = 0
    for keyData in range(4**L):
    #for keyData, valueData in distributionData.items():
        valueData = distributionData[keyData] if keyData in distributionData else 0
        for keySim in range(4**L):
        #for keySim, valueSim in ditributionSim.items():
           valueSim = ditributionSim[keySim] if keySim in ditributionSim else 0
           weight = w(keyData, keySim)/L
           dist += weight * (valueData - valueSim+1)**2
    return dist

def w(L, keyData, keySim):
        """
        computes a weight relative to two patterns
        :param keyData: pattern from file
        :param keySim: pattern from simulation
        :return: weight w
        """
        w = 0.0
        #iterations = self.L
        #while keyData != 0 or keySim != 0:
        for i in range(L):
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

def balance(L, pattern):
    mid = L/2
    r = 0
    l = 0
    for i in range(L):
        modP = pattern % 4
        for j in range(0, modP+2, 2):
            if i <= mid-1:
                l += 1
            elif i >= mid:
               r += 1
        pattern = pattern // 4
    b = (r-l)/abs(r-l) if r-l != 0 else r-l
    return b


def mahalonisDist(distributionData, distributionSim):
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

def buildWeights(L):
        #w = [[1, 2, 2, 4], [2, 1, 3, 2], [2, 3, 1, 2], [4, 2, 2, 1]]
        cov = np.zeros((4**L, 4**L))
        for i in range(0, 4**L):
            balanceI = balance(L, i)
            for j in range(i, 4**L):
                #if j != 6:
                    #continue
                balanceJ = balance(L, j)
                b = abs((balanceI-balanceJ)/2)
                weight = (w(L, i, j) + b) / L
                cov[i][j] = weight
                cov[j][i] = weight
        return cov


#data = {0:1}
#sim = {3:1}
#print(distanceFunction(data, sim))
#print(dist(data, sim))
#print(mahalonisDist(data, sim))
print(buildWeights(5))
