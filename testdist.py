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

def w(keyData, keySim):
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
            w += 3.0
        elif modData != modSim: #one position of methylation pattern different
            w += 2.0
        else:
            w += 1.0
        keyData = (keyData-modData) / 4
        keySim = (keySim-modSim) / 4
    return w

def mahalonisDist(distributionData, distributionSim, L=1):
        X = [distributionData[k] if k in distributionData else 0 for k in range(4**L)]
        Y = [distributionSim[k] if k in distributionSim else 0 for k in range(4**L)]
        mhy = np.average(X, weights=range(0,4**L))
        ny = np.average(Y, weights=range(0,4**L))
        z = list(zip([x-mhy for x in X], [y-ny for y in Y]))
        cov = np.cov(z)
        subtr = np.subtract(X, Y)
        d = np.inner(np.dot(subtr, cov), subtr)
        return math.sqrt(d)


data = {0:1}
sim = {3:1}
print(distanceFunction(data, sim))
print(dist(data, sim))
print(mahalonisDist(data, sim))
