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


def w(keyData, keySim, L):
        """
        computes a weight relative to two patterns
        :param keyData: pattern from file
        :param keySim: pattern from simulation
        :return: weight w
        """
        w = []
        blocksData = []
        blocksSim = []
        mismatches = []
        #iterations = self.L
        #while keyData != 0 or keySim != 0:
        for i in range(L):
            #iterations -= 1
            modData = keyData % 4
            modSim = keySim % 4
            if modData == modSim:
                w.append(1.0)
            else:
                mismatches.append(i)
                if modData+modSim == 3: #then the methylation state at both strands is complementary
                    w.append(3.0)
                elif modData != modSim: #one position of methylation pattern different
                    w.append(2.0)
            if modData != 0:
                blocksData.append(i)
            if modSim != 0:
                blocksSim.append(i)
            keyData = (keyData-modData) / 4
            keySim = (keySim-modSim) / 4
        #w += iterations
        return w, blocksData, blocksSim, mismatches


def balance(pattern, L):
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


def blocks(blocksData, blocksSim, L):
        c = 0
        if len(blocksData) != 0:
            for i in blocksSim:
                mini = min([abs(j - i) for j in blocksData])
                c += mini**2
        return c/L


def c(mismatches, L):
        if len(mismatches) == L:
            return [L]*L
        else:
            indices = [i for i in range(L)]
            c = []
            matches = [j for j in indices if j not in mismatches]
            for i in range(L-1, -1, -1):
                if i in mismatches:
                    c.append(min([abs(i - j) for j in matches])+1)
                else:
                    c.append(1)
            return c


def mahalonisDist(distributionData, distributionSim, L):
        X = [distributionData[k] if k in distributionData else 0 for k in range(4**L)]
        Y = [distributionSim[k] if k in distributionSim else 0 for k in range(4**L)]
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


def invert(cov):
        try:
            inv = np.linalg.inv(cov)
        except np.linalg.LinAlgError:
            #cov[0][0] += 0.001
            #self.invert(cov)
            pass
        else:
            pass
        return cov


def buildWeights(L):
        #w = [[1, 2, 2, 4], [2, 1, 3, 2], [2, 3, 1, 2], [4, 2, 2, 1]]
        cov = np.zeros((4**L, 4**L))
        for i in range(0, 4**L):
            #balanceI = self.balance(i)
            for j in range(i, 4**L):
                #balanceJ = self.balance(j)
                #balance = abs((balanceI-balanceJ)/2)
                W, blocksData, blocksSim, mismatches = w(i, j, L)
                C = c(mismatches, L)
                W = [x*y for x, y in zip(W, C)]
                W = sum(W) / L
                #W = (W + blocks(blocksData, blocksSim, L)) / L
                cov[i][j] = W
                cov[j][i] = W
        return cov


#data = {0:1}
#sim = {3:1}
#print(distanceFunction(data, sim))
#print(dist(data, sim))
#print(mahalonisDist(data, sim))
print(buildWeights(5))
