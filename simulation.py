import math
import random

from matplotlib import pyplot as plt
from scipy.optimize import minimize, basinhopping
from distribution import createDistri
from ABC import ABC

class Simulation:
    """
    simulates one cell division and methylation of DNMT1 + DNMT3A/B
    """
    upperStrand = []
    lowerStrand = []
    # if DNMT true, no knock out
    L = 0 #number of CpGs
    distances = [] #distances between CpGs
    distributionKO = [] #pattern distribution in KO-file
    patternDistriKO = [] # number of occurrences of each pattern in DNMT KO data for comparison with outcome of
    # simulation
    DNMT1KO = True
    #               rho             tau             mhy          delta
    #           (disassociation   (association  (maintenance   (de novo methyl
    #               prob)           prob)        methyl prob)       prob)
    probabilities = [[0.1,          0.8,          0.8,              0],  #DNMT1
                     [0.5,          0.5,          0,              1],   #DNMT3d (daughter strand)
                     [0.5,          0.5,          0,                1]]     #DNMT3p (parent strand)

    def __init__(self, WTfile, KOfile, distances, DNMT1KO=True):
        distributionWT, self.L, numOfPatterns = createDistri(WTfile) # distribution of methylation patterns for wildtype
        self.distributionKO, self.L, numOfPatterns = createDistri(KOfile) # distribution of methylation patterns after DNMT KO
        self.create_initial_distr(distributionWT) # create initial distribution from WT
        # number of occurrences of each pattern in DNMT KO data for comparison with outcome of simulation
        self.patternDistriKO = [i*numOfPatterns for i in self.distributionKO]
        self.distances = distances # distances between CpGs
        self.DNMT1KO = DNMT1KO

    def create_initial_distr(self, distribution):
        """
        choose an initial distribution
        :param distribution: methylation pattern distribution
        :return:
        """
        upperStrand = []
        lowerStrand = []
        #select pattern randomly
        randomdistribution = random.random()
        #print(randomdistribution)
        #start with last pattern with highest probability
        for index, item in reversed(list(enumerate(distribution))):
            if randomdistribution <= item:
                randomdistribution = index
                break;
            else:
                randomdistribution -= distribution[index]
        for i in range(self.L):
            mod = randomdistribution % 4
            if mod == 0:
                upperStrand.append(0)
                lowerStrand.append(0)
            elif mod == 1:
                upperStrand.append(1)
                lowerStrand.append(0)
            elif mod == 2:
                upperStrand.append(0)
                lowerStrand.append(1)
            elif mod == 3:
                upperStrand.append(1)
                lowerStrand.append(1)
            randomdistribution = (randomdistribution-mod) / 4
        self.upperStrand = upperStrand
        self.lowerStrand = lowerStrand

    def simulate(self, allProbs, DNMT1=True, DNMT3=True):
        '''
        simulate celldivision with methylation
        :param allProbs: methylation probabilities
        :param DNMT1: True if no DNMT1KO
        :param DNMT3: True if no DNMT3KO
        :return: strsnds after celldivisions
        '''
        #select random strand
        upperStrand = self.upperStrand if random.random() <= 0.5 else self.lowerStrand
        #decide about number of cell divisions
        celldivisions = 1 if DNMT1 and DNMT3 else 41 if DNMT3 else 26

        for c in range(celldivisions):
            lowerStrand = [0]*self.L
            if DNMT1:
                upperStrand, lowerStrand = self.simulateDNMT(allProbs[0], upperStrand, lowerStrand)
            if DNMT3:
                upperStrand, lowerStrand = self.simulateDNMT(allProbs[1], upperStrand, lowerStrand)
                lowerStrand, upperStrand = self.simulateDNMT(allProbs[2], lowerStrand, upperStrand)
            upperStrand = upperStrand if random.random() <= 0.5 else lowerStrand
        #print(upperStrand, lowerStrand)
        return upperStrand, lowerStrand

    def simulateDNMT(self, probs, parentS, daughterS):
        '''
        simulate one run of methylation enzyme
        :param probs: probability matrix
        :param parentS: template strand
        :param daughterS: strand to which enzyme is bound
        :return: both strands after methylation
        '''
        rho = probs[0]
        tau = probs[1]
        mhy = probs[2]
        delta = probs[3]
        bound = False
        CpG = True #if current position is a CpG
        length = sum(self.distances) + len(self.distances) + 1  # total number of positions
        cpgPos = 0 #cpg positions
        cpg = -1 #counting current CpG
        for pos in range(length):
            #if current position is a CpG
            if pos == cpgPos:
                cpg += 1
                cpgPos = cpgPos + self.distances[cpg] + 1 if cpg <= self.L-2 else 0 #compute next CpG
                CpG = True
            else:
                CpG = False
            #if enzyme stays associated
            if (not bound and random.random() <= tau) or bound:
                #probabilty for de novo or maintenance methylation
                if CpG and ((parentS[cpg] and random.random() <= mhy) or (not parentS[cpg] and random.random() <= delta)):
                   daughterS[cpg] = 1
                bound = True if random.random() <= 1-rho else False
        return parentS, daughterS

    def computePatternDistribution(self, probabilities):
        '''
        perform multiple iterations of simulation and compute pattern distribution of simulation
        :param probabilities: parameters for simulation
        :return: pattern distribution of simulation
        '''
        allProbs = [] * 3
        if self.DNMT1KO:
            #prob DNMT3p = DNMT3d???
            allProbs.append([])
            allProbs.append(probabilities)
            allProbs.append(probabilities)
        else:
           allProbs.append(probabilities)
        patterns = dict()
        #perform multiple iterations and store resulting patterns
        for i in range(10000):
            upperStrand, lowerStrand = self.simulate(allProbs, DNMT1=not self.DNMT1KO, DNMT3=self.DNMT1KO)
            pattern = 0
            for l in range(self.L):
                pattern += 4 ** (self.L - l - 1) * (upperStrand[l] + lowerStrand[l]*2)
            patterns[pattern] = patterns.get(pattern, 0) + 1
        patterns = {k: float(v/10000) for k, v in patterns.items()}
        return patterns

    def computeLH(self, probabilities):
        '''
        compute likelihood
        :param probabilities: parameters for simulation
        :return: negative likelihood
        '''
        #compute likelihood
        #likelihood = 1.0
        likelihood = 0.0 #for log-likelihood
        epsilon = 0.0000001  # add epsilon to all distribution values which are 0
        patterns = self.computePatternDistribution(probabilities)
        for k, v in enumerate(self.patternDistriKO):
            #if v != 0:
                simDistri = patterns[k] if k in patterns else epsilon
                likelihood += v * math.log(simDistri) #for log-likelihood
                #likelihood *= simDisrtri ** v
        print(likelihood, probabilities)
        return -likelihood

    def distanceFunction(self, patterns):
        '''
        A distance function for an ABC-algorithm, computes the distance of two distributions as sum of differences
        between all counts of patterns
        :param patterns: pattern distribution of simulation
        :return: distance between pattern distribution of simulation and pattern distribution of KO-file
        '''
        dist = 0
        for k, v in enumerate(self.distributionKO):
            dist += abs(patterns[k]-v) if k in patterns else v
        return dist

    def minimizeLH(self, probabilities):
        '''
        use parameter optimization to minimize likelihood
        :param probabilities: initial parameters
        '''
        # extra arguments passed to the objective function and its derivatives
        # bounds for parameter space
        bnds = ((0, 1), (0, 1), (0, 1), (0, 1))  # here 4 parameters bound between 0 and 1
        # use method L-BFGS-B because the problem is smooth and bounded
        #sol = minimize(self.computeLH, probabilities, method='L-BFGS-B', bounds=bnds, options={'disp': True})
        minimizer_kwargs = dict(method="L-BFGS-B", bounds=bnds, options={'disp': True})
        sol = basinhopping(self.computeLH, probabilities, minimizer_kwargs=minimizer_kwargs)
        print(sol)


#DNMT1KO:
sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT1KO.txt", [13, 14])
#DNMT3KO:
#sim = Simulation("Daten/ySatWTJ1C.txt", "Daten/ySatDNMT3abKO.txt", [13, 14], False)

#likelihood computation
#sim.minimizeLH(sim.probabilities[0])
#sim.minimizeLH(sim.probabilities[1])
#sim.minimizeLH([0,0.5,0.5,0])

#ABC
abc = ABC()
abc.abc(sim.computePatternDistribution, sim.distanceFunction, 0.8)
plt.plot(sim.distributionKO)
plt.xlabel('pattern value')
plt.ylabel('distribution value')
plt.title('pattern distribution of data')
plt.show()
