import math
import random

from scipy.optimize import minimize
from distribution import createDistri


class Simulation:
    """
    simulates one cell division and methylation of DNMT1 + DNMT3A/B
    """
    upperStrand = []
    lowerStrand = []
    # if DNMT true, no knock out
    L = 0 #number of CpGs
    #               rho             tau             mhy          delta
    #           (disassociation   (association  (maintenance   (de novo methyl
    #               prob)           prob)        methyl prob)       prob)
    probabilities = [[0.1,          0.8,          0.8,              0],  #DNMT1
                     [0.5,          0.5,          0,              1],   #DNMT3d (daughter strand)
                     [0.5,          0.5,          0,                1]]     #DNMT3p (parent strand)

    def __init__(self, filename, DNMT1=True, DNMT3=True):
        distribution, self.L, numOfPatterns = createDistri(filename) #distribution of methylation patterns
        self.create_initial_distr(distribution)
        dist = [i*numOfPatterns for i in distribution]

        self.minimizeLH(self.probabilities[0], dist, False)
        self.minimizeLH(self.probabilities[1], dist, True)

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
        bound = True if random.random() <= 1-rho else False
        for pos in range(self.L):
            #if enzyme stays associated
            if bound:
                #probabilty for de novo or maintenance methylation
                if (parentS[pos] and random.random() <= mhy) or (not parentS[pos] and random.random() <= delta):
                    daughterS[pos] = 1
                bound = True if random.random() <= 1-rho else False
            else:
                bound = True if random.random() <= tau else False
        return parentS, daughterS

    def computeLH(self, probabilities, distribution, DNMT1KO, iterations=10000):
        '''
        perform multiple iterations of simulation and compute likelihood
        :param probabilities: parameters for simulation
        :param distribution: original pattern distribution
        :param DNMT1KO: True if DNMT1KO, False if DNMT3KO
        :param iterations: number of iterations of simulation
        :return: negative likelihood
        '''
        allProbs = [] * 3
        if DNMT1KO:
            #prob DNMT3p = DNMT3d???
            allProbs.append([])
            allProbs.append(probabilities)
            allProbs.append(probabilities)
        else:
           allProbs.append(probabilities)
        patterns = dict()
        likelihood = 0

        #perform multiple iterations and store resulting patterns
        for i in range(iterations):
            upperStrand, lowerStrand = self.simulate(allProbs, DNMT1=not DNMT1KO, DNMT3=DNMT1KO)
            pattern = 0
            for l in range(self.L):
                pattern += 4 ** (self.L - l - 1) * (upperStrand[l] + lowerStrand[l]*2)
            patterns[pattern] = patterns.get(pattern, 0) + 1
        patterns = {k: float(v/iterations) for k, v in patterns.items()}

        #compute likelihood
        for k, v in patterns.items():
        #if v is 1.0 here, likelihood is 0?!
            likelihood += distribution[k] * math.log(v)
        #print(-likelihood, probabilities)
        return -likelihood #minimize negative log-Likelihood <=> maximize (log)Likelihood

    def minimizeLH(self, probabilities, distribution, DNMT1KO):
        '''
        use parameter optimization to minimize likelihood
        :param probabilities: initial parameters
        :param distribution: original pattern distribution
        :param DNMT1KO: True if DNMT1KO, False if DNMT3KO
        '''
        #extra arguments passed to the objective function and its derivatives
        args = (distribution, DNMT1KO)
        #bounds for parameter space
        bnds = ((0, 1), (0, 1), (0, 1), (0, 1))  # here 4 parameters bound between 0 and 1
        sol = minimize(self.computeLH, probabilities, args=args, bounds=bnds, options={'disp': True})
        #sol = minimize(self.computeLH, probabilities, method='L-BFGS-B', args=args, bounds=bnds, options={'disp': True})
        #sol = minimize(self.computeLH, probabilities, method='SLSQP', args=args, bounds=bnds, options={'disp': True})
        print(sol)


Simulation("Daten/IAPDnmt1KO.txt", DNMT1=False)
