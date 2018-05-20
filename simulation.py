import random

from distribution import createDistri


class Simulation:
    """
    simulates one cell division and methylation of DNMT1 + DNMT3A/B
    """
    upperStrand = []
    lowerStrand = []
    # if DNMT true, no knock out
    DNMT1 = True
    DNMT3 = True
    L = 0 #number of CpGs
    #               rho             tau             mhy          delta
    #           (disassociation   (association  (maintenance   (de novo methyl
    #               prob)           prob)        methyl prob)       prob)
    probabilities = [[0.1,          0.8,          0.8,              0.05],  #DNMT1
                     [0.5,          0.5,          0.5,              0.5],   #DNMT3d (daughter strand)
                     [0.5,          0.5,          1,                1]]     #DNMT3p (parent strand)

    def __init__(self, filename, DNMT1=True, DNMT3=True):
        self.DNMT1 = DNMT1
        self.DNMT3 = DNMT3

        distribution, self.L = createDistri(filename) #distribution of methylation patterns
        self.create_initial_distr(distribution)

        self.simulate(self.probabilities)

    def create_initial_distr(self, distribution):
        """
        choose an initial distribution
        :param distribution: methylation pattern distribution
        :return:
        """
        #select pattern randomly
        randomdistribution = random.random()
        print(randomdistribution)
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
                self.upperStrand.append(0)
                self.lowerStrand.append(0)
            elif mod == 1:
                self.upperStrand.append(1)
                self.lowerStrand.append(0)
            elif mod == 2:
                self.upperStrand.append(0)
                self.lowerStrand.append(1)
            elif mod == 3:
                self.upperStrand.append(1)
                self.lowerStrand.append(1)
            randomdistribution = (randomdistribution-mod) / 4
        print(self.upperStrand)
        print(self.lowerStrand)

    def simulate(self, allProbs):
        '''
        simulate one cell division with methylation
        :param allProbs: methylation probabilities
        :return: 
        '''
        #select random strand
        self.upperStrand = self.upperStrand if random.random() <= 0.5 else self.lowerStrand
        self.lowerStrand = [0]*self.L
        print(self.upperStrand)
        print(self.lowerStrand)
        # which one to do first?
        if self.DNMT1:
            self.upperStrand, self.lowerStrand = self.simulateDNMT(allProbs[0], self.upperStrand, self.lowerStrand)
        if self.DNMT3:
            self.upperStrand, self.lowerStrand = self.simulateDNMT(allProbs[1], self.upperStrand, self.lowerStrand)
            self.lowerStrand, self.upperStrand = self.simulateDNMT(allProbs[2], self.lowerStrand, self.upperStrand)
        print(self.upperStrand)
        print(self.lowerStrand)

    def simulateDNMT(self, probs, parentS, daughterS):
        '''
        simulate one run of methylation enzyme
        :param probs: probability matrix
        :param parentS: template strand
        :param daughterS: strand to which enzyme is bound
        :return: both strands after methylation
        '''
        bound = True
        rho = probs[0]
        tau = probs[1]
        mhy = probs[2]
        delta = probs[3]
        for pos in range(self.L):
            # which first - methylation or disassociation? can DNMT disassociate before methylation
            #if enzyme stays associated
            if bound and random.random() <= 1-rho:
                #probabilty for de novo or maintenance methylation
                if (parentS[pos] and random.random() <= mhy) or (not parentS[pos] and random.random() <= delta):
                    daughterS[pos] = 1
            else:
                bound = True if random.random() <= tau else False
        return parentS, daughterS


Simulation("Daten\IAPWTJ1C.txt")
