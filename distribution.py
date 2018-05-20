import numpy as np


def dataconv_general(file):
    #open file
    fobj = open(file, 'r')
    
    first_line = fobj.readline()
    L=len(first_line.split())
    #initialize emtpy array with correct dimension
    arr = np.array([]).reshape(0,L)
    # go through file line by line
    for line in fobj:
        #reshaping lines
        line = line.strip()
        columns = line.split()
        #omit rows with x entries (missing read)
        if ('X' in columns) ==False and ('x' in columns) ==False:
            #convert strings to int
            columns=list(map(int,columns))
            #attach to output vector
            arr=np.vstack([arr,columns])
    fobj.close()
    #convert entries to use base-4 numbers
    arr[arr==3]=2
    arr[arr==4]=3
    
    #calculate index of pattern
    cod=np.zeros(len(arr))
    for i in range(0,L):
        cod=cod+arr[:,i]*4**(L-i-1)
    cod=cod.astype(int)            
    
    return [arr, cod]
    
    
def createDistri(file1):

    #convert data
    [ig1, cod]  = dataconv_general(file1)   

    L=len(ig1[0])
    #initialize distribution
    distribution=np.zeros(4**L)

    #calculate distribution from data
    for i in range(0,len(cod)):
        aux=cod[i];
        distribution[aux]+=1
    
    
    #normalization of distribution
    distribution=distribution/len(cod)
    
    return distribution, L
