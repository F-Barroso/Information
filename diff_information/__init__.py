def diff_entropy(data):
    '''Computes differential entropy from one dimensional data array'''
    
    #bin width estimated from Freedman-Diaconis rule
    dx = 2*iqr(data)/(len(data)**(1/3))
    counts,bins=np.histogram(data,bins=np.round((np.max(data)-np.min(data))/dx).astype(int))
    dx = bins[1]-bins[0] #redefine bin width
    
    pdf = counts/np.sum(counts*dx)
    return -np.sum(pdf[pdf!=0]*np.log(pdf[pdf!=0]))*dx

def diff_entropydd(data):
    '''Computes differential entropy from d-dimensional data array, with d greater than 1'''
    
    #bin widths estimated from Freedman-Diaconis rule  
    dx = 2*iqr(data,axis=1)/(len(data[0])**(1/3))
    counts,bins=np.histogramdd(data.T, bins=np.round((np.max(data,axis=1)-np.min(data,axis=1))/dx).astype(int))
    dx = np.array([bins[i][1]-bins[i][0] for i in range(len(dx))])  #redefine bin widths
    
    pdf = counts/np.sum(counts*np.prod(dx))
    return -np.sum(pdf[pdf!=0]*np.log(pdf[pdf!=0]))*np.prod(dx)
    
def diff_mi(A,B,pos_save=False):
    '''Computes the differential mutual information between two 1-dimensional data arrays, A and B. - I(A,B).'''
    
    #bin widths estimated from Freedman-Diaconis rule  
    dA = 2*iqr(A)/(len(A)**(1/3))
    dB = 2*iqr(B)/(len(B)**(1/3))
    
    if dA==0 or dB==0: #integral of 0 size
        return 0.

    countA,binA=np.histogram(A, bins=np.round((np.max(A)-np.min(A))/dA).astype(int))
    if pos_save:
        posA, countA = np.where(countA!=0)[0], countA[countA!=0]
    dA = binA[1] - binA[0] #redefine bin widths
    pA = countA/np.sum(countA*dA)
    
    countB,binB=np.histogram(B, bins=np.round((np.max(B)-np.min(B))/dB).astype(int))
    if pos_save:
        posB, countB = np.where(countB!=0)[0], countB[countB!=0]
    dB = binB[1] - binB[0] #redefine bin widths
    pB = countB/np.sum(countB*dB)
    
    countAB=np.histogramdd([A,B], bins=[binA,binB])[0]
    if pos_save:
        posAB, countAB = np.array(np.where(countAB!=0)).T, countAB[countAB!=0]
    pAB = countAB/np.sum(countAB*dA*dB)
    
    if pos_save:
        return np.sum(pAB*np.log([pAB[i]/(pA[posA==posAB[i][0]][0]*pB[posB==posAB[i][1]][0]) for i in range(len(posAB))]))*dA*dB
    del countAB,countA,countB
    
    nonzero = pAB!=0
    return np.sum(pAB[nonzero]*np.log(pAB[nonzero]/np.tensordot(pA,pB,axes=0)[nonzero]))*dA*dB
