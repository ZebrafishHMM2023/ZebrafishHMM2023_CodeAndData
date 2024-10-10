import numpy as np
from tqdm import tqdm


def MSR_q(seq, q):
    """Compute the Mean Squared Reorientation of a sequence.
    
    Parameters : 
    ------------
    seq : 1D array of floats
        sequence of re-orientations angles.
    q : integer
        the lag.

    Return :
    --------
    m_q : float
        the computed M_q (see math function above)
    """
    b = []
    for n in range(len(seq)-q+1):
        a = 0
        for i in range(q):
            a += seq[n+i]
        b.append(a**2)
    m_q = np.mean(b)
    return m_q

def MSR(seq, qs):
    """Compute the MSR of a sequence at multiple lags.
    
    Parameters : 
    ------------
    seq : 1D array of floats
        sequence of re-orientations angles.
    qs : 1D array of integers
        the lags.

    Return :
    --------
    m_qs : 1D array
        the computed M_q's (see math function above)
    """
    m_qs = np.empty_like(qs, dtype=np.float_)
    for i,q in enumerate(qs):
        m_qs[i] = MSR_q(seq, q)

    return m_qs

def MSRs(seqs, qs):
    """Compute the MSR of multiple sequences at multiple lags.
    
    Parameters : 
    ------------
    seqs : liste of 1D array of floats
        sequences of re-orientations angles.
    qs : 1D array of integers
        the lags.

    Return :
    --------
    m_qs : 2D array
        the computed M_q's (see math function above)
    """
    m_qs = np.empty((len(seqs),len(qs)), dtype=np.float_)
    for i in tqdm(range(len(seqs)), desc="Computing MSRs"):
        m_qs[i] = MSR(seqs[i], qs)

    return m_qs

def randomize_sequences(X):
    lenghts = [len(x) for x in X]
    long_seq = np.concatenate(X)
    Y = np.random.choice(long_seq, size=(len(X),np.mean(lenghts).astype(int)))
    #y = np.random.permutation(np.concatenate(X))
    #cumlens = np.insert(np.cumsum(lenghts), 0,0)
    #Y = []
    #for i in range(1,len(cumlens)):
    #    m,n = cumlens[i-1], cumlens[i]
    #    Y.append(y[m:n])
    return Y