import numpy as np
from tqdm import tqdm

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

def MSR_q(trajs, q):
    return np.mean(np.concatenate([[ np.abs(np.sum(traj[i:i+q+1]))**2 for i in range(len(traj)-q)] for traj in trajs]))
def MSR(trajs, qs):
    mu = np.mean(np.concatenate(trajs))
    trajs_mu = [traj - mu for traj in trajs]
    return np.array([MSR_q(trajs_mu, q) for q in tqdm(qs)])