import numpy as np


def threshold_classifier(dtheta_seqs, left_threshold, right_threshold):
    """Classify re-orientation angles into 3 states based on 2 thresholds.

    Re-orientation angle δθ are classifed as followed :
        - δθ < :left_threshold:  → 1
        - δθ > :right_threshold: → 2
        - :left_threshold: < δθ < :right_threshold: → 0

    Parameters :
    ------------
    dtheta_seqs : list of 1d arrays of floats
        sequences of re-orientation angles
    left_threshold : number
        lower threshold
    right_threshold : number
        upper threshold

    Returns :
    ---------
    bout_seqs : list of 1d arrays of ints
        sequence of states.
    """

    def single_classifer(seq):
        bout_seq = np.zeros_like(seq, dtype=np.int8)
        idx_left = np.where(seq < left_threshold)
        idx_right = np.where(seq > right_threshold)
        bout_seq[idx_left] = 1
        bout_seq[idx_right] = 2
        return bout_seq

    bout_seqs = [single_classifer(seq) for seq in dtheta_seqs]
    return bout_seqs


def bout_proba(bout_seqs, n_state=3):
    """Compute the a priori probability of each state.

    Parameters :
    ------------
    bout_seqs : list of 1d arrays of ints
        sequence of states. 0='Forward', 1='Left', 2='Right'
    n_state : int
        number of states

    Returns :
    ---------
    p_bout : array of 3 float
        probability of each state 'Forward', 'Left', 'Right'
    """

    all = np.concatenate(bout_seqs)
    h, _ = np.histogram(all, bins=np.arange(n_state + 1))
    p_bout = h / h.sum()
    return p_bout


def bout_transitions(bout_seqs, n_state=3):
    """Compute the transition matrix between states.

    Parameters :
    ------------
    bout_seqs : list of 1d arrays of ints
        sequence of states. 0='Forward', 1='Left', 2='Right'
    n_state : int
        number of states

    Returns :
    ---------
    P : array of shape (n_state, n_state)
        transition matrix between states
    """

    P = np.zeros((n_state, n_state))
    for bout_seq in bout_seqs:
        for i in range(len(bout_seq) - 1):
            P[bout_seq[i], bout_seq[i + 1]] += 1
    P = P / P.sum(axis=1, keepdims=True)
    return P


def steady_state(T):
    """taken from https://stackoverflow.com/questions/52137856/steady-state-probabilities-markov-chain-python-implementation"""
    n = T.shape[0]
    q = np.c_[T - np.eye(n), np.ones(n)]
    QTQ = np.dot(q, q.T)
    return np.linalg.solve(QTQ, np.ones(n))


def plot_transition_matrix(
    ax,
    T,
    labels=["forward", "left", "right"],
    xlabel="Current step",
    ylabel="Next step",
):
    """Routine to plot transition matrices.

    Parameters :
    ------------
    ax : matplotlib.Axes
        matplotlib axis where to plot the transition matrix.
    T : 2D square array
        transition probability matrix.
    labels : list of strings
        labels for each state of the transition matrix.
        Must have the same lenght as :T:.
        default : ["forward", "left", "right"]
    xlabel : str
        X axis label
    ylabel : str
        Y axis label
    """
    assert T.ndim == 2, f":T: must be 2D. You gave {T.ndim}D."
    assert T.shape[0] == T.shape[1], f":T: must be a square matrix. You gave {T.shape}."
    assert (
        (T >= 0) * (T <= 1)
    ).all(), f":T: must be a matrix of probabilities. Some values you gave were not between 0 and 1."
    h = ax.imshow(T.T, vmin=0, vmax=1, cmap="plasma", origin="lower", aspect="equal")
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            ax.text(i, j, f"{T[i,j]:0.2f}", ha="center", va="center")
    l = np.arange(len(labels))
    ax.set_xticks(l, labels)
    ax.set_yticks(l, labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return h


def center_bouts(bout_seqs):
    """Center the sequence of states."""
    return [np.where(seq == 2, -1, seq) for seq in bout_seqs]


def rolling_window(a, window):
    """Create a rolling window view of the input 1D array.

    taken from https://stackoverflow.com/a/6811241"""
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def find_stubborn(centered, q=1):
    """Find stubborn sequences in a centered sequence of states."""
    rolled = rolling_window(centered, q + 2)
    cond_turn_bigin = rolled[:, 0] != 0
    cond_turn_end = rolled[:, -1] != 0
    cond_forward_middle = (rolled[:, 1:-1] == 0).all(axis=1)
    cond_stubborn = cond_turn_bigin * cond_turn_end * cond_forward_middle
    cond_switched = rolled.sum(axis=1) == 0
    inds_switched = np.where(cond_switched * cond_stubborn)[0]
    inds_stayed = np.where(~cond_switched * cond_stubborn)[0]
    return inds_switched, inds_stayed


def stubbornness_factor(centereds, q=1):
    """Compute the stubbornness factor of a sequence of states."""
    changed = 0
    stayed = 0
    for seq in centereds:
        i_sw, i_st = find_stubborn(seq, q=q)
        changed += len(i_sw)
        stayed += len(i_st)

    N = changed + stayed
    try:
        f = stayed / changed
        err = (
            1 / np.sqrt(N) * f * (np.sqrt(changed / stayed) + np.sqrt(stayed / changed))
        )
    except ZeroDivisionError:
        f = np.nan
        err = np.nan
    return f, err
