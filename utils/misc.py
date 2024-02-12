import matplotlib.pyplot as pl
import numpy as np


def confusion(seqs1, seqs2, n1, n2):
    """Compute the confusion matrix between two labbelled sequences."""
    C = np.zeros((n1, n2))
    for i in range(len(seqs1)):
        b3, b4 = seqs1[i], seqs2[i]
        C[b3, b4] += 1
    C = C / C.sum(axis=1)[:, np.newaxis]
    return C


def plot_confusion(ax, C, labs1, labs2):
    """Plot the confusion matrix."""
    h = ax.imshow(C.T, cmap="cividis", vmin=0, vmax=1)
    ax.set_xticks(range(len(labs1)), labs1)
    ax.set_yticks(range(len(labs2)), labs2)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            ax.text(i, j, f"{C[i,j]:0.2f}", ha="center", va="center")
    ax.set_aspect("equal")
    return h


def find_single_streaks(bs, states, lenghts):
    """Find the streaks of a single sequence.

    Parameters
    ----------
    bs : array
        The sequence of states.
    states : array
        The states to look for.
    lenghts : array
        The lenghts of streaks to look for.

    Returns
    -------
    Ns : array
        The number of streaks of each length.
    """
    Ns = np.empty((len(states), len(lenghts) - 1))
    for i, s in enumerate(states):
        mask = bs == s
        streak_lengths = np.diff(
            np.where(np.concatenate(([mask[0]], mask[:-1] != mask[1:], [True])))[0]
        )[::2]
        n, _ = np.histogram(streak_lengths, lenghts)
        Ns[i, :] = n
    return Ns


def find_streaks(bs, states, lenghts):
    """Find the streaks of a sets of sequences.

    Parameters
    ----------
    bs : list of arrays
        The sequences of states.
    states : array
        The states to look for.
    lenghts : array
        The lenghts of streaks to look for.

    Returns
    -------
    Ns : array
        The number of streaks of each length.
    """
    Ns = np.c_[[find_single_streaks(b, states, lenghts) for b in bs]]
    Ns = Ns.sum(axis=0)
    Ns = Ns / Ns.sum()
    return Ns


def streak_decay_const(streak_lenghts, lenghts, first_n=7):
    """Estimate the decay constant of the streaks.

    Parameters
    ----------
    streak_lenghts : array
        The number of streaks of each length.
    lenghts : array
        The lenghts of streaks to look for.
    first_n : int
        The number of points to use for the fit.

    Returns
    -------
    L0F : float
        The decay constant of the forward streaks.
    L0T : float
        The decay constant of the turning streaks.
    bF : float
        The intercept of the forward streaks.
    bT : float
        The intercept of the turning streaks.
    """
    # forward
    AF, BF = np.polyfit(lenghts[:first_n], np.log(streak_lenghts[0, :first_n]), 1)
    L0F, bF = -1 / AF, np.exp(BF)

    # turns
    AT, BT = np.polyfit(
        lenghts[:first_n], np.log(streak_lenghts[1:, :first_n].mean(axis=0)), 1
    )
    L0T, bT = -1 / AT, np.exp(BT)
    return L0F, L0T, bF, bT
