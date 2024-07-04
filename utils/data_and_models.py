import tarfile
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat


def extract_models(path="../Models/hmms_20240125.tar.gz"):
    """
    Extracts the models from the tar.gz file
    :param path: path to the tar.gz file
    :return: path to extracted models
    """
    path = Path(path)
    to = path.parent.joinpath(Path(path.stem).stem)
    if "mixture" in path.name:
        to = to.with_suffix(".h5")
    if to.exists():
        return to
    else:
        print(f"Extracting models from {path} to {to}... ")
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(path=path.parent)
        return to


def extract_data(path="../Data/behavior_free_swimming.tar.gz"):
    """
    Extracts the data from the tar.gz file
    :param path: path to the tar.gz file
    :return: path to extracted data
    """
    path = Path(path)
    to = path.parent.joinpath(Path(path.stem).stem).with_suffix(".h5")
    if to.exists():
        return to
    else:
        print(f"Extracting data from {path} to {to}... ")
        with tarfile.open(path, "r:gz") as tar:
            tar.extractall(path=path.parent)
        return to


def format_sequences(sequences):
    """Re-format data of reorientation sequences to remove NaNs.

    In the data, sequences are postpended with NaNs when the tacking was lost.
    This function removes the NaNs.

    Parameters :
    ------------
    :sequences: 2d array ( N_sequences x N_bouts )
        sequences of reorientation angles, padded with NaNs

    Return :
    --------
    :seqs: list of 1d arrays (N_sequences) (n_bouts_per_sequence,)
        list of sequences. Sequence have different number of bouts.
    """
    seqs = [seq[~np.isnan(seq)] for seq in sequences]
    return seqs


def load_sequences(path, temp):
    """Load reorientation sequences from hdf5 file at specific temperature.

    Parameters :
    ------------
    :path: str
        path to the hdf5 file
    :temp: int
        temperature at which the sequences were recorded

    Return :
    --------
    :X: list of 1d arrays (N_sequences) (n_bouts_per_sequence,)
        list of sequences. Sequence have different number of bouts.

    """
    file = h5py.File(path, "r")
    X = file[f"/behaviour/{temp}/dtheta"][:]
    file.close()
    return format_sequences(X)


def load_full_traj(path, temp, traj_i):
    """Load full trajectory from hdf5 file at specific temperature.

    Parameters :
    ------------
    :path: str
        path to the hdf5 file
    :temp: int
        temperature at which the sequences were recorded
    :traj_i: int
        index of the trajectory to load

    Return :
    --------
    """
    file = h5py.File(path, "r")
    dtheta = file[f"/behaviour/{temp}/dtheta"][traj_i, :]
    xpos = file[f"/behaviour/{temp}/xpos"][traj_i, :]
    ypos = file[f"/behaviour/{temp}/ypos"][traj_i, :]
    file.close()
    return dtheta[np.isfinite(dtheta)], xpos[np.isfinite(xpos)], ypos[np.isfinite(ypos)]


def model_name(model, temp):
    return f"hmm-{model}-sym-T{temp}"


def load_transmat(models_path, model, temp):
    """Load transition matrix from hdf5 file at specific temperature."""
    path = models_path.joinpath(model_name(model, temp) + ".hdf5")
    f = h5py.File(path, "r")
    T = f["transition_matrix"][:]
    f.close()
    return T.T


def load_transmat_multi(models_path, model, temp):
    """Load transition matrix from hdf5 file at specific temperature."""
    paths = models_path.glob(model_name(model, temp) + "-rep*.hdf5")
    Ts = []
    for path in paths:
        f = h5py.File(path, "r")
        T = f["transition_matrix"][:]
        Ts.append(T.T)
        f.close()
    return np.c_[Ts]


def load_emission(models_path, model, temp):
    """Load emission parameters from hdf5 file at specific temperature."""
    path = models_path.joinpath(model_name(model, temp) + ".hdf5")
    f = h5py.File(path, "r")
    try:
        F = f["σforw"][:]
    except KeyError:
        F = f["forw"][:]
    return {"forward": F, "turn": f["turn"][:]}


def load_emission_multi(models_path, model, temp):
    """Load emission parameters from hdf5 file at specific temperature."""
    paths = models_path.glob(model_name(model, temp) + "-rep*.hdf5")
    Fs = []
    Ts = []
    for path in paths:
        f = h5py.File(path, "r")
        try:
            Fs.append(f["σforw"][:])
        except KeyError:
            Fs.append(f["forw"][:])
        Ts.append(f["turn"][:])
    return {"forward": np.c_[Fs], "turn": np.c_[Ts]}


def load_viterbi(models_path, model, temp):
    """Load viterbi sequences from txt file at specific temperature."""
    path = models_path.joinpath(model_name(model, temp) + "-viterbi.txt")
    seqs = []
    with open(path, "r") as txt:
        lines = txt.readlines()
        for line in lines:
            x = (np.fromstring(line, sep="\t") - 1).astype(np.int_)
            seqs.append(x)
    return seqs


def load_genbouts(models_path, model, temp):
    """Load generated bouts from txt file at specific temperature."""
    path = models_path.joinpath(model_name(model, temp) + "-simulated-bouts.txt")
    seqs = []
    with open(path, "r") as txt:
        lines = txt.readlines()
        for line in lines:
            seqs.append(np.fromstring(line, sep="\t"))
    return seqs


def load_LLHs(models_path, model, temp, set="tests"):
    """Load log-likelihoods from txt file at specific temperature."""
    paths = models_path.glob(model_name(model, temp) + f"-rep*-{set}-lls.txt")
    LLHs = []
    for path in paths:
        LLHs.append(np.loadtxt(path))
    return np.concatenate(LLHs)


def ARTR_all_fish_temp(path):
    file = h5py.File(path, "r")
    all_fish_temp = file["fish_temp_combinations"][()]
    file.close()
    return all_fish_temp


def load_ARTR(datapath, fish, temp):
    """Load spikes and magnetization of ARTR."""
    file = h5py.File(datapath, "r")
    grp = file[f"T{temp}-fish{fish}"]
    dt = grp["dt"][()]
    L = grp["L"][()]
    R = grp["R"][()]
    file.close()
    mL = L.mean(axis=1)
    mR = R.mean(axis=1)
    return mL, mR, L, R, dt


def load_ARTR_magnet(datapath, fish, temp):
    """Load magnetization of ARTR."""
    mL, mR, _, _, _ = load_ARTR(datapath, fish, temp)
    return np.c_[mL, mR]


def load_ARTR_viterbi(path, fish, temp):
    """Load viterbi state sequence of ARTR."""
    file = h5py.File(path, "r")
    grp = file[f"T{temp}-fish{fish}"]
    state_seq = grp["viterbi"][()]
    file.close()
    return state_seq


def load_ARTR_gen(path, fish, temp):
    """Load generated state sequence of ARTR."""
    file = h5py.File(path, "r")
    grp = file[f"T{temp}-fish{fish}"]
    gen_seq = grp["gen_states"][()]
    gen_L = grp["gen_L"][()]
    gen_R = grp["gen_R"][()]
    file.close()
    return gen_seq, gen_L, gen_R


def load_ARTR_fields(path, fish, temp):
    """Load the field of ARTR HMMs."""
    name = f"artr_temperature={temp}-fish={fish}.hdf5"
    path = path.joinpath(name)
    file = h5py.File(path, "r")
    h = file["h"][()]
    file.close()
    return h


def load_ARTR_dt(path, fish, temp):
    """Load the field of ARTR HMMs."""
    name = f"artr_temperature={temp}-fish={fish}.hdf5"
    path = path.joinpath(name)
    file = h5py.File(path, "r")
    dt = file["time_unit"][()]
    file.close()
    return dt


def load_ARTR_sejourn(path, fish, temp):
    """Load the field of ARTR HMMs."""
    name = f"artr_temperature={temp}-fish={fish}.hdf5"
    path = path.joinpath(name)
    file = h5py.File(path, "r")
    stF = file["sojourn_times_F"][()]
    stL = file["sojourn_times_L"][()]
    stR = file["sojourn_times_R"][()]
    file.close()
    return stF, stL, stR


def load_ARTR_transmat(path, fish, temp):
    """Load the field of ARTR HMMs."""
    name = f"artr_temperature={temp}-fish={fish}.hdf5"
    path = path.joinpath(name)
    file = h5py.File(path, "r")
    T = file["transition_matrix"][()]
    file.close()
    return T.T
