# This file contains the setup for the figures in the paper

# Graphical setup
import matplotlib.pyplot as plt

plt.style.use("../style.mplstyle")
# plt.rcParams.update({"text.usetex": True})

# importing the necessary libraries
import sys
from pathlib import Path

import cmcrameri.cm as cmc
import h5py
import numpy as np

# importing local tools
sys.path.append("../")
import utils.MarkovChains as MC
# from utils.confusion import *
from utils.data_and_models import (extract_data, extract_models, load_emission,
                                   load_emission_multi, load_full_traj,
                                   load_genbouts, load_LLHs, load_sequences,
                                   load_transmat, load_transmat_multi,
                                   load_viterbi)

# from utils.format_data import format_sequences
# from utils.h5tree.hdf5_tree_viewer import h5tree_view
# from utils.load_models import *
# from utils.stubbornness import load_stubb_as_DF

# Extracting the data
extract_data("../Data/behavior_free_swimming.tar.gz")

# Extracting the models
extract_models("../Models/hmms_20240125.tar.gz")
extract_models("../Models/longtrajectories_20240202.tar.gz")

# Global variables
DATAPATH = Path("../Data/behavior_free_swimming.h5")
MODELPATH = Path("../Models/hmms_20240125/")
DTHETA_LIM = 100
DTHETA_CMAP = cmc.vanimo
FLR_colors = DTHETA_CMAP([0.5, 0.9, 0.1])
TEMPS = np.array([18, 22, 26, 30, 33])
TEMPS_COLS = {t: c for t, c in zip(TEMPS, cmc.roma_r.resampled(len(TEMPS)).colors)}
DTHETA_THRESH = 10
