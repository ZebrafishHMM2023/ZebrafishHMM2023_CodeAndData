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
from utils.data_and_models import *

# from utils.format_data import format_sequences
# from utils.h5tree.hdf5_tree_viewer import h5tree_view
# from utils.load_models import *
# from utils.stubbornness import load_stubb_as_DF

# Extracting the data
extract_data("../Data/behavior_free_swimming.tar.gz")
extract_data("../Data/behavior_single_fish.tar.gz")
extract_data("../Data/neuro.tar.gz")

# Extracting the models
extract_models("../Models/hmms_20240125.tar.gz")
extract_models("../Models/longtrajectories_20240202.tar.gz")
extract_models("../Models/hmms_ARTR_20240620.tar.gz")
extract_models("../Models/gen_behavior.tar.gz")
extract_models("../Models/gen_neuro.tar.gz")


# Extracting others
extract_data("../Data/generated_neuro_MSR.tar.gz")

# Global variables

## Behavior
bDATAPATH = Path("../Data/behavior_free_swimming.h5")
blDATAPATH = Path("../Data/behavior_single_fish.h5")
bMODELPATH = Path("../Models/hmms_20240125/")
blMODELPATH = Path("../Models/longtrajectories_20240202/")
bGENPATH = Path("../Models/gen_behavior/")
DTHETA_LIM = 100
DTHETA_CMAP = cmc.vanimo
FLR_colors = DTHETA_CMAP([0.5, 0.9, 0.1])
TEMPS = np.array([18, 22, 26, 30, 33])
TEMPS_COLS = {t: c for t, c in zip(TEMPS, cmc.roma_r.resampled(len(TEMPS)).colors)}
DTHETA_THRESH = 10

## Neuro
nDATAPATH = Path("../Data/neuro.h5")
nMODELPATH = Path("../Models/hmms_ARTR_20240620/")
nOUTPATH = Path("../Models/hmms_ARTR_outputs.h5")
nGENPATH = Path("../Models/gen_neuro/")
nMSRPATH = Path("../Data/generated_neuro_MSR/")
ALL_ARTRs = ARTR_all_fish_temp(nDATAPATH)
ARTR_CMAP = plt.cm.coolwarm
COLOR_Lartr = ARTR_CMAP(0)
COLOR_Rartr = ARTR_CMAP(256)
