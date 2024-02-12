# Structure and individuality of navigation in zebrafish larvae

This repository contains the data and code to reproduce the figures from :   
*Mattéo Dommanget-Kott, Jorge Fernandez-De-Cossio-Diaz, Monica Coraggioso, Volker Bormuth, Rémi Monasson, Georges Debrégeas, and Simona Cocco*. ‘**Structure and Individuality of Navigation in Zebrafish Larvae**’, February 2024. [https://hal.science/hal-04445557](https://hal.science/hal-04445557).

![Fig2c](Figures/panels/Fig2/example_labeling_part2.svg)

## Layout of this Repository

This repo is structured as follows :

```
.                                          descriptions
├── 📁Figures                              ├── 📁 Making and storing figures
│   ├── 📄setup.py                         │   ├── 📄 basic setup and imports for the notebooks
│   ├── 📔Fig1.ipynb                       │   ├── 📔 ipython notebook to make panels of Fig1
│   ├── 📔Fig2.ipynb                       │   ├── 📔 ipython notebook to make panels of Fig2
│   ├── 📔Fig3.ipynb                       │   ├── 📔 ipython notebook to make panels of Fig3
│   ├── 📔Fig4.ipynb                       │   ├── 📔 ipython notebook to make panels of Fig4
│   ├── 📔Fig5.ipynb                       │   ├── 📔 ipython notebook to make panels of Fig5
│   ├── 📔MixtureModel.ipynb               │   ├── 📔 Julia ipython notebook to compute mixture models of Fig2
│   ├── 📁panels                           │   ├── 📁 storing figure panels in svg format
│   │   ├── 📁Fig1                         │   │   ├── 📁 storing panels for Fig1
│   │   │   ├── 📊panel_name.svg           │   │   │   ├── 📊 ...
│   │   │   ├── 📊...                      │   │   │   ├── 📊 ...
│   │   ├── 📁Fig2                         │   │   ├── 📁 storing panels for Fig2
│   │   │   ├── 📊panel_name.svg           │   │   │   ├── 📊 ...
│   │   │   ├── 📊...                      │   │   │   ├── 📊 ...
│   │   ...
├── 📁Data                                 ├── 📁 storing datasets
│   ├── 💾behavior_free_swimming.tar.gz    │   ├── 💾 freely swimming trajectories (multiple fish)
├── 📁Models                               ├── 📁 storing HMM models and associated info
│   ├── 💾hmms_20240125.tar.gz             │   ├── 💾 models for multi-fish swiming data
│   ├── 💾longtrajectories_20240202.tar.gz │   ├── 💾 models for single-fish swiming data
│   ├── 💾mixtures.tar.gz                  │   ├── 💾 mixture models
├── 📁utils                                ├── 📁 usefull functions and routines
│   ├── 📄data_and_models.py               │   ├── 📄 functions for loading data + models related stuff
│   ├── 📄MarkovChains.py                  │   ├── 📄 Markov Chains related functions
│   ├── 📄misc.py                          │   ├── 📄 random bits of usefull stuff
├── 📄style.mplstyle                       ├── 📄 matplotlib style
├── 📄README.md                            ├── 📄 the repo readme
├── 📄LICENSE                              ├── 📄 terms of the GNU GPL v3 license
├── 📄ZebrafishHMM2024.env                 ├── 📄 anaconda explicit environnement
```

## Data availability

The data comes from the paper [Thermal modulation of Zebrafish exploratory statistics reveals constraints on individual behavioral variability, Le Goc et al. 2021, BMC Biol](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01126-w).  
The original dataset can be directly downloaded [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.3r2280ggw).

A re-organised version of this dataset is included in this repo in archive at `./Data/behavior_free_swimming.tar.gz`.  
To load and view the content of this file from this repository you can use :

```python
from utils.data_and_models import extract_data
import h5py
extract_data("./Data/behavior_free_swimming.tar.gz")
file = h5py.File("./Data/behavior_free_swimming.h5", "r")
#h5tree_view(file)
```

Here is the structure of this file. For each temperature ([18°C, 22°C, 36°C, 30°C, 33°C]) it contains many swimming parameters, But importantly $\delta\theta_{f,n}$ (=`dtheta`) with $f$ the fish trajectory and $n$ the bout number, the reorientation angles of the fish at each bouts.

```
. /path/to/behaviour_free_swimming.h5
├── 📁behaviour
│   ├── 📁18
│   │   └── 🏷️temperature = `18°C`
│   │   ├── ...
│   │   ├── 🔢dtheta ⚙️(532, 641)float64
│   │   │   └── 🏷️unit = `degree`
│   │   ├── ...
│   ├── 📁 ...
```

⚠️ `dtheta` contain `NaNs` ! During the experiments, larvae regularly left the field of view of the camera, therefore the trajectories don't all last the same number of bouts. In order to save the trajectories in arrays, they were post-padded with `np.nan`s. In order to load the sequences of `dtheta` for all the trajectories and remove the `NaNs`, you can use :

```python
from utils.data_and_models import load_sequences

temp = 18 # for 18°C
dthetas = load_sequences(
    "./Data/behavior_free_swimming.h5",
    temp,
)
```

## Hidden Markov Models
Hidden Markov Models were computed from a custom implementation which can be found at [here](https://github.com/ZebrafishHMM2023/ZebrafishHMM2023.jl).  
All models infered from the data and used in this work are included in the present repo at `./Models/hmms_20240125.tar.gz` (for multi-fish experiements), and `./Models/longtrajectories_20240202.tar.gz` (for single-fish experiments).


## Tutorial

- LJP : [tutorial with holes](https://github.com/EmeEmu/IBIO-Banyuls2023-Python/blob/main/day4_HMMs.ipynb), [solution](https://github.com/EmeEmu/IBIO-Banyuls2023-Python/blob/main/corrections/day4_HMMs_correction.ipynb)
- ENS : ?????
