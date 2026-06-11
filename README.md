# Multi-area Decision Dynamics Across Human Cortex Shape Confidence

This repository contains the MATLAB code used to reproduce the analyses and figures reported in:

**Toso A., Arazi A., de la Rocha J., Tsetsos K., & Donner T.H.**
**Multi-area Decision Dynamics Across Human Cortex Shape Confidence**
bioRxiv (2025)

The repository contains code for:

* Behavioral analyses.
* Computational model fitting and model comparison.
* Model simulations of choice and confidence.
* Neural decoding analyses from source-reconstructed MEG data.
* Brain-map visualizations.
* Reproduction of main and supplementary figures reported in the manuscript.

---

# Requirements

The analyses were developed and tested using:
* MATLAB R2020b
The repository relies on several external toolboxes listed below.

---

# Repository Structure

```text
├── Figure_1/
├── Figure_2/
├── Figure_3/
├── Figure_4/
├── functions/
└── data/
```

## Figure folders

Each figure directory contains the MATLAB scripts used to reproduce the corresponding figures reported in the manuscript.

---

# Data Availability

The behavioral and processed MEG datasets required to reproduce the analyses can be downloaded from:

https://figshare.com/s/96510af31ee742b4498e

After downloading, place the data into the repository `data/` directory while preserving the original folder structure.

The full raw MEG dataset will be made available through the University of Hamburg Center for Sustainable Data Management.

The repository includes all processed data required to reproduce the published figures. Recomputing decoder outputs directly from source-level MEG activity requires access to the raw/source-reconstructed MEG data.

---

# Data Organization

```text
data/
├── behavior/
├── model/
├── best_model_fit/
├── decoder_output/
├── brain_map_decoder/
└── generalization_matrix/
```

### behavior/
Contains behavioral responses and task-related variables used throughout the analyses.

### model/

Contains outputs from fitting the candidate computational models to behavioral data, including trial-wise predictions and model likelihoods.

### best_model_fit/

Contains simulations and outputs associated with the best-fitting model architecture.

### decoder_output/

Contains processed neural decoder outputs used throughout the manuscript.

### brain_map_decoder/

Contains ROI-based and whole-cortex decoder outputs used for brain-map analyses.

### generalization_matrix/

Contains temporal and cross-condition generalization analyses used in Figures 2–3 and related supplementary figures.

---

# Computational Modeling

The repository evaluates multiple computational models of perceptual decision-making and confidence.

The model fitting framework is located in:

```text
functions/model_fit/
```

### Model fitting

The central fitting routine is:

```text
model_fit.m
```

This function estimates model parameters by maximizing the likelihood of the observed behavioral responses.

The fitting framework includes:

* Model-specific cost functions.
* Model-specific simulation functions.
* Parameter optimization routines.

### Model simulations

The function:

```text
model_simulation.m
```

uses fitted model parameters to generate simulated behavioral responses. These simulations are used to generate the model predictions shown throughout Figure 1 and associated supplementary figures.

### Decision-variable simulations

The function:

```text
leaky_integrator_simulations.m
```

simulates model decision variables. Unlike the fitting procedure, these simulations use a temporal resolution that is ten times finer than the original stimulus sampling rate. This higher-resolution simulation allows detailed characterization of accumulation dynamics throughout the decision process.

The resulting decision variables are used for the model-based analyses presented in:

* Figure 2
* Figure 3
* Supplementary Figure S5
* Supplementary Figure S6

---

# Neural Decoding

The neural decoding analyses are implemented in:

```text
functions/decoder/
```

These functions operate on source-reconstructed MEG activity and generate decoder outputs used throughout the manuscript.
The decoding framework supports both Unilateral and Bilateral ROI analyses.
The processed decoder outputs necessary to reproduce the published figures are included in the repository (under data/).

---

# External Dependencies

## Bayesian Adaptive Direct Search (BADS)

Model fitting was performed using Bayesian Adaptive Direct Search (BADS):

Acerbi, L., & Ma, W. J. (2017).

Repository:

https://github.com/acerbilab/bads

## VBA Toolbox

Protected exceedance probabilities used for model comparison (Supplementary Figure S10) were computed using the VBA toolbox:

Daunizeau, J., Adam, V., & Rigoux, L. (2014).

Repository:

https://github.com/MBB-team/VBA-toolbox

## MATLAB GIfTI Library

Brain surface visualizations shown in Figure 2 were generated using the MATLAB GIfTI Library:

https://www.nitrc.org/projects/gifti/

---

# Reproducing the Figures

1. Clone this repository.
2. Download the datasets from Figshare.
3. Place the downloaded data in the `data/` directory.
4. Install the required external toolboxes.
5. Add the `functions/` directory to the MATLAB path.
6. Run the scripts contained in the corresponding `Figure_*` directories.

The figure scripts are organized to reproduce the analyses reported in the manuscript.

---

# Citation

If you use this code or data, please cite:

Toso, A., Arazi, A., de la Rocha, J., Tsetsos, K., & Donner, T. H. (2025).

*Multi-area Decision Dynamics Across Human Cortex Shape Confidence.*

bioRxiv.
