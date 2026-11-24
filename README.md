# Multi-area Decision Dynamics Across Human Cortex Shape Confidence

MATLAB code to reproduce the figures of:

**Toso A., Arazi A., de la Rocha J., Tsetsos K., & Donner T.H.**
**Multi-area Decision Dynamics Across Human Cortex Shape Confidence**

Every script reproduces the panels in its name (`F` = main figure, `FS` = supplementary figure), e.g. `F1B_FS3A_B.m` = Figure 1B and Figure S3A-B. Scripts are run from anywhere: paths are relative to the repository root (the folder above the script).

---

## Requirements

* MATLAB R2020b, Statistics and Machine Learning Toolbox
* Included in `functions/`: VBA toolbox (Figure S2A), BADS (model fitting and recovery), gifti toolbox, brewermap, shadedErrorBar, permutationTest, fdr_bh

## Data

The behavioral data, model fits and decoder outputs can be downloaded from (uLink will be available upon publication) and go into `data/`:

```text
data/
├── behavior/data.mat                         All_Data: subjects x drug conditions
├── model/
│   ├── results/                              All_data_model_<m>_sbj_<s>.mat, Parameters_model_<m>_sbj_<s>.mat
│   └── recovery/                             rec_g09_d<dataset>_f<model>.mat (Model_recovery/)
├── brain_map_decoder/                        to reproduce Figure 2A brainmaps
├── generalization_matrix/                    temporal cross-decoding matrices of PMd and M1
└── decoder_output/
    ├── ROI_bilateral/                        bilateral decoders of the selected ROIs (Figures 4, 5)
    ├── PMd_M1_hemi/                          single-hemisphere decoders of 6a, 6d, 4 (Figure 2)
    ├── PMd_M1_hemi_behavior/                 single-hemisphere decoders of 6A6D, 4, with confidence and RT (Figure 3)
    ├── all_areas/                            choice and confidence decoders on all areas (Figure S4E-G)
    ├── all_VX/                               decoder weights of the PMd/M1 vertices (Figure 2E)
    ├── bilateral/stimulus/                   PMd/M1 decoder output, 100 Hz (Figure S7F)
    └── v1/                                   contrast decoding from V1 (Figure S7B-C)
```

---

## Scripts

### Figure_1 (behavior and model)

| Script | Panels |
|---|---|
| `F1B_FS3A_B.m` | 1B confidence vs evidence; S3A confidence-RT correlation; S3B RT vs evidence |
| `F1C_FS1E_F.m` | 1C psychometric function; S1E psychophysical kernel; S1F primacy index |
| `F1D_E.m` | 1D confidence vs evidence for correct and error trials; 1E their correlations |
| `F1F_G.m` | 1F confidence kernels; 1G kernel asymmetry |
| `FS1C_D.m` | S1C-D correlation of left and right contrast fluctuations |
| `FS2_A.m` | S2A Bayesian model comparison (exceedance probability) |
| `FS2_B.m` | S2B model recovery |
| `FS2_C.m` | S2C mean squared error of the confidence kernels |
| `FS3C_D.m` | S3C RT vs evidence for correct and error trials; S3D their correlations |
| `FS3E_F.m` | S3E-F kernels split by RT |

### Figure_2 (neural DVs in PMd)

| Script | Panels |
|---|---|
| `F2A.m` | 2A whole-cortex map of choice decoding |
| `F2B_FS4A.m` | 2B (PMd), S4A (M1) bilateral choice decoding |
| `F2C_FS4B.m` | 2C (PMd), S4B (M1) temporal cross-decoding matrix |
| `F2D_FS4C.m` | 2D (PMd), S4C (M1) variance of the decoder output |
| `F2E.m` | 2E decoder weights per vertex |
| `F2F.m` | 2F example trials |
| `F2G.m` | 2G winning and losing neural DVs |
| `F2H_FS6D.m` | 2H across-trial correlation of neural (PMd) and model DVs; S6D same for M1 |
| `F2I_FS6F.m` | 2I (PMd), S6F (M1) within-trial cross-correlograms |
| `FS4E_G.m` | S4E-G choice and confidence decoders on all areas |
| `FS6A_B.m` | S6A-B correlation of left and right sample means |
| `FS6C_FS7A_D_E_G_H_I.m` | S6C, S7A, D, E, G, H, I model parameters |
| `FS6E.m` | S6E within-trial cross-correlograms of model DVs |
| `FS7B_C.m` | S7B-C contrast decoding from V1 |
| `FS7F.m` | S7F leak vs neural DV slope |

### Figure_3 (neural DVs, confidence and RT)

| Script | Panels |
|---|---|
| `F3B_C_FS9B_C_F_FS10.m` | 3B-C, S9B-C, S9F, S10 regressions of confidence and RT on winning and losing neural DVs (PMd, M1) |
| `F3D_FS8B.m` | 3D best-fitting model, S8B model with rectification: neural vs model DVs |
| `FS8A.m` | S8A kernel asymmetry of the model with rectification |
| `FS9D_E.m` | S9D-E model DVs and decision time, collapsing vs fixed bound |

### Figure_4 (association cortex)

| Script | Panels |
|---|---|
| `ROI_selection_TableS1.m` | selection of the ROIs (three criteria), Table S1 |
| `FS5.m` | S5A-B characterization of the selected ROIs |
| `F4_A_D_F5_A_B.m` | 4A-D and 5A-B decoding maps, normalized DVs, ramp onset, correlation matrices |

### Figure_5 (residual DVs)

| Script | Panels |
|---|---|
| `F5_D_FS11_A_D_F.m` | 5D, S11A-D, S11F residual DVs of winning and losing areas |

---

## Functions

| Function | Use |
|---|---|
| `leaky_integrator_simulations_fine.m` | accumulator model simulated in 10 ms steps with the fitted parameters |
| `model_simulation.m` | accumulator model with the 10-sample time step of the fits (Figure 1D-E) |
| `absorb_bound_dvs.m` | absorbing bound on simulated accumulators |
| `cluster_perm_1samp.m` | cluster-based permutation test against zero |
| `decoder/` | decoders run on the source-reconstructed MEG data (not included) |
| `model_fit/` | model fitting (`model_fit.m`) and the 19 models |

`Model_recovery/` generates the synthetic datasets of the model recovery and fits them; settings in `recovery_setup.m`.

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
---

# Citation

If you use this code or data, please cite:

Toso, A., Arazi, A., de la Rocha, J., Tsetsos, K., & Donner, T. H. (2025).

*Multi-area Decision Dynamics Across Human Cortex Shape Confidence.*

bioRxiv.