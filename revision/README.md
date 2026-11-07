# Revision code: Neubrandt, Lenkey & Vervaeke, Nature Communications (NCOMMS-25-18616)

MATLAB code added during revision (MATLAB R2025b with the Statistics and Machine Learning Toolbox and the Signal Processing Toolbox).
The preprocessing pipeline that creates the session files (sData) is in the repository root (`pipeline.m`, `+codes/`).
Several scripts are named after the figure numbers used during revision: "Figure3…" scripts make panels of final **Figure 4**, "Figure4…" scripts make final **Figure 5**.

## glm_figure2 — final Figure 2 (VIP encoding model)

| Script | Purpose |
|---|---|
| `vipGLM_fitCell.m` | Per-cell ridge-regression encoding model (GLM): cross-validated fit and unique ΔR² per feature |
| `vipGLM_session.m` | Fits all cells of one session |
| `vipGLM_cohort.m` | Main entry point: runs the GLM on all sessions and saves the cohort results (`cohortGLM_*.mat`) |
| `vipGLM_paperFigure.m` | Builds Figure 2a–g from the cohort results and the session files |
| `vipGLM_paperFigure_dots.m` | Final version of Figure 2 (adds per-session dots to panel b) |
| `vipGLM_drawExampleTrace.m`, `vipGLM_exampleCellTrace.m` | Example cell, observed vs held-out prediction (Figure 2c) |
| `vipGLM_positionTuning.m`, `vipGLM_featureTrialTrace.m` | Tuning and trial-by-trial activity per feature group (Figure 2e, f) |
| `vipGLM_recovery.m` | Validation of the GLM on surrogate data with known ground truth |
| `vipGLM_clusterValidation.m`, `vipGLM_clusterVsGLM.m`, `vipGLM_pcaOverlay.m` | Number of k-means clusters and comparison of clusters with GLM groups (review) |
| `vipGLM_primarySecondary.m`, `vipGLM_featureSwitchMap.m`, `vipGLM_novelOnsetDiag.m`, `vipGLM_inspect.m` | Additional analyses and an interactive single-cell inspector |

## figures_4_5_S7 — final Figures 4, 5 and Supplementary Figure 7 (statistics and plots)

| Script | Final panel |
|---|---|
| `Figure4a.m` | Figure 4a (place field emergence, LME) |
| `Figure3b.m` | Figure 4b (PV correlation, Novel 1 window) |
| `Figure3ef.m`, `NewFig4f.m` | Figure 4e, f, h, i (place field drift; per-mouse comparison) |
| `Figure3lm.m` | Figure 4l and Supplementary Figure 7a–c (backward shift) |
| `Figure3o.m` | Figure 4m (backward shift vs distance to reward) |
| `Figure3n.m` | Supplementary Figure 7d (place field density vs distance to reward) |
| `Figure3p.m` | Supplementary Figure 7e–h (place field width vs running speed) |
| `Figure4ef.m` | Figure 5e, f (per-block speed and lick statistics) |

## figure1_and_tools

| Script | Purpose |
|---|---|
| `figure1f.m` | Figure 1f, g (activity of all VIP neurons across the familiar–novel–familiar session) |
| `pupilTrackerRaw.m` | Pupil diameter and position from the eye-camera video (Figure 1h, Supplementary Figure 2) |
| `sessionInspectorGUI.m`, `vipCellInspector.m` | Interactive inspection of sessions and single VIP cells |

## Inputs

- Session files (sData) of the VIP imaging experiments; the default folder in the scripts is `~/Desktop/VIP-GCaMP-sData`.
- Cohort GLM results (`cohortGLM_*.mat`), produced by `vipGLM_cohort.m`.
- Per-session and per-place-field tables for Figures 4, 5 and Supplementary Figures 6–7 (`Fig3*.mat`, `Fig4Data.mat`, …).
