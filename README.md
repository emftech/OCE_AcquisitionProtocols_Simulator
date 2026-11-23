# OCE_AcquisitionProtocols_Simulator

A MATLAB simulation framework for designing and comparing **acquisition protocols for wave-based Optical Coherence Elastography (OCE)**. It covers the conventional and **Multimeridian MB mode** and the novel **Synchronized BM-mode (SBM)** protocol proposed in the thesis. It answers one question: *for a given scan pattern and timing, how well can the phase-speed (stiffness) map be recovered?*

> **About SBM:** the Synchronized BM-mode protocol is a **novel acquisition protocol** introduced in this thesis. It has **not yet been published** in a peer-reviewed journal. Please cite the thesis if you use or build on it.

Developed in the undergraduate thesis of **Estéfano Machado Ferrer**, under the supervision of Prof. Fernando Zvietcovich, PhD, at the Biophotonics and Biomedical Optics Research Group (GIBIO), Pontificia Universidad Católica del Perú.

> **Thesis:** E. Machado Ferrer, *Development of an Acquisition Protocol for Non-Cylindrical Wave Propagation in Wave-Based OCE*, PUCP, 2026. PDF in [`docs/`](docs/).
>
> This repository is **simulation only**. The hardware implementation on a real OCT system (scan-pattern generator, LabVIEW integration and processing of experimental data) is in the companion repository **[OCE_Experimental_Protocols](https://github.com/emftech/OCE_Experimental_Protocols)**.

## What the simulator does

```
 1. Wave propagation     k-Wave simulation of a harmonic or quasi-harmonic source in a custom 2D medium
                         (homogeneous / cylindrical, or with stiff inclusions / non-cylindrical)
          │
 2. Acquisition vectors  exact galvo X/Y positions and trigger timing of the chosen protocol,
                         on the OCT A-line clock (200 kHz)
          │
 3. Sampling             the wavefield is sampled where and when the OCT beam actually is,
                         optionally with eye motion (saccades)
          │
 4. Recovery             MB: per-meridian M-B frames  ·  SBM: equivalent-time recovery with
                         per-pixel acquisition-phase compensation
          │
 5. Speed estimation     1D PG along meridians, 2D phase gradient (PG), TV-regularized 2D PG (PG-TV)
          │
 6. Metrics              accuracy error Ea, precision error Ep and CNR in background and inclusion ROIs
```

### Protocols

| Protocol | Scan pattern | How the wave is captured |
|---|---|---|
| **Multimeridian MB mode** | M-B scans along N meridians through the excitation source | each meridian is an M-mode B-scan; the meridians are interpolated into an en-face map |
| **Synchronized BM mode (SBM)** (novel, unpublished) | ordinary raster BM scan, repeated over N_F frames | the scan is locked to the excitation so each frame samples every pixel at a slightly later phase; stacking frames reconstructs the 2D wavefield in equivalent time |

### Simulation cases

| Case | Medium | Excitation | Live scripts |
|---|---|---|---|
| Ideal | non-cylindrical (inclusions) | harmonic | `Sampling_Acq_Sim_MB_v4_Ideal.mlx`, `Sampling_Acq_Sim_SBM_v4_Ideal.mlx` |
| Cylindrical | homogeneous | harmonic (SBM), quasi-harmonic (MB) | `Sampling_Acq_Sim_SBM_v4_Cylindrical_Ideal.mlx`, `Sampling_Acq_Sim_MB_v4_Ideal_Cylindrical_QH.mlx` |
| Saccades | as above, plus simulated eye motion | harmonic / quasi-harmonic | `Sampling_Acq_Sim_MB_v4_Saccades.mlx`, `Sampling_Acq_Sim_SBM_v4_Saccades.mlx`, `Sampling_Acq_Sim_MB_v4_Saccades_Cylindrical_QH.mlx` |

`SBM_Aux.mlx` holds auxiliary SBM timing calculations.

## Repository layout

```
├── Sampling_Acq_Sim_*_v4_*.mlx   live scripts, one per protocol and case (documented step by step)
├── SBM_Aux.mlx
├── Codes/                        framework functions
│   ├── Wave_Prop_Simulation_v3.m, CustomMedium_{Harmonic,QuasiHarmonic}_kWave_Acoustic.m   wave propagation
│   ├── AcqVectors_MultimeridianMB.m, AcqVectors_SBM*.m                                     acquisition vectors
│   ├── SamplingSim_fromAcqVecs*.m, Sampling_fromAcqVecs_Generalized.m                      sampling
│   ├── Recovery_fromSampledFrames_SBM*.m                                                   SBM recovery
│   ├── Prepare_{MB,SBM}_2DPhaseGradient*.m, PhaseGradient*.m, TVDenoise_ChambollePock2D.m speed estimators
│   ├── Compute_{AccuracyError,PrecisionError,CNR}_Masks*.m, Get_*Mask*.m                   metrics and ROIs
│   └── Add_Saccades_*.m                                                                    eye-motion models
├── Batch_Simulations/            batch study matched to the experiments (18 x 18 mm, c = 10.54 m/s)
│   ├── AcqSim_Config.m           all parameters and run lists
│   ├── Run_AcqSim_Batch.m        runs every MB and SBM configuration
│   ├── Analyze_AcqSim_Batch.m    summary tables, metrics vs. meridians / B-scans, dispersion
│   └── display_helpers/
├── figures/                      representative thesis results
└── docs/                         thesis manuscript (PDF)
```

## Requirements

- MATLAB R2026a (developed and tested).
- Signal Processing, Image Processing, Statistics and Machine Learning, and Curve Fitting toolboxes. Parallel Computing Toolbox is optional.
- [k-Wave](http://www.k-wave.org) (LGPL-3.0, not included): download it and place it in `Codes/k-Wave/`.

## Quick start

1. Clone the repository and put k-Wave in `Codes/k-Wave/`.
2. Open a live script, for example `Sampling_Acq_Sim_SBM_v4_Ideal.mlx`, and run it section by section. The first section finds the repository folder, adds `Codes/` and its subfolders (including k-Wave) to the path, and results are saved under `Results/` (ignored by git). The current folder is never changed. If you run a script outside the Live Editor, set the current folder to the repository root first.
3. Scripts with `sim_opts.useGUI = true` open a GUI to define the medium (background speed, inclusions). The chosen medium is saved in your MATLAB preferences, and headless runs (`useGUI = false`) reuse it, or the default medium if none was saved.
4. For the batch study, run `Batch_Simulations/Run_AcqSim_Batch.m` and then `Batch_Simulations/Analyze_AcqSim_Batch.m`.

## Representative results

| | |
|---|---|
| ![Acquisition phase maps](figures/Acq_Phase_Maps.png) | ![Non-cylindrical wave propagation](figures/NonCylindrical_Elastic_WaveProp_Snapshots.png) |
| Acquisition phase of every pixel for the scan protocol | Simulated non-cylindrical wave propagation (ground truth) |
| ![Raw sampled frames](figures/Sampled_Frames_NonCylindrical_Raw_Snapshots.png) | ![Compensated sampled frames](figures/Sampled_Frames_NonCylindrical_Compensated_Snapshots.png) |
| Sampled frames before phase compensation | Sampled frames after phase compensation |
| ![Recovered frames](figures/Recovered_Frames_Snapshots_NonCylindrical.png) | ![Phase comparison](figures/Phase_Comparison_NonCylindrical.png) |
| Recovered wavefield with SBM | Recovered vs. ideal phase |
| ![Phase difference](figures/Phase_Diff_Maps_NonCylindrical.png) | ![Recovery without synchronization](figures/Recovery_NoSync.png) |
| Phase-difference maps | Recovery without synchronization, for comparison |

More figures, the MB vs. SBM comparison and the full metrics are in the thesis.

## Citation

If you use this code, please cite the thesis (and the corresponding journal articles once published):

```bibtex
@thesis{MachadoFerrer2026SBM,
  author = {Machado Ferrer, Estéfano},
  title  = {Development of an Acquisition Protocol for Non-Cylindrical Wave Propagation in Wave-Based OCE},
  school = {Pontificia Universidad Católica del Perú},
  year   = {2026},
  type   = {Undergraduate thesis}
}
```

## License

Code: MIT (see [LICENSE](LICENSE)). The thesis manuscript and figures are © the author, all rights reserved. k-Wave is distributed separately under LGPL-3.0.
