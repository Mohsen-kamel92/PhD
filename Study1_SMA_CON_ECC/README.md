# Spatial muscle activity during concentric and eccentric contractions

MATLAB analysis pipeline for a study comparing the spatial distribution of
tibialis anterior activity between torque-matched concentric and eccentric
dorsiflexion contractions, recorded with high-density surface EMG.

The code takes raw dynamometer and HDsEMG recordings through synchronization,
passive torque correction, EMG normalization, spatial analysis and group
statistics, and reproduces every figure and reported value in the manuscript.

---

## What is here

```
run_pipeline.m      driver; runs the stages in order from one root selection
scripts/            the analysis stages
functions/          helper functions and their licences
data/               processed group dataset
```

The raw recordings are too large for this repository. They are archived
separately; see **Data** below.

---

## Quick start

Reproducing the reported statistics and figures needs only the processed group
dataset, which ships here.

1. Clone the repository.
2. Install [spm1d for MATLAB](https://spm1d.org).
3. Open `run_pipeline.m` and set `PIPE.spm1dPath` to where spm1d lives.
4. Edit `PIPE.stages` to keep only the Stage 06B onward entries.
5. Press Run and select the `data` folder when prompted.

Every stage can also be run on its own. Open it, press Run, and choose the same
folder when the dialog appears.

---

## The pipeline

Each stage reads what the previous one wrote. Stages 00 to 05B walk the raw
participant folders; Stage 06A assembles the group dataset; everything after it
reads only that file.

| Stage | Script | Does |
|---|---|---|
| 00 | `Stage00_sync.m` | Aligns the dynamometer and HDsEMG recordings on a common time base |
| 01 | `Stage01_Matching_Validity.m` | Locates repetitions, finds the constant velocity window, checks torque against the prescribed ramp |
| 02A | `Stage02A_active_torque.m` | Fits and subtracts the passive torque, converts the windows to sample indices |
| 02B | `Stage02B_active_torqueMVC.m` | The same for the maximal voluntary contractions |
| 03 | `Stage03_bipolar_SignalNormalization.m` | Bipolar derivation, band-pass filtering, RMS envelope, normalization to %MVA |
| 04 | `Stage04_WeightedCentroid_processing.m` | Amplitude-weighted centroid in millimetres on the grid |
| 05A | `Stage05A_setLevel_aggregation_acrossReps.m` | Averages repetitions into set-level curves, with early and late bins |
| 05B | `Stage05B_ConditionLevel_aggregation_acrossSets.m` | Averages sets into condition-level curves |
| 06A | `Stage06A_buildGroup_dataset.m` | Assembles the group dataset the statistics read |
| 06B | `Stage06B_scalar_RM_ANOVA.m` | Two-way repeated measures ANOVA on window means, with paired follow-ups |
| 06C1 | `Stage06C1_SPM_TimeSeries_inference.m` | Time-series inference across the analysis window |
| 06C2 | `Stage06C2_Early_vs_Late.m` | Early against late repetitions |
| 06C3 | `Stage06C3_heat_maps.m` | Group-mean spatial maps |
| 06C4 | `Stage06C4_entropy_spm.m` | Shannon spatial entropy |

Every stage writes a CSV summary and a folder of quality control figures next to
its outputs, so any step can be checked without rerunning the analysis.

### Analysis decisions worth knowing

Each script's header explains its own choices in full. Four apply throughout:

**The analysis window.** Only the middle 20 to 80% of the constant velocity
phase is analysed. This is where the two contraction modes align most closely in
joint angle, and where between-participant variability is lowest.

**Eccentric time reversal.** The two modes traverse the range of motion in
opposite directions, so eccentric series are reversed in time before comparison.
This aligns them by joint angle rather than by elapsed time. Scalar means are
unaffected.

**Normalization.** Amplitude is expressed relative to a single scalar per
participant, the 95th percentile of the grid-averaged MVC RMS. Normalizing each
channel to its own MVC value would re-weight the spatial map and displace the
centroid; one scalar leaves the spatial pattern intact.

**Torque.** Torque is a manipulation check rather than an outcome. It is tested
so that the quality of matching can be judged, not so that a difference in it
can be interpreted.

---

## Data

**In this repository.** `data/groupData_stage6A.mat` holds every outcome as
participant by condition by time arrays, which is what Stages 06B onward read.

**Archived separately.** The raw SPIKE and `.sig` recordings, and the
intermediate files from Stages 00 to 05B, are too large to host here. See the
archive linked in the manuscript's data availability statement.

### Expected layout for the raw stages

```
<ROOT>/<participant>/<CON_75|CON_90|ECC_75|ECC_90>/<Set_N>/   dynamometer .mat
<ROOT>/<participant>/zEMG/<condition>/<Set_N>/                .sig files
<ROOT>/<participant>/MVC/                                     MVC recording
<ROOT>/<participant>/Ramps/                                   prescribed ramps
```

Every `Set_*` folder present is treated as retained data. Attempts that failed
the live torque matching criterion were deleted during collection rather than
kept, so no exclusion rule is applied in the code.

---

## Requirements

MATLAB R2023a or later, for `xregion` in the quality control figures. Earlier
releases run the analysis but will error on some plots.

Signal Processing Toolbox, for `butter` and `filtfilt`. Statistics and Machine
Learning Toolbox, for `fitrm`, `ranova` and `ttest`.

[spm1d for MATLAB](https://spm1d.org), for the time-series inference in Stages
06C1, 06C2 and 06C4. Everything else runs without it.

### Functions in this repository

| Function | Used for | Author |
|---|---|---|
| `Wfilt`, `Wvel` | Zero-phase filtering and velocity | B. J. Raiteri |
| `Wcu` | Winter cutoff correction for dual passing | W. van den Hoorn |
| `Auto_RestTorque`, `Normalize_Torque` | Passive torque fit and subtraction | Human Movement Science, RUB |
| `distinguishable_colors` | Participant colours | Tim Holy |
| `swtest` | Shapiro-Wilk test | Ahmed Ben Saïda |

`distinguishable_colors` and `swtest` are redistributed under their BSD 2-clause
licences, included in `functions/` as `LICENSE_distinguishable_colors.txt` and
`LICENSE_swtest_BenSaida.txt`.

### Not included

`MAECS_read`, which reads the raw `.sig` recordings, was developed at LISiN,
Politecnico di Torino, and is available from the system developers. Only Stage
00 requires it.

---

## Filtering

Torque and joint angle are low-pass filtered with a second-order Butterworth
applied in both directions, giving a fourth-order zero-phase response, with the
cutoff corrected for dual passing (Winter, 2009). EMG is band-pass filtered
between 20 and 450 Hz with a fourth-order Butterworth applied in both
directions. The Winter correction does not apply to band-pass filters and is not
used there.

---

## Citation

If you use this code, please cite the associated paper. Details will be added
here on publication.

## Licence

See `LICENSE`. Third-party functions keep their own licences, included in
`functions/`.

## Contact

Abdelmohsen Eldhma, Department of Human Movement Science, Ruhr University
Bochum.
Abdelmohsen.Eldhma@ruhr-uni-bochum.de
