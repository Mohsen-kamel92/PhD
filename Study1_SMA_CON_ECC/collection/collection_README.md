# Data collection

The scripts used to run the experiment and to build the torque traces the
participants followed. These run at collection time and are not part of the
analysis pipeline, which begins with the recordings they produce.

---

## Files

| File | Runs in | Does |
|---|---|---|
| `collection_protocol.s2s` | Spike2 | The experiment. Draws the target trace in an XY view, provides live torque feedback against it, and triggers the dynamometer rotation and the EMG system |
| `fs_sequencer_adaptable.pls` | Spike2 | Output sequencer called by the protocol script. Drives the dynamometer, sets the DAC1 trigger the analysis uses for synchronization, and controls the digital lines |
| `Sampling_Config.s2cx` | Spike2 | Sampling configuration: channels, rates and scaling |
| `SaveMAT.s2s` | Spike2 | Exports a selected window of a Spike2 recording to `.mat`, which is what Stage 00 of the analysis reads |
| `Flip_Contraction.m` | MATLAB | Builds the target traces from a participant's maximal concentric contraction |
| `dlg_SetTrace.m` | MATLAB | Dialog used by `Flip_Contraction` to set the ramp duration, the scaling factors and the contraction mode |

---

## How the target traces were made

Each participant first performed a maximal concentric dorsiflexion, recorded in
Spike2 and exported with `SaveMAT.s2s`. `Flip_Contraction.m` then took that
recording and produced the traces they followed for the rest of the session:

1. Torque and joint angle are low-pass filtered, and the rotation phase is
   located from the points where the angle stops changing.
2. The passive torque-angle relation is fitted and subtracted, leaving the
   active torque produced during the rotation.
3. The active component is scaled to each target intensity and the passive
   component added back, so the trace the participant sees is a real torque
   they could produce rather than a fraction of a total.
4. The concentric trace is reversed in time to give the eccentric trace, so both
   modes follow the same torque at the same joint angle.
5. Resting torque at the short and the long position is appended at either end,
   giving a trace that starts and ends where the limb rests.

The result is written as a two-column text file, time against torque, which the
Spike2 protocol reads into its XY view.

### The file names

Traces are saved as `ConMatch_P1_75.txt`, `EccMatch_P1_90.txt` and so on. The
number is the scaling factor as a percentage. The `P1` or `P2` reflects how many
scaling factors were entered in one run and does not identify the participant,
which is why the analysis matches ramp files on contraction mode and intensity
alone.

---

## Requirements

Spike2 version 8 or later, with a CED Power1401 for the sequencer.

MATLAB with the Signal Processing Toolbox, for `findchangepts` and `filtfilt`.
`Flip_Contraction.m` also uses `Wcu`, `Auto_RestTorque` and `Normalize_Torque`
from the `functions` folder one level up, which it adds to the path itself.

---

## Authorship

`collection_protocol.s2s` was written by Paolo Tecchio and Tobias Weingarten and
adapted for this study. `SaveMAT.s2s`, `Flip_Contraction.m` and `dlg_SetTrace.m`
were written in the Human Movement Science lab at Ruhr University Bochum. All
are included here with the authors' agreement and are released under the licence
in the folder above.
