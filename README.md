# PACE

**Parameter ensemble, Adaptive thresholds, Clinical optimization, and MIC Extrapolation**

Companion code and data for:

> Dewaele K, André E, Tzfadia O, Rigouts L, de Jong B, Lagrou K. PACE: three- to seven-fold faster antifungal MICs from live-cell imaging of standard broth microdilution plates in *Aspergillus fumigatus* and *Candida* species. *Manuscript submitted for publication*, 2026.

---

## Overview

PACE predicts minimum inhibitory concentrations (MICs) early from time-lapse microscopy of unmodified EUCAST broth
microdilution plates. Per well, it derives growth features from background-corrected absorption (BCA), an image-based
biomass measure, and optimises proportional growth/no-growth thresholds (relative to each isolate's growth control)
against the reference MIC at every timepoint. A dynamic ensemble then switches between the best-performing
parameter-threshold-shift combinations ("triplets") over the incubation time.

Results reported in the paper (first timepoint at which the ensemble reaches at least 90 % essential agreement, EA):
- *A. fumigatus*: isavuconazole 7 h, amphotericin B 9 h, voriconazole 9 h, posaconazole 12 h (reference read-out at 48 h).
- *Candida* spp.: fluconazole 8.5 h (reference read-out at 24 h).
- For isavuconazole, amphotericin B, voriconazole and fluconazole these times changed by at most 0.5 h in
  leave-one-isolate-out cross-validation (Supplementary Table S2).

---

## Repository structure

```
.
├── set_thresh.py                     # PACE pipeline: phase-1 threshold grid and phase-2 single-triplet runs
├── commands_phase1.txt               # phase 1: EA- and CA-optimised grids, all drugs
├── commands_phase2_aspergillus.txt   # phase 2: selected triplets, A. fumigatus
├── commands_phase2_candida.txt       # phase 2: selected triplets, Candida
├── commands_phase1_sesa_candida.txt  # method comparison (SESA feature), Candida
├── commands_phase2_sesa_candida.txt
├── cleanup_temporary_files.py        # removes intermediate .pkl files (see the warning under step 4)
├── requirements.txt
├── data/
│   ├── aspergillus/                  # 30 A. fumigatus plates (UniExplorer Excel exports)
│   └── candida/                      # 20 Candida plates
├── job_templ/                        # SLURM job templates (placeholders for account and e-mail)
├── pace_tools/
│   ├── data_extraction.py            # feature extraction from the Excel files
│   ├── dictionaries.py               # reference MICs, EUCAST v12 breakpoints, species, plate layout
│   ├── evaluate.py                   # EA, CA, VME, ME and threshold selection
│   ├── predict.py                    # growth/no-growth pattern detection and MIC prediction
│   └── plot.py, setup.py, phot_dict.py, examine_pickles.py
├── analysis/
│   ├── select_triplets.py            # phase 2: triplets from the phase-1 runs -> phase-2 commands
│   ├── run_ensembles.py              # phase 3: ensembles, times to result, Figure 5 values, Text S2, Tables S1 and S3
│   ├── run_loocv.py                  # leave-one-isolate-out cross-validation (Table S2)
│   └── pace_analysis/                # shared code (ensemble rule, triplet rule, EUCAST ATU rule, scoring)
└── supplementary_tables/             # Supplementary Tables S1-S3 of the paper (legend on the first sheet)
    ├── Table_S1_reference_vs_PACE_MICs.xlsx
    ├── Table_S2_crossvalidated_time_to_result.xlsx
    └── Table_S3_azole_class_calls.xlsx
```

---

## Data

Each file in `data/` is a UniExplorer (oCelloScope) Excel export of one plate, named by its study code. Sheet
"Job configuration" holds the acquisition settings (including the imaging interval); sheet "Feature values" holds,
per well and timepoint, the image-analysis features BCA, TA and SESAfungi (exported as log10 values) and their
normalised versions. PACE evaluates every plate at 30-minute timepoints.

*Aspergillus fumigatus* (30 isolates, `data/aspergillus/`):
- QC strain ATCC 204305 (`ASFU_ATCC_204305_S`)
- 13 azole-resistant isolates with TR34/L98H (`ASFU_R_TR34_*`)
- 8 isolates with TR46/Y121F/T289A (`ASFU_R_TR46_*`)
- 1 azole-resistant isolate without a cyp51A mutation (`ASFU_R_unknown_22`)
- 7 wild-type susceptible isolates (`ASFU_S_*`)

*Candida* (20 isolates, `data/candida/`):
- *C. albicans* ATCC 90029; *C. auris*: CDC AR Bank strains (clades I, II, V), 5 clinical isolates and 2 clinical
  isolates of clades III and IV; *C. dubliniensis*, *C. glabrata*, *C. krusei*, *C. parapsilosis* and *C. tropicalis*
  (ATCC and clinical isolates).

Reference MICs are in `pace_tools/dictionaries.py` (well codes; `Z0` = no valid reference MIC). Drug rows:
`antimycotics_dict_A` and `antimycotics_dict_C`.

---

## Installation

Python 3.11 with the package versions in `requirements.txt` (NumPy must be < 2):

```bash
pip install -r requirements.txt
```

`set_thresh.py` imports `pace_tools/` from its own directory; run all commands from the repository root.

---

## Reproducing the paper

All commands write their session folders to the working directory (or to `-o`). The analysis scripts find their inputs
below the folders you pass them.

**0. Data dictionaries.** The first run of each genus extracts the features from `data/` and saves
`genus_A_d_window_6.pkl` (or `genus_C_...`) in the working directory. Move it to `output/` so that later runs, which
point `-d` there, load it instead of extracting the data again:

```bash
mkdir -p output
# run the first command of commands_phase2_aspergillus.txt (or of commands_phase1.txt), then:
mv genus_A_d_window_6.pkl output/
```

**1. Phase 1, threshold grid** (`commands_phase1.txt`, 27 runs). For every timepoint, shift (-2 to +2) and parameter,
all proportional thresholds 0.00-0.99 are scored, with two criteria: EA (`-n EA_EA...`) and CA >= 0.9 (`-n CA_0.9...`).
These are the long runs: about 10 min (*Candida*, 0-24 h) to 25 min (*A. fumigatus*, 0-48 h) each on a recent
laptop; `job_templ/` has SLURM templates.

**2. Phase 2, stable triplets.** Select the triplets from the phase-1 runs and run each one over all timepoints:

```bash
python3 analysis/select_triplets.py PHASE1_DIR --commands   # prints the lines of commands_phase2_*.txt
bash commands_phase2_aspergillus.txt                        # 43 runs, about 1 min and 2.5 GB memory each
bash commands_phase2_candida.txt                            # 60 runs
```

**3. Phase 3, ensembles and per-isolate results:**

```bash
python3 analysis/run_ensembles.py --phase2 PHASE2_DIR --out results/
```

| Output (`results/`) | Paper |
|---|---|
| `figure5_schedules.txt` | Supplementary Text S2 (ensemble per drug) |
| `figure5_metrics.csv`, `figure5.png` | Figure 5 (EA, CA, VME and ME per timepoint) |
| `figure5_time_to_result.csv` | times to result, with EA and categorical errors at that time |
| `figure5_predicted_mics.csv` | Table S1: reference and PACE MIC, category and EA per isolate, here at every timepoint (Table S1 shows the time to result, 24 h and 48 h); also sheet Predicted MICs of Table S3 |
| `figure5_azole_class_calls.csv` | Table S3, sheet Calls: predicted azole categories, integrated call and check per isolate and timepoint |
| `figure5_azole_class_summary.csv` | Table S3, sheet Summary: counts per timepoint for the integrated call |

The tables themselves are in `supplementary_tables/`. The comparison rows on sheet Summary of Table S3 (one-R rule,
each azole alone) follow from the per-azole categories in `figure5_azole_class_calls.csv` and
`figure5_predicted_mics.csv`.

`--which s10` gives the CA-optimised ensembles of Supplementary Figure S10 instead.

Figure 5 uses the EA-optimised triplets for *A. fumigatus* and the EA- and CA-optimised triplets for *Candida*. The
ensemble keeps the active triplet as long as it meets the EA, CA, VME and ME targets; when no triplet meets them, a
switch requires an EA gain of more than 2 percentage points. Triplets that are not active for at least 2 consecutive
hours in a first pass are dropped, and the selection is repeated. Categories in
the paper follow the EUCAST v12 notes for the areas of technical uncertainty (isavuconazole 2 mg/L, posaconazole
0.25 mg/L; `analysis/pace_analysis/eucast_atu.py`), and VME and ME rates use resistant and susceptible isolates as
denominators. During threshold optimisation, `set_thresh.py` counts the ATU values as R.

**4. Leave-one-isolate-out cross-validation** (Supplementary Table S2):

```bash
python3 analysis/run_loocv.py --phase1 PHASE1_DIR --phase2 PHASE2_DIR --out results/
```

For every isolate, phase 1, phase 2 and the ensemble are repeated without it, and the left-out isolate is predicted
with the resulting ensemble. No PACE re-runs are needed: a phase-1 run stores every isolate's prediction for every
threshold, and with proportional thresholds an isolate's prediction does not depend on the other isolates. The script
reads the per-timepoint files `performance_dict_<row>_<time>_EA.pkl` of the EA phase-1 runs; do not run
`cleanup_temporary_files.py` on those folders before this step.

`table_s2_loocv.csv` holds the values of Table S2 (`supplementary_tables/Table_S2_crossvalidated_time_to_result.xlsx`).
`loocv_by_time.csv` gives the EA per timepoint, in-sample and cross-validated. `loocv_checks.txt` reports, per drug,
whether the triplets re-derived from phase 1 equal the phase-2 runs and whether the in-sample schedule equals the
Figure 5 ensemble.

---

## Key arguments of `set_thresh.py`

| Flag | Description |
|---|---|
| `-g` | Genus: `A` (Aspergillus) or `C` (Candida) |
| `-a` | Antimycotic row label (e.g. `A` = amphotericin B, `B` = voriconazole; see `dictionaries.py`) |
| `-i` | Input directory containing the Excel files |
| `-d` | Saved data dictionary (see step 0) |
| `-t` | Timepoint range as `(start, stop, step)` in seconds |
| `-s` | Threshold(s): range `[min, max, step]` or a single value |
| `-p` | Features to evaluate (list) |
| `-c` | Optimisation criterion, e.g. `[("EA", "EA")]` or `[("CA", 0.9)]` |
| `-m` | Shift (bias correction, in wells): an integer, or `[2]` for the range -2 to +2 |
| `-n` | Run label (used in output file names) |
| `-l` | Log level (`info`, `debug`, `warning`) |
| `-o` | Output directory (defaults to the current working directory) |

Each run creates a time-stamped session folder with a performance table (`performance_<label>_<time stamp>.pkl`;
inspect it with `python pace_tools/examine_pickles.py`), plots and a log file.

---

## Citation

If you use this code or data, please cite the paper above (see also `CITATION.cff`).

## License

MIT License. See [LICENSE](LICENSE).
