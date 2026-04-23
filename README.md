# PACE

**Parameter ensemble, Adaptive thresholds, Clinical optimization, and MIC Extrapolation**

Companion code and data for:

> Dewaele K *et al.* Same-day MIC results in drug-resistant *Aspergillus fumigatus* and *Candida auris* using live-cell imaging of a broth microdilution assay. *Manuscript submitted for publication*, 2026.

---

## Overview

PACE is an algorithmic framework for accelerated minimum inhibitory concentration (MIC) prediction from live-cell imaging of standard EUCAST broth microdilution plates. It extracts kinetic features from time-lapse absorbance measurements and systematically optimises growth/no-growth thresholds to predict MIC values hours before conventional read-out, while meeting CLSI diagnostic performance standards.

The framework is platform-agnostic: it operates on per-well absorbance time-series and does not depend on a specific microscopic imager.

Key results:
- **>5-6-fold** acceleration in time-to-result compared to standard EUCAST read-out
- *A. fumigatus*: isavuconazole (8 h), amphotericin B (9 h), voriconazole (9.5 h)
- *Candida* spp.: amphotericin B (5.5 h), fluconazole (8.5 h), micafungin (10.5 h)

---

## How it works

Broth microdilution plates are imaged at 30-minute intervals using a microscopic timelapse imager. For each well, background-corrected absorbance (BCA) time-series are used to derive kinetic features:

| Feature | Description |
|---|---|
| `BCA` | Raw background-corrected absorbance |
| `BCANormalized` | BCA normalised to growth-control well |
| `BCA_rate` | Instantaneous rate of change |
| `BCA_mean_rate` | Running mean rate |
| `BCA_cumulative` | Cumulative BCA over time |

For each antimycotic, a threshold grid search identifies the feature value and time window at which growth/no-growth classification across dilution wells best predicts the reference MIC. The best threshold is selected by maximising Essential Agreement (EA) and/or Categorical Agreement (CA) against EUCAST v12 clinical breakpoints.

An optional **bias correction** step detects and corrects systematic offsets between predicted and reference MIC distributions.

---

## Repository structure

```
.
├── set_thresh.py                    # Main pipeline script
├── cleanup_temporary_files.py       # Utility to remove intermediate .pkl files
├── commands_phase1.txt              # Example commands: threshold grid search
├── commands_phase2_aspergillus.txt  # Example commands: single-configuration runs (Aspergillus)
├── commands_phase2_candida.txt      # Example commands: single-configuration runs (Candida)
├── commands_phase1_sesa_candida.txt # Example commands: Candida (SESA panel)
├── commands_phase2_sesa_candida.txt # Example commands: Candida single-config (SESA)
├── data/
│   ├── aspergillus/                 # 30 A. fumigatus isolates (Excel input files)
│   └── candida/                     # 20 Candida spp. isolates (Excel input files)
├── job_templ/                       # SLURM job templates for HPC execution
└── pace_tools/
    ├── __init__.py
    ├── data_extraction.py           # Kinetic feature extraction from Excel files
    ├── dictionaries.py              # Reference MIC dicts, breakpoint tables, species info
    ├── evaluate.py                  # Agreement metrics (EA, CA, VME, ME) and threshold selection
    ├── predict.py                   # Growth/no-growth pattern detection and MIC prediction
    ├── plot.py                      # Result visualisation
    ├── setup.py                     # I/O directory management
    ├── phot_dict.py                 # Photometric MIC comparison utility
    └── examine_pickles.py           # Inspect pipeline result pickle files
```

---

## Data

The `data/` directory contains Excel files with raw time-lapse absorbance measurements.

**Aspergillus fumigatus** (30 isolates, `data/aspergillus/`):
- 1 ATCC quality-control strain (`ASFU_ATCC_204305_S`)
- 13 azole-resistant isolates harbouring the TR34/L98H mutation (`ASFU_R_TR34_*`)
- 8 isolates harbouring TR46/Y121F/T289A (`ASFU_R_TR46_*`)
- 1 isolate with unknown resistance mechanism (`ASFU_R_unknown_22`)
- 7 wild-type susceptible isolates (`ASFU_S_*`)

**Candida** (20 isolates, `data/candida/`):
- *C. albicans* -- 1 ATCC strain (`CDAL_ATCC_90029`)
- *C. auris* -- 3 CDC reference strains (clades I, II, V) + 5 clinical + 2 clinical (clades III, IV)
- *C. dubliniensis*, *C. glabrata*, *C. krusei*, *C. parapsilosis*, *C. tropicalis* -- ATCC and clinical isolates

Antifungal rows per plate are defined in `pace_tools/dictionaries.py` (`antimycotics_dict_A`, `antimycotics_dict_C`).

---

## Installation

Python >= 3.8 is required. Install dependencies:

```bash
pip install numpy pandas matplotlib seaborn openpyxl
```

No additional installation step is needed; `set_thresh.py` imports from the `pace_tools/` package in the same directory.

---

## Usage

### Phase 1 -- threshold grid search

Sweeps a range of thresholds across all timepoints to identify optimal parameter-threshold combinations:

```bash
# Aspergillus fumigatus, amphotericin B (row A), optimise for CA >= 0.9
./set_thresh.py \
  -g A -a A \
  -i data/aspergillus \
  -d output/genus_A_d_window_6.pkl \
  -t '(0, 172800, 1800)' \
  -s '[0.0, 1.0, 0.01]' \
  -c '[("CA", 0.9)]' \
  -p '["BCA", "BCANormalized", "BCA_rate", "BCA_mean_rate", "BCA_cumulative"]' \
  -m '[2]' \
  -n CA_0.9 \
  -l info
```

See `commands_phase1.txt` for all genus/antimycotic combinations used in the paper.

### Phase 2 -- single-configuration run

Runs a specific parameter-threshold-timewindow combination identified in phase 1:

```bash
# Aspergillus fumigatus, voriconazole (row B), bias -2, BCA threshold 0.80
./set_thresh.py \
  -g A -a B \
  -i data/aspergillus \
  -d output/genus_A_d_window_6.pkl \
  -t '(0, 172801, 1800)' \
  -s '0.80' \
  -p "['BCA']" \
  -m -2 \
  -n voriconazole_BCA_-2_0p80 \
  -l info
```

See `commands_phase2_aspergillus.txt` and `commands_phase2_candida.txt` for all configurations.

### Key arguments

| Flag | Description |
|---|---|
| `-g` | Genus: `A` (Aspergillus) or `C` (Candida) |
| `-a` | Antimycotic row label (e.g. `A`=amphotericin B, `B`=voriconazole; see `dictionaries.py`) |
| `-i` | Input directory containing Excel files |
| `-d` | Output pickle file for performance summaries |
| `-t` | Timepoint range as `(start, stop, step)` in seconds |
| `-s` | Threshold(s) to test: range `[min, max, step]` or single value |
| `-p` | Kinetic features to evaluate (list) |
| `-c` | Optimisation criterion, e.g. `[("CA", 0.9)]` or `[("EA", "EA")]` |
| `-m` | Bias correction offset(s) to apply (integer well offsets) |
| `-n` | Run label (used in output filenames) |
| `-l` | Log level (`info`, `debug`, `warning`) |
| `-o` | Output directory (defaults to current working directory) |

### HPC / SLURM

Template SLURM job scripts for parallel execution on an HPC cluster are provided in `job_templ/`. Edit the account name, email address, and environment source line for your cluster before use.

---

## Output

Each run creates a timestamped session folder containing:
- Performance summary pickle (`.pkl`) -- inspectable with `examine_pickles.py`
- Per-isolate MIC prediction plots
- Overview heatmaps of EA/CA across the threshold-timepoint grid
- Log file

```bash
# Inspect a performance summary
python pace_tools/examine_pickles.py output/my_session/performance_summary.pkl
```

---

## Citation

If you use this code or data, please cite:

> Dewaele K *et al.* Same-day MIC results in drug-resistant *Aspergillus fumigatus* and *Candida auris* using live-cell imaging of a broth microdilution assay. *Manuscript submitted for publication*, 2026.

---

## License

MIT License. See [LICENSE](LICENSE).
