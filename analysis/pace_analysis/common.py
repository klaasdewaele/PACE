"""Shared definitions: drugs, plate layout, reference MICs and categories, and discovery of PACE output files.

Reference MICs, species and EUCAST v12 breakpoints come from pace_tools/dictionaries.py, the same tables that
set_thresh.py uses. Plate column 1 holds the highest concentration; a reference code 'Z0' means that the
isolate-drug pair has no valid reference MIC and is not scored.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from pace_tools import dictionaries as D  # noqa: E402

GENUS_NAME = {'A': 'A. fumigatus', 'C': 'Candida'}
NAMES = {'A': D.antimycotics_dict_A, 'C': D.antimycotics_dict_C}
DRUGS = [('A', r) for r in 'ABCDEF'] + [('C', r) for r in 'ABCDEFGH']
SHORT = {('A', 'A'): 'Af AmB', ('A', 'B'): 'Af VOR', ('A', 'C'): 'Af POS', ('A', 'D'): 'Af ITR', ('A', 'E'): 'Af ISA',
         ('A', 'F'): 'Af OLO', ('C', 'A'): 'Ca AmB', ('C', 'B'): 'Ca 5FC', ('C', 'C'): 'Ca VOR', ('C', 'D'): 'Ca POS',
         ('C', 'E'): 'Ca ANF', ('C', 'F'): 'Ca ITR', ('C', 'G'): 'Ca FLC', ('C', 'H'): 'Ca MCF'}
# highest concentration (mg/L) per plate row = plate column 1
TOP = {('A', 'A'): 16, ('A', 'B'): 16, ('A', 'C'): 8, ('A', 'D'): 16, ('A', 'E'): 16, ('A', 'F'): 8,
       ('C', 'A'): 8, ('C', 'B'): 64, ('C', 'C'): 8, ('C', 'D'): 8, ('C', 'E'): 8, ('C', 'F'): 8,
       ('C', 'G'): 128, ('C', 'H'): 8}
MIC = {'A': D.MIC_dict_A, 'C': D.MIC_dict_C}
SPECIES = {'A': D.species_dict_A, 'C': D.species_dict_C}
BP = {'A': D.breakpoint_dict_A, 'C': D.breakpoint_dict_C}

PARAMS = ['BCA', 'BCANormalized', 'BCA_rate', 'BCA_mean_rate', 'BCA_cumulative']   # phase-1 order (-p argument)
SHIFTS = [-2, -1, 0, 1, 2]                                                          # phase-1 '-m [2]'
THRS = np.round(np.arange(0.0, 1.0, 0.01), 2)                                       # phase-1 '-s [0.0, 1.0, 0.01]'
FAILED = 13                 # distance >= 13: no unique inhibition boundary (failed prediction)
CE_CODE = {'VME': 1, 'ME': 2, 'ATU_VME': 3, 'ATU_ME': 4}
# targets hard-coded in the ensemble rule (rates over all categorisable isolates, as stored by the PACE code)
T_EA, T_CA, T_VME, T_ME = 0.90, 0.90, 0.015, 0.03

# repeat inoculations of ASFU_R_TR46_17 (repeatability plates) are not part of the 30-isolate panel
REPEATS = ('ASFU_R_TR46_17_',)


# conventional labels of the two-fold series below 0.125 mg/L
_LOW = {0.0625: '0.06', 0.03125: '0.03', 0.015625: '0.016', 0.0078125: '0.008', 0.00390625: '0.004',
        0.001953125: '0.002', 0.0009765625: '0.001'}


def fmt(x):
    return _LOW.get(x, f'{x:g}')


def mic_label(genus, row, col):
    """Plate column -> MIC label (mg/L), off-scale values as '>top' and '<=lowest'."""
    if col is None or (isinstance(col, float) and np.isnan(col)):
        return ''
    col = int(col)
    top = TOP[(genus, row)]
    if col <= 0:
        return '>' + fmt(top)
    if col >= 12:
        return '<=' + fmt(top / 2 ** 10)
    return fmt(top / 2 ** (col - 1))


def col_to_mg_l(genus, row, col):
    """Plate column -> MIC in mg/L; off-scale MICs are placed one dilution beyond the tested range."""
    top = TOP[(genus, row)]
    if col <= 0:
        return 2.0 * top
    if col >= 12:
        return top / 2 ** 11
    return top / 2 ** (col - 1)


def ref_col(genus, iso, row):
    """Reference plate column, or None for a pair without a valid reference MIC ('Z0')."""
    code = MIC[genus][iso][row]
    return None if code.startswith('Z') else int(code[1:])


def breakpoints(genus, iso, row):
    bp = BP[genus].get(f'{row}_{SPECIES[genus][iso]}', [np.nan, np.nan])
    return None if np.isnan(bp[0]) else bp


def ref_table(genus, row):
    """Per isolate: reference column, validity, breakpoints and the reference category as the PACE code assigns it
    (S / R / 'ATU' = any column between the S and R breakpoint; for Candida that is the I category)."""
    recs = []
    for iso in MIC[genus]:
        if iso.startswith(REPEATS):
            continue
        code = MIC[genus][iso][row]
        valid = not code.startswith('Z')
        col = int(code[1:])
        bp = BP[genus].get(f'{row}_{SPECIES[genus][iso]}', [np.nan, np.nan])
        has_bp = not np.isnan(bp[0])
        cat = ''
        if valid and has_bp:
            cat = 'S' if col >= bp[0] else ('R' if col <= bp[1] else 'ATU')
        recs.append({'isolate': iso, 'species': SPECIES[genus][iso], 'ref_valid': valid,
                     'ref_col': col if valid else np.nan, 'ref_col_internal': col,
                     'ref_MIC': mic_label(genus, row, col) if valid else 'n.a.',
                     'has_bp': has_bp, 'bp_S_col': bp[0], 'bp_R_col': bp[1], 'ref_cat': cat})
    return pd.DataFrame(recs).set_index('isolate')


# ---------------------------------------------------------------------------------------------------------------------
# PACE output files
# ---------------------------------------------------------------------------------------------------------------------
P2_PATTERN = re.compile(r'performance_(A|C)_([A-H])_.*?_single_parameter_(.+?)_(-?\d+)_(\d+p\d+)_(\d+)_(\d+)_(.+?)'
                        r'_(\d{2}-\d{2}-\d{4}-\d{2}-\d{2}-\d{2})\.pkl$')


def find_phase2(folders):
    """All phase-2 performance pickles below `folders` (set_thresh.py writes one per run, with the session time stamp
    in its name). Returns {(genus, row): [record, ...]} with records sorted by file name, the order in which the
    original analysis read them. If the same run label occurs more than once, the most recent run is kept."""
    found = {}
    for top in folders:
        for dirpath, _, files in os.walk(top):
            for fn in files:
                m = P2_PATTERN.search(fn)
                if not m:
                    continue
                g, row, par, shift, thr, w0, w1, opt, stamp = m.groups()
                label = fn[:m.start(9) - 1]
                rec = {'file': fn, 'path': os.path.join(dirpath, fn), 'genus': g, 'row': row, 'parameter': par,
                       'shift': int(shift), 'threshold': float(thr.replace('p', '.')),
                       'window_s': (int(w0), int(w1)), 'opt': 'CA' if 'CA_0.9' in opt else 'EA',
                       'stamp': pd.to_datetime(stamp, format='%d-%m-%Y-%H-%M-%S')}
                old = found.get(label)
                if old is None or rec['stamp'] > old['stamp']:
                    found[label] = rec
    out = {}
    for rec in sorted(found.values(), key=lambda r: r['file']):
        out.setdefault((rec['genus'], rec['row']), []).append(rec)
    return out


def find_phase1_grid(folders, genus, row):
    """Session folder of the EA-optimised phase-1 run of one drug (commands_phase1.txt, -n EA_EA...), or None."""
    hits = []
    for top in folders:
        for dirpath, dirnames, _ in os.walk(top):
            for d in dirnames:
                if d.startswith(f'{genus}_{row}_EA_time_0-') and '_bias_[' in d:
                    hits.append(os.path.join(dirpath, d))
            dirnames[:] = [d for d in dirnames if '_time_' not in d]     # do not descend into run folders
    return max(hits, key=os.path.getmtime) if hits else None           # most recent run
