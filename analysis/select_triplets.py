#!/usr/bin/env python3
"""Phase 2: select the stable triplets from phase-1 runs and print the phase-2 commands.

Input: the phase-1 performance tables that set_thresh.py writes into each phase-1 session folder
(performance_<genus>_<row>_EA_time_..._EA_EA_<time stamp>.pkl, or ..._0.9_time_..._CA_0.9_... for the CA >= 0.9 runs),
given as files or as folders to search.

  python3 analysis/select_triplets.py PHASE1_DIR            # all phase-1 runs below PHASE1_DIR
  python3 analysis/select_triplets.py run.pkl --commands    # print the command lines only

The rule (pace_analysis/triplets.py) takes, in a sliding window of consecutive timepoints, the most frequent best
parameter|shift and merges consecutive windows into periods; the threshold of a period is the mean of its distinct best
thresholds. With the paper's phase-1 runs this gives the triplets of commands_phase2_aspergillus.txt and
commands_phase2_candida.txt.
"""
import argparse
import os
import re
import sys

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pace_analysis.common import NAMES                                  # noqa: E402
from pace_analysis.triplets import command_line, run_label, triplet_rule   # noqa: E402

P1 = re.compile(r'^performance_(A|C)_([A-H])_(EA|0\.9)_time_\d+-\d+_bias_\[\d\]_.*_\d{2}-\d{2}-\d{4}-\d{2}-\d{2}-\d{2}\.pkl$')
P1_LOOSE = re.compile(r'^performance_(A|C)_([A-H])_(EA|0\.9)_time_\d+-\d+_bias_.*\.pkl$')   # renamed copies


def find(paths):
    out = []
    for p in paths:
        if os.path.isfile(p):
            out.append(p)
            continue
        for dirpath, _, files in os.walk(p):
            out += [os.path.join(dirpath, f) for f in files if P1.match(f)]
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--commands', action='store_true', help='print only the phase-2 command lines')
    args = ap.parse_args()
    files = find(args.paths)
    if not files:
        sys.exit('no phase-1 performance tables found')
    for path in files:
        m = P1.match(os.path.basename(path)) or P1_LOOSE.match(os.path.basename(path))
        if not m:
            sys.exit(f'not a phase-1 performance table: {path}')
        genus, row, crit = m.groups()
        opt = 'EA' if crit == 'EA' else 'CA'
        trips = triplet_rule(pd.read_pickle(path), genus)
        if not args.commands:
            print(f"\n# {NAMES[genus][row]} ({'A. fumigatus' if genus == 'A' else 'Candida'}, row {row}), "
                  f"{'EA' if opt == 'EA' else 'CA >= 0.9'} phase-1 run: {os.path.basename(path)}")
            for t in trips:
                print(f'#   {t[0]} shift {t[1]:+d} threshold {t[2]:.2f}, selected {t[3] / 3600:g}-{t[4] / 3600:g} h  '
                      f'-> {run_label(genus, row, t, opt)}')
        for t in trips:
            print(command_line(genus, row, t, opt))


if __name__ == '__main__':
    main()
