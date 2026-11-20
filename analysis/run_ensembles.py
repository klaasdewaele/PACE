#!/usr/bin/env python3
"""Phase 3: ensembles, times to result and the per-isolate predictions reported in the paper.

Reads the phase-2 runs (commands_phase2_aspergillus.txt, commands_phase2_candida.txt) found below --phase2 and writes
to --out:
  <fig>_schedules.txt      triplet switches per drug (Supplementary Text S2 for Figure 5)
  <fig>_metrics.csv        per drug and timepoint: active triplet, EA, CA, VME and ME (Figure 5 / S10 curves; EUCAST ATU
                           rule and CLSI denominators, see pace_analysis/scoring.py), plus the PACE code's own scoring
  <fig>_time_to_result.csv first timepoint with EA >= 90 % (within 0-24 h and over the whole run), with EA and
                           categorical errors at that time
  <fig>_predicted_mics.csv per drug, isolate and timepoint: reference and predicted MIC and category (Table S1)
  <fig>_azole_class_calls.csv, <fig>_azole_class_summary.csv
                           A. fumigatus class-level azole call per isolate and timepoint, 0-24 h (Table S3)
  <fig>.png                the ensemble curves per drug
where <fig> is figure5 (default: EA-optimised ensembles; A. fumigatus from the EA runs, Candida from all runs) or s10
(--which s10: CA_0.9 runs only).

  python3 analysis/run_ensembles.py --phase2 PHASE2_DIR --out results/
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pace_analysis import common as K                                           # noqa: E402
from pace_analysis.ensemble import build_ensembles, schedule_lines, time_to_result   # noqa: E402
from pace_analysis.scoring import azole_class_calls, score                      # noqa: E402


def schedules_text(ens):
    keys = sorted(((f"{'Aspergillus' if g == 'A' else 'Candida'} - {K.NAMES[g][r]}", e) for (g, r), e in ens.items()),
                  key=lambda kv: kv[0])
    lines = ['ENSEMBLE MODEL SUMMARIES PER ANTIMYCOTIC', '']
    for key, e in keys:
        lines.append(f'{key}:')
        lines += [f'  - {x}' for x in schedule_lines(e)]
        lines.append('')
    return '\n'.join(lines)


def ttr_table(ens, met):
    rows = []
    for (g, r), e in ens.items():
        m = met[(met.genus == K.GENUS_NAME[g]) & (met.row == r)].set_index('time_h')
        rec = {'genus': K.GENUS_NAME[g], 'antifungal': K.NAMES[g][r],
               'time_to_result_0_24h': time_to_result(e, 24.0), 'time_to_result_any_time': time_to_result(e)}
        t = rec['time_to_result_0_24h']
        if not np.isnan(t):
            x = m.loc[t]
            rec.update({'n': int(x.n), 'EA': f'{int(x.EA_count)}/{int(x.n)}', 'CA': f'{int(x.CA_count)}/{int(x.n_CA)}',
                        'VME': f'{int(x.VME_count)}/{int(x.n_R)}', 'ME': f'{int(x.ME_count)}/{int(x.n_S)}',
                        'minor': int(x.minor_count), 'triplet': x.triplet})
        rows.append(rec)
    return pd.DataFrame(rows)


def plot(ens, met, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    keys = sorted(ens)
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.8 * nrow), sharey=True, squeeze=False)
    for ax, (g, r) in zip(axs.flat, keys):
        m = met[(met.genus == K.GENUS_NAME[g]) & (met.row == r)].sort_values('time_h')
        for col, lab, c in [('EA_pct', 'EA', 'tab:blue'), ('CA_pct', 'CA', 'tab:green'),
                            ('VME_pct', 'VME (% of R)', 'tab:red'), ('ME_pct', 'ME (% of S)', 'tab:orange')]:
            if m[col].notna().any():
                ax.plot(m.time_h, m[col], label=lab, color=c, lw=1.2)
        ax.axhline(90, color='grey', lw=0.6, ls='--')
        t = time_to_result(ens[(g, r)])
        if not np.isnan(t):
            ax.axvline(t, color='black', lw=0.8, ls=':')
        ax.set_title(f'{K.GENUS_NAME[g]} {K.NAMES[g][r]}', fontsize=9)
        ax.set_xlabel('incubation time (h)', fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in list(axs.flat)[len(keys):]:
        ax.axis('off')
    axs.flat[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phase2', nargs='+', required=True, help='folder(s) containing the phase-2 runs')
    ap.add_argument('--out', required=True)
    ap.add_argument('--which', choices=['figure5', 's10'], default='figure5')
    ap.add_argument('--no-plot', action='store_true')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    p2 = K.find_phase2(args.phase2)
    if not p2:
        sys.exit(f'no phase-2 runs found below {args.phase2}')
    ens = build_ensembles(p2, args.which)
    met, preds = score(ens)
    f = lambda name: os.path.join(args.out, f'{args.which}_{name}')
    with open(f('schedules.txt'), 'w') as fh:
        fh.write(schedules_text(ens))
    met.to_csv(f('metrics.csv'), index=False)
    ttr = ttr_table(ens, met)
    ttr.to_csv(f('time_to_result.csv'), index=False)
    preds.to_csv(f('predicted_mics.csv'), index=False)
    if args.which == 'figure5' and ('A', 'B') in ens:
        calls = azole_class_calls(preds)
        calls.to_csv(f('azole_class_calls.csv'), index=False)
        summ = calls.pivot_table(index='time_h', columns='check', values='isolate', aggfunc='count', fill_value=0)
        summ.to_csv(f('azole_class_summary.csv'))
    if not args.no_plot:
        plot(ens, met, f'{args.out}/{args.which}.png')
    mism = int(met.code_CE_mismatches_vs_stored.sum())
    print(ttr.to_string(index=False))
    print(f'\nconsistency: categorical errors recounted with the code rules differ from the stored ones at {mism} '
          f'isolate-timepoints (expected 0)')


if __name__ == '__main__':
    main()
