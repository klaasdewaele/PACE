#!/usr/bin/env python3
"""Leave-one-isolate-out cross-validation of the whole selection procedure (Supplementary Table S2).

For each drug and each isolate with a valid reference MIC, the complete selection is repeated on the other isolates:
phase-1 thresholds (EA and CA >= 0.9 rules), the phase-2 triplet rule and the phase-3 ensemble (minimum 2 h per triplet,
2 % switching tolerance). The left-out isolate is then predicted at every timepoint with the triplet that this ensemble
selects. Pooled over the left-out isolates, this gives a cross-validated EA per timepoint and a cross-validated time to
result (first timepoint, 0-24 h, with EA >= 90 %). The same procedure on all isolates gives the in-sample values.

Input: the EA-optimised phase-1 runs (commands_phase1.txt, '-n EA_EA...'), one session folder per drug, found below the
--phase1 folders. Optionally --phase2 (the phase-2 runs): the in-sample triplets and schedules are then checked against
the phase-2 based ensembles of Figure 5.

  python3 analysis/run_loocv.py --phase1 PHASE1_DIR --phase2 PHASE2_DIR --out results/
Writes table_s2_loocv.csv, loocv_by_time.csv and loocv_checks.txt to --out. Takes a few minutes.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pace_analysis import common as K                            # noqa: E402
from pace_analysis.ensemble import build_ensembles, schedule, select_runs   # noqa: E402
from pace_analysis.grid import load_grid, metrics, phase1_best, phase1_best_CA   # noqa: E402
from pace_analysis.triplets import triplet_rule                  # noqa: E402

T24 = 24.0


def triplet_id(par, thr, shift):
    """Same string as construct_ensemble_schedule's Triplet_ID."""
    return f'{par}_{thr}_{shift}'


class GridDrug:
    """The candidate triplets of one drug, evaluated on the phase-1 grid: arrays (candidate, time, isolate)."""

    def __init__(self, grid, ref, g, r, trips):
        self.iso, self.times = list(grid['iso']), grid['times']
        drug = K.NAMES[g][r].replace(' ', '_') if g == 'A' else K.NAMES[g][r]

        def file_order(x):          # the order of the phase-2 file names, in which the original analysis read them
            return (f"{x[1]}_{drug}_single_parameter_{x[0]}_{x[1]}_" + f'{x[2]:.2f}'.replace('.', 'p') +
                    f"_{x[3]}_{x[4]}_{'CA_0.9' if x[5] == 'CA' else 'EA_EA'}")
        rows, seen = [], set()
        for p, s, thr, t0, t1, opt in sorted(trips, key=file_order):
            tid = triplet_id(p, thr, s)
            if tid in seen:         # the same triplet from both criteria (identical data): keep the first
                continue
            seen.add(tid)
            rows.append({'parameter': p, 'shift': s, 'threshold': thr, 'opt': opt, 'tid': tid,
                         'si': K.SHIFTS.index(s), 'pi': K.PARAMS.index(p), 'ki': int(round(thr * 100))})
        self.cands = pd.DataFrame(rows)
        si, pi, ki = self.cands.si.values, self.cands.pi.values, self.cands.ki.values
        self.dist = np.transpose(grid['dist'][:, si, pi, ki, :], (1, 0, 2))
        self.ce = np.transpose(grid['ce'][:, si, pi, ki, :], (1, 0, 2))
        self.avail = np.ones(self.dist.shape[:2], bool)
        self.valid_ref = ref['ref_valid'].values
        self.has_bp = ref['has_bp'].values

    def metrics(self, include):
        return metrics(self.dist, self.ce, include, self.valid_ref, self.has_bp, None, None)

    def schedule(self, include):
        return schedule(self.metrics(include), self.avail, list(self.cands.tid), self.times)


def derive(grid, ref, g, include):
    """Phase-1 rules and the triplet rule on the included isolates -> candidate triplets as used in Figure 5
    (A. fumigatus: EA-optimised only; Candida: EA-optimised and CA_0.9)."""
    trips = [(p, s, t, a, b, 'EA') for p, s, t, a, b in triplet_rule(phase1_best(grid, ref, include), g)]
    if g == 'C' and ref.has_bp.any():
        trips += [(p, s, t, a, b, 'CA') for p, s, t, a, b in triplet_rule(phase1_best_CA(grid, ref, include), g)]
    return trips


def first_time(times_h, ea):
    hit = [t for t, e in zip(times_h, ea) if e >= K.T_EA and t <= T24 + 1e-9]
    return min(hit) if hit else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phase1', nargs='+', required=True, help='folder(s) containing the phase-1 EA session folders')
    ap.add_argument('--phase2', nargs='*', default=[], help='folder(s) containing the phase-2 runs (optional check)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--drugs', default='', help="subset, e.g. 'A:A,C:G' (genus:row); default all with a phase-1 run")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    drugs = [tuple(x.split(':')) for x in args.drugs.split(',')] if args.drugs else K.DRUGS
    p2 = K.find_phase2(args.phase2) if args.phase2 else {}
    fig5 = build_ensembles(p2, 'figure5') if p2 else {}

    t0 = time.time()
    table, by_time, checks = [], [], []
    for g, r in drugs:
        folder = K.find_phase1_grid(args.phase1, g, r)
        if folder is None:
            checks.append(f'{K.SHORT[(g, r)]}: no phase-1 EA run found, skipped')
            print(checks[-1], flush=True)
            continue
        grid = load_grid(folder, g, r)
        ref = K.ref_table(g, r).loc[grid['iso']]
        I = len(grid['iso'])
        th = grid['times'] / 3600.0
        tsel = np.where(th <= T24 + 1e-9)[0]

        # in-sample: the same procedure on all isolates
        all_inc = np.ones(I, bool)
        gd_all = GridDrug(grid, ref, g, r, derive(grid, ref, g, None))
        ch_all = gd_all.schedule(all_inc)
        m_all = gd_all.metrics(all_inc)
        ea_in = [m_all['EA'][ch_all[t], t] for t in tsel]
        n_in = [int(m_all['n'][ch_all[t], t]) for t in tsel]
        k_in = [int(m_all['n_EA'][ch_all[t], t]) for t in tsel]
        if (g, r) in fig5:
            e = fig5[(g, r)]
            e = e[e['time'] / 3600.0 <= T24 + 1e-9]
            tid = dict(zip((e['time'] / 3600.0).round(2), e['Triplet_ID']))
            diff = sum(tid.get(round(th[t], 2)) != gd_all.cands.tid[ch_all[t]] for t in tsel)
            p2_tids = {triplet_id(x['parameter'], x['threshold'], x['shift'])
                       for x in select_runs(p2[(g, r)], g, 'figure5')}
            same = set(gd_all.cands.tid) == p2_tids
            checks.append(f'{K.SHORT[(g, r)]}: triplets re-derived from the phase-1 grid are '
                          f"{'identical to' if same else 'DIFFERENT from'} the phase-2 runs ({len(p2_tids)}); "
                          f'in-sample schedule vs Figure 5 ensemble, 0-24 h: {diff} of {len(tsel)} timepoints differ')

        # leave one isolate out
        n = np.zeros(len(tsel), int)
        k = np.zeros(len(tsel), int)
        n_held = 0
        for i in range(I):
            if not ref.ref_valid.values[i] or np.all(np.isnan(grid['dist'][tsel][..., i])):
                continue
            n_held += 1
            inc = all_inc.copy()
            inc[i] = False
            gd = GridDrug(grid, ref, g, r, derive(grid, ref, g, inc))
            ch = gd.schedule(inc)
            for j, t in enumerate(tsel):
                c = ch[t]
                if c < 0:
                    continue
                dd = gd.dist[c, t, i]
                if np.isnan(dd):
                    continue
                n[j] += 1
                k[j] += abs(dd) <= 1
        ea_loo = [k[j] / n[j] if n[j] else np.nan for j in range(len(tsel))]
        for j, t in enumerate(tsel):
            by_time.append({'genus': K.GENUS_NAME[g], 'drug': K.NAMES[g][r], 'time_h': th[t],
                            'triplet_in_sample': gd_all.cands.tid[ch_all[t]],
                            'n_in_sample': n_in[j], 'EA_in_sample_count': k_in[j], 'EA_in_sample': ea_in[j],
                            'n_loocv': int(n[j]), 'EA_loocv_count': int(k[j]), 'EA_loocv': ea_loo[j]})
        t_in, t_loo = first_time(th[tsel], ea_in), first_time(th[tsel], ea_loo)
        rec = {'genus': K.GENUS_NAME[g], 'antifungal': K.NAMES[g][r], 'isolates': n_held,
               'time_to_result_in_sample_h': t_in, 'time_to_result_cross_validated_h': t_loo,
               'EA_at_in_sample_time_in_sample': '', 'EA_at_in_sample_time_cross_validated': ''}
        if not np.isnan(t_in):
            j = int(np.where(np.isclose(th[tsel], t_in))[0][0])
            rec['EA_at_in_sample_time_in_sample'] = f'{k_in[j]}/{n_in[j]}'
            rec['EA_at_in_sample_time_cross_validated'] = f'{k[j]}/{n[j]}'
        table.append(rec)
        print(f"{K.SHORT[(g, r)]}: in-sample {t_in} h, cross-validated {t_loo} h, held out {n_held} "
              f"({time.time() - t0:.0f} s)", flush=True)

    pd.DataFrame(table).to_csv(os.path.join(args.out, 'table_s2_loocv.csv'), index=False)
    pd.DataFrame(by_time).to_csv(os.path.join(args.out, 'loocv_by_time.csv'), index=False)
    with open(os.path.join(args.out, 'loocv_checks.txt'), 'w') as fh:
        fh.write('\n'.join(checks) + '\n')
    print('\n'.join(checks))


if __name__ == '__main__':
    main()
