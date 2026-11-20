"""The phase-1 grid and the phase-1 selection rules, for the leave-one-isolate-out cross-validation.

A phase-1 run (commands_phase1.txt) scores every parameter x threshold (0.00-0.99) x shift (-2..+2) at every timepoint
and stores each isolate's prediction (performance_dict_<row>_<time>_EA.pkl in the run's per-timepoint folders).
Because a proportional threshold is applied to each isolate's own wells (threshold x the isolate's growth control),
an isolate's prediction does not depend on the other isolates. The whole selection procedure (phase-1 thresholds,
phase-2 triplets and the ensemble) can therefore be repeated without any one isolate, and the left-out isolate scored,
without re-running PACE.

phase1_best() and phase1_best_CA() re-implement the phase-1 selection of set_thresh.py / evaluate.py for the criteria
('EA', 'EA') and ('CA', 0.9) on any subset of isolates.
"""
import os
import pickle
import re

import numpy as np
import pandas as pd

from .common import CE_CODE, PARAMS, SHIFTS, THRS, ref_table


def load_grid(folder, genus, row):
    """Phase-1 grid of one drug as arrays dist[t, s, p, k, i] (float; NaN = no image) and ce[...] (int8 codes),
    with times in seconds and the isolate order of ref_table()."""
    iso = list(ref_table(genus, row).index)
    ix = {i: k for k, i in enumerate(iso)}
    subs = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d)) and '_bias_' in d
            and re.search(r'_time_(\d+)_bias_', d)]
    times = sorted({int(re.search(r'_time_(\d+)_bias_', d).group(1)) for d in subs})
    T, S, P, K, I = len(times), len(SHIFTS), len(PARAMS), len(THRS), len(iso)
    dist = np.full((T, S, P, K, I), np.nan, dtype=np.float32)
    ce = np.zeros((T, S, P, K, I), dtype=np.int8)
    found = np.zeros((T, S), bool)
    for d in subs:
        t = int(re.search(r'_time_(\d+)_bias_', d).group(1))
        s = int(re.search(r'_bias_(-?\d+)_', d).group(1))
        pk = os.path.join(folder, d, f'performance_dict_{row}_{t}_EA.pkl')
        if not os.path.exists(pk):
            continue
        with open(pk, 'rb') as fh:
            perf = pickle.load(fh)
        ti, si = times.index(t), SHIFTS.index(s)
        found[ti, si] = True
        for pi, par in enumerate(PARAMS):
            for thr, a in perf[par].items():
                ki = int(round(thr * 100))
                for f, v in a['distance_per_file'].items():
                    if f in ix:
                        dist[ti, si, pi, ki, ix[f]] = v
                for f, v in a['CE_per_file'].items():
                    if f in ix:
                        ce[ti, si, pi, ki, ix[f]] = CE_CODE[v]
    if not found.all():
        missing = [(times[t], SHIFTS[s]) for t, s in zip(*np.where(~found))]
        raise RuntimeError(f'incomplete phase-1 run in {folder}: {len(missing)} time x shift results missing, '
                           f'e.g. {missing[:3]}')
    return {'times': np.array(times), 'dist': dist.astype(float), 'ce': ce, 'iso': iso}


def metrics(dist, ce, include, valid_ref, has_bp, refR, refS):
    """EA / CA / VME / ME over the last axis (isolates), for the isolates in `include`, as the PACE code computes them:
    EA over isolates with data, CA and error rates over those with a breakpoint; pairs without a reference excluded."""
    data = ~np.isnan(dist) & include & valid_ref
    nEA = data.sum(-1)
    eaok = (np.abs(np.nan_to_num(dist, nan=99)) <= 1) & data
    cas = data & has_bp
    nCA = cas.sum(-1)
    out = {'n': nEA, 'n_EA': eaok.sum(-1), 'n_CA': nCA,
           'VME': ((ce == 1) & cas).sum(-1), 'ME': ((ce == 2) & cas).sum(-1),
           'ATU_VME': ((ce == 3) & cas).sum(-1), 'ATU_ME': ((ce == 4) & cas).sum(-1)}
    with np.errstate(invalid='ignore', divide='ignore'):
        out['EA'] = np.where(nEA > 0, out['n_EA'] / nEA, np.nan)
        err = out['VME'] + out['ME'] + out['ATU_VME'] + out['ATU_ME']
        out['CA'] = np.where(nCA > 0, (nCA - err) / nCA, np.nan)
        out['VME_rate'] = np.where(nCA > 0, out['VME'] / nCA, np.nan)
        out['ME_rate'] = np.where(nCA > 0, out['ME'] / nCA, np.nan)
    return out


def phase1_best(grid, ref, include=None):
    """Phase-1 per-timepoint selection, criterion ('EA', 'EA'):
    per parameter: highest EA, ties -> fewest weighted errors (2*VME + ME + 2*ATU_VME + ATU_ME), first threshold;
    across parameters: EA desc, total |distance| asc, weighted errors asc, CA desc (stable, -p order);
    across shifts -2..+2: highest EA, first shift on ties.
    Returns the columns the triplet rule needs (time, parameter, phase, threshold, best_shift_corr)."""
    inc = np.ones(len(grid['iso']), bool) if include is None else include
    r = ref.loc[grid['iso']]
    data = ~np.isnan(grid['dist']) & inc & r['ref_valid'].values
    nEA = data.sum(-1)
    EA = (np.abs(np.nan_to_num(grid['dist'], nan=99)) <= 1) & data
    with np.errstate(invalid='ignore', divide='ignore'):
        EA = EA.sum(-1) / nEA
    tot = np.where(data, np.abs(np.nan_to_num(grid['dist'])), 0).sum(-1)
    cas = data & r['has_bp'].values
    nCA = cas.sum(-1)
    ce = grid['ce']
    W = 2 * ((ce == 1) & cas).sum(-1) + ((ce == 2) & cas).sum(-1) + 2 * ((ce == 3) & cas).sum(-1) + \
        ((ce == 4) & cas).sum(-1)
    err = ((ce > 0) & cas).sum(-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        CA = np.where(nCA > 0, (nCA - err) / nCA, np.nan)
    T, S, P, K = EA.shape
    rows = []
    for t in range(T):
        best_shift, best_val = None, -np.inf
        per_shift = []
        for s in range(S):
            cand = []
            for p in range(P):
                ea = EA[t, s, p]
                mx = np.nanmax(ea)
                ks = np.where(ea == mx)[0]
                k = ks[np.argmin(W[t, s, p, ks])]
                cand.append((p, k))
            keyed = sorted(range(P), key=lambda q: (-EA[t, s, cand[q][0], cand[q][1]],
                                                    tot[t, s, cand[q][0], cand[q][1]],
                                                    W[t, s, cand[q][0], cand[q][1]],
                                                    -(np.nan_to_num(CA[t, s, cand[q][0], cand[q][1]], nan=-9))))
            p, k = cand[keyed[0]]
            per_shift.append((s, p, k))
            if EA[t, s, p, k] > best_val:
                best_val, best_shift = EA[t, s, p, k], (s, p, k)
        for j, (s, p, k) in enumerate(per_shift + [best_shift]):
            rows.append({'time': int(grid['times'][t]), 'parameter': PARAMS[p], 'threshold': float(THRS[k]),
                         'phase': f'shift_{SHIFTS[s]}', 'shift': SHIFTS[s], 'k': k, 'p': p, 's': s,
                         'best_shift_corr': 'best' if j == S else np.nan,
                         'EA': EA[t, s, p, k], 'CA': CA[t, s, p, k], 'W': W[t, s, p, k]})
    return pd.DataFrame(rows)


def _grid_arrays(grid, ref, include=None):
    inc = np.ones(len(grid['iso']), bool) if include is None else include
    r = ref.loc[grid['iso']]
    data = ~np.isnan(grid['dist']) & inc & r['ref_valid'].values
    nEA = data.sum(-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        EA = ((np.abs(np.nan_to_num(grid['dist'], nan=99)) <= 1) & data).sum(-1) / nEA
    tot = np.where(data, np.abs(np.nan_to_num(grid['dist'])), 0).sum(-1)
    cas = data & r['has_bp'].values
    nCA = cas.sum(-1)
    ce = grid['ce']
    a = {'EA': EA, 'tot': tot, 'nCA': nCA}
    for name, code in CE_CODE.items():
        a[name] = ((ce == code) & cas).sum(-1)
    a['W'] = 2 * a['VME'] + a['ME'] + 2 * a['ATU_VME'] + a['ATU_ME']
    with np.errstate(invalid='ignore', divide='ignore'):
        a['CA'] = np.where(nCA > 0, (nCA - a['VME'] - a['ME'] - a['ATU_VME'] - a['ATU_ME']) / nCA, np.nan)
    return a


def phase1_best_CA(grid, ref, include=None, target0=0.9):
    """Phase-1 per-timepoint selection, criterion ('CA', 0.9):
    target = 0.90, lowered by 0.01 until at least one threshold of any parameter has CA >= target; per parameter
    among thresholds with CA >= target, then across parameters: weighted errors, VME, ATU_VME, ME, ATU_ME,
    total |distance| (ascending), EA, CA (descending), stable; across shifts: fewest weighted errors, first shift."""
    a = _grid_arrays(grid, ref, include)
    T, S, P, K = a['EA'].shape
    rows = []
    for t in range(T):
        per_shift, best, best_w = [], None, np.inf
        for s in range(S):
            target = target0
            ca = a['CA'][t, s]
            while not (ca >= target).any():
                target = round(target - 0.01, 2)
                if target < 0:
                    break
            if target < 0:
                continue

            def keys(p, k):
                return (a['W'][t, s, p, k], a['VME'][t, s, p, k], a['ATU_VME'][t, s, p, k], a['ME'][t, s, p, k],
                        a['ATU_ME'][t, s, p, k], a['tot'][t, s, p, k], -a['EA'][t, s, p, k], -a['CA'][t, s, p, k])
            cand = []
            for p in range(P):
                ks = np.where(ca[p] >= target)[0]
                if len(ks):
                    k = min(ks, key=lambda kk: keys(p, kk))
                    cand.append((p, k))
            p, k = min(cand, key=lambda pk: keys(*pk))
            per_shift.append((s, p, k, target))
            if a['W'][t, s, p, k] < best_w:
                best_w, best = a['W'][t, s, p, k], (s, p, k, target)
        if best is None:          # no categorisable isolate at this time (possible with few isolates left out)
            continue
        for j, (s, p, k, target) in enumerate(per_shift + [best]):
            rows.append({'time': int(grid['times'][t]), 'parameter': PARAMS[p], 'threshold': float(THRS[k]),
                         'phase': f'shift_{SHIFTS[s]}', 'shift': SHIFTS[s], 'k': k, 'p': p, 's': s,
                         'best_shift_corr': 'best' if j == len(per_shift) else np.nan, 'criterium': target,
                         'EA': a['EA'][t, s, p, k], 'CA': a['CA'][t, s, p, k], 'W': a['W'][t, s, p, k]})
    return pd.DataFrame(rows)
