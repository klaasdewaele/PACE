"""Loading phase-2 runs and building the phase-3 ensembles.

build_ensembles() reproduces how the paper's figures were made:
  - Figure 5 ('figure5'): A. fumigatus from the EA-optimised phase-2 runs only (CA_0.9 runs excluded); Candida from
    all phase-2 runs (EA and CA_0.9); construct_ensemble_schedule(min_duration_hours=2).
  - Supplementary Figure S10 ('s10'): CA_0.9 runs only, both genera.

schedule() is a numpy re-implementation of the same rule for the leave-one-isolate-out cross-validation, where the
rule is re-run thousands of times. On all isolates it gives the same schedule as construct_ensemble_schedule (checked
by run_loocv.py).
"""
import numpy as np
import pandas as pd

from .common import T_CA, T_EA, T_ME, T_VME
from .ensemble_rule import construct_ensemble_schedule, summarize_ensemble_schedule

REQUIRED = ['time', 'EA', 'CA', 'ME_rate', 'VME_rate', 'ATU_VME_rate', 'ATU_ME_rate', 'distance_per_file']


def load_drug_frame(records):
    """Concatenate the phase-2 tables of one drug, in the order of `records` (sorted file names)."""
    dfs = []
    for rec in records:
        df = pd.read_pickle(rec['path'])
        missing = [c for c in REQUIRED if c not in df.columns]
        if missing:
            raise ValueError(f"{rec['file']}: missing columns {missing}")
        df = df.copy()
        df['parameter'] = rec['parameter']
        df['threshold'] = rec['threshold']
        if 'best_shift_corr' not in df.columns:
            df['best_shift_corr'] = f"shift_{rec['shift']}"
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def select_runs(records, genus, which):
    if which == 'figure5':
        return [r for r in records if genus == 'C' or r['opt'] == 'EA']
    if which == 's10':
        return [r for r in records if r['opt'] == 'CA']
    raise ValueError(which)


def build_ensembles(phase2, which='figure5'):
    """{(genus, row): ensemble DataFrame (one row per timepoint)} from find_phase2() output."""
    out = {}
    for (g, row), records in sorted(phase2.items()):
        runs = select_runs(records, g, which)
        if not runs:
            continue
        ens = construct_ensemble_schedule(load_drug_frame(runs), min_duration_hours=2)
        if not ens.empty:
            out[(g, row)] = ens.sort_values('time').reset_index(drop=True)
    return out


def schedule_lines(ens):
    """Supplementary Text S2 format: 'At 9.0h: Switch to BCA (thr=0.93, shift_-2)'."""
    return summarize_ensemble_schedule(ens)


def time_to_result(ens, t_max_h=None):
    """First timepoint with ensemble EA >= 90 % (the Figure 5 marker); NaN if never reached."""
    e = ens if t_max_h is None else ens[ens['time'] / 3600.0 <= t_max_h + 1e-9]
    ok = np.where(e['EA'].values >= T_EA)[0]
    return float(e['time'].values[ok[0]] / 3600.0) if len(ok) else np.nan


# ---------------------------------------------------------------------------------------------------------------------
# numpy re-implementation of construct_ensemble_schedule (used by the cross-validation)
# ---------------------------------------------------------------------------------------------------------------------
def _first_argmax(v):
    """Index of the first maximum, NaN treated as -inf (= pandas sort_values(ascending=False).iloc[0] on the
    candidate sets of this study)."""
    w = np.where(np.isnan(v), -np.inf, v)
    return int(np.argmax(w))


def schedule(m, avail, tids, times, min_duration_hours=2.0, tol=0.02):
    """Chosen candidate index per timepoint (-1 where no candidate has data).

    m: dict of (candidate, time) arrays 'EA', 'CA', 'VME_rate', 'ME_rate' (rates over all categorisable isolates, as
    the PACE code stores them); avail: (candidate, time) bool; tids: triplet ID per candidate; times in seconds.
    Pass 1 selects with hysteresis; triplets whose longest streak in pass 1 is shorter than min_duration_hours are
    then excluded (unless no other candidate has data) and the selection is repeated (pass 2)."""
    C, T = avail.shape
    EA, CA, VR, MR = m['EA'], m['CA'], m['VME_rate'], m['ME_rate']
    passing = (EA >= T_EA) & (CA >= T_CA) & (VR <= T_VME) & (MR <= T_ME)    # NaN -> False, as in pandas

    def select(idx, t, cur):
        s = EA[idx, t]
        p = passing[idx, t]
        best_overall = idx[_first_argmax(s)]
        if cur is not None and cur in idx:
            if passing[cur, t]:
                return cur
            if not p.any() and EA[cur, t] >= EA[best_overall, t] - tol:
                return cur
        if p.any():
            ip = idx[p]
            return ip[_first_argmax(EA[ip, t])]
        return best_overall

    def run(allowed):
        cur, ch = None, np.full(T, -1)
        for t in range(T):
            idx = np.where(avail[:, t])[0]
            if idx.size == 0:
                continue
            fidx = idx[allowed[idx]]
            if fidx.size == 0:
                fidx = idx
            cur = select(fidx, t, cur)
            ch[t] = cur
        return ch

    ch1 = run(np.ones(C, bool))
    if min_duration_hours <= 0:
        return ch1
    ok = ch1 >= 0
    tt, cc = times[ok], ch1[ok]
    nxt = np.append(tt[1:], tt[-1] + (tt[-1] - tt[-2] if len(tt) > 1 else 0))
    dur = (nxt - tt) / 3600.0
    ids = np.array([tids[c] for c in cc])
    best = {}
    k = 0
    while k < len(ids):
        j = k
        while j + 1 < len(ids) and ids[j + 1] == ids[k]:
            j += 1
        best[ids[k]] = max(best.get(ids[k], 0.0), dur[k:j + 1].sum())
        k = j + 1
    valid = {i for i, d in best.items() if d >= min_duration_hours}
    allowed = np.array([tids[c] in valid for c in range(C)])
    return run(allowed)
