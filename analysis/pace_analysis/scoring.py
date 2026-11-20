"""Per-isolate predictions of an ensemble, and the categorical metrics reported in the paper.

Each ensemble row is the full phase-2 row of the active triplet. Its 'distance_per_file' gives, per isolate imaged at
that timepoint, the predicted minus the reference MIC in two-fold dilutions; as in pace_tools/predict.py:determine_CA,
the predicted plate column is the reference column minus the distance (column 1 = highest concentration).

Two scorings:
  - 'code': what set_thresh.py stores (dictionaries.py breakpoints, no ATU resolution, error rates over all
    categorisable isolates). Recounted here only as a consistency check against the stored values.
  - 'paper' (Figure 5, Tables S1 and S3):
      * EUCAST ATU for A. fumigatus (eucast_atu.py): isavuconazole 2 mg/L is S if voriconazole <= 1 mg/L, otherwise R;
        posaconazole 0.25 mg/L is S if itraconazole <= 1 mg/L, otherwise R. A reference MIC is resolved with the
        reference voriconazole/itraconazole MIC, a predicted MIC with the predicted voriconazole/itraconazole MIC of
        the same isolate at the same timepoint, from the same set of ensembles.
      * CLSI denominators: VME % = VME / reference-R isolates, ME % = ME / reference-S isolates, per timepoint, over
        isolates with a prediction and a valid reference.
      * CA: agreement over all categorisable isolates; VMEs, MEs and minor errors (disagreements involving the Candida
        I category) count as disagreements.
"""
import numpy as np
import pandas as pd

from . import eucast_atu as ATU
from .common import (FAILED, GENUS_NAME, NAMES, SPECIES, breakpoints, col_to_mg_l, mic_label, ref_col)


def code_category(g, iso, row, col):
    """Category as predict.py:determine_CA assigns it: S / R / 'ATU' (any column between the S and R breakpoint;
    never occurs for A. fumigatus, the I category for Candida)."""
    bp = breakpoints(g, iso, row)
    if bp is None:
        return None
    return 'S' if col >= bp[0] else ('R' if col <= bp[1] else 'ATU')


def paper_category(g, iso, row, col, decider_col):
    """EUCAST category with the ATU rule; 'I' for the Candida I category. decider_col: plate column of the deciding
    drug (voriconazole for isavuconazole, itraconazole for posaconazole) from the same source, or None."""
    bp = breakpoints(g, iso, row)
    if bp is None:
        return None
    rule = ATU.EUCAST_ATU.get((SPECIES[g][iso], row))
    if rule is not None:
        dec = None if decider_col is None else col_to_mg_l(g, rule['decided_by'], decider_col)
        resolved = ATU.resolve(rule, col_to_mg_l(g, row, col), dec)
        if resolved is not None:
            return resolved
    return 'S' if col >= bp[0] else ('R' if col <= bp[1] else 'I')


def classify(ref, pred):
    """'VME' (R -> S), 'ME' (S -> R), 'minor' (any other disagreement), '' (agreement)."""
    if ref == pred:
        return ''
    if ref == 'R' and pred == 'S':
        return 'VME'
    if ref == 'S' and pred == 'R':
        return 'ME'
    return 'minor'


def code_error(ref, pred):
    """predict.py:determine_CA mapping (reference, predicted category) -> VME / ME / ATU_VME / ATU_ME."""
    return {('S', 'R'): 'ME', ('R', 'S'): 'VME', ('S', 'ATU'): 'ATU_ME', ('ATU', 'R'): 'ATU_ME',
            ('ATU', 'S'): 'ATU_VME', ('R', 'ATU'): 'ATU_VME'}.get((ref, pred), '')


def predicted_columns(ensembles):
    """{(genus, row, time_s, isolate): predicted plate column, or None for a failed prediction}."""
    out = {}
    for (g, row), e in ensembles.items():
        for t, dist in zip(e['time'].astype(int), e['distance_per_file']):
            for iso, d in dist.items():
                rc = ref_col(g, iso, row)
                if rc is None:
                    raise RuntimeError(f'{iso} {NAMES[g][row]}: a pair without a reference MIC was scored; '
                                       'set_thresh.py must exclude it (get_denominators)')
                out[(g, row, t, iso)] = None if d >= FAILED else rc - int(d)
    return out


def score(ensembles):
    """Per drug x timepoint metrics (paper and code scoring) and per isolate x timepoint predictions."""
    pred = predicted_columns(ensembles)
    rows, preds = [], []
    for (g, row), e in ensembles.items():
        for _, r in e.iterrows():
            t = int(r['time'])
            ce_stored = r['CE_per_file'] if isinstance(r.get('CE_per_file'), dict) else {}
            c = dict.fromkeys(['n', 'n_EA', 'n_CA', 'n_R', 'n_S', 'n_I', 'VME', 'ME', 'minor', 'failed',
                               'code_VME', 'code_ME', 'code_ATU_VME', 'code_ATU_ME', 'ce_mismatch'], 0)
            for iso, d in sorted(r['distance_per_file'].items()):
                rc = ref_col(g, iso, row)
                pc = pred[(g, row, t, iso)]
                c['n'] += 1
                c['n_EA'] += abs(int(d)) <= 1
                rule = ATU.EUCAST_ATU.get((SPECIES[g][iso], row))
                ref_dec = ref_col(g, iso, rule['decided_by']) if rule else None
                pred_dec = pred.get((g, rule['decided_by'], t, iso)) if rule else None
                has_bp = breakpoints(g, iso, row) is not None
                rcat = paper_category(g, iso, row, rc, ref_dec) if has_bp else None
                pcat = paper_category(g, iso, row, pc, pred_dec) if (has_bp and pc is not None) else None
                err = ''
                if has_bp:
                    c['n_CA'] += 1
                    c['n_R'] += rcat == 'R'
                    c['n_S'] += rcat == 'S'
                    c['n_I'] += rcat == 'I'
                    if pc is None:
                        err = 'failed'
                        c['failed'] += 1
                    else:
                        err = classify(rcat, pcat)
                        if err:
                            c[err] += 1
                        err_code = code_error(code_category(g, iso, row, rc), code_category(g, iso, row, pc))
                        if err_code:
                            c['code_' + err_code] += 1
                        c['ce_mismatch'] += err_code != ce_stored.get(iso, '')
                preds.append({'genus': GENUS_NAME[g], 'row': row, 'drug': NAMES[g][row], 'time_h': t / 3600.0,
                              'isolate': iso, 'species': SPECIES[g][iso],
                              'ref_MIC_mg_L': mic_label(g, row, rc), 'ref_category': rcat or '',
                              'pred_MIC_mg_L': 'failed' if pc is None else mic_label(g, row, pc),
                              'pred_category': pcat or '', 'distance': int(d), 'within_EA': abs(int(d)) <= 1,
                              'error': err, 'triplet': f"{r['parameter']} thr={r['threshold']} {r['best_shift_corr']}"})
            n, nca = c['n'], c['n_CA']
            ok = nca - c['VME'] - c['ME'] - c['minor'] - c['failed']
            code_err = c['code_VME'] + c['code_ME'] + c['code_ATU_VME'] + c['code_ATU_ME']
            rows.append({
                'genus': GENUS_NAME[g], 'row': row, 'drug': NAMES[g][row], 'time_h': t / 3600.0,
                'triplet': f"{r['parameter']} thr={r['threshold']} {r['best_shift_corr']}",
                'n': n, 'EA_count': c['n_EA'], 'EA_pct': 100.0 * c['n_EA'] / n if n else np.nan,
                'n_CA': nca, 'CA_count': ok, 'CA_pct': 100.0 * ok / nca if nca else np.nan,
                'n_R': c['n_R'], 'n_S': c['n_S'], 'n_I': c['n_I'],
                'VME_count': c['VME'], 'VME_pct': 100.0 * c['VME'] / c['n_R'] if c['n_R'] else np.nan,
                'ME_count': c['ME'], 'ME_pct': 100.0 * c['ME'] / c['n_S'] if c['n_S'] else np.nan,
                'minor_count': c['minor'], 'failed_count': c['failed'],
                'code_CA_pct': 100.0 * (nca - code_err) / nca if nca else np.nan,
                'code_VME_count': c['code_VME'], 'code_ME_count': c['code_ME'],
                'code_ATU_VME_count': c['code_ATU_VME'], 'code_ATU_ME_count': c['code_ATU_ME'],
                'code_CE_mismatches_vs_stored': c['ce_mismatch'],
                'stored_EA': r['EA'], 'stored_CA': r['CA'], 'stored_VME_rate': r['VME_rate'],
                'stored_ME_rate': r['ME_rate']})
    return pd.DataFrame(rows), pd.DataFrame(preds)


# ---------------------------------------------------------------------------------------------------------------------
# Class-level azole call (A. fumigatus; Supplementary Table S3)
# ---------------------------------------------------------------------------------------------------------------------
AZOLES = {'B': 'voriconazole', 'E': 'isavuconazole', 'C': 'posaconazole'}


def combine(cats):
    """Reference class category: 'R' if any azole is R, 'S' if all are S, otherwise 'IND'."""
    if 'R' in cats:
        return 'R'
    if all(c == 'S' for c in cats):
        return 'S'
    return 'IND'


def azole_class_calls(preds, t_max_h=24.0):
    """Per A. fumigatus isolate x timepoint: S if voriconazole, isavuconazole and posaconazole are all predicted S; R if
    at least two are predicted R; otherwise no call ('nc'). Checked against the reference class category ('R' if the
    reference is R to any of the three, 'S' if S to all three; categories after ATU resolution)."""
    a = preds[(preds.genus == GENUS_NAME['A']) & preds.row.isin(list(AZOLES)) & (preds.time_h <= t_max_h + 1e-9)]
    out = []
    for (iso, t), x in a.groupby(['isolate', 'time_h']):
        if set(x.row) != set(AZOLES):
            continue
        p = dict(zip(x.row, x.pred_category))
        rf = dict(zip(x.row, x.ref_category))
        v = [p[r] for r in 'BEC']
        call = 'S' if v.count('S') == 3 else 'R' if v.count('R') >= 2 else 'nc'
        ref_class = combine([rf[r] for r in 'BEC'])
        check = ('VME' if call == 'S' and ref_class == 'R' else 'ME' if call == 'R' and ref_class == 'S'
                 else 'nc' if call == 'nc' else 'ok')
        out.append({'isolate': iso, 'time_h': t, 'voriconazole': p['B'], 'isavuconazole': p['E'],
                    'posaconazole': p['C'], 'call': call, 'reference_class': ref_class, 'check': check})
    return pd.DataFrame(out)
