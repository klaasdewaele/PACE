"""Phase 2: selection of stable parameter-threshold-shift triplets from a phase-1 run.

analyze_best_params_over_time() and consolidate_periods() are copied verbatim from the analysis notebook used for the
paper (2_Parameter_triplet_selection_and_evaluation.ipynb). The phase-1 output lists, per timepoint, the best
(parameter, threshold, shift) ('best_shift_corr' == 'best'). In a sliding window of consecutive timepoints the most
frequent parameter|shift combination is taken; consecutive windows with the same combination are merged into one period,
whose threshold is the mean of the distinct best thresholds in that period (rounded to two decimals). Window: 24
timepoints (12 h) for A. fumigatus, 12 timepoints (6 h) for Candida.
"""
import numpy as np
import pandas as pd

from .common import NAMES

WINDOW = {'A': 24, 'C': 12}


# --- This function is correct and unchanged ---
def analyze_best_params_over_time(df: pd.DataFrame, window_size: int = 6, time_col: str = 'time') -> pd.DataFrame:
    """
    Analyzes the most frequent 'best' combination and also finds the most frequent 'c' label
    for that combination within the window.
    """
    best_df = df[df['best_shift_corr'] == 'best'].copy()
    #display(best_df)
    if best_df.empty:
        # This is not an error, just a state. Reduce verbose output.
        # print("No rows found with 'best_shift_corr' == 'best'. Cannot perform analysis.")
        return pd.DataFrame()

    unique_times = sorted(best_df[time_col].unique())
    if len(unique_times) < window_size:
        window_size = len(unique_times)
    if window_size == 0:
        return pd.DataFrame()

    results = []
    for i in range(len(unique_times) - window_size + 1):
        time_window_values = unique_times[i : i + window_size]
        window_df = best_df[best_df[time_col].isin(time_window_values)]
        
        if not window_df.empty:
            counts = window_df.groupby(['parameter', 'phase']).size()
            if not counts.empty:
                param, phase = counts.idxmax()
                
                dominant_rows = window_df[
                    (window_df['parameter'] == param) & (window_df['phase'] == phase)
                ]
                
                thresholds_list = dominant_rows['threshold'].tolist()
                mean_threshold = np.mean(thresholds_list) if thresholds_list else np.nan
                
                if 'c' in dominant_rows.columns and not dominant_rows['c'].empty:
                    c_label = dominant_rows['c'].astype(str).mode()[0]
                else:
                    c_label = "N/A"

                results.append({
                    'start_time': time_window_values[0],
                    'end_time': time_window_values[-1],
                    'best_parameter': param,
                    'best_phase': phase,
                    'thresholds_in_window': thresholds_list,
                    'mean_threshold': mean_threshold,
                    'c_label': c_label
                })

    return pd.DataFrame(results) if results else pd.DataFrame()

# --- This function is correct and unchanged ---
def consolidate_periods(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Consolidates consecutive windows and finds the most frequent 'c' label for the entire period.
    """
    if results_df.empty:
        return results_df

    results_df['param_phase'] = results_df['best_parameter'] + ' | ' + results_df['best_phase']
    results_df['block_id'] = (results_df['param_phase'] != results_df['param_phase'].shift()).cumsum()

    agg_rules = {
        'start_time': 'min', 'end_time': 'max', 'best_parameter': 'first',
        'best_phase': 'first', 'param_phase': 'first',
        'thresholds_in_window': 'sum', 'c_label': lambda x: list(x)
    }
    consolidated_df = results_df.groupby('block_id').agg(agg_rules).reset_index(drop=True)
    consolidated_df['thresholds_in_window'] = consolidated_df['thresholds_in_window'].apply(lambda x: list(dict.fromkeys(x)))
    consolidated_df['mean_threshold'] = consolidated_df['thresholds_in_window'].apply(np.mean)
    consolidated_df['final_c_label'] = consolidated_df['c_label'].apply(lambda x: pd.Series(x).mode()[0])
    consolidated_df = consolidated_df.drop(columns=['c_label', 'thresholds_in_window'])
    return consolidated_df


def triplet_rule(phase1_df, genus):
    """[(parameter, shift, threshold, start_s, end_s), ...] for one phase-1 table."""
    res = consolidate_periods(analyze_best_params_over_time(phase1_df, window_size=WINDOW[genus]))
    out = []
    for _, r in res.iterrows():
        shift = int(r['best_phase'].replace('shift_', ''))
        out.append((r['best_parameter'], shift, float(f"{r['mean_threshold']:.2f}"), int(r['start_time']),
                    int(r['end_time'])))
    return out


def run_label(genus, row, triplet, opt):
    """The -n label of the phase-2 run, as in commands_phase2_*.txt (opt: 'EA' or 'CA')."""
    par, shift, thr, t0, t1 = triplet
    drug = NAMES[genus][row].replace(' ', '_') if genus == 'A' else NAMES[genus][row]
    crit = 'EA_EA' if opt == 'EA' else 'CA_0.9'
    return f"{drug}_single_parameter_{par}_{shift}_{thr:.2f}".replace('.', 'p') + f"_{t0}_{t1}_{crit}"


def command_line(genus, row, triplet, opt):
    """Phase-2 command for one triplet (run from the repository root)."""
    par, shift, thr, _, _ = triplet
    data, end = ('data/aspergillus', 172801) if genus == 'A' else ('data/candida', 86400)
    return (f'./set_thresh.py -g {genus} -a {row} -i {data} -d output/genus_{genus}_d_window_6.pkl '
            f'-t "(0, {end}, 1800)" -s "{thr:.2f}" -p "[\'{par}\']" -m {shift} '
            f'-n "{run_label(genus, row, triplet, opt)}" -l info')
