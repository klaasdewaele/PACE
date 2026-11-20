"""Phase 3: the ensemble rule of the paper, copied verbatim from the analysis code that produced Figure 5 and
Supplementary Text S2 (statistical_utils_v2.construct_ensemble_schedule and summarize_ensemble_schedule).

construct_ensemble_schedule(df, min_duration_hours=2) is the published setting (minimum 2 h per triplet, 2 % EA
switching tolerance). Its input is the concatenation of the phase-2 performance tables of one drug (one row per
triplet and timepoint); see ensemble.load_drug_frame().
"""
import pandas as pd


def construct_ensemble_schedule(df, min_duration_hours=2, hysteresis_tolerance=0.02):
    """
    Constructs a dynamic ensemble schedule based on CLSI criteria with Hysteresis.
         
    Args:
        df: DataFrame containing performance data.
        min_duration_hours: Minimum continuous duration (hours) to keep a triplet.
        hysteresis_tolerance: The EA margin required to force a switch when failing.
    """
    df = df.copy()
    df['time'] = pd.to_numeric(df['time'])
    unique_times = sorted(df['time'].unique())
         
    # Add Triplet_ID upfront to make tracking easier
    def get_triplet_id(row):
        b_val = row.get('best_shift_corr', 'N/A')
        if 'shift_' in str(b_val): b_val = str(b_val).replace('shift_', '')
        return f"{row['parameter']}_{row['threshold']}_{b_val}"
             
    df['Triplet_ID'] = df.apply(get_triplet_id, axis=1)
         
    def select_best_triplet(candidate_df, current_id=None):
        if candidate_df.empty:
            return None, None, None
                     
        pass_mask = (
            (candidate_df['EA'] >= 0.90) &
            (candidate_df['CA'] >= 0.90) &
            (candidate_df['VME_rate'] <= 0.015) &
            (candidate_df['ME_rate'] <= 0.03)
        )
                 
        passing_df = candidate_df[pass_mask]
        best_overall = candidate_df.sort_values('EA', ascending=False).iloc[0]
                 
        # --- HYSTERESIS LOGIC (STICKINESS) ---
        if current_id is not None:
            current_rows = candidate_df[candidate_df['Triplet_ID'] == current_id]
            if not current_rows.empty:
                current_row = current_rows.iloc[0]
                                 
                # Rule 1: If current parameter is STILL passing CLSI criteria, KEEP IT.
                if current_id in passing_df['Triplet_ID'].values:
                    return current_row, 'PASS', ''
                                     
                # Rule 2: If NO parameters are passing, but current parameter is within 
                # tolerance of the absolute best, KEEP IT to avoid jitter.
                if passing_df.empty:
                    if current_row['EA'] >= (best_overall['EA'] - hysteresis_tolerance):
                        reason = 'FAIL (Kept via Hysteresis)'
                        return current_row, 'FAIL', reason
         
        # --- NORMAL SELECTION (If we couldn't stick to the current one) ---
        if not passing_df.empty:
            # Pick the best passing one
            best_passing = passing_df.sort_values('EA', ascending=False).iloc[0]
            return best_passing, 'PASS', ''
        else:
            # Pick the best overall and log failure reasons
            reasons = []
            if best_overall['EA'] < 0.90: reasons.append('EA<90%')
            if best_overall['CA'] < 0.90: reasons.append('CA<90%')
            if best_overall['VME_rate'] > 0.015: reasons.append('VME>1.5%')
            if best_overall['ME_rate'] > 0.03: reasons.append('ME>3%')
            return best_overall, 'FAIL', ', '.join(reasons)
     
    # --- Pass 1: Initial Selection with Hysteresis ---
    raw_results = []
    current_id = None
         
    for t in unique_times:
        current_df = df[df['time'] == t]
        if current_df.empty: continue
                     
        best_row, status, reason = select_best_triplet(current_df, current_id)
        current_id = best_row['Triplet_ID'] # Update currently active triplet
                 
        result = best_row.to_dict()
        result['Ensemble_Status'] = status
        result['Failure_Reason'] = reason
        result['Bias_Value'] = result['Triplet_ID'].split('_')[-1]
        raw_results.append(result)
             
    df_raw = pd.DataFrame(raw_results)
         
    if min_duration_hours <= 0 or df_raw.empty:
        df_raw['Ensemble_Phase'] = df_raw.apply(
            lambda r: f"{r['parameter']} (t={r['threshold']}, b={r['Bias_Value']})", axis=1)
        return df_raw
     
    # --- Duration Calculation & Filtering ---
    df_raw = df_raw.sort_values('time')
    df_raw['next_time'] = df_raw['time'].shift(-1)
    if len(df_raw) > 1:
        last_interval = df_raw['time'].iloc[-1] - df_raw['time'].iloc[-2]
        df_raw.loc[df_raw.index[-1], 'next_time'] = df_raw['time'].iloc[-1] + last_interval
    else:
        df_raw['next_time'] = df_raw['time']
     
    df_raw['duration_hours'] = (df_raw['next_time'] - df_raw['time']) / 3600.0
    df_raw['block'] = (df_raw['Triplet_ID'] != df_raw['Triplet_ID'].shift()).cumsum()
    block_durations = df_raw.groupby(['Triplet_ID', 'block'])['duration_hours'].sum().reset_index()
    max_streak_durations = block_durations.groupby('Triplet_ID')['duration_hours'].max()
    valid_triplets = max_streak_durations[max_streak_durations >= min_duration_hours].index.tolist()
         
    # --- Pass 2: Re-selection with Filter AND Hysteresis ---
    final_results = []
    current_id = None # Reset state for pass 2
         
    for t in unique_times:
        current_df = df[df['time'] == t].copy()
        if current_df.empty: continue
                 
        filtered_df = current_df[current_df['Triplet_ID'].isin(valid_triplets)]
                 
        # If no valid triplets exist at this time, fall back to all candidates
        target_df = filtered_df if not filtered_df.empty else current_df
                 
        best_row, status, reason = select_best_triplet(target_df, current_id)
        current_id = best_row['Triplet_ID'] # Update currently active triplet
                     
        result = best_row.to_dict()
        result['Ensemble_Status'] = status
        result['Failure_Reason'] = reason
        result['Bias_Value'] = result['Triplet_ID'].split('_')[-1]
        result['Ensemble_Phase'] = f"{result['parameter']} (t={result['threshold']}, b={result['Bias_Value']})"
                 
        final_results.append(result)
             
    return pd.DataFrame(final_results)


def summarize_ensemble_schedule(ensemble_df):
    """
    Returns a human-readable list of strings describing the ensemble model transitions.
    """
    if 'time' not in ensemble_df.columns or 'parameter' not in ensemble_df.columns:
        return ["No valid ensemble data."]
        
    summary = []
    current_model = None
    
    # Sort by time just in case
    df = ensemble_df.sort_values('time')
    
    for i, row in df.iterrows():
        # Identify the model by param, threshold, and shift
        model_str = f"{row['parameter']} (thr={row['threshold']}, {row.get('best_shift_corr', 'no_shift')})"
        
        if model_str != current_model:
            summary.append(f"At {row['time'] / 3600.0:.1f}h: Switch to {model_str}")
            current_model = model_str
            
    return summary
