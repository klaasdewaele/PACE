import numpy as np

# ---------------------------------------------------------------------------
# Antimycotic identifier mappings
# Keys correspond to well row labels on EUCAST broth microdilution plates.
# ---------------------------------------------------------------------------

antimycotics_dict_A = {
    'A': 'amphotericin B',
    'B': 'voriconazole',
    'C': 'posaconazole',
    'D': 'itraconazole',
    'E': 'isavuconazole',
    'F': 'olorofim',
}

antimycotics_dict_C = {
    'A': 'amphotericin B',
    'B': '5-flucytosine',
    'C': 'voriconazole',
    'D': 'posaconazole',
    'E': 'anidulafungin',
    'F': 'itraconazole',
    'G': 'fluconazole',
    'H': 'micafungin',
}

# ---------------------------------------------------------------------------
# EUCAST proportional inhibition thresholds
# Fraction of growth-control OD used as the MIC endpoint in the
# photometric (EUCAST) algorithm.
# ---------------------------------------------------------------------------

eucast_A = {'A': 0.10, 'B': 0.10, 'C': 0.10, 'D': 0.10, 'E': 0.10, 'F': 0.10}
eucast_C = {'A': 0.10, 'B': 0.50, 'C': 0.50, 'D': 0.50, 'E': 0.50, 'F': 0.50, 'G': 0.50, 'H': 0.50}

# ---------------------------------------------------------------------------
# EUCAST clinical breakpoint tables (EUCAST v12, 2026)
# Values are EUCAST plate well numbers (integer indices on a 12-well row).
#
# Encoding:
#   breakpoint_list[0]  → susceptible (S) if MIC well number >= this value
#   breakpoint_list[1]  → resistant (R)   if MIC well number <= this value
#   Wells between the two values correspond to ATU / intermediate category.
#
# Note for voriconazole (B) against A. fumigatus:
#   MIC 2 mg/L is formally ATU in EUCAST, but overlaps with R.
#   Here it is kept as a separate ATU category; ATU_VME (predicted S
#   instead of ATU) can be penalised as VME when optimising thresholds.
# ---------------------------------------------------------------------------

breakpoint_dict_A = {
    'A_ASFU': [5, 4],   # amphotericin B:  S ≤1, R ≥2
    'B_ASFU': [5, 4],   # voriconazole:    S ≤1, R ≥2
    'C_ASFU': [7, 6],   # posaconazole:    S ≤0.125, R ≥0.25
    'D_ASFU': [5, 4],   # itraconazole:    S ≤1, R ≥2
    'E_ASFU': [5, 4],   # isavuconazole:   S ≤1, R ≥2
    'F_ASFU': [np.nan, np.nan],  # olorofim: no EUCAST breakpoint
}

breakpoint_dict_C = {
    # C. albicans
    'A_CDAL': [4, 3],
    'B_CDAL': [np.nan, np.nan],
    'C_CDAL': [8, 5],
    'D_CDAL': [8, 7],
    'E_CDAL': [10, 9],
    'F_CDAL': [8, 7],
    'G_CDAL': [7, 5],
    'H_CDAL': [9, 8],

    # C. dubliniensis
    'A_CDDU': [4, 3],
    'B_CDDU': [np.nan, np.nan],
    'C_CDDU': [8, 5],
    'D_CDDU': [8, 7],
    'E_CDDU': [9, 8],
    'F_CDDU': [8, 7],
    'G_CDDU': [7, 5],
    'H_CDDU': [8, 7],

    # C. glabrata (H updated vs. v11)
    'A_CDGL': [4, 3],
    'B_CDGL': [np.nan, np.nan],
    'C_CDGL': [np.nan, np.nan],
    'D_CDGL': [np.nan, np.nan],
    'E_CDGL': [8, 7],
    'F_CDGL': [np.nan, np.nan],
    # EUCAST S guideline (0.001 mg/L) is below the lowest plate concentration;
    # G11 corresponds to 0.125 mg/L and is used as the S boundary.
    'G_CDGL': [11, 3],
    'H_CDGL': [8, 7],

    # C. krusei
    'A_CDKR': [4, 3],
    'B_CDKR': [np.nan, np.nan],
    'C_CDKR': [np.nan, np.nan],
    'D_CDKR': [np.nan, np.nan],
    'E_CDKR': [8, 7],
    'F_CDKR': [np.nan, np.nan],
    'G_CDKR': [np.nan, np.nan],
    'H_CDKR': [np.nan, np.nan],

    # C. parapsilosis
    'A_CDPA': [4, 3],
    'B_CDPA': [np.nan, np.nan],
    'C_CDPA': [7, 5],
    'D_CDPA': [8, 7],
    'E_CDPA': [2, 1],
    'F_CDPA': [7, 6],
    'G_CDPA': [7, 5],
    'H_CDPA': [2, 1],

    # C. tropicalis (C updated vs. v11)
    'A_CDTR': [4, 3],
    'B_CDTR': [np.nan, np.nan],
    'C_CDTR': [7, 5],
    'D_CDTR': [8, 7],
    'E_CDTR': [8, 7],
    'F_CDTR': [7, 6],
    'G_CDTR': [7, 5],
    'H_CDTR': [8, 7],

    # C. auris (E, G, H updated vs. v11)
    # EUCAST S guideline for amfotericin B (0.001 mg/L) is below the lowest
    # plate concentration; A11 (0.008 mg/L) is used as the S boundary.
    'A_CDAU': [11, 2],
    'B_CDAU': [np.nan, np.nan],
    'C_CDAU': [np.nan, np.nan],
    'D_CDAU': [np.nan, np.nan],
    'E_CDAU': [6, 5],
    'F_CDAU': [np.nan, np.nan],
    'G_CDAU': [np.nan, np.nan],  # IE in v12 guideline
    'H_CDAU': [6, 5],

    # Other Candida spp.
    'A_other': [np.nan, np.nan],
    'B_other': [np.nan, np.nan],
    'C_other': [np.nan, np.nan],
    'D_other': [np.nan, np.nan],
    'E_other': [np.nan, np.nan],
    'F_other': [np.nan, np.nan],
    'G_other': [7, 5],
    'H_other': [np.nan, np.nan],
}

# ---------------------------------------------------------------------------
# Reference MIC tables — Aspergillus fumigatus
#
# Values are well-number codes in the format <row><column> (e.g. 'B5' means
# row B, column 5 on the EUCAST plate).  'Z0' indicates that no valid MIC
# was obtained for that drug/isolate combination (excluded from analysis).
#
# Strain naming:
#   ASFU_ATCC_204305_S  — ATCC quality-control strain (susceptible)
#   ASFU_S_1 … ASFU_S_7 — clinical susceptible isolates
#   ASFU_R_TR34_*       — azole-resistant isolates with TR34/L98H mechanism
#   ASFU_R_TR46_*       — azole-resistant isolates with TR46/Y121F/M172I mechanism
#   ASFU_R_unknown_22   — azole-resistant isolate, resistance mechanism unknown
# ---------------------------------------------------------------------------

MIC_dict_A = {
    'ASFU_ATCC_204305_S': {'A': 'A6', 'B': 'B5', 'C': 'C8', 'D': 'D7', 'E': 'E4', 'F': 'F3'},
    'ASFU_R_TR34_1':  {'A': 'Z0', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'Z0'},  # amphotericin B excluded (suspect run quality)
    'ASFU_R_TR34_2':  {'A': 'A7', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E3', 'F': 'F8'},
    'ASFU_R_TR34_3':  {'A': 'A6', 'B': 'B1', 'C': 'C4', 'D': 'D0', 'E': 'E0', 'F': 'F5'},
    'ASFU_R_TR34_4':  {'A': 'A7', 'B': 'B2', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F8'},
    'ASFU_R_TR34_5':  {'A': 'A8', 'B': 'B3', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F8'},
    'ASFU_R_TR34_6':  {'A': 'A6', 'B': 'B3', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F8'},
    'ASFU_R_TR34_7':  {'A': 'A8', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E3', 'F': 'F8'},
    'ASFU_R_TR34_8':  {'A': 'A8', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F8'},
    'ASFU_R_TR34_9':  {'A': 'A7', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F9'},
    'ASFU_R_TR34_10': {'A': 'A8', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E3', 'F': 'Z0'},
    'ASFU_R_TR34_11': {'A': 'A7', 'B': 'B2', 'C': 'C5', 'D': 'D0', 'E': 'E1', 'F': 'F9'},
    'ASFU_R_TR34_12': {'A': 'A7', 'B': 'B2', 'C': 'C4', 'D': 'D0', 'E': 'E1', 'F': 'F6'},
    'ASFU_R_TR34_13': {'A': 'A6', 'B': 'B2', 'C': 'C6', 'D': 'D0', 'E': 'E2', 'F': 'Z0'},
    'ASFU_R_TR46_14': {'A': 'A6', 'B': 'B0', 'C': 'C0', 'D': 'D0', 'E': 'E0', 'F': 'Z0'},
    'ASFU_R_TR46_15': {'A': 'A6', 'B': 'B0', 'C': 'C5', 'D': 'D4', 'E': 'E0', 'F': 'F8'},
    'ASFU_R_TR46_16': {'A': 'A7', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F6'},
    'ASFU_R_TR46_17': {'A': 'A6', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'F7'},
    'ASFU_R_TR46_17_2': {'A': 'A7', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'Z0'},
    'ASFU_R_TR46_17_3': {'A': 'A7', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'Z0'},
    'ASFU_R_TR46_18': {'A': 'A6', 'B': 'B0', 'C': 'C0', 'D': 'D0', 'E': 'E0', 'F': 'Z0'},
    'ASFU_R_TR46_19': {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D5', 'E': 'E0', 'F': 'F9'},
    'ASFU_R_TR46_20': {'A': 'A6', 'B': 'B0', 'C': 'C3', 'D': 'D0', 'E': 'E0', 'F': 'Z0'},
    'ASFU_R_TR46_21': {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D0', 'E': 'E0', 'F': 'F8'},
    'ASFU_R_unknown_22': {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D6', 'E': 'E0', 'F': 'Z0'},
    'ASFU_S_1': {'A': 'A8', 'B': 'B6', 'C': 'C7', 'D': 'D7', 'E': 'E6', 'F': 'F9'},
    'ASFU_S_2': {'A': 'A6', 'B': 'B6', 'C': 'C7', 'D': 'D6', 'E': 'E6', 'F': 'Z0'},
    'ASFU_S_3': {'A': 'A7', 'B': 'B6', 'C': 'C7', 'D': 'D6', 'E': 'E4', 'F': 'F7'},
    'ASFU_S_4': {'A': 'A6', 'B': 'B6', 'C': 'C8', 'D': 'D7', 'E': 'E5', 'F': 'Z0'},
    'ASFU_S_5': {'A': 'A7', 'B': 'B6', 'C': 'C8', 'D': 'D8', 'E': 'E6', 'F': 'F8'},
    'ASFU_S_6': {'A': 'A7', 'B': 'B5', 'C': 'C6', 'D': 'D6', 'E': 'E4', 'F': 'F7'},
    'ASFU_S_7': {'A': 'A7', 'B': 'B6', 'C': 'C8', 'D': 'D8', 'E': 'E5', 'F': 'F8'},
}

# ---------------------------------------------------------------------------
# Reference MIC tables — Candida spp.
#
# Values follow the same well-number encoding as MIC_dict_A.
# Clinical isolates of C. auris (CDAU_clinical_*) are anonymous patient
# samples; all other C. auris entries are CDC reference strains.
#
# Strain naming:
#   CDAL  C. albicans   CDDU  C. dubliniensis   CDGL  C. glabrata
#   CDKR  C. krusei     CDPA  C. parapsilosis   CDTR  C. tropicalis
#   CDAU  C. auris
#   _ATCC_*  ATCC/IQC quality-control strains
#   _IQC     internal quality-control strain
#   _CDC_*   CDC reference strains
#   _clinical_*  de-identified clinical isolates
# ---------------------------------------------------------------------------

MIC_dict_C = {
    'CDAL_ATCC_90029':      {'A': 'A5', 'B': 'B6',  'C': 'C13', 'D': 'D10', 'E': 'E13', 'F': 'F10', 'G': 'G10', 'H': 'H13'},
    'CDAU_CDC_clade_I':     {'A': 'A6', 'B': 'B10', 'C': 'C9',  'D': 'D9',  'E': 'E3',  'F': 'F8',  'G': 'G7',  'H': 'H5'},
    'CDAU_CDC_clade_II':    {'A': 'A6', 'B': 'B9',  'C': 'C9',  'D': 'D10', 'E': 'E9',  'F': 'F10', 'G': 'G6',  'H': 'H9'},
    'CDAU_CDC_clade_V':     {'A': 'A5', 'B': 'B9',  'C': 'C4',  'D': 'D6',  'E': 'E4',  'F': 'F6',  'G': 'G3',  'H': 'H0'},
    'CDAU_clinical_1':      {'A': 'A5', 'B': 'B9',  'C': 'C0',  'D': 'D9',  'E': 'E5',  'F': 'F6',  'G': 'G3',  'H': 'H7'},
    'CDAU_clinical_2':      {'A': 'A4', 'B': 'B9',  'C': 'C3',  'D': 'D6',  'E': 'E3',  'F': 'F5',  'G': 'G0',  'H': 'H6'},
    'CDAU_clinical_3':      {'A': 'A4', 'B': 'B9',  'C': 'C5',  'D': 'D7',  'E': 'E1',  'F': 'F6',  'G': 'G2',  'H': 'H0'},
    'CDAU_clinical_4':      {'A': 'A6', 'B': 'B8',  'C': 'C3',  'D': 'D6',  'E': 'E5',  'F': 'F6',  'G': 'G1',  'H': 'H7'},
    'CDAU_clinical_5':      {'A': 'A5', 'B': 'B9',  'C': 'C2',  'D': 'D4',  'E': 'E5',  'F': 'F4',  'G': 'G1',  'H': 'H8'},
    'CDAU_clinical_clade_III': {'A': 'A5', 'B': 'B8', 'C': 'C0', 'D': 'D0', 'E': 'E3',  'F': 'F0',  'G': 'G0',  'H': 'H2'},
    'CDAU_clinical_clade_IV':  {'A': 'A5', 'B': 'B9', 'C': 'C9', 'D': 'D9', 'E': 'E6',  'F': 'F7',  'G': 'G7',  'H': 'H8'},
    'CDDU_clinical':        {'A': 'A6', 'B': 'B13', 'C': 'C10', 'D': 'D10', 'E': 'E13', 'F': 'F13', 'G': 'G8',  'H': 'H13'},
    'CDGL_ATCC_90030':      {'A': 'A5', 'B': 'B13', 'C': 'C6',  'D': 'D4',  'E': 'E9',  'F': 'F4',  'G': 'G4',  'H': 'H10'},
    'CDGL_clinical_1':      {'A': 'A5', 'B': 'B13', 'C': 'C6',  'D': 'D4',  'E': 'E8',  'F': 'F5',  'G': 'G5',  'H': 'H10'},
    'CDGL_clinical_2':      {'A': 'A6', 'B': 'B10', 'C': 'C6',  'D': 'D5',  'E': 'E9',  'F': 'F6',  'G': 'G5',  'H': 'H13'},
    'CDKR_ATCC_6258':       {'A': 'A5', 'B': 'B4',  'C': 'C7',  'D': 'D7',  'E': 'E9',  'F': 'F7',  'G': 'G4',  'H': 'H8'},
    'CDKR_clinical':        {'A': 'A4', 'B': 'B5',  'C': 'C5',  'D': 'D7',  'E': 'E9',  'F': 'F7',  'G': 'G3',  'H': 'H7'},
    'CDPA_ATCC_22019':      {'A': 'A5', 'B': 'B9',  'C': 'C9',  'D': 'D9',  'E': 'E4',  'F': 'F8',  'G': 'G7',  'H': 'H4'},
    'CDTR_IQC':             {'A': 'A5', 'B': 'B10', 'C': 'C8',  'D': 'D8',  'E': 'E8',  'F': 'F8',  'G': 'G9',  'H': 'H10'},
    'CDTR_clinical':        {'A': 'A5', 'B': 'B13', 'C': 'C7',  'D': 'D8',  'E': 'E9',  'F': 'F8',  'G': 'G7',  'H': 'H9'},
}

# ---------------------------------------------------------------------------
# Species lookup tables
# Map sample identifier → species code, used by the pipeline to look up the
# correct breakpoint entry in breakpoint_dict_A / breakpoint_dict_C.
# ---------------------------------------------------------------------------

species_dict_A = {
    'ASFU_ATCC_204305_S': 'ASFU',
    'ASFU_R_TR34_1':  'ASFU', 'ASFU_R_TR34_2':  'ASFU', 'ASFU_R_TR34_3':  'ASFU',
    'ASFU_R_TR34_4':  'ASFU', 'ASFU_R_TR34_5':  'ASFU', 'ASFU_R_TR34_6':  'ASFU',
    'ASFU_R_TR34_7':  'ASFU', 'ASFU_R_TR34_8':  'ASFU', 'ASFU_R_TR34_9':  'ASFU',
    'ASFU_R_TR34_10': 'ASFU', 'ASFU_R_TR34_11': 'ASFU', 'ASFU_R_TR34_12': 'ASFU',
    'ASFU_R_TR34_13': 'ASFU',
    'ASFU_R_TR46_14': 'ASFU', 'ASFU_R_TR46_15': 'ASFU', 'ASFU_R_TR46_16': 'ASFU',
    'ASFU_R_TR46_17': 'ASFU', 'ASFU_R_TR46_17_2': 'ASFU', 'ASFU_R_TR46_17_3': 'ASFU',
    'ASFU_R_TR46_18': 'ASFU', 'ASFU_R_TR46_19': 'ASFU', 'ASFU_R_TR46_20': 'ASFU',
    'ASFU_R_TR46_21': 'ASFU',
    'ASFU_R_unknown_22': 'ASFU',
    'ASFU_S_1': 'ASFU', 'ASFU_S_2': 'ASFU', 'ASFU_S_3': 'ASFU', 'ASFU_S_4': 'ASFU',
    'ASFU_S_5': 'ASFU', 'ASFU_S_6': 'ASFU', 'ASFU_S_7': 'ASFU',
}

species_dict_C = {
    'CDAL_ATCC_90029':       'CDAL',
    'CDAU_CDC_clade_I':      'CDAU',
    'CDAU_CDC_clade_II':     'CDAU',
    'CDAU_CDC_clade_V':      'CDAU',
    'CDAU_clinical_1':       'CDAU',
    'CDAU_clinical_2':       'CDAU',
    'CDAU_clinical_3':       'CDAU',
    'CDAU_clinical_4':       'CDAU',
    'CDAU_clinical_5':       'CDAU',
    'CDAU_clinical_clade_III': 'CDAU',
    'CDAU_clinical_clade_IV':  'CDAU',
    'CDDU_clinical':         'CDDU',
    'CDGL_ATCC_90030':       'CDGL',
    'CDGL_clinical_1':       'CDGL',
    'CDGL_clinical_2':       'CDGL',
    'CDKR_ATCC_6258':        'CDKR',
    'CDKR_clinical':         'CDKR',
    'CDPA_ATCC_22019':       'CDPA',
    'CDTR_IQC':              'CDTR',
    'CDTR_clinical':         'CDTR',
}

# ---------------------------------------------------------------------------
# Full list of derived parameters available for threshold optimisation.
# Base measurements (BCA, BCANormalized, SESAfungi, SESAfungiNormalized,
# TA, TANormalized) are extended with temporal derivatives (rate, mean_rate,
# cumulative, cumulative_rate, cumulative_mean_rate, cumulative_norm) and
# within-row ratios relative to neighbouring wells (_ratio suffix).
# ---------------------------------------------------------------------------

all_parameters = [
    'BCA', 'BCA_rate', 'BCA_mean_rate', 'BCA_cumulative', 'BCA_cumulative_rate',
    'BCA_cumulative_mean_rate', 'BCA_cumulative_norm',
    'BCA_ratio', 'BCA_rate_ratio', 'BCA_mean_rate_ratio', 'BCA_cumulative_ratio',
    'BCA_cumulative_rate_ratio', 'BCA_cumulative_mean_rate_ratio', 'BCA_cumulative_norm_ratio',
    'BCANormalized', 'BCANormalized_rate', 'BCANormalized_mean_rate',
    'BCANormalized_cumulative', 'BCANormalized_cumulative_rate',
    'BCANormalized_cumulative_mean_rate', 'BCANormalized_cumulative_norm',
    'BCANormalized_ratio', 'BCANormalized_rate_ratio', 'BCANormalized_mean_rate_ratio',
    'BCANormalized_cumulative_ratio', 'BCANormalized_cumulative_rate_ratio',
    'BCANormalized_cumulative_mean_rate_ratio', 'BCANormalized_cumulative_norm_ratio',
    'SESAfungi', 'SESAfungi_rate', 'SESAfungi_mean_rate', 'SESAfungi_cumulative',
    'SESAfungi_cumulative_rate', 'SESAfungi_cumulative_mean_rate', 'SESAfungi_cumulative_norm',
    'SESAfungi_ratio', 'SESAfungi_rate_ratio', 'SESAfungi_mean_rate_ratio',
    'SESAfungi_cumulative_ratio', 'SESAfungi_cumulative_rate_ratio',
    'SESAfungi_cumulative_mean_rate_ratio', 'SESAfungi_cumulative_norm_ratio',
    'SESAfungiNormalized', 'SESAfungiNormalized_rate', 'SESAfungiNormalized_mean_rate',
    'SESAfungiNormalized_cumulative', 'SESAfungiNormalized_cumulative_rate',
    'SESAfungiNormalized_cumulative_mean_rate', 'SESAfungiNormalized_cumulative_norm',
    'SESAfungiNormalized_ratio', 'SESAfungiNormalized_rate_ratio',
    'SESAfungiNormalized_mean_rate_ratio', 'SESAfungiNormalized_cumulative_ratio',
    'SESAfungiNormalized_cumulative_rate_ratio', 'SESAfungiNormalized_cumulative_mean_rate_ratio',
    'SESAfungiNormalized_cumulative_norm_ratio',
    'TA', 'TA_rate', 'TA_mean_rate', 'TA_cumulative', 'TA_cumulative_rate',
    'TA_cumulative_mean_rate', 'TA_cumulative_norm',
    'TA_ratio', 'TA_rate_ratio', 'TA_mean_rate_ratio', 'TA_cumulative_ratio',
    'TA_cumulative_rate_ratio', 'TA_cumulative_mean_rate_ratio', 'TA_cumulative_norm_ratio',
    'TANormalized', 'TANormalized_rate', 'TANormalized_mean_rate',
    'TANormalized_cumulative', 'TANormalized_cumulative_rate',
    'TANormalized_cumulative_mean_rate', 'TANormalized_cumulative_norm',
    'TANormalized_ratio', 'TANormalized_rate_ratio', 'TANormalized_mean_rate_ratio',
    'TANormalized_cumulative_ratio', 'TANormalized_cumulative_rate_ratio',
    'TANormalized_cumulative_mean_rate_ratio', 'TANormalized_cumulative_norm_ratio',
]
