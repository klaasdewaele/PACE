#!/usr/bin/env python3
"""
phot_dict.py — photometric MIC derivation utility.

Reads time-lapse absorbance Excel files from the oCelloscope instrument,
applies the EUCAST proportional-inhibition endpoint, and prints a table
comparing visually read reference MICs (from dictionaries.MIC_dict_*)
against photometrically derived MIC values.

Usage:
    Run interactively and enter 'candida' or 'aspergillus' when prompted
    for the genus, then supply the path to the folder containing the raw
    absorbance Excel files.
"""

import glob
import os

import pandas as pd

import dictionaries

genus = input("Genus (candida / aspergillus): ").strip().lower()

if genus == "candida":
    eucast = dictionaries.eucast_C
    ref_dict = dictionaries.MIC_dict_C
    no_of_rows = 8
    index_row = 3
    path = input("Path to Candida photometric data folder: ").strip()
    path_olo = None
    olo_files = []
else:
    eucast = dictionaries.eucast_A
    ref_dict = dictionaries.MIC_dict_A
    no_of_rows = 6
    index_row = 4
    path = input("Path to Aspergillus photometric data folder: ").strip()
    path_olo = input("Path to olorofim photometric data folder (leave blank to skip): ").strip()
    olo_files = os.listdir(path_olo) if path_olo else []

# ---------------------------------------------------------------------------
# Read absorbance data and derive photometric MIC per file and drug row
# ---------------------------------------------------------------------------
results = {}
excel_files = glob.glob(os.path.join(path, '*'))

for file in excel_files:
    df = pd.read_excel(file, header=index_row, usecols='A:M', index_col=0)
    average_blank = float(df.iat[38, 2])
    file_results = {}
    df = df.head(no_of_rows)

    # For Aspergillus: optionally replace the olorofim row (F) with data
    # from a dedicated olorofim plate when available.
    filename = file.split('/')[-1]
    if genus == 'aspergillus' and filename in olo_files and path_olo:
        olo_path = os.path.join(path_olo, filename)
        df_olo = pd.read_excel(olo_path, header=4, usecols='A:M', index_col=0)
        df_olo = df_olo.head(8)
        df.iloc[5] = df_olo.iloc[6]

    for index, row in df.iterrows():
        index_key = index.strip()   # drug row label used as dict key
        index_out = index_key       # row label used in output value (may become 'Z' for olorofim)
        column_m_value = float(row[12])
        mic_value = 0

        for col in row.index[:-1]:
            # Identify the lowest concentration where OD is suppressed to
            # <= eucast threshold fraction of the growth-control OD.
            if (row[col] - average_blank) <= (eucast[index_key] * (column_m_value - average_blank)):
                mic_value = col
                if mic_value == 11:
                    mic_value = 13  # well 11 re-mapped to 13 (plate layout)
            elif col == 1 and index_key != 'A':
                # Amphotericin B (row A) is excluded from early-break logic
                # because intrinsic absorbance prevents full 90% suppression
                # at well 1 for non-susceptible isolates.
                break

        # Mark olorofim (row F) as not determined when no dedicated plate is used.
        if index_key == 'F' and filename not in olo_files and genus == 'aspergillus':
            index_out = 'Z'
            mic_value = 0

        file_results[index_key] = f'{index_out}{mic_value}'

    results[filename] = file_results

# ---------------------------------------------------------------------------
# Compare photometric MICs against visual reference MICs
# ---------------------------------------------------------------------------
table = []
for file, value in ref_dict.items():
    for drug, vis_mic in value.items():
        phot_mic = results.get(file, {}).get(drug, 'N/A')
        row_data = {'file': file, 'drug': drug, 'vis': vis_mic, 'phot': phot_mic}
        if phot_mic != 'N/A':
            diff = int(vis_mic[1:]) - int(phot_mic[1:])
            row_data['d'] = diff
            row_data['flag'] = '!!!' if abs(diff) > 1 else ('*' if abs(diff) == 1 else '')
        else:
            row_data['d'] = 'N/A'
            row_data['flag'] = ''
        table.append(row_data)

df_out = pd.DataFrame(table)
for file in ref_dict.keys():
    print(df_out.loc[df_out['file'] == file].to_string())
