#!/data/leuven/348/vsc34807/miniconda3/bin/python
import os
import glob
import pandas as pd
import pprint
import dictionaries


input_dir = input("Genus: ")

#path = input_dir.split(',')[0].strip()
#genus = input_dir.split(',')[1].strip()
genus = str(input_dir)

# EUCAST dictionaries: threshold for growth inhibition compared to growth control
# !!! SOP aflezing Aspergillus zegt: "De visueel bepaalde MIC waarde stemt overeen met de laagste concentratie van antifungaal die 95% groei inhibitie veroorzaakt (OD <5% cut-off)."
# !!! Recent EUCAST document zegt: "Spectrophotometric MIC endpoint for amphotericin and azoles against A. fumigatus: The lowest concentration of drug leading to ≥90% reduction of optical density (OD) (the following wavelengths have been used: 405, 490 or 540 nm) of that of the drug-free control is the MIC value."
eucast_A_old = {'A': 0.05, 'B': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05, 'F': 0.05}
eucast_A = {'A': 0.10, 'B': 0.10, 'C': 0.10, 'D': 0.10, 'E': 0.10, 'F': 0.10}
eucast_C = {'A': 0.10, 'B': 0.5, 'C': 0.5, 'D': 0.5, 'E': 0.5, 'F': 0.5, 'G': 0.5, 'H': 0.5}


if genus == "candida":
    eucast = eucast_C
    ref_dict = dictionaries.MIC_dict_C
    no_of_rows = 8
    index_row = 3 
    print(f'Genus is Candida')
    path = "/data/leuven/348/vsc34807/phot_can"
    rita_dict = {}
    olo_files = []
else:
    eucast = eucast_A
    ref_dict = dictionaries.MIC_dict_A
    no_of_rows = 6
    index_row = 4
    print(f'Genus is Aspergillus')
    path_olo = "/data/leuven/348/vsc34807/olorofim_ruwe_data"
    rita_dict = dictionaries.rita_dict
    path = "/data/leuven/348/vsc34807/phot_asp"
    # Use glob to get a list of all Excel files in the specified path
    olo_files = glob.glob(os.path.join(path_olo, '*'))
    olo_files = os.listdir(path_olo)

# Create a dictionary to store the results
results = {}
#print(f'rita dict: {rita_dict}')

excel_files = glob.glob(os.path.join(path, '*'))
files = pprint.pformat(excel_files)
#print(f'Excel files: {files}')
#print(f'Olo files: {olo_files}')
i = 0
# Iterate through each Excel file
for file in excel_files:
    #print(f'\nfile {file.split("/")[-1]}')
    # Read the Excel file into a dataframe
    df = pd.read_excel(file, header=index_row, usecols='A:M', index_col=0)
    #print(df)
    average_blank = float(df.iat[38, 2])
    # Create a dictionary to store the results for this file
    file_results = {}
    df = df.head(no_of_rows)

    if file.split('/')[-1] in olo_files:
    
        olo_path = os.path.join(path_olo, file.split('/')[-1])
        df_olo = pd.read_excel(olo_path, header=4, usecols='A:M', index_col=0)
        df_olo = df_olo.head(8)
        #print(f'df olo: {df_olo}')
        df_olo.set_index(df_olo.index, inplace=True)    
        #print(f"Replaced {df.iloc[5]} with {df_olo.iloc[6]}")
        df.iloc[5] = df_olo.iloc[6]
        #print(f'replacement nr {i}')
        i += 1
        #print(f'df post: {df}')

    # Iterate through each row in the dataframe
    for index, row in df.iterrows():
        #print(f'index {index}')
        index_olo = index
        # Get the value in column M for this row
        column_m_value = float(row[12])
        mic_value = 0   
        # Iterate through each column in the row (excluding column M)
        for col in row.index[:-1]:
            #print(f'col {col}')
            # If the value in the column is at least 90% less than the value in column M, update the mic_value variable
            if (row[col] - average_blank) <= (eucast[index.strip()] * (column_m_value - average_blank)):
                if index.strip() == 'D':
                    pass #print(f'value {(row[col] - average_blank)} is lower than {(eucast[index.strip()] * (column_m_value - average_blank))}: mic_value = {col}')
                mic_value = col
                if mic_value == 11:
                   mic_value = 13
            # For Aspergillus, given crystallisation in multiple drugs, comment out the following:
            elif col == 1 and index.strip() != 'A': #if we're in the first column and suppression less than 90%, break out for loop and leave mic_value 0
                # Do not apply this for amphotericin B, where "suppression" is always lower than 90% because of intrinsic absorbance
                #print(f'value {(row[col] - average_blank)} is HIGHER than {(eucast[index.strip()] * (column_m_value - average_blank))}: and column 0: mic_value remains {mic_value}')
                break
        if (index.strip() == 'F') and not (file.split('/')[-1] in olo_files) and genus == 'aspergillus':
            #print(f'skipped olorofim for file {file}')   
            index_olo = "Z"
            mic_value = 0
        # Add the mic_value to the file_results dictionary using the row index as the key
        file_results[index.strip()] = f'{index_olo.strip()}{mic_value}'
  
    # Add the file_results dictionary to the results dictionary using the file name as the key
    results[file.split('/')[-1]] = file_results

# Print the results dictionary
results_formatted = pprint.pformat(results)
print(results_formatted)
#print(type(rita_dict))
row = {}
r_number = '!'
table = []
#print(ref_dict)
for file, value in ref_dict.items():
    for a, m in value.items():
        row = {}
        print(f'{file}, {a}, {m}')
        row['file'] = file
        row['drug'] = a
        row['vis'] = m
        r = rita_dict.get(file)
        if r:
            r_number = int(r[a][1:])
            diff_rit = int(m[1:])-r_number 
            row['rit'] = f'{a+str(r_number)}'
            row['d_rit'] = diff_rit
            if abs(diff_rit) > 1:
                warning = "!!!"
            elif abs(diff_rit) == 1:
                warning = "*"
            else:
                warning = ""
            row['r*'] = warning
            if a == 'F':
                row['rit'] = 'NA'
                row['d_rit'] = 'NA'
                row['r*'] = ''
        else:
            row['rit'] = 'NA'
            row['d_rit'] = 'NA'
            row['r*'] = ''
        diff_phot = int(m[1:])-int(results[file][a][1:])
        row['fot'] = f'{results[file][a]}'
        if abs(diff_phot) > 1:
            warning = "!!!"
        elif abs(diff_phot) == 1:
            warning = "*"
        else:
            warning = ""
        row['d_fot'] = diff_phot
        row['f*'] = warning
        table.append(row)
        #if m == results[file][a]:
            #print(f'\t {a}: ok')
        #else:
            #print(f'\t {a} !!! discrepancy: {m} (visual) = {results[file][a]} (photometer)')
df = pd.DataFrame(table)

for file in ref_dict.keys():
    df_per_file = df.loc[df['file'] == file] 
    dis = df_per_file.to_string()
    print(dis)
