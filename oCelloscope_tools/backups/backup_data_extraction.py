import pandas as pd
import numpy as np
import os
import math
import pprint
from datetime import datetime
import statistics
import logging
import itertools
import pickle


def get_add_par(y, window, x):
    # Rate of change in parameter: y_rate
    y_rate = [] # np.array needed for plotting ('.empty()' function); here lists are created that are converted to arrays: https://stackoverflow.com/questions/10121926/initialise-numpy-array-of-unknown-length
    y_rate_old = []
    y_rate_old_time_corrected = []

    if not y[0] == 0:
        previous = y[0]
    else:
        previous = y[1]

    previous_old = y[0]
    previous_x = x[0]
    for index, row in enumerate(y[1:]): # y: start with second value in dataframe for calculating difference with first value
        # for x (time), similarly take difference between second and first value for normalising time
        y_rate_old.append(row - previous_old)
        # probably better to use rate normalised for start value and for time units:
        #logging.debug(f'y {row}, previous {previous}, diff {row - previous}, time {x[index + 1]}, previous {previous_x}: diff {x[index + 1] - previous_x}')
        y_rate.append(float(((row - previous) / previous) / (x[index + 1] - previous_x)))
        y_rate_old_time_corrected.append(float((row - previous_old) / (x[index + 1] - previous_x))) 
        previous = row
        previous_old = row
        previous_x = x[index + 1]
    # y_rate contains one entry less than y (since calculated between 1st and 2nd entry and not between 1th and 0th entry)
    # for plotting purposes duplicate first entry of y_rate (alternative 0 or np.NaN give trouble with further calculations)
    y_rate[0:0] = [y_rate[0]]
    y_rate_old[0:0] = [y_rate_old[0]]
    y_rate_old_time_corrected[0:0] = [y_rate_old_time_corrected[0]]

    # Mean rate of change in a window of n values:
    y_mean_rate = []
    for index in range(len(y_rate)):
        values_in_window = y_rate[max(index - (window - 1), 0):index + 1]
        mean = statistics.mean(values_in_window)
        y_mean_rate.append(mean)
        #logging.debug(f'y is {y[index]}, y_rate is {y_rate[index]}, window {index + 1}/{len(y_rate)}: {values_in_window}; average: {mean}')#; y_mean_rate length {len(y_mean_rate)}')
    #logging.debug(f'y mean rate: {y_mean_rate}')
        
    """ Example with window size 3
    y: [1, 2, 4, 7, 11, 16, 22, 29, 37]; length 9
    y_rate: [1, 2, 3, 4, 5, 6, 7, 8]; length 8
    y_rate extended: [1, 1, 2, 3, 4, 5, 6, 7, 8]; length 9
    y is 1, y_rate is 1, window 1/9: [1]; average: 1
    y is 2, y_rate is 1, window 2/9: [1, 1]; average: 1
    y is 4, y_rate is 2, window 3/9: [1, 1, 2]; average: 1.33
    y is 7, y_rate is 3, window 4/9: [1, 2, 3]; average: 1.75
    y is 11, y_rate is 4, window 5/9: [2, 3, 4]; average: 2.5
    y is 16, y_rate is 5, window 6/9: [3, 4, 5]; average: 3.5
    y is 22, y_rate is 6, window 7/9: [4, 5, 6]; average: 4.5
    y is 29, y_rate is 7, window 8/9: [5, 6, 7]; average: 5.5
    y is 37, y_rate is 8, window 9/9: [6, 7, 8]; average: 6.5
    """
    # Cumulative values for parameter
    y_cumulative = []
    previous = 0
    for row in y:
        sum_ = row + previous#np.nan_to_num(previous)
        y_cumulative.append(sum_)
        previous = sum_
    
    # Cumulative value for rate of change in parameter
    y_cumulative_rate = []
    previous = 0
    for row in y_rate:
        sum_ = row + previous#np.nan_to_num(previous)
        y_cumulative_rate.append(sum_)
        previous = sum_
    
    # Cumulative value for mean rate of change
    y_cumulative_mean_rate = []
    previous = 0
    for row in y_mean_rate:
        sum_ = row + previous#np.nan_to_num(previous)
        y_cumulative_mean_rate.append(sum_)
        previous = sum_

    # Cumulative value for parameter, normalized according to first measurement (normalized as a ratio of first measurement)
    # Determine first non-zero value (BCANormalized, SESAfungiNormalized, … are normalized by subtracting first value from all other values)
    for x in y_cumulative:
        if x != 0: 
            first_non_zero = x
            break
    y_cumulative_norm = [value / first_non_zero for value in y_cumulative]

    # For reviewing obtained data in columns besides each other
    """juxta = pd.DataFrame()
    juxta['y'] = y
    juxta['y_cumulative'] = y_cumulative
    juxta['y_cumulative_norm'] = y_cumulative_norm
    juxta['y_rate'] = y_rate
    juxta['y_cumulative_rate'] = y_cumulative_rate
    juxta['y_mean_rate'] = y_mean_rate
    juxta['y_cumulative_mean_rate'] = y_cumulative_mean_rate
    pprint.pprint(juxta)"""
    
    return np.array(y_rate), np.array(y_mean_rate), np.array(y_cumulative), np.array(y_cumulative_rate), np.array(y_cumulative_mean_rate), np.array(y_cumulative_norm), np.array(y_rate_old), np.array(y_rate_old_time_corrected)

def extract_rate_cum(d, antimycotic, file, well, parameter, data, window):

    # Get time column (x) for each file in order to calculate rate per time unit
    x = d[file]

    (d[antimycotic][file][well][f'{parameter}_rate'],
        d[antimycotic][file][well][f'{parameter}_mean_rate'],
        d[antimycotic][file][well][f'{parameter}_cumulative'],
        d[antimycotic][file][well][f'{parameter}_cumulative_rate'],
        d[antimycotic][file][well][f'{parameter}_cumulative_mean_rate'],
        d[antimycotic][file][well][f'{parameter}_cumulative_norm'],
        d[antimycotic][file][well][f'{parameter}_rate_old'],
        d[antimycotic][file][well][f'{parameter}_rate_old_time_corrected']) = get_add_par(data, window, x)
    
def calculate_ratios(d, antimycotic, file, well, parameter):
    global filtered_parameters
    logging.debug(f'--- --- --- --- calculate_ratios() for {antimycotic}, {file}, {well}, {parameter}')
    # for files in training_files, well is 'MIC±x' - for others this is just a number [1-12]
    # At this point we have initialized all derived parameters, including the ones with ratios, as empty np.arrays at d[antimycotic][file][well] for each well
    # While populating these ratios per parameter, avoid creating BCA_ratio_ratio: therefore filter out keys with 'ratio'
    filtered_parameters = [p for p in list(d[antimycotic][file][well].items()) if (parameter == p[0].split('_')[0]) and ('ratio' not in p[0])] 
    # Note: filtered_parameters both contains filtered keys AND values in the format ([key,value], [key,value])
    for [par, par_values] in filtered_parameters: # note: iterating through d as list, otherwise error: https://www.delftstack.com/howto/python/python-dictionary-changed-size-during-iteration/        
        logging.debug(f'--- --- --- --- --- par {par}, par_val[:3] {par_values[:3]}')
        # extract integer value of well number (e.g. MIC-10: -10; well name 10, must remain 10)
        # wells are in format [1, 2, 3, 4, …] or [MIC-2, MIC-1, MIC, MIC+1, …]: isolate prefix ('MIC' or '') from well number (1, 2, 3, 4 or -2, -1, '', 1, …)
        prefix_list = []
        [prefix_list.append(char) for char in well if (not char.isnumeric()) and char != '-']
        prefix = "".join(prefix_list)
        suffix_list = []
        [suffix_list.append(char) for char in well if (char.isnumeric()) or char == '-']
        suffix = int("".join(suffix_list) or 0) # trouble when int('') in case of 'MIC': https://stackoverflow.com/questions/25388166/convert-empty-string-to-zero
        # Determine which neighbouring well has a concentration one log2 dilution higher than current well (leftside well in EUCAST configuration)
        # For wells expressed relative to MIC, this is the well number plus one (higher concentrations are expressed as positive values compared to MIC)
        # For wells with an absolute well number (A1, A2, A3, …) this is the well number minus one
        if prefix == 'MIC': 
            leftside_well = f'MIC{suffix + 1}'
            if leftside_well == 'MIC0': 
                leftside_well = 'MIC'
            # Problem when suffix is -1, you get 'MIC0': change to 'MIC'
        else: 
            leftside_well = str(suffix - 1)
 
        df_ratio = pd.DataFrame()
        logging.debug(f'--- --- --- --- --- --- well {well}, leftside: {leftside_well}')
        if (d[antimycotic][file][leftside_well][par].size == 0):
            denominator = np.zeros_like(par_values)
        else:
            denominator = d[antimycotic][file][leftside_well][par]
        if (par_values.size != 0):# and (d[antimycotic][file][leftside_well][par].size != 0):
            df_ratio = np.divide(par_values, denominator, out=np.full_like(par_values, 1), where=(denominator!=0))
            #https://www.semicolonworld.com/question/54586/how-to-return-0-with-divide-by-zero
            #juxta = pd.DataFrame()
            #juxta['par_values'] = par_values
            #juxta['leftside well'] = denominator 
            #juxta['ratio'] = df_ratio
            #logging.debug(f'juxta: {juxta}')
            
        par_without_ratio = str(par.split('_ratio')[0]) 
        # At this point we have initialized all derived parameters, including the ones with ratios, as empty np.arrays
        # Now we iterate through all par and par_values for this well in d and therefore also encounter e.g. BCA_ratio, etc
        d[antimycotic][file][well][f'{par_without_ratio}_ratio'] = df_ratio
        #print(f'{par_without_ratio}_ratio')
        #pprint.pprint(d[antimycotic][file][well][f'{par_without_ratio}_ratio'][0:3])
                     
def extract_parameter(d, df, parameter, file, antimycotics, window, well_names, MIC_dict, training_files):
    # Data per parameter spans over a number of rows in the source Excel file: determine start and stop row
    start = 0
    stop = 0
    for index, row in df.iterrows():   
        if (row['Features'] == parameter):  
            start = index + 2 # Values start two rows below index of parameter name
            for i, r in df.iloc[index + 2:].iterrows(): # iterating through rest of dataframe, stop is index of first empty row
                # Iterrows() yields Series object with i as index and r as values per row WITH column headers
                if (isinstance(r['Features'], float)):
                    if np.isnan(r['Features']): # Empty cell detected after series of values
                        stop = i
                        break
    
    # CREATE DATAFRAME FOR PRESENT PARAMETER WITH APPROPRIATE HEADERS
    header_row = 2 # the row index where the headers are in the source dataframe
    df_parameter = df.copy()[start:stop] 
    df_parameter.columns = df.iloc[header_row] # Create column names from row index 2 in dataframe
        
    # Add time column to d at highest level if not yet done (this has to be done only once for each file)
    if d[file].size == 0:
        d[file] = df_parameter[['Time (seconds)']].values

    # ITERATE THROUGH COLUMNS; IDENTIFY MIC-COLUMN AND NAME COLUMNS RELATIVE TO MIC-COLUMN PER ANTIMYCOTIC
    for antimycotic in antimycotics.keys():
        logging.info(f'--- --- Antimycotic: {antimycotic}')
        
        # Create entries for all antimycotic-parameter combinations, initialize as empty np.array
        derived_parameters = [f'{parameter}_rate', f'{parameter}_mean_rate', f'{parameter}_cumulative', 
                              f'{parameter}_cumulative_rate', f'{parameter}_cumulative_mean_rate', 
                              f'{parameter}_cumulative_norm', f'{parameter}_ratio', f'{parameter}_rate_old', f'{parameter}_rate_old_time_corrected',
                              f'{parameter}_rate_ratio', f'{parameter}_mean_rate_ratio', 
                              f'{parameter}_cumulative_ratio', f'{parameter}_cumulative_rate_ratio', 
                              f'{parameter}_cumulative_mean_rate_ratio', f'{parameter}_cumulative_norm_ratio',
                              f'{parameter}_rate_old_ratio', f'{parameter}_rate_old_time_corrected_ratio']
        for well in well_names:
            d[antimycotic][file][well][parameter] = np.array([])
            for derived_par in derived_parameters:
                d[antimycotic][file][well][derived_par] = np.array([])
        
        # Create dataframe per antimycotic/parameter: filter columns that contain letter of current antimycotic (e.g. 'A', 'B', ...)
        # Note: error may occur if other column names would also contain one of these capitalized letters
        df_parameter_antimycotic = df_parameter.filter(regex=antimycotic)
                
        # Define column corresponding to MIC:
        MIC_column = MIC_dict[file][antimycotic]
    
        if file in training_files:
            # Iterate through dataframe and define MIC and other columns
            for name, content in df_parameter_antimycotic.items():
                logging.debug(f'--- --- --- Extracting well {name}')
                if name == MIC_column:
                    d[antimycotic][file]['MIC'][parameter] = content.values # .values returns a numpy array
                    extract_rate_cum(d, antimycotic, file, 'MIC', parameter, content.values, window)  
                    calculate_ratios(d, antimycotic, file, 'MIC', parameter)
                else: 
                    difference = int(MIC_column[1:]) - int(name[1:])
                    d[antimycotic][file][f'MIC{difference}'][parameter] = content.values
                    extract_rate_cum(d, antimycotic, file, f'MIC{difference}', parameter, content.values, window)
                    calculate_ratios(d, antimycotic, file, f'MIC{difference}', parameter)
                    if int(name[1:]) == 12:
                        d[antimycotic][file]['GC'][parameter] = content.values
                        extract_rate_cum(d, antimycotic, file, f'GC', parameter, content.values, window)
        else:
            # In case you want to extract data without reference to MIC
            for name, content in df_parameter_antimycotic.items():
                d[antimycotic][file][name[1:]][parameter] = content.values # .values returns a numpy array
                extract_rate_cum(d, antimycotic, file, name[1:], parameter, content.values, window)  
                calculate_ratios(d, antimycotic, file, name[1:], parameter)  

        
        if df_parameter_antimycotic.empty: # This is the case when no (olorofim ('F')) data is available
            for i in well_names: # Replace empty dictionaries as values with empty np.array
                d[antimycotic][file][i][parameter] = np.array([])
                for par in derived_parameters:
                    d[antimycotic][file][i][par] = np.array([]) # use empty np.array([]) instead of empty pd.DataFrame: plotting uses empty() function

def extract_derive_parameters(rootdir, filelist, parameterlist, antimycotics, window, ref_wells, training_files, MIC_dict):
    logging.debug(f'rootdir = {rootdir} \nfilelist = {filelist} \n parameterlist = {parameterlist} \n antimycotics = {antimycotics} \n window = {window} \n ref_wells = {ref_wells} \n training_files = {training_files} \n MIC_dict = {MIC_dict}')  

    # Initialize d with antimycotics at highest level: it will contain all data extracted from files per antimycotic + derived parameters
    d = {}
    for antimycotic in antimycotics:
        d[antimycotic] = {}
        
    for index, file in enumerate(filelist):
        logging.info(f'Reading {file} ({index + 1}/{len(filelist)})')
        df = pd.read_excel(os.path.join(rootdir, file), sheet_name=1)

        # Determine how wells will be named: relative to MIC (['MIC-2', 'MIC-1', 'MIC']) or absolute well numbers [1, 2, 3, 4, …]])
        if file in training_files:
            well_names = ref_wells
        else:
            well_names = [str(name) for name in [*range(13)]] # yields: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            # Why include 0? Because when ratios with leftside well are calculated, you need to initialize well 0 as containing no values

        # dictionary d contains raw data per file, well, parameter: 
        # [antimycotic] > [file] > [well] > [parameter] > data-arrays
        for antimycotic in antimycotics.keys():
            d[antimycotic][file] = {}
            for well in well_names:
                d[antimycotic][file][well] = {}
                
        # At highest level add per file an entry to d that will contain all time values for that file
        d[file] = np.array([])

        for parameter in parameterlist:       
            # Create entry for each parameter per file in d for each antimycotic
            logging.info(f'--- Extracting parameter {parameter}')
            
            extract_parameter(d, df, parameter, file, antimycotics, window, well_names, MIC_dict, training_files)
    return d
