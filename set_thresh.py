#!/home/kdewaele/.conda/envs/ocelloscope_env/bin/python 
#!/Users/kdewaele/miniconda3/envs/ocelloscope/bin/python

import argparse
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
import ast
import glob
from oCelloscope_tools import dictionaries
from oCelloscope_tools import setup
from oCelloscope_tools import data_extraction
from oCelloscope_tools import plot
from oCelloscope_tools import predict
from oCelloscope_tools import evaluate
import subprocess
import sys
import json
import random
import re
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')

home = '/scratch/leuven/348/vsc34807/oCelloscope' # Default folder for storing pipeline output 
home = '/home/kdewaele/ocelloscope/output' 
# copy_path used in evaluate.get_best_parameter() function for copying summary data for all runs to a single folder
copy_path = '/scratch/leuven/348/vsc34807/summary'
copy_path = '/home/kdewaele/ocelloscope/output/summary'
scripts_dir = '/home/kdewaele/ocelloscope/scripts'

def cleanup_temporary_pkl_files(output_dir, timepoint=None, antimycotic=None, keep_performance_files=True, recursive=False):
    """
    Remove temporary pkl files to save disk space.
    
    Args:
        output_dir: Directory to clean up
        timepoint: Specific timepoint to clean up (if None, cleans all)
        antimycotic: Specific antimycotic to clean up (if None, cleans all)
        keep_performance_files: Whether to preserve performance_*.pkl files (default: True)
        recursive: Whether to clean up subdirectories recursively
    """
    if not os.path.exists(output_dir):
        return
    
    files_removed = 0
    total_size_removed = 0
    
    # If recursive, clean up subdirectories first
    if recursive:
        try:
            for item in os.listdir(output_dir):
                item_path = os.path.join(output_dir, item)
                if os.path.isdir(item_path) and not item.startswith('.'):
                    # Recursively clean subdirectories
                    cleanup_temporary_pkl_files(item_path, timepoint, antimycotic, keep_performance_files, recursive=False)
        except Exception as e:
            logging.warning(f'Could not list directory {output_dir} for recursive cleanup: {e}')
    
    # Clean up threshold_dict files
    pattern = os.path.join(output_dir, "*_threshold_dict.pkl")
    for file_path in glob.glob(pattern):
        try:
            file_size = os.path.getsize(file_path)
            os.remove(file_path)
            files_removed += 1
            total_size_removed += file_size
            logging.debug(f'Removed threshold_dict file: {os.path.basename(file_path)}')
        except Exception as e:
            logging.warning(f'Could not remove {file_path}: {e}')
    
    # Clean up di_*.pkl files (use os.listdir to avoid glob pattern issues with special characters)
    try:
        all_files = os.listdir(output_dir)
        for filename in all_files:
            if filename.startswith('di_') and filename.endswith('.pkl'):
                # Check if this file matches our criteria
                should_remove = False
                
                if timepoint and antimycotic:
                    # Clean specific timepoint/antimycotic
                    if f"di_{antimycotic}_timepoint_{timepoint}_" in filename:
                        should_remove = True
                elif timepoint:
                    # Clean specific timepoint for all antimycotics
                    if f"_timepoint_{timepoint}_" in filename:
                        should_remove = True
                elif antimycotic:
                    # Clean specific antimycotic for all timepoints
                    if f"di_{antimycotic}_timepoint_" in filename:
                        should_remove = True
                else:
                    # Clean all di files
                    should_remove = True
                
                if should_remove:
                    file_path = os.path.join(output_dir, filename)
                    try:
                        file_size = os.path.getsize(file_path)
                        os.remove(file_path)
                        files_removed += 1
                        total_size_removed += file_size
                        logging.debug(f'Removed di file: {filename}')
                    except Exception as e:
                        logging.warning(f'Could not remove {file_path}: {e}')
    except Exception as e:
        logging.warning(f'Could not list directory {output_dir}: {e}')
    
    # Clean up performance_dict files (also temporary)
    pattern = os.path.join(output_dir, "performance_dict_*.pkl")
    for file_path in glob.glob(pattern):
        try:
            file_size = os.path.getsize(file_path)
            os.remove(file_path)
            files_removed += 1
            total_size_removed += file_size
            logging.debug(f'Removed performance_dict file: {os.path.basename(file_path)}')
        except Exception as e:
            logging.warning(f'Could not remove {file_path}: {e}')
    
    # Clean up all_par files (also temporary)
    pattern = os.path.join(output_dir, "all_par_*.pkl")
    for file_path in glob.glob(pattern):
        try:
            file_size = os.path.getsize(file_path)
            os.remove(file_path)
            files_removed += 1
            total_size_removed += file_size
            logging.debug(f'Removed all_par file: {os.path.basename(file_path)}')
        except Exception as e:
            logging.warning(f'Could not remove {file_path}: {e}')
    
    if files_removed > 0:
        size_mb = total_size_removed / (1024 * 1024)
        logging.info(f'Cleanup: Removed {files_removed} temporary pkl files, freed {size_mb:.1f} MB')
    else:
        logging.debug('Cleanup: No temporary pkl files found to remove')

def is_single_configuration(parameters_to_predict, thresholds, timepoint, antimycotic, sec_bias_corr):
    """
    Detect if we're testing a single configuration (1 parameter, 1 threshold, 1 bias).
    In this case, optimization is unnecessary - we just need performance calculation.
    """
    # Check if we have exactly one parameter
    if len(parameters_to_predict) != 1:
        return False
    
    parameter = parameters_to_predict[0]
    
    # Check if this parameter has exactly one threshold
    if len(thresholds[timepoint][antimycotic][parameter]) != 1:
        return False
    
    # Check if we're testing a single bias (not a range)
    if isinstance(sec_bias_corr, list):
        return False
    
    return True

def calculate_single_configuration_performance(timepoint, antimycotic, parameter, threshold, di, files_with_data, files_with_breakpoint, bias, args_dict, session_time, output_dir, di_path):
    """
    Direct performance calculation for single configuration scenarios.
    Bypasses optimization overhead when testing only one threshold/parameter/bias combination.
    """
    logging.info(f'--- --- Single configuration detected: calculating performance directly for {parameter} with threshold {threshold}')
    
    # Initialize the di structure that get_agreement_parameters_ext expects
    # This mimics what get_best_thresholds_bis does but for a single threshold
    if antimycotic not in di:
        di[antimycotic] = {}
    if timepoint not in di[antimycotic]:
        di[antimycotic][timepoint] = {}
    if parameter not in di[antimycotic][timepoint]:
        di[antimycotic][timepoint][parameter] = {}
    
    # Initialize the performance metrics list [total_MIC_distance, VME_count, ME_count, ATU_VME_count, ATU_ME_count]
    di[antimycotic][timepoint][parameter][threshold] = [0, 0, 0, 0, 0]
    
    # We need to calculate the actual performance metrics by iterating through files
    # This is what normally happens in get_best_thresholds_bis but we do it directly here
    for file in files_with_data:
        try:
            # Get the MIC distance and categorical errors for this file
            # Based on predict.py, the data is stored under 'MIC'
            distance = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['MIC_distance']
            VME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['VME']
            ME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ME']
            ATU_VME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_VME']
            ATU_ME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_ME']
            
            # Sum up the metrics across all files (matching evaluate.py line 169-175)
            di[antimycotic][timepoint][parameter][threshold][0] += abs(distance) if not np.isnan(distance) else 0  # total_MIC_distance
            if file in files_with_breakpoint:
                # if no breakpoint, these values remain 0
                di[antimycotic][timepoint][parameter][threshold][1] += VME if not np.isnan(VME) else 0       # VME_count
                di[antimycotic][timepoint][parameter][threshold][2] += ME if not np.isnan(ME) else 0        # ME_count  
                di[antimycotic][timepoint][parameter][threshold][3] += ATU_VME if not np.isnan(ATU_VME) else 0   # ATU_VME_count
                di[antimycotic][timepoint][parameter][threshold][4] += ATU_ME if not np.isnan(ATU_ME) else 0    # ATU_ME_count
                
        except KeyError as e:
            logging.debug(f'No prediction data for file {file}: {e}')
            # If no prediction data, treat as a major error (distance = 12 is used for failed predictions)
            di[antimycotic][timepoint][parameter][threshold][0] += 12
            if file in files_with_breakpoint:
                di[antimycotic][timepoint][parameter][threshold][2] += 1  # Count as ME
    
    # Now calculate performance metrics directly using the existing evaluation function
    agreement_dict = evaluate.get_agreement_parameters_ext(antimycotic, timepoint, parameter, threshold, di, files_with_data, files_with_breakpoint)
    
    # Create the performance_per_timepoint dictionary in the same format as the optimization functions
    performance_per_timepoint = {}
    
    # Copy all argument information (same as get_best_parameter does)
    for key, value in args_dict.items():
        performance_per_timepoint[key] = value
    
    # Add performance metrics
    performance_per_timepoint['time'] = timepoint
    performance_per_timepoint['parameter'] = parameter
    performance_per_timepoint['threshold'] = threshold
    performance_per_timepoint['session_time'] = session_time
    performance_per_timepoint['di_path'] = di_path
    performance_per_timepoint['output_dir'] = output_dir
    
    # Add all the agreement metrics
    for key, value in agreement_dict.items():
        performance_per_timepoint[key] = value
    
    # Map 'distance_per_file' to 'dist_per_file' for consistency with optimization path
    # (get_agreement_parameters_ext returns 'distance_per_file', but rest of code expects 'dist_per_file')
    performance_per_timepoint['dist_per_file'] = agreement_dict.get('distance_per_file', {})
    
    # Calculate total_errors (sum of weighted categorical errors)
    performance_per_timepoint['total_errors'] = agreement_dict.get('min_errors_weighted', 0)
    
    logging.info(f'--- --- Single configuration performance: EA={agreement_dict.get("EA", "N/A"):.3f}, CA={agreement_dict.get("CA", "N/A"):.3f}, Total errors={performance_per_timepoint["total_errors"]}')
    
    return performance_per_timepoint

parser = argparse.ArgumentParser()
parser.add_argument('-g', type=str, required=True, help="Genus: Aspergillus (A) or Candida (C).")
parser.add_argument('-a', type=str, required=True, help="Antimycotic acccording to row letter on EUCAST broth microdilution plates layout.")
parser.add_argument('-d', type=str, required=False, help="Optional: path to dictionary d with extracted raw data and derived parameters. If not supplied data will be extracted with default window 6. If integer provided, create d with custom window size.")
parser.add_argument('-t', type=str, required=True, help="Provide start time, end time and interval in seconds, as a tuple enclosed by quotation marks.")
parser.add_argument('-p', type=str, required=False, help="Provide parameters to predict as a list enclosed by quotation marks. If none supplied, all parameters will be predicted.")
parser.add_argument('-s', type=str, required=False, help="The number of thresholds (steps) between minimum and maximum values at a certain timepoint. If '1' supplied, EUCAST approach will be used (threshold of percent inhibition compared to growth control using default dictionaries. If float supplied from 0.0 to 1.0, this will be selected as the proportion inhibition compared to growth control set as threshold per file. If string representation of list of [start, end, step] is provided, a range of proportional values between start and end will be created.")
parser.add_argument('-i', type=str, required=True, help="Input directory containing source Excel files.")
parser.add_argument('-b', type=str, required=False, help="Number of files used for bagging.")
parser.add_argument('-o', type=str, required=False, help="Output path for creation of session folder(s). If none supplied, home folder will be used." )
parser.add_argument('-l', type=str, required=False, help="Desired logging level.")
parser.add_argument('-r', type=str, required=False, help="Files in input directory not to be included in analysis. Supply as list enclosed by quotation marks.")
parser.add_argument('-m', type=str, required=False, help="Provide integer to move detected pattern with specific bias correction, or provide string representation of list with one integer, to create range of biases for -i to +i.")
parser.add_argument('-rerun', type=str, required=False, help="Flag used for recursive execution of scripts")
parser.add_argument('-c', type=str, required=False, help="Criterium for threshold optimization: minimise total_MIC_distance for essential agreement (EA) or minimise min_errors_weighted for categorical agreement (CA), or optimise EA or CA with a target value so that selected thresholds will meet at least supplied threshold, in order to avoid CA-EA 'trade-off'. Supply as string representation of list and tuple, e.g. \"[('EA', 'total_MIC_distance')]\". When supplying threshold for minimum target value (e.g. >= 0.9) for, for example, CA, tool will secondarily select threshold with best EA. Problem with absolute CA optimization approach is that these will always primarily select maximum categorical agreement (say 1), and only select thresholds with best EA among these (while EA of these may be uniformly very poor if CA is maximised at expense of EA). Therefore, this option balances CA en EA. For absolute maximisation of EA, there's usually a better correlation with maximal CA - the opposite is not true.")
parser.add_argument('-n', type=str, required=False, help="An informative string to be added to the session directory name.")
parser.add_argument('-x', type=str, required=False, help="Pattern to be detected. If not supplied, use default pattern [False, True, True, True]")
parser.add_argument('--min_CA', type=str, required=False, help="When selecting best threshold, tool will select thresholds with at least supplied categorical agreement. Among these thresholds, tool will secondarily select threshold with best EA. Difference with threshold optimization criterion [('CA', 'min_errors_weighted')] is that latter will always primarily select maximum categorical agreement (say 1), and only select thresholds with best EA among these (while EA of these may be uniformly very poor if CA is maximised at expense of EA). Therefore, this option balances CA en EA.")
parser.add_argument('--min_EA', type=str, required=False, help="See --min_CA option: here thresholds will be preselected that at least meet supplied essential agreement.")


# Argument parsing
args=parser.parse_args()
args_overview = pprint.pformat(vars(args))
arguments = vars(args)

# Time
start, end, step  = ast.literal_eval(args.t)
timepoints = [*range(start, end, step)]
if (start + step) > end: # Run for one timepoint, show only start label, for range, show start and end_label in '{start}{-end_label}' format
    end_label = ''
else:
    end_label = f'-{end}'

# Optimization criterium
if args.c:
    criteria = ast.literal_eval(args.c)
    criteria_label = criteria[0][1]
else:
    criteria = [('EA', 'total_MIC_distance'), ('CA', 'min_errors_weighted')] 
    criteria_label = f'{criteria[0][1]}-{criteria[1][1]}'

# Create label for naming output directory and log file
label = f"{str(args.g)}_{str(args.a)}_{criteria_label}_time_{start}{end_label}_bias_{str(args.m)}_{str(args.n)}"

# Set-up output directory
output_dir, session_time  = setup.output_setup(str(args.o), home, label)

# Set-up log file handler: will be saved in home directory - add output_dir if it must be in session folder
handler = logging.FileHandler(os.path.join(home, output_dir, f"{label}.log"))
logger = logging.getLogger()
logger.addHandler(handler)

# Logging level
if args.l:
    logging.getLogger().setLevel(level=args.l.upper())
    handler.setLevel(level=args.l.upper())

# Rerun
if args.rerun:
    logging.info(f'\n\n<<< <<< <<< RERUN ({args.rerun}) >>> >>> >>>\n')

# Print arguments
logging.info(f'arguments: \n{args_overview}')
logging.info(f'Timepoints: start {start}, end {end}, step {step}: {timepoints}')
logging.info(f'criterium/a = {criteria}')

# Number of thresholds between min-max
eucast = False
if args.s:
    steps = args.s
    if not isinstance(ast.literal_eval(steps), list):
    # In case a single value is provided: if string repr of float > will be used as a prop thresh, if int will be used as an abs number of steps (except if 1 - then default EUCAST relative thresh are used): in case of integer threshold_dict is left empty
        # Use regex to differentiate integers (1, 2 ...) from float (0.0, 1.0, ...) 
        if re.match("^[+-]?\d+?\.\d+?$", steps) is not None:
            logging.info(f'Steps is {steps}. EUCAST algorithm applied with custom threshold_dict.')
            eucast = True
            prop = float(steps)
            threshold_dict = {'A': [prop], 'B': [prop], 'C': [prop], 'D': [prop], 'E': [prop], 'F': [prop], 'G': [prop], 'H': [prop]}
        elif re.match("^[+-]?\d+?$", steps) is not None:
            steps = int(steps)
            logging.info(f'Steps is {steps}.')
            if steps == 1:
                eucast = True
                logging.info(f'Steps is {steps}. EUCAST algorithm applied with routine dictionaries.')
            threshold_dict = {}
        else:
            logging.info(f'Invalid args.s argument provided. Quitting.')
            quit()
    else: # in case a list is provided, a number of relative thresholds will be set; threshold dict will contain for present antimycotic this range of thresholds [0.09, 0.1, 0.11, 0.12 ...]
        range_prop = ast.literal_eval(steps)
        logging.debug(f'range_prop {range_prop}, type {type(range_prop)}')
        try:
            threshold_dict = {}
            threshold_dict[str(args.a)] = [*np.round(np.arange(range_prop[0], range_prop[1], range_prop[2]), 2)] 
        except Exception as e:
            logging.warning(f'Creating range of proportion values failed: {e}. Quitting.')
            quit()
        logging.debug(f'Args.s is list: creating range of proportional values: {threshold_dict.items()}')
        eucast = True
else:
    steps = 10
    logging.info(f'No steps argument provided; using default {steps} steps')
    threshold_dict = {}

# Pattern
if args.x:
    pattern = ast.literal_eval(args.x)
else:
    pattern = [False, True, True, True]
logging.info(f'Detect pattern: {pattern}')

# Genus
if str(args.g) == 'A':
    antimycotics_dict = dictionaries.antimycotics_dict_A
    breakpoint_dict = dictionaries.breakpoint_dict_A
    MIC_dict = dictionaries.MIC_dict_A
    species_dict = dictionaries.species_dict_A
    if eucast and steps == 1:
        threshold_dict = dictionaries.eucast_A
        logging.debug(f'EUCAST threshold dict: {threshold_dict}')
    logging.info(f'Aspergillus dictionaries selected')
else:
    antimycotics_dict = dictionaries.antimycotics_dict_C
    breakpoint_dict = dictionaries.breakpoint_dict_C
    MIC_dict = dictionaries.MIC_dict_C
    species_dict = dictionaries.species_dict_C
    if eucast and steps == 1:
        threshold_dict = dictionaries.eucast_C
        logging.debug(f'EUCAST threshold dict: {threshold_dict}')
    logging.info(f'Candida dictionaries selected')
# Only when using proportional thresholds should threshold_dict not be empty: it contains a proportion or list of proportions (float between 0.0 to 1)
# It better should be called proportional_thresholds_dict, to differentiate from thresholds{}, that has structure thresholds[timepoint][antimycotic][parameter] = []

# Input files
input_dir = setup.input_setup(str(args.i))
filelist = os.listdir(input_dir)
logging.debug(f'filelist: {filelist}')
logging.info(f'{len(filelist)} file(s) found')
if args.r:
    files_to_remove = ast.literal_eval(args.r)
    logging.info(f'{len(files_to_remove)} file(s) left out')
    [filelist.remove(file_to_remove) for file_to_remove in files_to_remove if file_to_remove in filelist]
training_files = filelist 
files_pDST = training_files + filelist
# If path to extracted data dictionary d provided, load, if not provided or unable to read, create d
# If integer provided, create d with custom window size
create_d = True
d = {}
if not args.d:
    logging.info(f'No dictionary with extracted data provided - will be created with default window size 6')
    window = 6
elif os.path.isfile(str(args.d)):
    logging.info(f'Dictionary with extracted data provided at {str(args.d)}. Loading ...')
    with open(str(args.d), 'rb') as f:
        try: 
            d = pickle.load(f)
            create_d = False
            window = 'user'
        except: 
            logging.error('Could not read dictionary at {str(args.d)} - will be created with default window size 6')
            window = 6
else:
    try:
        window = int(args.d)
        logging.info(f'Extracting data, using window size {window}')
    except ValueError:
        logging.warning(f'No valid path or window value provided - data extraction dictionary will be created with default window size.')
        window = 6

# Parameters for data extraction: these must encompass those that will be predicted
parameters_to_extract = ['BCA', 'BCANormalized', 'SESAfungi', 'SESAfungiNormalized', 'TA', 'TANormalized']
logging.info(f'Parameters to be extracted: {parameters_to_extract}')

# Plate layout: init_wells is range of wells per antimycotic used to construct the name of wells relative to MIC (MIC-2, MIC-3, etc): put in ref_wells list
wells_per_row = 12
init_wells = [*range(-(wells_per_row + 1), wells_per_row + 2)] # yields: [-13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
# we need well +/-13 compared to MIC for calculating ratios of neighboring wells: in case rMIC is limit value, and pMIC +/-12, then neighboring wells may become +/-13 (although obviously there's no data for such wells)
init_wells.remove(0) # to avoid 'MIC0'
init_wells.append('') # instead get 'MIC'
ref_wells = [f'MIC{i}' for i in init_wells]
ref_wells.reverse() # to put 'negative' wells at right side (only relevant for plotting), cfr. layout EUCAST plates
ref_wells.append('GC') # GC = Growth control
ref_wells.insert(wells_per_row + 1, ref_wells.pop(0)) # to put 'MIC' in middle of the list

if create_d:
    d = data_extraction.extract_derive_parameters(input_dir, filelist, parameters_to_extract, antimycotics_dict, window, ref_wells, training_files, MIC_dict)
    d_path = os.path.join(home, f'genus_{str(args.g)}_d_window_{window}.pkl')
    # add d_path to current args, to prevent repeated data extraction for nested reruns of pipeline (that take a copy of parent args.d)
    args.d = d_path
    with open(d_path, 'wb') as f:
        pickle.dump(d, f)
        logging.info(f'Succesfully saved data extraction dictionary at {d_path}')

# Systematic error correction
sec_bias_corr = False # Secondary bias corr: allows recursive execution of pipeline with (a range of) bias correction factor(s) 
if not args.m:
    bias = 0
elif str(args.m) == 'average':
    # average bias will be calculated and pipeline rerun to correct this 
    sec_bias_corr = 'average' # variable not used for now, bias correction always on
    logging.info(f'Bias correction for average systematic error.')
else:
    args_m = ast.literal_eval(args.m)
    if isinstance(args_m, list):
        # sec_bias_corr will be a range of biases for nested execution of pipeline for each of these biases
        bias = 0
        sec_bias_corr = [i for i in range(-args_m[0], args_m[0]+1)]
        logging.info(f'Multiple bias correction requested for range {sec_bias_corr}')
    elif isinstance(args_m, int):
        # Single bias integer supplied, secondary bias corr remains False
        bias = args_m
        logging.info(f'Bias is {bias}')
    else:
        logging.warning(f'Bias argument not recognized: args.m {args.m}, args_m {args_m}, type args_m {type(args_m)}')

# Antimycotics to predict: code optimized for one at a time
antimycotics_to_predict = [str(args.a)]

# Parameters to predict
if args.p:
    parameters_to_predict = ast.literal_eval(args.p)
else:
    parameters_to_predict = dictionaries.all_parameters
if eucast: # Ratio parameters cannot be used for proportional prediction: exclude them
    parameters_to_predict = [para for para in parameters_to_predict if 'ratio' not in para]
logging.info(f'Predicting {len(parameters_to_predict)} derived parameter(s)')

# In case no pMIC data is available: use supplied threshold in threshold_dict to predict files
if args.b:
    bagging_args = ast.literal_eval(args.b)
    cycles = int(bagging_args[0])
    cycles_start = cycles
    bootstrap_size = int(bagging_args[1])
    threshold_dict_path = str(bagging_args[2])
    #threshold_dict = # path to dictionariy with training thresholds dictionaries.threshold_dict = {}
    training = True
    if os.path.isfile(threshold_dict_path):
        training = False
        logging.info(f'{cycles}  bagging cycle(s) remain. Bootstrap size {bootstrap_size}, threshold_dict path provided: {threshold_dict_path}.')
        with open(threshold_dict_path, 'rb') as f:
            logging.info(f'--- Reading threshold dictionary ...')
            threshold_dict = pickle.load(f)
    else:
        training = True
        logging.info(f'{cycles} bagging cycle(s) remain. Bootstrap size {bootstrap_size}, no valid threshold_dict path provided. Setting thresholds.')
else:
    training = True
logging.info(f'Training? {training}')

def get_denominators(antimycotic):
    no_data = []
    no_breakpoint = []
    for file, items in d[antimycotic].items():
        # get species for this file:
        species = species_dict[file]
        logging.debug(f'--- --- --- Species for file {file} is {species}')
        if (file in filelist) and ((not items['MIC']['BCA'].any()) and (not items['MIC1']['BCA'].any()) and (not items['MIC-1']['BCA'].any())): # in case for this file no data is available for the present antimycotic
            no_data.append(file)
            logging.debug(f'--- --- --- --- No data for file {file}')
        try:
            if (np.isnan(breakpoint_dict[f'{antimycotic}_{species}'][0])): # no breakpoint found
                no_breakpoint.append(file)
                logging.debug(f'--- --- --- --- No breakpoint found for file {file} with species {species}')
        except KeyError:
                no_breakpoint.append(file)
                logging.debug(f'--- --- --- --- No entry for species {species} (file {file}) in breakpoint_dict')
            
    data = list(set(filelist) - set(no_data))
    breakpoint = list(set(data) - set(no_breakpoint))
    logging.info(f'--- --- Files with available data ({len(data)})')
    logging.info(f'--- --- No data available for {no_data}')
    logging.info(f'--- --- Files with available data and breakpoint ({len(breakpoint)})')
    return data, breakpoint 


def finish_rerun(rerun_label, performance_per_timepoint):
    beautiful_dict = pprint.pformat(performance_per_timepoint)
    logging.info(f'--- --- --- Performance_per_timepoint after rerun: \n\n {beautiful_dict}\n')
    logging.warning(f'\n\n<<< <<< <<< RERUN FINISHED ({rerun_label}) >>> >>> >>> Passing dict to parent script\n')
    sys.stdout.write(json.dumps(performance_per_timepoint))
    sys.stdout.flush()
    sys.exit()

def pipeline_rerun(**kwargs):
    arguments = kwargs.get('args')
    kwargs.pop('args')
    arg_dict = vars(arguments).copy()
    logging.debug(f'args: {arguments} \n arg_dict {arg_dict}')
    for arg_name, arg_value in kwargs.items():
        arg_dict[arg_name] = arg_value     

    # Convert args to a list of command line arguments
    arg_list = ["-{}={}".format(key, value) for key, value in arg_dict.items() if value != None]
    # Run nested rerun of pipeline specifically for this criterium and with correction factor
    process = subprocess.run([os.path.join(scripts_dir, "set_thresh.py")] + arg_list, stdout=subprocess.PIPE)

    performance_per_timepoint = json.loads(process.stdout)
    # Get correct output_dir and di dict for plot()
    output_dir_rerun = performance_per_timepoint['output_dir']
    with open(performance_per_timepoint['di_path'], 'rb') as f:
        di_plot = pickle.load(f)
    logging.info(f'--- --- --- Returning performance_per_timepoint that results from rerun \nModified di at {performance_per_timepoint["di_path"]}\nOutput dir is {performance_per_timepoint["output_dir"]}.\n')

    return performance_per_timepoint, output_dir_rerun, di_plot


performance_chron = []

for timepoint in timepoints:
    logging.info(f'Timepoint {timepoint}')
    di = {}
    thresholds = {}
    for antimycotic in antimycotics_to_predict:
        logging.info(f'--- Antimycotic: {antimycotic}')
        # Differentiate files with lacking data for certain antimycotics (e.g. Olorofim) and files for species that do not have breakpoints
        files_with_data, files_with_breakpoint = get_denominators(antimycotic)

        if args.b:
            files_with_data = random.sample(files_with_data, min(bootstrap_size, len(files_with_data)))

        thresholds, parameters_to_predict, bias_from_dict = predict.define_thresholds(timepoint, antimycotic, parameters_to_predict, steps, filelist, d, thresholds, eucast, threshold_dict, training)
        # predict() assigns True/False to data points according to whether threshold was passed (categorise()) and predicts MIC depending on pattern of growth in neighboring wells (predict_MIC())
        # it also initializes thresholds according to steps between min/max value - this function could be separated from predict() function 

        multiple_unpenalized = True # In case multiple pattern matches, the rightmost is chosen. If multiple_unpenalized False, then in case of multiple matches this will be penalized as a ME.
        # This function is executed for all files with any data for that antimycotic, regardless of whether there's data for present timepoint - among these files the function identifies files lacking data for specific timepoint, and returns these as a set
        no_data_for_timepoint, thresh_dict = predict.predict(timepoint, antimycotic, files_with_data, parameters_to_predict, ref_wells, steps, thresholds, di, MIC_dict, breakpoint_dict, threshold_dict, d, species_dict, output_dir, bias, eucast, pattern, multiple_unpenalized)

        # Now modify files_with_data and files_with_breakpoint as only to include files that have data at this timepoint
        files_with_data = list(set(files_with_data).difference(no_data_for_timepoint))
        files_with_breakpoint = list(set(files_with_breakpoint).difference(no_data_for_timepoint)) 

        # get_best_thresholds_bis() finds which thresholds PER PARAMETER results in best EA (and secondarily best CA) and best CA (with secondarily best EA) (or optimisation criteria can also apply - see get_best_thresholds_bis())
        # CA can only be calculated for files with species that have a breakpoint
        # July 17: added as argument 'criteria' variable to enable custom CA/EA optimisation threshold selection, for example, (CA, 0.90), means at least CA threshold of 90%
        
        # Check if this is a single configuration scenario - bypass optimization if so
        if is_single_configuration(parameters_to_predict, thresholds, timepoint, antimycotic, sec_bias_corr):
            logging.info(f'--- --- Single configuration detected - bypassing optimization')
            parameter = parameters_to_predict[0]
            threshold = thresholds[timepoint][antimycotic][parameter][0]
            
            # Save di (same structure as optimization path)
            criteria_label_for_timepoint = criteria[0][1]
            di_path = os.path.join(output_dir, f'di_{antimycotic}_timepoint_{timepoint}_steps_{steps}_bias_{str(bias)}_{criteria_label_for_timepoint}.pkl')
            with open(di_path, 'wb') as f:
                pickle.dump(di, f)
            logging.info(f'--- --- Succesfully saved di dictionary at {di_path}')
            
            # Calculate performance directly without optimization
            performance_per_timepoint = calculate_single_configuration_performance(
                timepoint, antimycotic, parameter, threshold, di, files_with_data, 
                files_with_breakpoint, bias, arguments, session_time, output_dir, di_path
            )
            
            # Set criteria info to match expected structure
            criteria_for_timepoint = [(criteria[0][0], criteria[0][1])]
            
        else:
            # Standard optimization path
            skip_flag = False
            modified_target, skip_flag = evaluate.get_best_thresholds_bis(timepoint, antimycotic, parameters_to_predict, di, thresholds, files_with_data, files_with_breakpoint, criteria[0], output_dir, skip_flag)
            logging.info(f'modified_target: {modified_target}')
            if isinstance(modified_target, (int, float)) and np.isnan(modified_target):
                logging.warning(f'No targets found that meet specified criterion - skipping timepoint.')
                break
            elif modified_target != criteria[0][1]:
                logging.warning(f'Note: criterium ({criteria[0][1]}) was modified: {modified_target}, moving on with new criterium')
                criteria_for_timepoint = [(criteria[0][0], modified_target)]
                criteria_label_for_timepoint = criteria[0][1]
            else: # if modified target is the same as original target
                criteria_for_timepoint = [(criteria[0][0], criteria[0][1])]
                criteria_label_for_timepoint = criteria[0][1]

            # Save di
            di_path = os.path.join(output_dir, f'di_{antimycotic}_timepoint_{timepoint}_steps_{steps}_bias_{str(bias)}_{criteria_label_for_timepoint}.pkl')
            with open(di_path, 'wb') as f:
                pickle.dump(di, f)
            logging.info(f'--- --- Succesfully saved di dictionary at {di_path}')

            for EA_CA, criterium in criteria_for_timepoint:
                logging.info(f'--- --- EA_CA is {EA_CA}, criterium is {criterium}')
                logging.info(f'--- --- Finding best parameter-threshold pair for criterium: {criterium}')
                # Performance_per_timepoint also contains di_path and output_dir, threshold_dict path contains the path to best parameter and threshold for time and antimycotic
                performance_per_timepoint, threshold_dict_path = evaluate.get_best_parameter(timepoint, antimycotic, parameters_to_predict, EA_CA, criterium, di, output_dir, di_path, arguments, session_time, bias)#, steps, window, pattern)
        
        # Common post-processing for both single configuration and optimization paths
        # Add files for which no data was available for this timepoint
        performance_per_timepoint['no_data'] = list(no_data_for_timepoint)
        # Add training status
        performance_per_timepoint['training'] = training

        # The path to the new di and output_dir are entries in the dictionary
        if args.rerun == 'bias': # 'bias' flag set in case of rerun of pipeline for bias correction - finish_rerun() to prevent endless bias correction
            # finish_rerun writes performance_per_timepoint to parent script and does sys.flush() and sys.exit()
            finish_rerun(args.rerun, performance_per_timepoint)
        elif args.rerun == 'bagging' and not training: # training False here because threshold_dict path provided as pipeline argument (args.b)
            # When this clause? For 'test' runs while bagging: bagging True and training False: prevent bias correction and nested bagging 
            logging.info(f'--- --- --- Non-training run: no bias correction.')
            finish_rerun(args.rerun, performance_per_timepoint)
        #elif criterium == "total_MIC_distance": # Only in case of minimizing total_MIC_distance, do bias correction - not for min_weighted_errors
        else: # removed requirement of min_distance criterium: what is rationale behind this? Also, this would not work when using "CA_0.90" optimization.
        # If setting threshold, detect systematic error: average MIC distances across files of best threshold: if different from 0, skip plotting, and rerun pipeline with bias correction
            dist_per_file = performance_per_timepoint['dist_per_file']
            # calculate median
            average = int(round(statistics.median(dist_per_file.values()), 0))
            #total = 0
            #for file, distance in dist_per_file.items():
                #total += distance
            #average = int(round(total/len(dist_per_file.keys()), 0))
            beautiful_dict = pprint.pformat(performance_per_timepoint)
            logging.info(f'--- --- --- Bias of {average} log2 dilutions. Performance_per_timepoint before rerun: \n\n {beautiful_dict}\n')

            # Document pre-bias correction performance
            performance_per_timepoint['phase'] = 'pre_bias_corr'
            performance_chron.append(performance_per_timepoint)

            
            if isinstance(sec_bias_corr, list):
                # Initialize based on whether we want to minimize or maximize the criterion
                if criterium in ['EA', 'CA']:
                    # For EA and CA, we want the HIGHEST values (best performance)
                    best_total_errors = -np.inf
                    comparison_func = lambda new_val, best_val: new_val > best_val
                    logging.info(f'Optimizing for HIGHEST {criterium} values')
                else:
                    # For min_errors_weighted, total_MIC_distance, etc., we want the LOWEST values
                    best_total_errors = np.inf
                    comparison_func = lambda new_val, best_val: new_val < best_val
                    logging.info(f'Optimizing for LOWEST {criterium} values')
                
                for shift in sec_bias_corr:
                    logging.info(f'\n\n shift {shift} in {sec_bias_corr}')
                    performance_per_timepoint_shift, output_dir_rerun_shift, di_plot_shift = pipeline_rerun(m=shift, rerun='bias', c=[(criteria[0][0], criteria[0][1])], t=f'({timepoint}, {timepoint + 1}, 2)', o=output_dir, args=args)
                    performance_per_timepoint_shift['phase'] = f'shift_{shift}'
                    performance_chron.append(performance_per_timepoint_shift)
                    if comparison_func(performance_per_timepoint_shift['total_errors'], best_total_errors):
                        logging.info(f'New best bias correction: {shift} - total errors is {performance_per_timepoint_shift["total_errors"]}, previous {best_total_errors}')
                        best_total_errors = performance_per_timepoint_shift['total_errors']
                        performance_per_timepoint = performance_per_timepoint_shift.copy()
                        output_dir_rerun = output_dir_rerun_shift
                        di_plot = di_plot_shift
                performance_per_timepoint['best_shift_corr'] = 'best'
                logging.warning(f'added dict with best total errors {best_total_errors} and shift {performance_per_timepoint["phase"]} to performance_chron list')
                performance_chron.append(performance_per_timepoint) 

            # This clause enables automated bias dtection and bias correction, in practice empiric bias correction over range works better than one derived from average error
            elif average != 0:
                logging.info(f'Systematic bias detected of {average} log2 dilutions; plotting skipped and rerun pipeline with correction factor')
                # Rerun pipeline with 'bias' flag for this one timepoint and add to performance_chronological list
                # pipeline_rerun() @provides modified di (called 'di_plot') and output directory of the nested run, for the plotting function
                performance_per_timepoint, output_dir_rerun, di_plot = pipeline_rerun(m=average, rerun='bias', c=[(criteria[0][0], criteria[0][1])], t=f'({timepoint}, {timepoint + 1}, 2)', o=output_dir, args=args)
                bias_corr_cycle_number = 0
                while bias_corr_cycle_number <= 2:
                    average = int(round(statistics.median(performance_per_timepoint['dist_per_file'].values()), 0))
                    logging.info(f'Remaining systematic error: {average}')
                    if average != 0:
                        bias_corr_cycle_number += 1
                        logging.info(f'Rerun bias correction number {bias_corr_cycle_number}')
                        # increase bias correction: isolate sign of average and add bias_corr_cycle_number: -1 becomes -2, becomes -3 ...
                        average = int(math.copysign(1, average) * (abs(average) + bias_corr_cycle_number))
                        logging.info(f'New bias correction: {average}')
                        performance_per_timepoint, output_dir_rerun, di_plot = pipeline_rerun(m=average, rerun='bias', c=[(criteria[0][0], criteria[0][1])], t=f'({timepoint}, {timepoint + 1}, 2)', o=output_dir, args=args)
                    else:
                        break
                performance_per_timepoint['phase'] = 'bias_corrected_cycle_{bias_corr_cycle_number}'
                performance_chron.append(performance_per_timepoint)
            else: # In case of no bias document performance_per_timepoint
                di_plot = di
                output_dir_rerun = output_dir  # Use main output directory when no bias correction
                performance_per_timepoint['phase'] = 'no_bias_corr'
                performance_chron.append(performance_per_timepoint)

            if args.rerun == 'bagging': # To prevent endless recursive bagging (see below), terminate this training cycle here
                logging.warning(f'Terminating training cycle. Writing performance_per_timepoint to parent process.')
                finish_rerun(args.rerun, performance_per_timepoint)

            # For bagging, we need here to rerun pipeline n cycles times
            while args.b and cycles != 0:
                logging.warning(f'Cycle {cycles}/{cycles_start}. Now TEST RUN.')
                # rerun pipeline in non-training mode: this is flagged by providing path to threshold_dict, which was returned from evaluate.get_best_parameter()
                # What does training False mean? define_threshold() will define one parameter-threshold pair in thresholds dict (thresholds[timepoint][antimycotic][parameter]) 
                # The pipeline will terminate after evaluate.get_best_parameter() (which returns performance per timepoint dataframe) 
                training_data = files_with_data + files_to_remove # These files need to be removed in a test run: supply for -r argument
                logging.info(f'These files will be removed from test run: {training_data}')
                performance_per_timepoint, output_dir_rerun, di_plot = pipeline_rerun(args=args, b=f'({cycles}, {bootstrap_size}, "{str(threshold_dict_path)}")', rerun='bagging', c=[(EA_CA, criterium)], t=f'({timepoint}, {timepoint + 1}, 2)', o=output_dir, r=f"{training_data}")
                performance_per_timepoint['phase'] = 'no_bias_corr'
                performance_chron.append(performance_per_timepoint)
                # do a fresh training cycle with cycles-1 - no dict path provided, which triggers training True  
                cycles -= 1
                logging.warning(f'Starting cycle {cycles}/{cycles_start}. Now TRAINING RUN.')
                performance_per_timepoint, output_dir_rerun, di_plot = pipeline_rerun(args=args, b=f'({cycles}, {bootstrap_size}, {None})', rerun='bagging', c=[(EA_CA, criterium)], t=f'({timepoint}, {timepoint + 1}, 2)', o=output_dir)
                # Still to remove files used in training: f'\"{performance_per_timepoint['dist_per_file'].keys()}\"'
                performance_chron.append(performance_per_timepoint)                 

            # highlight has (here redundant) structure with par and threshold, but this is needed for plotting specific thresholds from within predict() 
            highlight = {} 
            for file, distance in performance_per_timepoint['dist_per_file'].items(): 
                highlight[file] = {} # highlight[file][parameter][threshold]
                highlight[file][performance_per_timepoint['parameter']] = {}
                highlight[file][performance_per_timepoint['parameter']][performance_per_timepoint['threshold']] = distance
                
            #plot.plot(output_path=output_dir, d=d, di=di_plot, params=[performance_per_timepoint['parameter']], am=[antimycotic], files=files_with_data, 
            #         timeframe=((20000, 90000)), criterium=criterium, thresh_dict = thresh_dict,# timepoint - 0.2 * timepoint, timepoint + 0.2 * timepoint
            #         well_legend=ref_wells, thresh_y=performance_per_timepoint['threshold'], timepoint=timepoint, highlight=highlight,
            #        label=performance_per_timepoint, y_lim = threshold_dict)
            

            beautiful_dict = pprint.pformat(performance_per_timepoint)
            logging.info(f'--- --- --- Best performance_per_timepoint dict: \n\n {beautiful_dict}\n')

            performance_chron_df = pd.DataFrame(performance_chron)
            performance_chron_df.to_pickle(os.path.join(output_dir, f'performance_{label}.pkl'))  
            performance_chron_df.to_pickle(os.path.join(copy_path, f'performance_{label}_{session_time}.pkl'))  
            logging.info(f'Saved performance dataframe.')

            # Clean up temporary files after bias correction cycles complete
            cleanup_temporary_pkl_files(output_dir_rerun, timepoint, antimycotic)

    # Clean up remaining temporary files after timepoint completes
    logging.info(f'Timepoint {timepoint} completed. Cleaning up remaining temporary files.')
    cleanup_temporary_pkl_files(output_dir, timepoint, recursive=True)

plot.plot_overview(performance_chron_df, *antimycotics_to_predict, len(timepoints), len(parameters_to_predict), steps, output_dir, label, copy_path)

# Final cleanup of any remaining temporary files
cleanup_temporary_pkl_files(output_dir, recursive=True)
logging.info('Pipeline completed. All temporary pkl files cleaned up.')




