import logging
import numpy as np
from oCelloscope_tools import plot
import os
import pickle

def determine_CA(file, distance, antimycotic, MIC_dict, breakpoint_dict, species_dict):
    logging.debug(f'--- --- --- Determining CA for {file} and {antimycotic}')
    failed_pred = []
    ME, VME, ATU_ME, ATU_VME = 0, 0, 0, 0

    if distance >= 13: # When pattern has been found multiple times or was not found, distance becomes 13 - (-1) (depending on len(ref_wells))
        failed_pred.append(file)
        logging.debug(f'--- --- --- --- Prediction failed for {file} and {antimycotic}')
        # In this case, we return 1 for VME: this adds counts failed prediction as a VME
        # You could also return 1, 0, 0, 0, to penalize failed prediction as VME
        return 0, 1, 0, 0
    else:
        try: 
            breakpoint_list = breakpoint_dict[f'{antimycotic}_{species_dict[file]}']
        except KeyError:
            logging.warning(f'--- --- --- --- No entry in breakpoint_dict for {antimycotic} and {species_dict[file]}; breakpoint_list set to [np.NaN, np.NaN]')
            breakpoint_list = [np.NaN, np.NaN]

        pred_well_number = int(MIC_dict[file][antimycotic][1:]) - distance # Use the well of ref MIC to calculate well of pred MIC (expressed relative to ref MIC)
        if np.isnan(breakpoint_list[0]):
            pred_cat = np.NaN
        elif pred_well_number >= breakpoint_list[0]: 
            pred_cat = 'S'
        elif pred_well_number <= breakpoint_list[1]: 
            pred_cat = 'R'
        else: 
            pred_cat = 'ATU'
        logging.debug(f'--- --- --- --- pred_well_number {pred_well_number}; breakpoint_dict {antimycotic}_{species_dict[file]}')
        logging.debug(f'--- --- --- --- --- pred_cat: {pred_cat}')

        MIC_well_number = int(MIC_dict[file][antimycotic][1:])
        if np.isnan(breakpoint_list[0]): 
            ref_cat = np.NaN
        elif MIC_well_number >= breakpoint_list[0]: 
            ref_cat = 'S'
        elif MIC_well_number <= breakpoint_list[1]: 
            ref_cat = 'R'
        else: 
            ref_cat = 'ATU'
        logging.debug(f'--- --- --- --- MIC_well_number {MIC_well_number}; breakpoint_dict {antimycotic}_{species_dict[file]}')
        logging.debug(f'--- --- --- --- --- pred_cat {pred_cat}, ref_cat: {ref_cat}')
        if (pred_cat == 'R') and (ref_cat == 'S'): 
            ME = 1
        elif (pred_cat == 'S') and (ref_cat == 'R'): 
            VME = 1
        elif (pred_cat == 'ATU') and (ref_cat == 'S'): 
            ATU_ME = 1
        elif (pred_cat == 'R') and (ref_cat == 'ATU'):
            ATU_ME = 1
        elif (pred_cat == 'S') and (ref_cat == 'ATU'): 
            ATU_VME = 1
        elif (pred_cat == 'ATU') and (ref_cat == 'R'): 
            ATU_VME = 1
        elif pred_cat == ref_cat:
            pass
        elif np.isnan(ref_cat): 
            VME, ME, ATU_VME, ATU_ME = np.NaN, np.NaN, np.NaN, np.NaN
        logging.debug(f'--- --- --- --- VME, ME, ATU_VME, ATU_ME: {VME, ME, ATU_VME, ATU_ME}')
        
        return VME, ME, ATU_VME, ATU_ME

def predict_MIC(antimycotic, file, ref_wells, parameter, timepoint, threshold, di, bias, eucast, pat):
    # detect desired pattern: e.g. [False, True, True, True]
    if not eucast:
        pattern = pat #[False, True, True, True]
    else:
        pattern = [False, True]
    logging.debug(f'--- --- --- --- Pattern: {pattern}')
    # for the present threshold, 'bool_list' will contain True/False in a specific order
    bool_list = []
    #logging.debug(f'--- --- --- --- predict() ref_wells is {ref_wells}')
    error = 'no_error'
    # ref_wells is: ['MIC12', 'MIC11', 'MIC10', 'MIC9', 'MIC8', 'MIC7' … 'MIC-12', 'GC']
    # Avoid looking for patterns in growth control, therefore loop through ref_wells[:-1]
    for ref_well in ref_wells[:-1]: # omit GC
        bool_list.append(di[antimycotic][file][ref_well][parameter][timepoint][threshold]['bool']) # bool_list yields e.g.: [{}, {}, {}, {}, {}, {}, {}, True, False, False, False, False, True, True, True, True]
    logging.debug(f'--- --- --- --- bool list {bool_list}')
    # >>> First detect predicted MICs that are limit values: that is, growth in highest log2 dilution or absence of growth in lowest log2 dilution
    # growth in highest log2 dilution with a certain threshold would translate to bool_list of [{}, {}, True, …]
    # absence of growth in lowest log2 dilution with a certain threshold: bool_list of [{}, {}, False, True, False, {}]       
    # previous_item is to check value of list item previous to the one considered in the loop: must be {} or beginning of list
    previous_item = None
    # First find the first occurence in bool_list of 'True': this implies growth at that cutoff in leftmost well
    for index, item in enumerate(bool_list):
        # item must be True and previous item cannot be True or False: only detect [{}, True] pattern and [True, …, …]
        if item and (type(previous_item) != bool):
            pattern_position = [(index + 0) - 1]  # - 1 because if reference MIC is limit-high value, it will be the well NEXT to the well with highest log2 dilution and growth
            # This is because a reference MIC that is a limit-high value implies that no measurements are available for this well
            # So we detect the pattern position as the extreme well where measurements are still available
            # And we extrapolate that the well next to it is the MIC: that's why we shift one well to the left using - 1
            # INDEX: 0 1 2 3 4 5 6 7 8 9 10 11 12 13
            # An index of 7 is placed where on a scale that goes from 1 to 12 to MIC to -1 to -12 (25 elements in total, starting at 1)?
            # It is placed at well 8, which is +5 log2 dilutions from MIC
            # WELLS: 1 2 3 4 5 6 7 8 9 1011 12 MIC
            logging.debug(f'--- --- --- --- --- Limit high value detected: index {index}, pattern_pos {pattern_position}')
            break # as soon as True detected, break out of loop and retain index of True: highest log2 dilution with growth
        elif item == False and (type(previous_item) != bool):
            # Patterns that include multiple False (like False, False, False, False, True), will never be found for high MICs (bool_list None None False True True)
            # Therefore, add False in the beginning according to the number of False in pattern  
            pass
        previous_item = item
    # Paradoxically, often when growth at highest log2 dilution, at lower dilutions, growth is LESS (Eagle type resistance?)
    # This means that, if we have detected growth at highest log2 dilution, we should neglect eventual [False, {}] or […, False] patterns in lower log2 dilution wells
    # That's why we put what follows in else clause: only if NO limit high value detected, look for a limit low value
    else: # if no [{}, True] or [True, …] pattern detected: go on to detect [False, {}] and […, False] pattern
        bool_list_extended = bool_list.copy()
        # For this, reverse order of list and detect first occurence of False, following {} or as first item in the list
        previous_item = None
        for index, item in enumerate(bool_list[::-1]):
            if item == False and (type(previous_item) != bool):
                pattern_position = [len(bool_list) - (index + 1) + 1] # (index + 1) because 0 is included in index
                logging.debug(f'--- --- --- --- --- Limit low value detected: index {index}, pattern_pos {pattern_position}')   
                break # as soon as False detected, break out of loop and retain index of lowest False: lowest log2 dilution where no growth

            # Now we're iterating through reverse bool_list, detect trailing True's (see explanation below)
            if item == True and (type(previous_item) != bool):
                # First determine how many True's should be added, according to the number of trailing True's in pattern
                number_Trues = pattern.count(True) - 1 # minus one because only if a trailing True is detected, will we attach the remaining number of required True's
                bool_list_extended[(len(bool_list_extended) - index):len(bool_list_extended) - index + number_Trues] = [True] * number_Trues
                logging.debug(f'--- --- --- --- --- Trailing True detected: add {number_Trues} Trues to bool_list')   
            previous_item = item

        # >>> DETECT PATTERN
        else: # if no limit high/low value detected, go on to detect [False, True, True, True] pattern or variant of it
            # Determine the position of this pattern in the range of wells
            # This position corresponds to the position of the predicted MIC
            # One cave at: what if growth only in well +11 or +12: […, {}, {}, False, True, {}, {}]? Pattern will not be detected
            # That's why, if there's a part of the pattern occurring at end of bool_list: append True's

            pattern_position = [x for x in range(len(bool_list)) if bool_list_extended[x:x+len(pattern)] == pattern]
            # 'Pattern position' contains the position of a pattern across ALL WELLS belonging to a certain antimycotic/file/parameter/timepoint
            if len(pattern_position) == 1 and not eucast: # Exclude EUCAST approach: here we don't differentiate multiple from single detection
                # one position where the pattern is detected
                logging.debug(f'--- --- --- --- --- Pattern detected at position {pattern_position}')
                # apply bias
                pattern_position = [pattern_position[0] + bias]
                logging.debug(f'--- --- --- --- --- Pattern position after bias correction {pattern_position}')
            else:
                if len(pattern_position) > 1 and not eucast: # when required pattern is detected multiple times: set pattern outside range of valid wells, at -1
                    pattern_position = [-1] # putting pattern position at -1 maximizes the distance to the MIC, that is at position 12
                    error = 'multiple'

                elif not pattern_position: # if pattern not found, similarly, set position at -1
                    if all(item is None for item in bool_list_extended):
                        logging.debug(f'--- --- --- --- --- No data for pattern-detection: all values None')
                        error = 'no_data'
                        pattern_position = [np.NaN]
                    else:
                        pattern_position = [-1]
                        error = 'not_found'
                else: # This is case of EUCAST approach: both multiple or single pattern detections included: select rightmost pattern position + apply bias correction 
                    logging.debug(f'--- --- --- --- --- Pattern detected at position {pattern_position}')
                    pattern_position = [pattern_position[-1] + bias]
                    logging.debug(f'--- --- --- --- --- Pattern position after bias correction {pattern_position}')


                logging.debug(f'--- --- --- --- --- Pattern error: {error}')

    return pattern_position, error

def get_index(timepoint, file, d):
    margin = 2000
    time_column = d[file]
    index = np.where((time_column > (timepoint - margin)) & (time_column < (timepoint + margin))) # is empty array in case run was shorter than timepoint!
    index = np.asarray(index)
    if index.size == 0: return np.NaN
    else: return index[0][0] # to convert multidimensional arrays to single values: [[23 24][ 0  0]] becomes 23

def define_thresholds(timepoint, antimycotic, parameters, steps, filelist, d, thresholds, eucast, threshold_dict, training):

    logging.debug(f'--- Predict(timepoint: {timepoint}, {len(antimycotic)} antimycotic(s), files: {filelist}, {len(parameters)} parameter(s))')
    thresholds[timepoint] = {}
    thresholds[timepoint][antimycotic] = {}
    if not eucast:
    # 'thresholds' dictionary contains a range of n thresholds between min/max values across all wells/files per antimycotic and timepoint
    # for training, a random subsample of filelist will have been selected, and supplied to predict as "filelist" (after excluding files without data?)
        if training:
            logging.info(f'--- --- --- Proceeding to setting thresholds based on {steps} steps.')
        else: # if not, fetch validation thresholds from a dictionary
            parameter = '_'.join(threshold_dict[f'{antimycotic}_{timepoint}'].split('_')[:-2]) # format threshold_dict, e.g. {'A_46800': 'SESAfungi_ratio_-2_0.9886754950370688'}
            threshold = [float(threshold_dict[f'{antimycotic}_{timepoint}'].split('_')[-1])]
            bias = int(threshold_dict[f'{antimycotic}_{timepoint}'].split('_')[-2])
            logging.debug(f'--- --- parameter {parameter}, threshold {threshold}, bias {bias}')
            thresholds[timepoint][antimycotic][parameter] = threshold
            return thresholds, [parameter]
    else: # For EUCAST approach, also initialize full thresholds dict hierarchy (a bit redundant since same float applied to all parameters - this is for uniformity in downstream analyses)
        for parameter in parameters: 
            thresholds[timepoint][antimycotic][parameter] = [threshold_dict[antimycotic]] 
        return thresholds, parameters
    logging.debug(f'--- --- Thresholds dictionary: {thresholds}')

    # First, get range of min-max values for this antimycotic and timepoint across all files and across all wells
    # Put all values for this antimycotic, for this timepoint, across all files in thresholds_min_max dictionary: 
    # Initialize thresholds_min_max as a dictionary of lists per parameter
    thresholds_min_max = {}
    for parameter in parameters:
        thresholds_min_max[parameter] = []

    # Populate thresholds_min_max: min/max values across all files and all wells
    for file in filelist:
        index = get_index(timepoint, file, d)
        # Consider for this antimycotic and for all files and all wells and all parameters the values at timepoint for present antimycotic
        for well, values in d[antimycotic][file].items():
            for parameter in parameters:
                #print(f'{antimycotic}, {file}, {well}, {parameter}: \n {d[antimycotic][file][well][parameter][0:2]}')
                try: 
                    value = d[antimycotic][file][well][parameter][index] # to capture KeyErrors
                    thresholds_min_max[parameter].append(value)
                except Exception as e: pass #print(f'for par {parameter} this threshold exception: {e}')

    # Define min-max values in this list that contains all values of all wells of all files for each parameter for antimycotic/timepoint
    # Then define define range of cutoffs accordingly: split range between min/max values in n parts ('steps')
    # Put a list of n cutoffs in thresholds dictionary: THIS IS SPECIFIC PER ANTIMYCOTIC/TIMEPOINT AND IS RE-INITIALIZED FOR EACH ANTIMYCOTIC/TIMEPOINT
    for parameter in parameters:
        min_ = min(thresholds_min_max[parameter])
        max_ = max(thresholds_min_max[parameter])
        # Alternatively, define min/max as a smaller range to exclude outliers, after sorting from lowest to highest
        lowest_to_highest = []
        lowest_to_highest = thresholds_min_max[parameter].sort()
        outlier_exclusion = int(0.1 * len(filelist))
        min_ = thresholds_min_max[parameter][outlier_exclusion]
        max_ = thresholds_min_max[parameter][-(outlier_exclusion + 1)]
        try: 
            thresholds[timepoint][antimycotic][parameter] = [*np.linspace(min_, max_, steps)]
        except Exception as e:
            logging.warning(f'Error creating thresholds dict entry: {e}')
            thresholds[timepoint][antimycotic][parameter] = []  
    logging.info(f'--- --- {steps} threshold(s) defined')

    return thresholds, parameters

def categorise(antimycotic, filelist, thresholds, parameters, timepoint, di, d, eucast):
    # Iterate through all files
    eucast_threshold_values = {}
    for file in filelist:
        logging.info(f'--- --- --- Categorise {file}')
        if eucast: eucast_threshold_values[file] = {} # leave empty dict if non-EUCAST, so plot() knows threshold is thresh_y and NOT specific to file
        di[antimycotic][file] = {}
        # Get the row index of datapoints for this timepoint
        index = get_index(timepoint, file, d)

        # >>> First populate ['bool'] per well, timepoint, threshold with True or False, according to whether threshold was passed or not
        # Iterate through wells in d with data arrays per parameter
        for well, values in d[antimycotic][file].items(): # e.g. well: MIC-12, value: {'BCA': array([], dtype=float64), 'BCA_rate': array([], dtype=float64)
            di[antimycotic][file][well] = {} # Copy well names to dict di
            logging.debug(f'--- --- --- --- Categorising {len(parameters)} parameter(s) for well {well}')
            for parameter in parameters: # Copy parameter names to dict di
                di[antimycotic][file][well][parameter] = {}
                di[antimycotic][file][well][parameter][timepoint] = {} # where in d timepoints and data are in a single array, di differentiates different timepoints
                for threshold in thresholds[timepoint][antimycotic][parameter]: # Iterate over all pre-defined thresholds and define them in dict di
                    di[antimycotic][file][well][parameter][timepoint][threshold] = {}
                    if (values[parameter].size == 0) or (np.isnan(index)):
                        di[antimycotic][file][well][parameter][timepoint][threshold]['bool'] = None
                    # If data is available: assign True if datapoint surpasses threshold and else False
                    elif not eucast: # If non EUCAST approach, bool value assigned according to absolute threshold
                        if values[parameter][index] > threshold: 
                            pos = True
                        else: 
                            pos = False
                        di[antimycotic][file][well][parameter][timepoint][threshold]['bool'] = pos
                    else: # EUCAST approach: bool value assigned according to proportion compared to GC
                        # this is only valid for parameters with min 0
                        eucast_value = threshold * d[antimycotic][file]['GC'][parameter][index]
                        if values[parameter][index] > eucast_value: 
                            pos = True
                        else: 
                            pos = False
                        di[antimycotic][file][well][parameter][timepoint][threshold]['bool'] = pos
                        eucast_threshold_values[file][parameter] = eucast_value
                        logging.debug(f"--- --- --- --- --- {pos} for EUCAST threshold {eucast_value} ({threshold} * {d[antimycotic][file]['GC'][parameter][index]})")
    # This dictionary is only used for plotting of thresholds
    logging.debug(f"--- --- --- --- --- EUCAST threshold dict: {eucast_threshold_values}")
    return eucast_threshold_values

def predict(timepoint, antimycotic, filelist, parameters, ref_wells, steps, thresholds, di, MIC_dict, breakpoint_dict, threshold_dict, d, species_dict, output_dir, bias, eucast, pattern):
    di[antimycotic] = {}
    highlight = {}

    # categorise() assigns True/False in 'bool' entry per well, timepoint, parameter
    # e.g. di['A']['A-93624221 72 scan areas.xlsx']['MIC3']['BCA'][43200].items() yields: dict_items([(7.010701392915223, {'bool': True}), (8.37132208988804, {'bool': False})])
    # If eucast True, categorise() determines thresholds per file based on proportions intialized in threshold dict: we need time index for this and it is already used in categorise() 
    # Therefore, categorise() also returns a dict that per file contains the EUCAST criteria set threshold - eucast_threshold_dict is empty for non-EUCAST mode  
    eucast_threshold_dict = categorise(antimycotic, filelist, thresholds, parameters, timepoint, di, d, eucast)
    no_data = []

    # Now zoom out, and per timepoint and threshold detect patterns across wells to predict MIC
    for file in filelist:
        highlight[file] = {}
        logging.debug(f'--- Predicting MIC for file {file} and {len(parameters)} parameter(s)')
        for parameter in parameters:
            highlight[file][parameter] = {}
            # Iterate over all thresholds; arbitrarily use 'MIC' entry to get the thresholds ([timepoint.keys()])
            logging.debug(f'--- --- CA for parameter: {parameter}')
            random_well_key = list(di[antimycotic][file].keys())[0] # Pick any well because thresholds are defined for each
            for threshold in di[antimycotic][file][random_well_key][parameter][timepoint].keys():
                logging.debug(f'--- --- --- Threshold {threshold}')
                pattern_position, error = predict_MIC(antimycotic, file, ref_wells, parameter, timepoint, threshold, di, bias, eucast, pattern)
                # Calculate distance between predicted MIC and reference MIC
                # position of MIC is always in the middle for training files, that is position 12/24 - calculate distance by: MIC_position - pattern_position
                #distance = 12 - int(*pattern_position)

                if not np.isnan(pattern_position[0]):
                    distance = int((len(ref_wells) - 2) / 2) - int(*pattern_position) # -2 because 'GC' and ... included
                    # Determine categorical agreement
                    VME, ME, ATU_VME, ATU_ME = determine_CA(file, distance, antimycotic, MIC_dict, breakpoint_dict, species_dict)
                    highlight[file][parameter][threshold] = distance
                else: # pattern_position is np.NaN when all values are None, i.e. when for this timepoint file has no data. You can only conveniently find this out by iterating through ref_wells for each file, as predict_MIC() does
                    distance = np.NaN
                    VME, ME, ATU_VME, ATU_ME = np.NaN, np.NaN, np.NaN, np.NaN
                    # when this is the case, record files without data for this specific timepoint in list
                    no_data.append(file)
                # now put the distance between the pattern and the MIC and associated categorical errors for this particular threshold under the 'MIC' entry
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['MIC_distance'] = distance
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['VME'] = VME
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ME'] = ME
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_VME'] = ATU_VME
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_ME'] = ATU_ME
                di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['error'] = error
                logging.debug(f'--- --- --- - Dist {distance}: categorical errors: {VME} VME, {ME} ME, {ATU_VME} ATU_VME, {ATU_ME} ATU_ME, prediction errors: {error}')

    """
    for parameter in parameters:
        for threshold in di[antimycotic][file][random_well_key][parameter][timepoint].keys():
            logging.debug(f'Plotting for {parameter} and {threshold}')
            plot.plot(output_path=output_dir, d=d, di=di, params=[parameter], am=[antimycotic], files=filelist,
                     timeframe=((0, 60000)), # timepoint - 0.2 * timepoint, timepoint + 0.2 * timepoint
                     well_legend=ref_wells, thresh_y=threshold, timepoint=timepoint, highlight=highlight,
                     thresh_dict = eucast_threshold_dict,
                    label=f"distance: {distance}")#; thresh_y {performance_per_timepoint['threshold']}; distances {performance_per_timepoint['dist_per_file']}; CE per file: {performance_per_timepoint['CE_per_file']}")#, y_lim=(8,9))
    """
    
    return set(no_data)

