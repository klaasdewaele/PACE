"""
evaluate.py — threshold optimisation and performance evaluation.

Provides functions to:

  - ``get_agreement_parameters_short()`` — fast proxy metrics (total MIC
    distance, weighted categorical error count) for inner optimisation loops.
  - ``get_agreement_parameters_ext()`` — full performance metrics including
    essential agreement (EA), categorical agreement (CA), VME/ME/ATU rates,
    and per-file distance breakdowns.
  - ``get_best_thresholds_bis()`` — selects the optimal threshold per
    kinetic parameter according to one of several criteria (minimise total
    MIC distance, minimise weighted errors, maximise EA or CA, or enforce
    a minimum EA/CA target).
  - ``get_best_parameter()`` — across all parameters with their optimal
    thresholds, selects the single best parameter-threshold pair and writes
    the result to disk for subsequent bias correction and plotting.
"""

import logging
import os
import pickle
import statistics

import numpy as np
import pandas as pd

from . import plot

def get_agreement_parameters_short(antimycotic, timepoint, parameter, threshold, di):
    logging.debug(f'--- --- --- --- Getting MIC-distance and weighted CA for {threshold}')
    agreement_dict = {}
    # Now all MIC_distances and errors have been summed across all files per threshold, check which threshold is the best performing, both for EA and CA
    # Because VME is worse than ME, calculate weighted sum of number of errors
    # weighted_CA will be np.NaN in case no breakpoints were available
    agreement_dict['total_MIC_distance'] = di[antimycotic][timepoint][parameter][threshold][0]
    VME_tot = di[antimycotic][timepoint][parameter][threshold][1]
    ME_tot = di[antimycotic][timepoint][parameter][threshold][2]
    ATU_VME_tot = di[antimycotic][timepoint][parameter][threshold][3] 
    ATU_ME_tot = di[antimycotic][timepoint][parameter][threshold][4]
    weighting_factor = 2 # VMEs are counted as double errors - also use factor for ATU_VMEs?
    # Note: confusingly, here agreement_dict['weighted_CA'] is not a categorical agreement, but rather the sum of categorical errors
    # In get_agreement_parameters_ext, agreement_dict['CA'] IS the categorical agreement
    agreement_dict['weighted_CA'] = (VME_tot * weighting_factor) + ME_tot + (ATU_VME_tot) + ATU_ME_tot
    agreement_dict['unweighted_CA'] = VME_tot + ME_tot + ATU_VME_tot + ATU_ME_tot
    # If prediction failed for a sample for a threshold, this is counted as an ME; see determine_CA() with distance 12
    logging.debug(f'--- --- --- --- MIC distance is {agreement_dict["total_MIC_distance"]}; weighted_CA is {agreement_dict["weighted_CA"]}')
    
    # In case no breakpoints are available, all values become 0
    return agreement_dict

def get_agreement_parameters_ext(antimycotic, timepoint, parameter, threshold, di, files_with_data, files_with_breakpoint):
    logging.debug(f'--- --- Getting extended agreement parameters for {threshold}, {antimycotic}, {timepoint}, {parameter}')

    agreement_dict = {}
    list_multiple = []
    list_not_found = []
    CE_per_file = {} # categorical errors per file
    dist_per_file = {}
    # apart from MIC distance and categorical errors, also collect the occurence of prediction errors for that threshold per file
    # These prediction errors are collected in a list that contains the file-names where each of the errors occurred
    for file in files_with_data:
        error = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['error']
        if error == 'multiple':
            list_multiple.append(file)
        elif error == 'not_found':
            list_not_found.append(file)                 
        # Similarly, dist_per_file is dictionary that contains distances for this threshold per file, CE_per_file the categorical errors
        dist_per_file[file] = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['MIC_distance']
        # CE_per_file contains categorical errors per file: loop through keys of di[antimycotic][file]['MIC'][parameter][timepoint][threshold]
        # key, value = {'VME': 0, 'ME': '0', 'ATU_VME': 1, etc} — CE_per_file = {'filename.xlsx': 'VME', 'filename2.xlsx': 'ME', etc}
        CE_present_file = dict([(file, key) for key, value in di[antimycotic][file]['MIC'][parameter][timepoint][threshold].items() if ('ME' in key) and (value == 1)])
        # {di[antimycotic][file]['MIC'][parameter][timepoint][threshold].items()} yields:
        # dict_items([('bool', True), ('MIC_distance', 12), ('VME', 0), ('ME', 1), ('ATU_VME', 0), ('ATU_ME', 0), ('error', 'no_error')])
        CE_per_file = {**CE_per_file, **CE_present_file}
        
    agreement_dict['multiple'] = list_multiple
    agreement_dict['not_found'] = list_not_found
    agreement_dict['distance_per_file'] = dist_per_file
    agreement_dict['CE_per_file'] = CE_per_file
      
    # Create di[antimycotic][timepoint]["best-threshold-per-parameter"] entry that describes summed MIC distances and multiple-/not_found-errors
    agreement_dict['threshold'] = threshold
    # In case no breakpoint available, di[antimycotic][timepoint][parameter][threshold][0-4] will be np.NaN
    agreement_dict['total_MIC_distance'] = di[antimycotic][timepoint][parameter][threshold][0]
    VME_tot = di[antimycotic][timepoint][parameter][threshold][1]
    ME_tot = di[antimycotic][timepoint][parameter][threshold][2]
    ATU_VME_tot = di[antimycotic][timepoint][parameter][threshold][3] 
    ATU_ME_tot = di[antimycotic][timepoint][parameter][threshold][4]
    agreement_dict['VME'] = VME_tot
    agreement_dict['ME'] = ME_tot
    agreement_dict['ATU_VME'] = ATU_VME_tot
    agreement_dict['ATU_ME'] = ATU_ME_tot
    weighting_factor = 2 # VMEs are counted as double errors - also use factor for ATU_VMEs?
    agreement_dict['min_errors_weighted'] = (VME_tot * weighting_factor) + ME_tot + (ATU_VME_tot * weighting_factor) + ATU_ME_tot
    agreement_dict['min_errors_unweighted'] = VME_tot + ME_tot + ATU_VME_tot + ATU_ME_tot

    denominator_EA = len(files_with_data)
    denominator_CA = len(files_with_breakpoint)
    logging.debug(f'--- --- --- denominator EA {denominator_EA}')
    logging.debug(f'--- --- --- denominator CA {denominator_CA}')

    #denominator = len(filelist) # use as denominator the total number of files used for valiation
    # An alternative denominator would be the number of succesful analyses, excluding failed analysis due to multiple/not_found errors
    # To calculate this denominator use denominator - (multiple_EA + not_found_EA)

    if denominator_CA != 0:
        agreement_dict['VME_rate'] = VME_tot / denominator_CA
        agreement_dict['ME_rate'] = ME_tot / denominator_CA
        agreement_dict['ATU_VME_rate'] = ATU_VME_tot / denominator_CA
        agreement_dict['ATU_ME_rate'] = ATU_ME_tot / denominator_CA
        agreement_dict['CA'] = (denominator_CA - VME_tot - ME_tot - ATU_VME_tot - ATU_ME_tot) / denominator_CA
    else:
        agreement_dict['VME_rate'] = np.NaN
        agreement_dict['ME_rate'] =  np.NaN
        agreement_dict['ATU_VME_rate'] = np.NaN
        agreement_dict['ATU_ME_rate'] = np.NaN
        agreement_dict['CA'] = np.NaN

    agreement_dict['failed_rate'] = (len(list_multiple) + len(list_not_found)) / denominator_EA
    CE = 0
    for value in dist_per_file.values():
        if (value > 1) or (value < -1):
            CE += 1
    agreement_dict['EA'] = (denominator_EA - CE) / denominator_EA
    agreement_dict['d_EA'] = denominator_EA
    agreement_dict['d_CA'] = denominator_CA
    agreement_dict['f_EA'] = files_with_data
    agreement_dict['f_CA'] = files_with_breakpoint
    #logging.debug(f'agreement_dict {agreement_dict["EA"]}')

    return agreement_dict


def get_best_thresholds_bis(timepoint, antimycotic, parameters, di, thresholds, files_with_data, files_with_breakpoint, criteria, output_dir, skip_flag):
    logging.info(f'--- --- --- Determining best threshold for {len(parameters)} parameter(s) for antimycotic {antimycotic} and timepoint {timepoint}')
    """
    We have calculated the distance between pattern and MIC per antimycotic and per file individually for each threshold
    
    Now we can sum these distances to get sum of MIC distances - we can calculate categorical errors, and also 'higher level' metrics, such as essential and categorical agreeement
    In this function, We'll select the best performing threshold PER PARAMETER, according to the criterion that is supplied (min_distance, min_errors_weighted, EA, CA, or a combination of errors)
    In the function get_best_parameter(), we'll loop over all parameters (with their best threshold) and select the best parameter    

    Note, there are three ways of selecting 'optimal' threshold for either EA or CA:
    1) Approach 1: Optimise EA/CA with as proxy parameters the minimal number of weighted categorical errors (if optimising CA) or minimal total sum of MIC distances (if optimising EA)
       When there are multiple matching thresholds, secondarily select threshold with lowest total MIC distance and lowest number of categorical errors, respectively
       Note: this approach is faster since there's no need to calculate for every of many thresholds the essential and categorical agreement and all associated categorical errors. It is implemented in get_agreement_parameters_short()). 
       However, when optimising EA or CA, it is better to maximise EA and CA directly instead of using total MIC distance and total number of categorical errors as proxies)
        Approach 2 does this and when computational factors are not limited, should be the default approach.
    2) Approach 2: Optimise EA/CA with criterion the best essential agreement (if optimising EA), or best categorical agreement (optimising CA)
       If multiple matches, secondarily select threshold with best CA and EA respectively
    3) Approach 3: Optimise EA/CA with a fixed threshold for essential agreement or categorical agreement (e.g. at least 0.90 EA) - this is necessary for antimycotics where CA/EA tradeoff is prominent (itraconazole for example)
    """

    # Initialize new hierarchy for storing sum of errors per threshold
    di[antimycotic][timepoint] = {}
    EA_CA, criterium = criteria
    # dict for each threshold with associated errors: will be used for multi-objective optimisation algorithm
    performance_dict = {}
    counter = 0
    for parameter in parameters:
        performance_dict[parameter] = {}

        logging.debug(f'\n--- --- Par {parameter}')
        # Initialize per antimycotic, timepoint, parameter and threshold, a list that contains [total MIC distance, #VME, #ME, #ATU_VME, #ATU_ME]
        di[antimycotic][timepoint][parameter] = {}
        for threshold in thresholds[timepoint][antimycotic][parameter]:
            di[antimycotic][timepoint][parameter][threshold] = [0, 0, 0, 0, 0]
        # Initialize per antimycotic, timepoint and parameter a dict that will contain performance of best threshold for that parameter
        di[antimycotic][timepoint][f'{parameter}_best_EA'] = {}
        di[antimycotic][timepoint][f'{parameter}_best_CA'] = {}

        # Approach 1: variables to track the lowest MIC distance and total weighted errors
        minimum_EA = np.inf # minimum_EA and _CA store the lowest essential/categorical error rate (i.e. MIC distance and weighted errors total, respectively) across all thresholds = approach 1
        minimum_CA = np.inf
        # Approach 2: variables to track the highest EA and CA
        maximum_EA = 0 # idem but store essential agreement and categorical agreement (not total MIC distance and sum weighted errors) = approach 2
        maximum_CA = 0
        list_best_EA = [] # lists that contain the thresholds meeting criterion
        list_best_CA = []

        for threshold in thresholds[timepoint][antimycotic][parameter]:
            performance_dict[parameter][threshold] = {}

            logging.debug(f'--- --- --- Counting errors across {len(files_with_data)} files for threshold {threshold}')
           
            # This loop stores MIC_distance and categorical errors in a list of five items (MIC_distance, VME, ME, ATU_VME, ATU_ME)
            ### IS IT SPECIFIC TO PURPOSE 1?
            for file in files_with_data:
                di[antimycotic][timepoint][parameter][threshold][0] += abs(di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['MIC_distance'])
                if file in files_with_breakpoint:
                    # if no breakpoint, these values remain 0
                    # Fix: check for NaN before adding, as Z-prefixed MICs result in NaN errors
                    val_VME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['VME']
                    val_ME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ME']
                    val_ATU_VME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_VME']
                    val_ATU_ME = di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_ME']

                    di[antimycotic][timepoint][parameter][threshold][1] += val_VME if not np.isnan(val_VME) else 0
                    di[antimycotic][timepoint][parameter][threshold][2] += val_ME if not np.isnan(val_ME) else 0
                    di[antimycotic][timepoint][parameter][threshold][3] += val_ATU_VME if not np.isnan(val_ATU_VME) else 0
                    di[antimycotic][timepoint][parameter][threshold][4] += val_ATU_ME if not np.isnan(val_ATU_ME) else 0                   
            
            # if approach 2 or 3: slower approaches that involve calculation of all performance metrics for each threshold
            if type(criterium) == float or criterium == "EA" or criterium == "CA": 
                logging.debug(f'Approach 2 or 3: getting extended agreement_dict for {parameter} and {threshold}')
                # in these cases, for each threshold, perform extended calculation of essential agreement, categorical agreement, VME, etc
                agreement_dict = get_agreement_parameters_ext(antimycotic, timepoint, parameter, threshold, di, files_with_data, files_with_breakpoint)
                present_EA = agreement_dict['EA'] # here, "present_EA" and "CA" are truly essential/categorical agreement - however, in below else clause, "present_EA" actually means total "MIC_distance", and "present_CA" means sum_weighted_errors
                present_CA = agreement_dict['CA']
                present_weighted_cat_errors = agreement_dict['min_errors_weighted'] # this is specifically used in approach 2, where for thresholds with best EA, we secondarily select thresholds with best weighted errors, rather than best CA
                # We'll put this information also in a separate dict, for multi-objective optimisation algorithm (independent from this algorithm)
                logging.debug(f'writing agreement dict for {threshold} for {parameter} to performance_dict')
                for metric, performance in agreement_dict.items():
                    logging.debug(f'>>> for EA/CA {EA_CA} and criterium {criterium} metric {metric} has performance: {performance}')
                    performance_dict[parameter][threshold][metric] = performance
                
            # if approach 1: fastest approach
            else:
                # For this threshold, only determine total MIC distance and weighted categorical agreement
                agreement_dict = get_agreement_parameters_short(antimycotic, timepoint, parameter, threshold, di)
                present_EA = agreement_dict['total_MIC_distance']
                present_CA = agreement_dict['weighted_CA'] # reminder: weighted_CA is the weighted sum of categorical errors, not categorical agreement!
                # in case no breakpoints available, such as for Olorofim, present_CA will be 0

            # if the threshold under considerations meets specified criteria, add it to a list of all thresholds that meet those criteria, named 'list_best_EA' or 'list_best_CA', according to what we're optimising
            # approach 3:
            if type(criterium) == float:
                # Regardless of whether we're optimising EA or CA, always check for this threshold if it is above target for either EA and CA
                # Note, here we do not only record threshold with its EA/CA performance, but also the number of categorical errors, which can be useful for further subselecting thresholds
                if present_EA >= criterium:
                    list_best_EA.append((threshold, present_EA, present_CA, agreement_dict['VME'], agreement_dict['ME'], agreement_dict['ATU_VME'], agreement_dict['ATU_ME'], agreement_dict['total_MIC_distance'], agreement_dict['min_errors_weighted'], agreement_dict))
                if present_CA >= criterium:
                    list_best_CA.append((threshold, present_EA, present_CA, agreement_dict['VME'], agreement_dict['ME'], agreement_dict['ATU_VME'], agreement_dict['ATU_ME'], agreement_dict['total_MIC_distance'], agreement_dict['min_errors_weighted'], agreement_dict))

            # in approach 1 and 2, select best threshold for EA/CA and secondarily for CA/EA - these approaches are vulnerable to CA/EA trade-off
            # approach 1:
            elif criterium != 'CA' and criterium != 'EA': # approach 1: present_EA will be total MIC_distance or weighted_sum_errors - always subselect thresholds optimised for EA secondarily for CA, and reverse
                if present_EA <= minimum_EA:
                    # If present_EA at least equal to minimum_EA: add to list_best
                    # No need to reset minimum_EA
                    list_best_EA.append(f'{threshold}_{present_CA}') # Perhaps better to use tuple instead of string
                    logging.debug(f'--- --- --- --- --- Added this threshold to list_best_EA: {list_best_EA}')
                    if present_EA < minimum_EA:
                        # If in addition it is lower than minimum_EA, reset minimum_EA, reset list_best and append first entry
                        list_best_EA = []
                        list_best_EA.append(f'{threshold}_{present_CA}')
                        logging.debug(f'--- --- --- --- --- --- New best threshold. New minimum: {present_EA}. Reset list_best_EA: {list_best_EA}')
                        minimum_EA = present_EA
                    
                if present_CA <= minimum_CA:
                    list_best_CA.append(f'{threshold}_{present_EA}')
                    logging.debug(f'--- --- --- --- --- Added this threshold to list_best_CA: {list_best_CA}')
                    if present_CA < minimum_CA:
                        # If in addition it is lower than minimum_EA, reset minimum_EA, reset list_best and append first entry
                        list_best_CA = []
                        list_best_CA.append(f'{threshold}_{present_EA}')
                        logging.debug(f'--- --- --- --- --- --- New best threshold. New minimum: {present_CA}. Reset list_best_CA: {list_best_CA}')
                        minimum_CA = present_CA
            else: # approach 2: criterium is EA or CA: maximise EA and CA
                if present_EA >= maximum_EA:
                    list_best_EA.append(f'{threshold}_{present_weighted_cat_errors}') # Perhaps better to use tuple instead of string
                    logging.debug(f'--- --- --- --- --- Added this threshold to list_best_EA: {list_best_EA}')
                    if present_EA > maximum_EA:
                        list_best_EA = []
                        list_best_EA.append(f'{threshold}_{present_weighted_cat_errors}')
                        logging.debug(f'--- --- --- --- --- --- New best threshold. New maximum: {present_EA}. Reset list_best_EA: {list_best_EA}')
                        maximum_EA = present_EA 
                if present_CA >= maximum_CA:
                    list_best_CA.append(f'{threshold}_{present_EA}')
                    logging.debug(f'--- --- --- --- --- Added this threshold to list_best_CA: {list_best_CA}')
                    if present_CA > maximum_CA:
                        list_best_CA = []
                        list_best_CA.append(f'{threshold}_{present_EA}')
                        logging.debug(f'--- --- --- --- --- --- New best threshold. New maximum: {present_CA}. Reset list_best_CA: {list_best_CA}')
                        maximum_CA = present_CA
        
        # approach 3 
        # If thresholds not selected based on minimizing EA/CA or min_distance/weighted errors, but have been added to list_best_EA/CA based on them exceeding a target EA/CA
        # Then we still have to determine which threshold is best
        final_target = criterium
        if type(criterium) == float:
             # initialise agreement dictionaries for EA and CA
             agreement_dict_EA = {}
             agreement_dict_CA = {}
             # at this point, you have list of thresholds meeting target for EA and for CA - first check if any thresholds meet this criterion, if not try with lower target
             if (len(list_best_EA) == 0 and EA_CA == "EA") or (len(list_best_CA) == 0 and EA_CA == "CA"):
                 counter += 1 # count every time we skip a parameter; if it equals the length of parameters, then no parameters meet target
                 # in that case, rerun get_best_thresholds_bis with a lower target
                 if counter == len(parameters):
                     new_target = round(criterium - 0.01, 2)
                     if new_target < 0: # no targets have been found that meet criterion
                         return np.NaN, True 
                     criteria = (EA_CA, new_target)
                     logging.warning(f'--- --- --- No thresholds found that meet target of {criterium} for {EA_CA} for any parameter - retry with lower target: {new_target}')
                     final_target, skip_flag = get_best_thresholds_bis(timepoint, antimycotic, parameters, di, thresholds, files_with_data, files_with_breakpoint, criteria, output_dir, skip_flag)
                 else: # until we have reached last parameter, just skip parameter if it doesn't have thresholds that meet target
                     #logging.info(f'list best ea {list_best_EA}, {EA_CA}, list best ca {list_best_CA}')
                     logging.debug(f'For {parameter} no thresholds found above target of {criterium}')
                     continue
             else:
                 # an entry in list_best_EA looks like this: ((threshold, present_true_EA, present_true_CA, agreement_dict['VME'], agreement_dict['ME'], agreement_dict['ATU_VME'], agreement_dict['ATU_ME']))
                 if len(list_best_EA) != 0:
                     # for thresholds with EA above target, sort according to desired specifications
                     df = pd.DataFrame(list_best_EA, columns=['threshold', 'EA', 'CA', 'VME', 'ME', 'ATU_VME', 'ATU_ME', 'total_MIC_distance', 'min_errors_weighted', 'dict'])
                     # define the other metrics (apart from EA) for which you want to optimise
                     # e.g. the below list is optimised for voriconazole: among thresholds above 90%, optimise weighted CA and VME/VME_ATU first (that is, get rid of VMEs) 
                     # then MEs and ATU_MEs, then CA, and finally total MIC_distance and EA
                     priority_list = (['min_errors_weighted', 'VME', 'ATU_VME', 'ME', 'ATU_ME', 'CA', 'total_MIC_distance', 'EA'], [True, True, True, True, True, False, True, False])
                     df.sort_values(priority_list[0], ascending=priority_list[1], inplace=True)
                     agreement_dict_EA = df.iloc[0]['dict']
                     logging.info(f'\nThese are the thresholds to choose from for EA:\n{df}')
                     logging.info(f'This one we eventually chose:\n{df.iloc[0]}\n')
                 if len(list_best_CA) != 0:
                     df = pd.DataFrame(list_best_CA, columns=['threshold', 'EA', 'CA', 'VME', 'ME', 'ATU_VME', 'ATU_ME', 'total_MIC_distance', 'min_errors_weighted', 'dict'])
                     # in case we have an antimycotic with good EA but lesser CA (e.g. when reference MICs are close to breakpoint, so that 1 log2 dilution errors result in categorical errors)
                     # here, we have preselected thresholds with CA above for example 90%: for these we may now subselect to optimise 1) for grave categorical errors VMEs, 2) for EA
                     priority_list = (['min_errors_weighted', 'VME', 'ATU_VME', 'ME', 'ATU_ME', 'total_MIC_distance', 'EA', 'CA'], [True, True, True, True, True, True, False, False])
                     df.sort_values(priority_list[0], ascending=priority_list[1], inplace=True)
                     agreement_dict_CA = df.iloc[0]['dict'] 
                     logging.info(f'\nThese are the thresholds to choose from for CA:\n{df}')
                     logging.info(f'This one we eventually chose:\n{df.iloc[0]}\n')

        elif criterium != 'CA' and criterium != 'EA':
            logging.info(f'--- --- --- For {parameter}, list_best_EA contains {len(list_best_EA)} entries: {list_best_EA}') 
            # All thresholds have been evaluated, now do an extra selection among those selected to get the best for EA and CA:
            # If no breakpoints (such as for olorofim) list_best_CA contains all evaluated thresholds, with associated EA
            # list_best_EA contains entries and present_CA will be 0
            # So if no breakpoints, first entry in list_best_EA will be selected as best for EA and secondarily CA
            # And in list_best_CA, threshold with lowest EA will be selected
            best_weighted_CA_for_best_EA = np.inf
            best_thresh_EA = None
            for threshold_CA in list_best_EA:
                if float(threshold_CA.split('_')[1]) < best_weighted_CA_for_best_EA:
                    best_weighted_CA_for_best_EA = float(threshold_CA.split('_')[1])
                    best_thresh_EA = float(threshold_CA.split('_')[0])
            logging.debug(f'--- --- --- --- thresh {best_thresh_EA} has best CA: {best_weighted_CA_for_best_EA}')
            logging.debug(f'--- --- --- For {parameter}, list_best_CA contains {len(list_best_CA)} entries: {list_best_CA}')
            # Select among thresholds with similar CA best performing with EA
            best_EA_for_best_CA = np.inf
            best_thresh_CA = None
            for threshold_EA in list_best_CA:
                if float(threshold_EA.split('_')[1]) < best_EA_for_best_CA:
                    best_EA_for_best_CA = float(threshold_EA.split('_')[1])
                    best_thresh_CA = float(threshold_EA.split('_')[0])
            logging.debug(f'--- --- --- --- thresh {best_thresh_CA} has best EA: {best_EA_for_best_CA}')
        
            # Agreement dict will need to lookup in di including the bias, and will need to store parameter including bias
            # Read data from parameter with best EA primarily and best weighted CA secondarily
            agreement_dict_EA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_EA, di, files_with_data, files_with_breakpoint)
            agreement_dict_CA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_CA, di, files_with_data, files_with_breakpoint)  

        else: # approach 2
            logging.debug(f'--- --- --- --- approach 2')
            # NOTE: among thresholds optimised for EA, we secondarily use min_weighted_errors and NOT CA, because there is no advantage in using an unweighted measure of CA
            # For thresholds optimised for CA, we secondarily use EA - however, this is probably inferior to the weighted categorical measure of approach 1 
            logging.debug(f'--- --- --- For {parameter}, list_best_EA contains {len(list_best_EA)} entries: {list_best_EA}')
            best_weighted_CA_for_best_EA = np.inf # we want to minimise weighted categorical errors
            best_thresh_EA = None
            for threshold_CA in list_best_EA:
                if float(threshold_CA.split('_')[1]) < best_weighted_CA_for_best_EA:
                    best_weighted_CA_for_best_EA = float(threshold_CA.split('_')[1])
                    best_thresh_EA = float(threshold_CA.split('_')[0])
            logging.debug(f'--- --- --- --- thresh {best_thresh_EA} has best CA: {best_weighted_CA_for_best_EA}')
            logging.debug(f'--- --- --- For {parameter}, list_best_CA contains {len(list_best_CA)} entries: {list_best_CA}')
            # Select among thresholds with similar CA best performing with EA
            best_EA_for_best_CA = 0 # we want to MAXIMISE EA
            best_thresh_CA = None
            for threshold_EA in list_best_CA:
                if float(threshold_EA.split('_')[1]) > best_EA_for_best_CA:
                    best_EA_for_best_CA = float(threshold_EA.split('_')[1])
                    best_thresh_CA = float(threshold_EA.split('_')[0])
            logging.debug(f'--- --- --- --- thresh {best_thresh_CA} has best EA: {best_EA_for_best_CA}')

            # Agreement dict will need to lookup in di including the bias, and will need to store parameter including bias
            # Read data from parameter with best EA primarily and best weighted CA secondarily
            if best_thresh_EA is not None:
                agreement_dict_EA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_EA, di, files_with_data, files_with_breakpoint)
            else:
                logging.warning(f'No valid threshold found for EA optimization for parameter {parameter}')
                agreement_dict_EA = {}
            
            if best_thresh_CA is not None:
                agreement_dict_CA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_CA, di, files_with_data, files_with_breakpoint)
            else:
                logging.warning(f'No valid threshold found for CA optimization for parameter {parameter}')
                agreement_dict_CA = {}
        
        # if skip_flag is True, it means in one of the recursive functions we had threshold that resulted in performance that matched
        # it means we can skip the present function call and all parent calls
        if skip_flag:
            logging.info(f'skip_flag is True for target {criterium} - final_target is {final_target} - move on')
            # final_target will be criterium that was successful
            return final_target, True

        for key, value in agreement_dict_EA.items():
            di[antimycotic][timepoint][f'{parameter}_best_EA'][key] = value
        for key, value in agreement_dict_CA.items():
            di[antimycotic][timepoint][f'{parameter}_best_CA'][key] = value
        #logging.info(f"parameters best threshold: {di[antimycotic][timepoint][f'{parameter}_best_EA'].items()}")

    with open(os.path.join(output_dir, f"performance_dict_{antimycotic}_{timepoint}_{criterium}.pkl"), "wb") as f:
        logging.info(f'saving performance dict at:')
        logging.info(os.path.join(output_dir, f"performance_dict_{antimycotic}_{timepoint}_{criterium}.pkl"))
        pickle.dump(performance_dict, f)

    logging.info(f'will return criterium {criterium} (final target is {final_target}) and skip flag True')
    return criterium, True

def get_best_parameter(timepoint, antimycotic, parameters, EA_CA, criterium, di, output_dir, di_path, args_dict, session_time, bias):#, steps, window, pattern):    
    # Now for each parameter an optimal EA and CA threshold has been determined, check which parameter is best for EA and CA
    logging.info(f'--- --- --- Determining best parameter among {len(parameters)} parameter(s) for criterion {EA_CA} and {criterium}')

    # create a dictionary to contain the best thresholds of all parameters 
    ranking_dict = {}
    for parameter in parameters:
        ranking_dict[parameter] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']
    
    df = pd.DataFrame.from_dict(ranking_dict).T
    #logging.info(f'columns {df.columns}: df\n{df}')
    # approach 2        
    if criterium == "EA":
        priority_list = (['EA', 'total_MIC_distance', 'min_errors_weighted', 'CA'], [False, True, True, False])
    elif criterium == "CA":
        priority_list = (['CA', 'min_errors_weighted', 'total_MIC_distance', 'EA'], [False, True, True, False])
    # approach 1
    elif criterium == "min_errors_weighted":
        priority_list = (['min_errors_weighted', 'CA', 'total_MIC_distance', 'EA'], [True, False, True, False])
    elif criterium == "total_MIC_distance":
        priority_list = (['total_MIC_distance', 'EA', 'min_errors_weighted', 'CA'], [True, False, True, False])
        priority_list = (['total_MIC_distance'], [True])
    # approach 3
    elif type(criterium) == float:
        if EA_CA == "EA":
            priority_list = (['min_errors_weighted', 'VME', 'ATU_VME', 'ME', 'ATU_ME', 'CA', 'total_MIC_distance', 'EA'], [True, True, True, True, True, False, True, False])
        elif EA_CA == "CA":
            priority_list = (['min_errors_weighted', 'VME', 'ATU_VME', 'ME', 'ATU_ME', 'total_MIC_distance', 'EA', 'CA'], [True, True, True, True, True, True, False, False])
    logging.info(f'EA_CA {EA_CA}, criterium {criterium}')

    # sort_list determines the order by which we sort the dataframe that contains for all parameters the best performing ones
    logging.info(f'--- --- --- --- Sorting parameters and thresholds - main criterium {criterium}, priority list {priority_list}')
    df.sort_values(priority_list[0], ascending=priority_list[1], inplace=True)
    logging.info(f"parameters df\n{df[['threshold', 'EA', 'CA', 'VME', 'ATU_VME', 'ME', 'ATU_ME', 'min_errors_weighted', 'total_MIC_distance']]}\n")
    parameter = df.index[0]
    threshold = df.iloc[0]['threshold']
    logging.info(f'--- --- --- Best parameter is {parameter}, best threshold {threshold} (bias correction is {bias})')

    best_par_dict = {}
    best_par_dict['time'] = timepoint
    best_par_dict['antimycotic'] = antimycotic
    best_par_dict['parameter'] = parameter
    best_par_dict['criterium'] = criterium
    best_par_dict['threshold'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['threshold']
    best_par_dict['EA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['EA']
    best_par_dict['CA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['CA']
    best_par_dict['VME_rate'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['VME_rate']
    best_par_dict['ME_rate'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['ME_rate']
    best_par_dict['ATU_VME_rate'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['ATU_VME_rate']
    best_par_dict['ATU_ME_rate'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['ATU_ME_rate']
    best_par_dict['CE_per_file'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['CE_per_file']
    best_par_dict['dist_per_file'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['distance_per_file']
    best_par_dict['d_EA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['d_EA']
    best_par_dict['f_EA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['f_EA']
    best_par_dict['d_CA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['d_CA']
    best_par_dict['f_CA'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['f_CA']
    best_par_dict['output_dir'] = output_dir
    best_par_dict['di_path'] = di_path
    
    # Fix: Set total_errors based on the actual optimization criterion, not always EA_CA
    if criterium == "EA" or criterium == "CA":
        # For direct EA/CA optimization, use the optimized metric
        best_par_dict['total_errors'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium]
        logging.info(f'Setting total_errors to {criterium} value: {best_par_dict["total_errors"]}')
    elif criterium == "min_errors_weighted":
        # For weighted errors optimization, use min_errors_weighted
        best_par_dict['total_errors'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['min_errors_weighted']
        logging.info(f'Setting total_errors to min_errors_weighted value: {best_par_dict["total_errors"]}')
    elif criterium == "total_MIC_distance":
        # For MIC distance optimization, use total_MIC_distance
        best_par_dict['total_errors'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['total_MIC_distance']
        logging.info(f'Setting total_errors to total_MIC_distance value: {best_par_dict["total_errors"]}')
    elif type(criterium) == float:
        # For threshold-based optimization, use min_errors_weighted as the primary optimization metric
        best_par_dict['total_errors'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}']['min_errors_weighted']
        logging.info(f'Setting total_errors to min_errors_weighted value for threshold-based optimization: {best_par_dict["total_errors"]}')
    else:
        # Fallback to EA_CA if criterium is not recognized
        best_par_dict['total_errors'] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][EA_CA]
        logging.warning(f'Unrecognized criterium {criterium}, falling back to {EA_CA} value: {best_par_dict["total_errors"]}')
    
    #best_par_dict['bias'] = bias
    #best_par_dict['steps'] = steps
    #best_par_dict['window'] = window
    #best_par_dict['pattern'] = pattern
    best_par_dict['session_time'] = session_time
    best_par_dict.update(args_dict)
    logging.info(f'best_par_dict: {best_par_dict}')

    logging.info(f'--- --- --- Saving sorted dataframe.')
    with open(os.path.join(output_dir, f"all_par_{antimycotic}_{timepoint}_{criterium}.pkl"), "wb") as f:
        pickle.dump(df, f)
    # plot parameters clustered by substring
    #plot.par_ranking(sorted_dict, EA_CA, criterium, output_dir, timepoint, antimycotic) 

    logging.debug(f'--- --- --- --- best_par_dict: {best_par_dict}')
    threshold_dict = {}
    threshold_dict[f'{antimycotic}_{timepoint}'] = f'{best_par_dict["parameter"]}_{bias}_{best_par_dict["threshold"]}' 
    threshold_dict_path = os.path.join(output_dir, f"{antimycotic}_{timepoint}_{criterium}_threshold_dict.pkl")
    with open(threshold_dict_path, "wb") as f:
        pickle.dump(threshold_dict, f)
    logging.debug(f'--- --- --- --- Saved threshold_dict in {output_dir}')

    return best_par_dict, threshold_dict_path
