import logging
import numpy as np
import os
import pickle
from . import plot

def get_agreement_parameters_short(antimycotic, timepoint, parameter, threshold, di):
    logging.debug(f'--- --- --- --- Getting MIC-distance and weighted CA for {threshold}')
    agreement_dict = {}
    # Now all MIC_distances and errors have been summed across all files per threshold, check which threshold is the best performing, both for EA and CA
    # Because VME is worse than ME, calculate weighted sum of number of errors
    # weighted_CA will be np.NaN in case no breakpoints were available
    agreement_dict['MIC_distance'] = di[antimycotic][timepoint][parameter][threshold][0]
    VME_tot = di[antimycotic][timepoint][parameter][threshold][1]
    ME_tot = di[antimycotic][timepoint][parameter][threshold][2]
    ATU_VME_tot = di[antimycotic][timepoint][parameter][threshold][3] 
    ATU_ME_tot = di[antimycotic][timepoint][parameter][threshold][4]
    weighting_factor = 2 # VMEs are counted as double errors - also use factor for ATU_VMEs?
    agreement_dict['weighted_CA'] = (VME_tot * weighting_factor) + ME_tot + (ATU_VME_tot) + ATU_ME_tot
    agreement_dict['unweighted_CA'] = VME_tot + ME_tot + ATU_VME_tot + ATU_ME_tot
    # If prediction failed for a sample for a threshold, this is counted as an ME; see determine_CA() with distance 12
    logging.debug(f'--- --- --- --- MIC distance is {agreement_dict["MIC_distance"]}; weighted_CA is {agreement_dict["weighted_CA"]}')
    
    # In case no breakpoints are available, all values become 0
    return agreement_dict

def get_agreement_parameters_ext(antimycotic, timepoint, parameter, threshold, di, files_with_data, files_with_breakpoint):
    logging.debug(f'--- --- Getting extended agreement parameters for {threshold}')

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
    agreement_dict['min_distance'] = di[antimycotic][timepoint][parameter][threshold][0]
    VME_tot = di[antimycotic][timepoint][parameter][threshold][1]
    ME_tot = di[antimycotic][timepoint][parameter][threshold][2]
    ATU_VME_tot = di[antimycotic][timepoint][parameter][threshold][3] 
    ATU_ME_tot = di[antimycotic][timepoint][parameter][threshold][4]
    agreement_dict['VME'] = VME_tot
    agreement_dict['ME'] = ME_tot
    agreement_dict['ATU_VME'] = ATU_VME_tot
    agreement_dict['ATU_ME'] = ATU_ME_tot
    weighting_factor = 2 # VMEs are counted as double errors - also use factor for ATU_VMEs?
    agreement_dict['min_errors_weighted'] = (VME_tot * weighting_factor) + ME_tot + (ATU_VME_tot) + ATU_ME_tot
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

    return agreement_dict


def get_best_thresholds_bis(timepoint, antimycotic, parameters, di, thresholds, files_with_data, files_with_breakpoint):
    logging.info(f'--- --- --- Determining best threshold for {len(parameters)} parameter(s) for antimycotic {antimycotic} and timepoint {timepoint}')
    # We have calculated the distance between pattern and MIC per antimycotic and per file individually for each threshold
    # Now we need to sum these distances across all files
    # Initialize new hierarchy for storing sum of errors per antimycotic, timepoint, parameter and threshold
    di[antimycotic][timepoint] = {}
    for parameter in parameters:
        logging.debug(f'\n--- --- Par {parameter}')
        # Initialize per antimycotic, timepoint, parameter and threshold, a list that contains [total MIC distance, #VME, #ME, #ATU_VME, #ATU_ME]
        di[antimycotic][timepoint][parameter] = {}
        for threshold in thresholds[timepoint][antimycotic][parameter]:
            di[antimycotic][timepoint][parameter][threshold] = [0, 0, 0, 0, 0]
        # Initialize per antimycotic, timepoint and parameter a dict that will contain performance of best threshold for that parameter
        di[antimycotic][timepoint][f'{parameter}_best_EA'] = {}
        di[antimycotic][timepoint][f'{parameter}_best_CA'] = {}

        # Among thresholds,sum errors over all files, and then select threshold with lowest errors; 
        # list_best_EA/CA contains best thresholds for this timepoint and parameter
        minimum_EA = np.inf # minimum_EA and _CA store the lowest essential/categorical error rate across all thresholds
        minimum_CA = np.inf
        list_best_EA = []
        list_best_CA = []

        for threshold in thresholds[timepoint][antimycotic][parameter]:
            logging.debug(f'--- --- --- Counting errors across {len(files_with_data)} files for threshold {threshold}')
            # Determine across all files the MIC_distance and number of errors for each threshold: store in list of five items
            # Additional level: di[antimycotic][timepoint][parameter][threshold][bias][1]: for bias in bias_list (-2, -1, 0, +1, +2)
            for file in files_with_data:
                di[antimycotic][timepoint][parameter][threshold][0] += abs(di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['MIC_distance'])
                if file in files_with_breakpoint:
                    # if no breakpoint, these values remain 0
                    di[antimycotic][timepoint][parameter][threshold][1] += (di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['VME'])
                    di[antimycotic][timepoint][parameter][threshold][2] += (di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ME'])
                    di[antimycotic][timepoint][parameter][threshold][3] += (di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_VME'])
                    di[antimycotic][timepoint][parameter][threshold][4] += (di[antimycotic][file]['MIC'][parameter][timepoint][threshold]['ATU_ME'])                   
            # For this threshold, first determine total MIC distance and weighted categorical agreement
            agreement_dict = get_agreement_parameters_short(antimycotic, timepoint, parameter, threshold, di)            
            # For EA, create di[antimycotic][timepoint]["best-threshold-per-parameter_EA"] entry that describes summed MIC distances and multiple-/not_found-errors
            present_EA = agreement_dict['MIC_distance']
            present_CA = agreement_dict['weighted_CA']
            # in case no breakpoints available, such as for Olorofim, present_CA will be 0 - present_CA is the sum of all errors (not categorical agreement!)
            # Adapt {threshold} string to include bias: {threshold_bias}
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
                # If present_EA at least equal to minimum_EA: add to list_best
                # No need to reset minimum_EA
                list_best_CA.append(f'{threshold}_{present_EA}')
                logging.debug(f'--- --- --- --- --- Added this threshold to list_best_CA: {list_best_CA}')
                if present_CA < minimum_CA:
                    # If in addition it is lower than minimum_EA, reset minimum_EA, reset list_best and append first entry
                    list_best_CA = []
                    list_best_CA.append(f'{threshold}_{present_EA}')
                    logging.debug(f'--- --- --- --- --- --- New best threshold. New minimum: {present_CA}. Reset list_best_CA: {list_best_CA}')
                    minimum_CA = present_CA

            
        logging.debug(f'--- --- --- For {parameter}, list_best_EA contains {len(list_best_EA)} entries: {list_best_EA}')
 
        # if bias added to string, you'll have to select last item after split in order to get performance EA or CA

        # All thresholds have been evaluated, now do an extra selection among those selected to get the best for EA and CA:
        # If no breakpoints (such as for olorofim) list_best_CA contains all evaluated thresholds, with associated EA
        # list_best_EA contains entries and present_CA will be 0
        # So if no breakpoints, first entry in list_best_EA will be selected as best for EA and secondarily CA
        # And in list_best_CA, threshold with lowest EA will be selected
        best_weighted_CA_for_best_EA = np.inf
        best_thresh_EA = None
        for threshold_CA in list_best_EA:
            if int(threshold_CA.split('_')[1]) < best_weighted_CA_for_best_EA:
                best_weighted_CA_for_best_EA = int(threshold_CA.split('_')[1])
                best_thresh_EA = float(threshold_CA.split('_')[0])
        logging.debug(f'--- --- --- --- thresh {best_thresh_EA} has best CA: {best_weighted_CA_for_best_EA}')
        logging.debug(f'--- --- --- For {parameter}, list_best_CA contains {len(list_best_CA)} entries: {list_best_CA}')
        # Select among thresholds with similar CA best performing with EA
        best_EA_for_best_CA = np.inf
        best_thresh_CA = None
        for threshold_EA in list_best_CA:
            if int(threshold_EA.split('_')[1]) < best_EA_for_best_CA:
                best_EA_for_best_CA = int(threshold_EA.split('_')[1])
                best_thresh_CA = float(threshold_EA.split('_')[0])
        logging.debug(f'--- --- --- --- thresh {best_thresh_CA} has best EA: {best_EA_for_best_CA}')
        
        # Agreement dict will need to lookup in di including the bias, and will need to store parameter including bias
        # Read data from parameter with best EA primarily and best weighted CA secondarily
        agreement_dict_EA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_EA, di, files_with_data, files_with_breakpoint)
        agreement_dict_CA = get_agreement_parameters_ext(antimycotic, timepoint, parameter, best_thresh_CA, di, files_with_data, files_with_breakpoint)  
            
        for key, value in agreement_dict_EA.items():
            di[antimycotic][timepoint][f'{parameter}_best_EA'][key] = value
        for key, value in agreement_dict_CA.items():
            di[antimycotic][timepoint][f'{parameter}_best_CA'][key] = value
        #logging.info(f"parameters best threshold: {di[antimycotic][timepoint][f'{parameter}_best_EA'].items()}")

def get_best_parameter(timepoint, antimycotic, parameters, EA_CA, criterium, di, output_dir, di_path, bias):    
    # Now for each parameter an optimal EA and CA threshold has been determined, check which parameter is best for EA and CA
    minimum = np.inf
    best_par_dict = {}
    logging.info(f'--- --- --- Determining best parameter among {len(parameters)} parameter(s) for criterion {EA_CA} and {criterium}')
    ranking_dict = {}
    for parameter in parameters:
        ranking_dict[parameter] = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium]
        logging.debug(f"--- --- --- --- Par {parameter}: {di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium]}")
    logging.info(f'--- --- --- --- Sorting parameters and thresholds according to performance ...')
    sorted_dict = sorted(ranking_dict.items(), key=lambda x: x[1])
         
        if di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium] < minimum:
            logging.debug(f"--- --- --- --- --- New minimum: par {parameter} has {di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium]}, previous min {minimum}")
            minimum = di[antimycotic][timepoint][f'{parameter}_best_{EA_CA}'][criterium]
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
            best_par_dict['total_dist'] = minimum
            best_par_dict['bias'] = bias

    with open(os.path.join(output_dir, f"all_par_{antimycotic}_{timepoint}_{criterium}.pkl"), "wb") as f:
        pickle.dump(sorted_dict, f)
    # plot parameters clustered by substring
    plot.par_ranking(sorted_dict, EA_CA, criterium, output_dir, timepoint, antimycotic) 

    logging.debug(f'--- --- --- --- best_par_dict: {best_par_dict}')
    threshold_dict = {}
    threshold_dict[f'{antimycotic}_{timepoint}'] = f'{best_par_dict["parameter"]}_{bias}_{best_par_dict["threshold"]}' 
    threshold_dict_path = os.path.join(output_dir, f"{antimycotic}_{timepoint}_{criterium}_threshold_dict.pkl")
    with open(threshold_dict_path, "wb") as f:
        pickle.dump(threshold_dict, f)
    logging.debug(f'--- --- --- --- Saved threshold_dict in {output_dir}')

    return best_par_dict, threshold_dict_path
