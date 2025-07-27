import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
import itertools
from datetime import datetime
import os
import seaborn as sns 
from matplotlib import colors
from matplotlib import figure 
import pprint
import pandas as pd

def plot(**kwargs):
    parameters = kwargs.get('params')
    antimycotics = kwargs.get('am')
    samples = kwargs.get('files')
    timeframe = kwargs.get('timeframe')
    well_order = kwargs.get('well_legend')
    twin_parameters = kwargs.get('par_twin')
    thresh_y = kwargs.get('thresh_y')
    thresh_dict = kwargs.get('thresh_dict')
    y_lim = kwargs.get('y_lim')
    y_lim_twin = kwargs.get('y_lim_twin')
    timepoint = kwargs.get('timepoint')
    label = kwargs.get('label')
    dist_per_file = kwargs.get('highlight')
    d = kwargs.get('d')
    di = kwargs.get('di')
    output_path = kwargs.get('output_path')
    criterium = kwargs.get('criterium')

    if not twin_parameters: twin_parameters = [] # twin_parameters should be an empty list and not just None, to avoid zipping problems
    ## still to define globally instead of here
    files_pDST = samples
    ##    

    if samples[0] in files_pDST: # Check first file in samples: if it is a training sample, construct ref_wells for plotting wells relative to MIC in right order
        # well_order_dict: used for creating titles per subplot
        well_order_dict = {}
        well_order_dict['MIC'] = 'MIC'
        well_order_dict['GC'] = 'GC'
        well_order.remove('MIC') # well_order is a copy of ref_wells
        well_order.remove('GC')
        for i in well_order:
            well_order_dict[i] = '{} log{}'.format(i.split('MIC')[1], '\u2082')
        well_order.insert(int((len(well_order)) / 2), 'MIC') # re-insert 'MIC' in the middle
        well_order.append('GC') # append GC at end
        # If needed, decrease no. of plots per row: cut off 'crop' number of subplots from both sides and re-append 'GC' at the end
        crop = 4
        well_order_cropped = well_order[(crop - 1):(len(well_order) - crop)]
        well_order_cropped.append('GC')
        logging.debug(f'>>> >>> >>> Well order: \n {well_order} \n well order cropped: \n {well_order_cropped}')        

        # linestyle dict will make MIC and GC solid linestyles, wells with lower log dilutions dashed, and containing higher log dilutions dotted
        linestyle_dict = {}
        for i in well_order:
            if (i == 'MIC') or (i == 'GC'):
                linestyle_dict[i] = 'solid'
            elif (int(i[3:]) > 0):
                linestyle_dict[i] = 'dotted'
            else:
                linestyle_dict[i] = 'dashed'
    else:
        logging.info(f'No phenotypic DST data available for {samples}')
        well_order = list(d[antimycotics[0]][samples[0]].keys()) # Get the well names/order from the keys of the first sample of the first antimycotic in d
        well_order_cropped = well_order[1:]
        linestyle_dict = {}
        well_order_dict = {}
        for i in well_order_cropped:
            well_order_dict[i] = i
            linestyle_dict[i] = 'solid'


    # For plotting of antimycotics per row, change len(samples) to len(antimycotics)
    fig, axs = plt.subplots(len(samples)+1, len(well_order_cropped), sharey = True, sharex = True, squeeze=False, constrained_layout=True)           
    fig.suptitle(f'Antimycotic: {antimycotics}, parameter: {parameters}, time {timepoint}, criterium {criterium}, threshold {thresh_y}')
    # For plotting of antimycotics per row, change len(samples) to len(antimycotics)
    fig.set_size_inches(30, len(samples) * 3.5)
    twin_axes = []
    
    for par, twin_par in itertools.zip_longest(parameters, twin_parameters):
        # create highlight dictionary that per file contains predicted well number
        # dist_per_file: highlight[file][parameter][threshold] contains distance
        if dist_per_file:
            logging.debug(f'>>> >>> >>> dist_per_file: {dist_per_file}')
            for antimycotic in antimycotics: 
                highlight = {}
                for file, par_thresh in dist_per_file.items():
                    highlight[file] = int((len(well_order_cropped) - 2) / 2) - par_thresh[par][thresh_y]
                    #logging.debug(f'--- \n well order cropped: {well_order_cropped} \n length: {len(well_order_cropped)} \n total {int((len(well_order_cropped) - 2) / 2)} \n distance {distance}')
        logging.debug(f'>>> >>> >>> Highlight dictionary: \n {highlight.items()}')

        row = 0
        # For plotting antimycotics per row instead of samples, put antimycotics on highest level instead of samples
        for sample in samples:
            for antimycotic in antimycotics:   
                axs[row, 0].set_ylabel(sample)
                twin_axes = []
                logging.info(f'>>> >>> >>> Plotting {antimycotic} for {sample}')
                for column, well in enumerate(well_order_cropped):
                    # Fetch data for plotting
                    y = d[antimycotic][sample][well][par]
                    x = d[sample]
                    if (y.size == 0): 
                        y = [np.NaN]
                        x = [np.NaN]
                    elif thresh_y: # only plot threshold line in subplots where data is
                        if thresh_dict: # in case threshold line is different per file, fetch from dictionary 
                            try: y_thresh_line = thresh_dict[sample][par][thresh_y]
                            except KeyError: logging.debug(f'>>> >>> >>> >>> No data available for file {key} at this timepoint') 
                        else:
                            y_thresh_line = thresh_y
                        axs[row, column].plot((timeframe[0], timepoint), (y_thresh_line, y_thresh_line), linestyle='dotted', color='grey')
                        axs[row, column].plot(timepoint, y_thresh_line, '+', color='black')

                    axs[row, column].plot(x, y, linestyle=linestyle_dict[well], label=f"{sample.split(' ')[0]} {par}")
                    # y_lim can be user-supplied tuple of custom y_lim values - other option is y_lim as a list of multiple values
                    # if thresh_dict is non-empty (proportional thresholds), for each file a custom threshold is used: all thresholds are read and interval between min/max used for setting y limits
                    # if thresh_dict is empty (absolute thresholds), the interval between all evenly spread thresholds is used for setting y limits
                    # thresh_dict contains per file and parameter the relative thresholds set according to proportional thresholds in threshold_dict: thresh_dict[file]['BCA'][0.12] = 1.83
                    # threshold_dict contains per antimycotic the proportional thresholds; threshold_dict['A'] = [0.1, 0.11, 0.12 ...]
                    # thresholds is the equivalent of thresh_dict for absolute thresholds: thresholds[timepoint][antimycotic][parameter] = [1.83, 1.94, ...]
                    logging.debug(f'>>> >>> >>> >>> y_lim is {y_lim}; type {type(y_lim)}')
                    # for proportional thresholds y_lim will be 
                    # Setting y limits: if tuple provided as y_lim: use this to set window
                    if y_lim and isinstance(y_lim, tuple):
                        logging.debug(f'>>> >>> >>> >>> >>> y_lim provided as tuple.')
                        axs[row, column].set_ylim(y_lim)
                    # Alternatively, y_lim is a list of threshold values in case absolute thresholds are used: these are evenly spaced and the interval between values is used to define ideal window for plotting y values
                    # In case relative thresholds are used, y_lim is a list of thresholds of format y_limk['A'] = [0.02, 0.03, 0.04]
                    elif y_lim and isinstance(y_lim[antimycotic], list) and len(y_lim[antimycotic]) > 1:
                        range_list = []
                        if thresh_dict: # add all values for this threshold across all samples to range_list, sort
                            for key, value in thresh_dict.items():
                                try: range_list.append(value[par][thresh_y])
                                except KeyError: logging.debug(f'>>> >>> >>> >>> No data available for file {key} at this timepoint')
                            range_list.sort()
                            interval = range_list[-1] - range_list[0] 
                            bottom = range_list[0] - 5 * interval
                            top = range_list[-1] + 20 * interval
                        else:
                            interval = y_lim[antimycotic][1] - y_lim[antimycotic][0] # if not thresh_dict, all values are evenly spread, take interval between first and second value
                            bottom = y_lim[antimycotic][0] - 5 * interval
                            top = y_lim[antimycotic][-1] + 5 * interval
                        axs[row, column].set_ylim(bottom, top)
                        logging.debug(f'>>> >>> >>> >>> >>> List of thresholds provided: y_lim will be set according to min/max threshold values: {interval}, {bottom} to {top}')  
                    
                    axs[row, column].set_xlim(timeframe)
                    axs[row, column].title.set_text(f"{well_order_dict[well]}")
                    handles, labels = axs[0, int((len(well_order_cropped)) / 2) - 1].get_legend_handles_labels()

                    # For plotting second parameter on twin y-axis
                    if twin_par:
                        y_twin = d[antimycotic][sample][well][twin_par]
                        if (y_twin.size == 0): 
                            y_twin = [np.NaN]
                            x = [np.NaN]
                        ax_twin = axs[row, column].twinx()
                        ax_twin.plot(x, y_twin, linestyle=linestyle_dict[well], label=f"{sample.split(' ')[0]} {twin_par}", color='red')
                        if y_lim_twin: 
                            ax_twin.set_ylim(y_lim_twin)
                        twin_axes.append(ax_twin)
                        twin_axes[0].get_shared_y_axes().join(*twin_axes) # https://stackoverflow.com/questions/53831482/share-secondary-y-axis-in-looped-seaborn-plots?noredirect=1&lq=1
                        for twin_ax in twin_axes[:-1]: # Only have ticks at rightmost subplot
                            twin_ax.yaxis.set_tick_params(labelright=False)
                        handles, labels = [(a + b) for a, b in zip(axs[0, int((len(well_order_cropped)) / 2) - 1].get_legend_handles_labels(), ax_twin.get_legend_handles_labels())] #https://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
                    
                    # For assigning colors to subplots according to whether threshold was passed or not
                    pass_threshold = None
                    try: 
                        pass_threshold = di[antimycotic][sample][well][par][timepoint][thresh_y]['bool']
                        logging.debug(f'>>> >>> >>> >>> threshold = {round(thresh_y, 2)}, pass: {pass_threshold}') 
                    except KeyError:
                        logging.warning('KeyError: growth/no growth entry not found in di dict')
                    if pass_threshold:
                        axs[row, column].set_facecolor('lightpink')
                    elif pass_threshold == False:
                        axs[row, column].set_facecolor('honeydew')
    
            axs[row, (len(well_order_cropped) - 1)].set_facecolor('gainsboro')
            if samples[0] in files_pDST:
                for spine in axs[row, int((len(well_order_cropped) - 1) / 2)].spines.values():
                    spine.set_edgecolor('grey')
                [x.set_linewidth(1) for x in axs[row, int((len(well_order_cropped) - 1) / 2)].spines.values()]
            else:
                [x.set_linewidth(1) for x in axs[row, int(MIC_dict[sample][antimycotic][1:])].spines.values()]
                for spine in axs[row, int(MIC_dict[sample][antimycotic][1:])].spines.values():
                    spine.set_edgecolor('grey')
            if dist_per_file:
               #logging.debug(f'>>> >>> >>> >>> dist per file: \n {dist_per_file}, \n highlight sample: {highlight.items()}')
                try:
                    if not highlight[sample] <= -1: # this will be -1 when no cropping, it will be more negative if plotting is cropped
                        #axs[row, highlight[sample]].set_facecolor('lemonchiffon')
                        axs[row, highlight[sample]].annotate("*", xy=(0.5, 0.7), xycoords="axes fraction", ha="center", va="center", fontsize=30)
                except Exception as e:
                    logging.warning(f'Error {e} - Could not highlight subplot in row {row}; {file} not found in dictionary: {highlight}')
            row += 1
    #fig.legend(handles, labels, loc='upper left')
    #fig.legend(loc='upper left')
    #plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0)
    #beautiful_dict = pprint.pformat(label)
    #plt.figtext(0.5, 0.01, f"{beautiful_dict}", ha="center")
    logging.info(f'label {label}')
    if label:
        CE_per_file = pd.DataFrame(label['CE_per_file'], index=[0]).transpose()
        dist_per_file = pd.DataFrame(label['dist_per_file'], index=[0]).transpose()
        tab = pd.concat([CE_per_file, dist_per_file], axis=0,  keys=['CE_per_file', 'dist_per_file'])
        tab = pd.DataFrame(tab[0])
        tab.reset_index(level=[0, 1], inplace=True, names=['Parameter', 'Sample'])
        plt.subplot(len(samples)+1, 1, len(samples)+1)
        no_rows = tab.shape[0]
        # bbox: first coordinate is a shift on the x-axis, second coordinate is a gap between plot and text box (table in your case), third coordinate is a width of the text box, fourth coordinate is a height of text box.
        table = plt.table(cellText=tab.values, colLabels=tab.columns, loc='center', cellLoc='right', bbox=[0, -no_rows*0.8, 0.4, 3])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        plt.axis('off')

    now = datetime.now().strftime('%H-%M-%S')
    logging.debug(f'output path {output_path}, antimycotics {antimycotics}, timepoint {timepoint}, par {par}, thresh {thresh_y}, now {now}')
    plt.savefig(os.path.join(output_path, f'{antimycotics[0]}_{timepoint}_{criterium}_{par}_thresh_{round(thresh_y, 2)}_{now}.png'), bbox_inches="tight")
    logging.info(f'>>> >>> >>> File saved at {antimycotics[0]}_{timepoint}_{criterium}_{par}_thresh_{round(thresh_y, 2)}_{now}.png')
    plt.close()

def get_color_dict(df):
    # Assign a color to parameters, clustered per substring
    substrings = ['BCANormalized', 'SESAfungiNormalized', 'TANormalized', 'BCA', 'SESAfungi', 'TA']
    color_list = list(sns.color_palette("hls", 8))
    color_dict = {}
    for substring in substrings:
        color_dict[substring] = color_list[0]
        color_list.pop(0)
    colors_plot = {}
    parameter_names = df['parameter'].unique()
    for parameter in parameter_names:
        for substring, color in color_dict.items():
            if substring in parameter:
                colors_plot[parameter] = color
                break
    return colors_plot

def plot_overview(df, a, t, p, s, output_path, label, copy_path):

    # Extract unique criteria from the dataframe to determine what to plot
    unique_criteria = df['criterium'].unique() if 'criterium' in df.columns else []
    
    # Initialize dataframes for EA and CA data
    df_EA = pd.DataFrame()
    df_CA = pd.DataFrame()
    
    # Determine which criteria are present and filter data accordingly
    if 'EA' in unique_criteria:
        df_EA = df.loc[(df['criterium'] == 'EA')]
    elif 'min_distance' in unique_criteria:
        df_EA = df.loc[(df['criterium'] == 'min_distance')]
    
    if 'CA' in unique_criteria:
        df_CA = df.loc[(df['criterium'] == 'CA')]
    elif 'min_errors_weighted' in unique_criteria:
        df_CA = df.loc[(df['criterium'] == 'min_errors_weighted')]
    
    # Also check for the 'c' column which might contain the criteria information
    if 'c' in df.columns:
        # Extract EA and CA data based on the 'c' column
        df_EA_from_c = df.loc[(df['c'].str.contains('EA', na=False))]
        df_CA_from_c = df.loc[(df['c'].str.contains('CA', na=False))]
        
        # Also include single configuration results for both EA and CA plotting
        df_single_config = df.loc[(df['c'].str.contains('single_config', na=False))]
        if not df_single_config.empty:
            # Reset indices to avoid alignment issues when concatenating
            df_EA_from_c = df_EA_from_c.reset_index(drop=True)
            df_CA_from_c = df_CA_from_c.reset_index(drop=True)
            df_single_config = df_single_config.reset_index(drop=True)
            df_EA_from_c = pd.concat([df_EA_from_c, df_single_config], ignore_index=True)
            df_CA_from_c = pd.concat([df_CA_from_c, df_single_config], ignore_index=True)
        
        # Use the data from 'c' column if the criterium-based filtering didn't work
        if df_EA.empty and not df_EA_from_c.empty:
            df_EA = df_EA_from_c
        if df_CA.empty and not df_CA_from_c.empty:
            df_CA = df_CA_from_c

    # Handle the case where phase column might not exist - create filters specific to each subset
    if not df_EA.empty and 'phase' in df_EA.columns:
        phase_filter_EA = df_EA['phase'] != 'pre_bias_corr'
        df_EA_without_pre = df_EA.loc[phase_filter_EA]
    elif not df_EA.empty:
        df_EA_without_pre = df_EA
    else:
        df_EA_without_pre = pd.DataFrame()
    
    if not df_CA.empty and 'phase' in df_CA.columns:
        phase_filter_CA = df_CA['phase'] != 'pre_bias_corr'
        df_CA_without_pre = df_CA.loc[phase_filter_CA]
    elif not df_CA.empty:
        df_CA_without_pre = df_CA
    else:
        df_CA_without_pre = pd.DataFrame()
    
    # Handle best_shift_corr filtering - check each subset individually
    if not df_EA.empty and 'best_shift_corr' in df_EA.columns:
        df_EA_per = df_EA.loc[(df_EA['best_shift_corr'] == 'best')][['EA', 'CA', 'VME_rate', 'ME_rate', 'time', 'ATU_VME_rate', 'ATU_ME_rate']]
    elif not df_EA_without_pre.empty:
        df_EA_per = df_EA_without_pre[['EA', 'CA', 'VME_rate', 'ME_rate', 'time', 'ATU_VME_rate', 'ATU_ME_rate']]
    else:
        df_EA_per = pd.DataFrame()
    
    if not df_CA.empty and 'best_shift_corr' in df_CA.columns:
        df_CA_per = df_CA.loc[(df_CA['best_shift_corr'] == 'best')][['EA', 'CA', 'VME_rate', 'ME_rate', 'time', 'ATU_VME_rate', 'ATU_ME_rate']]
    elif not df_CA_without_pre.empty:
        df_CA_per = df_CA_without_pre[['EA', 'CA', 'VME_rate', 'ME_rate', 'time', 'ATU_VME_rate', 'ATU_ME_rate']]
    else:
        df_CA_per = pd.DataFrame()
    
    dfm_EA = df_EA_per.melt('time', var_name='cols', value_name='vals') if not df_EA_per.empty else pd.DataFrame()
    dfm_CA = df_CA_per.melt('time', var_name='cols', value_name='vals') if not df_CA_per.empty else pd.DataFrame()

    timepoints = df["time"].unique()
    fig, axes = plt.subplots(2, 3, sharey=True)
    fig.set_size_inches(max(5 + len(timepoints) * 0.6, 17), max(len(timepoints) * 0.4, 10))
    fig.suptitle(f'Antimycotic: {a}, {t} timepoints, {p} parameters, {s} steps')
    #fig.autofmt_xdate()

    colors_plot = get_color_dict(df)
    
    # Get the actual criteria names for the titles
    ea_criterion = df_EA['criterium'].iloc[0] if not df_EA.empty and 'criterium' in df_EA.columns else 'EA'
    ca_criterion = df_CA['criterium'].iloc[0] if not df_CA.empty and 'criterium' in df_CA.columns else 'CA'
    
    if not df_EA.empty:
        scatter_ess_EA = sns.scatterplot(ax=axes[0,0], x='time', y='EA', data=df_EA, hue='parameter', style="phase", palette=colors_plot)
        timepoints = df_EA['time'].unique()
        for point in timepoints:
            min_ = df_EA.loc[(df_EA['time']) == point]['EA'].min()
            max_ = df_EA.loc[(df_EA['time']) == point]['EA'].max()
            scatter_ess_EA.plot((point, point), (min_, max_), color='whitesmoke', zorder=0)
        axes[0,0].set_title(f'EA; {ea_criterion}')
        scatter_ess_EA.plot((timepoints.min(), timepoints.max()), (0.9, 0.9), linestyle = 'dotted', color='grey', zorder=0)
        box = axes[0,0].get_position()
        axes[0,0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        axes[0,0].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=5)
        
        scatter_ess_CA = sns.scatterplot(ax=axes[0,1], x='time', y='CA', data=df_EA, hue='parameter', style="phase", palette=colors_plot)
        for point in timepoints:
            min_ = df_EA.loc[(df_EA['time']) == point]['CA'].min()
            max_ = df_EA.loc[(df_EA['time']) == point]['CA'].max()
            scatter_ess_CA.plot((point, point), (min_, max_), color='whitesmoke', zorder=0)
        axes[0,1].set_title(f'CA; {ea_criterion}')
        scatter_ess_CA.plot((timepoints.min(), timepoints.max()), (0.9, 0.9), linestyle = 'dotted', color='grey', zorder=0)
        box = axes[0,1].get_position()
        axes[0,1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        axes[0,1].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=5)       
        
        sns.pointplot (ax=axes[0,2], x="time", y="vals", hue='cols', data=dfm_EA, scale=0.4, errorbar=None, palette=sns.color_palette("hls", 8))
        axes[0,2].set_title(f'Error rates; {ea_criterion}')
        axes[0,2].set_xticklabels(labels = timepoints, rotation=90, fontsize=8)

    if not df_CA.empty:
        scatter_err_EA = sns.scatterplot(ax=axes[1,0], x='time', y='EA', data=df_CA, hue='parameter', style="phase",  palette=colors_plot)
        axes[1,0].set_title(f'EA; {ca_criterion}')
        timepoints = df_CA['time'].unique()
        for point in timepoints:
            min_ = df_CA.loc[(df_CA['time']) == point]['EA'].min()
            max_ = df_CA.loc[(df_CA['time']) == point]['EA'].max()
            scatter_err_EA.plot((point, point), (min_, max_), color='whitesmoke', zorder=0)        
        scatter_err_EA.plot((timepoints.min(), timepoints.max()), (0.9, 0.9), linestyle = 'dotted', color='grey', zorder=0)
        box = axes[1,0].get_position()
        axes[1,0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        axes[1,0].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=5)
        
        axes[1,1].set_title(f'CA; {ca_criterion}')
        scatter_err_CA = sns.scatterplot(ax=axes[1,1], x='time', y='CA', data=df_CA, hue='parameter', style="phase", palette=colors_plot)
        for point in timepoints:
            min_ = df_CA.loc[(df_CA['time']) == point]['CA'].min()
            max_ = df_CA.loc[(df_CA['time']) == point]['CA'].max()
            scatter_err_CA.plot((point, point), (min_, max_), color='whitesmoke', zorder=0)
        scatter_err_CA.plot((timepoints.min(), timepoints.max()), (0.9, 0.9), linestyle = 'dotted', color='grey', zorder=0)
        box = axes[1,1].get_position()
        axes[1,1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
        axes[1,1].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=5)       
         
        pointplot = sns.pointplot(ax=axes[1,2], x="time", y="vals", hue='cols', data=dfm_CA, scale=0.4, errorbar=None, palette=sns.color_palette("hls", 8))
        axes[1,2].set_title(f'Error rates; {ca_criterion}')
        axes[1,2].set_xticklabels(labels = timepoints, rotation=90, fontsize=8)

    now = datetime.now().strftime('%d-%m-%Y-%H-%M-%S')
    plt.savefig(os.path.join(output_path, f'{label}.png'), dpi=400)
    plt.savefig(os.path.join(copy_path, f'{label}_{now}.png'), dpi=400)
    plt.close()

def par_ranking(sorted_dict, EA_CA, criterium, output_dir, timepoint, antimycotic):
    parameter_names = [x[0] for x in sorted_dict]
    parameter_values = [x[1] for x in sorted_dict]

    # Assign a color to parameters, clustered per substring
    substrings = ['BCANormalized', 'SESAfungiNormalized', 'TANormalized', 'BCA', 'SESAfungi', 'TA'] 
    color_list = list(colors.TABLEAU_COLORS)
    color_list = list(sns.color_palette("hls", 8))
    color_dict = {}

    for substring in substrings:
        color_dict[substring] = color_list[0]
        color_list.pop(0)    
    colors_plot = []

    for parameter in parameter_names:
        for substring, color in color_dict.items():
            if substring in parameter:
                colors_plot.append(color)
                break
        else:
            colors_plot.append('gray')

    # Determine the number of labels for adjusting figure width
    num_labels = len(parameter_names)
    figure_width = num_labels * 0.6
    fig = figure.Figure(figsize=(figure_width, 20))
    ax = fig.add_subplot(111)
    ax.bar(parameter_names, parameter_values, color=colors_plot)

    ax.set_title(f'Parameters with best {criterium}')
    ax.set_xticks(range(num_labels))
    ax.set_xticklabels(parameter_names, rotation=90)
    fig.savefig(os.path.join(output_dir, f'all_par_{antimycotic}_{timepoint}_{criterium}.png'))
    plt.close()
    logging.info(f'--- --- --- Saved parameter ranking chart at {output_dir}')
