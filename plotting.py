import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
from statistics import mean, median, stdev
import pandas as pd
import numpy as np

def plot_OD(df, variable, subtract_background = False, errorbars=False, yscale='log', append_title='', save=False, pdf=None):
    '''
    Plots OD measurements from DataFrame that is returned from od_query function.

    Args:
        df (pandas.DataFrame): Dataframe returned from od_query() function.
        subtract_background (bool): Whether to subtract background reading 
            from all measurements.
        yscale (str): 'log' or 'linear'
        append_title (str): Additional text to add to the figure title
    Returns:
        None    
    '''
        
    # Prepare data for plotting
    df = df.sort_values('datetime')
    df['od_background_subtracted'] = df['od'] - df['background']
    
    if subtract_background:
        value = 'od_background_subtracted'
    else:
        value = 'od'

    # Define different combinations of conditions that will be plotted
    conditions = df[[variable, 'strain_name']].drop_duplicates().dropna()
    conditions['label'] = conditions.apply(
        lambda x: f'{x["strain_name"]} - {x[variable]}', axis=1
    )
    conditions['colors'] = colormaps['tab20'].colors[:len(conditions)]

    # For the legend
    handles = []
    labels = []

    # Create a figure and a set of subplots
    fig, ax = plt.subplots()

    # Figure will need to be stretched horizontally for readability
    fig_width, fig_height = fig.get_size_inches() # Get the current figure size
    fig.set_size_inches(fig_width * 3.5, fig_height) # Doubling the width
    
    total_transfers = df['passage'].max()

    # For the defined conditions and at each transfer,
    # plot the OD readings
    for i, row in conditions.iterrows():
        handle_line = Line2D([0], [0], label=row['label'], color=row['colors'])
        handles.append(handle_line)
        label = row['label']
        labels.append(label)
        
        for t in range(1, total_transfers+1):
            this_condition = df.loc[
                (df[variable] == row[variable]) &
                (df['strain_name'] == row['strain_name']) &
                (df['passage'] == t)
            ]
            if errorbars == False:
                
                this_condition = this_condition.groupby('datetime'
                                                       )[value].agg(mean).to_frame().reset_index()
                plt.plot(
                    this_condition['datetime'],
                    this_condition[value],
                    color=row['colors'],
                    marker='o',
                    markersize=2
                )
                
            else:
                y = this_condition.groupby('datetime')[value].agg(mean)
                yerr = this_condition.groupby('datetime')[value].agg(stdev)
                x = this_condition.groupby('datetime')[value].agg(mean).index
                plt.errorbar(
                    x=x,
                    y=y,
                    yerr=yerr,
                    color=row['colors'],
                    marker='o',
                    markersize=2
                )
    
    # Configure and label the axes and tickmarks
    plt.yscale(yscale)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    plt.xlabel('datetime')
    plt.ylabel('OD') 

    # early_datetimes = df.groupby('passage')['datetime'].agg(
    #     lambda x: sorted(list(set(x)))[3])
    # secax = ax.secondary_xaxis('top')
    # secax.set_xticks(early_datetimes, np.arange(1, total_transfers+1, 1))
    # secax.set_xlabel('transfer')

    # Set title and legend
    plt.title(append_title + " " + value, fontsize=16)
    # plt.legend(handles=handles, labels=labels, loc='lower center')
    plt.legend(handles=handles, labels=labels, loc='lower center', bbox_to_anchor=(.5, -.4))

    if save:
        plt.savefig(pdf, format='pdf', bbox_inches='tight')
    
    plt.show()


def plot_OD_separate(df, variable, subtract_background = False, errorbars=False, yscale='log', append_title='', save=False, pdf=None):
    '''  
    '''
        
    # Prepare data for plotting
    df = df.sort_values('datetime')
    df['od_background_subtracted'] = df['od'] - df['background']
    
    if subtract_background:
        value = 'od_background_subtracted'
    else:
        value = 'od'

    # Define different combinations of conditions that will be plotted
    conditions = df[[variable, 'strain_name']].drop_duplicates().dropna()
    conditions['label'] = conditions.apply(
        lambda x: f'{x["strain_name"]} - {x[variable]}', axis=1
    )
    conditions['colors'] = colormaps['tab20'].colors[:len(conditions)]

    # For the legend
    handles = []
    labels = []

    total_transfers = df['passage'].max()
    
    # Create a figure and a set of subplots
    fig, ax = plt.subplots(nrows=1, ncols=total_transfers, sharey=True)

    # Figure will need to be stretched horizontally for readability
    fig_width, fig_height = fig.get_size_inches() # Get the current figure size
    fig.set_size_inches(fig_width * 3.5, fig_height) # Doubling the width

    # For the defined conditions and at each transfer,
    # plot the OD readings
    
    for i, row in conditions.iterrows():
        handle_line = Line2D([0], [0], label=row['label'], color=row['colors'])
        handles.append(handle_line)
        label = row['label']
        labels.append(label)
        
        for t in range(1, total_transfers+1):
            this_condition = df.loc[
                (df[variable] == row[variable]) &
                (df['strain_name'] == row['strain_name']) &
                (df['passage'] == t)
            ]
            if errorbars == False:
                
                this_condition = this_condition.groupby('timepoint'
                                                       )[value].agg(mean).to_frame().reset_index()
                ax[t-1].plot(
                    this_condition['timepoint'],
                    this_condition[value],
                    color=row['colors'],
                    marker='o',
                    markersize=2
                )
                
            else:
                y = this_condition.groupby('timepoint')[value].agg(mean)
                yerr = this_condition.groupby('timepoint')[value].agg(stdev)
                x = this_condition.groupby('timepoint')[value].agg(mean).index
                ax[t-1].errorbar(
                    x=x,
                    y=y,
                    yerr=yerr,
                    color=row['colors'],
                    marker='o',
                    markersize=2
                )
            
    # Configure and label the axes and tickmarks
    plt.yscale(yscale)
    # ax.yaxis.set_major_formatter(ScalarFormatter())
    
    fig.supxlabel('timepoint')
    fig.supylabel('OD')
    plt.subplots_adjust(wspace=0)

    # early_timepoints = df.groupby('passage')['datetime'].agg(
    #     lambda x: sorted(list(set(x)))[3])
    # secax = ax.secondary_xaxis('top')
    # secax.set_xticks(early_datetimes, np.arange(1, total_transfers+1, 1))
    # secax.set_xlabel('transfer')

    # Set title and legend
    fig.suptitle(append_title + " " + value, fontsize=16)
    # plt.legend(handles=handles, labels=labels, loc='lower center')
    fig.legend(handles=handles, labels=labels, loc='lower center', bbox_to_anchor=(.5, -.4))
    
    if save:
        plt.savefig(pdf, format='pdf', bbox_inches='tight')
    
    plt.show()
    
def plot_growth_metric(df, value, variable, errorbars=False, append_title=''):
    '''
    Plots growthrates from DataFrame that is returned from query_growth_rate function.

    Args:
        df (pandas.DataFrame): Dataframe returned from od_query() function.
        append_title (str): Additional text to add to the figure title
    Returns:
        None    
    '''

    # Define different combinations of conditions that will be plotted
    conditions = df[[variable, 'strain_name']].drop_duplicates().dropna()
    conditions['label'] = conditions.apply(
        lambda x: f'{x["strain_name"]} - {x[variable]}', axis=1
    )
    conditions['colors'] = colormaps['tab20'].colors[:len(conditions)]

    # For the legend
    handles = []
    labels = []

    # Create a figure and a set of subplots
    fig, ax = plt.subplots()

    # Figure will need to be stretched horizontally for readability
    fig_width, fig_height = fig.get_size_inches() # Get the current figure size
    fig.set_size_inches(fig_width * 3.5, fig_height) # Doubling the width


    # For the defined conditions and at each passage,
    # plot the growth metric
    for i, row in conditions.iterrows():
        handle_line = Line2D([0], [0], label=row['label'], color=row['colors'])
        handles.append(handle_line)
        label = row['label']
        labels.append(label)
        
        this_condition = df.loc[
            (df[variable] == row[variable]) &
            (df['strain_name'] == row['strain_name'])
        ]
        
        if errorbars == False:

            plt.plot(
                this_condition.groupby('passage')['passage'].agg(mean),
                this_condition.groupby('passage')[value].agg(mean),
                color=row['colors'],
                marker='o',
                markersize=2
            )
            
        else:
            y = this_condition.groupby('passage')[value].agg(mean)
            yerr = this_condition.groupby('passage')[value].agg(stdev)
            x = this_condition.groupby('passage')['passage'].agg(mean)
            plt.errorbar(
                x=x,
                y=y,
                yerr=yerr,
                color=row['colors'],
                marker='o',
                markersize=2
            )
        
    # # Configure and label the axes and tickmarks
    # plt.yscale(yscale)
    # ax.yaxis.set_major_formatter(ScalarFormatter())
    plt.xlabel('transfer')
    plt.ylabel(value) 

    # early_datetimes = df.groupby('passage')['datetime'].agg(
    #     lambda x: sorted(list(set(x)))[3])
    # secax = ax.secondary_xaxis('top')
    # secax.set_xticks(early_datetimes, np.arange(1, total_transfers+1, 1))
    # secax.set_xlabel('transfer')

    # Set title and legend
    plt.title(append_title + " " + value, fontsize=16)
    # plt.legend(handles=handles, labels=labels, loc='lower center')
    plt.legend(handles=handles, labels=labels, loc='lower center', bbox_to_anchor=(.27, -.5, .47, .102))
    
    plt.show()