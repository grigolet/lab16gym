import pandas as pd
import numpy as np
import math

def read_data(path, use_cols, col_names):
    '''
    This function reads the data from a txt file.
    It needs list of the columns to read and their names
    '''
    #reading the data fram
    df = pd.read_csv(path, header=None, sep="\s+", usecols=use_cols, dtype="float64")
    #assigning columns' names
    df.columns = col_names
    return df

def remove_offset(df, offset):
    '''This function remove offset from the res
    signal BEFORE amplification. '''
    df.v_res = df.v_res - offset

def amplify_data(df, ampl):
    '''
    This function receives a list of tuples with
    columns to amplify and the factor
    [("v_scheda",7.43),('v_r',2)].
    It created new columns in the dataframe
    appending the _ampl string to the name.'''
    for a in ampl:
        df[a[0]+'_ampl'] = df[a[0]]*a[1]

def calculate_data(df, R):
    '''This function calculate and add to dataframe
    the current and v_probe given the resistance R'''
    df['v_probe'] = df.v_scheda_ampl - df.v_res_ampl
    df['i_raw'] = (df.v_res_ampl / R)


def capacity_correction(df, cap):
    I_cap = (np.gradient(df.v_probe)/np.gradient(df.t)) * cap
    df['i'] = df['i_raw'] + I_cap


def bin_data(df, intervals, n_bin):
    '''This funcions bins the data in zones and bins.
    -Intervals is a list of the extremes of the zones.
    -n_bin is a list the number of bin in each zone.
    It returns a dictionary with groups created
    and the lists of means and stds of every column in df,
    for each zone'''
    zones = []
    #reading zones
    for a in range(len(intervals)-1):
        #np.linspace generates bins and .append group them in a list
        zones.append(np.linspace(intervals[a],
                                 intervals[a+1],
                                 n_bin[a]))
    #creating data for each zone
    groups = []
    means = []
    stds = []
    #binning for each zone
    for a in range(len(intervals)-1):
        group = df.groupby(pd.cut(df["v_probe"],zones[a])) #split the group
        groups.append(group) #add the group to the groups list
        means.append(group.mean()) #the same
        stds.append(group.std())
    return {"groups": groups, #dictionary of the groups
            "means": means,
            "stds" : stds}

def concat_data(data):
    '''This function concatenate the data of different zones
    into x,y,xerr,yerr. It returns a numpy tuple with (x,y,xerr,yerr)'''
    #Here I'm getting the means and stds for every zone and
    #concatenating them!
    x = np.concatenate([m.v_probe.values for m in data['means']])
    #for each zone get the array values and concatenate them
    y = np.concatenate([m.i.values for m in data['means']])
    xerr = np.concatenate([s.v_probe.values for s in data['stds']])
    yerr = np.concatenate([s.i.values for s in data['stds']])
    return (x,y,xerr,yerr)


def get_plot_data(df, intervals, n_bin):
    '''This funciont gets a dataframe and the intervals
    definition and returns the groups information
    and plots data (x,y,xerr,yerr). It simply
    unifies the previous two functions.'''
    bins_data = bin_data(df, intervals, n_bin)
    plt_data = concat_data(bins_data)
    return (plt_data, bins_data)
