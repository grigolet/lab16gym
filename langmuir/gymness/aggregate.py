import os
import json
import sys
# import pandas as pd
from pandas import DataFrame
import numpy as np

A = 12 * 10**-6

def aggregate_data(path):
    '''
    This function reads all data from a directory
    and puts it in a dataframe.
    '''
    # list for all data
    items = [x for x in os.listdir(path)]
    print("Reading dataset:" + str(items))
    data = []
    # reading all json files
    for item in items:
        data += json.load(open(path + "/" + item + "/fit_results.json", 'r'))
    # creating dataframe
    df = DataFrame(data, columns=['id', 'probe', 'x', 'I_B', 'P', 'W_em', 'P_off',
                        'T_e', 'err_T_e', 'n_e', 'err_n_e', 'f_ep', 'err_f_ep',
                        'v_plasma', 'err_v_plasma', 'alpha', 'err_alpha',
                        'v_floating', 'err_v_floating',
                        'i_ionic_sat', 'err_i_ionic_sat'])
    # post_processing on the whole data
    post_processing(df)
    return df


def post_processing(df):
    df['j_sat'] = df['i_ionic_sat'] / A


def remove_bias_ne(df):
    '''This function removes the bias getting
    the scale in the middle point of the selected columns'''
    j1 = df[(df.x == 0) & (df.probe == 1)]['j_sat'].values[0]
    j2 = df[(df.x == 0) & (df.probe == 2)]['j_sat'].values[0]
    medio = (j2 + j1) / 2
    delta1 = (medio - j1) / j1
    delta2 = (medio - j2) / j2
    df.ix[df.probe == 1, 'n_e'] *= (1 + delta1)
    df.ix[df.probe == 2, 'n_e'] *= (1 + delta2)


def remove_bias(df, columns):
    '''This function removes the bias in the central
    point for the column specified'''
    deltas = []
    for column in columns:
        s1 = df[(df.probe == 1) & (df.x == 0)][column].values[0]
        s2 = df[(df.probe == 2) & (df.x == 0)][column].values[0]
        delta = (s2 - s1) / 2
        df.ix[df.probe == 1, column] += delta
        df.ix[df.probe == 2, column] -= delta
        deltas.append((round(((delta * 2) / s1) * 100, 4),
                       round(((delta * 2) / s2) * 100, 4)))
    return deltas

def get_fuh(df, B):
    '''This funcion adds to the dataframe the f_upperhybrid'''
    f_ce = np.multiply(2.80e6, B)
    df["B"] = B
    df["f_uh"] = (df.f_ep**2 + f_ce**2)**0.5
    df["err_f_uh"] = (df.f_ep / ((df.f_ep**2 + df.B**2)**(3 / 2))) * df.err_f_ep


if __name__ == "__main__":
    # reading from the conf file
    path = sys.argv[1]
    df = aggregate_data(path)
