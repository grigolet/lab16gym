# importing our lib for binning and fitting
import os
import json 
import sys
import multiprocessing
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import binning
import fitting
from scipy import stats  # need for chi2 statistics


def plot_raw_data(df, path):
    #plotting raw data
    fig, ax = plt.subplots(1, 1)
    ax.errorbar(df['v_probe'], df['i'], fmt='o')
    ax.set_xlabel("V (Volt)")
    ax.set_ylabel("I probe (A)")
    ax.grid(True)
    fig.set_size_inches(18, 12)
    fig.savefig(path + "/raw_data.png")


def plot_binned_data(plt_data, path):
    fig, ax = plt.subplots(1, 1)
    ax.grid(True)
    ax.plot(plt_data[0], plt_data[1], ':b')
    ax.errorbar(x=plt_data[0],
                y=plt_data[1],
                xerr=plt_data[2],
                yerr=plt_data[3])
    ax.set_xlabel("V (Volt)")
    ax.set_ylabel("I probe (A)")
    fig.set_size_inches(18, 12)
    fig.savefig(path + "/binned_data.png")



def elaborate_dataset(fit_param, print_fig, input_dir, output_dir, capacity_1,
                      capacity_2,offset_1,offset_2, R_1, R_2, probes, columns):
    fit_results = []
    ############## Fitting FIRST probe
    if 1 in probes:
        print('\n####### Fitting ' + fit_param['id'] + " - Probe 1 ######")
        fit_path = output_dir + "/" + fit_param['id'] + "/probe1"
        os.makedirs(fit_path, exist_ok=True)
        file_path = input_dir + "/" + fit_param['id'] + ".txt"
        #read the first three columns
        df = binning.read_data(file_path, columns[1],
                               ['t', 'v_res', 'v_scheda'])
        #remoing offset_1
        binning.remove_offset(df, offset_1)
        #We have to change sign of V_r for first probe
        binning.amplify_data(df, [('v_scheda', fit_param['ampl_scheda']),
                                  ('v_res', fit_param['ampl_res'])])
        #compute the current and v_probe
        #I_gain is the correction for the first probe
        binning.calculate_data(df, R_1)
        print('Loaded dataframe')
        #calculate capacity correction
        binning.capacity_correction(df, capacity_1)

        #elaborating data into bins and getting plot data
        plt_data, bins_data = binning.get_plot_data(df, fit_param['zones'],
                                    fit_param['bins'])

        #fitting initial data
        initial_data = fitting.fit_initial_param(bins_data, plt_data)

        output_file = open(fit_path + "/output.txt", 'w')
        output_file.write('Analysis params:'+\
                          '\nampl_scheda: {}'.format(fit_param['ampl_scheda'])+\
                          '\nampl_res: {}'.format(fit_param['ampl_res'])+\
                          '\nzones: {}'.format(fit_param['zones'])+\
                          '\nbins: {}'.format(fit_param['bins'])+\
                          '\nN. iterations: {}'.format(fit_param['n_fit_iter']))
        output_file.write('\n\nInitial params:\n' + initial_data[4])

        #fitting all data
        N_iterations = fit_param['n_fit_iter']
        data_fitted = fitting.fit_data(plt_data, initial_data[:3],
                                       N_iterations)
        #printing results
        print('Fit exec time: {}'.format(round(data_fitted['exec_time'], 5)))
        output_file.write('\n\nFitted data:'+\
              '\nN. iterations: {}'.format(N_iterations) +'\n'+\
                data_fitted['output_str'] +\
                '\n\nExec time: {} s'.format(
                round(data_fitted['exec_time'],5)))
        output_file.close()
        #saving results for exporting
        fit_output = {
            "id": fit_param['id'],
            "probe": 1,
            "W_em": fit_param['W_em'],
            "P": fit_param['P'],
            "I_B": fit_param['I_B'],
            "x": fit_param['x'] - 10  #the short probe is
            #10 cm before the long probe that is the reference
        }
        fit_output.update(data_fitted['data'])
        fit_results.append(fit_output)

        if print_fig:
            #plotting raw data
            plot_raw_data(df, fit_path)
            #plotting binned data
            plot_binned_data(plt_data, fit_path)
            #plotting fit results
            initial_data[3].savefig(fit_path + "/fit_i_ionic.png")
            data_fitted['temp_fig'].savefig(fit_path + "/temp_plot.png")
            data_fitted['fit_fig'].savefig(fit_path + "/fit_plot.png")

    if 2 in probes:
        ############## Fitting SECOND probe
        print('\n####### Fitting ' + fit_param['id'] + " - Probe 2 ######")
        fit_path = output_dir + "/" + fit_param['id'] + "/probe2"
        os.makedirs(fit_path, exist_ok=True)
        file_path = input_dir + "/" + fit_param['id'] + ".txt"
        #read the first three columns
        df = binning.read_data(file_path, columns[2],
                             ['t', 'v_res', 'v_scheda'])
        #remoing offset_2
        binning.remove_offset(df, offset_2)
        #the signs are correct
        binning.amplify_data(df, [('v_scheda', fit_param['ampl_scheda']),
                                  ('v_res', fit_param['ampl_res'])])
        #compute the current and v_probe
        binning.calculate_data(df, R_2)
        print('Loaded dataframe\n')
        #calculate capacity correction
        binning.capacity_correction(df, capacity_2)

        #elaborating data into bins and getting plot data
        plt_data, bins_data = binning.get_plot_data(df, fit_param['zones'],
                                        fit_param['bins'])

        #fitting initial data
        initial_data = fitting.fit_initial_param(bins_data, plt_data)

        output_file = open(fit_path + "/output.txt", 'w')
        output_file.write('Analysis params:'+\
                          '\nampl_scheda: {}'.format(fit_param['ampl_scheda'])+\
                          '\nampl_res: {}'.format(fit_param['ampl_res'])+\
                          '\nzones: {}'.format(fit_param['zones'])+\
                          '\nbins: {}'.format(fit_param['bins'])+\
                          '\nN. iterations: {}'.format(fit_param['n_fit_iter']))
        output_file.write('\n\nInitial params:\n' + initial_data[4])

        #fitting all data
        N_iterations = fit_param['n_fit_iter']
        data_fitted = fitting.fit_data(plt_data, initial_data[:3],
                                       N_iterations)
        #printing results
        print('Fit exec time: {}'.format(round(data_fitted['exec_time'], 5)))
        output_file.write('\n\nFitted data:'+\
              '\nN. iterations: {}'.format(N_iterations) +'\n'+\
                data_fitted['output_str'] +\
                '\n\nExec time: {} s'.format(
                round(data_fitted['exec_time'],5)))
        output_file.close()
        #saving results for exporting
        fit_output = {
            "id": fit_param['id'],
            "probe": 2,
            "W_em": fit_param['W_em'],
            "P": fit_param['P'],
            "I_B": fit_param['I_B'],
            "x": fit_param['x']  #the second probe is the reference for x
        }
        fit_output.update(data_fitted['data'])
        fit_results.append(fit_output)
        if print_fig:
            #plotting raw data
            plot_raw_data(df, fit_path)
            #plotting binned data
            plot_binned_data(plt_data, fit_path)
            #plotting fit results
            initial_data[3].savefig(fit_path + "/fit_i_ionic.png")
            data_fitted['temp_fig'].savefig(fit_path + "/temp_plot.png")
            data_fitted['fit_fig'].savefig(fit_path + "/fit_plot.png")

    #closing all the figures
    plt.close('all')
    #printing fit results
    json.dump(fit_results, open(output_dir+ "/"+ \
                        fit_param['id']+"/fit_results.json",'w'),
                        indent=4)


if __name__ == "__main__":
    #reading config url
    f = sys.argv[1]
    print_fig = sys.argv[2] == "y"

    #reading configs
    conf = json.load(open(f))
    input_dir = conf['input_dir']
    output_dir = conf['output_dir']
    os.makedirs(output_dir, exist_ok=True)
    capacity_1 = conf['capacity_1']
    capacity_2 = conf['capacity_2']
    offset_1 = conf['offset_1']
    offset_2 = conf['offset_2']
    R_1 = conf['R_1']
    R_2 = conf['R_2']
    probes = conf['probes']
    columns = {1: conf['columns1'],
               2: conf['columns2'] }
    fit_param_generic = conf.get('params')


    for fit_param in conf['fits']:
        #adding the generic param
        for p in fit_param_generic:
            if not p in fit_param:
                fit_param[p] = fit_param_generic[p]
        '''
        p = multiprocessing.Process(
            target=elaborate_dataset,
            args=(fit_param, print_fig, input_dir, output_dir, capacity_1,
                  capacity_2, offset_1, offset_2, R_1, R_2, probes,columns))
        p.start()
        '''
        elaborate_dataset(fit_param, print_fig, input_dir, output_dir, capacity_1,
                  capacity_2, offset_1, offset_2, R_1, R_2, probes,columns)