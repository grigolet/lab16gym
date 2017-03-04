import numpy as np
import matplotlib.pyplot as plt
import time
import math
from scipy.odr import Model, RealData, ODR

A = 15 * 10**-6
# m_i = (938 + 939) * 10**6
c = 3 * 10**8
e = - 1.60217653 * 10**-19

def fit_initial_param(bins_data, plt_data):
    '''
    This function computes the inizial parameter for the fit
    from the data. It search V_floating, alpha and I ionic saturation.

    It returns a tuple with (V_floating, alpha, i_ionic, figure, output_str).
    -Figure is the plot of the linear fit
    -output_str is a preformatted output string result.
    '''

    #1) Extracting V_floating
    #selecting the data from the right zone (where V_floating is located)
    right_range = bins_data['means'][1]
    #we found the min I value, but we interested
    #in the corresponding x value that is V_floating
    i_floating = right_range['i'].abs().min()
    #searching data row corresponding to i_floating and accessing the v index
    v_floating = right_range[right_range['i'].abs()==i_floating]['v_probe'].values[0]

    #2) Fitting the ionization current
    # We start by defining our fit function
    def ionic_saturation_function(p, v_probe):
        return p[0] + p[1]*(v_probe -v_floating)

    #we select the data from the first range (the left one)
    x_left  = bins_data['means'][0].v_probe.values
    y_left  = bins_data['means'][0].i.values
    sx_left = bins_data['stds'][0].v_probe.values
    sy_left = bins_data['stds'][0].i.values
    #next, we define our type of model
    linear_fit = Model(ionic_saturation_function)
    odr_data   = RealData(x_left, y_left, sx=sx_left, sy=sy_left)
    #beta0 represents and array with the parameter guessing
    odr = ODR(odr_data, linear_fit, beta0=[0.0002, 0.001])
    #running fit and getting result
    output = odr.run()
    #saving the parameter obtained by the fit
    q                 = output.beta[0]
    sigma_q           = output.sd_beta[0]
    alpha             =  - output.beta[1] / q
    #the ionic current is equal to q with this fitting function
    i_ionic_sat       = q


    #creating the figure for linear fitting
    y_fit = q - q*alpha*(x_left -v_floating)
    fig1 ,ax = plt.subplots(1,1)
    ax.plot(x_left, y_fit, 'r', label='Fit curve', linewidth=2.)
    #fmt doesn't plot the line between data points
    ax.errorbar(x_left,y_left, xerr=sx_left,
                yerr=sy_left, fmt='o', label='Data')
    ax.set_xlabel('$ V\;(V) $', fontsize=18)
    ax.set_ylabel('$ I\;(mA) $', fontsize=18)
    ax.legend(bbox_to_anchor=(0.15,1))
    fig1.set_size_inches(18,12)
    ax.grid()

    output = []
    output.append('Ionic saturation current = {0} A'.format(round(i_ionic_sat,5)))
    output.append('Alpha                    = {0}'.format(round(alpha, 5)))
    output.append('V floating               = {0} V'.format(round(v_floating,5)))
    output_str = '\n'.join(output)

    return (v_floating, i_ionic_sat, alpha, fig1, output_str)

def fit_data(plt_data, initial_data, n_iterations, m_i):
    t1 = time.time()
    #unpacking variables
    x, y, xerr, yerr = plt_data
    v_floating, i_ionic_sat, alpha = initial_data
    #splitting data for iterations
    iterations_data = split_data(n_iterations, plt_data, v_floating)
    # we need to save every fit result in order
    # to use it with the next fitting.
    # We define our first guessing parameters as fit_results.
    #We define 5 the initial electronic Temperature
    fit_results = [(i_ionic_sat, alpha, v_floating, 10.)]
    #separate list for the errors results from fitting
    error_results = [(0,0,0,0)]
    #list for reduced chi2 results
    chi2r_results = []

    #iterating o
    for i in range(n_iterations):
        #getting data for i index
        x_fit, y_fit, sx_fit, sy_fit = iterations_data[i]

        #next, we define our type of model
        langmuir_fit   = Model(fpff)
        odr_data   = RealData(x_fit, y_fit, sx=sx_fit, sy=sy_fit)
        #every cycle we use the previous fit
        #result as beta 0 parameters
        odr = ODR(odr_data, langmuir_fit,
                  beta0=fit_results[i], maxit=1000)
        #running fit and getting result
        output = odr.run()
        #saving data
        fit_results.append(tuple(output.beta))
        error_results.append(tuple(output.sd_beta))
        chi2r_results.append(output.res_var)

    #now we have to search for the best fit
    #using the principle of min temperature
    #unpacking results for temperature
    temp_fit = [b[3] for b in fit_results]
    # discard the first temp value because it has been arbitrarly
    # set
    temp_fit = temp_fit[1:]
    iterations_data = iterations_data[1:]
    fit_results = fit_results[1:]
    error_results = error_results[1:]
    chi2r_results = chi2r_results[1:]
    #find the minimum of the temperature to get all associated
    #parameters and use them as the best parameters for fitting
    min_temp = min(temp_fit)
    min_index = temp_fit.index(min_temp)
    last_v_probe = iterations_data[min_index][0][-1]
    best_fit_results =  fit_results[min_index]
    best_fit_errors  = error_results[min_index]
    chi2r_best       = chi2r_results[min_index]

    #creating text output
    output = []
    output.append('Minimum temperature found: {0}'.format(min_temp))
    output.append('Correspondent index of iteration: {0}'.format(min_index))
    output.append('V_probe limit of fit: {} V'.format(last_v_probe))
    output.append('\nParameters associated with the lowest temperature:')
    output.append('Ionic saturation current = {0} ± {1} mA'.format(
            best_fit_results[0], best_fit_errors[0]))
    output.append('Alpha                    =  {0} ± {1}'.format(
            best_fit_results[1], best_fit_errors[1]))
    output.append('V floating               =  {0} ± {1} V'.format(
            best_fit_results[2], best_fit_errors[2]))
    output.append('Temp                     =  {0} ± {1}'.format(
            best_fit_results[3], best_fit_errors[3]))
    output.append('Reduced Chi2             =  {0}'.format(
            round(chi2r_best,3)))
    output_str = '\n'.join(output)

    #Creating plots
    #1) Plot of temperatures
    #plot of temperatures with respect to the number of iterations.
    temp_fig, temp_plot = plt.subplots(1,1)
    temp_plot.grid(True)
    temp_plot.set_xlabel('Number of iterations', fontsize=18)
    temp_plot.set_ylabel('Temperature (eV)', fontsize=18)
    temp_plot.plot(temp_fit,'go', label='Temperature')
    temp_plot.axhline(y= best_fit_results[3],linewidth=2, color = 'r')
    temp_fig.set_size_inches(18,12)

    #2) Final plot of the fit
    plot_text = ('Parameters obtained with the \nminimum'+
        ' temperature fitting: \n\n' +
        'last $V_{probe} = %.3f \, V $\n' +
        '$ I_{i,sat} = %.5f \pm %.5f\, A $\n'  +
        '$\\alpha = %.5f \pm %.5f $\n' +
        '$V_{floating} = %.3f \pm %.3f \,V $\n' +
        '$T_e = %.3f \pm %.3f \,eV $\n' +
        "$\\chi^2 reduced = %.3f $") % (last_v_probe,
            best_fit_results[0], best_fit_errors[0],
            best_fit_results[1], best_fit_errors[1],
            best_fit_results[2], best_fit_errors[2],
            best_fit_results[3], best_fit_errors[3] ,chi2r_best)
    #calculating i from fitted data
    #We use x because we want the i in all the x range
    i_fit = fpff(best_fit_results, x)
    fit_fig, fit_plot = plt.subplots(1,1)
    fit_plot.grid(True)
    fit_plot.plot(x, i_fit, 'r', label='Fit curve', linewidth=2.)
    #fmt is the key to avoid line that connects markers
    fit_plot.errorbar(x,y, xerr=xerr, yerr=yerr,
                fmt='o', label='Data')
    fit_plot.set_xlabel('$ V\;(V) $', fontsize=18)
    fit_plot.set_ylabel('$ I\;(A) $', fontsize=18)
    fit_plot.set_xlim([x[0]-5, x[-1]+5 ])
    fit_plot.set_ylim([-0.003, 0.035])
    fit_plot.legend(bbox_to_anchor=(0.11,1))
    fit_plot.text(-50., 0.02, plot_text, fontsize=18,
                bbox=dict(facecolor='white', edgecolor='black', ))
    fit_plot.axvspan(last_v_probe, x[-1]+5 ,
                facecolor='blue', alpha=0.2)
    fit_fig.set_size_inches(18,12)
    t2 = time.time()
    #calculate addictional data
    ne_data = calculate_ne(
        best_fit_results[0],best_fit_errors[0],
        best_fit_results[3], best_fit_errors[3], m_i)
    vplasma_data = calculate_vplasma(
        best_fit_results[2], best_fit_errors[2],
        best_fit_results[3], best_fit_errors[3])
    f_ep_data = calculate_fep(ne_data[0], ne_data[1])
    #creating dictionary with data
    data = {
        "i_ionic_sat": best_fit_results[0],
        "err_i_ionic_sat" : best_fit_errors[0],
        "alpha" : best_fit_results[1],
        "err_alpha" : best_fit_errors[1],
        "v_floating" : best_fit_results[2],
        "err_v_floating" : best_fit_errors[2],
        "T_e" : best_fit_results[3],
        "err_T_e": best_fit_errors[3],
        "n_e": ne_data[0],
        "err_n_e": ne_data[1],
        "v_plasma": vplasma_data[0],
        "err_v_plasma": vplasma_data[1],
        "f_ep": f_ep_data[0],
        "err_f_ep": f_ep_data[1]}
    #returning a dictionary with all data
    return {
            'best_iteration': min_index,
            "data": data,
            'chi2r_results': chi2r_results,
            'last_v_probe': last_v_probe,
            'output_str': output_str,
            'temp_fig': temp_fig,
            'fit_fig': fit_fig,
            'exec_time': t2-t1
    }


def fpff(p, x):
    '''Our Four Parameters Fitting Function'''
    return p[0]*(1-p[1]*(x-p[2])-np.exp((x-p[2])/p[3]))


def split_data(n, plt_data, v_float):
    '''
    This function split the x,y,xerr,yerr data arrays in blocks.
    The first block is from the beginning to v_float.
    The reamining data is splitted in n blocks.
    Then the blocks are concatenated recursively to create
    consecutive iterations blocks.
    The function return a list of tuples with (x,y,xerr,yerr) arrays
    for every block.
    '''
    x,y,xerr,yerr = plt_data
    start_index = np.where(x==v_float)[0][0]#the right end starts from v_floating

    #defining all the data. Starting from the v_floating point
    x_iter = x[start_index:]
    y_iter = y[start_index:]
    sx_iter = xerr[start_index:]
    sy_iter = yerr[start_index:]
    result = []
    x_blocks = np.array_split(x_iter, n) #split the data in n parts
    y_blocks = np.array_split(y_iter, n)
    sx_blocks = np.array_split(sx_iter, n)
    sy_blocks = np.array_split(sy_iter, n)
    result = [(x[:start_index],y[:start_index], #saving a tuple of intervals where we will fit
              xerr[:start_index], yerr[:start_index])]
    for i in range(n):
        result.append((
                #append to the tuple the i-th frame of data
                np.concatenate([result[i][0],x_blocks[i]]),
                np.concatenate([result[i][1], y_blocks[i]]),
                np.concatenate([result[i][2], sx_blocks[i]]),
                np.concatenate([result[i][3], sy_blocks[i]])
            ))
    return result

def calculate_ne(I, sigma_I, T_e, sigma_Te, m_i):
        '''This function return (n_e, sigma_ne)'''
        n_e = I / (0.6 * e * A * c) *  math.sqrt(m_i * 10**6 / (T_e))
        sigma_ne = ( I/(e*c*A))*math.sqrt(((m_i*sigma_I**2)/(I**2 * T_e)) +\
                    ((m_i* (sigma_Te**2))/(T_e**3)))
        return (n_e,sigma_ne)

def calculate_vplasma(V_floating, sigma_V_floating, Te, sigma_te):
    V_plasma =  V_floating + 3*Te
    sigma_V_plasma = sigma_V_floating + 3 * sigma_te
    return (V_plasma, sigma_V_plasma)

def calculate_fep(n_e, sigma_ne):
    return (8.98 *(n_e)**0.5, (8.98 /(2* (n_e)**(3/2) ) )*sigma_ne  )
