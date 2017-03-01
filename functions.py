# -*- coding: utf-8 -*-
import numpy as np
import peakutils
import json
from scipy import interpolate


def get_peak_ratio(file_name):
    wavelenght, intensity = np.genfromtxt(file_name, unpack=True)

    peakind = peakutils.indexes(intensity, thres=0.1, min_dist=0.5)
    
    x_fit=np.zeros(21)
    y_fit=np.zeros(21)

    for i in range(0,21):
        x_fit[i] = wavelenght[peakind[np.abs(wavelenght[peakind] - 480).argmin()] + i - 10]
        y_fit[i] = intensity[peakind[np.abs(wavelenght[peakind] - 480).argmin()] + i - 10]
    I480, x480, s480 = peakutils.peak.gaussian_fit(x_fit, y_fit, center_only=False)
            
    for i in range(0,21):
        x_fit[i] = wavelenght[peakind[np.abs(wavelenght[peakind] - 488).argmin()] + i - 10]
        y_fit[i] = intensity[peakind[np.abs(wavelenght[peakind] - 488).argmin()] + i - 10]
    I488, x488, s488 = peakutils.peak.gaussian_fit(x_fit, y_fit, center_only=False)

    ratio = I480 / I488
    ratio_err = ratio*(s480/I480 + s488/I488)

    return ratio

        
# definisco la funzione irradianza
def irradiance(x):
    A = 43.1162499557604
    B = -4456.27408561124
    C = 0.958841675471763
    D = 128.862702098102
    E = -127865.310911202
    F = 49796783.0638663
    G = -7243608027.29156
    H = 0
    return np.power(x,-5)*np.exp(A+B/x)*(C+D/x+E/np.power(x,2)+F/np.power(x,3)+G/np.power(x,4)+H/np.power(x,5))


#definisco L_lambda=radianza di diffusione
def diffuser_radiance(x, rho):
    return np.sqrt(2)/np.pi*irradiance(x)*rho   


#definisco Intensity assoluta corretta (viene 1/J)
def I_correct_abs(I, I_lamp, rho, x, x_lamp):
    L=2e8                       #nm, cosi va via con i nm di x
    h_c=1.9864e-25                  #J
    diff_r=diffuser_radiance(x, rho)
    I_lamp_correct=np.interp(x, x_lamp, I_lamp)
    I_abs=(I*diff_r)/I_lamp_correct
    return I_abs*(x/(h_c))*((4*np.pi)/L)


def get_absolute_intensity(I_lamp_file, rho_file, I_argon_file, acq_time):
    #leggo il file riferito alla I_lamp e normalizzo I_lamp rispetto al tempo (3sec di acquisizione)
    x_lamp, I_lamp=np.genfromtxt(I_lamp_file, unpack=True)
    I_lamp=I_lamp/3

    #leggo il file riferito alla rho=riflettivita di diffusione a 45gradi
    x_rho, rho=np.genfromtxt(rho_file, unpack=True, skip_header=2, delimiter=',')

    #leggo il file riferito alla I_argon e la normalizzo a un secondo (file di dati presi ogni 0.01s)
    x_I, I_argon=np.genfromtxt(I_argon_file, unpack=True)
    I_argon=I_argon/acq_time

    #creo corrispondenza tra i lambda del file (x_rho, rho) e I_argon e I_correct attraverso interpolazione
    I_argon_correct=np.interp(x_rho, x_I, I_argon)
    I_correct=I_correct_abs(I_argon_correct, I_lamp, rho, x_rho, x_lamp)

    #troviamo la temperatura
    n_e = 2e11              #cm-3
    # usando la legge dei gas perfetti
    n_0 = ((495)/(1.60e-19*2.57e-2))*10e-6  #cm-3 (P=6.38e-5/1.29 mbar=495 Pa, Kx_B=1.38e-23/8.62e-5=1.60e19 J/eV, T=T_amb=2.57e-2eV, 10e-6=fattore per udm
    X750 = I_correct[abs(x_rho - 750).argmin()]/(n_0*n_e)*10e-4 #trova il valore di I + vicino al lamdba=750
    I_abs_750 = I_correct[abs(x_rho - 750).argmin()]*10e-4

    return [np.array([x_rho, I_correct]), n_0, X750, I_abs_750]


def get_profile_ratio(config_file, data_path):
    """This file reads the config file, stored as a .json and the
    data path where the data is stored and returns the following:
    an array of average intensity line ratios,
    the power and 
    a dict with the results grouped by power"""
    # open the config file
    with open(config_file) as config_file_data:
        config_file = json.load(config_file_data)

    ratios = {}
    results = {}
    ratio_mean = np.array([])
    w_em = np.array([])

    # iterate over the items in teh config file to 
    # get the powers and file names
    for key, value in config_file.items():
        peak_ratio = get_peak_ratio(data_path + key + '.dat')
        if value['w'] not in results.keys():
            results[value['w']] = []
        if value['w'] not in ratios.keys():
            ratios[value['w']] = []

        results[value['w']].append({
            'id': key,
            'ratio': peak_ratio,
            'p': value['p']
        })

        ratios[value['w']].append(peak_ratio)

    for w, values_ratio in ratios.items():
        ratio_mean = np.append(ratio_mean, np.mean(values_ratio))
        w_em =np.append(w_em, w)

    return ratio_mean, w_em, results


def get_density(peak_ratio, fants_file='fants.txt'):
    """
    Take an array of line ratios (I480/I488) and interpolate the fants.txt file
    in order to get the corresponding density
    """
    
    if (not isinstance (peak_ratio, np.ndarray)):
        # read data
        n, r_1, r_2, r_3, r_4, r_5, r_6, r_7 = np.genfromtxt(fants_file, unpack=True)
        #ratio = np.array([r_1, r_2, r_3, r_4, r_5, r_6, r_7])
        lista = [r_1, r_2, r_3, r_4, r_5, r_6, r_7]
        n_e = np.empty([7])
        i=0
        for r_i in lista:
            f_interp = interpolate.interp1d(r_i, n)
            np.vectorize(f_interp)
            n_e[i] = f_interp(peak_ratio)
            i = i + 1
        n_e_mean = np.mean(n_e)
        n_e_min = min(n_e)
        n_e_max = max(n_e)
        return n_e_mean, n_e_min, n_e_max

    """
    Per capire le modifiche nel ciclo sotto secondo me è meglio guardare il ciclo sopra, che è quello "vecchio" con la sola aggiunta
    della dondizione. Se il peak_ratio che prende in ingresso la funzione è un vettore (cosa che avviene in in peak_analysis.ipynb 
    nella cella immediatamente successiva a quella che stampa i risultati, quando operiamo su tutti i file specificati nel json) quando 
    chiamo f_interp(peak_ratio) nel ciclo for il risultato non è un numero, ma un array di numeri, uno per ogni file analizzato; 
    cerco quindi di mettere in n_e[i] un vettore, da cui l'errore "setting array element with a sequence". Per risolvere il problema
    se il peak_ratio in ingresso è un vettore faccio eseguire il ciclo sopra per ogni elemento di peak_ratio; 
    il ciclo sopra ritorna tre numeri, quello sotto tre array di numeri;
    """               
                                                
    if (isinstance (peak_ratio, np.ndarray)):
        n, r_1, r_2, r_3, r_4, r_5, r_6, r_7 = np.genfromtxt(fants_file, unpack=True)
        lista = [r_1, r_2, r_3, r_4, r_5, r_6, r_7]
        n_e = np.empty([7])
        l = len(peak_ratio)
        n_e_mean_array = np.empty([l])
        n_e_min_array = np.empty([l])
        n_e_max_array = np.empty([l])
        
        for j in range (0,l):
            i=0
            for r_i in lista:
                f_interp = interpolate.interp1d(r_i, n)
                np.vectorize(f_interp)
                n_e[i] = f_interp(peak_ratio[j])
                i = i + 1
            n_e_mean_array[j] = np.mean(n_e)
            n_e_min_array[j] = min(n_e)
            n_e_max_array[j] = max(n_e)   
        
        return n_e_mean_array, n_e_min_array, n_e_max_array