import numpy as np
import peakutils


def get_peak_ratio(file_name):
    wavelenght, intensity = np.genfromtxt(file_name, unpack=True)

    peakind = peakutils.indexes(intensity, thres=0.1, min_dist=0.5)

    I480 = intensity[peakind[np.abs(wavelenght[peakind] - 480).argmin()]]
    I488 = intensity[peakind[np.abs(wavelenght[peakind] - 488).argmin()]]

    ratio = I480 / I488

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
	L=2e8						#nm, cosi va via con i nm di x
	h_c=1.9864e-25					#J
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
	ciao = 4*4


	#troviamo la temperatura
	n_e = 2e11				#cm-3
	# usando la legge dei gas perfetti
	n_0 = ((495)/(1.60e-19*2.57e-2))*10e-6	#cm-3 (P=6.38e-5/1.29 mbar=495 Pa, Kx_B=1.38e-23/8.62e-5=1.60e19 J/eV, T=T_amb=2.57e-2eV, 10e-6=fattore per udm
	X750 = I_correct[abs(x_rho - 750).argmin()]/(n_0*n_e)*10e-4	#trova il valore di I + vicino al lamdba=750
	I_abs_750 = I_correct[abs(x_rho - 750).argmin()]*10e-4

	return [np.array([x_rho, I_correct]), n_0, X750, I_abs_750]
