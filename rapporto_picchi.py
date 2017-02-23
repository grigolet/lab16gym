import numpy as np
#from scipy import signal
import peakutils

deme_directory = "/home/demetrio/Dropbox/Laboratorio_Plasmi/lab16gym/OES_Studenti_2017/"
file = "Ar15022017_a1_acq_time_0.5s.txt"

wavelenght, intensity = np.genfromtxt(deme_directory + file, unpack = True)

#peakind = signal.find_peaks_cwt(intensity, np.arange(1,5,0.01))
peakind = peakutils.indexes(intensity, thres = 0.1, min_dist = 0.5)

I480 = intensity[peakind[np.abs(wavelenght[peakind] - 480).argmin()]]
I488 = intensity[peakind[np.abs(wavelenght[peakind] - 488).argmin()]]

ratio = I480/I488