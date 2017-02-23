from scipy.interpolate import interp1d
import numpy as np

# read position z and B at 1000A
x, y = np.genfromtxt('data/B_z_1000A.txt', unpack=True)
    
x = (x + 0.64) * 100
y = y * 1e4  # in gauss
# now we have the field at 1000A 

def getB(current):
    y_scaled = y *  (current/ 1000)
    #returning the interpolation funct
    return interp1d(x, y_scaled)
    

