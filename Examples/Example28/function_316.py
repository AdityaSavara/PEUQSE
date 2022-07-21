# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:27:52 2022

@author: NWT

This file contains the simulated function for modeling the concentration of Cr
in a sample relative to depth. A wrapper function which takes all parameters as
a single nested array is also implemented, as PEUQSE requires to run properly
"""

import numpy as np
import scipy.special as sp
import sys
import data_316

# This simulation function takes a single independant variable, depth (x) and returns the concentration predicted
# at any particular depth given time. a, b, and c are the effective diffusion coefficient, initial chromium concentration,
# and the horizontal offset for the function, respectively.
def simulationFunction(x,a,b,c):
    Dcr = a
    C0 = b # 16.825 is the literature value
    Cs = c # Surface concentraion
    t = data_316.time # 3000 hr for this dataset
    x = np.array(x) # distances must be a numpy array to run through PEUQSE properly
    if Dcr < 0.0:
        return [float("nan")]
    else:
        y =  (Cs-C0)*(1-sp.erf(((x*10**(-6)))/(2*np.sqrt(Dcr*t*3600)))) + C0
        return y

x_values_for_data = []  # This is just initializing the global value to avoid confusion
x_values_for_data = data_316.distances

def simulation_function_wrapper(parametersArray): 
    global x_values_for_data
    a_given = parametersArray[0] # a "given" just means this wrapper will simulate using whatever a value it receives.
    b_given = parametersArray[1] # b "given" just means this wrapper will simulate using whatever b value it receives.
    c_given = parametersArray[2] # c "given" just means this wrapper will simulate using whatever c value it receives.
    
    # Call the simulation function with all the parameters passed in, y is returned as a numpy array
    y = simulationFunction(x_values_for_data, a_given, b_given, c_given) 
    # This is a check to determine if y contains any NaN values (NOTE: THIS CODE MAY SOON BE IMPLEMENTED IN PEUQSE AND BE REDUNDANT)
    y_array = np.array(y)
    nans_in_array = np.isnan(y_array)
    
    # If block checks for any NaN values in the result from the simulation function. If so, the wrapper returns a None object and the result is
    # treated as a 0-probability point
    if True in nans_in_array:
        # print("NaN detected in result array for simulation function.")
        sys.stdout.flush()
        return None
    else:
        # Only reached if the results contain no NaN values, hence a valid entry
        return y

# This wrapper takes a logarithm parameter a, and calls the simulation function with a converted to its anti-log, then returns the result
def simulation_function_using_log_a_wrapper(parametersWithAntiLogAArray):
    parametersArray = [10**parametersWithAntiLogAArray[0], parametersWithAntiLogAArray[1], parametersWithAntiLogAArray[2]]
    y = simulation_function_wrapper(parametersArray)
    return y                         
                                 

if __name__ == "__main__":
    print("function_316.py: Running independently.")
    print(simulation_function_using_log_a_wrapper([ np.log10(4.2000e-19),  1.6825e+01, 5]))