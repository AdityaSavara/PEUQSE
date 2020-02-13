# These functions are for using the odeint-based temperature-programmed desorption model.  These functions set up the odeint call and process the data for obtaining the rates and the output values feed to the likelihood probability computation in 'run_me.py'.  These functions are imported in UserInput and the names are set there to be accessed here by 'UserInput. ' syntax.  EAW moved code here originally written by AS 2020/01/27.
import numpy as np
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
import scipy
from scipy.integrate import odeint
import pandas as pd
from tprmodel import tprequation

#Need to define beta directly, or define dt and dT.
dT = 0.77 #Set this to 0 for an isothermal experiment.
dt = 0.385
beta_dTdt = dt/dT #This beta is heating rate. This will be set to 0 if somebody sets TPR to false. Not to be confused with 1/(T*k_b) which is often also called beta. User can put beta in manually.
T_0 = 152.96 #this is the starting temperature.
initial_concentrations_array = [0.5, 0.5]
temp_points = np.array([0,49,99,149])
temp_points = (np.linspace(0,224, num=225)).astype(int)


def rate_tot_summing_func(rate):
    rate_tot = -np.sum(rate, axis=0)
    return rate_tot

# def log10_wrapper_func(rate):
    # rate_tot_four_points = rate_tot_summing_func(rate)
    # loggedRateValues = np.log10(rate_tot_four_points)
    # return loggedRateValues

def no_log_wrapper_func(rate):
    rate = rate_tot_summing_func(rate)
    return rate

def observedResponsesFunc():
    global experiment
    values = experiment[temp_points]
    print("line 31, processing", values)
    return values

# def observedResponsesProxyFunc():
    # global experiment
    # log10Values = np.log10(experiment[temp_points])
    # print("line 37, processing", log10Values)
    # return log10Values

def TPR_simulationFunctionWrapper(discreteParameterVector): 
    global times
    sample_list = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [tprequation, initial_concentrations_array, times, (*sample_list,beta_dTdt,T_0) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *sample_list, beta_dTdt,T_0] 
    simulationFunction = tprequation
    simulationOutput = tprequation(*simulationInputArguments) # EAW 2020/01/08
    return simulationOutput

def import_experimental_settings(Filename): 
    global times
    global experiment
    global errors
    
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000 
    print(len(experiment))
    #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    global fourPoints
    fourPoints = True
    if fourPoints == True:
        print(temp_points)
        print(len(temp_points))
        errors = errors[temp_points]
        times = times[temp_points]
    return times, experiment, errors
