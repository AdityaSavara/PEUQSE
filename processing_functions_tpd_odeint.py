# These functions are for using the odeint-based temperature-programmed desorption model.  These functions set up the odeint call and process the data for obtaining the rates and the output values feed to the likelihood probability computation in 'run_me.py'.  These functions are imported in UserInput and the names are set there to be accessed here by 'UserInput. ' syntax.  EAW moved code here originally written by AS 2020/01/27.
import numpy as np
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
import scipy
from scipy.integrate import odeint
import pandas as pd

#Need to define beta directly, or define dt and dT.
dT = 0.77 #Set this to 0 for an isothermal experiment.
dt = 0.385
beta_dTdt = dt/dT #This beta is heating rate. This will be set to 0 if somebody sets TPR to false. Not to be confused with 1/(T*k_b) which is often also called beta. User can put beta in manually.
T_0 = 152.96 #this is the starting temperature.

def rate_tot_summing_func(rate):
    rate_tot = -np.sum(rate, axis=0)
    return rate_tot

def rate_tot_four_points_func(rate): #Multiple layers of wrapper functions are fine.
    rate_tot = rate_tot_summing_func(rate)
    #rate_tot_four_points = np.array(rate_tot[temp_points])
    return rate_tot

def log10_wrapper_func(rate, fourpoints = False):
    if fourpoints == True:
        temp_points = UserInput.temp_points #range(225) #FIXME: There should be nothing hard coded here. You can hard code it in userinput if you want.
        rate = rate[temp_points]
    loggedRateValues = np.log10(rate)
    return loggedRateValues

def observedResponses():
    global experiment
    #return np.log10(experiment[UserInput.temp_points])
    return np.log10(experiment)

def TPR_simulationFunctionWrapper(discreteParameterVector): 
    global times
    sample_list = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [UserInput.model_function_name, UserInput.initial_concentrations_array, times, (*sample_list,beta_dTdt,T_0) ] #FIXME: Times needs to occur in UserInput. This needs to all occur in somekind of UserFunctions module called from UserInput, should not be passed in here. 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. #FIXME: initialArgs and equation should come from UserInput, not be hardcoded here.            
    simulationInputArguments = [tpr_theta, times, *sample_list, beta_dTdt,T_0] #FIXME:To be passed in from userInput
    simulationFunction = UserInput.model_function_name #FIXME: To be passed in from userInput
    simulationOutput = UserInput.model_function_name(*simulationInputArguments) # EAW 2020/01/08
    return simulationOutput

def import_experimental_settings(Filename): #FIXME: This is obviously not very general. Though actually, we don't need it here. This should just go into UserInput as code rather than a function. These variables will become something like: UserInput.times, UserInput.observedResponse, UserInput.responseUncertainties.
    global times
    global experiment
    global errors
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000 
    print(len(experiment))
    #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    return times, experiment, errors
