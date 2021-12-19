# These functions are for using the odeint-based temperature-programmed desorption model.  These functions set up the odeint call and process the data for obtaining the rates and the output values feed to the likelihood probability computation in 'run_me.py'.  These functions are imported in UserInput and the names are set there to be accessed here by 'UserInput. ' syntax.  EAW moved code here originally written by AS 2020/01/27.
import numpy as np
import scipy
from scipy.integrate import odeint
import pandas as pd
#from tprmodel import tprequation, tprequationPiecewise

#Need to define beta directly, or define dt and dT.
dT = 0.77 #Set this to 0 for an isothermal experiment.
dt = 0.385
beta_dTdt = dt/dT #This beta is heating rate. This will be set to 0 if somebody sets TPR to false. Not to be confused with 1/(T*k_b) which is often also called beta. User can put beta in manually.
T_0 = 152.96 #this is the starting temperature.
initial_concentrations_array = [0.5, 0.5]
temp_points = np.array([0,49,99,149])
temp_points = (np.linspace(0,224, num=225)).astype(int)


def no_log_wrapper_func(individual_rates_vector):
    rate_tot = -np.sum(individual_rates_vector, axis=0)
    return rate_tot

def observedResponsesFunc():
    global experiment_rates
    values = experiment_rates[temp_points]
    #print("line 31, processing", values)
    return values


def TPR_internalPiecewiseSimulationFunctionWrapper(discreteParameterVector): 
    global times
    discreteParameterVectorList = list(discreteParameterVector) #converting to list so can use list expansion in arguments. 
    regularParams = discreteParameterVector[0:4] #This is Ea1, A1, gamma1
    Ea_modParams = discreteParameterVector[4:] #This is Ea_mod1, Ea_mod2, etc. from 3rd index to end.
    tpr_theta_Arguments = [tprequationPiecewise, initial_concentrations_array, times, (*regularParams,beta_dTdt,T_0,*Ea_modParams) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *regularParams, beta_dTdt,T_0, *Ea_modParams] 
    simulationOutput = tprequationPiecewise(*simulationInputArguments)
    return simulationOutput

def TPR_simulationFunctionWrapper(discreteParameterVector): 
    global times
    discreteParameterVectorList = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [tprequation, initial_concentrations_array, times, (*discreteParameterVectorList,beta_dTdt,T_0) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *discreteParameterVectorList, beta_dTdt,T_0] 
    simulationOutput = tprequation(*simulationInputArguments) # EAW 2020/01/08
    return simulationOutput

def TPR_simulationFunctionWrapperPiecewise(discreteParameterVector): 
    global times
    discreteParameterVectorList = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [tprequation, initial_concentrations_array, times, (*discreteParameterVectorList,beta_dTdt,T_0) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *discreteParameterVectorList, beta_dTdt,T_0] 
    simulationOutput = tprequation(*simulationInputArguments) # EAW 2020/01/08
    return simulationOutput 

def TPR_integerated_simulationFunctionWrapper(discreteParameterVector): 
    simulationOutput = TPR_simulationFunctionWrapper(discreteParameterVector)
    rate = no_log_wrapper_func(simulationOutput)
    global times
    from PEUQSE import littleEulerGivenArray
    times, integrated_desorption, rate = littleEulerGivenArray(0, times, rate)
    return integrated_desorption

def TPR_integerated_simulationFunctionWrapperPiecewise(discreteParameterVector): 
    simulationOutput = TPR_simulationFunctionWrapperPiecewise(discreteParameterVector)
    rate = no_log_wrapper_func(simulationOutput)
    global times
    from PEUQSE import littleEulerGivenArray
    times, integrated_desorption, rate = littleEulerGivenArray(0, times, rate)
    return integrated_desorption  
    
def import_experimental_settings(Filename): 
    global times
    global experiment_rates
    global errors
    
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment_rates = np.array(experiments_df['AcHBackgroundSubtracted'])/2000 
    #print(len(experiment))
    #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    global fourPoints
    fourPoints = False #<-- found this was true, somewhere along the way.
    if fourPoints == True:
        #print(temp_points)
        #print(len(temp_points))
        errors = errors[temp_points]
        times = times[temp_points]
    return times, experiment_rates, errors


def import_integrals_settings(Filename): 
    times, experiment_rates, experiment_rates_uncertainties = import_experimental_settings(Filename)
    from PEUQSE import littleEulerGivenArray, littleEulerUncertaintyPropagation
    times, integrated_desorption, experiment_rates = littleEulerGivenArray(0, times, experiment_rates)
    integrated_desorption_uncertainties = littleEulerUncertaintyPropagation(experiment_rates_uncertainties, times, 0.0)#The 0.0 is an initial coverage uncertainty.
    return times, integrated_desorption, integrated_desorption_uncertainties
        
#This is made for Example 8 which is using the Constant Errors and single site.
def import_experimental_settings_single(Filename): 
    global times
    global experiment_rates
    global errors
    
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment_rates = np.array(experiments_df['AcHBackgroundSubtracted'])/1000 
    #print(len(experiment))
    #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    global fourPoints
    return times, experiment_rates, errors

    
dT = 0.77 #Set this to 0 for an isothermal experiment.
dt = 0.385
beta_dTdt = dt/dT #This beta is heating rate. This will be set to 0 if somebody sets TPR to false. Not to be confused with 1/(T*k_b) which is often also called beta. User can put beta in manually.
T_0 = 152.96 #this is the starting temperature.
initial_concentrations_array = [0.5, 0.5]
temp_points = np.array([0,49,99,149])
temp_points = (np.linspace(0,224, num=225)).astype(int)
