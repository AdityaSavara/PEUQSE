# These functions are for using the odeint-based temperature-programmed desorption model.  These functions set up the odeint call and process the data for obtaining the rates and the output values feed to the likelihood probability computation in 'run_me.py'.  These functions are imported in UserInput and the names are set there to be accessed here by 'UserInput. ' syntax.  EAW moved code here originally written by AS 2020/01/27.
import numpy as np
import scipy
from scipy.integrate import odeint
import pandas as pd
from tprmodel import tprequation, tprequationPiecewise

#Since beta is not one of the parameters for estimation, we need to define beta directly, or define dt and dT.
#This can then be changed in the namespace later if desired. The same for the inital temperature and initial conncentrations.
dT = 0.77 #Set this to 0 for an isothermal experiment.
dt = 0.385
beta_dTdt = dt/dT #This beta is heating rate. This will be set to 0 if somebody sets TPR to false. Not to be confused with 1/(T*k_b) which is often also called beta. User can put beta in manually.
T_0 = 152.96 #this is the starting temperature.
initial_concentrations_array = [0.5, 0.5]

#This is intended to be a post processing function, since we usually (in TPR) use the negative of the rate. 
def neg_sum_of_all_rates(individual_rates_vector):
    rate_tot = -np.sum(individual_rates_vector, axis=0)
    return rate_tot

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


#this is what CheKiPEUQ considers the "Basic" simulation function. It takes the discrete parameter vector, and then simulates.
def TPR_simulationFunctionWrapperRatioFirst(discreteParameterVectorRatioFirst): 
    global times
    discreteParameterVectorList = list(discreteParameterVectorRatioFirst) #converting to list so can use pop and also list expansion in arguments.  
    initial_concentrations_array[0] = 1-discreteParameterVectorRatioFirst[0]
    initial_concentrations_array[1] = discreteParameterVectorRatioFirst[0]
    discreteParameterVectorList.pop(0)
    tpr_theta_Arguments = [tprequation, initial_concentrations_array, times, (*discreteParameterVectorList,beta_dTdt,T_0) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *discreteParameterVectorList, beta_dTdt,T_0] 
    simulationOutput = tprequation(*simulationInputArguments) 
    return simulationOutput

#Below is like above, but for doing a piecewise coverage dependent simulation.
def TPR_simulationFunctionWrapperPiecewise(discreteParameterVector): 
    global times
    discreteParameterVectorList = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [tprequation, initial_concentrations_array, times, (*discreteParameterVectorList,beta_dTdt,T_0) ] 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. 
    simulationInputArguments = [tpr_theta, times, *discreteParameterVectorList, beta_dTdt,T_0] 
    simulationOutput = tprequation(*simulationInputArguments)
    return simulationOutput 

#below converts the simulation into an integral, so it is just a wrapper around the "Base" case above.
def TPR_integerated_simulationFunctionWrapperRatioFirst(discreteParameterVectorRatioFirst): 
    simulationOutput = TPR_simulationFunctionWrapperRatioFirst(discreteParameterVectorRatioFirst)
    rate = neg_sum_of_all_rates(simulationOutput)
    global times
    from CheKiPEUQ import littleEulerGivenArray
    times, integrated_desorption, rate = littleEulerGivenArray(0, times, rate)
    return integrated_desorption

#Below is like above, but for doing a piecewise coverage dependent simulation.
def TPR_integerated_simulationFunctionWrapperPiecewise(discreteParameterVector): 
    simulationOutput = TPR_simulationFunctionWrapperPiecewise(discreteParameterVector)
    rate = neg_sum_of_all_rates(simulationOutput)
    global times
    from CheKiPEUQ import littleEulerGivenArray
    times, integrated_desorption, rate = littleEulerGivenArray(0, times, rate)
    return integrated_desorption  

#This is a support function that just returns the abscissa ("times") the observed values ("experiment_rates") and the uncertainties for those responses ("errors")    
def import_experimental_settings(Filename): 
    global times
    global experiment_rates
    global errors
    
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment_rates = np.array(experiments_df['AcHBackgroundSubtracted']) 
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    return times, experiment_rates, errors


#This is a support function like above only it returns the integrated amount of species desorbed. 
def import_integrals_settings(Filename): 
    times, experiment_rates, experiment_rates_uncertainties = import_experimental_settings(Filename)
    from CheKiPEUQ import littleEulerGivenArray, littleEulerUncertaintyPropagation
    times, integrated_desorption, experiment_rates = littleEulerGivenArray(0, times, experiment_rates)
    integrated_desorption_uncertainties = littleEulerUncertaintyPropagation(experiment_rates_uncertainties, times, 0)#The 0.2 is an initial coverage uncertainty.
    return times, integrated_desorption, integrated_desorption_uncertainties
        
#This is made for Example 8 which is using the Constant Errors and single site.
def import_experimental_settings_single(Filename): 
    global times
    global experiment_rates
    global errors
    
    experiments_df = pd.read_csv(Filename)
    times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    experiment_rates = np.array(experiments_df['AcHBackgroundSubtracted']) 
    #print(len(experiment))
    #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    errors = np.array(experiments_df['Errors']) #.to_numpy()
    return times, experiment_rates, errors

    