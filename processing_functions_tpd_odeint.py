# These functions are for using the odeint-based temperature-programmed desorption model.  These functions set up the odeint call and process the data for obtaining the rates and the output values feed to the likelihood probability computation in 'run_me.py'.  These functions are imported in UserInput and the names are set there to be accessed here by 'UserInput. ' syntax.  EAW moved code here originally written by AS 2020/01/27.
import numpy as np
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
import scipy
from scipy.integrate import odeint
import pandas as pd

def rate_tot_summing_func(rate):
    rate_tot = -np.sum(rate, axis=0)
    return rate_tot

def rate_tot_four_points_func(rate): #Multiple layers of wrapper functions are fine.
    rate_tot = rate_tot_summing_func(rate)
    temp_points = UserInput.temp_points #range(225) #FIXME: There should be nothing hard coded here. You can hard code it in userinput if you want.
    rate_tot_four_points = np.array(rate_tot[temp_points])
    return rate_tot_four_points

def log10_wrapper_func(rate):
    rate_tot_four_points = rate_tot_four_points_func(rate)
    loggedRateValues = np.log10(rate_tot_four_points)
    return loggedRateValues

def observedResponses(experiment):
    return np.log10(experiment[UserInput.temp_points])

def simulationFunctionWrapper(self, discreteParameterVector): #FIXME: This should be defined in UserInput and passed in. User is responsible for it.
    # from tprmodel import tprequation # This is moved to the beginning of this file EAW 2020/01/08
    sample_list = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
    tpr_theta_Arguments = [UserInput.model_function_name, UserInput.initial_concentrations_array, self.times, (*sample_list,self.UserInput.beta_dTdt,self.UserInput.T_0) ] #FIXME: Times needs to occur in UserInput. This needs to all occur in somekind of UserFunctions module called from UserInput, should not be passed in here. 
    tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. #FIXME: initialArgs and equation should come from UserInput, not be hardcoded here.            
    simulationInputArguments = [tpr_theta, self.times, *sample_list, self.UserInput.beta_dTdt,self.UserInput.T_0] #FIXME:To be passed in from userInput
    simulationFunction = UserInput.model_function_name #FIXME: To be passed in from userInput
    simulationOutput = UserInput.model_function_name(*simulationInputArguments) # EAW 2020/01/08
    return simulationOutput

def import_experimental_settings(self): #FIXME: This is obviously not very general. Though actually, we don't need it here. This should just go into UserInput as code rather than a function. These variables will become something like: UserInput.times, UserInput.observedResponse, UserInput.responseUncertainties.
    experiments_df = pd.read_csv(UserInput.Filename)
    self.times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    self.experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    self.errors = np.array(experiments_df['Errors']) #.to_numpy()
    #return self.times, self.experiment, self.errors
