import numpy as np
###User sets their model equation####
from tprmodel import tprequation # EAW 2020/01/13 
model_function_name = tprequation # EAW 2020/01/08
from processing_functions_tpd_odeint import rate_tot_summing_func, observedResponsesFunc,  TPR_simulationFunctionWrapper, import_experimental_settings, no_log_wrapper_func #,observedResponsesProxyFunc, log10_wrapper_func 
#####Temperature Programmed Reaction Settings#####
TPR = True #Set to false if doing an isothermal experiment.

import processing_functions_tpd_odeint
num_data_points = len(processing_functions_tpd_odeint.temp_points) #FIXME: This is temporary hardcoding to get it out of the run_me file.

####BELOW ARE MODEL PARAMETERS, WE WILL WANT TO COMBINE THESE INTO A LIST OF PARAMETERS###
InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
InputConstants= [] #TODO: ERIC, WE SHOULD EITHER DESIGN YOUR CODE TO ALLOW CONSTANTS SEPARATELY, OR TO HAVE UNCERTAINTIES OF ZERO TO MAKE THINGS INTO A CONSTANT. THAT IS UP TO YOU AT THIS STAGE.

#####Experimental Data Input Files#####
Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedLargerErrors.csv'
times, experiment, observedResponses_uncertainties = import_experimental_settings(Filename)

#####Chemical Kinetic Model Input Files#####
parameterNamesAndMathTypeExpressionsDict = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}

#####Chemical Kinetic Model Initial Concentrations#####
#initial_concentrations_dict = {}
#initial_concentrations_array = [0.5, 0.5]


#####Bayesian Probability Parameters#####
verbose = False
#TODO: I think that we should replace below with mu_prior = np.array(InputParameterValues). That will make it easier for most users.
mu_prior = np.array(InputParameterInitialValues)
#mu_prior = np.array([41.5, 41.5, 13.0, 13.0, 0.1, 0.1]) 

if len(np.shape(InputParametersInitialValuesUncertainties)) == 1 and (len(InputParametersInitialValuesUncertainties) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
    cov_prior = np.diagflat(InputParametersInitialValuesUncertainties) 
elif len(np.shape(InputParametersInitialValuesUncertainties)) > 1: #If it's non-1D, we assume it's already a covariance matrix.
    cov_prior = InputParametersInitialValuesUncertainties
else: #If a blank list is received, that means the user
    print("The covariance of the priors is undefined because InputParametersInitialValuesUncertainties is blank.")
#    cov_prior = np.array([[200.0, 0., 0., 0., 0., 0.], 
#                          [0., 200.0, 0., 0., 0., 0.],
#                          [0., 0., 13.0, 0., 0., 0.],
#                          [0., 0., 0., 13.0, 0., 0.],
#                          [0., 0., 0., 0., 0.1, 0.],
#                          [0., 0., 0., 0., 0., 0.1]])
#

				  
######MCMC settings:#####
mcmc = True 
mcmc_random_seed = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
mcmc_length = 1000
mcmc_burn_in = 500 
modulate_accept_probability = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.

######gridSamplingSettings#####
gridSampling = False    


######processing functions for odeint-based temperature-programmed desorption model##########
simulationFunctionWrapper = TPR_simulationFunctionWrapper
import_experimental_settings = import_experimental_settings

simulationOutputProcessingFunction = no_log_wrapper_func
observedResponses = observedResponsesFunc()

######mumpce plots#####
model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': 'Parameter 0', 'parameter_value': 1.0},
 {'parameter_number': 1, 'parameter_name': 'Parameter 1'},
 {'parameter_number': 2, 'parameter_name': 'Parameter 2'},
 {'parameter_number': 3, 'parameter_name': 'Parameter 3'},
 {'parameter_number': 4, 'parameter_name': 'Parameter 4'},
 {'parameter_number': 5, 'parameter_name': 'Parameter 5'},
 {'parameter_number': 6, 'parameter_name': 'Parameter 6'}])
active_parameters = np.array([0, 1, 2, 4, 6])
posterior_mu_vector = np.array([-0.58888733,1.1200355, 0.00704044, -1.62385888,0.80439847])
posterior_cov_matrix = np.array([[ 0.0148872,-0.01894579, -0.01047339,0.01325883,0.04734254],
 [-0.01894579,0.04284732, -0.00131389, -0.04801795, -0.04545703],
 [-0.01047339, -0.00131389,0.02343653,0.01588293, -0.05618226],
 [ 0.01325883, -0.04801795,0.01588293,0.08171972,0.00875017],
 [ 0.04734254, -0.04545703, -0.05618226,0.00875017,0.20669273]])
prior_mu_vector = np.array([-0.98888733,0.8200355, 0.01204044, -7.02385888,0.40439847])
prior_cov_matrix = 10*posterior_cov_matrix
pairs_of_parameter_indices = [[0, 1], [1, 2],[3, 4]]
# Other options include: figure_name, fontsize, num_y_ticks, num_x_ticks, colormap_posterior_customized, colormap_prior_customized, contours_normalized, colorbars
