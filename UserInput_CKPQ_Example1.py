import numpy as np

#####Temperature Programmed Reaction Settings#####
TPR = True #Set to false if doing an isothermal experiment.


#####Chemical Kinetic Model Input Info#####
from tprmodel import tprequation # EAW 2020/01/13 
model_function_name = tprequation # EAW 2020/01/08
from processing_functions_tpd_odeint import observedResponsesFunc,  TPR_simulationFunctionWrapper, import_experimental_settings, no_log_wrapper_func 
import processing_functions_tpd_odeint


#THIS is not part of a normal UserInput file. This is just part of our toy model.
Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedLargerErrors.csv'
times, experiment, observedResponses_uncertainties = import_experimental_settings(Filename)
import pandas as pd
experiments_df = pd.read_csv(Filename)


num_data_points = len(processing_functions_tpd_odeint.temp_points) #FIXME: This is temporary hardcoding to get it out of the run_me file.

simulated_response_plot_settings = {}
simulated_response_plot_settings['x_label'] = 'T (K)'
simulated_response_plot_settings['y_label'] = r'$rate (s^{-1})$'
#simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
simulated_response_plot_settings['figure_name'] = 'tprposterior'

######processing functions for odeint-based temperature-programmed desorption model##########
import_experimental_settings 

####BELOW ARE MODEL PARAMETERS, WE WILL WANT TO COMBINE THESE INTO A LIST OF PARAMETERS###
model = {} 
model['InputParameterInitialValues'] = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
#InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
model['InputParametersInitialValuesUncertainties'] = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
prep_cov_mat = np.zeros((6,6))
np.fill_diagonal(prep_cov_mat,model['InputParametersInitialValuesUncertainties'])
model['PriorCovarianceMatrix'] = prep_cov_mat
model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
#InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
model['simulateByInputParametersOnlyFunction'] = TPR_simulationFunctionWrapper
model['simulationOutputProcessingFunction'] = None

InputConstants= [] #TODO: ERIC, WE SHOULD EITHER DESIGN YOUR CODE TO ALLOW CONSTANTS SEPARATELY, OR TO HAVE UNCERTAINTIES OF ZERO TO MAKE THINGS INTO A CONSTANT. THAT IS UP TO YOU AT THIS STAGE.

#####Experimental Data Input Files#####
responses = {}
responses['responses_abscissa'] = np.array(experiments_df['AcH - T'])
responses['responses_observed'] = np.array(experiments_df['AcHBackgroundSubtracted'])/2000
responses['responses_observed_uncertainties'] = observedResponses_uncertainties
#responses_abscissa = np.array(experiments_df['AcH - T'])
#responses_observed = np.array(experiments_df['AcHBackgroundSubtracted'])/2000
observedResponses = observedResponsesFunc()


#####Parameter Estimation Inputs#####
parameter_estimation_settings = {}
parameter_estimation_settings['verbose'] = False
parameter_estimation_settings['exportLog'] = True
parameter_estimation_settings['exportAllSimulatedOutputs'] = False
parameter_estimation_settings['checkPointFrequency'] = None
parameter_estimation_settings['gridsearch'] = False
parameter_estimation_settings['verbose']
parameter_estimation_settings['scaling_uncertainties_type'] = "std"
				  
######MCMC settings:#####
parameter_estimation_settings['mcmc'] = True 
parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
parameter_estimation_settings['mcmc_mode'] = 'MAP_finding' #can be 'unbiased', 'MAP_finding', or 'HPD_exploring', the exploring one should take the MAP as an initial guess.
parameter_estimation_settings['mcmc_length'] = 10000
parameter_estimation_settings['mcmc_burn_in'] = 500 
parameter_estimation_settings['mcmc_relative_step_length'] = 0.1
parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.

######gridSamplingSettings#####
parameter_estimation_settings['gridSampling'] = False    


######mumpce plots#####
model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$E_{a1}$'},
{'parameter_number': 1, 'parameter_name': r'$E_{a2}$'},
{'parameter_number': 2, 'parameter_name': r'$log(A_{1})$'},
{'parameter_number': 3, 'parameter_name': r'$log(A_{2})$'},
{'parameter_number': 4, 'parameter_name': r'$\gamma_{1}$'},
{'parameter_number': 5, 'parameter_name': r'$\gamma_{2}$'}])
active_parameters = np.array([0, 1, 2, 3, 4, 5])
pairs_of_parameter_indices = [[0, 1], [1, 2],[2, 3],[3, 4],[4, 5]]
contour_settings_custom = {'figure_name': 'Example_1_plots_mumpce','fontsize':'auto' ,'num_y_ticks': 'auto','num_x_ticks':'auto','contours_normalized':False,'center_on':'all','colorbars':True} #'colormap_posterior_customized':'Oranges','colormap_prior_customized':'Greens'
parameter_pairs_for_contour_plots = [] #This will accept either strings (for variable names) or integers for positions.