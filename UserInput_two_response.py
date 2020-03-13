import numpy as np
import processing_function_two_response as fun

simulated_response_plot_settings = {}

####BELOW ARE MODEL PARAMETERS, WE WILL WANT TO COMBINE THESE INTO A LIST OF PARAMETERS###
model = {} 
model['InputParameterInitialValues'] = [1, 5] # theta_1, theta_2
model['InputParametersInitialValuesUncertainties'] = [1, 1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
model['simulationOutputProcessingFunction'] = fun.direct_parameters_to_observations
prep_cov_mat = np.zeros((2,2))
np.fill_diagonal(prep_cov_mat,model['InputParametersInitialValuesUncertainties'])
model['PriorCovarianceMatrix'] = prep_cov_mat
model['parameterNamesAndMathTypeExpressionsDict'] = {'theta_1':r'$/theta_{1}$','theta_2':r'$theta_{2}$'}

InputConstants= [] #TODO: ERIC, WE SHOULD EITHER DESIGN YOUR CODE TO ALLOW CONSTANTS SEPARATELY, OR TO HAVE UNCERTAINTIES OF ZERO TO MAKE THINGS INTO A CONSTANT. THAT IS UP TO YOU AT THIS STAGE.

#####Experimental Data Input Files#####
responses = {}
responses['responses_abscissa'] = np.array([0, 1]) # These represent theta_1 and theta_2
responses['responses_observed'] = np.array([[2], [3]]) # [2, 3], a 1-D array, may be also tested.  It should yield a similar result, apart from the randomness involved in every new simulation instance.
responses['responses_observed_uncertainties'] = np.array([1, 1])
#####Parameter Estimation Inputs#####
parameter_estimation_settings = {}
parameter_estimation_settings['verbose'] = False
parameter_estimation_settings['exportAllSimulatedOutputs'] = False
parameter_estimation_settings['checkPointFrequency'] = None
parameter_estimation_settings['exportLog'] = False
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
model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$\theta_{1}$'},
{'parameter_number': 1, 'parameter_name': r'$\theta_{2}$'}])
active_parameters = np.array([0, 1])
pairs_of_parameter_indices = [[0, 1], [1, 0]]
contour_settings_custom = {'figure_name': 'Example_1_plots_mumpce','fontsize':'auto' ,'num_y_ticks': 'auto','num_x_ticks':'auto','contours_normalized':False,'center_on':'all','colorbars':True} #'colormap_posterior_customized':'Oranges','colormap_prior_customized':'Greens'
parameter_pairs_for_contour_plots = [] #This will accept either strings (for variable names) or integers for positions.
