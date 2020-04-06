import numpy as np

#####Temperature Programmed Reaction Settings#####
TPR = True #Set to false if doing an isothermal experiment.


#Response Plot Settings
simulated_response_plot_settings = {}
simulated_response_plot_settings['x_label'] = ''
simulated_response_plot_settings['y_label'] = ''
#simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
simulated_response_plot_settings['figure_name'] = 'Posterior_Simulated'

####BELOW ARE MODEL PARAMETERS, WE WILL WANT TO COMBINE THESE INTO A LIST OF PARAMETERS###
model = {} 
model['InputParameterInitialValues'] = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
#InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
model['InputParametersInitialValuesUncertainties'] = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
prep_cov_mat = np.zeros((6,6))
np.fill_diagonal(prep_cov_mat,model['InputParametersInitialValuesUncertainties'])
model['PriorCovarianceMatrix'] = prep_cov_mat
model['parameterNamesAndMathTypeExpressionsDict'] = {} #This must be provided. It can be as simple as {"Param1":"1"} etc. but it must be a dictionary with strings as keys and as values. The next line is a comment with a more complicated example.
#model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
#InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
model['simulateByInputParametersOnlyFunction'] = None #A function must be provided! This cannot be left as None.
model['simulationOutputProcessingFunction'] = None
model['reducedParameterSpace'] = [] #This is to keep parameters as 'constants'. Any parameter index in this list will be allowed to change, the rest will be held as constants.

#####Experimental Data Input Files#####
responses = {}
responses['responses_abscissa'] = []
responses['responses_observed'] = []
responses['responses_observed_uncertainties'] = []
responses['reducedResponseSpace'] = []

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
#model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$E_{a1}$'},
#{'parameter_number': 1, 'parameter_name': r'$E_{a2}$'},
#{'parameter_number': 2, 'parameter_name': r'$log(A_{1})$'},
#{'parameter_number': 3, 'parameter_name': r'$log(A_{2})$'},
#{'parameter_number': 4, 'parameter_name': r'$\gamma_{1}$'},
#{'parameter_number': 5, 'parameter_name': r'$\gamma_{2}$'}])
active_parameters = []
#pairs_of_parameter_indices = [[0, 1], [1, 2],[2, 3],[3, 4],[4, 5]]
contour_settings_custom = {'figure_name': 'PosteriorContourPlots','fontsize':'auto' ,'num_y_ticks': 'auto','num_x_ticks':'auto','contours_normalized':False,'center_on':'all','colorbars':True} #'colormap_posterior_customized':'Oranges','colormap_prior_customized':'Greens'
parameter_pairs_for_contour_plots = [] #This will accept either strings (for variable names) or integers for positions.