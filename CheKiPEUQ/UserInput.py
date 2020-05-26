import numpy as np

#####Temperature Programmed Reaction Settings#####

#Response Plot Settings
simulated_response_plot_settings = {}
simulated_response_plot_settings['x_label'] = ''
simulated_response_plot_settings['y_label'] = ''
#simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
simulated_response_plot_settings['figure_name'] = 'Posterior_Simulated' #This is the default name for simulated response plots.
simulated_response_plot_settings['legend'] = True #Can be changed to false to turn off the legend.
#simulated_response_plot_settings['legendLabels'] = ['experiment', 'mu_guess', 'MAP'] here is an example of how to change the legend labels.
#simulated_response_plot_settings['fontdict']= {'size':20} A font dictionary can be passed in, this will be used for the axes and axes labels.

####BELOW ARE MODEL PARAMETERS, WE WILL WANT TO COMBINE THESE INTO A LIST OF PARAMETERS###
model = {} 
model['InputParameterPriorValues'] =  [] #Should be like: [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
model['InputParametersPriorValuesUncertainties'] = []# Should be like: [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
model['parameterNamesAndMathTypeExpressionsDict'] = {} #This must be provided. It can be as simple as {"Param1":"1"} etc. but it must be a dictionary with strings as keys and as values. The next line is a comment with a more complicated example.
#Example: model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
model['simulateByInputParametersOnlyFunction'] = None #A function must be provided! This cannot be left as None.
model['simulationOutputProcessingFunction'] = None
model['reducedParameterSpace'] = [] #This is to keep parameters as 'constants'. Any parameter index in this list will be allowed to change, the rest will be held as constants.
model['responses_simulation_uncertainties'] = None #Can be none, a list/vector, or can be a function that returns the uncertainties after each simulation is done. The easiest way would be to have a function that extracts a list that gets updated in another namespace after each simulation.
model['custom_logLikelihood'] = None #This should point to a function that takes the discrete parameter values as an argument and returns "logLikelihood, simulatedResponses". So the function returns a value for the logLikelihood (or proportional to it). The function must *also* return the simulated response output, though technically can just return the number 0 as the ssecond return.  The function can be a simple as accessing a global dictionary. This feature is intended for cases where the likelihood cannot be described by a normal/gaussian distribution.
model['custom_logPrior'] = None  #This feature has not been implemented yet, but is intended for cases where the prior distribution is not described by a normal distribution. The user will provide a function that takes in the parameters and returns a logPrior (or something proportional to a logPrior).
model['InputParameterPriorValues_upperBounds'] = [] #This should be a list/array of the same shape as InputParameterPriorValues. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then truncates any distribution to have a probability of ~0 when any of the parameters go outside of their bounds. ##As of May 4th 2020, this only has been checked for scaling_uncertainties_type = 'off'
model['InputParameterPriorValues_lowerBounds'] = []#This should be a list/array of the same shape as InputParameterPriorValues. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then truncates any distribution to have a probability of ~0 when any of the parameters go outside of their bounds. ##As of May 4th 2020, this only has been checked for scaling_uncertainties_type = 'off'

#####Experimental Data Input Files#####
responses = {}
responses['responses_abscissa'] = []
responses['responses_observed'] = []
responses['responses_observed_uncertainties'] = []
responses['reducedResponseSpace'] = []
responses['independent_variables_values'] = []
responses['independent_variables_names'] = []

#####Parameter Estimation Inputs#####
parameter_estimation_settings = {}
parameter_estimation_settings['verbose'] = False
parameter_estimation_settings['exportLog'] = True
parameter_estimation_settings['exportAllSimulatedOutputs'] = False
parameter_estimation_settings['checkPointFrequency'] = None
parameter_estimation_settings['scaling_uncertainties_type'] = "std" #"std" is for standard deviation. there is also "off" and the option of "mu" for using the absolute values of the mean(s) of the prior distribution(s). If a scalar is entered (a float) then that fixed value will be used for all scalings.
parameter_estimation_settings['undo_scaling_uncertainties_type'] = False #This undoing can be set to True but presently only works for the case of fixed scaling (a single scalar).
				  
######MCMC settings:#####
parameter_estimation_settings['mcmc_random_seed'] = None #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
parameter_estimation_settings['mcmc_mode'] = 'unbiased' #can be 'unbiased', 'MAP_finding', or 'HPD_exploring', the exploring one should take the MAP as an initial guess.
parameter_estimation_settings['mcmc_length'] = 10000
parameter_estimation_settings['mcmc_burn_in'] = 500 
parameter_estimation_settings['mcmc_relative_step_length'] = 0.1 #Default value is of 0.1, but values such as 1 are also quite reasonable. This is the step length relative to the covmat of the prior. So it is relative to the variance, not relative to the standard deviation.
parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
parameter_estimation_settings['mcmc_info_gain_cutoff'] = 0  #A typical value is 1E-5. Use 0 to turn this setting off. Allowing values that are too small will cause numerical errors, this serves as a highpass filter.
parameter_estimation_settings['mcmc_info_gain_returned'] = 'KL_divergence' # #current options are 'log_ratio' and 'KL_divergence' where 'KL' stands for Kullback-Leibler


#####Plot Settings#####
#possible dictionary fields include: dpi, figure_name, fontsize, x_label, y_label, figure_name, x_range, y_range
samplingScatterMatrixPlotsSettings ={}

######gridSamplingSettings##### 
#At present, all gridSampling settings are fed as arguments directly into the doGridSearch function.
#Perhaps that should be changed in the future so that wrapper functions can pass arguments to doGridSearch.
#parameter_estimation_settings['gridSampling'] = False    
#doGridSearch(self, searchType='doMetropolisHastings', export = True, verbose = False, gridSamplingIntervalSize = [], gridSamplingRadii = [], passThroughArgs = {}):

######mumpce plots#####
#model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$E_{a1}$'},
#{'parameter_number': 1, 'parameter_name': r'$E_{a2}$'},
#{'parameter_number': 2, 'parameter_name': r'$log(A_{1})$'},
#{'parameter_number': 3, 'parameter_name': r'$log(A_{2})$'},
#{'parameter_number': 4, 'parameter_name': r'$\gamma_{1}$'},
#{'parameter_number': 5, 'parameter_name': r'$\gamma_{2}$'}])
active_parameters = [] #Blank by default: gets populated with all parameters (or reduced parameters) if left blank. Warning: trying to set this manually while using the reduced parameters feature is not supported as of April 2020.
#pairs_of_parameter_indices = [[0, 1], [1, 2],[2, 3],[3, 4],[4, 5]] #This sets which parameters to plot contours for. By default, all pairs are plotted.
contour_settings_custom = {'figure_name': 'PosteriorContourPlots','fontsize':'auto' ,'num_y_ticks': 'auto','num_x_ticks':'auto','contours_normalized':True,'center_on':'all','colorbars':True} #'colormap_posterior_customized':'Oranges','colormap_prior_customized':'Greens'
#num_y_ticks and num_x_ticks must be either a string ('auto') or an integer (such as 4, either without string or with integer casting like int('5')).
parameter_pairs_for_contour_plots = [] #This will accept either strings (for variable names) or integers for positions.

####Design Of Experiments####
doe_settings = {} #To use this automated design of experiments the independent variables feature **must** be used.
doe_settings['info_gains_matrices_array_format'] = 'xyz' #options are 'xyz' and 'meshgrid'.  Images are only ouput when scanning two independent variables. If using more than two, it is probably better to use the 'xyz' format and inspect the final info_gains_matrices_array directly. Note that this setting must be set *before* running the doe command. You cannot change the format of the info_gains_matrices_array afterwards because the way the sampling is conducted will change based on this setting.
doe_settings['independent_variable_grid_center'] = [] #This must be a 1D array/list with length of number of independent variables.  
doe_settings['independent_variable_grid_interval_size'] = [] #This must be a 1D array/list with length of number of independent variables.  
doe_settings['independent_variable_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of independent variables.

doe_settings['on_the_fly_conditions_grids'] = True #This makes the independent variable grid each time. This costs more processing but less memory. As of April 2020 the other option has not been implemented but would just require making the combinations into a list the first time and then operating on a copy of that list.

#doe_settings['parameter_modulation_grid_center'] #We do NOT create such a variable. The initial guess variable is used, which is the center of the prior if not filled by the user.
doe_settings['parameter_modulation_grid_interval_size'] = [] #This must be 1D array/list with length of number of parameters.  These are all relative to the standard deviation of the prior of that parmaeter. 
doe_settings['parameter_modulation_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of paramaeters.