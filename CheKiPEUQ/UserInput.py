import numpy as np


#Response Plot Settings
simulated_response_plot_settings = {}
simulated_response_plot_settings['x_label'] = ''
simulated_response_plot_settings['y_label'] = ''
#simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
simulated_response_plot_settings['figure_name'] = 'Posterior_Simulated' #This is the default name for simulated response plots.
simulated_response_plot_settings['legend'] = True #Can be changed to false to turn off the legend.
#simulated_response_plot_settings['legendLabels'] = ['experiment', 'mu_guess', 'MAP'] here is an example of how to change the legend labels.
simulated_response_plot_settings['error_linewidth'] = 'auto' #Integer. Using "auto" or "None" sets to "20" when there is only 1 point, 1 when number of points is > 10, and "4" when number of points is between 1 and 10 and. Using '0' or 'none' will hide the error bars.
simulated_response_plot_settings['fontdict']= {'size':16} #A font dictionary can be passed in, this will be used for the axes and axes labels.

#### Model Paramerters Variables ###
model = {} 
model['InputParameterPriorValues'] =  [] #Should be like: [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
model['InputParametersPriorValuesUncertainties'] = []# Should be like: [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #A value of -1 in the Uncertainties indicates the parameter in question is described a univorm distribution In this case, the InputParameterPriorValues_upperBounds and InputParameterPriorValues_lowerBounds must be defined for each parmeter (can be defined as None for the non-uniform parameters).
model['InputParameterInitialGuess'] = [] #This is optional. An initial guess changes where the search is centered without changing the prior. If no initial guess is proided, the InputParameterPriorValues are taken as an initial guess.
model['parameterNamesAndMathTypeExpressionsDict'] = {} #This must be provided. It can be as simple as {"Param1":"1"} etc. but it must be a dictionary with strings as keys and as values. The next line is a comment with a more complicated example.
#Example: model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
model['populateIndependentVariablesFunction'] = None  #Not normally used. Mainly for design of experiments.
model['simulateByInputParametersOnlyFunction'] = None #A function must be provided! This cannot be left as None.
model['simulationOutputProcessingFunction'] = None
model['reducedParameterSpace'] = [] #This is to keep parameters as 'constants'. Any parameter index in this list will be allowed to change, the rest will be held as constants. For example, using [0,3,4] would allow the first, fourth, and fifth parameters to vary and would keep the rest as constants (note that array indexing is used).
model['responses_simulation_uncertainties'] = None #Can be none, a list/vector, or can be a function that returns the uncertainties after each simulation is done. The easiest way would be to have a function that extracts a list that gets updated in another namespace after each simulation.
model['custom_logLikelihood'] = None #This should point to a function that takes the discrete parameter values as an argument and returns "logLikelihood, simulatedResponses". So the function returns a value for the logLikelihood (or proportional to it). The function must *also* return the simulated response output, though technically can just return the number 0 as the ssecond return.  The function can be a simple as accessing a global dictionary. This feature is intended for cases where the likelihood cannot be described by a normal/gaussian distribution.
model['custom_logPrior'] = None  #This feature has  been implemented but not tested, it is intended for cases where the prior distribution is not described by a normal distribution. The user will provide a function that takes in the parameters and returns a logPrior (or something proportional to a logPrior). If MCMC will be performed, the user will still need to fill out InputParametersPriorValuesUncertainties with std deviations or a covariance matrix since that is used to decide the mcmc steps.
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
parameter_estimation_settings['checkPointFrequency'] = None #Deprecated. It will override all other checkpoint choices if it is changed from None. The user should use the similar variables below.
parameter_estimation_settings['scaling_uncertainties_type'] = "std" #"std" is for standard deviation. there is also "off" and the option of "mu" for using the absolute values of the mean(s) of the prior distribution(s). If a scalar is entered (a float) then that fixed value will be used for all scalings.
parameter_estimation_settings['undo_scaling_uncertainties_type'] = False #This undoing can be set to True but presently only works for the case of fixed scaling (a single scalar).
				  
######MCMC settings:#####
parameter_estimation_settings['mcmc_random_seed'] = None #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
parameter_estimation_settings['mcmc_mode'] = 'unbiased' #can be 'unbiased', 'MAP_finding', or 'HPD_exploring', the exploring one should take the MAP as an initial guess.
parameter_estimation_settings['mcmc_length'] = 10000   #This is the number of mcmc steps to take.
parameter_estimation_settings['mcmc_burn_in'] = 'auto' #This must be an integer or Auto. When it is set to auto it will be 10% of the mcmc_length (as of Oct 2020). 
parameter_estimation_settings['mcmc_relative_step_length'] = 0.1 #Default value is of 0.1, but values such as 1 are also quite reasonable. This is the step length relative to the covmat of the prior. So it is relative to the variance, not relative to the standard deviation.  As of Oct 2020, this only accepts the MetropolisHastings step size and not the EnsembleSliceSampling step size.
parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior during MetropolisHastings sampling. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
parameter_estimation_settings['mcmc_info_gain_cutoff'] = 0  #A typical value is 1E-5. Use 0 to turn this setting off. The purpose of this is that allowing values that are too small will cause numerical errors, this serves as a highpass filter.
parameter_estimation_settings['mcmc_info_gain_returned'] = 'KL_divergence' # #current options are 'log_ratio' and 'KL_divergence' where 'KL' stands for Kullback-Leibler
parameter_estimation_settings['mcmc_threshold_filter_samples'] = True #This feature removes low probability tails from the posterior. This can be important for getting mu_AP, especially when using ESS. Default is true.
parameter_estimation_settings['mcmc_threshold_filter_coefficient'] = 'auto' #This can be a float or the string 'auto'. Currently (Oct 2020), 'auto' sets the value is 2.0.  The smaller the value the more aggressive the filtering.
#The below settings are for ESS and/or parallel sampling#
parameter_estimation_settings['mcmc_nwalkers'] = 'auto'  #The number of walkers to use.  By default, if doing ESS, this is 4*numParameters. As of Oct 2020, this has no effect for MetropolisHastings.
parameter_estimation_settings['mcmc_maxiter'] = 1E6 #This is related to the expansions and contractions in ESS. It has a role similar to limiting the number of iterations in conventional regression. The ESS backend has a default value of 1E4, but in initial testing that was violated too often so 1E6 has been used now.
parameter_estimation_settings['mcmc_maxiter'] = 1E6 
parameter_estimation_settings['mcmc_walkerInitialDistribution'] = 'auto' #Can be 'uniform', 'gaussian', or 'identical'.  Auto will use 'uniform' during gridsearch and 'uniform' for most other cases.
parameter_estimation_settings['mcmc_checkPointFrequency'] = None #This is only for MH, not ESS. (as of Oct 2020)
parameter_estimation_settings['mcmc_parallel_sampling'] = False #This makes completely parallelized sampling of a single sampling. syntax to use is like "mpiexec -n 5 python runfile.py" where 5 is the number of processors. Currently, the first processor's results are thrown away.  In the future, this may change.


#####Plot Settings#####
#possible dictionary fields include: dpi, figure_name, fontsize, x_label, y_label, figure_name, x_range, y_range
samplingScatterMatrixPlotsSettings ={}

######multiStart (including gridsearch)##### 
parameter_estimation_settings['multistart_checkPointFrequency'] = None #Note: this setting does not work perfectly with ESS.
parameter_estimation_settings['multistart_parallel_sampling'] = False
parameter_estimation_settings['multistart_centerPoint'] = None #With None the centerPoint will be taken as model['InputParameterInitialGuess'] 
parameter_estimation_settings['multistart_numStartPoints'] = 0 #If this is left as zero it will be set as 3 times the number of active parameters.
parameter_estimation_settings['multistart_initialDistributionType'] = 'uniform' #Can be 'uniform', 'gaussian', 'identical', or 'grid'.
parameter_estimation_settings['multistart_relativeInitialDistributionSpread'] = 1.0 #This settting is for non-grid multistarts. The default value is 1.0. This scales the distribution's spread. By default, the uniform distribution, the points are sampled from a 2 sigma interval in each direction from the initial guess. This value then scales that range.
parameter_estimation_settings['multistart_gridsearchSamplingInterval'] = [] #This is for gridsearches and is in units of absolute intervals. By default, these intervals will be set to 1 standard deviaion each.  To changefrom the default, make a comma separated list equal to the number of parameters.
parameter_estimation_settings['multistart_gridsearchSamplingRadii'] = [] #This is for gridsearches and refers to the number of points (or intervals) in each direction to check from the center. For example, 3 would check 3 points in each direction plus the centerpointn for a total of 7 points along that dimension. For a 3 parameter problem, [3,7,2] would check radii of 3, 7, and 2 for those parameters.
parameter_estimation_settings['multistart_passThroughArgs'] = {}

parameter_estimation_settings['multistart_calculatePostBurnInStatistics'] = True
parameter_estimation_settings['multistart_keep_cumulative_post_burn_in_data'] = False
parameter_estimation_settings['multistart_exportLog'] = []
parameter_estimation_settings['multistart_passThroughArgs'] = {}


######mumpce plots##### #####contour plots#### 
#model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$E_{a1}$'},
#{'parameter_number': 1, 'parameter_name': r'$E_{a2}$'},
#{'parameter_number': 2, 'parameter_name': r'$log(A_{1})$'},
#{'parameter_number': 3, 'parameter_name': r'$log(A_{2})$'},
#{'parameter_number': 4, 'parameter_name': r'$\gamma_{1}$'},
#{'parameter_number': 5, 'parameter_name': r'$\gamma_{2}$'}])

contour_plot_settings = {} #TODO: change below to contour_plot_settings. Make "parameter_pairs_for_contour_plots" and "active_parameters" to be within that.
active_parameters = [] #Blank by default: gets populated with all parameters (or reduced parameters) if left blank. Warning: trying to set this manually while using the reduced parameters feature is not supported as of April 2020.
#pairs_of_parameter_indices = [[0, 1], [1, 2],[2, 3],[3, 4],[4, 5]] #This sets which parameters to plot contours for. By default, all pairs are plotted.
contour_settings_custom = {'figure_name': 'PosteriorContourPlots','fontsize':16 ,'num_y_ticks': 'auto','num_x_ticks':'auto','contours_normalized':True,'center_on':'all','colorbars':True}  #Fontsize can also be set to "auto"
#See the file mumpce_custom_plotting_example.py for the full set of arguments that can be provided inside contour_settings_custom.
#'colormap_posterior_customized':'Oranges','colormap_prior_customized':'Greens'
#num_y_ticks and num_x_ticks must be either a string ('auto') or an integer (such as 4, either without string or with integer casting like int('5')).
parameter_pairs_for_contour_plots = [] #This will accept either strings (for variable names) or integers for positions.

####Design Of Experiments####
#The design of experiments feature is used with syntax like PE_object.designOfExperiments()   
#This feature modulates the parameters to see how much info gain there would be in different parts of condition space (using synthetic data).
#The normal usage as of May 26 2020 is to first use the indpendent variables feature, which must be used, then fill the doe_settings below.
#Then call PE_object.doeParameterModulationCombinationsScanner()
#If a single parameter modulation grid is going to be used, one can instead define the independent variable grid and call like this: PE_object.createInfoGainPlots(plot_suffix="manual")
#For a real usage, see Example14doeFunctionExample
#Note: after calling either of these functions, the variables populated are PE_object.info_gains_matrices_array and PE_object.info_gain_matrix.  So if somebody wants to export these after using the functions, one can cycle across each info gains matrix inside PE_object.info_gains_matrices_array and export to csv.
#A key is printed out inside of Info_gain_parModulationGridCombinations.csv

doe_settings = {} #To use the design of experiments feature the independent variables feature **must** be used.
doe_settings['info_gains_matrices_array_format'] = 'xyz' #options are 'xyz' and 'meshgrid'.  Images are only ouput when scanning two independent variables. If using more than two independent variables, you will need to use the 'xyz' format and will need to analyze the final info_gains_matrices_array written to file directly. Note this variable must be set before running the doe command. You cannot change the format of the info_gains_matrices_array afterwards because the way the sampling is stored during a run is change based on this setting.

doe_settings['info_gains_matrices_multiple_parameters'] = 'sum' #The possible values are 'sum' or 'each'.  The 'each' choice exports info_gains for **each** parameter (and also exports the sums). #This feature is (for now) only for KL_divergence.

#Now we will define how big of a modulation / offset we want to apply to each parameter.
#doe_settings['parameter_modulation_grid_center'] #We do NOT create such a variable. The initial guess variable is used, which is the center of the prior if not filled by the user.
doe_settings['parameter_modulation_grid_interval_size'] = [] #This must be 1D array/list with length of number of parameters.  These are all relative to the standard deviation of the prior of that parmaeter. Such as [1,1]  #If you wish to manually change this setting, you should use all non-zero values for this setting, even if you plan to not modulate that parmeter (which parameters will be modulated are in the num_intervals variable below).
doe_settings['parameter_modulation_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of parameters. #Such as [1,1].  #This is the number of steps in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the parameter to vary.

#Now we will define the *conditions space* to explore (to find the highest info gain), which will be done for *each* modulation.
#Note that this means that responses['independent_variables_values'] must be used, AND it must be fed into the model simulation file as a connected variables array. (See  Example13doeFunctionExample directory runfiles).
doe_settings['independent_variable_grid_center'] = [] #This must be a 1D array/list with length of number of independent variables.  It's a central condition a grid will be made around. Example: [500, 0.5]
doe_settings['independent_variable_grid_interval_size'] = [] #This must be a 1D array/list with length of number of independent variables.  this how big of each step will be in each direction/dimension (it is the grid spacing).  You should always use non-zero values for this setting, even if you plan to not vary that independent variable. Like [100,0.2]
doe_settings['independent_variable_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of independent variables. #This is the number of steps in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the condition to vary. Example: [2,2] 

doe_settings['on_the_fly_conditions_grids'] = True #Oct 2020: do not change this setting. #Values are True or False. This makes the independent variable grid each time. This costs more processing but less memory. As of April 2020 the other option has not been implemented but would just require making the combinations into a list the first time and then operating on a copy of that list.