import numpy as np

#####Directory Settings####
directories = {} # lower case letters are recommended when creating directories to avoid possible errors
directories['graphs'] = "./graphs/"  #To include the graphs in the main runfile directory, just make this "./" 
directories['logs_and_csvs'] = "./logs_and_csvs/" #to include the logs and csvs in the main runfile directory, just make this "./" 
directories['pickles'] = "./pickles/" #to include the pickles and dills the main runfile directory, just make this "./" 

#####Experimental Data Input Files#####
responses = {}
responses['responses_abscissa'] = [] #Make 1 or more list or array within a list.
responses['responses_observed'] = [] #Make 1 list/array for each response.
responses['responses_observed_uncertainties'] = [] #Normally should not be blank, and should be provided with the same structure as responses_observed. One standard deviation of uncertainty should be provided for each response value. To set the responses_observed_uncertainties to zero, this variable or the values inside must really be set equal to 0. A blank list will not result in zeros and will autogenerate uncertainties relative to the responses_observed. A full covariance matrix can alternatively be used, but not all features are compatible with a full covariance matrix. If using a full covariance matrix, like in example 7j, please note that only the 'bottom left' of the covariance matrix is used, such that np.array([[2,1],[1,3]]) is actually the same as np.array([[2,-500],[1,3]]).  Here, the -500 would be ignored.
responses['responses_observed_max_covmat_size'] = 100 #The user should nor normally change this. if any response has more datapoints than this, that response will have variances evaluated separately (only the diagonal of the covmat) in a way which changes the computational cost to linear scaling.  For most regular computers around Jan 25 2020, the crossover happens after a few hundred points, so this variable has been set to have a default value of 100.                                                                                                                                                                                           
responses['responses_observed_weighting'] = [] #This feature is not recommended for normal use. If used, the input should be the same shape as responses_observed_uncertainties. This adds coefficients to responses_observed_uncertainties based on 1/(weighting)^0.5 to 'back propagate' any additional weighting terms (in analogy to variance weighted SSR).  If the responses_observed_uncertainties are appropriately defined, this should generally not be needed. This feature is only compatible when responses_observed_uncertainties consists of standard deviations rather than a covariance matrix.
responses['reducedResponseSpace'] = [] #If there are multiple responses, the user can ignore some of the responses during parameter estimation.
responses['independent_variables_values'] = []  #This field is mainly for design of experiments, but can be used as a type of connected variables for other cases also.
responses['independent_variables_names'] = [] #names associated with independent_variables_values.
responses['num_responses'] = 'auto' #'auto' is recommended, though an integer can be put in directly.

#(Optional) data transforms  This is for transforming the responses to improve the objective function.  Will be applied on simulated data also. 
#This feature is not compatible with simulatedResponses_upperBounds and simulatedResponses_lowerBounds as of Dec 2020. Contact the developer if this is needed
responses['data_overcategory'] = '' #Choices are currently 'transient_kinetics' and 'steady_state_kinetics'.  If this is used, then one also needs to provide response_types ( One for each response dimension). Additional features are welcome.
responses['response_types'] = [] #Response types can currently  be 'P' (product), 'I' (intermediate), 'R' (reactant), 'T' (temperature), 'O' (other)
responses['response_data_types'] = [] #Response data types can be 'c' (concentration), 'r' (rate), 'o' (other)


#### Model Paramerters Variables ###
model = {} 
model['InputParameterPriorValues'] =  [] #Should be like: [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
model['InputParametersPriorValuesUncertainties'] = []# Should be like: [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #A value of -1 in the Uncertainties indicates the parameter in question is described a univorm distribution In this case, the InputParameterPriorValues_upperBounds and InputParameterPriorValues_lowerBounds must be defined for each parmeter (can be defined as None for the non-uniform parameters). 
model['InputParameterInitialGuess'] = [] #This is optional. An initial guess changes where the search is centered without changing the prior. If no initial guess is proided, the InputParameterPriorValues are taken as an initial guess.
model['parameterNamesAndMathTypeExpressionsDict'] = {} #This should be a dictionary with strings as keys and as values, or it can just be a list. If it is not provided, the parmeters will get names like ParInd_0, ParInd_1, etc. The next line is a comment with a more complicated example.
#Example: model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
#A simpler list example would be ['Ea1', 'Ea2', 'logA1', 'logA2', 'gamma1', 'gamma2'], but then the graphs would have no subscripts.
model['populateIndependentVariablesFunction'] = None  #Not normally used. Mainly for design of experiments.
model['simulateByInputParametersOnlyFunction'] = None #A function must be provided! This cannot be left as None. The function should normally return an array the same size and shape as responses_observed, or should return a None object when the simulation fails or the result is considered non-physical. Alternatively, the function can written an object that needs to be processed further by SimulationOutputProcessingFunction.
model['simulationOutputProcessingFunction'] = None #An optional function may be provided which takes the outputs from simulateByInputParametersOnlyFunction and then processes them to match the size, shape, and scale of responses_observed. A None object can be returned when the simulation fails or the result is considered non-physical.
model['reducedParameterSpace'] = [] #This is to keep parameters as 'constants'. Any parameter index in this list will be allowed to change, the rest will be held as constants. For example, using [0,3,4] would allow the first, fourth, and fifth parameters to vary and would keep the rest as constants (note that array indexing is used).
model['responses_simulation_uncertainties'] = None #Optional. Can be none, a list/vector, or can be a function that returns the uncertainties after each simulation is done. The easiest way would be to have a function that extracts a list that gets updated in another namespace after each simulation.
model['custom_logLikelihood'] = None #Optional. This should point to a function that takes the discrete parameter values as an argument and returns "logLikelihood, simulatedResponses". So the function returns a value for the logLikelihood (or proportional to it). The function must *also* return the simulated response output, though technically can just return the number 0 as the ssecond return.  The function can be a simple as accessing a global dictionary. This feature is intended for cases where the likelihood cannot be described by a normal/gaussian distribution.
model['custom_logPrior'] = None  #Optional. This feature has  been implemented but not tested, it is intended for cases where the prior distribution is not described by a normal distribution. The user will provide a function that takes in the parameters and returns a logPrior (or something proportional to a logPrior). If MCMC will be performed, the user will still need to fill out InputParametersPriorValuesUncertainties with std deviations or a covariance matrix since that is used to decide the mcmc steps.
model['InputParameterPriorValues_upperBounds'] = [] #Optional. This should be a list/array of the same shape as InputParameterPriorValues. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then truncates any distribution to have a probability of ~0 when any of the parameters go outside of their bounds. ##As of May 4th 2020, this only has been checked for scaling_uncertainties_type = 'off'
model['InputParameterPriorValues_lowerBounds'] = []#Optional. This should be a list/array of the same shape as InputParameterPriorValues. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then truncates any distribution to have a probability of ~0 when any of the parameters go outside of their bounds. ##As of May 4th 2020, this only has been checked for scaling_uncertainties_type = 'off'
model['simulatedResponses_upperBounds'] = [] #Optional. Disallows responses outside of provided bounds. This should be a list/array of the same shape as responses_observed. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then sets the likelihood (and posterior) to ~0 when any of the responses go outside of their bounds.  Not compatible with data_overcategory feature.
model['simulatedResponses_lowerBounds'] = [] #Optional. Disallows responses outside of provided bounds. This should be a list/array of the same shape as responses_observed. Use a value of "None" for any parameter that should not be bounded in this direction.  The code then sets the likelihood (and posterior) to ~0 when any of the responses go outside of their bounds.  Not compatible with data_overcategory feature.
model['PosteriorParameterBounds'] = {} #Optional. Allows user to alter the parameter space after sampling has occurred. This should be a dictionary with the parameter name and a tuple with the lower bound first. 


#####Parameter Estimation Inputs#####
parameter_estimation_settings = {}
parameter_estimation_settings['verbose'] = False
parameter_estimation_settings['exportLog'] = True
parameter_estimation_settings['exportAllSimulatedOutputs'] = False #This feature (when set to true) behaves differently for multi-start and for mcmc. For mutli-start, all of the simulated responses for the final maps will be exported. For mcmc, all of the post-burn-in simulated outputs will be stored and exported.  Even if filtering is on, all of the simulated outputs will be exported, not just the filtered ones. This feature is presently not compatible with continueSampling. It will only export the simulatedOutputs from the most recent run. The feature has not been implemented for ESS.
parameter_estimation_settings['checkPointFrequency'] = None #Deprecated. It will override all other checkpoint choices if it is changed from None. The user should use the similar variables below.
parameter_estimation_settings['scaling_uncertainties_type'] = "std" #"std" is for standard deviation. there is also "off" and the option of "mu" for using the absolute values of the mean(s) of the prior distribution(s). If a scalar is entered (a float) then that fixed value will be used for all scalings.
parameter_estimation_settings['undo_scaling_uncertainties_type'] = False #This undoing can be set to True but presently only works for the case of fixed scaling (a single scalar).
parameter_estimation_settings['convergence_diagnostics'] = True #This feature when set to True will guide convergence analysis of the sampler. This is specific to the sampler.

######MCMC settings:#####
parameter_estimation_settings['mcmc_exportLog'] = True #exports additional information during the mcmc.
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
parameter_estimation_settings['mcmc_threshold_filter_benchmark'] = 'auto' # can be MAP, mu_AP, or auto to determine the filtering mechanism. The filtering with keyword MAP will filter points too far below the MAP. The keyword auto will filter samples that fall outside 2 std based on the distribution of the log of the log of the posterior distribution (the auto choice should be considered a convenient but not rigorous method of filtering, and essentially only removes probablities that are outliers).
##The below settings are for ESS and/or parallel sampling##
parameter_estimation_settings['mcmc_nwalkers'] = 'auto'  #The number of walkers to use.  By default, if doing ESS, this is 4*numParameters. As of Oct 2020, this has no effect for MetropolisHastings.
parameter_estimation_settings['mcmc_maxiter'] = 1E6 #This is related to the expansions and contractions in ESS. It has a role similar to limiting the number of iterations in conventional regression. The ESS backend has a default value of 1E4, but in initial testing that was violated too often so 1E6 has been used now.
parameter_estimation_settings['mcmc_maxiter'] = 1E6 
parameter_estimation_settings['mcmc_walkerInitialDistribution'] = 'auto' #Can be 'uniform', 'gaussian', or 'identical'.  Auto will use 'uniform' for most cases, including gridsearch.
parameter_estimation_settings['mcmc_walkerInitialDistributionSpread'] = 'auto' #Primarily for ensemble slice sampling. A number. The walkerInitialDistributionSpread is in relative units (relative to standard deviations). In the case of a uniform inital distribution the default level of spread is actually across two standard deviations, so the walkerInitialDistributionSpread is relative to that (that is, a value of 2 would give 2*2 = 4 for the full spread in each direction from the initial guess).
parameter_estimation_settings['mcmc_checkPointFrequency'] = None #This is only for MH, not ESS. (as of Oct 2020)
parameter_estimation_settings['mcmc_parallel_sampling'] = False #This makes completely parallelized sampling of a single sampling. syntax to use is like "mpiexec -n 5 python runfile.py" where 5 is the number of processors. Currently, the first processor's results are thrown away.  In the future, this may change.
parameter_estimation_settings['mcmc_continueSampling']  = 'auto' #This can be set to True if user would like to continue sampling from a previous result in the directory.  The mcmc_logP_and_parameter_samples.pkl file will be used.  Note that if one calls the same PE_object after mcmc sampling within a given python instance then continued sampling will also occur in that situation.


######multistart (including gridsearch)##### 
#To do a gridsearch, make multistart_initialPointsDistributionType into 'grid' and then set the two 'gridsearchSampling' variables.
#The multistart feature exports the best values to permutations_log_file.txt, and relevant outputs to permutations_initial_points_parameters_values.csv and permutations_MAP_logP_and_parameters_values.csv
parameter_estimation_settings['multistart_searchType'] = 'getLogP' #Possible searchTypes are: 'getLogP', 'doEnsembleSliceSampling', 'doMetropolisHastings', 'doOptimizeLogP', 'doOptimizeNegLogP', 'doOptimizeSSR'.  These can also be called by syntax like PE_object.doMultiStart('doEnsembleSliceSampling') in the runfile
parameter_estimation_settings['multistart_checkPointFrequency'] = None #Note: this setting does not work perfectly with ESS.
parameter_estimation_settings['multistart_parallel_sampling'] = False
parameter_estimation_settings['multistart_centerPoint'] = None #With None the centerPoint will be taken as model['InputParameterInitialGuess'] 
parameter_estimation_settings['multistart_numStartPoints'] = 0 #If this is left as zero it will be set as 3 times the number of active parameters.
parameter_estimation_settings['multistart_initialPointsDistributionType'] = 'sobol' #Can be 'uniform', 'gaussian', 'identical', 'grid', 'astroidal', 'sobol', and 'shell'. Astroidal option has a high density of points towards the center while shell has a higher density at the edges of the prior distribution.
parameter_estimation_settings['multistart_relativeInitialDistributionSpread'] = 0.866 # Normally, a user should not change this. This setting is for non-grid multistarts in sampling parameters. This scales the initial point sampling distribution's spread. By default, this setting is at 0.866. This means that non-grid multistarts will by default have initial points chosen in the interval of 2*(0.866) sigma initial distribution in each direction from the mu of the centerPoint. The sigma initial distribution of the parameters should not be confused with the sigma prior distribution of the parameters. Recognize that the initial points distribution will in general not match the prior distribution, and thus the initial points distribution has a separate sigma. However, it will be common to use the sigma of the prior distribution to define the sigma of the initial points distribution, and PEUQSE does so by default. From statistics, in the case of a uniform distribution, the standard deviation is given by (b-a)/(12^0.5).  Thus, the initial point sampling distribution's spread will be given by bounds of 2*(b-a)/3.464. If all (or any) parameters have uniform distributions, and if one desires to have an initial point distribution with the same boundaries as that uniform distribution [i.e., no rejections], then one must set relativeInitialDistributionSpread to less than or equal to 3.464/4 which is 0.866. Hence, the default is 0.866 to avoid rejections in high parameter spaces that have uniform distributions among some of the parameters. If all of the parameters have normal distributions, then 1.0 may be appropriate. Note that 0.866 means that we explore SD +/- 1.73, such that for normal distributions the 0.866 choice still results in a space encompassing over 90% of the values for a normal distribution (because a normal distribution has 90% at 1.64 and 95% at 1.96). Thus, even for normal distributions, there is little advantage to increasing this choice to 1.0
parameter_estimation_settings['multistart_gridsearchSamplingInterval'] = [] #This is for gridsearches and is in units of absolute intervals. By default, these intervals will be set to 1 standard deviaion each.  To changefrom the default, make a comma separated list equal to the number of parameters.
parameter_estimation_settings['multistart_gridsearchSamplingRadii'] = [] #This is for gridsearches and refers to the number of points (or intervals) in each direction to check from the center. For example, 3 would check 3 points in each direction plus the centerpointn for a total of 7 points along that dimension. For a 3 parameter problem, [3,7,2] would check radii of 3, 7, and 2 for those parameters.
parameter_estimation_settings['multistart_permutationsToSamples'] = True #if this is set to True, then when 'getLogP' is selected then the sampling results will be converted into a statistical sampling distribution so that posterior distribution plots and statistics can be generated. It is presently intended for use with gridsearch or uniform distribution search, for which the pseudo-sampling created is directly proportional to the posterior distribution. The pseudo-sampling created is not directly  proportional to the posterior with other multistart searches.
parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_samples'] = True #This setting does nothing. Currently the filtering is always on. This feature removes low probability tails from the posterior. This can be important for getting mu_AP. Default is true. This only has an effect if multistart_permutationsToSamples is set to True.
parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient'] = 'auto' #This can be a float or the string 'auto'. Currently (Oct 2020), 'auto' sets the value at 2.0.  The smaller the value the more aggressive the filtering. This only has an effect if multistart_permutationsToSamples is set to True. Not recommended to set this higher than 4 or 5 as the computer may run out of memory.
parameter_estimation_settings['multistart_continueSampling']  = 'auto' #This only works with multistart_permutationsToSamples. This can be set to True if user would like to continue sampling from a previous result in the directory.  The permutations_MAP_logP_and_parameters_values.pkl file will be used.  Note that if one calls the same PE_object after multistart_permutationsToSamples sampling within a given python instance then continued sampling will also occur in that situation.
parameter_estimation_settings['multistart_passThroughArgs'] = {}  #Typically, one would put here the arguments for doOptimizeNegLogP:  {'method':"Nelder-Mead", 'printOptimum':True, 'verbose':False}  To find additional details of which arguments can be used with doOptimizeNegLogP, see the function doOptimizeNegLogP in PEUQSE\InverseProblem.py
parameter_estimation_settings['multistart_calculatePostBurnInStatistics'] = True #This is mainly for multistart_permutationsToSamples. This calculates the mu_AP, among other values.
parameter_estimation_settings['multistart_keep_cumulative_post_burn_in_data'] = False
parameter_estimation_settings['multistart_exportLog'] = False #In the future, this will cause more information to be exported.

#####Plot Settings#####
plotting_ouput_settings={}
plotting_ouput_settings['setMatPlotLibAgg'] = 'auto' #This setting controls whether the matplot lib aggregator is turned on. #by default, on Windows machines this setting will become False during runtime and the plots will be generated 'normally'. By default, on Linux machines, this setting well be set to True during run time, and the matplotlib plot aggregator will be turned on.  The reason this setting exists is that on most supercomputers and clusters, which are usally linux based, the graphs will not be generated unless the aggregator is turned on.


#####Histogram Plot Settings#####
histogram_plot_settings = {}
histogram_plot_settings['histograms_as_density']= False # When true, the histograms will be shown as density with integral normalized to 1. It's good to keep this as False during analysis, then to set to True to make final graphs for presentations.  By default this is False. When this is False, for ESS or MH the histograms will show sampling frequency. When this is False, for uniform or random samplings, the histograms show a pseudo-sampling frequency with an order of magnitude proportional to multistart_gridsearch_threshold_filter_coefficient.
histogram_plot_settings['y_label'] = '' #will use defaults if ''
histogram_plot_settings['dpi'] = 220 #Dots per Inch specification. This changes the resolution of the image created. Larger values are recommended for important pictures.
histogram_plot_settings['show_initial_value'] = True #This option shows the initial value as a vertical line in the histogram.
histogram_plot_settings['show_MAP_value'] = True #This option shows the MAP as a vertical line in the histogram.
histogram_plot_settings['show_muAP_value'] = True #This option shows the muAP as a vertical line in the histogram.
histogram_plot_settings['vertical_linewidth'] = 1.5 #This value changes the linewidth of the initial value, MAP, and muAP. The units are points, where 1 point is dpi/72 pixels. Typical values are 1 for thin and 4 for thick.
histogram_plot_settings['x_label_size'] = 16
histogram_plot_settings['y_label_size'] = 16
histogram_plot_settings['axis_font_size'] = 16


#Response Plot Settings
simulated_response_plot_settings = {}
simulated_response_plot_settings['x_label'] = '' #will use defaults if ''
simulated_response_plot_settings['y_label'] = '' #will use defaults if ''
#simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
simulated_response_plot_settings['figure_name'] = 'Posterior_Simulated' #This is the default name for simulated response plots.
simulated_response_plot_settings['legend'] = True #Can be changed to false to turn off the legend.
#simulated_response_plot_settings['legendLabels'] = ['experiment', 'mu_guess', 'MAP'] here is an example of how to change the legend labels.
simulated_response_plot_settings['error_linewidth'] = 'auto' #Integer. Using "auto" or "None" sets to "20" when there is only 1 point, 1 when number of points is > 10, and "4" when number of points is between 1 and 10 and. Using '0' or 'none' will hide the error bars.
simulated_response_plot_settings['fontdict']= {'size':16} #A font dictionary can be passed in, this will be used for the axes and axes labels.

#Scatter Matrix Plot Settings
#possible dictionary fields include: dpi, figure_name, fontsize, x_label, y_label, figure_name, x_range, y_range
scatter_matrix_plots_settings ={}
scatter_matrix_plots_settings['combined_plots'] = 'auto' #True, False, or  'auto'. With 'auto', the combined plots are only created if there are 5 parameters or less.
scatter_matrix_plots_settings['individual_plots'] = 'auto' #True, False, or 'auto'. With 'auto', the individual_plots will only be created when there are more than 5 parameters.
scatter_matrix_plots_settings['dpi'] = 220 #Dots per Inch specification. This changes the resolution of the image created. Larger values are recommended for important pictures.
scatter_matrix_plots_settings['figure_name'] = 'scatter_matrix_posterior'
scatter_matrix_plots_settings['sampled_point_sizes'] = 1 #Size of sampled points in the scatter plot. Lower values are recommended (under 6).
scatter_matrix_plots_settings['cross_marker_size'] = 10 #Size of special points (MAP, muAP, and initial point) denoted by a cross. Values should be larger than sampled_point_size to increase clarity. 
scatter_matrix_plots_settings['sampled_point_transparency'] = 0.5 #Transparency of sampled points in scatter plot. This value ranges from 0 to 1. Values closer to 1 can make plots with many points look crowded and confusing. Mid values are recommended. 
scatter_matrix_plots_settings['cross_marker_transparency'] = 0.8 #Transparency of special points (MAP, muAP, and initial point) denoted by a cross. This value ranges from 0 to 1. This should be larger than sampled_point_transparency to emphasize the special point locations.
scatter_matrix_plots_settings['max_num_x_ticks'] = 'auto' #Maximum number of tick marks on x axis. This value should be an integer or the string 'auto' to let matplotlib choose the number of tick marks.
scatter_matrix_plots_settings['max_num_y_ticks'] = 'auto' #Maximum number of tick marks on y axis. This value should be an integer or the string 'auto' to let matplotlib choose the number of tick marks.
scatter_matrix_plots_settings['fontsize'] = 16 #Fontsize of labels on x and y axis. This value should be an integer. 
scatter_matrix_plots_settings['all_pair_permutations'] = True #By default, all possible permutations are created [(a,b) and (b,a)]. To get only combinations rather than all permutations, set this option to False.

#Scatter Heatmap Plot Settings
scatter_heatmap_plots_settings ={}
scatter_heatmap_plots_settings['figure_name'] = 'scatter_heatmap_posterior'
scatter_heatmap_plots_settings['dpi'] = 220 #Dots per Inch specification. This changes the resolution of the image created. Values between 220 and 600 dpi are typical.
scatter_heatmap_plots_settings['sampled_point_sizes'] = 1 #Size of sampled points in the scatter plot. Lower values are recommended (under 6).
scatter_heatmap_plots_settings['cross_marker_size'] = 10 #Size of special points (MAP, muAP, and initial point) denoted by a cross. Values should be larger than sampled_point_size to increase clarity. 
scatter_heatmap_plots_settings['sampled_point_transparency'] = 0.5 #Transparency of sampled points in scatter plot. This value ranges from 0 to 1. Values closer to 1 can make plots with many points look crowded and confusing. Mid values are recommended. 
scatter_heatmap_plots_settings['cross_marker_transparency'] = 0.8 #Transparency of special points (MAP, muAP, and initial point) denoted by a cross. This value ranges from 0 to 1. This should be larger than sampled_point_transparency to emphasize the special point locations.
scatter_heatmap_plots_settings['max_num_x_ticks'] = 'auto' #Maximum number of tick marks on x axis. This value should be an integer or the string 'auto' to let matplotlib choose the number of tick marks.
scatter_heatmap_plots_settings['max_num_y_ticks'] = 'auto' #Maximum number of tick marks on y axis. This value should be an integer or the string 'auto' to let matplotlib choose the number of tick marks.
scatter_heatmap_plots_settings['fontsize'] = 16 #Fontsize of labels on x and y axis. This value should be an integer. 
scatter_heatmap_plots_settings['all_pair_permutations'] = True #By default, all possible permutations are created [(a,b) and (b,a)]. To get only combinations rather than all permutations, set this option to False.

#contour plots# / #mumpce plots#
#model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': r'$E_{a1}$'},
#{'parameter_number': 1, 'parameter_name': r'$E_{a2}$'},
#{'parameter_number': 2, 'parameter_name': r'$log(A_{1})$'},
#{'parameter_number': 3, 'parameter_name': r'$log(A_{2})$'},
#{'parameter_number': 4, 'parameter_name': r'$\gamma_{1}$'},
#{'parameter_number': 5, 'parameter_name': r'$\gamma_{2}$'}])
contour_plot_settings = {}
contour_plot_settings['active_parameters'] = [] #Blank by default: gets populated with all parameters (or reduced parameters) if left blank. Warning: trying to set this manually while using the reduced parameters feature is not supported as of April 2020.
contour_plot_settings['parameter_pairs'] = [] #This will accept either strings (for variable names) or integers for positions. #This sets which parameters to plot contours for. By default, all pairs are plotted. For example,  [[0, 1], [1, 2],[2, 3],[3, 4],[4, 5]] 
contour_plot_settings['figure_name'] = 'PosteriorContourPlots'
contour_plot_settings['individual_plots'] = 'auto' #True, False, or 'auto'. With 'auto', the individual_plots will always be created.
contour_plot_settings['combined_plots'] = 'auto' #True, False, or  'auto'. With 'auto', the combined plots are only created if there are 5 pairs or less.
contour_plot_settings['zoom_std_devs'] = 2.5 #how zoomed in the image is.
contour_plot_settings['fontsize']=16  #sets the fontsize for everything except the colorbars. Can be an integer or the word 'auto', or the word "None". Should change space_between_subplots if fontsize is changed. 
contour_plot_settings['space_between_subplots'] = 0.40 #Typically a value between 0.20 and 5.0. Set to 0.40 by default. Should be changed when font size is changed. Fontsize 'auto' tends to make small fonts which needs smaller values like 0.20.
contour_plot_settings['cmap_levels'] = 4   #This is the number of contour levels.
contour_plot_settings['max_num_y_ticks'] = 'auto'  #adjusts number of y ticks (actually sets a maximum number of them). #max_num_y_ticks and max_num_x_ticks must be either a string ('auto') or an integer (such as 4, either without string or with integer casting like int('5')). This feature is recommended.  #Note that this is a *request* When it's not fulfilled exactly, the user can play with the number.
contour_plot_settings['max_num_x_ticks'] = 'auto'  #adjusts number of x ticks (actually sets a maximum number of them). #max_num_y_ticks and max_num_x_ticks must be either a string ('auto') or an integer (such as 4, either without string or with integer casting like int('5')).This feature is recommended. #Note that this is a *request* When it's not fulfilled exactly, the user can play with the number.
contour_plot_settings['num_pts_per_axis'] = 500 #This sets the resolution of the contours.
contour_plot_settings['dpi'] = 220 #Dots per Inch specification. This changes the resolution of the image created. Values between 220 and 600 dpi are typical.
contour_plot_settings['x_ticks'] = 'auto' #feed in an array of numbers directly. Not recommended to change.
contour_plot_settings['y_ticks'] = 'auto' #feed in an array of numbers directly. Not recommended to change.
contour_plot_settings['axis_limits'] = 'auto' #Feed in list of [x_min, x_max, y_min, y_max]. This is appropriate to use. If a list of lists is provided, then the individual_plots will each receive the appropriate axis_limits.
contour_plot_settings['contours_normalized']=True #This sets the scales on the color bars to 1.0.  Changing to False shows absolute density values for the posterior and prior. With all default settings, shows contours at 0.2, 0.4, 0.6., 0.8
contour_plot_settings['center_on']='all' # #can be 'all', 'prior' or 'posterior'. 
contour_plot_settings['colorbars']=True #can be true or false.
contour_plot_settings['colormap_posterior_customized'] = 'auto' #can also be 'Oranges' for example. #accepts a string (matplotlib colormap names, like 'Greens') or a list of tuples with 0-to-1 and colornames to interpolate between. For example, the default right now is:  [(0,    '#00FFFF'),(1,    '#0000FF')]. The tuple could have 0, 0.7, and 1, for example. #colors can be obtained from: https://www.htmlcsscolor.com/hex/244162  
contour_plot_settings['colormap_prior_customized'] = 'auto' #can also be 'Greens' for example. #accepts a string (matplotlib colormap names, like 'Oranges') or a list of tuples with 0-to-1 and colornames to interpolate between. For example, the default right now is:  [(0,    '#FFFF00'),(1,    '#FF0000')]. The tuple could have 0, 0.7, and 1, for example. #colors can be obtained from: https://www.htmlcsscolor.com/hex/244162  
#See the file mumpce_custom_plotting_example.py for the full set of arguments that can be provided inside contour_plot_settings.

####Design Of Experiments####
#The design of experiments feature is used with syntax like PE_object.designOfExperiments()   
#This feature modulates the parameters to see how much info gain there would be in different parts of condition space (using synthetic data).
#The normal usage as of May 26 2020 is to first use the indpendent variables feature, which must be used, then fill the doe_settings below.
#Then call PE_object.doeParameterModulationPermutationsScanner()
#If a single parameter modulation grid is going to be used, one can instead define the independent variable grid and call like this: PE_object.createInfoGainPlots(plot_suffix="manual")
#For a real usage, see Example13doeFunctionExample
#Note: after calling either of these functions, the variables populated are PE_object.info_gains_matrices_array and PE_object.info_gain_matrix.  So if somebody wants to export these after using the functions, one can cycle across each info gains matrix inside PE_object.info_gains_matrices_array and export to csv.
#A key is printed out inside of Info_gain_parModulationGridCombinations.csv

doe_settings = {} #To use the design of experiments feature the independent variables feature **must** be used.
doe_settings['info_gains_matrices_array_format'] = 'xyz' #options are 'xyz' and 'meshgrid'.  Images are only ouput when scanning two independent variables. If using more than two independent variables, you will need to use the 'xyz' format and will need to analyze the final info_gains_matrices_array written to file directly. Note this variable must be set before running the doe command. You cannot change the format of the info_gains_matrices_array afterwards because the way the sampling is stored during a run is change based on this setting.

doe_settings['info_gains_matrices_multiple_parameters'] = 'sum' #The possible values are 'sum' or 'each'. 'sum' is the default such that there is one averaged infogain matrix exported per modulation. The 'each' choice exports info_gains for **each** parameter per modulation (and also exports the sums). #This feature is (for now) only for KL_divergence.



#Now we will define how big of a modulation / offset we want to apply to each parameter.
#doe_settings['parameter_modulation_grid_center'] #We do NOT create such a variable. The initial guess variable is used, which is the center of the prior if not filled by the user.
doe_settings['parameter_modulation_grid_interval_size'] = [] #This must be 1D array/list with length of number of parameters.  These are all relative to the standard deviation of the prior of that parmaeter. Such as [1,1]. The default will set everything to 1.  #If you wish to manually change this setting, you should use all non-zero values for this setting, even if you plan to not modulate that parmeter (which parameters will be modulated are in the num_intervals variable below).
doe_settings['parameter_modulation_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of parameters. #Such as [1,1]. The default will be set to all 1. #This is the number of steps in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the parameter to vary.
doe_settings['parameter_modulation_grid_checkPointFrequency'] = None #None means no checkpoints. Recommended values are None or 1.
doe_settings['parallel_parameter_modulation'] = False  #this parallelizes the modulation of parameters. It is not compatible with parallel_conditions_exploration, so only one should be used at a time.

#Now we will define the *conditions space* to explore (to find the highest info gain), which will be done for *each* modulation.
#Note that this means that responses['independent_variables_values'] must be used, AND it must be fed into the model simulation file as a connected variables array. (See  Example13doeFunctionExample directory runfiles).
doe_settings['independent_variable_grid_center'] = [] #This must be a 1D array/list with length of number of independent variables.  It's a central condition a grid will be made around. Example: [500, 0.5]
doe_settings['independent_variable_grid_interval_size'] = [] #This must be a 1D array/list with length of number of independent variables.  this how big of each step will be in each direction/dimension (it is the grid spacing).  You should always use non-zero values for this setting, even if you plan to not vary that independent variable. Like [100,0.2]
doe_settings['independent_variable_grid_num_intervals'] = [] #This must be a 1D array/list with length of number of independent variables. #This is the number of steps in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the condition to vary. Example: [2,2] 
doe_settings['independent_variable_grid_checkPointFrequency'] = None #None means no checkpoints. Recommended values are None or 1.
doe_settings['parallel_conditions_exploration'] = False  #this parallelizes the modulation of the conditions exploration. It is not compatible with parallel_parameter_modulation, so only one should be used at a time.


doe_settings['on_the_fly_conditions_grids'] = True #Oct 2020: do not change this setting. #Values are True or False. This makes the independent variable grid each time. This costs more processing but less memory. As of April 2020 the other option has not been implemented but would just require making the combinations into a list the first time and then operating on a copy of that list.