import sys; sys.path.append('../../')
sys.path.append('../../../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    #####This file is like the other Example 15 runfiles, and has the two_site_ratio as giving the percent of sites that is site 2 (as a decimal) and this is the first parameter in the parameter vector.
    
    import processing_functions_tpd_odeint_two_site_NineParameters
    observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    times, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint_two_site_NineParameters.import_experimental_settings_single(observed_data_Filename)
    #experiments_datarame = pd.read_csv(observed_data_Filename)    
    
    
    UserInput.responses['responses_abscissa'] = times
    UserInput.responses['responses_observed'] = responses_observed
    UserInput.responses['responses_observed_uncertainties'] = observedResponses_uncertainties
    
    #We are going to use the built in transform to improve the optimization.
    UserInput.responses['data_overcategory'] = 'transient_kinetics'
    UserInput.responses['response_types']=['P'] #need a categorization for each response dimension.
    UserInput.responses['response_data_types']=['r'] #need a categorization for each response dimension.
    
    
    UserInput.simulated_response_plot_settings['x_label'] = 'time (s)'
    UserInput.simulated_response_plot_settings['y_label'] = r'Rate of Desorption (ML/s)'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example16a' #This creates the filename, also.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'scalingFactor':'scalingFactor', 'backgroundOffset':'backgroundOffset', 'site2Ratio':'site2Ratio','Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$'}
    UserInput.model['InputParameterPriorValues'] = [ 1.0, 0.0, 0.50, 41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    UserInput.model['InputParametersPriorValuesUncertainties'] = [ 0.10, 0.005, 0.50/3, 20, 20, 2, 2, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterInitialGuess'] = [1.31127639e+00,2.41324344e-03,0.395877298,4.03158778e+01,2.84225571e+01,2.59928937e+01,1.54584948e+01,9.48684766e-09,8.29500179e-02]
    #UserInput.model['InputParameterPriorValues_upperBounds'] = [ None, None, 1.0, None, None, None, None, None, None]
    #UserInput.model['InputParameterPriorValues_lowerBounds'] = [ 0, None, 0, 0, 0, None, None, 0, 0]
    
    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    
    #InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
    UserInput.model['simulateByInputParametersOnlyFunction'] = processing_functions_tpd_odeint_two_site_NineParameters.TPR_simulationFunctionWrapperNineParameters #This must simulate with *only* the parameters listed above, and no other arguments.
    UserInput.model['simulationOutputProcessingFunction'] = None  #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    
    
    UserInput.model['reducedParameterSpace'] = [3,5]
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "std"
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['exportLog'] = True
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 1000
     
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 20 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 1#Normally use 100 for example 3a
    UserInput.parameter_estimation_settings['mcmc_length'] = 100000#Normally use 10000 for example 3a.
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.01
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'] = 0
    UserInput.contour_plot_settings['contours_normalized'] = True
    UserInput.contour_plot_settings['colorbars'] = False
    UserInput.contour_plot_settings['fontsize']=20
    #UserInput.contour_plot_settings['center_on'] = 'prior'
    UserInput.contour_plot_settings['axis_limits'] = [0,80, 0, 60]
    UserInput.contour_plot_settings['max_num_x_ticks']= 4

    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
#    #Now we do parameter estimation.
#    PE_object.doMetropolisHastings()
   # PE_object.doSinglePoint()
#    PE_object.createAllPlots() #This function calls each of the below functions.
    
#    PE_object.doOptimizeNegLogP(method="Nelder-Mead", printOptimum=True, verbose=True)

#    PE_object.doGridSearch('getLogP')
    #PE_object.doGridSearch('doMetropolisHastings')
#    PE_object.doGridSearch('doOptimizeNegLogP', verbose = True,gridSamplingNumOfIntervals = [], passThroughArgs={'method':'BFGS'})
    #PE_object.doOptimizeNegLogP(method="Nelder-Mead", printOptimum=True, verbose=False, maxiter=5000)
    PE_object.doSinglePoint()
    #PE_object.doGridSearch('doOptimizeNegLogP', gridSamplingAbsoluteIntervalSize=UserInput.model['InputParametersPriorValuesUncertainties'], gridSamplingNumOfIntervals=[0,0,1,1,0, 1,0,1,0,1,0], passThroughArgs={"method":"Nelder-Mead", "maxiter":100, "verbose":False})#, "maxiter":1000, "verbose":False})
    print(PE_object.map_parameter_set, PE_object.map_logP)
    UserInput.simulated_response_plot_settings['fontdict']= {'size':15}
    UserInput.simulated_response_plot_settings['legend']=False
    #UserInput.simulated_response_plot_settings['legendLabels'] = [' f', 'g ', 'r ']
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
    #PE_object.createSimulatedResponsesPlots()
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.