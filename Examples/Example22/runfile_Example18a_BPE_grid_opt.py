import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    import processing_functions_tpd_odeint
    observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    times, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint.import_experimental_settings_single(observed_data_Filename)
    #experiments_datarame = pd.read_csv(observed_data_Filename)    
    processing_functions_tpd_odeint.initial_concentrations_array = [1.0] #Trying to force a single coverage.
    
    UserInput.responses['responses_abscissa'] = times
    UserInput.responses['responses_observed'] = responses_observed
    UserInput.responses['responses_observed_uncertainties'] = observedResponses_uncertainties
    
    #We are going to use the built in transform to improve the optimization.
    UserInput.responses['data_overcategory'] = 'transient_kinetics'
    UserInput.responses['response_types']=['P'] #need a categorization for each response dimension.
    UserInput.responses['response_data_type']=['r'] #need a categorization for each response dimension.    

    
    UserInput.simulated_response_plot_settings['x_label'] = 'time (s)'
    UserInput.simulated_response_plot_settings['y_label'] = r'Rate of Desorption (ML/s)'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example' #This creates the filename, also.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'scalingFactor':'scalingFactor', 'verticalOffset':'verticalOffset', 'Ea_1':r'$E_{a1}$','log_A1':r'$log(A_{1})$','gamma1':r'$\gamma_{1}$',  'modEa1':'modEa1', 'modEa2':'modEa2', 'modEa3':'modEa3', 'modEa4':'modEa4', 'modEa5':'modEa5', 'modEa6':'modEa6' }#, 'modEa5':'modEa5', 'modEa6':'modEa6', 'modEa7':'modEa7', 'modEa8':'modEa8', 'modEa9':'modEa9', 'modEa10':'modEa10'}
    UserInput.model['InputParameterPriorValues'] = [1.0, 0.0, 40.0, 13.0, 0.1, 
                   0.0, 0.0,   0.0, 0.0, 0.0, 0.0]  #0, 0.2, 0.4, 0.6, 0.8, 1.0
    UserInput.model['InputParametersPriorValuesUncertainties'] = [.1, 0.005, 20, 2, 0.3, 
                   0.1, 0.1,      0.1, 0.1, 0.1, 0.1] #,        10, 10,      10, 10,          10, 10] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterInitialGuess'] =[1.0, 0.0, 40.0, 13.0, 0.0, 
                   0.3, 0.0,   0.0, 0.0, 0.0, 0.0] 
    
    UserInput.model['InputParameterPriorValues_upperBounds'] = [ None, None,  None, None, None,     None, None, None, None, None, None]
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [ 0, None, 0, 0, None,     None, None, None, None, None, None]
    
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = processing_functions_tpd_odeint.TPR_internalPiecewiseSimulationFunctionWrapperScaledAndOffset #This must simulate with *only* the parameters listed above, and no other arguments.
    UserInput.model['simulationOutputProcessingFunction'] = processing_functions_tpd_odeint.no_log_wrapper_func  #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    UserInput.parameter_pairs_for_contour_plots=[[2,3],[2,4]]
        
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "off"                                                                                                 
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['exportLog'] = False
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 1
     
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 1
    UserInput.parameter_estimation_settings['mcmc_length'] = 1001
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.05
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 1000 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.

    UserInput.contour_settings_custom['contours_normalized'] = True
    

    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = CKPQ.parameter_estimation(UserInput)
    
#    #Now we do parameter estimation.
#    PE_object.doMetropolisHastings()
   # PE_object.doSinglePoint()
#    PE_object.createAllPlots() #This function calls each of the below functions.
    
#    PE_object.doOptimizeNegLogP(method="Nelder-Mead", printOptimum=True, verbose=True)

#    PE_object.doGridSearch('getLogP')
    #PE_object.doGridSearch('doMetropolisHastings')
#    PE_object.doGridSearch('doOptimizeNegLogP', verbose = True,gridSamplingRadii = [], passThroughArgs={'method':'BFGS'})
    
    PE_object.doGridSearch('doOptimizeNegLogP', gridSamplingAbsoluteIntervalSize=UserInput.model['InputParametersPriorValuesUncertainties'], gridSamplingNumOfIntervals=[0,0,1,1,0, 1,0,1,0,1,0], passThroughArgs={"method":"Nelder-Mead", "maxiter":100, "verbose":False})#, "maxiter":1000, "verbose":False})
    print(PE_object.map_parameter_set, PE_object.map_logP)
    
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
    #PE_object.createSimulatedResponsesPlots()
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.