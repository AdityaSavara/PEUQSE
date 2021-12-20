import sys; sys.path.append('../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    import observed_values_00  #Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    import simulation_model_00 #Simple example.
        
    UserInput.responses['responses_abscissa'] = observed_values_00.observed_data_x_values
    UserInput.responses['responses_observed'] = observed_values_00.observed_data_y_values
    UserInput.responses['responses_observed_uncertainties'] = observed_values_00.observed_data_y_values_uncertainties

    
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$time (s)$'
    

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [-1, -1] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterPriorValues_upperBounds'] = [1E6, 1E6] 
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [-1E6, -1E6]
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "off"
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0
        #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.
    #UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10000
    #UserInput.parameter_estimation_settings['mcmc_length'] = 100000 #The uninformed prior int his example has a "bad" MCMC walker so requires lots of sampling to converge.
    UserInput.parameter_estimation_settings['multistart_checkPointFrequency'] = 10000 #This example is long enough that it's good to get updates.
    
    UserInput.parameter_estimation_settings['multistart_searchType'] = 'getLogP'
    UserInput.parameter_estimation_settings['multistart_initialPointsDistributionType'] = 'uniform'
    UserInput.parameter_estimation_settings['multistart_exportLog'] = True
    UserInput.parameter_estimation_settings['multistart_gridsearch_threshold_filter_coefficient'] = 2.0 #The lower this is, the more the points become filtered. It is not recommended to go below 2.0.
    UserInput.parameter_estimation_settings['multistart_numStartPoints'] = 1000000
    UserInput.parameter_estimation_settings['multistart_relativeInitialDistributionSpread'] = 0.50
    
    
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    PE_object.doMultiStart()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()