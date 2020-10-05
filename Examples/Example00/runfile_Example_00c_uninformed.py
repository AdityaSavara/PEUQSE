import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    import simulation_model_00
    simulation_model_00.x_values_for_data = [600,1100,1400] #Setting the x_values_for_data
    
    
    UserInput.responses['responses_abscissa'] = simulation_model_00.x_values_for_data
    UserInput.responses['responses_observed'] = [360500, 580500, 1620500]
    UserInput.responses['responses_observed_uncertainties'] = [200000, 300000, 200000]    
    
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$time (s)$'
    UserInput.simulated_response_plot_settings['error_linewidth'] = 4

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [-1, -1] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterPriorValues_upperBounds'] = [1E6, 1E6] 
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [-1E6, -1E6]
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "off"
        #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10000
    UserInput.parameter_estimation_settings['mcmc_length'] = 1000000 #The uninformed prior int his example has a "bad" MCMC walker so requires lots of sampling to converge.
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 10000 #This example is long enough that it's good to get updates.
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = CKPQ.parameter_estimation(UserInput)
    PE_object.doMetropolisHastings()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlot()