import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    import observed_values_00  #Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    import simulation_model_00 #Simple example.
        
    UserInput.responses['responses_abscissa'] = observed_values_00.observed_data_x_values
    UserInput.responses['responses_observed'] = observed_values_00.observed_data_y_values
    UserInput.responses['responses_observed_uncertainties'] = observed_values_00.observed_data_y_values_uncertainties

    
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$time (s)$'
    UserInput.simulated_response_plot_settings['fontdict'] = {'size':16}
    
    
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [100, 200] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    
    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.

    UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] = 1000
    UserInput.parameter_estimation_settings['mcmc_length'] = 10000 #10000 is the default.


    #UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 This can be useful for testing.
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = CKPQ.parameter_estimation(UserInput)
    PE_object.doMetropolisHastings()
    
    
    ####NORMALLY WE WOULD NOW CREATE ALL PLOTS, BUT IN THIS EXAMPLE WE WILL DELAY DOING SO IN ORDER TO SHOW HOW TO SAVE AND LOAD A PE_OBJECT TO AND FROM A FILE.
    #PE_object.createAllPlots()

    #to save a PE_object for later
    PE_object.save_to_dill("PE_object_00a1")
    
    #to load a PE_object after some time, first one has to put (any) UserInput to create a PE_object, then to load from file.
    #these two steps can be done in a different python file. A different PE_object name is being used to emphasize that this process can be done from a different python file.
    PE_object2 = CKPQ.parameter_estimation(UserInput)
    PE_object2 = PE_object2.load_from_dill("PE_object_00a1")
    
    
    PE_object2.createAllPlots() #This function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()
