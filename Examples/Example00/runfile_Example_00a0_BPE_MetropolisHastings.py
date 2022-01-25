import sys; sys.path.append('../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    import observed_values_00  #Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    import simulation_model_00 #Simple example.
        
    #Provide the observed X values and Y values and uncertainties -- all should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
    UserInput.responses['responses_abscissa'] = observed_values_00.observed_data_x_values
    UserInput.responses['responses_observed'] = observed_values_00.observed_data_y_values
    UserInput.responses['responses_observed_uncertainties'] = observed_values_00.observed_data_y_values_uncertainties
   
    #Optional: provide labels for the responses axes and parameter names.
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$time (s)$'
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    
    #Provide the prior distribution and uncertainties of the individual parameters.
    UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [100, 200] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    #Provide a function that returns simulated values -- must of the same form as observed values, should be arrays or lists with nesting like [[1,2,3]] or [[1,2,3,4],[4,5,6,6]]
    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.

    #mcmc length should typically be on the order of 10,000 per parameter. By default, the burn in will be the first 10% of the mcmc length.
    UserInput.parameter_estimation_settings['mcmc_length'] = 10000 #10000 is the default.

    #After filinlg the variables of the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    #Now we can do the mcmc!
    PE_object.doMetropolisHastings()
    #Another option would be PE_object.doEnsembleSliceSampling(), one can also do grid search or an astroidal distribution search.
    
    #Finally, create all plots!
    PE_object.createAllPlots()
    ####NORMALLY WE WOULD NOW CREATE ALL PLOTS, BUT IN THIS EXAMPLE WE WILL DELAY DOING SO IN ORDER TO SHOW HOW TO SAVE AND LOAD A PE_OBJECT TO AND FROM A FILE.
    
    #Optionally, one can save a PE_object for later
    try:
        PE_object.save_to_dill("PE_object_00a0")
    except:
        pass
    
    #to load a PE_object after some time, first one has to put (any) UserInput to create a PE_object, then to load from file.
    #these two steps can be done in a different python file. A different PE_object name is being used to emphasize that this process can be done from a different python file.
    PE_object2 = PEUQSE.parameter_estimation(UserInput)
    PE_object2 = PE_object2.load_from_dill("PE_object_00a0")
    
    
    PE_object2.createAllPlots() #The createAllPlots function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()
