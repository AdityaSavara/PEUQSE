import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    import simulation_model_00
    simulation_model_00.x_values_for_data = [600,1100,1400] #Setting the x_values_for_data
    
    
    UserInput.responses['responses_abscissa'] = simulation_model_00.x_values_for_data
    UserInput.responses['responses_observed'] = [360500, 580500, 1620500]
    UserInput.responses['responses_observed_uncertainties'] = [200000, 300000, 200000]    
    
    # UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    # UserInput.simulated_response_plot_settings['y_label'] = r'$timee (s)$'


    #UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [100, 200] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.
   
    UserInput.parameter_estimation_settings['mcmc_length'] = 100
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 1
    UserInput.parameter_estimation_settings['scaling_uncertainties_type']="std"
    
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = CKPQ.parameter_estimation(UserInput)
    PE_object.doMetropolisHastings()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlot()
    
    
    
    #Now we do parameter estimation.
    #PE_object.doGridSearch('getLogP', verbose = False)
    PE_object.doMetropolisHastings()
    #PE_object.doOptimizeNegLogP(method="BFGS", printOptimum=True, verbose=True)
    #[map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, samples_simulatedOutputs, logP] = PE_object.doMetropolisHastings()
    
    
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.