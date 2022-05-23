import sys; sys.path.append('../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    import model_functions_example5 #We will need both the functions and the experimental data.
    
    observed_x_values, responses_observed, observedResponses_uncertainties = model_functions_example5.observed_x_values, model_functions_example5.responses_observed, model_functions_example5.observedResponses_uncertainties
    
    UserInput.responses['responses_abscissa'] = observed_x_values
    UserInput.responses['responses_observed'] = responses_observed
    UserInput.responses['responses_observed_uncertainties'] = observedResponses_uncertainties


    #UserInput.model['parameterNamesAndMathTypeExpressionsDict'] #We will just take the default names for the parameters.
    UserInput.model['InputParameterPriorValues'] = [0,0,0,0,0,0,0] #E modifiers will start at 0.
    UserInput.model['InputParametersPriorValuesUncertainties'] = [10000,10000,10000,10000,10000,10000,10000] #E modifiers uncertainties. To assume no covariance, a 1D standard deviation is acceptable.
    
    #Ideally, we should use half-normal distributions or similar. However, we will simply make a lower bound of 0 for the parameters.
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [0,0,0,0,0,0,0]
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = model_functions_example5.cantera_simulation_wrapper_example5 #This must simulate with *only* the parameters listed above, and no other arguments.

    UserInput.simulated_response_plot_settings['x_label']= 'time (s)'
    UserInput.simulated_response_plot_settings['y_label'] = 'Integral'
    
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] = 10
     
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 100
    UserInput.parameter_estimation_settings['mcmc_length'] = 300
    # UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.


    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    #Now we do parameter estimation.
    PE_object.doMetropolisHastings()
    #[map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, samples_simulatedOutputs, logP] = PE_object.doMetropolisHastings()
    
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.