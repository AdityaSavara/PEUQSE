import sys; sys.path.insert(0, '../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":    
    import processing_function_two_response #This will provide the "simulation" function.
    import numpy as np

    UserInput.responses['responses_abscissa'] = np.array([[0], [1]]) # These represent theta_1 and theta_2 abscissa points. Note that this is *nested* it is a multi-dimensional response, and produces posterior distributions for each response.
    UserInput.responses['responses_observed'] = np.array([[2], [3]]) # Note that this is *nested* it is a multi-dimensional response, and produces two posteriors.
    UserInput.responses['responses_observed_uncertainties'] = np.array([[1], [1]])

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'theta_1':r'$\theta_{1}$','theta_2':r'$\theta_{2}$'}
    UserInput.model['InputParameterPriorValues'] = [1, 5] 
    UserInput.model['InputParametersPriorValuesUncertainties'] = np.array([[1,0.0], [0.0,1]]) #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterInitialGuess'] = [1, 5] #This is where the mcmc chain will start.
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = processing_function_two_response.split_to_separated_lists #This must simulate with *only* the parameters listed above, and no other arguments.
    
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example_two_response' #This creates the filename, also.
    

    
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] = False
    UserInput.parameter_estimation_settings['checkPointFrequency'] = None
    UserInput.parameter_estimation_settings['exportLog'] = True
     
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'MAP_finding'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 100
    UserInput.parameter_estimation_settings['mcmc_length'] = 10000
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 1.0
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "std"
    UserInput.contour_plot_settings['parameter_pairs'] = [[0, 1]]
    UserInput.contour_plot_settings['contours_normalized'] = False

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
