import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":    
    import simulationFunctionExample14Tp #This will provide the "simulation" function.
    import numpy as np

    UserInput.responses['responses_abscissa'] = np.array([[1]]) # 
    UserInput.responses['responses_observed'] = np.array([[230]]) #
    UserInput.responses['responses_observed_uncertainties'] = np.array([[10]])

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'logA':r'log(A/s^-1)','Ea':r'Ea'}#,'Theta0':'Theta0', 'beta_H':'beta_H', 'n':'n'}
    UserInput.model['InputParameterPriorValues'] = [13, (37.63*1000)] 
    UserInput.model['InputParametersPriorValuesUncertainties'] = np.array([2, 25.30*1000]) #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    UserInput.model['InputParameterInitialGuess'] = [12, 45*1000] #This is where the mcmc chain will start.
    UserInput.model['responses_simulation_uncertainties'] = np.array([[25]])
    
    UserInput.model['simulateByInputParametersOnlyFunction'] = simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper #This must simulate with *only* the parameters listed above, and no other arguments.
        
    
    #UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example_two_response' #This creates the filename, also.
    

    
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] = False
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 1000
    UserInput.parameter_estimation_settings['exportLog'] = True
    UserInput.parameter_estimation_settings['mcmc'] = True 
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 1000
    UserInput.parameter_estimation_settings['mcmc_length'] = 20000
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 1.0
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "off"
    UserInput.parameter_pairs_for_contour_plots = [[0, 1]]
    #UserInput.contour_settings_custom['contours_normalized'] = True

    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = CKPQ.parameter_estimation(UserInput)
    
    #Now we do parameter estimation.
    PE_object.doMetropolisHastings()
    #[map_parameter_set, mu_AP_parameter_set, stdap_parameter_set, evidence, info_gain, samples, samples_simulatedOutputs, logP] = PE_object.doMetropolisHastings()
    
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlot()
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.

    # PE_object = CKPQ.parameter_estimation(UserInput)
    # PE_object.doOptimizeNegLogP(method="BFGS", printOptimum=True, verbose=False)
    # PE_object.createAllPlots()
    print("MAP parameters:", PE_object.map_parameter_set)
    print("mu_AP muap_parameter_set:", PE_object.mu_AP_parameter_set)
    print("std_AP stdap_parameter_set:", PE_object.stdap_parameter_set)
    print("posterior_cov_matrix:", PE_object.posterior_cov_matrix)
    
    
    ####BELOW WE TRY TO FIND THE PRIOR T_p AND ALSO THE POSTERIOR T_p#####
    import simulationFunctionExample14Tp
    print(simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(UserInput.model['InputParameterPriorValues']))
    upperRateAddition = np.array(-1*PE_object.stdap_parameter_set[0],PE_object.stdap_parameter_set[0])
    lowerRateAddition = np.array(1*PE_object.stdap_parameter_set[0],-1*PE_object.stdap_parameter_set[0])
    upperRateTp = (simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(UserInput.model['InputParameterPriorValues'])+upperRateAddition)
    lowerRateTp = (simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(UserInput.model['InputParameterPriorValues'])+lowerRateAddition)
    #Take the average of the difference and make it orthogonal to the model uncertainty of 25K.
    uncertaintyForTp = (((upperRateTp-lowerRateTp)/2)**2+25**2)**0.5
    print(uncertaintyForTp)
    
    print(simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(PE_object.mu_AP_parameter_set))
    upperRateTp=(simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(PE_object.mu_AP_parameter_set)+upperRateAddition)
    lowerRateTp=(simulationFunctionExample14Tp.getTpFromKineticParametersAndInitialCoverageWrapper(PE_object.mu_AP_parameter_set)+lowerRateAddition)
    uncertaintyForTp = (((upperRateTp-lowerRateTp)/2)**2+25**2)**0.5
    print(uncertaintyForTp)