import pandas as pd
import sys; sys.path.append('../../');  import PEUQSE as PEUQSE
import cantera as ct
import cantera.ck2cti as ck2cti
import PEUQSE.simulationDriver.canteraSimulate
import PEUQSE.simulationDriver.canteraKineticsParametersParser 
import numpy as np




if __name__ == "__main__":
    import PEUQSE.UserInput as UserInput
    
    import processing_functions_tpd_odeint #We just want to import the experimental data.
    observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    #times, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint.import_experimental_settings(observed_data_Filename) <-- this is non-integral version.
    times, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint.import_integrals_settings(observed_data_Filename)
    
    
    
    UserInput.responses['responses_abscissa'] = times
    UserInput.responses['responses_observed'] = responses_observed
    UserInput.responses['responses_observed_uncertainties'] = observedResponses_uncertainties


    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'Ea_1':r'$E_{a1}$','Ea_2':r'$E_{a2}$','log_A1':r'$log(A_{1})$','log_A2':r'$log(A_{2})$','gamma1':r'$\gamma_{1}$','gamma2':r'$\gamma_{2}$', 'modA1':'modA1', 'modA2':'modA2', 'modA3':'modA3', 'modA4':'modA4', 'modA5':'modA5', 'modA6':'modA6', 'modA7':'modA7' ,  'modE1':'modE1', 'modE2':'modE2', 'modE3':'modE3', 'modE4':'modE4', 'modE5':'modE5', 'modE6':'modE6', 'modE7':'modE7' }
    UserInput.model['InputParameterPriorValues'] = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1, # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
                                                     -1,-1,-1,-1,-1,-1,0 , #A modifiers.
                                                     60000,50000,40000,30000,20000,10000,0] #E modifiers.

    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    UserInput.model['InputParametersPriorValuesUncertainties'] = [20, 20, 13, 13, 0.1, 0.1,
                                                     1,1,1,1,1,1,1 , #A modifiers.
                                                     10000,10000,10000,10000,10000,10000,10000 #E modifiers.
                                                     ] 
                                                                   #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D variance is acceptable.

    UserInput.model['InputParameterInitialGuess'] = [20, 20, 13, 13, 0.0, 0.0,
                                                     0,0,0,0,0,0,0 , #A modifiers.
                                                     0,0,0,0,0,0,0 #E modifiers.
                                                     ] 
                                                                   #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D variance is acceptable.

    
    from model_functions_example6 import integrated_cantera_simulation_wrapper_example6
    UserInput.model['simulateByInputParametersOnlyFunction'] = integrated_cantera_simulation_wrapper_example6 #This must simulate with *only* the parameters listed above, and no other arguments.
    UserInput.model['simulationOutputProcessingFunction'] = None #Optional: a function to process what comes out of the simulation Function and then return an observable vector.

    UserInput.simulated_response_plot_settings['x_label']= 'time (s)'
    UserInput.simulated_response_plot_settings['y_label'] = 'Integral'
    
    UserInput.parameter_estimation_settings['verbose'] = False 
    
     
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10000 
    UserInput.parameter_estimation_settings['mcmc_length'] = 30000
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.


    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    #Now we do parameter estimation.
    #PE_object.doMetropolisHastings()
    PE_object.doOptimizeNegLogP()
    #[map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, samples_simulatedOutputs, logP] = PE_object.doMetropolisHastings()
    
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()
    #TODO: call the mum_pce plotting objects, which will be PE_object.createContourGraphs() or something like that.