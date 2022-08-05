# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:22:42 2022

@author: NWT

This program uses PEUQSE to obtain realistic parameters with uncertainties to
better fit data obtained from a paper studying the corrosion of 316L stainless
steel 

"""

import sys; sys.path.insert(0, '../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    import function_316
    import data_316 
    import numpy as np
        
    # Provided distances and Cr concetrations in the samples.
    UserInput.responses['responses_abscissa'] = data_316.distances
    UserInput.responses['responses_observed'] = data_316.concentrations
    UserInput.responses['responses_observed_uncertainties'] = data_316.errors
   
    # Labels for x and y axis, as well as parameter names provided.
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (um)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$Concentration (wt\%)$'
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'d_eff_cr','b':'init_cr_conc','c':'surface_conc'}
    
    # Provided the prior distribution and uncertainties of the individual parameters.
    UserInput.model['InputParameterPriorValues'] = [np.log10(4.2E-19), 16.825, data_316.concentrations[0]]
    UserInput.model['InputParametersPriorValuesUncertainties'] = [10, 1.0, data_316.errors[0]]
    
    
    # Optional bound setting lines for finding uninformed parameters.
    UserInput.model['InputParameterPriorValues_upperBounds'] = [np.log10(1.0E-17), 18.0, 18.0] 
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [np.log10(1.0E-19), 16.0, 0.0]
    

    # Provides simulation function for Cr concentration throughout a sample
    UserInput.model['simulateByInputParametersOnlyFunction'] = function_316.simulation_function_using_log_a_wrapper

    # Reduced sample size needed for EnsembleSliceSampling() due to single-mode data
    UserInput.parameter_estimation_settings['mcmc_length'] = 10000 # 10000 is the default.
    numWalkers = 24
    UserInput.parameter_estimation_settings['mcmc_nwalkers'] = numWalkers

    
    # take the MAP directly from the mcmc_log_file.txt after an earlier run or by eye
    estimatedModeParameterSet = [-17.95284098, 17.40621784, 7.21967115]

    # run a new simulation with a new PE_object that generates points around the MAP. 
    # This assumes that the MAP exists in the highest posterior density region
    # We take the array directly output from the getPointsNearExistingSample function and put this in initial guess
    path_to_previous_run_log = 'logs_and_csvs/mcmc_logP_and_parameter_samples.csv' # User specified path for logs with the log file name
    # getPointsNearExistingSample returns multiple outputs, we only care about the first, which are the starting point parameter sets
    extracted_parameter_samples, extracted_logP_values, extracted_objective_values = PEUQSE.getPointsNearExistingSample(numPointsToGet=numWalkers, existingSamples=path_to_previous_run_log, parameters_values=estimatedModeParameterSet)
    # PEUQSE creates a pickled object from the getPointsNearExistingSample, but for this case we will pickle the output from the function
    PEUQSE.pickleAnObject(extracted_parameter_samples, 'extracted_parameter_samples.pkl')
    UserInput.model['InputParameterInitialGuess'] = 'extracted_parameter_samples.pkl' # input the starting points as initial guess, this can handle a list of arrays for multiple walkers but only as a pickle

    # create new PE_object
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    # Run the program with EJS
    PE_object.doEnsembleJumpSampling()
    
    PE_object.createAllPlots()
