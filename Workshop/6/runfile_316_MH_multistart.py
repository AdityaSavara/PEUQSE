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
    
    # Multistart checkpoints set 
    UserInput.parameter_estimation_settings['multistart_checkPointFrequency'] = 1
    
    # Set walkers 
    UserInput.parameter_estimation_settings['multistart_numStartPoints'] = 12
    
    # Do multistart
    UserInput.parameter_estimation_settings['multistart_searchType'] = 'doMetropolisHastings'
    
    # Optional bound setting lines for finding uninformed parameters.
    UserInput.model['InputParameterPriorValues_upperBounds'] = [np.log10(1.0E-17), 18.0, 18.0] 
    UserInput.model['InputParameterPriorValues_lowerBounds'] = [np.log10(1.0E-19), 16.0, 0.0]

    # Provides simulation function for Cr concentration throughout a sample
    UserInput.model['simulateByInputParametersOnlyFunction'] = function_316.simulation_function_using_log_a_wrapper
    
    # UserInput.model['walkerInitialDistributionSpread'] = 0.25 # [Optional] line to reduce initial distribution if neccesary

    # Enable checkpoints
    # UserInput.parameter_estimation_settings['checkPointFrequency'] = 1

    # Reduced sample size needed for EnsembleSliceSampling() due to single-mode data
    UserInput.parameter_estimation_settings['mcmc_length'] = 100000 # 10000 is the default.
    
    # UserInput.parameter_estimation_settings['mcmc_threshold_filter_coefficient'] = 2.0 
       
    # UserInput.parameter_estimation_settings['mcmc_walkerInitialDistribution'] = 'identical'
    # After filinlg the variables of the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    # Run the program with Multistart MH
    PE_object.doMultiStart()
    
    PE_object.createAllPlots()