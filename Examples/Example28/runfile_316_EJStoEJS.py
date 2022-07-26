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
    
    # Guesses are provided, since the posteriors deviate significntly from literature values.
    # UserInput.model['InputParameterInitialGuess'] = [4.2E-19, 17.0, 5.0]

    # Provides simulation function for Cr concentration throughout a sample
    UserInput.model['simulateByInputParametersOnlyFunction'] = function_316.simulation_function_using_log_a_wrapper
    
    # UserInput.model['walkerInitialDistributionSpread'] = 0.25 # [Optional] line to reduce initial distribution if neccesary

    # Enable checkpoints
    # UserInput.parameter_estimation_settings['checkPointFrequency'] = 1

    # Reduced sample size since this mcmc run will essentially be a burn in 
    UserInput.parameter_estimation_settings['mcmc_length'] = 10000 # 10000 is the default.
    
    # UserInput.parameter_estimation_settings['mcmc_threshold_filter_coefficient'] = 2.0 
       
    # After filinlg the variables of the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object_burnin = PEUQSE.parameter_estimation(UserInput)
    
    # Run the program with EJS
    PE_object_burnin.doEnsembleJumpSampling()

    # return the MAP
    MAP_burnin = PE_object_burnin.map_parameter_set

    # run a new simulation with a new PE_object that generates points around the MAP. 
    # This assumes that the MAP exists in the highest posterior density region
    # A pickle file of the starting points are created in the pickles directory as specified
    numWalkers = 24 # run next mcmc run with 24 walkers
    path_to_previous_run_log = UserInput.directories['logs_and_csvs'] + 'mcmc_logP_and_parameter_samples.csv' # concatenate User specified path for logs with the log file name
    pickleFilePath = UserInput.directories['pickles'] + 'pointsNearExistingSample.pkl'
    PEUQSE.getPointsNearExistingSample(numPointsToGet=numWalkers, existingSamples=path_to_previous_run_log, parameters_values=MAP_burnin, pickleFileName=pickleFilePath)

    # set up next simulation
    UserInput.parameter_estimation_settings['mcmc_length'] = 100000
    UserInput.parameter_estimation_settings['mcmc_nwalkers'] = numWalkers
    UserInput.model['InputParameterInitialGuess'] = pickleFilePath # by inputting a string, the input is assumed to be a pickled file

    # create new PE_object **VERY IMPORTANT**
    PE_object = PEUQSE.parameter_estimation(UserInput)

    # run the next EJS sampling
    PE_object.doEnsembleJumpSampling()


    
    
    '''
    UserInput.parameter_estimation_settings['mcmc_threshold_filter_coefficient'] = 0.10 
    
    UserInput.parameter_estimation_settings['mcmc_length'] = 1000000
    UserInput.parameter_estimation_settings['mcmc_continueSampling']
    PE_object.doEnsembleJumpSampling()
    
    '''
    '''
    #use 1,000,000 steps for this second step.
    UserInput.parameter_estimation_settings['mcmc_length'] = 100000
    UserInput.parameter_estimation_settings['mcmc_continueSampling']
    PE_object.doEnsembleSliceSampling(movesType='global') #change sampling type.
    '''
    # PE_object.doMetropolisHastings()    

    # PE_object.doOptimizeSSR()

    # PE_object.doEnsembleSliceSampling()    

    # PE_object.doSinglePoint()
    
    # PE_object.doOptimizeLogP(method="BFGS", printOptimum=True, verbose=True)
    
    # Another option would be PE_object.doEnsembleSliceSampling(), one can also do grid search or an astroidal distribution search.
    
    # Create histograms and reponses only. PE_object.makeSamplingScatterMatrixPlot() does not seem to work
    # PE_object.createSimulatedResponsesPlots()
    # PE_object.makeHistogramsForEachParameter()
    
    PE_object_burnin.createAllPlots()
""" 
    #########Optional example of saving and loading PE_objects after running the mcmc.
    #########This feature requires having dill installed (pip install dill, https://pypi.org/project/dill/)
    try:
        import dill
        dillModuleExists = True
    except:
        dillModuleExists = False

    
    #Optionally, one can save a PE_object for later,if the dill module has been installed.
    if dillModuleExists == True: 
        PE_object.save_to_dill("PE_object_00a0")
        #to load a PE_object after some time, first one has to put (any) UserInput to create a PE_object, then to load from file.
        
        
        #Normally, we would do the loading and plotting in another python file, but for this example the syntax is being demonstrated below within the same file.
        PE_object2 = PEUQSE.parameter_estimation(UserInput)
        PE_object2 = PE_object2.load_from_dill("PE_object_00a0")
        PE_object2.createAllPlots()
"""