import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.integrate import odeint
#import dill
import sys; sys.path.insert(0, '../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    
    #observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    
    import simulation_functions_Langmuir_CO_H2O_four_parameters as simulation_functions

    UserInput.responses['responses_abscissa'] = [simulation_functions.T-25,simulation_functions.T,simulation_functions.T+25]
    UserInput.responses['responses_observed'] = [0.00, 0.00, 0.00] #The initial values won't be used during the DOE.
    UserInput.responses['responses_observed_uncertainties'] = [np.log10(1.5),np.log10(1.5),np.log10(1.5)]
    
    UserInput.responses['independent_variables_values'] = [400 , 0.1]
    UserInput.responses['independent_variables_names'] = ['T(K)', 'P_A(bar)']
    simulation_functions.connected_variables_values = UserInput.responses['independent_variables_values'] #It is important to push the list *into* the other module.
    
    
    UserInput.simulated_response_plot_settings['x_label'] = r'$Temperature (K)$'
    UserInput.simulated_response_plot_settings['y_label'] = r'$\theta_A$'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example12' #This creates the filename, also.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'delta_H_rxn':'delta_H_rxn', 'delta_S_rxn':'delta_S_rxn'}
    #Below uses numbers that Eric converted for eV units.
    delta_H_rxn = -0.687 - (-0.858)
    delta_S_rxn = -1.50e-3 - (-1.45E-3)
    delta_H_rxn_uncertainty = (5.00e-2**2 + 1.25e-1**2)**0.5
    delta_S_rxn_uncertainty = (2.07e-4**2 + 2.07e-4**2)**0.5
    UserInput.model['InputParameterPriorValues'] = [delta_H_rxn, delta_S_rxn]                                                                               
    UserInput.model['InputParametersPriorValuesUncertainties'] = [delta_H_rxn_uncertainty, delta_S_rxn_uncertainty] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D  This is a standard deviation!
#    UserInput.model['InputParameterInitialGuess'] = [-0.687, 1.50e-3, -0.858, 1.50e-3] #This is where the mcmc chain will start.
    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    #InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
    UserInput.model['populateIndependentVariablesFunction'] = simulation_functions.populate_pA_and_T #This is needed for design of experiments.                          
    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_functions.Langmuir_replacement_three_temperatures_log #This must simulate with *only* the parameters listed above, and no other arguments.
    #UserInput.model['simulationOutputProcessingFunction'] = processing_functions_tpd_odeint.no_log_wrapper_func #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "std"
    UserInput.parameter_estimation_settings['exportLog'] = False
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] = None
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 100
    UserInput.parameter_estimation_settings['mcmc_length'] = 200
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.05
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] = 'KL_divergence' #obtains the information gain using the Kullback-Leibler divergence    
   
    UserInput.contour_plot_settings['parameter_pairs'] = [[0, 1]]
    UserInput.contour_plot_settings['contours_normalized'] = False
    UserInput.contour_plot_settings['figure_name'] = 'Mumpce_contour_plot_Langmuir_compete_ads'
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    
    
    #It's good to run a test before doing a design of experiments.
#    PE_object = PEUQSE.parameter_estimation(UserInput)
#    PE_object.doMetropolisHastings()
    #PE_object.createAllPlots()
    
    
    UserInput.doe_settings['info_gains_matrices_array_format'] = 'meshgrid'
    UserInput.doe_settings['info_gains_matrices_multiple_parameters'] = 'sum'
    UserInput.doe_settings['independent_variable_grid_center'] = [500, 0.5]
    UserInput.doe_settings['independent_variable_grid_interval_size'] = [100, 0.1]
    UserInput.doe_settings['independent_variable_grid_num_intervals'] = [2,2] #This is the number in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the parameter to vary.
    
    UserInput.doe_settings['parameter_modulation_grid_interval_size'] = [1,1] #use a non-zero value even for parameters that you will not vary.
    UserInput.doe_settings['parameter_modulation_grid_num_intervals'] = [1,1] #make the number of intervals zero for any parameter that you don't want to vary.
    UserInput.doe_settings['parallel_conditions_exploration'] = True
    
    PE_object = PEUQSE.parameter_estimation(UserInput)
    
    
    
    
    PE_object.doeParameterModulationPermutationsScanner()
    #print(PE_object.info_gains_matrices_array[0])
    PE_object.createInfoGainPlots()
    
    
    #To obtain a single info gain matrix, for a single set of indepependet variables, we would use the following syntax:
    # del PE_object
    # UserInput.doe_settings['info_gains_matrices_array_format'] = 'meshgrid'
    # #We *still* have to define an independent variable grid.
    # UserInput.doe_settings['independent_variable_grid_center'] = [500, 0.5]
    # UserInput.doe_settings['independent_variable_grid_interval_size'] = [100, 0.1]
    # UserInput.doe_settings['independent_variable_grid_num_intervals'] = [2,2] #This is the number in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the parameter to vary.
    # #Note that we *no longer* define intervals for the parametric space.
    # simulation_functions.connected_variables_values = UserInput.responses['independent_variables_values'] #It is important to push the list *into* the other module.
    # PE_object2 = PEUQSE.parameter_estimation(UserInput)    
    # PE_object2.doeGetInfoGainMatrix(UserInput.model['InputParameterPriorValues']+UserInput.model['InputParametersPriorValuesUncertainties']) #This is an example with a +1SD perturbation.
    # PE_object2.createInfoGainPlots(plot_suffix="manual")