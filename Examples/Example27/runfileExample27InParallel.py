import sys; sys.path.insert(0, '../../');
import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput
import numpy as np
from multiprocessing import Pool
import os

#%% functions
def run_simulation(inputs):
    
    folder_name, priors = inputs

    import observed_values_00  #Just a simple example. The user can also put the values in directly into the runfile or extract from a csv, for example.
    import simulation_model_00 #Simple example.
    
    try:
        os.mkdir(folder_name)
    except OSError:
        print('')

    os.chdir(folder_name)

    UserInput.responses['responses_abscissa'] = observed_values_00.observed_data_x_values
    UserInput.responses['responses_observed'] = observed_values_00.observed_data_y_values
    UserInput.responses['responses_observed_uncertainties'] = observed_values_00.observed_data_y_values_uncertainties

    
    UserInput.simulated_response_plot_settings['x_label'] = 'distance (m)'
    UserInput.simulated_response_plot_settings['y_label'] = r'$time (s)$'
    UserInput.simulated_response_plot_settings['fontdict'] = {'size':16}
    
    
    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'a':'a','b':'b'}
    UserInput.model['InputParameterPriorValues'] = priors #prior expected values for a and b
    UserInput.model['InputParametersPriorValuesUncertainties'] = [100, 200] #required. #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D
    #UserInput.model['InputParameterInitialGuess'] = [150,400] #Can optionally change the initial guess to be different from prior means.

    
    UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This must simulate with *only* the parameters listed above, and no other arguments.

    UserInput.parameter_estimation_settings['mcmc_length'] = 10000 #10000 is the default.


    #UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 This can be useful for testing.
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    PE_object = PEUQSE.parameter_estimation(UserInput)
    PE_object.doMetropolisHastings()
    PE_object.createAllPlots() #This function calls each of the below functions so that the user does not have to.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlots()


#%% main
if __name__ == "__main__":
    initial_priors = [200, 500] # initial priors for parameters a and b
    p = 0.1 # percent perturbation that will be applied on parameter prior set for multiple mcmc runs
    mesh_set = [[x, x*(1-p), x*(1+p)] for x in initial_priors] 
    perturbed_priors = list(np.array(np.meshgrid(mesh_set[0], mesh_set[1])).T.reshape(-1, 2)) # creating a grid of perturbed sets where n**3 sets are created with n being the number of parameters
    perturbed_priors_tuples = [(f'run_{x}', y) for x, y in enumerate(perturbed_priors)] # labels the parameter sets to seperate runs into seperate folders
#%% pool submission
    pool = Pool()
    pool.map(run_simulation, perturbed_priors_tuples) # run all processes simultaniously 