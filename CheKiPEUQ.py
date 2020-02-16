import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd
import sys
import timeit
import copy
#import mumce_py.Project as mumce_pyProject #FIXME: Eric to fix plotting/graphing issue described in issue 9 -- https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/9
#import mumce_py.solution mumce_pySolution

class parameter_estimation:
    #Ip for 'inverse problem'. Initialize prior chain starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
    def __init__(self, UserInput = None):
        self.UserInput = UserInput #Note that this is a pointer, so the later lines are within this object.
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.model['parameterNamesAndMathTypeExpressionsDict'].keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        if self.UserInput.parameter_estimation_settings['verbose']: 
            print("Bayes Model Initialized")
        #Leaving the original dictionary object intact, but making a new object to make cov_prior.
        UserInput.InputParametersInitialValuesUncertainties = UserInput.model['InputParametersInitialValuesUncertainties']     
        if len(np.shape(UserInput.InputParametersInitialValuesUncertainties)) == 1 and (len(UserInput.InputParametersInitialValuesUncertainties) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
            UserInput.cov_prior = np.diagflat(UserInput.InputParametersInitialValuesUncertainties) 
        elif len(np.shape(UserInput.InputParametersInitialValuesUncertainties)) > 1: #If it's non-1D, we assume it's already a covariance matrix.
            UserInput.cov_prior = UserInput.InputParametersInitialValuesUncertainties
        else: #If a blank list is received, that means the user
            print("The covariance of the priors is undefined because InputParametersInitialValuesUncertainties is blank.")
        #    cov_prior = np.array([[200.0, 0., 0., 0., 0., 0.], 
        #                          [0., 200.0, 0., 0., 0., 0.],
        #                          [0., 0., 13.0, 0., 0., 0.],
        #                          [0., 0., 0., 13.0, 0., 0.],
        #                          [0., 0., 0., 0., 0.1, 0.],
        #                          [0., 0., 0., 0., 0., 0.1]])
#

        self.UserInput.num_data_points = len(UserInput.responses['responses_abscissa'])
        self.UserInput.mu_prior = np.array(UserInput.model['InputParameterInitialValues']) 
        #self.cov_prior = UserInput.cov_prior
        self.Q_mu = self.UserInput.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_cov = self.UserInput.cov_prior/10 # Take small steps. <-- looks like this 20 should be a user defined variable.
#        self.initial_concentrations_array = UserInput.initial_concentrations_array
        #self.modulate_accept_probability = UserInput.modulate_accept_probability
        #self.UserInput.import_experimental_settings(UserInput.Filename) #FIXME: This needs to get out of this function.
    #main function to get samples
    def doMetropolisHastings(self):
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if type(self.UserInput.parameter_estimation_settings['mcmc_random_seed']) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        samples_simulatedOutputs = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],self.UserInput.num_data_points)) #TODO: Consider moving this out of this function.
        samples = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],len(self.UserInput.mu_prior)))
        samples[0,:] = self.UserInput.mu_prior # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
        likelihoods_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        posteriors_un_normed_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        priors_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        #Code to initialize checkpoints.
        if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
            print("Starting MCMC sampling.")
            import timeit
            timeOfFirstCheckpoint = timeit.time.clock()
            timeCheckpoint = timeOfFirstCheckpoint #First checkpoint at time 0.
            numCheckPoints = self.UserInput.parameter_estimation_settings['mcmc_length']/self.UserInput.parameter_estimation_settings['checkPointFrequency']
        for i in range(1,self.UserInput.parameter_estimation_settings['mcmc_length']):
            if self.UserInput.parameter_estimation_settings['verbose']: print("MCMC sample number", i)                  
            proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov)
            prior_proposal = self.getPrior(proposal_sample)
            [likelihood_proposal, simulationOutput_proposal] = self.getLikelihood(proposal_sample)
            prior_current_location = self.getPrior(samples[i-1,:]) 
            [likelihood_current_location, simulationOutput_current_location] = self.getLikelihood(samples[i-1,:]) #FIXME: the previous likelihood should be stored so that it doesn't need to be calculated again.
            accept_probability = (likelihood_proposal*prior_proposal)/(likelihood_current_location*prior_current_location) 
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current likelihood',likelihood_current_location, 'Proposed likelihood', likelihood_proposal, '\nAccept_probability (gauranteed if above 1)', accept_probability)
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current posterior',likelihood_current_location*prior_current_location, 'Proposed Posterior', likelihood_proposal*prior_proposal)
            if self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'] != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'])
                accept_probability = accept_probability**(1/N_flatten) #TODO: add code that unflattens the final histograms, that way even with more sampling we still get an accurate final posterior distribution. We can also then add a flag if the person wants to keep the posterior flattened.
            if accept_probability > np.random.uniform():  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                if self.UserInput.parameter_estimation_settings['verbose']:
                  print('accept', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_proposal)
                samples[i,:] = proposal_sample
                samples_simulatedOutputs[i,:] = simulationOutput_proposal
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal #FIXME: Separate this block of code out into a helper function in the class, that way I can create another helper function for non-MCMC sampling.
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                if self.UserInput.parameter_estimation_settings['verbose']:
                  print('reject', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_current_location)
                samples[i,:] = samples[i-1,:]
                samples_simulatedOutputs[i,:] = simulationOutput_current_location
#                print("line 121", simulationOutput_current_location)
                posteriors_un_normed_vec[i] = likelihood_current_location*prior_current_location
                likelihoods_vec[i] = likelihood_current_location
                priors_vec[i] = prior_current_location
            if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
                if i%self.UserInput.parameter_estimation_settings['checkPointFrequency'] == 0: #The % is a modulus function.
                    timeSinceLastCheckPoint = timeit.time.clock() -  timeCheckpoint
                    timeCheckpoint = timeit.time.clock() - timeOfFirstCheckpoint
                    checkPointNumber = i/self.UserInput.parameter_estimation_settings['checkPointFrequency']
                    averagetimePerSampling = timeCheckpoint/i
                    print("MCMC sample number ", i, "checkpoint", checkPointNumber, "out of", numCheckPoints) 
                    print("averagetimePerSampling", averagetimePerSampling, "seconds")
                    print("timeSinceLastCheckPoint", timeSinceLastCheckPoint, "seconds")
                    print("Estimated time remaining", averagetimePerSampling*(self.UserInput.parameter_estimation_settings['mcmc_length']-i), "seconds")
            ########################################
        self.burn_in_samples = samples[:self.UserInput.parameter_estimation_settings['mcmc_burn_in']]
        self.post_burn_in_samples = samples[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
            self.post_burn_in_samples_simulatedOutputs = samples_simulatedOutputs[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_posteriors_un_normed_vec = posteriors_un_normed_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_logP_un_normed_vec = np.log(self.post_burn_in_posteriors_un_normed_vec)
        self.post_burn_in_likelihoods_vec = likelihoods_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_priors_vec = priors_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        # posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        self.evidence = np.mean(self.post_burn_in_posteriors_un_normed_vec)*np.sqrt(2*np.pi*np.std(self.post_burn_in_samples)**2)
        post_burn_in_posteriors_vec = self.post_burn_in_posteriors_un_normed_vec/self.evidence
        log_ratios = np.log(post_burn_in_posteriors_vec/self.post_burn_in_priors_vec)
        log_ratios[np.isinf(log_ratios)] = 0
        log_ratios = np.nan_to_num(log_ratios)
        self.info_gain = np.mean(log_ratios)
        map_logP = max(self.post_burn_in_logP_un_normed_vec)
        self.map_index = list(self.post_burn_in_logP_un_normed_vec).index(map_logP)
        self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.
        self.muap_parameter_set = np.mean(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and muap will be the same.
        self.stdap_parameter_set = np.std(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and muap will be the same.
        #TODO: should return the variance of each sample in the post_burn_in
        if self.UserInput.parameter_estimation_settings['verbose'] == True:
            print(self.map_parameter_set)
            print(self.muap_parameter_set)
            print(self.stdap_parameter_set)
        if self.UserInput.parameter_estimation_settings['exportLog'] == True:
            #The self.post_burn_in_samples_simulatedOutputs has length of ALL sampling including burn_in
            #The self.post_burn_in_samples is only AFTER burn in.
            #The self.post_burn_in_posteriors_un_normed_vec is AFTER burn in.           
            #TODO: Make header for mcmc_output
            mcmc_output = np.hstack((self.post_burn_in_logP_un_normed_vec,self.post_burn_in_samples))
            np.savetxt('mcmc_logP_and_parameter_samples.csv',mcmc_output, delimiter=",")             
            if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True: #By default, we should not keep this, it's a little too large with large sampling.
                np.savetxt('mcmc_all_simulated_outputs.csv',self.post_burn_in_samples_simulatedOutputs, delimiter=",")             
            with open("mcmc_log_file.txt", 'w') as out_file:
                out_file.write("MAP_logP:" +  str(map_logP) + "\n")
                out_file.write("self.map_index:" +  str(self.map_index) + "\n")
                out_file.write("self.map_parameter_set:" + str( self.map_parameter_set) + "\n")
                out_file.write("self.muap_parameter_set:" + str( self.muap_parameter_set) + "\n")
                out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
                out_file.write("evidence:" + str(self.evidence) + "\n")
        return [self.map_parameter_set, self.muap_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] # EAW 2020/01/08
    def getPrior(self,discreteParameterVector):
        probability = multivariate_normal.pdf(x=discreteParameterVector,mean=self.UserInput.mu_prior,cov=self.UserInput.cov_prior)
        return probability
    def getLikelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        #This is more in accordance with https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/11. 
        simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']
        simulationOutput =simulationFunction(discreteParameterVector) 
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        observedResponses = self.UserInput.responses['responses_observed']
        #To find the relevant covariance, we take the errors from the points.
        responses_cov = self.UserInput.responses['responses_observed_uncertainties'] 
        probability_metric = multivariate_normal.pdf(x=simulatedResponses,mean=observedResponses,cov=responses_cov)
        return probability_metric, simulatedResponses

    def makeHistogramsForEachParameter(self):
        import plotting_functions
        parameterSamples = self.post_burn_in_samples
        parameterNamesAndMathTypeExpressionsDict = self.UserInput.model['parameterNamesAndMathTypeExpressionsDict']
        plotting_functions.makeHistogramsForEachParameter(parameterSamples,parameterNamesAndMathTypeExpressionsDict)

    def makeSamplingScatterMatrixPlot(self, parameterSamples = [], parameterNamesAndMathTypeExpressionsDict={}, parameterNamesList =[], plot_settings={}):
        if 'dpi' not in  plot_settings:  plot_settings['dpi'] = 220
        if 'figure_name' not in  plot_settings:  plot_settings['figure_name'] = 'scatter_matrix_posterior'
        if parameterSamples  ==[] : parameterSamples = self.post_burn_in_samples
        if parameterNamesAndMathTypeExpressionsDict == {}: parameterNamesAndMathTypeExpressionsDict = self.UserInput.model['parameterNamesAndMathTypeExpressionsDict']
        if parameterNamesList == []: parameterNamesList = self.UserInput.parameterNamesList #This is created when the parameter_estimation object is initialized.        
        posterior_df = pd.DataFrame(parameterSamples,columns=[parameterNamesAndMathTypeExpressionsDict[x] for x in parameterNamesList])
        pd.plotting.scatter_matrix(posterior_df)
        plt.savefig(plot_settings['figure_name'],dpi=plot_settings['dpi'])
        
    def createSimulatedResponsesPlot(self, x_values=[], listOfYArrays=[], plot_settings={}):            
        if x_values == []: x_values = self.UserInput.responses['responses_abscissa']       
        if listOfYArrays ==[]:
            #Get map simiulated output and simulated responses.
            self.map_SimulatedOutput = self.UserInput.model['simulateByInputParametersOnlyFunction'](self.map_parameter_set)           
            if type(self.UserInput.model['simulationOutputProcessingFunction']) == type(None):
                self.map_SimulatedResponses = self.map_SimulatedOutput #Is this the log of the rate? If so, Why?
            if type(self.UserInput.model['simulationOutputProcessingFunction']) != type(None):
                self.map_SimulatedResponses =  self.UserInput.model['simulationOutputProcessingFunction'](self.map_SimulatedOutput)                    
            #Get muap simiulated output and simulated responses.
            self.muap_SimulatedOutput = self.UserInput.model['simulateByInputParametersOnlyFunction'](self.muap_parameter_set)
            if type(self.UserInput.model['simulationOutputProcessingFunction']) == type(None):
                self.muap_SimulatedResponses = self.muap_SimulatedOutput #Is this the log of the rate? If so, Why?
            if type(self.UserInput.model['simulationOutputProcessingFunction']) != type(None):
                self.muap_SimulatedResponses =  self.UserInput.model['simulationOutputProcessingFunction'](self.muap_SimulatedOutput) 
            listOfYArrays = [self.UserInput.responses['responses_observed'],self.map_SimulatedResponses, self.muap_SimulatedResponses]        
        if plot_settings == {}: 
            plot_settings = self.UserInput.simulated_response_plot_settings
            plot_settings['legendLabels'] = ['experiments', 'MAP','mu_AP']
            #Other allowed settings are like this, but will be fed in as simulated_response_plot_settings keys rather than plot_settings keys.
            #plot_settings['x_label'] = 'T (K)'
            #plot_settings['y_label'] = r'$rate (s^{-1})$'
            #plot_settings['y_range'] = [0.00, 0.025] #optional.
            #plot_settings['figure_name'] = 'tprposterior'
            
        import plotting_functions
        figureObject = plotting_functions.createSimulatedResponsesPlot(x_values, listOfYArrays, plot_settings)
        return figureObject #This figure is a matplotlib.pyplot as plt object.

    def createAllPlots(self):
        self.makeHistogramsForEachParameter()    
        self.makeSamplingScatterMatrixPlot()
        self.createSimulatedResponsesPlot()
        
if __name__ == "__main__":
    pass
