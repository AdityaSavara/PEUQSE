import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd
import sys

import copy
#import mumce_py.Project as mumce_pyProject #FIXME: Eric to fix plotting/graphing issue described in issue 9 -- https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/9
#import mumce_py.solution mumce_pySolution

class parameter_estimation:
    #Ip for 'inverse problem'. Initialize prior chain starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
    def __init__(self, UserInput = None):
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.parameterNamesAndMathTypeExpressionsDict.keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        self.UserInput = UserInput #Note that this is a pointer, so the last two lines effects are within this object.
        if self.UserInput.verbose: 
            print("Bayes Model Initialized")
        #self.mcmc_length = UserInput.mcmc_length
        #self.mcmc_burn_in = UserInput.mcmc_burn_in # Number of samples trimmed off the beginning of the Markov chain.
        #self.mu_prior = UserInput.mu_prior
        #self.cov_prior = UserInput.cov_prior
        self.Q_mu = self.UserInput.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_cov = self.UserInput.cov_prior/10 # Take small steps. <-- looks like this 20 should be a user defined variable.
#        self.initial_concentrations_array = UserInput.initial_concentrations_array
        #self.modulate_accept_probability = UserInput.modulate_accept_probability
        #self.UserInput.import_experimental_settings(UserInput.Filename) #FIXME: This needs to get out of this function.

    #main function to get samples
    def doMetropolisHastings(self):
        #TODO: add a checkpoint counter or progress bar of some kind.
        if hasattr(self.UserInput, "mcmc_random_seed"):
            if type(self.UserInput.mcmc_random_seed) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.mcmc_random_seed)
        samples_simulatedOutputs = np.zeros((self.UserInput.mcmc_length,self.UserInput.num_data_points)) #TODO: Consider moving this out of this function.
        samples = np.zeros((self.UserInput.mcmc_length,len(self.UserInput.mu_prior)))
        samples[0,:] = self.UserInput.mu_prior # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
        likelihoods_vec = np.zeros((self.UserInput.mcmc_length,1))
        posteriors_un_normed_vec = np.zeros((self.UserInput.mcmc_length,1))
        priors_vec = np.zeros((self.UserInput.mcmc_length,1))
        for i in range(1,self.UserInput.mcmc_length):
            if self.UserInput.verbose: print("MCMC sample number", i)
            if type(self.UserInput.checkPointFrequency) != None:
                if i%self.UserInput.checkPointFrequency == 0: print("MCMC sample number", i) #The % is a modulus function.
            proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov)
            prior_proposal = self.getPrior(proposal_sample)
            [likelihood_proposal, simulationOutput_proposal] = self.getLikelihood(proposal_sample)
            prior_current_location = self.getPrior(samples[i-1,:]) 
            [likelihood_current_location, simulationOutput_current_location] = self.getLikelihood(samples[i-1,:]) #FIXME: the previous likelihood should be stored so that it doesn't need to be calculated again.
            accept_probability = (likelihood_proposal*prior_proposal)/(likelihood_current_location*prior_current_location) 
            if self.UserInput.verbose: print('Current likelihood',likelihood_current_location, 'Proposed likelihood', likelihood_proposal, '\nAccept_probability (gauranteed if above 1)', accept_probability)
            if self.UserInput.verbose: print('Current posterior',likelihood_current_location*prior_current_location, 'Proposed Posterior', likelihood_proposal*prior_proposal)
            if self.UserInput.modulate_accept_probability != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.UserInput.modulate_accept_probability)
                accept_probability = accept_probability**(1/N_flatten) #TODO: add code that unflattens the final histograms, that way even with more sampling we still get an accurate final posterior distribution. We can also then add a flag if the person wants to keep the posterior flattened.
            if accept_probability > np.random.uniform():  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                if self.UserInput.verbose:
                  print('accept', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_proposal)
                samples[i,:] = proposal_sample
                samples_simulatedOutputs[i,:] = simulationOutput_proposal
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal #FIXME: Separate this block of code out into a helper function in the class, that way I can create another helper function for non-MCMC sampling.
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                if self.UserInput.verbose:
                  print('reject', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_current_location)
                samples[i,:] = samples[i-1,:]
                samples_simulatedOutputs[i,:] = simulationOutput_current_location
#                print("line 121", simulationOutput_current_location)
                posteriors_un_normed_vec[i] = likelihood_current_location*prior_current_location
                likelihoods_vec[i] = likelihood_current_location
                priors_vec[i] = prior_current_location
            ########################################
        self.burn_in_samples = samples[:self.UserInput.mcmc_burn_in]
        self.post_burn_in_samples = samples[self.UserInput.mcmc_burn_in:]
        if self.UserInput.exportAllSimulatedOutputs == True:
            self.post_burn_in_samples_simulatedOutputs = samples_simulatedOutputs[self.UserInput.mcmc_burn_in:]
        self.post_burn_in_posteriors_un_normed_vec = posteriors_un_normed_vec[self.UserInput.mcmc_burn_in:]
        self.post_burn_in_logP_un_normed_vec = np.log(self.post_burn_in_posteriors_un_normed_vec)
        self.post_burn_in_likelihoods_vec = likelihoods_vec[self.UserInput.mcmc_burn_in:]
        self.post_burn_in_priors_vec = priors_vec[self.UserInput.mcmc_burn_in:]
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
        if self.UserInput.verbose == True:
            print(self.map_parameter_set)
            print(self.muap_parameter_set)
            print(self.stdap_parameter_set)
        if self.UserInput.exportResults == True:
            #The self.post_burn_in_samples_simulatedOutputs has length of ALL sampling including burn_in
            #The self.post_burn_in_samples is only AFTER burn in.
            #The self.post_burn_in_posteriors_un_normed_vec is AFTER burn in.           
            #TODO: Make header for mcmc_output
            mcmc_output = np.hstack((self.post_burn_in_logP_un_normed_vec,self.post_burn_in_samples))
            np.savetxt('mcmc_logP_and_parameter_samples.csv',mcmc_output, delimiter=",")             
            if self.UserInput.exportAllSimulatedOutputs == True: #By default, we should not keep this, it's a little too large with large sampling.
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
        simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction
        simulationFunctionWrapper =  self.UserInput.simulationFunctionWrapper
        simulationOutput =simulationFunctionWrapper(discreteParameterVector) #FIXME: code should look like simulationOutput = self.UserInput.simulationFunction(*self.UserInput.simulationInputArguments)       
        #simulationOutputProcessingFunction = self.UserInput.log10_wrapper_func #FIXME: this will become fed in as self.UserInput.simulationOutputProcessingFunction
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
#            print("line 168", simulatedResponsesProxy)
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = self.UserInput.simulationOutputProcessingFunction(simulationOutput) #This will become simulatedResponses = self.UserInput.simulationOutputProcessingFunction(simulationOutput)
#            print("line 170", simulatedResponsesProxy)
        observedResponses = self.UserInput.observedResponses
        #To find the relevant covariance, we take the errors from the points.
        cov = self.UserInput.observedResponses_uncertainties #FIXME: We should not be doing subset of points like this here. Should happen at user input level.
        probability_metric = multivariate_normal.pdf(x=simulatedResponses,mean=observedResponses,cov=cov)
#        print("line 178", simulatedResponses, simulatedResponsesProxy)
        return probability_metric, simulatedResponses #FIXME: This needs to say probability_metric, simulatedResponses or something like that, but right now the sizes of the arrays do not match.

    def makeHistogramsForEachParameter(self):
        import plotting_functions
        parameterSamples = self.post_burn_in_samples
        parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        plotting_functions.makeHistogramsForEachParameter(parameterSamples,parameterNamesAndMathTypeExpressionsDict)

    def makeSamplingScatterMatrixPlot(self, parameterSamples = [], parameterNamesAndMathTypeExpressionsDict={}, parameterNamesList =[], plot_settings={}):
        if 'dpi' not in  plot_settings:  plot_settings['dpi'] = 220
        if 'figure_name' not in  plot_settings:  plot_settings['figure_name'] = 'scatter_matrix_posterior'
        if parameterSamples  ==[] : parameterSamples = self.post_burn_in_samples
        if parameterNamesAndMathTypeExpressionsDict == {}: parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        if parameterNamesList == []: parameterNamesList = self.UserInput.parameterNamesList        
        posterior_df = pd.DataFrame(parameterSamples,columns=[parameterNamesAndMathTypeExpressionsDict[x] for x in parameterNamesList])
        pd.plotting.scatter_matrix(posterior_df)
        plt.savefig(plot_settings['figure_name'],dpi=plot_settings['dpi'])
        
    def createSimulatedResponsesPlot(self, x_values=[], listOfYArrays=[], plot_settings={}):            
        if x_values == []: x_values = self.UserInput.responses_abscissa       
        if listOfYArrays ==[]:
            self.map_SimulatedOutput = self.UserInput.simulationFunctionWrapper(self.map_parameter_set)
            self.map_SimulatedResponses = self.UserInput.simulatedOutputProcessingFunction(self.map_SimulatedOutput) #for this line, always no log wrapper because want actual response and not proxy.            
            self.muap_SimulatedOutput = self.UserInput.simulationFunctionWrapper(self.muap_parameter_set)
            self.muap_SimulatedResponses = self.UserInput.simulatedOutputProcessingFunction(self.muap_SimulatedOutput) #for this line, always no log wrapper because want actual response and not proxy.
            listOfYArrays = [self.UserInput.responses_observed,self.map_SimulatedResponses, self.muap_SimulatedResponses]        
        if plot_settings == {}: plot_settings = self.UserInput.simulated_response_plot_settings
        import plotting_functions
        figureObject = plotting_functions.createSimulatedResponsesPlot(x_values, listOfYArrays, plot_settings)
        return figureObject #This figure is a matplotlib.pyplot as plt object.

    def createAllPlots(self):
        makeHistogramsForEachParameter()    
        makeSamplingScatterMatrixPlot()
        createSimulatedResponsesPlot()
        
if __name__ == "__main__":
    pass
