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
    #Ip for 'inverse problem'. Initialize chain with initial guess (prior if not provided) as starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  
    def __init__(self, UserInput = None):
        self.UserInput = UserInput #Note that this is a pointer, so the later lines are within this object.
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.model['parameterNamesAndMathTypeExpressionsDict'].keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        if self.UserInput.parameter_estimation_settings['verbose']: 
            print("Bayes Model Initialized")
        #Leaving the original dictionary object intact, but making a new object to make cov_prior.
        UserInput.InputParametersPriorValuesUncertainties = UserInput.model['InputParametersPriorValuesUncertainties']     
        if len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) == 1 and (len(UserInput.InputParametersPriorValuesUncertainties) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
            UserInput.cov_prior = np.diagflat(UserInput.InputParametersPriorValuesUncertainties) 
        elif len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) > 1: #If it's non-1D, we assume it's already a covariance matrix.
            UserInput.cov_prior = UserInput.InputParametersPriorValuesUncertainties
        else: #If a blank list is received, that means the user
            print("The covariance of the priors is undefined because InputParametersPriorValuesUncertainties is blank.")
        #    cov_prior = np.array([[200.0, 0., 0., 0., 0., 0.], 
        #                          [0., 200.0, 0., 0., 0., 0.],
        #                          [0., 0., 13.0, 0., 0., 0.],
        #                          [0., 0., 0., 13.0, 0., 0.],
        #                          [0., 0., 0., 0., 0.1, 0.],
        #                          [0., 0., 0., 0., 0., 0.1]])
#

        self.UserInput.num_data_points = len(UserInput.responses['responses_abscissa'])
        self.UserInput.mu_prior = np.array(UserInput.model['InputParameterPriorValues']) 
        #self.cov_prior = UserInput.cov_prior
        self.Q_mu = self.UserInput.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_cov = self.UserInput.cov_prior # Take small steps. 
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
        mcmc_step_modulation_history = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: Make this optional for efficiency. #This allows the steps to be larger or smaller. Make this same length as samples. In future, should probably be same in other dimension also, but that would require 2D sampling with each step.                                                                          
        if 'InputParameterInitialGuess' not in self.UserInput.model: #if an initial guess is not provided, we use the prior.
            self.UserInput.model['InputParameterInitialGuess'] = self.UserInput.mu_prior
        samples[0,:]=self.UserInput.model['InputParameterInitialGuess']  # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean if an initial guess is not provided..
        samples_drawn = samples*1.0 #this includes points that were rejected. #TODO: make this optional for efficiency.               
        likelihoods_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        posteriors_un_normed_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        log_postereriors_drawn = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: make this optional for efficiency. We don't want this to be 2D, so we don't copy posteriors_un_normed_vec.
        priors_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        #Code to initialize checkpoints.
        if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
            print("Starting MCMC sampling.")
            import timeit
            timeOfFirstCheckpoint = timeit.time.clock()
            timeCheckpoint = timeit.time.clock() - timeOfFirstCheckpoint #First checkpoint at time 0.
            numCheckPoints = self.UserInput.parameter_estimation_settings['mcmc_length']/self.UserInput.parameter_estimation_settings['checkPointFrequency']
        for i in range(1,self.UserInput.parameter_estimation_settings['mcmc_length']): #FIXME: Don't we need to start with i of 0?
            if self.UserInput.parameter_estimation_settings['verbose']: print("MCMC sample number", i)                  
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'unbiased':
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'MAP_finding':
                if i == 1: mcmc_step_dynamic_coefficient = 1
                mcmc_step_modulation_coefficient = np.random.uniform() + 0.5 #TODO: make this a 2D array. One for each parameter.
                mcmc_step_modulation_history[i] = mcmc_step_modulation_coefficient
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov*mcmc_step_dynamic_coefficient*mcmc_step_modulation_coefficient*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
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
                samples_drawn[i,:] = proposal_sample
                log_postereriors_drawn[i] = np.log(likelihood_proposal*prior_proposal)
                samples_simulatedOutputs[i,:] = simulationOutput_proposal
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal #FIXME: Separate this block of code out into a helper function in the class, that way I can create another helper function for non-MCMC sampling.
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                if self.UserInput.parameter_estimation_settings['verbose']:
                  print('reject', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_current_location)
                samples[i,:] = samples[i-1,:] #the sample is not kept if it is rejected, though we still store it in the samples_drawn.
                samples_drawn[i,:] = proposal_sample
                log_postereriors_drawn[i] = np.log(likelihood_proposal*prior_proposal)
                samples_simulatedOutputs[i,:] = simulationOutput_current_location
#                print("line 121", simulationOutput_current_location)
                posteriors_un_normed_vec[i] = likelihood_current_location*prior_current_location
                likelihoods_vec[i] = likelihood_current_location
                priors_vec[i] = prior_current_location
            if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
                if i%self.UserInput.parameter_estimation_settings['checkPointFrequency'] == 0: #The % is a modulus function.
                    timeSinceLastCheckPoint = (timeit.time.clock() - timeOfFirstCheckpoint) -  timeCheckpoint
                    timeCheckpoint = timeit.time.clock() - timeOfFirstCheckpoint
                    checkPointNumber = i/self.UserInput.parameter_estimation_settings['checkPointFrequency']
                    averagetimePerSampling = timeCheckpoint/i
                    print("MCMC sample number ", i, "checkpoint", checkPointNumber, "out of", numCheckPoints) 
                    print("averagetimePerSampling", averagetimePerSampling, "seconds")
                    print("timeSinceLastCheckPoint", timeSinceLastCheckPoint, "seconds")
                    print("Estimated time remaining", averagetimePerSampling*(self.UserInput.parameter_estimation_settings['mcmc_length']-i), "seconds")
                    if self.UserInput.parameter_estimation_settings['mcmc_mode'] != 'unbiased':
                        print("Most recent mcmc_step_dynamic_coefficient:", mcmc_step_dynamic_coefficient)
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] != 'unbiased':
                if i%100== 0: #The % is a modulus function to change the modulation coefficient every n steps.
                    if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'MAP_finding':
                        recent_log_postereriors_drawn=log_postereriors_drawn[i-100:i] 
                        recent_mcmc_step_modulation_history=mcmc_step_modulation_history[i-100:i]
                        #Make a 2D array and remove anything that is not finite.
                        #let's find out where the posterior is not finite:
                        recent_log_postereriors_drawn_is_finite = np.isfinite(recent_log_postereriors_drawn) #gives 1 if is finite, 0 if not.
                        #Now let's find the cases that were not...
                        not_finite_indices = np.where(recent_log_postereriors_drawn_is_finite == 0)
                        #Now delete the indices we don't want.
                        recent_log_postereriors_drawn = np.delete(recent_log_postereriors_drawn, not_finite_indices)
                        recent_mcmc_step_modulation_history = np.delete(recent_mcmc_step_modulation_history, not_finite_indices)
#                        recent_stacked = np.vstack((recent_log_postereriors_drawn,recent_mcmc_step_modulation_history)).transpose()                                              
#                        print(recent_stacked)
#                        np.savetxt("recent_stacked.csv",recent_stacked, delimiter=',')
                        #Numpy polyfit uses "x, y, degree" for nomenclature. We want posterior as function of modulation history.
                        linearFit = np.polynomial.polynomial.polyfit(recent_mcmc_step_modulation_history, recent_log_postereriors_drawn, 1) #In future, use multidimensional and numpy.gradient or something like that? 
                        #The slope is in the 2nd index of linearFit, despite what the documentation says.
                        #A positive slope means that bigger steps have better outcomes, on average.
                        if linearFit[1] > 0:
                            if mcmc_step_dynamic_coefficient < 10:
                                mcmc_step_dynamic_coefficient = mcmc_step_dynamic_coefficient*1.05
                        if linearFit[1] < 0:
                            if mcmc_step_dynamic_coefficient > 0.1:
                                mcmc_step_dynamic_coefficient = mcmc_step_dynamic_coefficient*0.95
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
        self.mu_AP_parameter_set = np.mean(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.
        self.stdap_parameter_set = np.std(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.
        #TODO: should return the variance of each sample in the post_burn_in
        if self.UserInput.parameter_estimation_settings['verbose'] == True:
            print(self.map_parameter_set)
            print(self.mu_AP_parameter_set)
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
                out_file.write("self.mu_AP_parameter_set:" + str( self.mu_AP_parameter_set) + "\n")
                out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
                out_file.write("evidence:" + str(self.evidence) + "\n")
                if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.mu_AP_parameter_set).any() > 0.10:
                    out_file.write("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632")
        if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.mu_AP_parameter_set).any() > 0.10:
            print("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632  ")
        return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] # EAW 2020/01/08
    def getPrior(self,discreteParameterVector):
        probability = multivariate_normal.pdf(x=discreteParameterVector,mean=self.UserInput.mu_prior,cov=self.UserInput.cov_prior)
        return probability
    def getLikelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
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
        #If our likelihood is  “probability of Response given Theta”…  we have a continuous probability distribution for both the response and theta. That means the pdf  must use binning on both variables. Eric notes that the pdf returns a probability density, not a probability mass. So the pdf function here divides by the width of whatever small bin is being used and then returns the density accordingly. Because of this, our what we are calling likelihood is not actually probability (it’s not the actual likelihood) but is proportional to the likelihood.
        #This we call it a probability_metric and not a probability. #TODO: consider changing likelihood and get likelihood to "likelihoodMetric" and "getLikelihoodMetric"
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
            #Get mu_AP simiulated output and simulated responses.
            self.mu_AP_SimulatedOutput = self.UserInput.model['simulateByInputParametersOnlyFunction'](self.mu_AP_parameter_set)
            if type(self.UserInput.model['simulationOutputProcessingFunction']) == type(None):
                self.mu_AP_SimulatedResponses = self.mu_AP_SimulatedOutput #Is this the log of the rate? If so, Why?
            if type(self.UserInput.model['simulationOutputProcessingFunction']) != type(None):
                self.mu_AP_SimulatedResponses =  self.UserInput.model['simulationOutputProcessingFunction'](self.mu_AP_SimulatedOutput) 
            #Get mu_guess simulated output and responses. 
            self.mu_guess_SimulatedOutput = self.UserInput.model['simulateByInputParametersOnlyFunction']( self.UserInput.model['InputParameterInitialGuess'])
            if type(self.UserInput.model['simulationOutputProcessingFunction']) == type(None):
                self.mu_guess_SimulatedResponses = self.mu_guess_SimulatedOutput #Is this the log of the rate? If so, Why?
            if type(self.UserInput.model['simulationOutputProcessingFunction']) != type(None):
                self.mu_guess_SimulatedResponses =  self.UserInput.model['simulationOutputProcessingFunction'](self.mu_guess_SimulatedOutput)                               
            listOfYArrays = [self.UserInput.responses['responses_observed'],self.map_SimulatedResponses, self.mu_AP_SimulatedResponses, self.mu_guess_SimulatedResponses]        
        if plot_settings == {}: 
            plot_settings = self.UserInput.simulated_response_plot_settings
            plot_settings['legendLabels'] = ['experiments', 'MAP','mu_posterior', 'mu_guess']
            #Other allowed settings are like this, but will be fed in as simulated_response_plot_settings keys rather than plot_settings keys.
            #plot_settings['x_label'] = 'T (K)'
            #plot_settings['y_label'] = r'$rate (s^{-1})$'
            #plot_settings['y_range'] = [0.00, 0.025] #optional.
            #plot_settings['figure_name'] = 'tprposterior'
            
        import plotting_functions
        figureObject = plotting_functions.createSimulatedResponsesPlot(x_values, listOfYArrays, plot_settings)
        return figureObject #This figure is a matplotlib.pyplot as plt object.

    def createMumpcePlots(self):
        import plotting_functions
        from plotting_functions import plotting_functions_class
        figureObject_beta = plotting_functions_class(self.UserInput) # The "beta" is only to prevent namespace conflicts with 'figureObject'.
        parameterSamples = self.post_burn_in_samples
        
        #TODO: the posterior mu_vector and cov_matrix should be calculated elsewhere.
        posterior_mu_vector = np.mean(parameterSamples, axis=0)
        posterior_cov_matrix = np.cov(self.post_burn_in_samples.T)
        #TODO: In future, worry about whether there are constants or not, since then we will have to trim down the prior.
        #Make the model_parameter_info object that mumpce Project class needs.
        self.UserInput.model_parameter_info = []#This variable name is for mumpce definition of variable names. Not what we would choose otherwise.
        for parameterIndex, parameterName in enumerate(self.UserInput.model['parameterNamesAndMathTypeExpressionsDict']):
            individual_model_parameter_dictionary = {'parameter_number': parameterIndex, 'parameter_name': self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'][parameterName]} #we are actually putting the MathTypeExpression as the parameter name when feeding to mum_pce.
            self.UserInput.model_parameter_info.append(individual_model_parameter_dictionary)
        self.UserInput.model_parameter_info = np.array(self.UserInput.model_parameter_info)
        numParams = len(self.UserInput.model_parameter_info)
        active_parameters = np.linspace(0, numParams-1, numParams) #just a list of whole numbers.
        active_parameters = np.array(active_parameters, dtype='int')
        #TODO: reduce active_parameters by anything that has been set as a constant.
        pairs_of_parameter_indices = self.UserInput.parameter_pairs_for_contour_plots
        if pairs_of_parameter_indices == []:
            import itertools 
            all_pairs_iter = itertools.combinations(active_parameters, 2)
            all_pairs_list = list(all_pairs_iter)
            pairs_of_parameter_indices = all_pairs_list #right now these are tuples, and we need lists inside.
            for  pairIndex in range(len(pairs_of_parameter_indices)):
                pairs_of_parameter_indices[pairIndex] = list(pairs_of_parameter_indices[pairIndex])
        elif type(pairs_of_parameter_indices[0]) == type('string'):
            pairs_of_parameter_indices = self.UserInput.pairs_of_parameter_indices
            for  pairIndex in range(len(pairs_of_parameter_indices)):
                firstParameter = int(self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'][pairIndex[0]])
                secondParameter = int(self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'][pairIndex[0]])
                pairs_of_parameter_indices[pairIndex] = [firstParameter, secondParameter]        
        figureObject_beta.mumpce_plots(model_parameter_info = self.UserInput.model_parameter_info, active_parameters = self.UserInput.active_parameters, pairs_of_parameter_indices = pairs_of_parameter_indices, posterior_mu_vector = posterior_mu_vector, posterior_cov_matrix = posterior_cov_matrix, prior_mu_vector = np.array(self.UserInput.model['InputParameterInitialValues']), prior_cov_matrix = self.UserInput.cov_prior, contour_settings_custom = self.UserInput.contour_settings_custom)
        return figureObject_beta

    def createAllPlots(self):
        self.makeHistogramsForEachParameter()    
        self.makeSamplingScatterMatrixPlot()
        self.createSimulatedResponsesPlot()
        self.createMumpcePlots()

'''Below are a bunch of functions for Euler's Method.'''
#This takes an array of dydt values. #Note this is a local dydtArray, it is NOT a local deltaYArray.
def littleEulerGivenArray(y_initial, t_values, dydtArray): 
    #numPoints = len(t_values)
    simulated_t_values = t_values #we'll simulate at the t_values given.
    simulated_y_values = np.zeros(len(simulated_t_values)) #just initializing.
    simulated_y_values[0] = y_initial
    dydt_values = dydtArray #We already have them, just need to calculate the delta_y values.
    for y_index in range(len(simulated_y_values)-1):
        localSlope = dydtArray[y_index]
        deltat_resolution = t_values[y_index+1]-t_values[y_index]
        simulated_y_values[y_index+1] = simulated_y_values[y_index] + localSlope * deltat_resolution
#        print(simulated_t_values[y_index+1], simulated_y_values[y_index+1], localSlope, localSlope * deltat_resolution)
#        print(simulated_y_values[y_index], simulated_t_values[y_index]*10-(simulated_t_values[y_index]**2)/2 +2)
    return simulated_t_values, simulated_y_values, dydt_values

#The initial_y_uncertainty is a scalar, the dydt_uncertainties is an array. t_values is an arrray, so the npoints don't need to be evenly spaced.
def littleEulerUncertaintyPropagation(dydt_uncertainties, t_values, initial_y_uncertainty=0):
    y_uncertainties = dydt_uncertainties*0.0
    y_uncertainties[0] = initial_y_uncertainty #We have no way to make an uncertainty for point 0, so we just use the same formula.
    for index in range(len(dydt_uncertainties)-1): #The uncertainty for each next point is propagated through the uncertainty of the current value and the delta_t*(dy/dt uncertainty), since we are adding two values.
        deltat_resolution = t_values[index+1]-t_values[index]
        y_uncertainties[index+1] = ((y_uncertainties[index])**2+(dydt_uncertainties[index]*deltat_resolution)**2)**0.5
    return y_uncertainties

#for calculating y at time t from dy/dt.  
def littleEulerGivenFunction(y_initial, deltat_resolution, dydtFunction, t_initial, t_final):
    numPoints = int((t_final-t_initial)/deltat_resolution)+1
    simulated_t_values = np.linspace(t_initial, t_final, numPoints)
    simulated_y_values = np.zeros(len(simulated_t_values)) #just initializing.
    dydt_values = np.zeros(len(simulated_t_values)) #just initializing.
    simulated_y_values[0] = y_initial
    for y_index in range(len(simulated_y_values)-1):
        localSlope = dydtFunction(simulated_t_values[y_index] ) 
        dydt_values[y_index]=localSlope
        simulated_y_values[y_index+1] = simulated_y_values[y_index] + localSlope * deltat_resolution
#        print(simulated_t_values[y_index+1], simulated_y_values[y_index+1], localSlope, localSlope * deltat_resolution)
#        print(simulated_y_values[y_index], simulated_t_values[y_index]*10-(simulated_t_values[y_index]**2)/2 +2)
    return simulated_t_values, simulated_y_values, dydt_values

def dydtNumericalExtraction(t_values, y_values, last_point_derivative = 0):
    lastIndex = len(simulated_y_values)-1
    delta_y_numerical = np.diff(np.insert(simulated_y_values,lastIndex,simulated_y_values[lastIndex])) #The diff command gives one less than what is fed in, so we insert the last value again. This gives a final value derivative of 0.
    delta_y_numerical[lastIndex] = last_point_derivative #now we set that last point to the optional argument.
    #It is ASSUMED that the t_values are evenly spaced.
    delta_t = t_values[1]-t_values[0]
    dydtNumerical = delta_y_numerical/delta_t
    return dydtNumerical
        
if __name__ == "__main__":
    pass
