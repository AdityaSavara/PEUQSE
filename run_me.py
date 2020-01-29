import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd

import copy
#import tprmodel
#from tprmodel import tprequation # moved from likelihood so import is not repeated at every sample EAW 2020/01/08
#import mumce_py.Project as mumce_pyProject #FIXME: Eric to fix plotting/graphing issue described in issue 9 -- https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/9
#import mumce_py.solution mumce_pySolution





#TURNING OFF THE WRITING OF FUNCTION FEATURE RIGHT NOW, TO PREVENT ERIC FROM HAVING TO WORRY ABOUT IT IF IT AFFECTS HIS WORK
# def writeTPRModelFile():
    # with open("tprmodel.py", "w") as myfile:
        # myfile.write("\n\
# import numpy as np\n\
# \n\
# # The below function must return a vector of rates. \n\
# def tprequation(tpr_theta,t," + UserInput.stringOfParameterNames + ",beta_dTdt,start_T): #beta_dTdT is the heating rate. \n\
    # if tpr_theta.ndim == 1:  #for consistency, making tpr_theta a 2D array if it does not start as 2D. \n\
        # tpr_theta2D = np.atleast_2d(tpr_theta)  \n\
    # if tpr_theta.ndim == 2: \n\
        # tpr_theta2D = np.array(tpr_theta) \n\
    # #Now find out how many species concentrations there are from the data: \n\
    # num_of_concentrations = len(tpr_theta2D[0]) \n\
    # Ea_Array = [Ea_1,Ea_2] \n\
    # log_A_array = [log_A1, log_A2] \n\
    # gamma_array = [gamma1, gamma2] \n\
    # \n\
    # T = start_T + beta_dTdt*t  \n\
    # kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1  \n\
    # \n\
    # ratesList = [] \n\
    # for rateIndex in range(num_of_concentrations): \n\
        # rate = -tpr_theta2D[:,rateIndex]*np.exp(-(Ea_Array[rateIndex]-kB*T*log_A_array[rateIndex]-gamma_array[rateIndex]*tpr_theta2D[:,rateIndex])/(kB*T))  \n\
    # \n\
        # #Shortened below to one line (above) \n\
        # # theta_i = tpr_theta2D[:,rateIndex] \n\
        # # Ea_i = Ea_Array[rateIndex] \n\
        # # log_A_i = log_A_array[rateIndex] \n\
        # # gamma_i = gamma_array[rateIndex] \n\
        # # rate = -theta_i*np.exp(-(Ea_i-kB*T*log_A_i-gamma_i*theta_i)/(kB*T))  \n\

      # \n\
        # #The above expression is the general form of this: rate_2 = -theta_2*np.exp(-(Ea_2-kB*T*log_A2-gamma2*theta_2)/(kB*T))       \n\
        # ratesList.append(rate) \n\
    # \n\
    # if tpr_theta.ndim == 1: \n\
        # ratesList = list(np.array(ratesList).flatten()) #for some reason, needs to be flattened for the MCMC. \n\
    # return ratesList \n\

    # ")
    


class ip:
    #Ip for 'inverse problem'. Initialize prior chain starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
    def __init__(self, UserInput = None):
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.parameterNamesAndMathTypeExpressionsDict.keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        self.UserInput = UserInput
        if self.UserInput.verbose: 
            print("Bayes Model Initialized")
        self.mcmc_length = UserInput.mcmc_length
        self.mcmc_burn_in = UserInput.mcmc_burn_in # Number of samples trimmed off the beginning of the Markov chain.
        self.mu_prior = UserInput.mu_prior
        self.cov_prior = UserInput.cov_prior
        self.Q_mu = self.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_cov = self.cov_prior/10 # Take small steps. <-- looks like this 20 should be a user defined variable.
#        self.initial_concentrations_array = UserInput.initial_concentrations_array
        self.modulate_accept_probability = UserInput.modulate_accept_probability
        

    #main function to get samples
    def MetropolisHastings(self):
        if hasattr(self.UserInput, "mcmc_random_seed"):
            if type(UserInput.mcmc_random_seed) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(UserInput.mcmc_random_seed)
        self.UserInput.import_experimental_settings(UserInput.Filename) #FIXME: This needs to get out of this function.
        samples_simulatedOutputs = np.zeros((self.mcmc_length,self.UserInput.num_data_points))
        samples = np.zeros((self.mcmc_length,len(self.mu_prior)))
        samples[0,:] = self.mu_prior # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
        likelihoods_vec = np.zeros((self.mcmc_length,1))
        posteriors_un_normed_vec = np.zeros((self.mcmc_length,1))
        priors_vec = np.zeros((self.mcmc_length,1))
        for i in range(1,self.mcmc_length):
            if self.UserInput.verbose: print("MCMC sample number", i)
            proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov)
            prior_proposal = self.prior(proposal_sample)
            [likelihood_proposal, simulationOutput_proposal] = self.likelihood(proposal_sample)
            prior_current_location = self.prior(samples[i-1,:])
            [likelihood_current_location, rate_tot_current_location] = self.likelihood(samples[i-1,:])
            accept_probability = (likelihood_proposal*prior_proposal)/(likelihood_current_location*prior_current_location) 
            if self.modulate_accept_probability != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.flatten_accept_probability)
                accept_probability = accept_probability**(1/N_flatten) #TODO: add code that unflattens the final histograms, that way even with more sampling we still get an accurate final posterior distribution. We can also then add a flag if the person wants to keep the posterior flattened.
            if accept_probability > np.random.uniform():  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                if self.UserInput.verbose:
                  print('accept')
                  #print(simulationOutput_proposal)
                samples[i,:] = proposal_sample
                samples_simulatedOutputs[i,:] = simulationOutput_proposal
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal #FIXME: Separate this block of code out into a helper function in the class, that way I can create another helper function for non-MCMC sampling.
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                if self.UserInput.verbose:
                  print('reject')
                  #print(rate_tot_current_location)
                samples[i,:] = samples[i-1,:]
                samples_simulatedOutputs[i,:] = rate_tot_current_location
                posteriors_un_normed_vec[i] = likelihood_current_location*prior_current_location
                likelihoods_vec[i] = likelihood_current_location
                priors_vec[i] = prior_current_location
            ########################################
        samples = samples[self.mcmc_burn_in:]
        posteriors_un_normed_vec = posteriors_un_normed_vec[self.mcmc_burn_in:]
        likelihoods_vec = likelihoods_vec[self.mcmc_burn_in:]
        priors_vec = priors_vec[self.mcmc_burn_in:]
        # posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        evidence = np.mean(posteriors_un_normed_vec)*np.sqrt(2*np.pi*np.std(samples)**2)
        posteriors_vec = posteriors_un_normed_vec/evidence
        log_ratios = np.log(posteriors_vec/priors_vec)
        log_ratios[np.isinf(log_ratios)] = 0
        log_ratios = np.nan_to_num(log_ratios)
        info_gain = np.mean(log_ratios)
        if self.UserInput.verbose == True:
            print("line 146", (evidence), (info_gain), len(samples), len(samples_simulatedOutputs), len(np.log(posteriors_un_normed_vec)))
            #mcmc_output = np.vstack((evidence, info_gain, samples, samples_simulatedOutputs, np.log(posteriors_un_normed_vec)))
            #np.savetxt('mcmc_output.csv',mcmc_output) # EAW 2020/01/08
        return [evidence, info_gain, samples, samples_simulatedOutputs, np.log(posteriors_un_normed_vec)] # EAW 2020/01/08
    def prior(self,discreteParameterVector):
        probability = multivariate_normal.pdf(x=discreteParameterVector,mean=self.mu_prior,cov=self.cov_prior)
        return probability
    def likelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        #This is more in accordance with https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/11. 
        simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction
        simulationFunctionWrapper =  self.UserInput.simulationFunctionWrapper

        simulationOutput =simulationFunctionWrapper(discreteParameterVector) #FIXME: code should look like simulationOutput = self.UserInput.simulationFunction(*self.UserInput.simulationInputArguments)
        
        
        #simulationOutputProcessingFunction = self.UserInput.log10_wrapper_func #FIXME: this will become fed in as self.UserInput.simulationOutputProcessingFunction
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput
        if type(simulationOutputProcessingFunction) != type(None):
            print("getting here")
            simulatedResponses = self.UserInput.simulationOutputProcessingFunction(simulationOutput) #This will become simulatedResponses = self.UserInput.simulationOutputProcessingFunction(simulationOutput)

            
        #temp_points = np.array([0,49,99,149]) moved to UserInput EAW 2020/01/27
        observedResponses = self.UserInput.observedResponses()#np.log10(self.experiment[temp_points]) #FIXME: This should not be hard coded here. Should be self.UserInput.obseredResponses
        
        #To find the relevant covariance, we take the errors from the points.
        cov = UserInput.errors[UserInput.temp_points] #FIXME: We should not be doing subset of points like this here. Should happen at user input level.
        #FIXME should become:
        #cov = self.UserInput.observedResponses_uncertainties
        
        
        #probability_metric = multivariate_normal.pdf(x=np.log10(rate_tot[temp_points]),mean=np.log10(self.experiment[temp_points]),cov=self.errors[temp_points]) #ERIC, THIS WAS THE PREVIOUS LINE. I BELIEVE YOUR LOG10 TRANSFORM IS NOT SOMETHING WE WOULD NORMALLY DO. AM I CORRECT?
        probability_metric = multivariate_normal.pdf(x=simulatedResponses,mean=observedResponses,cov=cov) #FIXME: should become self.UserInput.responseUncertantiesCov or something like that.
        simulationOutput = self.UserInput.rate_tot_summing_func(simulationOutput)
        #temp_points = self.UserInput.temp_points
        if self.UserInput.verbose: print('likelihood probability',probability_metric)
        #if self.UserInput.verbose: print('likelihood probability',probability_metric,'log10(rate_tot)',np.log10(rate_tot[temp_points]), 'log10(experiment)', np.log10(self.experiment[temp_points]), 'error', self.errors[temp_points])
        return probability_metric, simulationOutput #FIXME: This needs to say probability_metric, simulatedResponses or something like that, but right now the sizes of the arrays do not match.


def sampledParameterHistogramMaker(parameterName,parameterNamesAndMathTypeExpressionsDict, sampledParameterFiguresDictionary, sampledParameterAxesDictionary):
        parameterIndex = list(parameterNamesAndMathTypeExpressionsDict).index(parameterName)
        sampledParameterFiguresDictionary['parameterName'], sampledParameterAxesDictionary['parameterName'] = plt.subplots()   #making plt objects    
        sampledParameterAxesDictionary['parameterName'].hist(samples[:,parameterIndex]) #filling the object with data
        #setting the labels etc. and then exporting.
        sampledParameterAxesDictionary['parameterName'].set_ylabel('frequency')
        sampledParameterAxesDictionary['parameterName'].set_xlabel(UserInput.parameterNamesAndMathTypeExpressionsDict[parameterName])
        sampledParameterFiguresDictionary['parameterName'].tight_layout()
        sampledParameterFiguresDictionary['parameterName'].savefig(parameterName+'.png', dpi=220)

        #The above block makes code kind of like this in a dynamic fashion. Since we know how many we will need, a dictionary is used to avoid the need for 'exec' statements when making new parameters.
        # fig2, ax2 = plt.subplots()
        # ax2.hist(samples[:,1])
        # ax2.set_ylabel('frequency')
        # ax2.set_xlabel(r'$E_{a2}$')
        # fig2.tight_layout()
        # fig2.savefig('Ea2.png', dpi=220)
        

if __name__ == "__main__":
    import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
    UserInput.verbose = True    
    UserInput.mcmc_burn_in = 500
    UserInput.mcmc_length = 1000
    ip_object = ip(UserInput)
    [evidence, info_gain, samples, samples_simulatedOutputs, logP] = ip_object.MetropolisHastings()
    ############################################# The computation portion is contained above.
#    print(samples)
    np.savetxt('logP.csv',logP) # EAW 2020/01/08
    #post_mean = np.mean(samples, axis=0)
    experiments_df = pd.read_csv(UserInput.Filename)
#    print(experiments_df)
#    sys.exit()
    #start_T = UserInput.T_0 #this is the starting temperature.
    #times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
    #experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
    #errors = np.array(experiments_df['Errors']) #.to_numpy()
    #post_mean_list = list(post_mean) #converting to list so can use list expansion in arguments.
    #The simulation outputs are passed in the 'samples_simulatedOutputs', and so calling the model again in the following lines is no longer necessary.
    #tpr_theta = odeint(tprequation, UserInput.initial_concentrations_array, times, args = (*post_mean_list,UserInput.beta_dTdt, start_T)) # [0.5, 0.5] are the initial theta's.
    #rate = tprequation(tpr_theta, times, *post_mean_list, UserInput.beta_dTdt, start_T)
    #rate_tot = -np.sum(rate, axis=0) 
    
    fig0, ax0 = plt.subplots()
    if UserInput.verbose:
      print(np.mean(samples_simulatedOutputs,axis = 0))
    ax0.plot(np.array(experiments_df['AcH - T']),np.mean(samples_simulatedOutputs,axis = 0), 'r')
    ax0.plot(np.array(experiments_df['AcH - T']),np.array(experiments_df['AcHBackgroundSubtracted'])/2000,'g')
    ax0.set_ylim([0.00, 0.025])
    ax0.set_xlabel('T (K)')
    ax0.set_ylabel(r'$rate (s^{-1})$')
    ax0.legend(['model posterior', 'experiments'])
    fig0.tight_layout()
    fig0.savefig('tprposterior.png', dpi=220)
    posterior_df = pd.DataFrame(samples,columns=[UserInput.parameterNamesAndMathTypeExpressionsDict[x] for x in UserInput.parameterNamesList])
    pd.plotting.scatter_matrix(posterior_df)
    plt.savefig('scatter_matrix_posterior.png',dpi=220)
    #######Below function and codee makes the Histograms for each parameter#####
    
    
    
    #Make histograms for each parameter. Need to make some dictionaries where relevant objects will be stored.
    sampledParameterFiguresDictionary = copy.deepcopy(UserInput.parameterNamesAndMathTypeExpressionsDict)
    sampledParameterAxesDictionary = copy.deepcopy(UserInput.parameterNamesAndMathTypeExpressionsDict)
    for key in UserInput.parameterNamesAndMathTypeExpressionsDict:
        parameterName = key
        sampledParameterHistogramMaker(parameterName,UserInput.parameterNamesAndMathTypeExpressionsDict, sampledParameterFiguresDictionary, sampledParameterAxesDictionary)
