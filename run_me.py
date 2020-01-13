import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
import copy
#import tprmodel
#from tprmodel import tprequation # moved from likelihood so import is not repeated at every sample EAW 2020/01/08
#import mumce_py.Project as mumce_pyProject #FIXME: Eric to fix plotting/graphing issue described in issue 9 -- https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/9
#import mumce_py.solution mumce_pySolution


#Now will automatically populate some variables from above.    
def parseUserInputParameters():
    UserInput.parameterNamesList = list(UserInput.parameterNamesAndMathTypeExpressionsDict.keys())
    UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
parseUserInputParameters()

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
    
# writeTPRModelFile() #####This line should be commented out if tprmodel.py is going to be edited manually.
from tprmodel import tprequation
    


class ip:
    #Ip for 'inverse problem'. Initialize prior chain starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
    def __init__(self, UserInput = UserInput):
        self.UserInput = UserInput
        self.verbose = UserInput.verbose
        if self.verbose: 
            print("Bayes Model Initialized")
        self.mcmc_length = UserInput.mcmc_length
        self.mcmc_burn_in = UserInput.mcmc_burn_in # Number of samples trimmed off the beginning of the Markov chain.
        self.mu_prior = UserInput.mu_prior
#        self.start_T = UserInput.T_0
#        self.beta_dTdt = UserInput.beta_dTdt
        self.cov_prior = UserInput.cov_prior
        self.Q_mu = self.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_cov = self.cov_prior/10 # Take small steps. <-- looks like this 20 should be a user defined variable.
#        self.initial_concentrations_array = UserInput.initial_concentrations_array
        self.modulate_accept_probability = UserInput.modulate_accept_probability
        
    def import_experimental_settings(self): #FIXME: This is obviously not very general. Though actually, we don't need it here. This should just go into UserInput as code rather than a function. These variables will become something like: UserInput.times, UserInput.observedResponse, UserInput.responseUncertainties.
        experiments_df = pd.read_csv(UserInput.Filename)
        self.times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
        self.experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
        self.errors = np.array(experiments_df['Errors']) #.to_numpy()

    #main function to get samples
    def MetropolisHastings(self):
        self.import_experimental_settings()
        rate_tot_array = np.zeros((self.mcmc_length,len(self.experiment)))
        samples = np.zeros((self.mcmc_length,len(self.mu_prior)))
        samples[0,:] = self.mu_prior # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
        likelihoods_vec = np.zeros((self.mcmc_length,1))
        posteriors_un_normed_vec = np.zeros((self.mcmc_length,1))
        priors_vec = np.zeros((self.mcmc_length,1))
        for i in range(1,self.mcmc_length):
            if self.verbose: print("MCMC sample number", i)
            proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov)
            prior_proposal = self.prior(proposal_sample)
            [likelihood_proposal, rate_tot_proposal] = self.likelihood(proposal_sample)
            prior_current_location = self.prior(samples[i-1,:])
            [likelihood_current_location, rate_tot_current_location] = self.likelihood(samples[i-1,:])
            accept_probability = (likelihood_proposal*prior_proposal)/(likelihood_current_location*prior_current_location) 
            if self.modulate_accept_probability != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.flatten_accept_probability)
                accept_probability = accept_probability**(1/N_flatten) #TODO: add code that unflattens the final histograms, that way even with more sampling we still get an accurate final posterior distribution. We can also then add a flag if the person wants to keep the posterior flattened.
            if accept_probability > np.random.uniform():  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                if self.verbose:
                  print('accept')
                  print(rate_tot_proposal)
                samples[i,:] = proposal_sample
                rate_tot_array[i,:] = rate_tot_proposal
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal #FIXME: Separate this block of code out into a helper function in the class, that way I can create another helper function for non-MCMC sampling.
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                if self.verbose:
                  print('reject')
                  print(rate_tot_current_location)
                samples[i,:] = samples[i-1,:]
                rate_tot_array[i,:] = rate_tot_current_location
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
        return [evidence, info_gain, samples, rate_tot_array, np.log(posteriors_un_normed_vec)] # EAW 2020/01/08
    def prior(self,discreteParameterVector):
        probability = multivariate_normal.pdf(x=discreteParameterVector,mean=self.mu_prior,cov=self.cov_prior)
        return probability
    def likelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        #Ashi has made a bunch of  FIXME below according to https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/11. 


        #The below few lines replace what used to look like this: #tpr_theta = odeint(tprequation, self.initial_concentrations_array, self.times, args = (*sample_list,self.beta_dTdt,self.start_T))      

        #tprEquationOutput = tprequation(tpr_theta, self.times, *sample_list, self.beta_dTdt,self.start_T) #FIXME: this will be in UserInput, only arguments and function and highest level function will be passed in.


#        simulationInputArguments = [tpr_theta, self.times, *sample_list, self.UserInput.beta_dTdt,self.UserInput.T_0] #FIXME: To be passed in from userInput
       
        #Want to achief: simulationOutput = simulationWrapper(discreteParameterVector)
        
        def simulationFunctionWrapper(discreteParameterVector): #FIXME: This should be defined in UserInput and passed in. User is responsible for it.
            # from tprmodel import tprequation # This is moved to the beginning of this file EAW 2020/01/08
            sample_list = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
            tpr_theta_Arguments = [tprequation, self.UserInput.initial_concentrations_array, self.times, (*sample_list,self.UserInput.beta_dTdt,self.UserInput.T_0) ] #FIXME: Times needs to occur in UserInput. This needs to all occur in some kind of UserFunctions module called from UserInput, should not be passed in here. 
            tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. #FIXME: initialArgs and equation should come from UserInput, not be hardcoded here.            
            simulationInputArguments = [tpr_theta, self.times, *sample_list, self.UserInput.beta_dTdt,self.UserInput.T_0] #FIXME: To be passed in from userInput
            simulationFunction = tprequation #FIXME: To be passed in from userInput
            simulationOutput = UserInput.model_function_name(*simulationInputArguments) # EAW 2020/01/08
            return simulationOutput
        
        simulationOutput = simulationFunctionWrapper(discreteParameterVector) #FIXME: code should look like simulationOutput = self.UserInput.simulationFunction(*self.UserInput.simulationInputArguments)

        #intermediate_metric = np.mean(np.square(rate_tot - self.experiment) / np.square(self.errors ))            
        
#        simulationOutput = simulationFunction(*simulationInputArguments) #FIXME: code should look like simulationOutput = self.UserInput.simulationFunction(*self.UserInput.simulationInputArguments)
        

        sample_list = list(discreteParameterVector) #converting to list so can use list expansion in arguments.        
        tpr_theta_Arguments = [tprequation, self.UserInput.initial_concentrations_array, self.times, (*sample_list,self.UserInput.beta_dTdt,self.UserInput.T_0) ] #FIXME: This needs to be passed in from UserInput, and called as self.UserInput.simulationInputArguments. right now it's hard coded here. The tuple is args. Right now, we will not support named arguments (maybe later I will add code to do it).
        tpr_theta = odeint(*tpr_theta_Arguments) # [0.5, 0.5] are the initial theta's. #FIXME: initialArgs and equation should come from UserInput, not be hardcoded here.
        rate = tprequation(tpr_theta, self.times, *sample_list, self.UserInput.beta_dTdt,self.UserInput.T_0) #This is old code for Eric to delete when he understands the code and has fixed the "return" part of the likelihood.
        rate_tot = -np.sum(rate, axis=0) #This is old code for Eric to delete when he understands the code and has fixed the "return" part of the likelihood.
        #temp_points = np.array([0,49,99,149]) #range(225)        #This was a temporary line that should be deleted once Eric understands code flow.
        #simulatedResponses = np.log10(rate_tot[temp_points]) #This was a temporary line that should be deleted once Eric understands code flow.
        
        #FIXME: the below wrapper functions can be combined, but either way should be in userinput as described below.
        def rate_tot_summing_func(rate):  
            rate_tot = -np.sum(rate, axis=0)   
            return rate_tot
        def rate_tot_four_points_func(rate): #Multiple layers of wrapper functions are fine.
            rate_tot = rate_tot_summing_func(rate)
            temp_points = np.array([0,49,99,149]) #range(225) #FIXME: There should be nothing hard coded here. You can hard code it in userinput if you want.
            rate_tot_four_points = np.array(rate_tot[temp_points])
            return rate_tot_four_points
        def log10_wrapper_func(rate):
            rate_tot_four_points = rate_tot_four_points_func(rate)
            loggedRateValues = np.log10(rate_tot_four_points)
            return loggedRateValues
            
        simulationOutputProcessingFunction = log10_wrapper_func #FIXME: this will become fed in as self.UserInput.simulationOutputProcessingFunction
        simulatedResponses = simulationOutputProcessingFunction(simulationOutput) #This will become simulatedResponses = self.UserInput.simulationOutputProcessingFunction(simulationOutput)
            
        temp_points = np.array([0,49,99,149])
        observedResponses = np.log10(self.experiment[temp_points]) #FIXME: This should not be hard coded here. Should be self.UserInput.obseredResponses
        cov = self.errors[temp_points]
        #probability_metric = multivariate_normal.pdf(x=np.log10(rate_tot[temp_points]),mean=np.log10(self.experiment[temp_points]),cov=self.errors[temp_points]) #ERIC, THIS WAS THE PREVIOUS LINE. I BELIEVE YOUR LOG10 TRANSFORM IS NOT SOMETHING WE WOULD NORMALLY DO. AM I CORRECT?
        probability_metric = multivariate_normal.pdf(x=simulatedResponses,mean=observedResponses,cov=cov) #FIXME: should become self.UserInput.responseUncertantiesCov or something like that.
        if self.verbose: print('likelihood probability',probability_metric,'log10(rate_tot)',np.log10(rate_tot[temp_points]), 'log10(experiment)', np.log10(self.experiment[temp_points]), 'error', self.errors[temp_points])
        return probability_metric, rate_tot #FIXME: This needs to say probability_metric, simulatedResponses or something like that, but right now the sizes of the arrays do not match.
        
    
ip_object = ip()
[evidence, info_gain, samples, rate_tot_array, logP] = ip_object.MetropolisHastings()
############################################# The computation portion is contained above.
np.savetxt('logP.csv',logP) # EAW 2020/01/08
post_mean = np.mean(samples, axis=0)
experiments_df = pd.read_csv(UserInput.Filename)
start_T = UserInput.T_0 #this is the starting temperature.
times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/2000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
errors = np.array(experiments_df['Errors']) #.to_numpy()
post_mean_list = list(post_mean) #converting to list so can use list expansion in arguments.
#The simulation outputs are passed in the 'rate_tot_array', and so calling the model again in the following lines is no longer necessary.
#tpr_theta = odeint(tprequation, UserInput.initial_concentrations_array, times, args = (*post_mean_list,UserInput.beta_dTdt, start_T)) # [0.5, 0.5] are the initial theta's.
#rate = tprequation(tpr_theta, times, *post_mean_list, UserInput.beta_dTdt, start_T)
#rate_tot = -np.sum(rate, axis=0) 

fig0, ax0 = plt.subplots()
if UserInput.verbose:
  print(np.mean(rate_tot_array,axis = 0))
ax0.plot(np.array(experiments_df['AcH - T']),np.mean(rate_tot_array,axis = 0), 'r')
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


#Make histograms for each parameter. Need to make some dictionaries where relevant objects will be stored.
sampledParameterFiguresDictionary = copy.deepcopy(UserInput.parameterNamesAndMathTypeExpressionsDict)
sampledParameterAxesDictionary = copy.deepcopy(UserInput.parameterNamesAndMathTypeExpressionsDict)
for key in UserInput.parameterNamesAndMathTypeExpressionsDict:
    parameterName = key
    sampledParameterHistogramMaker(parameterName,UserInput.parameterNamesAndMathTypeExpressionsDict, sampledParameterFiguresDictionary, sampledParameterAxesDictionary)
