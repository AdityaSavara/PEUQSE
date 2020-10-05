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
try:
    import CiteSoft
except:
    import os #The below lines are to allow CiteSoftLocal to be called regardless of user's working directory.
    lenOfFileName = len(os.path.basename(__file__)) #This is the name of **this** file.
    absPathWithoutFileName = os.path.abspath(__file__)[0:-1*lenOfFileName]
    sys.path.append(absPathWithoutFileName)
    import CiteSoftLocal as CiteSoft

class parameter_estimation:
    #Inside this class, a UserInput namespace is provided. This has dictionaries of UserInput choices.
    #However, the code initally parses those choices and then puts processed versions in the SAME name space, but no longer in the dictionaries.
    #So functions in this class should (when possible) call the namespace variables that are not in dictionaries, unless the original userinput is desired.
    #'inverse problem'. Initialize chain with initial guess (prior if not provided) as starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  
    
    
    software_name = "CheKiPEUQ Bayesian Parameter Estimation"
    software_version = "1.0.0"
    software_unique_id = "https://doi.org/10.1002/cctc.202000953"
    software_kwargs = {"version": software_version, "author": ["Aditya Savara", "Eric A. Walker"], "doi": "https://doi.org/10.1002/cctc.202000953", "cite": "Savara, A. and Walker, E.A. (2020), CheKiPEUQ Intro 1: Bayesian Parameter Estimation Considering Uncertainty or Error from both Experiments and Theory. ChemCatChem. Accepted. doi:10.1002/cctc.202000953"} 
    @CiteSoft.module_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def __init__(self, UserInput = None):
        self.UserInput = UserInput #Note that this is a pointer, so the later lines are within this object.
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.model['parameterNamesAndMathTypeExpressionsDict'].keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        UserInput.parameterNamesAndMathTypeExpressionsDict = UserInput.model['parameterNamesAndMathTypeExpressionsDict']
        if self.UserInput.parameter_estimation_settings['verbose']: 
            print("Paremeter Estimation Object Initialized")
        
        #Setting this object so that we can make changes to it below without changing userinput dictionaries.
        self.UserInput.mu_prior = np.array(UserInput.model['InputParameterPriorValues']) 
        
        #Below code is mainly for allowing uniform distributions in priors.
        UserInput.InputParametersPriorValuesUncertainties = np.array(UserInput.model['InputParametersPriorValuesUncertainties'],dtype='float32') #Doing this so that the -1.0 check below should work.
        if -1.0 in UserInput.InputParametersPriorValuesUncertainties: #This means that at least one of the uncertainties has been set to "-1" which means a uniform distribution. 
            UserInput.InputParametersPriorValuesUniformDistributionsIndices = [] #intializing.
            if len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) != 1:
                print("A value of '-1' in the uncertainties signifies a uniform distribution for CheKiPEUQ. As of July 1st 2020, the uniform distribution feature is only compatible with a 1D of array for uncertainties and not compatible with providing a full covariance matrix. If you need such a feature, contact the developers because it could be implemented. Eventually, a more sophisiticated back end may be used which would allow such a feature.")
            # If there is a uniform distribution, that means two actions need to be taken:
             #First, we will populate InputParametersPriorValuesUncertainties with the standard deviation of a uniform distribution. This is so that the MCMC steps can be taken of the right size.
             #Second, that we will need to make a custom calculation when calculating the prior probability that effectively excludes this variable.  So we'll create an array of indices to help us with that.        
            #We will do both in a loop.
            UserInput.InputParametersPriorValuesUniformDistributionsKey  = UserInput.InputParametersPriorValuesUncertainties *1.0 #Just initalizing
            for parameterIndex, uncertaintyValue in enumerate(UserInput.InputParametersPriorValuesUncertainties):
                if uncertaintyValue == -1.0:
                    UserInput.InputParametersPriorValuesUniformDistributionsKey[parameterIndex] = 1.0 #This is setting the parameter as "True" for having a uniform distribution. 
                    UserInput.InputParametersPriorValuesUniformDistributionsIndices.append(parameterIndex)
                    #In the case of a uniform distribution, the standard deviation and variance are given by sigma = (b−a)/ √12 :   
                    #See for example  https://www.quora.com/What-is-the-standard-deviation-of-a-uniform-distribution-How-is-this-formula-determined
                    std_prior_single_parameter = (UserInput.model['InputParameterPriorValues_upperBounds'][parameterIndex] - UserInput.model['InputParameterPriorValues_lowerBounds'][parameterIndex])/(12**0.5)
                    UserInput.InputParametersPriorValuesUncertainties[parameterIndex] = std_prior_single_parameter #Note that going forward the array InputParametersPriorValuesUncertainties cannot be checked to see if the parameter is from a uniform distribution. Instead, InputParametersPriorValuesUniformDistributionsKey must be checked. 
                    #We will also fill the model['InputParameterPriorValues'] to have the mean of the two bounds. This can matter for some of the scaling that occurs later.
                    self.UserInput.mu_prior[parameterIndex] = (UserInput.model['InputParameterPriorValues_upperBounds'][parameterIndex] + UserInput.model['InputParameterPriorValues_lowerBounds'][parameterIndex])/2
        
        #Now to make covmat. Leaving the original dictionary object intact, but making a new object to make covmat_prior.
        if len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) == 1 and (len(UserInput.InputParametersPriorValuesUncertainties) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
            UserInput.std_prior = np.array(UserInput.InputParametersPriorValuesUncertainties, dtype='float32') #using 32 since not everyone has 64.
            UserInput.var_prior = np.power(UserInput.InputParametersPriorValuesUncertainties,2)
            UserInput.covmat_prior = np.diagflat(self.UserInput.var_prior) 
        elif len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) > 1: #If it's non-1D, we assume it's already a covariance matrix.
            UserInput.covmat_prior = np.array(UserInput.InputParametersPriorValuesUncertainties, dtype='float32')
            UserInput.var_prior = np.diagonal(UserInput.covmat_prior)
            UserInput.std_prior = np.power(UserInput.covmat_prior,0.5)
        else: #If a blank list is received, that means the user
            print("The covariance matrix of the priors is undefined because InputParametersPriorValuesUncertainties is blank.")
        #    cov_prior = np.array([[200.0, 0., 0., 0., 0., 0.], 
        #                          [0., 200.0, 0., 0., 0., 0.],
        #                          [0., 0., 13.0, 0., 0., 0.],
        #                          [0., 0., 0., 13.0, 0., 0.],
        #                          [0., 0., 0., 0., 0.1, 0.],
        #                          [0., 0., 0., 0., 0., 0.1]])


        #Making things at least 2d.  Also changing it to a purely internal variable because that way we don't edit the user input dictionary going forward.
        self.samples_of_prior = np.random.multivariate_normal(self.UserInput.mu_prior,UserInput.covmat_prior,UserInput.parameter_estimation_settings['mcmc_length'])
        UserInput.responses_abscissa = np.atleast_2d(UserInput.responses['responses_abscissa'])
        UserInput.responses_observed = np.atleast_2d(UserInput.responses['responses_observed'])
        UserInput.responses_observed_uncertainties = np.atleast_2d(UserInput.responses['responses_observed_uncertainties'])

        #TODO: This currently only is programmed for if the uncertainties are uncorrelated standard deviaions (so is not compatible with a directly fed cov_mat). Also, we need to figure out what do when they are not gaussian/symmetric.
        UserInput.responses_observed_transformed, UserInput.responses_observed_transformed_uncertainties  = self.transform_responses(UserInput.responses_observed, UserInput.responses_observed_uncertainties) #This creates transforms for any data that we might need it. The same transforms will also be applied during parameter estimation.
            
        #The below unusual code is because during doeParameterModulationCombinationsScanner, populate synthetic data calls init again.
        #So we will only call populateIndependentVariablesFunction if we're not in the middle of design of experiments.
        if not hasattr(self, 'middle_of_doe_flag'): #We check of the middle_of_doe_flag exists. #If the flag is not there and the populate function exists, we call it.
            if UserInput.model['populateIndependentVariablesFunction'] != None:
                UserInput.model['populateIndependentVariablesFunction'](UserInput.responses['independent_variables_values']) 
        if hasattr(self, 'middle_of_doe_flag'): #We check of the middle_of_doe_flag exists. If it's there, no problem.
            if self.middle_of_doe_flag == False: #If the flag is there, we only proceed to call the function if the flag is set to false.
                if UserInput.model['populateIndependentVariablesFunction'] != None:
                    UserInput.model['populateIndependentVariablesFunction'](UserInput.responses['independent_variables_values']) 
        
        self.UserInput.num_data_points = len(UserInput.responses_observed.flatten()) #FIXME: This is only true for transient data.
        #Now scale things as needed:
        if UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "off":
            self.UserInput.mu_prior_scaled = UserInput.mu_prior*1.0
            self.UserInput.var_prior_scaled = UserInput.var_prior*1.0
            self.UserInput.covmat_prior_scaled = UserInput.covmat_prior*1.0
        else:            
            if UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "std":
                self.UserInput.scaling_uncertainties = UserInput.std_prior #Could also be by mu_prior.  The reason a separate variable is made is because this will be used in the getPrior function as well, and having a separate variable makes it easier to trace. This scaling helps prevent numerical errors in returning the pdf.
            elif UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "mu":
                self.UserInput.scaling_uncertainties = UserInput.mu_prior
            else: #Else we assume that UserInput.parameter_estimation_settings['scaling_uncertainties_type'] has been set to a fixed float or vector. For now, we'll just support float.
                scaling_factor = float(UserInput.parameter_estimation_settings['scaling_uncertainties_type'])
                self.UserInput.scaling_uncertainties = (UserInput.mu_prior/UserInput.mu_prior)*scaling_factor #This basically makes a vector of ones times the scaling factor.
            #TODO: consider a separate scaling for each variable, taking the greater of either mu_prior or std_prior.
            #TODO: Consider changing how self.UserInput.scaling_uncertainties is done to accommodate greater than 1D vector. Right now we use np.shape(self.UserInput.scaling_uncertainties)[0]==1, but we could use np.shape(self.UserInput.scaling_uncertainties)==np.shape(UserInput.mu_prior)
            if np.shape(np.atleast_2d(self.UserInput.scaling_uncertainties))[0]==1:  #In this case, the uncertainties is not a covariance matrix.
                pass
            elif np.shape(np.atleast_2d(self.UserInput.scaling_uncertainties))[0]==np.shape(np.atleast_2d(self.UserInput.scaling_uncertainties))[1]: #In his case, the uncertainties are a covariance matrix so we take the diagonal (which are variances) and the square root of them.
                self.UserInput.scaling_uncertainties = (np.diagonal(self.UserInput.scaling_uncertainties))**0.5 #Take the diagonal which is variances, and            
            else:
                print("There is an unsupported shape somewhere in the prior.  The prior is currently expected to be 1 dimensional.")
                print(np.shape(self.UserInput.scaling_uncertainties))
                sys.exit()
            self.UserInput.mu_prior_scaled = np.array(UserInput.mu_prior/UserInput.scaling_uncertainties)
            self.UserInput.var_prior_scaled = np.array(UserInput.var_prior/(UserInput.scaling_uncertainties*UserInput.scaling_uncertainties))
            self.UserInput.covmat_prior_scaled = self.UserInput.covmat_prior*1.0 #First initialize, then fill.
            for parameterIndex, parameterValue in enumerate(UserInput.scaling_uncertainties):
                UserInput.covmat_prior_scaled[parameterIndex,:] = UserInput.covmat_prior[parameterIndex,:]/parameterValue
                #The next line needs to be on UserInput.covmat_prior_scaled and not UserInput.covmat_prior, since we're stacking the divisions.
                UserInput.covmat_prior_scaled[:,parameterIndex] = UserInput.covmat_prior_scaled[:,parameterIndex]/parameterValue        
        
        #To find the *observed* responses covariance matrix, meaning based on the uncertainties reported by the users, we take the uncertainties from the points. This is needed for the likelihood. However, it will be transformed again at that time.
        self.UserInput.num_response_dimensions = np.shape(UserInput.responses_abscissa)[0] #The first index of shape is the num of responses, but has to be after at_least2d is performed.
        self.observed_responses_covmat_transformed = returnShapedResponseCovMat(self.UserInput.num_response_dimensions, self.UserInput.responses_observed_transformed_uncertainties)
                       
        #self.covmat_prior = UserInput.covmat_prior
        self.Q_mu = self.UserInput.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_covmat = self.UserInput.covmat_prior # Take small steps. 
        #Getting initial guess of parameters and populating the internal variable for it.
        if 'InputParameterInitialGuess' not in self.UserInput.model: #if an initial guess is not provided, we use the prior.
            self.UserInput.model['InputParameterInitialGuess'] = self.UserInput.mu_prior
        #From now, we switch to self.UserInput.InputParameterInitialGuess because this is needed in case we're going to do reducedParameterSpace.
        self.UserInput.InputParameterInitialGuess = self.UserInput.model['InputParameterInitialGuess']
        #Now populate the simulation Functions. #NOTE: These will be changed if a reduced parameter space is used.
        self.UserInput.simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        self.UserInput.simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']
    
        #Now reduce the parameter space if requested by the user. #Considered having this if statement as a function called outside of init.  However, using it in init is the best practice since it forces correct ordering of reduceParameterSpace and reduceResponseSpace
        if len(self.UserInput.model['reducedParameterSpace']) > 0:
            print("Important: the UserInput.model['reducedParameterSpace'] is not blank. That means the only parameters allowed to change will be the ones in the indices inside 'reducedParameterSpace'.   All others will be held constant.  The values inside  'InputParameterInitialGuess will be used', and 'InputParameterPriorValues' if an initial guess was not provided.")
            self.reduceParameterSpace()
    
        #Now reduce the parameter space if requested by the user. #Considered having this if statement as a function called outside of init.  However, using it in init is the best practice since it forces correct ordering of reduceParameterSpace and reduceResponseSpace
        #This code must be **after** the reduceParameterSpace because this makes a wrapper for the simulationOutputProcessingFunction
        if len(self.UserInput.responses['reducedResponseSpace']) > 0:
            print("Important: the UserInput.model['reducedResponseSpace'] is not blank. That means the only responses examined will be the ones in the indices inside 'reducedReponseSpace'.   The values of all others will be discarded during each simulation.")
            self.reduceResponseSpace()
            
    
    def reduceResponseSpace(self):
        #This function has no explicit arguments, but takes everything in self.UserInput as an implied argument.
        #In particular, self.UserInput.responses['reducedResponseSpace']
        #it has two implied returns: 1) self.UserInput.simulationOutputProcessingFunction, 2) self.responses_covmat becomes reduced in size.
        
        UserInput = self.UserInput
        #First, we need to make a function that is going to reduce the dimensionality of the outputs outputs when there are simulations.
        #Make a deep copy of the existing function, so that we can use it if needed.
        self.UserInput.beforeReducedResponseSpaceSimulationOutputProcessingFunction = copy.deepcopy(self.UserInput.simulationOutputProcessingFunction)
        self.UserInput.beforeReducedResponseSpace_num_response_dimensions = self.UserInput.num_response_dimensions
        def extractReducedResponsesOutputsWrapper(simulatedOutput):
            #The simulatedOuput is an exlicit argument, the self.UserInput.model['reducedResponseSpace'] is an implicit argument.    
            #First, check if there is an OutputProcessing function to use on the simulatedOutput.
            if type(self.UserInput.beforeReducedResponseSpaceSimulationOutputProcessingFunction) != type(None):
                fullResponseOutput = self.UserInput.beforeReducedResponseSpaceSimulationOutputProcessingFunction(simulatedOutput) #We use the processing function to convert the simulated output to the actual responses, then we trim them as above.
            elif type(self.UserInput.beforeReducedResponseSpaceSimulationOutputProcessingFunction) == type(None): #if not, we take the output directly.
                fullResponseOutput = simulatedOutput
                                    
            #We could calculate the number of responses from fullResponseOutput, but we use self.UserInput.beforeReducedResponseSpace_num_response_dimensions as an implicit argument.
            reducedResponseOutput = []#Just intializing, then will append to it.
            for responseDimIndex in range(self.UserInput.beforeReducedResponseSpace_num_response_dimensions):
                #We'll only keep a responsDim if the responseDimIndex is named in self.UserInput.model['reducedResponseSpace']
                if responseDimIndex in self.UserInput.responses['reducedResponseSpace']:
                    reducedResponseOutput.append(fullResponseOutput[responseDimIndex])
            return reducedResponseOutput

        #Now get our first "implied return" by using the above function as the processing function.
        self.UserInput.simulationOutputProcessingFunction = extractReducedResponsesOutputsWrapper
    
        #Now we get our second "implied return" by reducing the response_abscissa, transformed response values, and their uncertainties.
        #TODO: consider making a different variable so that the dictionary does not need to get overwritten.
        self.UserInput.responses_abscissa = returnReducedIterable(self.UserInput.responses_abscissa, self.UserInput.responses['reducedResponseSpace'])
        self.UserInput.responses_observed_transformed = returnReducedIterable(self.UserInput.responses_observed_transformed, self.UserInput.responses['reducedResponseSpace'])
        self.UserInput.responses_observed_transformed_uncertainties = returnReducedIterable(self.UserInput.responses_observed_transformed, self.UserInput.responses['reducedResponseSpace'])
        self.UserInput.num_response_dimensions = np.shape(UserInput.responses_abscissa)[0]
    
        #Now we get our third "implied return" by reducing the response_covmat.
        self.observed_responses_covmat_transformed = returnReducedIterable(self.observed_responses_covmat_transformed, self.UserInput.responses['reducedResponseSpace'])
        return

    
    #This function reduces the parameter space. The only parameters allowed to change will be the ones in the indices inside 'reducedParameterSpace'.   All others will be held constant.  The values inside  'InputParameterInitialGuess will be used', and 'InputParameterPriorValues' if an initial guess was not provided.")    
    #These lines of code started in __init__ was moved outside of initializing the class so that someday people can call it later on after making the class object, if desired.
    #That way people can change to a different reduced parameter space without making a new object by updating what is in UserInput.model['reducedParameterSpace']
    #However, that 'later changing' is currently not supported. The indices *at present* only work out correctly when this is called at end of initialization.
    def reduceParameterSpace(self): 
        UserInput = self.UserInput
        
        self.UserInput.simulationFunction = self.simulateWithSubsetOfParameters #Now simulateWithSubsetOfParameters will be called as the simulation function.
        self.UserInput.simulationOutputProcessingFunction = None #We will use self.UserInput.model['simulationOutputProcessingFunction'], but we'll do it inside subsetOfParameterSpaceWrapper. So during parameter estimation there will be no separate call to a simulation output processing function.
        #Now start reducing various inputs...
        reducedIndices = UserInput.model['reducedParameterSpace']
        UserInput.InputParameterInitialGuess = returnReducedIterable(UserInput.InputParameterInitialGuess, reducedIndices)
        UserInput.parameterNamesList = returnReducedIterable(UserInput.parameterNamesList, reducedIndices)
        #We need to reparse to populate UserInput.stringOfParameterNames, can't use return Reduced Iterable.
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        #To make UserInput.parameterNamesAndMathTypeExpressionsDict we use a for loop to remove keys that should not be there anymore.
        #need to trim the dictionary based on what is in the UserInput.parameterNamesList variable
        parameterNamesAndMathTypeExpressionsDict = copy.deepcopy(self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'])
        for keyIndex in range(len(parameterNamesAndMathTypeExpressionsDict)):
            key = list(self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'])[keyIndex] #Need to call it out separately from original dictionary due to loop making the new dictionary smaller.
            if key not in self.UserInput.parameterNamesList:
                del parameterNamesAndMathTypeExpressionsDict[key] #Remove any parameters that were not in reduced parameter space.
        UserInput.parameterNamesAndMathTypeExpressionsDict = parameterNamesAndMathTypeExpressionsDict
        UserInput.InputParametersPriorValuesUncertainties = returnReducedIterable(UserInput.InputParametersPriorValuesUncertainties, reducedIndices)
        UserInput.std_prior     = returnReducedIterable( UserInput.std_prior    , reducedIndices )
        UserInput.var_prior     = returnReducedIterable( UserInput.var_prior   , reducedIndices  )
        UserInput.covmat_prior     = returnReducedIterable( UserInput.covmat_prior    , reducedIndices )
        self.UserInput.scaling_uncertainties     = returnReducedIterable( self.UserInput.scaling_uncertainties    , reducedIndices )
        self.UserInput.mu_prior     = returnReducedIterable( self.UserInput.mu_prior    , reducedIndices )
        self.UserInput.mu_prior_scaled     = returnReducedIterable( self.UserInput.mu_prior_scaled    , reducedIndices )
        self.UserInput.var_prior_scaled     = returnReducedIterable( self.UserInput.var_prior_scaled    , reducedIndices )
        self.UserInput.covmat_prior_scaled     = returnReducedIterable( self.UserInput.covmat_prior_scaled    , reducedIndices )
        self.Q_mu     = returnReducedIterable( self.Q_mu    , reducedIndices )
        self.Q_covmat     = returnReducedIterable( self.Q_covmat    , reducedIndices )
        #There are no returns. Everything above is an implied return.
        return

    def get_responses_simulation_uncertainties(self, discreteParameterVector): #FIXME: Make sure this works with responses['reducedResponseSpace']  and model['reducedParameterSpace']. I don't think it does.
        if type(np.array(self.UserInput.model['responses_simulation_uncertainties'])) == type(np.array([0])): #If it's an array, we take it as is.
            responses_simulation_uncertainties = np.array(self.UserInput.model['responses_simulation_uncertainties'])*1.0
        else:  #Else we assume it's a function taking the discreteParameterVector.
            responses_simulation_uncertainties = self.UserInput.model['responses_simulation_uncertainties'](discreteParameterVector) #This is passing an argument to a function.
        return responses_simulation_uncertainties

    def simulateWithSubsetOfParameters(self,reducedParametersVector): #This is a wrapper.
        #This function has implied arguments of ...
        #self.UserInput.model['InputParameterInitialGuess'] for the parameters to start with
        #self.UserInput.model['reducedParameterSpace'] a list of indices for which parameters are the only ones to change.
        #simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        #simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']
        #When this wrapper is used, EVERYWHERE ELSE will call it to do the simulation, by calling self.UserInput.simulationFunction and self.UserInput.simulationOutputProcessingFunction
        simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction'] #This is making a local simulation function. The global will be set ot simulateWithSubsetOfParameters.
        simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction'] #This is making a local simulation function. The global will be set to None.
        
        #now populate the discreteParameterVector first with the initial guess, then with the new reducedParameters vector.
        discreteParameterVector = copy.deepcopy(self.UserInput.model['InputParameterInitialGuess']) #This is the original one from the user, before any reduction.
        for reducedParameterIndex, parameterValue in enumerate(reducedParametersVector):
            #we find which index to put things into from #self.UserInput.model['reducedParameterSpace'], which is a list of indices.
            regularParameterIndex = self.UserInput.model['reducedParameterSpace'][reducedParameterIndex]
            discreteParameterVector[regularParameterIndex] = parameterValue
        if type(simulationFunction) != type(None):#This is the normal case.
            simulationOutput = simulationFunction(discreteParameterVector) 
        elif type(simulationOutput) == type(None):
            return 0, None #This is for the case that the simulation fails. User can have simulationOutput return a None type in case of failure. Perhaps should be made better in future. 

            
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        
        simulatedResponses = np.atleast_2d(simulatedResponses)
        #This is not needed:
        #observedResponses = np.atleast_2d(self.UserInput.responses_observed)
        return simulatedResponses
    
    def transform_responses(self, nestedAllResponsesArray, nestedAllResponsesUncertainties = []):
        nestedAllResponsesArray_transformed = copy.deepcopy(nestedAllResponsesArray) #First make a copy to populate with transformed values.
        nestedAllResponsesUncertainties_transformed = copy.deepcopy(nestedAllResponsesUncertainties) #First make a copy to populate with transformed values. If blank, we won't populate it.        
        UserInput = self.UserInput
        
        #TODO: Make little function for interpolation in case it's necessary (see below).
#        def littleInterpolator():
#            abscissaRange = UserInput.responses_abscissa[responseIndex][-1] - UserInput.responses_abscissa[responseIndex][0] #Last value minus first value.
#            UserInput.responses_observed = np.atleast_2d(UserInput.responses_observed)
#            UserInput.responses_observed_uncertainties = np.atleast_2d(UserInput.responses_observed_uncertainties)
        if 'data_overcategory' not in UserInput.responses:  #To make backwards compatibility.
            UserInput.responses['data_overcategory'] = ''
        if UserInput.responses['data_overcategory'] == 'transient_kinetics': #This assumes that the abscissa is always time.
            for responseIndex, response in enumerate(UserInput.responses_observed):
                #We will need the abscissa also, so need to check if there are independent abscissa or not:
                if len(UserInput.responses_abscissa) == 1: #This means there is only one abscissa.
                    abscissaIndex = 0
                else:
                    abscissaIndex = responseIndex
                #Now to do the transforms.
                if UserInput.responses['response_types'][responseIndex] == 'I':	 #For intermediate
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        #Perform the littleEuler twice.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray_transformed[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties_transformed[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'R':	#For reactant
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        LittleEuler
                if UserInput.responses['response_types'][responseIndex] == 'P':	 #For product
                    
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        #TODO: use responses['points_if_transformed'] variable to interpolate the right number of points. This is for data that's not already evenly spaced.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'O':
                    if UserInput.responses['response_data_type'][responseIndex] == 'o':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        LittleEuler
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        LittleEulerTwice
        if UserInput.responses['data_overcategory'] == 'steady_state_kinetics': #TODO: so far, this does not do anything. It assumes that the abscissa is never time.
            for responseIndex, response in enumerate(UserInput.responses_observed):
                if UserInput.responses['response_types'][responseIndex] == 'T':	 #For abscissa of temperature dependence. Will probably do a log transform.
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        pass
                if UserInput.responses['response_types'][responseIndex] == 'c':	 #For abscissa of concentration dependence.
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        pass
        return nestedAllResponsesArray_transformed, nestedAllResponsesUncertainties_transformed  
  
    #This helper function has been made so that gridSearch and design of experiments can call it.
    #Although at first glance it may seem like it should be in the gridCombinations module, that is a misconception. This is just a wrapper setting defaults for calling that module, such as using the prior for the grid interval when none is provided.
    #note that a blank list is okay for gridSamplingAbsoluteIntervalSize if doing a parameter grid, but not for other types of grids.
    def getGridCombinations(self, gridCenterVector, gridSamplingAbsoluteIntervalSize, gridSamplingNumOfIntervals, SpreadType="Addition",toFile=False):
        import CheKiPEUQ.CombinationGeneratorModule as CombinationGeneratorModule
        numParameters = len(gridCenterVector)
        if len(gridSamplingNumOfIntervals) == 0:
            gridSamplingNumOfIntervals = np.ones(numParameters, dtype='int') #By default, will make ones.
            numGridPoints = 3**numParameters
        else: 
            gridSamplingNumOfIntervals = np.array(gridSamplingNumOfIntervals, dtype='int')
            numGridPoints = 1 #just initializing.
            for radius in gridSamplingNumOfIntervals:
                numGridPoints=numGridPoints*(2*radius+1)
        if len(gridSamplingAbsoluteIntervalSize) == 0:
            gridSamplingAbsoluteIntervalSize = self.UserInput.std_prior #By default, we use the standard deviations associated with the priors.
        else: gridSamplingAbsoluteIntervalSize = np.array(gridSamplingAbsoluteIntervalSize, dtype='float')
        gridCombinations = CombinationGeneratorModule.combinationGenerator(gridCenterVector, gridSamplingAbsoluteIntervalSize, gridSamplingNumOfIntervals, SpreadType=SpreadType,toFile=toFile)
        return gridCombinations, numGridPoints  
  
    @CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def doGridSearch(self, searchType='getLogP', exportLog = True, gridSamplingAbsoluteIntervalSize = [], gridSamplingNumOfIntervals = [], passThroughArgs = {}):
        # gridSamplingNumOfIntervals is the number of variations to check in units of variance for each parameter. Can be 0 if you don't want to vary a particular parameter in the grid search.
        #TODO: the upper part of the gridsearch may not be compatibile with reduced parameter space. Needs to be checked.
        verbose = self.UserInput.parameter_estimation_settings['verbose']

        gridCenter = self.UserInput.InputParameterInitialGuess #We take what is in the variable self.UserInput.InputParameterInitialGuess for the center of the grid.    
        gridCombinations, numGridPoints = self.getGridCombinations(gridCenter, gridSamplingAbsoluteIntervalSize, gridSamplingNumOfIntervals)
        allGridResults = []
        #Initialize some things before loop.
        if (type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None)) or (verbose == True):
                import timeit
                timeAtGridStart = timeit.time.clock()
                timeAtLastGridPoint = timeAtGridStart #just initializing
        highest_logP = float('-inf') #Just initializing.
        #Start grid search loop.
        for combinationIndex,combination in enumerate(gridCombinations):
            self.UserInput.InputParameterInitialGuess = combination #We need to fill the variable InputParameterInitialGuess with the combination being checked.
            if searchType == 'getLogP':
                self.map_logP = self.getLogP(combination) #The getLogP function does not fill map_logP by itself.
                self.map_parameter_set = combination
                thisResult = [self.map_logP, str(self.map_parameter_set).replace(",","|").replace("[","").replace('(','').replace(')',''), 'None', 'None', 'None', 'None', 'None', 'None']
            if searchType == 'doMetropolisHastings':
                thisResult = self.doMetropolisHastings()
                #self.map_logP gets done by itself in Metrpolos hastings.                
            if searchType == 'doOptimizeNegLogP':
                thisResult = self.doOptimizeNegLogP(**passThroughArgs)
                #FIXME: the column headings of "thisResult" are wrong for the case of doOptimizeNegLogP.
                #What we really need to do is have the log file's column headings generated based on the searchType.
            if searchType == 'doOptimizeSSR':
                thisResult = self.doOptimizeSSR(**passThroughArgs)
            if (type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None)) or (verbose == True):
                timeAtThisGridPoint = timeit.time.clock()
                timeOfThisGridPoint = timeAtThisGridPoint - timeAtLastGridPoint
                averageTimePerGridPoint = (timeAtThisGridPoint - timeAtGridStart)/(combinationIndex+1)
                numRemainingGridPoints = numGridPoints - combinationIndex+1
                timeAtLastGridPoint = timeAtThisGridPoint #Updating.
            if self.map_logP > highest_logP: #This is the grid point in space with the highest value found so far and will be kept.
                bestResultSoFar = thisResult
                highest_logP = self.map_logP
                highest_logP_parameter_set = self.map_parameter_set
            allGridResults.append(thisResult)
            if verbose == True:
                print("GridPoint", combination, "number", combinationIndex+1, "out of", numGridPoints, "timeOfThisGridPoint", timeOfThisGridPoint)
                print("GridPoint", combinationIndex+1, "averageTimePerGridPoint", "%.2f" % round(averageTimePerGridPoint,2), "estimated time remaining", "%.2f" % round( numRemainingGridPoints*averageTimePerGridPoint,2), "s" )
                print("GridPoint", combinationIndex+1, "current logP", self.map_logP, "highest logP", highest_logP)
            elif type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None): #If verbose off but checkpoint frequency is on.
                if (combinationIndex/self.UserInput.parameter_estimation_settings['checkPointFrequency']).is_integer():
                    print("GridPoint", combination, "number", combinationIndex+1, "out of", numGridPoints, "timeOfThisGridPoint", timeOfThisGridPoint)
                    print("GridPoint", combinationIndex+1, "averageTimePerGridPoint", "%.2f" % round(averageTimePerGridPoint,2), "estimated time remaining", "%.2f" % round( numRemainingGridPoints*averageTimePerGridPoint,2), "s" )
                    print("GridPoint", combinationIndex+1, "current logP", self.map_logP, "highest logP", highest_logP)
            
        #TODO: export the allGridResults to file at end of search in a nicer format.        
        #First set the initial guess back to the center of the grid.
        self.UserInput.InputParameterInitialGuess = gridCenter
        #Now populate the map etc. with those of the best result.
        self.map_logP = highest_logP 
        self.map_parameter_set = highest_logP_parameter_set 
        if exportLog == True:
            with open("gridsearch_log_file.txt", 'w') as out_file:
                out_file.write("result: " + "self.map_logP, self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec" + "\n")
                for resultIndex, result in enumerate(allGridResults):
                    out_file.write("result: " + str(resultIndex) + " " +  str(result) + "\n")
        print("Final map results from gridsearch:", self.map_parameter_set, "final logP:", self.map_logP)
        if searchType == 'doMetropolisHastings':
            #Metropolis hastings has other variables to populate.
            #[self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] =
            return bestResultSoFar # [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] 
        if searchType == 'doOptimizeNegLogP':            
            return bestResultSoFar# [self.map_parameter_set, self.map_logP]
        if searchType == 'getLogP':          
            return bestResultSoFar# [self.map_parameter_set, self.map_logP]

    #The below function is a helper function that is used during doeInfoGainMatrix. However, it can certainly be used for other purposes.
    def populateResponsesWithSyntheticData(self, parModulationCombination):
        #For each parameter Modulation Combination we are going to obtain a matrix of info_gains that is based on a grid of the independent_variables.
        #First we need to make some synthetic data using parModulationCombination for the discreteParameterVector
        discreteParameterVector = parModulationCombination
        simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.  
        simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
        simulationOutput =simulationFunction(discreteParameterVector)
        if type(simulationOutput)==type(None):
            return float('-inf'), None #This is intended for the case that the simulation fails. User can return "None" for the simulation output. Perhaps should be made better in future.
        if np.array(simulationOutput).any()==float('nan'):
            return float('-inf'), None #This is intended for the case that the simulation fails without returning "None".
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        simulatedResponses = np.atleast_2d(simulatedResponses)
        #need to check if there are any 'responses_simulation_uncertainties'. #TODO: This isn't really implemented yet.
        if type(self.UserInput.model['responses_simulation_uncertainties']) == type(None): #if it's a None type, we keep it as a None type
            responses_simulation_uncertainties = None
        else:  #Else we get it based on the the discreteParameterVector
            responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(discreteParameterVector)
        
        synthetic_data  = simulatedResponses
        synthetic_data_uncertainties = responses_simulation_uncertainties
        #We need to populate the "observed" responses in userinput with the synthetic data.
        self.UserInput.responses['responses_observed'] = simulatedResponses
        self.UserInput.responses['responses_observed_uncertainties'] = simulatedResponses
        #Now need to do something unusual: Need to call the __init__ function again so that the arrays get reshaped as needed etc.
        self.__init__(self.UserInput)
    
    #This function requires first populating the doe_settings dictionary in UserInput in order to know which conditions to explore.
    software_name = "CheKiPEUQ Bayesian Design of Experiments"
    software_version = "1.0.2"
    software_unique_id = "https://doi.org/10.1002/cctc.202000976"
    software_kwargs = {"version": software_version, "author": ["Eric A. Walker", "Kishore Ravisankar", "Aditya Savara"], "doi": "https://doi.org/10.1002/cctc.202000976", "cite": "Eric Alan Walker, Kishore Ravisankar, Aditya Savara. CheKiPEUQ Intro 2: Harnessing Uncertainties from Data Sets, Bayesian Design of Experiments in Chemical Kinetics. ChemCatChem. Accepted. doi:10.1002/cctc.202000976"} 
    @CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    @CiteSoft.module_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def doeGetInfoGainMatrix(self, parameterCombination):#Note: There is an implied argument of info_gains_matrices_array_format being 'xyz' or 'meshgrid'
        #At present, we *must* provide a parameterCombination because right now the only way to get an InfoGainMatrix is with synthetic data assuming a particular parameterCombination as the "real" or "actual" parameterCombination.
        doe_settings = self.UserInput.doe_settings
        self.middle_of_doe_flag = True  #This is a work around that is needed because right now the synthetic data creation has an __init__ call which is going to try to modify the independent variables back to their original values if we don't do this.
        info_gain_matrix = [] #Right now, if using KL_divergence, each item in here is a single array. It is a sum across all parameters. 
        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each':
            info_gain_matrices_each_parameter = [] #make a matrix ready to copy info_gain_matrix. 
            #need to make a list of lists (or similar) to fill it with the individual matrices necessary.
            numParameters = len(self.UserInput.InputParametersPriorValuesUncertainties)
            for parameterIndex in range(0,numParameters):#looping across number of parameters...
                info_gain_matrices_each_parameter.append([]) #These are empty lists create to indices and initialize each parameter's info_gain_matrix. They will be appended to later.
            self.info_gain_matrices_each_parameter = info_gain_matrices_each_parameter #Need to initialize this since it's nested so can't be initialized in a loop later.
        if self.UserInput.doe_settings['info_gains_matrices_array_format'] == 'xyz':
            self.info_gains_matrices_array_format = 'xyz'            
            #For the IndependentVariables the grid info must be defined ahead of time. On the fly conditions grid means it's generated again fresh for each parameter combination. (We are doing it this way out of convenience during the first programming of this feature).
            if doe_settings['on_the_fly_conditions_grids'] == True:
                conditionsGridCombinations, numGridPoints = self.getGridCombinations(doe_settings['independent_variable_grid_center'], doe_settings['independent_variable_grid_interval_size'], doe_settings['independent_variable_grid_num_intervals'])
            #Here is the loop across conditions.                
            for conditionsCombinationIndex,conditionsCombination in enumerate(conditionsGridCombinations):    
                #It is absolutely critical that we *do not* use syntax like self.UserInput.responses['independent_variables_values'] = xxxx
                #Because that would move where the pointer is going to. We need to instead populate the individual values in the simulation module's namespace.
                #This population Must occur here. It has to be after the indpendent variables have changed, before synthetic data is made, and before the MCMC is performed.
                self.UserInput.model['populateIndependentVariablesFunction'](conditionsCombination)
                self.populateResponsesWithSyntheticData(parameterCombination)
                [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doMetropolisHastings()
                conditionsCombination = np.array(conditionsCombination) #we're going to make this an array before adding to the info_gain matrix.
                conditionsCombinationAndInfoGain = np.hstack((conditionsCombination, info_gain))
                info_gain_matrix.append(conditionsCombinationAndInfoGain)
                if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above lines for the sum.
                    for parameterIndex in range(0,numParameters):#looping across number of parameters...
                        conditionsCombinationAndInfoGain = np.hstack((conditionsCombination, np.array(self.info_gain_each_parameter[parameterIndex]))) #Need to pull the info gain matrix from the nested objected named info_gain_each_parameter
                        #Below mimics the line above which reads info_gain_matrix.append(conditionsCombinationAndInfoGain)
                        info_gain_matrices_each_parameter[parameterIndex].append(conditionsCombinationAndInfoGain)
            self.info_gain_matrix = np.array(info_gain_matrix) #this is an implied return in addition to the real return.
            if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above line for the sum.
                for parameterIndex in range(0,numParameters):#looping across number of parameters...
                    self.info_gain_matrices_each_parameter[parameterIndex]= np.array(info_gain_matrices_each_parameter[parameterIndex])
            self.middle_of_doe_flag = False #Set this back to false once info gain matrix is ready.
            return np.array(info_gain_matrix)            
        if self.UserInput.doe_settings['info_gains_matrices_array_format'] == 'meshgrid':
            self.info_gains_matrices_array_format = 'meshgrid'  
            if len(doe_settings['independent_variable_grid_center']) !=2:
                print("CURRENTLY THE INFOGAIN MESHGRID OPTION IS ONLY SUPPORTED FOR TWO INDEPENDENT VARIABLES. Use doe_settings['independent_variable_grid_center'] = 'xyz' and run again.")
                sys.exit()
            #STEP 1 is just to append each info_gain matrix to info_gain_matrix, and step 2 is 
            #For loop to generate info_gains_matrix.
            #For the IndependentVariables the grid info must be defined ahead of time. On the fly conditions grid means it's generated again fresh for each parameter combination. (We are doing it this way out of convenience during the first programming of this feature).
            if doe_settings['on_the_fly_conditions_grids'] == True:
                independentVariable1CentralValue = doe_settings['independent_variable_grid_center'][0]
                independentVariable2CentralValue = doe_settings['independent_variable_grid_center'][1]
                independentVariable1UpperValue = independentVariable1CentralValue + doe_settings['independent_variable_grid_interval_size'][0]*doe_settings['independent_variable_grid_num_intervals'][0]
                independentVariable1LowerValue = independentVariable1CentralValue - doe_settings['independent_variable_grid_interval_size'][0]*doe_settings['independent_variable_grid_num_intervals'][0]
                independentVariable2UpperValue =  independentVariable2CentralValue + doe_settings['independent_variable_grid_interval_size'][1]*doe_settings['independent_variable_grid_num_intervals'][1]
                independentVariable2LowerValue =  independentVariable2CentralValue - doe_settings['independent_variable_grid_interval_size'][1]*doe_settings['independent_variable_grid_num_intervals'][1]
                independentVariable1ValuesArray = np.linspace(independentVariable1LowerValue,independentVariable1UpperValue,doe_settings['independent_variable_grid_num_intervals'][0]*2+1)
                independentVariable2ValuesArray = np.linspace(independentVariable2LowerValue,independentVariable2UpperValue,doe_settings['independent_variable_grid_num_intervals'][1]*2+1)
                self.meshGrid_independentVariable1ValuesArray = independentVariable1ValuesArray #This is sortof an implied return.
                self.meshGrid_independentVariable2ValuesArray = independentVariable2ValuesArray #This is sortof an implied return.
                #Here is the loop across conditions.
                for indValue2 in independentVariable2ValuesArray: #We know from experience that the outer loop should be over the YY variable.
                    for indValue1 in independentVariable1ValuesArray: #We know from experience that the inner loop should be over the XX variable.
                        #It is absolutely critical that we *do not* use syntax like self.UserInput.responses['independent_variables_values'] = xxxx
                        #Because that would move where the pointer is going to. We need to instead populate the individual values in the simulation module's namespace.
                        #This population Must occur here. It has to be after the indpendent variables have changed, before synthetic data is made, and before the MCMC is performed.
                        self.UserInput.model['populateIndependentVariablesFunction']([indValue1,indValue2])
                        self.populateResponsesWithSyntheticData(parameterCombination)
                        [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doMetropolisHastings()
                        conditionsCombination = np.array([indValue1,indValue2])
                        conditionsCombinationAndInfoGain = np.hstack((conditionsCombination, info_gain))
                        info_gain_matrix.append(conditionsCombinationAndInfoGain) #NOTE that the structure *includes* the combinations.
                        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above lines for the sum.
                            for parameterIndex in range(0,numParameters):#looping across number of parameters...
                                conditionsCombinationAndInfoGain = np.hstack((conditionsCombination, np.array(self.info_gain_each_parameter[parameterIndex]))) #Need to pull the info gain matrix from the nested objected named info_gain_each_parameter
                                #Below mimics the line above which reads info_gain_matrix.append(conditionsCombinationAndInfoGain)
                                info_gain_matrices_each_parameter[parameterIndex].append(conditionsCombinationAndInfoGain)
                self.info_gain_matrix = np.array(info_gain_matrix) #this is an implied return in addition to the real return.
                if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above line for the sum.
                    for parameterIndex in range(0,numParameters):#looping across number of parameters...
                        self.info_gain_matrices_each_parameter[parameterIndex]= np.array(info_gain_matrices_each_parameter[parameterIndex])
                self.middle_of_doe_flag = False #Set this back to false once info gain matrix is ready.
                return np.array(info_gain_matrix)
    
    #This function requires population of the UserInput doe_settings dictionary. It automatically scans many parameter modulation combinations.
    def doeParameterModulationCombinationsScanner(self):
        import CheKiPEUQ.CombinationGeneratorModule as CombinationGeneratorModule
        doe_settings = self.UserInput.doe_settings 
        #For the parameters, we are able to use a default one standard deviation grid if gridSamplingAbsoluteIntervalSize is a blank list.
        #doe_settings['parameter_modulation_grid_center'] #We do NOT create such a variable in user input. The initial guess variable is used, which is the center of the prior if no guess has been provided.
        parModulationGridCenterVector = self.UserInput.InputParameterInitialGuess
        numParameters = len(parModulationGridCenterVector)
        parModulationGridIntervalSizeAbsolute = doe_settings['parameter_modulation_grid_interval_size']*self.UserInput.std_prior
        parModulationGridCombinations, numGridPoints = self.getGridCombinations(parModulationGridCenterVector,parModulationGridIntervalSizeAbsolute, doe_settings['parameter_modulation_grid_num_intervals'])
        
        parModulationGridCombinations= np.array(parModulationGridCombinations)
        if len(self.UserInput.parameterNamesList) == len(self.UserInput.InputParametersPriorValuesUncertainties): #then we assume variable names have been provided.
            headerString = self.UserInput.stringOfParameterNames #This variable is a string, no brackets.
        else: #else no variable names have been provided.
            headerString = ''
        np.savetxt("Info_gain__parModulationGridCombinations.csv", parModulationGridCombinations, delimiter=",", encoding =None, header=headerString)
        
        #We will get a separate info gain matrix for each parModulationCombination, we'll store that in this variable.
        info_gains_matrices_list = []
        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #just making analogous structure which exists for sum.
            info_gains_matrices_lists_one_for_each_parameter = [] #make a matrix ready to copy info_gains_matrices_list. 
            #need to make a list of lists (or similar) to fill it with the individual matrices necessary.
            numParameters = len(self.UserInput.InputParametersPriorValuesUncertainties)
            for parameterIndex in range(0,numParameters):#looping across number of parameters...
                info_gains_matrices_lists_one_for_each_parameter.append([]) #These are empty lists create to indices and initialize each parameter's info_gain_matrix. They will be appended to later.
        for parModulationCombinationIndex,parModulationCombination in enumerate(parModulationGridCombinations):                
            #We will get separate info gain matrix for each parameter modulation combination.
            info_gain_matrix = self.doeGetInfoGainMatrix(parModulationCombination)
            #Append the info gain matrix obtainend.
            info_gains_matrices_list.append(np.array(info_gain_matrix))
            if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above lines which were for the sum.
                    for parameterIndex in range(0,numParameters):#looping across number of parameters...
                        info_gains_matrices_lists_one_for_each_parameter[parameterIndex].append(np.array(self.info_gain_matrices_each_parameter[parameterIndex]))
        self.info_gains_matrices_array=np.array(info_gains_matrices_list) #This is an implied return, but we will also return it.
        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above line for the sum.
            self.info_gains_matrices_arrays_one_for_each_parameter = list(self.UserInput.InputParametersPriorValuesUncertainties) #initializing it with right length, then will fill it.
            for parameterIndex in range(0,numParameters):#looping across number of parameters...
                self.info_gains_matrices_arrays_one_for_each_parameter[parameterIndex]= np.array(info_gains_matrices_lists_one_for_each_parameter[parameterIndex]) #make each an array like above.
            self.info_gains_matrices_arrays_one_for_each_parameter = np.array(self.info_gains_matrices_arrays_one_for_each_parameter)
        #TODO: write the self.info_gains_matrices_array individual elements to file.
        #for modulationIndex in range(len(self.info_gains_matrices_array)):
            #self.info_gains_matrices_array[modulationIndex]  #Write this to file. This is 'xyz' format regardless of whether self.info_gains_matrices_array_format == 'xyz'  or =='meshgrid' is used.
        return self.info_gains_matrices_array
    
    def createInfoGainPlots(self, parameterIndices=[], plot_suffix = ''):
        #parameterIndices should be a list of parameters if the user only wants as subset of parameters. The default, a blank list, will do all if the setting for doing each is on.
        
        #first make the modulation plots for the Sum.
        self.createInfoGainModulationPlots(parameterIndex=None, plot_suffix = plot_suffix)
        #now, by default, loop through and make plots fore each parameterIndex if the setting for that is on.
        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each':
            if len(parameterIndices) > 0: #if the user has provided a list of parameters, we will only make the plots for those parameters.
               for parameterIndex in parameterIndices:
                    plotSuffixString = "_par_" + str(parameterIndex) + plot_suffix
                    self.createInfoGainModulationPlots(parameterIndex=parameterIndex, plot_suffix = plotSuffixString)
            if len(parameterIndices) == 0: #This is the default case, and we'll make plots for each parameter.
                numParameters = len(self.UserInput.InputParametersPriorValuesUncertainties)
                for parameterIndex in range(0,numParameters):
                    plotSuffixString = "_par_" + str(parameterIndex) + plot_suffix
                    self.createInfoGainModulationPlots(parameterIndex=parameterIndex, plot_suffix = plotSuffixString)
    
    def createInfoGainModulationPlots(self, parameterIndex=None, plot_suffix = ''): 
        #Right now, when using KL_divergence and design of experiments there is an option of UserInput.doe_settings['info_gains_matrices_multiple_parameters'] = 'each' or 'sum'
        #the default is sum. But when it is 'each', then it is possible to plot separate info_gains for each parameter.
        #Note: the below code *does not* add a suffix to inidicate when a parameter Index has been fed.
        #TODO: The variable "parameterInfoGainIndex" is made with the presumption that later we'll have to add another index when we have info_gains for each parameter. In that case it will become like this:
        #xValues = self.info_gains_matrices_array[modulationIndex][:,0] will become xValues = self.info_gains_matrices_array[modulationIndex][parameterInfoGainIndex][:,0]
        #self.meshGrid_independentVariable1ValuesArray will remain unchanged.      
        
        import CheKiPEUQ.plotting_functions as plotting_functions
        
        #assess whether the function is called for the overall info_gain matrices or for a particular parameter.
        if parameterIndex==None:  #this means we're using the regular info gain, not the parameter specific case.
            #Normally, the info gain plots should be stored in self.info_gains_matrices_array.
            #However, in case it does not exist or there are none in there, then we assume the person is trying to make just one. So we take the most recent info gain matrix.
            try:
                if len(self.info_gains_matrices_array) >= 0: #normally, it should exist and be populated.
                    pass
                if len(self.info_gains_matrices_array) == 0:#in case it exists but is not populated, we'll populated.
                    self.info_gains_matrices_array = np.array([self.info_gain_matrix])
            except: #if it does not yet exist, we create it and populate it.
                    self.info_gains_matrices_array = np.array([self.info_gain_matrix])
            local_info_gains_matrices_array = self.info_gains_matrices_array #We have to switch to a local variable since that way below we can use the local variable whether we're doing the 'global' info_gains_matrices array or a parameter specific one.
        if parameterIndex!=None:
            if hasattr(self, 'info_gains_matrices_arrays_one_for_each_parameter'): #this structure will only exist if doeParameterModulationCombinationsScanner has been called.
                local_info_gains_matrices_array = np.array(self.info_gains_matrices_arrays_one_for_each_parameter)[:][parameterIndex] #each "row" is a modulation, and within that are structures for each parameter.  This is further described in the document InfoGainMatrixObjectsStructure.docx
            else: #if a modulation has not been run, and simply doeGetInfoGainMatrix was done, then the larger structure might not exist and we have to just pull out by the parameter index and then make it nested as for a regular info_gain sum.
                local_info_gains_matrices_array = np.array([self.info_gain_matrices_each_parameter[parameterIndex]])
        #At present, plots are only made if the number of independent variables is 2.
        if len(self.UserInput.doe_settings['independent_variable_grid_center']) == 2:
            if self.info_gains_matrices_array_format == 'xyz':                
                for modulationIndex in range(len(local_info_gains_matrices_array)):
                    xValues = local_info_gains_matrices_array[modulationIndex][:,0]
                    yValues = local_info_gains_matrices_array[modulationIndex][:,1]
                    zValues = local_info_gains_matrices_array[modulationIndex][:,2]
                    plotting_functions.makeTrisurfacePlot(xValues, yValues, zValues, figure_name = "Info_gain_TrisurfacePlot_modulation_"+str(modulationIndex)+plot_suffix)
            if self.info_gains_matrices_array_format == 'meshgrid':        
                for modulationIndex in range(len(local_info_gains_matrices_array)):
                    #Now need to get things prepared for the meshgrid.
                    #NOTE: we do not pull XX and YY from local_info_gains_matrices_array because that is 1D and these are 2D arrays made a different way.
                    #xValues = local_info_gains_matrices_array[modulationIndex][:,0] #Still correct, but not being used.
                    #yValues = local_info_gains_matrices_array[modulationIndex][:,1] #Still correct, but not being used.
                    XX, YY = np.meshgrid(self.meshGrid_independentVariable1ValuesArray, self.meshGrid_independentVariable2ValuesArray)
                    zValues = local_info_gains_matrices_array[modulationIndex][:,2]
                    ZZ = zValues.reshape(XX.shape) #We know from experience to reshape this way.
                    plotting_functions.makeMeshGridSurfacePlot(XX, YY, ZZ, figure_name = "Info_gain_Meshgrid_modulation_"+str(modulationIndex)+plot_suffix)
        else:
            print("At present, createInfoGainPlots and createInfoGainModulationPlots only create plots when the length of  independent_variable_grid_center is 2. We don't currently support creation of other dimensional plots.")
    def getLogP(self, proposal_sample): #The proposal sample is specific parameter vector.
        [log_likelihood_proposal, simulationOutput_proposal] = self.getLogLikelihood(proposal_sample)
        log_prior_proposal = self.getLogPrior(proposal_sample)
        log_numerator_or_denominator = log_likelihood_proposal+log_prior_proposal #Of the Metropolis-Hastings accept/reject ratio
        return log_numerator_or_denominator
        
    def getNegLogP(self, proposal_sample): #The proposal sample is specific parameter vector. We are using negative of log P because scipy optimize doesn't do maximizing. It's recommended minimize the negative in this situation.
        neg_log_postererior = -1*self.getLogP(proposal_sample)
        return neg_log_postererior

    def doOptimizeNegLogP(self, simulationFunctionAdditionalArgs = (), method = None, optimizationAdditionalArgs = {}, printOptimum = True, verbose=True, maxiter=0):
        #THe intention of the optional arguments is to pass them into the scipy.optimize.minimize function.
        # the 'method' argument is for Nelder-Mead, BFGS, SLSQP etc. https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        #Note that "maxiter=0" just means to use the default.
        initialGuess = self.UserInput.InputParameterInitialGuess
        import scipy.optimize
        if verbose == False:
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(self.getNegLogP, initialGuess, method = method)
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(self.getNegLogP, initialGuess, method = method, options={"maxiter": maxiter})
        if verbose == True:
            verbose_simulator = verbose_optimization_wrapper(self.getNegLogP)
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"disp": True})
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"maxiter": maxiter})
            #print(f"Number of calls to Simulator instance {verbose_simulator.num_calls}") <-- this is the same as the "Function evaluations" field that gets printed.
            
        self.map_parameter_set = optimizeResult.x #This is the map location.
        self.map_logP = -1.0*optimizeResult.fun #This is the map logP
        if printOptimum == True:
            print("Final results from doOptimizeNegLogP:", self.map_parameter_set, "final logP:", self.map_logP)
        return [self.map_parameter_set, self.map_logP]


    def getSSR(self, discreteParameterVector): #The proposal sample is specific parameter vector. 
        #First do a parameter bounds check. We'll return an inf if it fails.
        passedBoundsCheck = self.doInputParameterBoundsChecks(discreteParameterVector)
        if passedBoundsCheck == False:
            return float('inf')
        
        #If within bounds, proceed to get the simulated responses.
        simulatedResponses = self.getSimulatedResponses(discreteParameterVector)
        if type(simulatedResponses) == type(None):
            return float('inf') #This is intended for the case that the simulation fails, indicated by receiving an 'nan' or None type from user's simulation function.
        
        #now calculate the SSR if nothing has failed.
        Residuals = np.array(simulatedResponses) - np.array(self.UserInput.responses_observed)
        SSR = np.sum(Residuals**2)
        return SSR

    def doOptimizeSSR(self, simulationFunctionAdditionalArgs = (), method = None, optimizationAdditionalArgs = {}, printOptimum = True, verbose=True, maxiter=0):
        #THe intention of the optional arguments is to pass them into the scipy.optimize.minimize function.
        # the 'method' argument is for Nelder-Mead, BFGS, SLSQP etc. https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        #Note that "maxiter=0" just means to use the default.
        initialGuess = self.UserInput.InputParameterInitialGuess
        import scipy.optimize
        if verbose == False:
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(self.getSSR, initialGuess, method = method)
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(self.getSSR, initialGuess, method = method, options={"maxiter": maxiter})
        if verbose == True:
            verbose_simulator = verbose_optimization_wrapper(self.getSSR)
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"disp": True})
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"maxiter": maxiter})
            #print(f"Number of calls to Simulator instance {verbose_simulator.num_calls}") <-- this is the same as the "Function evaluations" field that gets printed.
            
        self.opt_parameter_set = optimizeResult.x #This is the best fit parameter set.
        self.opt_SSR = optimizeResult.fun #This is the best fit SSR.
        if printOptimum == True:
            print("Final results from doOptimizeSSR:", self.opt_parameter_set, "final SSR:", self.opt_SSR)
        #FIXME: Right now, the createAllPlots command will not work unless we populate the map parameter set, so that is what we are doing. But a better longterm solution needs to be made. In which the graph says "opt" rather than "MAP" and uses the appropriate variables.
        #TODO: Also need to add things like WSSR based on magnitude and variance weightings.
        self.map_parameter_set = self.opt_parameter_set
        return [self.opt_parameter_set, self.opt_SSR]

    
    #This function is meant to be called from the runfile when testing a new function etc. It allows a simulation plot to be created.
    #This is *not* recommended for use in other functions, where it is recommended that getLogP be called directly.
    def doSinglePoint(self, discreteParameterVector=None, objectiveFunction='logP'):
        #objectiveFunction can be 'logP' or 'SSR'
        if type(discreteParameterVector)==type(None): #If somebody did not feed a specific vector, we take the initial guess.
            discreteParameterVector = self.UserInput.InputParameterInitialGuess
        if objectiveFunction=='logP':
            self.map_parameter_set = discreteParameterVector
            self.map_logP = self.getLogP(discreteParameterVector)
            objectiveFunctionValue = self.map_logP
        if objectiveFunction=='SSR':
            self.opt_parameter_set = discreteParameterVector
            self.opt_SSR = self.getSSR(discreteParameterVector)
            objectiveFunctionValue = self.opt_SSR
        return [discreteParameterVector, objectiveFunctionValue]
    
    #main function to get samples #TODO: Maybe Should return map_log_P and mu_AP_log_P?
    @CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def doMetropolisHastings(self):
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if type(self.UserInput.parameter_estimation_settings['mcmc_random_seed']) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        samples_simulatedOutputs = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],self.UserInput.num_data_points)) #TODO: Consider moving this out of this function.
        samples = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],len(self.UserInput.mu_prior)))
        mcmc_step_modulation_history = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: Make this optional for efficiency. #This allows the steps to be larger or smaller. Make this same length as samples. In future, should probably be same in other dimension also, but that would require 2D sampling with each step.                                                                          
        samples[0,:]=self.UserInput.InputParameterInitialGuess  # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean if an initial guess is not provided.. Do not use self.UserInput.model['InputParameterInitialGuess']  because that doesn't work with reduced parameter space feature.
        samples_drawn = samples*1.0 #this includes points that were rejected. #TODO: make this optional for efficiency.               
        log_likelihoods_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        log_posteriors_un_normed_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        log_postereriors_drawn = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: make this optional for efficiency. We don't want this to be 2D, so we don't copy log_posteriors_un_normed_vec.
        log_priors_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        #Code to initialize checkpoints.
        if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
            print("Starting MCMC sampling.")
            import timeit
            timeOfFirstCheckpoint = timeit.time.clock()
            timeCheckpoint = timeit.time.clock() - timeOfFirstCheckpoint #First checkpoint at time 0.
            numCheckPoints = self.UserInput.parameter_estimation_settings['mcmc_length']/self.UserInput.parameter_estimation_settings['checkPointFrequency']
        #Before sampling should fill in the first entry for the posterior vector we have created. #FIXME: It would probably be better to start with i of 0 in below sampling loop. I believe that right now the "burn in" and "samples" arrays are actually off by an index of 1. But trying to change that alters their length relative to other arrays and causes problems. Since we always do many samples and this only affects the initial point being averaged in twice, it is not a major problem. It's also avoided if people use a burn in of at least 1.
        log_posteriors_un_normed_vec[0]= self.getLogP(samples[0])
        for i in range(1, self.UserInput.parameter_estimation_settings['mcmc_length']): #FIXME: Don't we need to start with i of 0?
            sampleNumber = i #This is so that later we can change it to i+1 if the loop starts from i of 0 in the future.
            if self.UserInput.parameter_estimation_settings['verbose']: print("MCMC sample number", sampleNumber)                  
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'unbiased':
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_covmat*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'MAP_finding':
                if i == 1: mcmc_step_dynamic_coefficient = 1
                mcmc_step_modulation_coefficient = np.random.uniform() + 0.5 #TODO: make this a 2D array. One for each parameter.
                mcmc_step_modulation_history[i] = mcmc_step_modulation_coefficient
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_covmat*mcmc_step_dynamic_coefficient*mcmc_step_modulation_coefficient*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
            log_prior_proposal = self.getLogPrior(proposal_sample)
            [log_likelihood_proposal, simulationOutput_proposal] = self.getLogLikelihood(proposal_sample)
            log_prior_current_location = self.getLogPrior(samples[i-1,:]) #"current" location is the most recent accepted location, because we haven't decided yet if we're going to move.
            [log_likelihood_current_location, simulationOutput_current_location] = self.getLogLikelihood(samples[i-1,:]) #FIXME: the previous likelihood should be stored so that it doesn't need to be calculated again.
            log_accept_probability = (log_likelihood_proposal + log_prior_proposal) - (log_likelihood_current_location + log_prior_current_location) 
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current log_likelihood',log_likelihood_current_location, 'Proposed log_likelihood', log_likelihood_proposal, '\nLog of Accept_probability (gauranteed if above 0)', log_accept_probability)
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current posterior',log_likelihood_current_location+log_prior_current_location, 'Proposed Posterior', log_likelihood_proposal+log_prior_proposal)
            if self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'] != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'])
                #Our logP are of the type e^logP = P. #This is base 'e' because the logpdf functions are base e. Ashi checked the sourcecode.
                #The flattening code works in part because P is always < 1, so logP is always negative. 1/N_flatten at front brings negative number closer to zero which is P closer to 1. If logP is already positive, it will stay positive which also causes no problem.
                #TODO: add code that unflattens the final histograms, that way even with more sampling we still get an accurate final posterior distribution. We can also then add a flag if the person wants to keep the posterior flattened.
                log_accept_probability = (1/N_flatten)*log_accept_probability
            randomNumber = np.random.uniform()
            log_randomNumber = np.log(randomNumber) #This is base 'e' because the logpdf functions are base e. Ashi checked the sourcecode.
            if log_accept_probability > log_randomNumber:  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                if self.UserInput.parameter_estimation_settings['verbose']:
                  print('accept', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_proposal)
                samples[i,:] = proposal_sample
                samples_drawn[i,:] = proposal_sample
                log_postereriors_drawn[i] = (log_likelihood_proposal+log_prior_proposal) #FIXME: should be using getlogP
                samples_simulatedOutputs[i,:] = simulationOutput_proposal
                log_posteriors_un_normed_vec[i] = log_likelihood_proposal+log_prior_proposal 
                log_likelihoods_vec[i] = log_likelihood_proposal
                log_priors_vec[i] = log_prior_proposal
            else:
                if self.UserInput.parameter_estimation_settings['verbose']:
                  print('reject', proposal_sample)
                  sys.stdout.flush()
                  #print(simulationOutput_current_location)
                samples[i,:] = samples[i-1,:] #the sample is not kept if it is rejected, though we still store it in the samples_drawn.
                samples_drawn[i,:] = proposal_sample
                log_postereriors_drawn[i] = (log_likelihood_proposal+log_prior_proposal)
                samples_simulatedOutputs[i,:] = simulationOutput_current_location
                log_posteriors_un_normed_vec[i] = log_likelihood_current_location+log_prior_current_location
                log_likelihoods_vec[i] = log_likelihood_current_location
                log_priors_vec[i] = log_prior_current_location
            if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
                if sampleNumber%self.UserInput.parameter_estimation_settings['checkPointFrequency'] == 0: #The % is a modulus function.
                    timeSinceLastCheckPoint = (timeit.time.clock() - timeOfFirstCheckpoint) -  timeCheckpoint
                    timeCheckpoint = timeit.time.clock() - timeOfFirstCheckpoint
                    checkPointNumber = sampleNumber/self.UserInput.parameter_estimation_settings['checkPointFrequency']
                    averagetimePerSampling = timeCheckpoint/(sampleNumber)
                    print("MCMC sample number ", sampleNumber, "checkpoint", checkPointNumber, "out of", numCheckPoints) 
                    print("averagetimePerSampling", averagetimePerSampling, "seconds")
                    print("timeSinceLastCheckPoint", timeSinceLastCheckPoint, "seconds")
                    print("Estimated time remaining", averagetimePerSampling*(self.UserInput.parameter_estimation_settings['mcmc_length']-sampleNumber), "seconds")
                    if self.UserInput.parameter_estimation_settings['mcmc_mode'] != 'unbiased':
                        print("Most recent mcmc_step_dynamic_coefficient:", mcmc_step_dynamic_coefficient)
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] != 'unbiased':
                if sampleNumber%100== 0: #The % is a modulus function to change the modulation coefficient every n steps.
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
        self.burn_in_samples = samples[0:self.UserInput.parameter_estimation_settings['mcmc_burn_in']] #FIXME: this line will have to change with i indexing changed.
        self.post_burn_in_samples = samples[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:] #FIXME: this line will have to change with i indexing changed.
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
            self.post_burn_in_samples_simulatedOutputs = samples_simulatedOutputs[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_log_posteriors_un_normed_vec = log_posteriors_un_normed_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_logP_un_normed_vec = (self.post_burn_in_log_posteriors_un_normed_vec)
        self.post_burn_in_log_likelihoods_vec = log_likelihoods_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_log_priors_vec = log_priors_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        
        #BELOW CALCULATE INFOGAIN RELATED QUANTITIES.
        # posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        self.evidence = np.mean(np.exp(self.post_burn_in_log_posteriors_un_normed_vec))/np.linalg.norm(self.post_burn_in_samples)# another variety:*np.sqrt(2*np.pi*np.std(self.post_burn_in_samples)**2)                    
        if self.UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] == 'log_ratio':
            if self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'] == 0:        
                #we have log A, and we want log(A/B).  #log (e^log(A) / B )  = log(A/B).  
                #But we could also do...  log(A) - log(B) = log(A/B). So changing to that.
                post_burn_in_log_posteriors_vec = self.post_burn_in_log_posteriors_un_normed_vec - np.log(self.evidence) 
                log_ratios = (post_burn_in_log_posteriors_vec-self.post_burn_in_log_priors_vec) #log10(a/b) = log10(a)-log10(b)
                log_ratios[np.isinf(log_ratios)] = 0
                log_ratios = np.nan_to_num(log_ratios)
                self.info_gain_log_ratio_each_parameter = None #TODO: create a list or array of arrays such that the index is the parameter number.
                self.info_gain_log_ratio = np.mean(log_ratios) #NOTE: The log_ratio info_gain is *always* calculated, at this line or below.
            elif self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'] != 0:        
                #Need to consider using a truncated evidence array as well, but for now will not worry about that.
                #First intialize the stacked array.
                #Surprisingly, the arrays going in haves shapes like 900,1 rather than 1,900 so now transposing them before stacking.
                stackedLogProbabilities = np.vstack((self.post_burn_in_log_priors_vec.transpose(), self.post_burn_in_log_posteriors_un_normed_vec.transpose()))
                #Now, we are going to make a list of abscissaIndices to remove, recognizing that numpy arrays are "transposed" relative to excel.
                abscissaIndicesToRemove = [] 
                #FIXME: Below there are some "if verbose", but those should not be printed, they should be collected and exported to a file at the end.
                for abscissaIndex in range(np.shape(stackedLogProbabilities)[1]):
                    if self.UserInput.parameter_estimation_settings['verbose']:
                        print("parameter set:", self.post_burn_in_samples[abscissaIndex])
                    ordinateValues = stackedLogProbabilities[:,abscissaIndex]
                    #We mark anything where there is a 'nan':
                    if np.isnan( ordinateValues ).any(): #A working numpy syntax is to have the any outside of the parenthesis, for this command, even though it's a bit strange.
                        abscissaIndicesToRemove.append(abscissaIndex)
                        if self.UserInput.parameter_estimation_settings['verbose']:
                            print(abscissaIndex, "removed nan (log_prior, log_posterior)", ordinateValues, np.log( self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff']))
                    elif (ordinateValues < np.log( self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'] ) ).any(): #again, working numpy syntax is to put "any" on the outside. We take the log since we're looking at log of probability. This is a natural log.
                        abscissaIndicesToRemove.append(abscissaIndex)
                        if self.UserInput.parameter_estimation_settings['verbose']:
                            print(abscissaIndex, "removed small (prior, posterior)",  np.exp(ordinateValues), self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'])
                    else:
                        if self.UserInput.parameter_estimation_settings['verbose']:
                            print(abscissaIndex, "kept (prior, posterior)",  np.exp(ordinateValues), self.UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'])
                        pass
                #Now that this is finshed, we're going to do the truncation using numpy delete.
                stackedLogProbabilities_truncated = stackedLogProbabilities*1.0 #just initializing.
                stackedLogProbabilities_truncated = np.delete(stackedLogProbabilities, abscissaIndicesToRemove, axis=1)
                post_burn_in_log_priors_vec_truncated = stackedLogProbabilities_truncated[0]
                post_burn_in_log_posteriors_un_normed_vec_truncated = stackedLogProbabilities_truncated[1] #We have to truncate with not normalized, so we add the normalization in here.
                post_burn_in_log_posteriors_vec_truncated = np.log  ( np.exp( post_burn_in_log_posteriors_un_normed_vec_truncated) /self.evidence)
                #Now copy the same lines that Eric had used above, only change to using log_ratios_truncated
                log_ratios_truncated = (post_burn_in_log_posteriors_vec_truncated-post_burn_in_log_priors_vec_truncated)
                log_ratios_truncated[np.isinf(log_ratios_truncated)] = 0
                log_ratios_truncated = np.nan_to_num(log_ratios_truncated)
                self.info_gain_log_ratio_each_parameter = None #TODO: create a list or array of arrays such that the index is the parameter number.
                self.info_gain_log_ratio = np.mean(log_ratios_truncated) #NOTE: The log_ratio info_gain is *always* calculated, at this line or earlier. 
                #TODO: Export the below things.
                #post_burn_in_log_posteriors_vec_non_truncated = self.post_burn_in_log_posteriors_un_normed_vec - np.log(self.evidence)
                #print(post_burn_in_log_posteriors_vec_truncated) #TODO: Export this
                #print(post_burn_in_log_priors_vec_truncated)  #TODO: Export this
        if self.UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] == 'log_ratio':
            self.info_gain = self.info_gain_log_ratio
        if self.UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] == 'KL_divergence':
            #Below is the KL_divergence info_gain calculation.
            length, width = self.post_burn_in_samples.shape
            self.info_gain_KL = 0
            self.info_gain_KL_each_parameter  = []
            for param in range(width):
                (density0,bins0,pathces0)=plt.hist([self.samples_of_prior[:,param].flatten(),self.post_burn_in_samples[:,param].flatten()],bins=100,density=True)
                current_info_gain_KL = density0[1]*np.log(density0[1]/density0[0])
                current_info_gain_KL = current_info_gain_KL[np.isfinite(current_info_gain_KL)]
                current_info_gain_KL = np.sum(current_info_gain_KL)
                self.info_gain_KL_each_parameter.append(current_info_gain_KL) #could make this optional, but normally shouldn't take much memory.
                self.info_gain_KL = self.info_gain_KL + current_info_gain_KL
            self.info_gain_each_parameter = self.info_gain_KL_each_parameter #could make this optional, but normally shouldn't take much memory.
            self.info_gain = self.info_gain_KL
        
        #BELOW calculate MAP and mu_AP related quantities.
        map_logP = max(self.post_burn_in_logP_un_normed_vec)
        self.map_logP = map_logP
        self.map_index = list(self.post_burn_in_logP_un_normed_vec).index(map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
        self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.
        self.mu_AP_parameter_set = np.mean(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.
        self.stdap_parameter_set = np.std(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.
        #TODO: Probably should return the variance of each sample in the post_burn_in
        if self.UserInput.parameter_estimation_settings['verbose'] == True:
            print(self.map_parameter_set)
            print(self.mu_AP_parameter_set)
            print(self.stdap_parameter_set)
        if self.UserInput.parameter_estimation_settings['exportLog'] == True:
            #The self.post_burn_in_samples_simulatedOutputs has length of ALL sampling including burn_in
            #The self.post_burn_in_samples is only AFTER burn in.
            #The self.post_burn_in_log_posteriors_un_normed_vec is AFTER burn in.           
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
                out_file.write("self.stdap_parameter_set:" + str( self.stdap_parameter_set) + "\n")
                out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
                out_file.write("evidence:" + str(self.evidence) + "\n")
                out_file.write("posterior_cov_matrix:" + "\n" + str(np.cov(self.post_burn_in_samples.T)) + "\n")
                if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.var_prior).any() > 0.10:
                    pass #Disabling below warning until if statement is fixed.
                    #out_file.write("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% of prior variance in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632")
        if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.var_prior).any() > 0.10:  
            pass #Disabling below warning until if statement it is fixed.
            #print("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% of prior variance in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632  ")
        return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] # EAW 2020/01/08
    def getLogPrior(self,discreteParameterVector):
        if type(self.UserInput.model['custom_logPrior']) != type(None):
            logPrior = self.UserInput.model['custom_logPrior'](discreteParameterVector)
            return logPrior
        boundsChecksPassed = self.doInputParameterBoundsChecks(discreteParameterVector)
        if boundsChecksPassed == False: #If false, return a 'zero probability' type result. Else, continue getting log of prior..
            return float('-inf') #This approximates zero probability.        
        if self.UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "off":
            discreteParameterVector_scaled = np.array(discreteParameterVector)*1.0
        elif self.UserInput.parameter_estimation_settings['scaling_uncertainties_type'] != "off":
            if np.shape(self.UserInput.scaling_uncertainties)==np.shape(discreteParameterVector):
                discreteParameterVector_scaled = np.array(discreteParameterVector)/self.UserInput.scaling_uncertainties
            else: #TODO: If we're in the else statemnt, then the scaling uncertainties is a covariance matrix, for which we plan to do row and column scaling, which has not yet been implemented. #We could pobably just use the diagonal in the short term.
                print("WARNING: There is an error in your self.UserInput.scaling_uncertainties. Contact the developers with a bug report.")
                discreteParameterVector_scaled = np.array(discreteParameterVector)*1.0

        if hasattr(self.UserInput, 'InputParametersPriorValuesUniformDistributionsIndices') == False: #this is the normal case, no uniform distributionns.
            logPrior = multivariate_normal.logpdf(x=discreteParameterVector_scaled,mean=self.UserInput.mu_prior_scaled,cov=self.UserInput.covmat_prior_scaled)
        elif hasattr(self.UserInput, 'InputParametersPriorValuesUniformDistributionsIndices') == True: #This means that at least one variable has a uniform prior distribution. So we need to remove that  parameter before doing the multivariate_normal.logpdf.
            #Note that this if-statement is intentionally after the scaling uncertainties because that feature can be compatible with the uniform distribution.
            discreteParameterVector_scaled_truncated = np.delete(discreteParameterVector_scaled, self.UserInput.InputParametersPriorValuesUniformDistributionsIndices) #delete does not change original array.
            mu_prior_scaled_truncated = np.delete(self.UserInput.mu_prior_scaled, self.UserInput.InputParametersPriorValuesUniformDistributionsIndices) #delete does not change original array.
            var_prior_scaled_truncated = np.delete(self.UserInput.var_prior_scaled, self.UserInput.InputParametersPriorValuesUniformDistributionsIndices) #delete does not change original array.
            #Presently, we don't have full covmat support with uniform distributions. In principle, it would be better to use covmat_prior_scaled and delete the rows and columns since then we might have covmat support.
            #For now, we just make the truncated covmat from the var_prior. We currently don't have full covmat support for the case of uniform distributions.
            covmat_prior_scaled_truncated = np.diagflat(var_prior_scaled_truncated) 
            if len(covmat_prior_scaled_truncated) == 0: #if all variables are uniform, then need to return log(1) which is 0.
                logPrior = 0
            else:
                logPrior = multivariate_normal.logpdf(x=discreteParameterVector_scaled_truncated,mean=mu_prior_scaled_truncated,cov=covmat_prior_scaled_truncated)
        #Note: Below code should be okay regardless of whether there are uniform distributions since it only adjusts logPrior by a scalar.
        if self.UserInput.parameter_estimation_settings['undo_scaling_uncertainties_type'] == True:
            try:
                scaling_factor = float(self.UserInput.parameter_estimation_settings['scaling_uncertainties_type'])
                logPrior = logPrior - np.log(scaling_factor)
            except:
                if self.UserInput.parameter_estimation_settings['scaling_uncertainties_type'] != "off":
                    print("Warning: undo_scaling_uncertainties_type is set to True, but can only be used with a fixed value for scaling_uncertainties_type.  Skipping the undo.")
        return logPrior
        
    def doInputParameterBoundsChecks(self, discreteParameterVector): #Bounds are considered part of the prior, so are set in InputParameterPriorValues_upperBounds & InputParameterPriorValues_lowerBounds
        if len(self.UserInput.model['InputParameterPriorValues_upperBounds']) > 0:
            upperCheck = boundsCheck(discreteParameterVector, self.UserInput.model['InputParameterPriorValues_upperBounds'], 'upper')
            if upperCheck == False:
                return False
        if len(self.UserInput.model['InputParameterPriorValues_lowerBounds']) > 0:
            lowerCheck = boundsCheck(discreteParameterVector, self.UserInput.model['InputParameterPriorValues_lowerBounds'], 'lower')
            if lowerCheck == False:
                return False
        return True #If the test has gotten here without failing any of the tests, we return true.

    #This helper function must be used because it allows for the output processing function etc. It has been separated from getLogLikelihood so that it can be used by doOptimizeSSR etc.
    def getSimulatedResponses(self, discreteParameterVector): 
        simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.  
        simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
        simulationOutput =simulationFunction(discreteParameterVector) 
        if type(simulationOutput)==type(None):
            return None #This is intended for the case that the simulation fails. User can return "None" for the simulation output.
        if np.array(simulationOutput).any()==float('nan'):
            print("WARNING: Your simulation output returned a 'nan' for parameter values " +str(discreteParameterVector) + ". 'nan' values cannot be processed by the CheKiPEUQ software and this set of Parameter Values is being assigned a probability of 0.")
            return None #This is intended for the case that the simulation fails in some way without returning "None". 
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput 
        elif type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        simulatedResponses = np.atleast_2d(simulatedResponses)
        return simulatedResponses
    
    def getLogLikelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        #First do upper and lower bounds checks, if such bounds have been provided.
        boundsChecksPassed = self.doInputParameterBoundsChecks(discreteParameterVector)
        if boundsChecksPassed == False: #If false, return a 'zero probability' type result. Else, continue getting log likelihood.
            return float('-inf'), None #This approximates zero probability.

        #Check if user has provided a custom log likelihood function.
        if type(self.UserInput.model['custom_logLikelihood']) != type(None):
            logLikelihood, simulatedResponses = self.UserInput.model['custom_logLikelihood'](discreteParameterVector)
            simulatedResponses = np.array(simulatedResponses).flatten()
            return logLikelihood, simulatedResponses
        #else pass is implied.
        
        #Now get the simulated responses.
        simulatedResponses = self.getSimulatedResponses(discreteParameterVector)
        if type(simulatedResponses) == type(None):
            return float('-inf'), None #This is intended for the case that the simulation fails, indicated by receiving an 'nan' or None type from user's simulation function.
        
        #need to check if there are any 'responses_simulation_uncertainties'.
        if type(self.UserInput.model['responses_simulation_uncertainties']) == type(None): #if it's a None type, we keep it as a None type
            responses_simulation_uncertainties = None
        else:  #Else we get it based on the the discreteParameterVector
            responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(discreteParameterVector)

        #Now need to do transforms. Transforms are only for calculating log likelihood. If responses_simulation_uncertainties is "None", then we need to have one less argument passed in and a blank list is returned along with the transformed simulated responses.
        if type(responses_simulation_uncertainties) == type(None):
            simulatedResponses_transformed, blank_list = self.transform_responses(simulatedResponses) #This creates transforms for any data that we might need it. The same transforms were also applied to the observed responses.
            responses_simulation_uncertainties_transformed = None
            simulated_responses_covmat_transformed = None
        else:
            simulatedResponses_transformed, responses_simulation_uncertainties_transformed = self.transform_responses(simulatedResponses, responses_simulation_uncertainties) #This creates transforms for any data that we might need it. The same transforms were also applied to the observed responses.
            simulated_responses_covmat_transformed = returnShapedResponseCovMat(self.UserInput.num_response_dimensions, responses_simulation_uncertainties_transformed)  #assume we got standard deviations back.
        observedResponses_transformed = self.UserInput.responses_observed_transformed
        simulatedResponses_transformed_flattened = np.array(simulatedResponses_transformed).flatten()
        observedResponses_transformed_flattened = np.array(observedResponses_transformed).flatten()
        #If our likelihood is  “probability of Response given Theta”…  we have a continuous probability distribution for both the response and theta. That means the pdf  must use binning on both variables. Eric notes that the pdf returns a probability density, not a probability mass. So the pdf function here divides by the width of whatever small bin is being used and then returns the density accordingly. Because of this, our what we are calling likelihood is not actually probability (it’s not the actual likelihood) but is proportional to the likelihood.
        #Thus we call it a probability_metric and not a probability. #TODO: consider changing names of likelihood and get likelihood to "likelihoodMetric" and "getLikelihoodMetric"
        #Now we need to make the comprehensive_responses_covmat.
        #First we will check whether observed_responses_covmat_transformed is square or not. The multivariate_normal.pdf function requires a diagonal values vector to be 1D.
        observed_responses_covmat_transformed = self.observed_responses_covmat_transformed
        observed_responses_covmat_transformed_shape = np.shape(observed_responses_covmat_transformed) 
                
        #In general, the covmat could be a function of the responses magnitude and independent variables. So eventually, we will use non-linear regression or something to estimate it. However, for now we simply take the observed_responses_covmat_transformed which will work for most cases.
        #TODO: use Ashi's nonlinear regression code (which  he used in this paper https://www.sciencedirect.com/science/article/abs/pii/S0920586118310344).  Put in the response magnitudes and the independent variables.
        #in future it will be something like: if self.UserInput.covmat_regression== True: comprehensive_responses_covmat = nonLinearCovmatPrediction(self.UserInput['independent_variable_values'], observed_responses_covmat_transformed)
        #And that covmat_regression will be on by default.  We will need to have an additional argument for people to specify whether magnitude weighting and independent variable values should both be considered, or just one.
        if type(simulated_responses_covmat_transformed) == type(None):
            comprehensive_responses_covmat = observed_responses_covmat_transformed
        else: #Else we add the uncertainties, assuming they are orthogonal. Note that these are already covmats so are already variances that can be added directly.
            comprehensive_responses_covmat = observed_responses_covmat_transformed + simulated_responses_covmat_transformed
        comprehensive_responses_covmat_shape = copy.deepcopy(observed_responses_covmat_transformed_shape) #no need to take the shape of the actual comprehensive_responses_covmat since they must be same. This is probably slightly less computation.
        if len(comprehensive_responses_covmat_shape) == 1: #Matrix is square because has only one value.
            log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed_flattened,x=observedResponses_transformed_flattened,cov=comprehensive_responses_covmat)
        elif comprehensive_responses_covmat_shape[0] == comprehensive_responses_covmat_shape[1]:  #Else it is 2D, check if it's square.
            log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed_flattened,x=observedResponses_transformed_flattened,cov=comprehensive_responses_covmat)
            #TODO: Put in near-diagonal solution described in github: https://github.com/AdityaSavara/CheKiPEUQ/issues/3
        else:  #If it is not square, it's a list of variances so we need to take the 1D vector version.
            try:
                log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed_flattened,x=observedResponses_transformed_flattened,cov=comprehensive_responses_covmat[0])                
            except:
                log_probability_metric = 1
            if log_probability_metric == 1:
                log_probability_metric = -1E100 #Just initializing, then will add each probability separately.
                for responseValueIndex in range(len(simulatedResponses_transformed_flattened)):
                    try:
                        current_log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed_flattened[responseValueIndex],x=observedResponses_transformed_flattened[responseValueIndex],cov=comprehensive_responses_covmat[0][responseValueIndex])    
                    except: #The above is to catch cases when the multivariate_normal fails.
                        current_log_probability_metric = float('-inf')
                    log_probability_metric = current_log_probability_metric + log_probability_metric
                    if float(current_log_probability_metric) == float('-inf'):
                        print("Warning: There are posterior points that have zero probability. If there are too many points like this, the MAP and mu_AP returned will not be meaningful.")
                        current_log_probability_metric = -1E100 #Just choosing an arbitrarily very severe penalty. I know that I have seen 1E-48 to -303 from the multivariate pdf, and values inbetween like -171, -217, -272. I found that -1000 seems to be worse, but I don't have a systematic testing. I think -1000 was causing numerical errors.
                        log_probability_metric = current_log_probability_metric + log_probability_metric
        return log_probability_metric, simulatedResponses.flatten()

    def makeHistogramsForEachParameter(self):
        import CheKiPEUQ.plotting_functions as plotting_functions 
        parameterSamples = self.post_burn_in_samples
        parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        plotting_functions.makeHistogramsForEachParameter(parameterSamples,parameterNamesAndMathTypeExpressionsDict)

    def makeSamplingScatterMatrixPlot(self, parameterSamples = [], parameterNamesAndMathTypeExpressionsDict={}, parameterNamesList =[], plot_settings={}):
        if 'dpi' not in  plot_settings:  plot_settings['dpi'] = 220
        if 'figure_name' not in  plot_settings:  plot_settings['figure_name'] = 'scatter_matrix_posterior'
        if parameterSamples  ==[] : parameterSamples = self.post_burn_in_samples
        if parameterNamesAndMathTypeExpressionsDict == {}: parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        if parameterNamesList == []: parameterNamesList = self.UserInput.parameterNamesList #This is created when the parameter_estimation object is initialized.        

        posterior_df = pd.DataFrame(parameterSamples,columns=[parameterNamesAndMathTypeExpressionsDict[x] for x in parameterNamesList])
        pd.plotting.scatter_matrix(posterior_df)
        plt.savefig(plot_settings['figure_name'],dpi=plot_settings['dpi'])
        
    def createSimulatedResponsesPlots(self, allResponses_x_values=[], allResponsesListsOfYArrays =[], plot_settings={},allResponsesListsOfYUncertaintiesArrays=[] ): 
        #allResponsesListsOfYArrays  is to have 3 layers of lists: Response > Responses Observed, mu_guess Simulated Responses, map_Simulated Responses, (mu_AP_simulatedResponses) > Values
        if allResponses_x_values == []: 
            allResponses_x_values = np.atleast_2d(self.UserInput.responses_abscissa)
        if allResponsesListsOfYArrays  ==[]: #In this case, we assume allResponsesListsOfYUncertaintiesArrays == [] also.
            allResponsesListsOfYUncertaintiesArrays = [] #Set accompanying uncertainties list to a blank list in case it is not already one. Otherwise appending would mess up indexing.
            simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.
            simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
            
            #We already have self.UserInput.responses_observed, and will use that below. So now we get the simulated responses for the guess, MAP, mu_ap etc.
            
            #Get mu_guess simulated output and responses. 
            self.mu_guess_SimulatedOutput = simulationFunction( self.UserInput.InputParameterInitialGuess) #Do NOT use self.UserInput.model['InputParameterInitialGuess'] because that won't work with reduced parameter space requests.
            if type(simulationOutputProcessingFunction) == type(None):
                self.mu_guess_SimulatedResponses = np.atleast_2d(self.mu_guess_SimulatedOutput)
            if type(simulationOutputProcessingFunction) != type(None):
                self.mu_guess_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.mu_guess_SimulatedOutput)     )
            #Check if we have simulation uncertainties, and populate if so.
            if type(self.UserInput.model['responses_simulation_uncertainties']) != type(None):
                self.mu_guess_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.UserInput.InputParameterInitialGuess)
                
            #Get map simiulated output and simulated responses.
            self.map_SimulatedOutput = simulationFunction(self.map_parameter_set)           
            if type(simulationOutputProcessingFunction) == type(None):
                self.map_SimulatedResponses = np.atleast_2d(self.map_SimulatedOutput)
            if type(simulationOutputProcessingFunction) != type(None):
                self.map_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.map_SimulatedOutput)     )
            #Check if we have simulation uncertainties, and populate if so.
            if type(self.UserInput.model['responses_simulation_uncertainties']) != type(None):
                self.map_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.map_parameter_set)
            
            if hasattr(self, 'mu_AP_parameter_set'): #Check if a mu_AP has been assigned. It is normally only assigned if mcmc was used.           
                #Get mu_AP simiulated output and simulated responses.
                self.mu_AP_SimulatedOutput = simulationFunction(self.mu_AP_parameter_set)
                if type(simulationOutputProcessingFunction) == type(None):
                    self.mu_AP_SimulatedResponses = np.atleast_2d(self.mu_AP_SimulatedOutput)
                if type(simulationOutputProcessingFunction) != type(None):
                    self.mu_AP_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.mu_AP_SimulatedOutput)      )
                #Check if we have simulation uncertainties, and populate if so.
                if type(self.UserInput.model['responses_simulation_uncertainties']) != type(None):
                    self.mu_AP_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.mu_AP_parameter_set)
            
            #Now to populate the allResponsesListsOfYArrays and the allResponsesListsOfYUncertaintiesArrays
            for responseDimIndex in range(self.UserInput.num_response_dimensions):
                if not hasattr(self, 'mu_AP_parameter_set'): #Check if a mu_AP has been assigned. It is normally only assigned if mcmc was used.                
                    listOfYArrays = [self.UserInput.responses_observed[responseDimIndex], self.mu_guess_SimulatedResponses[responseDimIndex], self.map_SimulatedResponses[responseDimIndex]]        
                    allResponsesListsOfYArrays.append(listOfYArrays)
                    #Now to do uncertainties, there are two cases. First case is with only observed uncertainties and no simulation ones.
                    if type(self.UserInput.model['responses_simulation_uncertainties']) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                        allResponsesListsOfYUncertaintiesArrays.append( [self.UserInput.responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                    else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                        allResponsesListsOfYUncertaintiesArrays.append([self.UserInput.responses_observed_uncertainties[responseDimIndex],self.mu_guess_responses_simulation_uncertainties,self.map_responses_simulation_uncertainties]) #We need to give a list for each response dimension.                    
                if hasattr(self, 'mu_AP_parameter_set'):
                    listOfYArrays = [self.UserInput.responses_observed[responseDimIndex], self.mu_guess_SimulatedResponses[responseDimIndex], self.map_SimulatedResponses[responseDimIndex], self.mu_AP_SimulatedResponses[responseDimIndex]]        
                    allResponsesListsOfYArrays.append(listOfYArrays)
                    if type(self.UserInput.model['responses_simulation_uncertainties']) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                        allResponsesListsOfYUncertaintiesArrays.append( [self.UserInput.responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                    else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                        allResponsesListsOfYUncertaintiesArrays.append([self.UserInput.responses_observed_uncertainties[responseDimIndex],self.mu_guess_responses_simulation_uncertainties,self.map_responses_simulation_uncertainties,self.mu_AP_responses_simulation_uncertainties]) #We need to give a list for each response dimension.                    

        if plot_settings == {}: 
            plot_settings = self.UserInput.simulated_response_plot_settings
            if 'legendLabels' not in plot_settings: #The normal case:
                if hasattr(self, 'mu_AP_parameter_set'): 
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'MAP','mu_AP']
                else: #Else there is no mu_AP.
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'MAP']
                if hasattr(self, "opt_SSR"): #This means we are actually doing an optimization, and self.opt_SSR exists.
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'CPE']
                    print("line 1201 got here!!!")
            #Other allowed settings are like this, but will be fed in as simulated_response_plot_settings keys rather than plot_settings keys.
            #plot_settings['x_label'] = 'T (K)'
            #plot_settings['y_label'] = r'$rate (s^{-1})$'
            #plot_settings['y_range'] = [0.00, 0.025] #optional.
            #plot_settings['figure_name'] = 'tprposterior'
        if 'figure_name' not in plot_settings:
            plot_settings['figurename'] = 'Posterior'
        import CheKiPEUQ.plotting_functions as plotting_functions
        allResponsesFigureObjectsList = []
        for responseDimIndex in range(self.UserInput.num_response_dimensions): #TODO: Move the exporting out of the plot creation and/or rename the function and possibly have options about whether exporting graph, data, or both.
            #Some code for setting up individual plot settings in case there are multiple response dimensions.
            individual_plot_settings = copy.deepcopy(plot_settings) #we need to edit the plot settings slightly for each plot.
            if self.UserInput.num_response_dimensions == 1:
                responseSuffix = '' #If there is only 1 dimension, we don't need to add a suffix to the files created. That would only confuse people.
            if self.UserInput.num_response_dimensions > 1:
                responseSuffix = "_"+str(responseDimIndex)
            individual_plot_settings['figure_name'] = individual_plot_settings['figure_name']+responseSuffix
            if 'x_label' in plot_settings:
                if type(plot_settings['x_label']) == type(['list']) and len(plot_settings['x_label']) > 1: #the  label can be a single string, or a list of multiple response's labels. If it's a list of greater than 1 length, then we need to use the response index.
                    individual_plot_settings['x_label'] = plot_settings['x_label'][responseDimIndex]
            if 'y_label' in plot_settings:
                if type(plot_settings['y_label']) == type(['list']) and len(plot_settings['y_label']) > 1: #the  label can be a single string, or a list of multiple response's labels. If it's a list of greater than 1 length, then we need to use the response index.
                    individual_plot_settings['y_label'] = plot_settings['y_label'][responseDimIndex]                
            #TODO, low priority: we can check if x_range and y_range are nested, and thereby allow individual response dimension values for those.                               
            if np.shape(allResponses_x_values)[0] == 1: #This means a single abscissa for all responses.
                figureObject = plotting_functions.createSimulatedResponsesPlot(allResponses_x_values[0], allResponsesListsOfYArrays[responseDimIndex], individual_plot_settings, listOfYUncertaintiesArrays=allResponsesListsOfYUncertaintiesArrays[responseDimIndex])
                np.savetxt(individual_plot_settings['figure_name']+".csv", np.vstack((allResponses_x_values[0], allResponsesListsOfYArrays[responseDimIndex])).transpose(), delimiter=",", header='x_values, observed, sim_initial_guess, sim_MAP, sim_mu_AP', comments='')
            if np.shape(allResponses_x_values)[0] > 1: #This means a separate abscissa for each response.
                figureObject = plotting_functions.createSimulatedResponsesPlot(allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex], individual_plot_settings, listOfYUncertaintiesArrays=allResponsesListsOfYUncertaintiesArrays[responseDimIndex])
                np.savetxt(individual_plot_settings['figure_name']+".csv", np.vstack((allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex])).transpose(), delimiter=",", header='x_values, observed, sim_initial_guess, sim_MAP, sim_mu_AP', comments='')
            allResponsesFigureObjectsList.append(figureObject)
        return allResponsesFigureObjectsList  #This is a list of matplotlib.pyplot as plt objects.

    def createMumpcePlots(self):
        import CheKiPEUQ.plotting_functions as plotting_functions
        from CheKiPEUQ.plotting_functions import plotting_functions_class
        figureObject_beta = plotting_functions_class(self.UserInput) # The "beta" is only to prevent namespace conflicts with 'figureObject'.
        parameterSamples = self.post_burn_in_samples
        
        #TODO: the posterior mu_vector and cov_matrix should be calculated elsewhere.
        posterior_mu_vector = np.mean(parameterSamples, axis=0)
        posterior_cov_matrix = np.cov(self.post_burn_in_samples.T)
        self.posterior_cov_matrix = posterior_cov_matrix
        #TODO: In future, worry about whether there are constants or not, since then we will have to trim down the prior.
        #Make the model_parameter_info object that mumpce Project class needs.
        self.UserInput.model_parameter_info = []#This variable name is for mumpce definition of variable names. Not what we would choose otherwise.
        for parameterIndex, parameterName in enumerate(self.UserInput.parameterNamesAndMathTypeExpressionsDict):
            individual_model_parameter_dictionary = {'parameter_number': parameterIndex, 'parameter_name': self.UserInput.parameterNamesAndMathTypeExpressionsDict[parameterName]} #we are actually putting the MathTypeExpression as the parameter name when feeding to mum_pce.
            self.UserInput.model_parameter_info.append(individual_model_parameter_dictionary)
        self.UserInput.model_parameter_info = np.array(self.UserInput.model_parameter_info)
        if len(self.UserInput.active_parameters) == 0:
            numParams = len(self.UserInput.model_parameter_info)
            active_parameters = np.linspace(0, numParams-1, numParams) #just a list of whole numbers.
            active_parameters = np.array(active_parameters, dtype='int')
        else:
            active_parameters = self.UserInput.active_parameters
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
                firstParameter = int(self.UserInput.parameterNamesAndMathTypeExpressionsDict[pairIndex[0]])
                secondParameter = int(self.UserInput.parameterNamesAndMathTypeExpressionsDict[pairIndex[0]])
                pairs_of_parameter_indices[pairIndex] = [firstParameter, secondParameter]        
        figureObject_beta.mumpce_plots(model_parameter_info = self.UserInput.model_parameter_info, active_parameters = active_parameters, pairs_of_parameter_indices = pairs_of_parameter_indices, posterior_mu_vector = posterior_mu_vector, posterior_cov_matrix = posterior_cov_matrix, prior_mu_vector = np.array(self.UserInput.mu_prior), prior_cov_matrix = self.UserInput.covmat_prior, contour_settings_custom = self.UserInput.contour_settings_custom)
        return figureObject_beta

    @CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def createAllPlots(self):
        try:
            self.makeHistogramsForEachParameter()    
            self.makeSamplingScatterMatrixPlot()
        except:
            pass

        try:        
            self.createMumpcePlots()
        except:
            pass

        try:
            self.createSimulatedResponsesPlots()
        except:
            pass

class verbose_optimization_wrapper: #Learned how to use callback from Henri's post https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
    def __init__(self, simulationFunction):
        self.simulationFunction = simulationFunction
        self.FirstCall = True # Just intializing.
        self.iterationNumber = 0 # Just intializing.
    
    def simulateAndStoreObjectiveFunction(self, discreteParameterVector):
        #This class function is what we feed to the optimizer. It mainly keeps track of what has been tried so far.
        simulationOutput = self.simulationFunction(discreteParameterVector) # the actual evaluation of the function
        self.lastTrialDiscreteParameterVector = discreteParameterVector
        self.lastTrialObjectiveFunction = simulationOutput
        return simulationOutput
    
    def callback(self, discreteParameterVector, *extraArgs):
        #This class function has to be passed in as the callback function argument to the optimizer.
        #basically, it gets 'called' between iterations of the optimizer.
        #Some optimizers give back extra args, so there is a *extraArgs argument above.
        if self.FirstCall == True:
            parameterNamesString = ""
            for parameterIndex in range(len(discreteParameterVector)):
                parameterName = f"Par-{parameterIndex+1}"
                parameterNamesString += f"{parameterName:10s}\t"
            headerString = "Iter  " + parameterNamesString + "ObjectiveF"
            print(headerString)
            self.FirstCall = False
        
        iterationNumberString = "{0:4d}  ".format(self.iterationNumber)
        discreteParameterVector = self.lastTrialDiscreteParameterVector #We take the stored one rather than the one provided to make sure that we're getting the same one as the stored objective function.
        parameterValuesString = ""
        for parameterValue in discreteParameterVector:
            parameterValuesString += f"{parameterValue:10.5e}\t"
        currentObjectiveFunctionValue = f"{self.lastTrialObjectiveFunction:10.5e}"
        iterationOutputString = iterationNumberString + parameterValuesString + currentObjectiveFunctionValue
        print(iterationOutputString)
        self.iterationNumber += 1 #In principle, could be done inside the simulateAndStoreObjectiveFunction, but this way it is after the itration number has been printed.

'''Below are a bunch of functions for Euler's Method.'''
#This takes an array of dydt values. #Note this is a local dydtArray, it is NOT a local deltaYArray.
software_name = "Integrated Production (Objective Function)"
software_version = "1.0.0"
software_unique_id = "https://doi.org/10.1016/j.susc.2016.07.001"
software_kwargs = {"version": software_version, "author": ["Aditya Savara"], "doi": "https://doi.org/10.1016/j.susc.2016.07.001", "cite": "Savara, Aditya. 'Simulation and fitting of complex reaction network TPR: The key is the objective function.' Surface Science 653 (2016): 169-180."} 
@CiteSoft.module_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
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
def littleEulerUncertaintyPropagation(dydt_uncertainties, t_values, initial_y_uncertainty=0, forceNonzeroInitialUncertainty=True):
    y_uncertainties = dydt_uncertainties*0.0
    y_uncertainties[0] = initial_y_uncertainty #We have no way to make an uncertainty for point 0.
    for index in range(len(dydt_uncertainties)-1): #The uncertainty for each next point is propagated through the uncertainty of the current value and the delta_t*(dy/dt uncertainty), since we are adding two values.
        deltat_resolution = t_values[index+1]-t_values[index]
        y_uncertainties[index+1] = ((y_uncertainties[index])**2+(dydt_uncertainties[index]*deltat_resolution)**2)**0.5
    if forceNonzeroInitialUncertainty==True:
        if initial_y_uncertainty == 0: #Errors are caused if initial_y_uncertainty is left as zero, so we take the next uncertainty as an assumption for a reasonable base estimate of the initial point uncertainty.
            y_uncertainties[0] = y_uncertainties[1]   
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


#TODO: move this into some kind of support module for parsing. Like XYYYDataFunctions or something like that.
def returnReducedIterable(iterableObjectToReduce, reducedIndices):
    #If a numpy array or list is provided, the same will be returned. Else, a list will be returned.
    #For arrays, only 1D and square 2D are supported. Anything else will only do the first axis.
    reducedIterable = copy.deepcopy(iterableObjectToReduce) #Doing this initially so that unsupported cases will still return something.
    
    #In most cases, we use a little function that makes a list to do the reduction.
    def returnReducedList(iterableObjectToReduce, reducedIndices):
        reducedList = [] #just initializing.
        for elementIndex,element in enumerate(iterableObjectToReduce):
            if elementIndex in reducedIndices:
                reducedList.append(element)
        return reducedList

    #Now to do the actual reduction.
    if type(iterableObjectToReduce)== type(np.array([0])):
        if len(np.shape(iterableObjectToReduce)) == 1: #If it's 1D, we can just use a list and convert back to numpy array.
            reducedIterableAsList = returnReducedList(iterableObjectToReduce, reducedIndices)
            reducedIterable = np.array(reducedIterableAsList)
        if len(np.shape(iterableObjectToReduce)) == 2: #If it's a 2D square matrix, then we will still support it.
            if np.shape(iterableObjectToReduce)[0] == np.shape(iterableObjectToReduce)[1]: #Make sure it is square before trying to do more:
                #FIRST GO ACROSS THE ROWS.
                reducedIterableAsList = returnReducedList(iterableObjectToReduce, reducedIndices)
                partiallyReducedIterable = np.array(reducedIterableAsList)
                #NOW TRANSPOSE, DO IT AGAIN, AND THEN TRANSPOSE BACK.
                partiallyReducedIterable = partiallyReducedIterable.transpose()
                reducedIterableAsList = returnReducedList(partiallyReducedIterable, reducedIndices)
                reducedIterable = np.array(reducedIterableAsList).transpose() #convert to array and transpose
            else: #If it's 2D but not square, we just reduce along the row axis (main axis)
                reducedIterableAsList = returnReducedList(iterableObjectToReduce, reducedIndices)
                reducedIterable = np.array(reducedIterableAsList)
    else: # the following is included in the else, type(iterableObjectToReduce)== type(['list']):
        reducedIterable = returnReducedList(iterableObjectToReduce, reducedIndices)
    if np.shape(reducedIterable) == np.shape(iterableObjectToReduce):
        print("returnReducedIterable received an object type or size that is not supported.")
    return reducedIterable

def returnShapedResponseCovMat(numResponseDimensions, uncertainties):
    #The uncertainties, whether transformed or not, must be one of the folllowing: a) for a single dimension response can be a 1D array of standard deviations, b) for as ingle dimension response can be a covmat already (so already variances), c) for a multidimensional response we *only* support standard deviations at this time.
    if numResponseDimensions == 1:
        shapedUncertainties = np.array(uncertainties, dtype="float") #Initializing variable. 
        if np.shape(shapedUncertainties)[0] == (1): #This means it's just a list of standard deviations and needs to be squared to become variances.
            shapedUncertainties = np.square(shapedUncertainties) # Need to square standard deviations to make them into variances.
        else:
            shapedUncertainties = shapedUncertainties
    elif numResponseDimensions > 1:  #if the dimensionality of responses is greater than 1, we only support providing standard deviations. Will flatten and square.
        shapedUncertainties = np.array(uncertainties, dtype="float") #Filling variable.  
        shapedUncertainties = shapedUncertainties.flatten() 
        shapedUncertainties = np.square(shapedUncertainties) #Need to square standard deviations to make them into variances.
    return shapedUncertainties

def boundsCheck(parameters, parametersBounds, boundsType):
    #Expects three arguments.
    #the first two are 1D array like arguments (parameters and a set of *either* upper bounds or lower bounds)
    #The third argumment is the type of bounds, either 'upper' or 'lower'
    #In practice, this means the function usually needs to be called twice.
    #A "None" type is expected for something that is not bounded in that direction. 
    
    #We first need to make arrays and remove anything that is None in the bounds.
    parameters = np.array(parameters).flatten()
    parametersBounds = np.array(parametersBounds).flatten()
    #to remove, we use brackets that pull out the indices where the comparison is not None. This is special numpy array syntax.
    parametersTruncated = parameters[parametersBounds != None]
    parametersBoundsTruncated = parametersBounds[parametersBounds != None]    
    if boundsType.lower() == 'upper': #we make the input into lower case before proceeding.
        upperCheck = parametersTruncated < parametersBoundsTruncated #Check if all are smaller.
        if False in upperCheck: #If any of them failed, we return False.
            return False
        else:
            pass #else we do the lower bounds check next.
    if boundsType.lower() == 'lower':
        lowerCheck = parametersTruncated > parametersBoundsTruncated #Check if all are smaller.
        if False in lowerCheck: #If any of them failed, we return False.
            return False
        else:
            pass
    return True #If we have gotten down to here without returning False, both checks have passed and we return true.
        
if __name__ == "__main__":
    pass

