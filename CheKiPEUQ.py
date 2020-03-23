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
    #Inside this class, a UserInput namespace is provided. This has dictionaries of UserInput choices.
    #However, the code initally parses those choices and then puts processed versions in the SAME name space, but no longer in the dictionaries.
    #So functions in this class should (when possible) call the namespace variables that are not in dictionaries, unless the original userinput is desired.
    #'inverse problem'. Initialize chain with initial guess (prior if not provided) as starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  
    def __init__(self, UserInput = None):
        self.UserInput = UserInput #Note that this is a pointer, so the later lines are within this object.
        #Now will automatically populate some variables from UserInput
        UserInput.parameterNamesList = list(UserInput.model['parameterNamesAndMathTypeExpressionsDict'].keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        UserInput.parameterNamesAndMathTypeExpressionsDict = UserInput.model['parameterNamesAndMathTypeExpressionsDict']
        if self.UserInput.parameter_estimation_settings['verbose']: 
            print("Bayes Model Initialized")
        #Leaving the original dictionary object intact, but making a new object to make covmat_prior.
        UserInput.InputParametersPriorValuesUncertainties = UserInput.model['InputParametersPriorValuesUncertainties']     
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

        self.UserInput.mu_prior = np.array(UserInput.model['InputParameterPriorValues']) 
        #Making things at least 2d.
        UserInput.responses['responses_abscissa'] = np.atleast_2d(UserInput.responses['responses_abscissa'])
        UserInput.responses['responses_observed'] = np.atleast_2d(UserInput.responses['responses_observed'])
        UserInput.responses['responses_observed_uncertainties'] = np.atleast_2d(UserInput.responses['responses_observed_uncertainties'])


        UserInput.responses['responses_observed_transformed'], UserInput.responses['responses_observed_transformed_uncertainties']  = self.transform_responses(UserInput.responses['responses_observed'], UserInput.responses['responses_observed_uncertainties']) #This creates transforms for any data that we might need it. The same transforms will also be applied during parameter estimation.
             
        self.UserInput.num_data_points = len(UserInput.responses['responses_observed'].flatten()) #FIXME: This is only true for transient data.
        #Now scale things as needed:
        if UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "std":
            self.UserInput.scaling_uncertainties = UserInput.std_prior #Could also be by mu_prior.  The reason a separate variable is made is because this will be used in the getPrior function as well, and having a separate variable makes it easier to trace. This scaling helps prevent numerical errors in returning the pdf.
        elif UserInput.parameter_estimation_settings['scaling_uncertainties_type'] == "mu":
            self.UserInput.scaling_uncertainties = UserInput.mu_prior
        #TODO: consider a separate scaling for each variable, taking the greater of either mu_prior or std_prior.
        self.UserInput.mu_prior_scaled = np.array(UserInput.mu_prior/UserInput.scaling_uncertainties)
        self.UserInput.var_prior_scaled = np.array(UserInput.var_prior/(UserInput.scaling_uncertainties*UserInput.scaling_uncertainties))
        self.UserInput.covmat_prior_scaled = self.UserInput.covmat_prior*1.0 #First initialize, then fill.
        for parameterIndex, parameterValue in enumerate(UserInput.scaling_uncertainties):
            UserInput.covmat_prior_scaled[parameterIndex,:] = UserInput.covmat_prior[parameterIndex,:]/parameterValue
            UserInput.covmat_prior_scaled[:,parameterIndex] = UserInput.covmat_prior[:,parameterIndex]/parameterValue           

        #To find the responses covariance matrix, we take the errors from the points. This is needed for the likelihood.
        self.UserInput.num_response_dimensions = np.shape(UserInput.responses['responses_abscissa'])[0] #The first index of shape is the num of responses, but has to be after at_least2d is performed.
        if self.UserInput.num_response_dimensions == 1:
            responses_covmat = np.array(self.UserInput.responses['responses_observed_transformed_uncertainties']) #Filling variable.
            if np.shape(responses_covmat)[0] == (1): #This means it's just a list of standard deviations and needs to be squared to become variances.
                responses_covmat = np.square(responses_covmat)
            else:
                responses_covmat = responses_covmat
        elif self.UserInput.num_response_dimensions > 1:  #if the dimensionality of responses is greater than 1, we only support providing standard deviations. Will flatten and square.
            responses_covmat = np.array(self.UserInput.responses['responses_observed_transformed_uncertainties']) #Filling variable.
            responses_covmat = responses_covmat.flatten() 
            responses_covmat = np.square(responses_covmat)
        self.responses_covmat = responses_covmat            
                       
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
        
        #Now reduce the parameter space if requested by the user. #TODO: consider having this as a function outside of init
        if len(self.UserInput.model['reducedParameterSpace']) > 0:
            print("Important: the UserInput.model['reducedParameterSpace'] is not blank. That means the only parameters allowed to change will be the ones in the indices inside 'reducedParameterSpace'.   All others will be held constant.  The values inside  'InputParameterInitialGuess will be used', and 'InputParameterPriorValues' if an initial guess was not provided.")
            self.reduceParameterSpace()
    
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

    def simulateWithSubsetOfParameters(self,reducedParametersVector): #This is a wrapper.
        #This function has implied arguments of ...
        #self.UserInput.model['InputParameterInitialGuess'] for the parameters to start with
        #self.UserInput.model['reducedParameterSpace'] a list of indices for which parameters are the only ones to change.
        #simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        #simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']
        #When this wrapper is used, EVERYWHERE ELSE will call it to do the simulation, by calling self.UserInput.simulationFunction and self.UserInput.simulationOutputProcessingFunction
        simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']
        
        #now populate the discreteParameterVector first with the initial guess, then with the new reducedParameters vector.
        discreteParameterVector = copy.deepcopy(self.UserInput.model['InputParameterInitialGuess']) #This is the original one from the user, before any reduction.
        for reducedParameterIndex, parameterValue in enumerate(reducedParametersVector):
            #we find which index to put things into from #self.UserInput.model['reducedParameterSpace'], which is a list of indices.
            regularParameterIndex = self.UserInput.model['reducedParameterSpace'][reducedParameterIndex]
            discreteParameterVector[regularParameterIndex] = parameterValue
        try:
            simulationOutput = simulationFunction(discreteParameterVector) 
        except:
            return 0, None #This is for the case that the simulation fails. Should be made better in future.
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        
        simulatedResponses = np.atleast_2d(simulatedResponses)
        #This is not needed:
        #observedResponses = np.atleast_2d(self.UserInput.responses['responses_observed'])
        return simulatedResponses
    
    def transform_responses(self, nestedAllResponsesArray, nestedAllResponsesUncertainties = []):
        nestedAllResponsesArray_transformed = copy.deepcopy(nestedAllResponsesArray) #First make a copy to populate with transformed values.
        nestedAllResponsesUncertainties_transformed = copy.deepcopy(nestedAllResponsesUncertainties) #First make a copy to populate with transformed values. If blank, we won't populate it.        
        UserInput = self.UserInput
        
        #TODO: Make little function for interpolation in case it's necessary (see below).
#        def littleInterpolator():
#            abscissaRange = UserInput.responses['responses_abscissa'][responseIndex][-1] - UserInput.responses['responses_abscissa'][responseIndex][0] #Last value minus first value.
#            UserInput.responses['responses_observed'] = np.atleast_2d(UserInput.responses['responses_observed'])
#            UserInput.responses['responses_observed_uncertainties'] = np.atleast_2d(UserInput.responses['responses_observed_uncertainties'])
        if 'kinetics_type' not in UserInput.model:  #To make backwards compatibility.
            UserInput.model['kinetics_type'] = ''
        if UserInput.model['kinetics_type'] == 'transient': #This assumes that the abscissa is always time.
            for responseIndex, response in enumerate(UserInput.responses['responses_observed']):
                #We will need the abscissa also, so need to check if there are independent abscissa or not:
                if len(UserInput.responses['responses_abscissa']) == 1: #This means there is only one abscissa.
                    abscissaIndex = 0
                else:
                    abscissaIndex = responseIndex
                #Now to do the transforms.
                if UserInput.responses['response_types'][responseIndex] == 'I':	 #For intermediate
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses['responses_abscissa'][abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses['responses_abscissa'][abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        #Perform the littleEuler twice.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses['responses_abscissa'][abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses['responses_abscissa'][abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses['responses_abscissa'][abscissaIndex], nestedAllResponsesArray_transformed[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties_transformed[responseIndex], UserInput.responses['responses_abscissa'][abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
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
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses['responses_abscissa'][abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses['responses_abscissa'][abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'o':
                    if UserInput.responses['response_data_type'][responseIndex] == 'o':
                        pass
                    if UserInput.responses['response_data_type'][responseIndex] == 'c':
                        LittleEuler
                    if UserInput.responses['response_data_type'][responseIndex] == 'r':
                        LittleEulerTwice
        if UserInput.model['kinetics_type'] == 'steady_state': #TODO: so far, this does not do anything. It assumes that the abscissa is never time.
            for responseIndex, response in enumerate(UserInput.responses['responses_observed']):
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
  
    def doGridSearch(self, searchType='doMetropolisHastings', export = True, verbose = False, gridSamplingIntervalSize = [], gridSamplingRadii = [], passThroughArgs = {}):
        # gridSamplingRadii is the number of variations to check in units of variance for each parameter. Can be 0 if you don't want to vary a particular parameter in the grid search.
        #TODO: the upper part of the gridsearch may not be compatibile with reduced parameter space. Needs to be checked.
        import CombinationGeneratorModule
        numParameters = len(self.UserInput.parameterNamesList)
        if len(gridSamplingRadii) == 0:
            gridSamplingRadii = np.ones(numParameters, dtype='int') #By default, will make ones.
            numGridPoints = 3**numParameters
        else: 
            gridSamplingRadii = np.array(gridSamplingRadii, dtype='int')
            numGridPoints = 1 #just initializing.
            for radius in gridSamplingRadii:
                numGridPoints=numGridPoints*(2*radius+1)
        if len(gridSamplingIntervalSize) == 0:
            gridSamplingIntervalSize = self.UserInput.std_prior #By default, we use the standard deviations associated with the priors.
        else: gridSamplingIntervalSize = np.array(gridSamplingRadii, dtype='float')
        gridCenter = self.UserInput.InputParameterInitialGuess #We take what is in the variable self.UserInput.InputParameterInitialGuess for the center of the grid.
        gridCombinations = CombinationGeneratorModule.combinationGenerator(gridCenter, gridSamplingIntervalSize, gridSamplingRadii, SpreadType="Addition",toFile=False)
        allGridResults = []
        #Initialize some things before loop.
        if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
                import timeit
                timeAtGridStart = timeit.time.clock()
                timeAtLastGridPoint = timeAtGridStart #just initializing
        highest_logP = float('-inf') #Just initializing.
        #Start grid search loop.
        for combinationIndex,combination in enumerate(gridCombinations):
            self.UserInput.InputParameterInitialGuess = combination #We need to fill the variable InputParameterInitialGuess with the combination being checked.
            if searchType == 'getLogP':
                thisResult = self.getLogP(combination)
                self.map_logP = thisResult #The getLogP function does not fill map_logP by itself.
                self.map_parameter_set = combination
            if searchType == 'doMetropolisHastings':
                thisResult = self.doMetropolisHastings()
            if searchType == 'doOptimizeNegLogP':
                thisResult = self.doOptimizeNegLogP(**passThroughArgs)
            if type(self.UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None):
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
                print("GridPoint", combination, "number", combinationIndex, "out of", numGridPoints, "timeOfThisGridPoint", timeOfThisGridPoint)
                print("GridPoint", combinationIndex, "averageTimePerGridPoint", "%.2f" % round(averageTimePerGridPoint,2), "estimated time remaining", "%.2f" % round( numRemainingGridPoints*averageTimePerGridPoint,2), "s" )
                print("GridPoint", combinationIndex, "current logP", self.map_logP, "highest logP", highest_logP)
        #TODO: export the allGridResults to file at end of search in a nicer format.        
        #First set the initial guess back to the center of the grid.
        self.UserInput.InputParameterInitialGuess = gridCenter
        #Now populate the map etc. with those of the best result.
        self.map_logP = highest_logP 
        self.map_parameter_set = highest_logP_parameter_set 
        with open("gridsearch_log_file.txt", 'w') as out_file:
            out_file.write("result: " + "self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec" + "\n")
            for resultIndex, result in enumerate(allGridResults):
                out_file.write("result:" + str(resultIndex) +  str(result) + "\n")
            print("Final map results from gridsearch:", self.map_parameter_set, "final logP:", self.map_logP)
        if searchType == 'doMetropolisHastings':
            #Metropolis hastings has other variables to populate.
            #[self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] =
            return bestResultSoFar # [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] 
        if searchType == 'doOptimizeNegLogP':            
            return bestResultSoFar# [self.map_parameter_set, self.map_logP]
        if searchType == 'simplegrid':          
            return bestResultSoFar# [self.map_parameter_set, self.map_logP]
            

    def getLogP(self, proposal_sample): #The proposal sample is specific parameter vector.
        [log_likelihood_proposal, simulationOutput_proposal] = self.getLogLikelihood(proposal_sample)
        log_prior_proposal = self.getLogPrior(proposal_sample)
        log_numerator_or_denominator = log_likelihood_proposal+log_prior_proposal #Of the Metropolis-Hastings accept/reject ratio
        return log_numerator_or_denominator
        
    def getNegLogP(self, proposal_sample): #The proposal sample is specific parameter vector. We are using negative of log P because scipy optimize doesn't do maximizing. It's recommended minimize the negative in this situation.
        neg_log_postererior = -1*self.getLogP(proposal_sample)
        return neg_log_postererior

    def doOptimizeNegLogP(self, simulationFunctionAdditionalArgs = (), method = None, optimizationAdditionalArgs = {}, printOptimum = True, verbose=True):
        #THe intention of the optional arguments is to pass them into the scipy.optimize.minimize function.
        # the 'method' argument is for Nelder-Mead, BFGS, SLSQP etc. https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize
        initialGuess = self.UserInput.InputParameterInitialGuess
        import scipy.optimize
        if verbose == False:
            optimizeResult = scipy.optimize.minimize(self.getNegLogP, initialGuess, method = method)
        if verbose == True:
            verbose_simulator = verbose_optimization_wrapper(self.getNegLogP)
            optimizeResult = scipy.optimize.minimize(verbose_simulator.simulate, initialGuess, method=method, callback=verbose_simulator.callback, options={"disp": True})
            #print(f"Number of calls to Simulator instance {verbose_simulator.num_calls}") <-- this is the same as the "Function evaluations" field that gets printed.
            
        self.map_parameter_set = optimizeResult.x #This is the map location.
        self.map_logP = -1.0*optimizeResult.fun #This is the map logP
        if printOptimum == True:
            print("Final results from doOptimizeNegLogP:", self.map_parameter_set, "final logP:", self.map_logP)
        return [self.map_parameter_set, self.map_logP]
    
    #main function to get samples #TODO: Maybe Should return map_log_P and mu_AP_log_P?
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
        for i in range(1,self.UserInput.parameter_estimation_settings['mcmc_length']): #FIXME: Don't we need to start with i of 0?
            if self.UserInput.parameter_estimation_settings['verbose']: print("MCMC sample number", i)                  
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'unbiased':
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_covmat*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
            if self.UserInput.parameter_estimation_settings['mcmc_mode'] == 'MAP_finding':
                if i == 1: mcmc_step_dynamic_coefficient = 1
                mcmc_step_modulation_coefficient = np.random.uniform() + 0.5 #TODO: make this a 2D array. One for each parameter.
                mcmc_step_modulation_history[i] = mcmc_step_modulation_coefficient
                proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_covmat*mcmc_step_dynamic_coefficient*mcmc_step_modulation_coefficient*self.UserInput.parameter_estimation_settings['mcmc_relative_step_length'])
            log_prior_proposal = self.getLogPrior(proposal_sample)
            [log_likelihood_proposal, simulationOutput_proposal] = self.getLogLikelihood(proposal_sample)
            log_prior_current_location = self.getLogPrior(samples[i-1,:]) 
            [log_likelihood_current_location, simulationOutput_current_location] = self.getLogLikelihood(samples[i-1,:]) #FIXME: the previous likelihood should be stored so that it doesn't need to be calculated again.
            log_accept_probability = (log_likelihood_proposal + log_prior_proposal) - (log_likelihood_current_location + log_prior_current_location) 
            #accept_probability = np.power(10, log_accept_probability) Don't use this!
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current log_likelihood',log_likelihood_current_location, 'Proposed log_likelihood', log_likelihood_proposal, '\nLog of Accept_probability (gauranteed if above 0)', log_accept_probability)
            if self.UserInput.parameter_estimation_settings['verbose']: print('Current posterior',log_likelihood_current_location+log_prior_current_location, 'Proposed Posterior', log_likelihood_proposal+log_prior_proposal)
            if self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'] != 0: #This flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy.
                N_flatten = float(self.UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability'])
                #This is e^logP = P. #This is base 'e' because the logpdf functions are base e. Ashi checked the sourcecode.
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
#                print("line 121", simulationOutput_current_location)
                log_posteriors_un_normed_vec[i] = log_likelihood_current_location+log_prior_current_location
                log_likelihoods_vec[i] = log_likelihood_current_location
                log_priors_vec[i] = log_prior_current_location
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
        self.post_burn_in_log_posteriors_un_normed_vec = log_posteriors_un_normed_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_logP_un_normed_vec = (self.post_burn_in_log_posteriors_un_normed_vec)
        self.post_burn_in_log_likelihoods_vec = log_likelihoods_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        self.post_burn_in_log_priors_vec = log_priors_vec[self.UserInput.parameter_estimation_settings['mcmc_burn_in']:]
        # posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        #FIXME: Log was not propagated correctly here. Below line used to be self.evidence = np.mean(self.post_burn_in_posteriors_un_normed_vec)*np.sqrt(2*np.pi*np.std(self.post_burn_in_samples)**2)
        #So either need to make post_burn_in_posteriors_un_normed_vec again before this step, or need to change below line.
        self.evidence = np.mean(self.post_burn_in_log_posteriors_un_normed_vec)*np.sqrt(2*np.pi*np.std(self.post_burn_in_samples)**2)
        post_burn_in_log_posteriors_vec = self.post_burn_in_log_posteriors_un_normed_vec/self.evidence
        log_ratios = (post_burn_in_log_posteriors_vec-self.post_burn_in_log_priors_vec) #log10(a/b) = log10(a)-log10(b)
        log_ratios[np.isinf(log_ratios)] = 0
        log_ratios = np.nan_to_num(log_ratios)
        self.info_gain = np.mean(log_ratios)
        map_logP = max(self.post_burn_in_logP_un_normed_vec)
        self.map_logP = map_logP
        self.map_index = list(self.post_burn_in_logP_un_normed_vec).index(map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
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
                out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
                out_file.write("evidence:" + str(self.evidence) + "\n")
                if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.var_prior).any() > 0.10:
                    out_file.write("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% of prior variance in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632")
        if abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.var_prior).any() > 0.10:  
            print("Warning: The MAP parameter set and mu_AP parameter set differ by more than 10% of prior variance in at least one parameter. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632  ")
        return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_logP_un_normed_vec] # EAW 2020/01/08
    def getLogPrior(self,discreteParameterVector):
        discreteParameterVector_scaled = np.array(discreteParameterVector/self.UserInput.scaling_uncertainties)
        log_probabilityPrior = multivariate_normal.logpdf(x=discreteParameterVector_scaled,mean=self.UserInput.mu_prior_scaled,cov=self.UserInput.covmat_prior_scaled)
        return log_probabilityPrior
    def getLogLikelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.  
        simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
        try:
            simulationOutput =simulationFunction(discreteParameterVector) 
        except:
            return 0, None #This is for the case that the simulation fails. Should be made better in future.
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        simulatedResponses = np.atleast_2d(simulatedResponses)
        simulatedResponses_transformed, blank_list = self.transform_responses(simulatedResponses) #This creates transforms for any data that we might need it. The same transforms were also applied to the observed responses.
        observedResponses_transformed = self.UserInput.responses['responses_observed_transformed']
        simulatedResponses_transformed_flattened = np.array(simulatedResponses_transformed).flatten()
        observedResponses_transformed_flattened = np.array(observedResponses_transformed).flatten()
 
        #If our likelihood is  “probability of Response given Theta”…  we have a continuous probability distribution for both the response and theta. That means the pdf  must use binning on both variables. Eric notes that the pdf returns a probability density, not a probability mass. So the pdf function here divides by the width of whatever small bin is being used and then returns the density accordingly. Because of this, our what we are calling likelihood is not actually probability (it’s not the actual likelihood) but is proportional to the likelihood.
        #This we call it a probability_metric and not a probability. #TODO: consider changing likelihood and get likelihood to "likelihoodMetric" and "getLikelihoodMetric"
        
        #Now we will check whether responses_covmat is square or not. If it's square, we take it as is. If it's not square, we take the nested object inside since the multivariate_normal.pdf function requires a diagonal values vector to be 1D.
        responses_covmat = self.responses_covmat
        responses_covmat_shape = np.shape(responses_covmat)
        if len(responses_covmat_shape) == 1: #Matrix is square because has only one value.
            log_probability_metric = multivariate_normal.logpdf(x=simulatedResponses_transformed_flattened,mean=observedResponses_transformed_flattened,cov=responses_covmat)
        elif responses_covmat_shape[0] == responses_covmat_shape[1]:  #Else it is 2D, check if it's square.
            probability_metric = multivariate_normal.logpdf(x=simulatedResponses_transformed_flattened,mean=observedResponses_transformed_flattened,cov=responses_covmat)
            #TODO: Put in near-diagonal solution described in github.
        else:  #If it is not square, it's a list of variances so we need to take the 1D vector version.
            try:
                log_probability_metric = multivariate_normal.logpdf(x=simulatedResponses_transformed_flattened,mean=observedResponses_transformed_flattened,cov=responses_covmat[0])                
            except:
                log_probability_metric = 1
            if log_probability_metric == 1:
                log_probability_metric = -1E100 #Just initializing, then will add each probability separately.
                for responseValueIndex in range(len(simulatedResponses_transformed_flattened)):
                    current_log_probability_metric = multivariate_normal.logpdf(x=simulatedResponses_transformed_flattened[responseValueIndex],mean=observedResponses_transformed_flattened[responseValueIndex],cov=responses_covmat[0][responseValueIndex])    
                    log_probability_metric = current_log_probability_metric + log_probability_metric
                    if float(current_log_probability_metric) == float('-inf'):
                        print("Warning: There are posterior points that have zero probability. If there are too many points like this, the MAP and mu_AP returned will not be meaningful.")
                        current_log_probability_metric = -1E100 #Just choosing an arbitrarily very severe penalty. I know that I have seen 1E-48 to -303 from the multivariate pdf, and values inbetween like -171, -217, -272. I found that -1000 seems to be worse, but I don't have a systematic testing. I think -1000 was causing numerical errors.
                        log_probability_metric = current_log_probability_metric + log_probability_metric
            #print(log_probability_metric)
        return log_probability_metric, simulatedResponses.flatten()

    def makeHistogramsForEachParameter(self):
        import plotting_functions #This is going to become import CheKIPEUQ.plotting_functions as plotting_functions
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
        
    def createSimulatedResponsesPlots(self, allResponses_x_values=[], allResponsesListsOfYArrays =[], plot_settings={},allResponsesListsOfYUncertainties=[] ): 
        #allResponsesListsOfYArrays  is to have 3 layers of lists: Response > Responses Observed, mu_guess Simulated Responses, map_Simulated Responses, (mu_AP_simulatedResponses) > Values
        if allResponses_x_values == []: allResponses_x_values = self.UserInput.responses['responses_abscissa']       
        if allResponsesListsOfYUncertainties == []: allResponsesListsOfYUncertainties = self.UserInput.responses['responses_observed_uncertainties']
        if allResponsesListsOfYArrays  ==[]:
            
            simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.
            simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
            
            #Get mu_guess simulated output and responses. 
            self.mu_guess_SimulatedOutput = simulationFunction( self.UserInput.InputParameterInitialGuess) #Do NOT use self.UserInput.model['InputParameterInitialGuess'] because that won't work with reduced parameter space requests.
            if type(simulationOutputProcessingFunction) == type(None):
                self.mu_guess_SimulatedResponses = np.atleast_2d(self.mu_guess_SimulatedOutput)
            if type(simulationOutputProcessingFunction) != type(None):
                self.mu_guess_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.mu_guess_SimulatedOutput)     )
                
            #Get map simiulated output and simulated responses.
            self.map_SimulatedOutput = simulationFunction(self.map_parameter_set)           
            if type(simulationOutputProcessingFunction) == type(None):
                self.map_SimulatedResponses = np.atleast_2d(self.map_SimulatedOutput)
            if type(simulationOutputProcessingFunction) != type(None):
                self.map_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.map_SimulatedOutput)     )
            
            if hasattr(self, 'mu_AP_parameter_set'): #Check if a mu_AP has been assigned. It is normally only assigned if mcmc was used.           
                #Get mu_AP simiulated output and simulated responses.
                self.mu_AP_SimulatedOutput = simulationFunction(self.mu_AP_parameter_set)
                if type(simulationOutputProcessingFunction) == type(None):
                    self.mu_AP_SimulatedResponses = np.atleast_2d(self.mu_AP_SimulatedOutput)
                if type(simulationOutputProcessingFunction) != type(None):
                    self.mu_AP_SimulatedResponses =  np.atleast_2d(     simulationOutputProcessingFunction(self.mu_AP_SimulatedOutput)      )
                for responseDimIndex in range(self.UserInput.num_response_dimensions):
                    listOfYArrays = [self.UserInput.responses['responses_observed'][responseDimIndex], self.mu_guess_SimulatedResponses[responseDimIndex], self.map_SimulatedResponses[responseDimIndex], self.mu_AP_SimulatedResponses[responseDimIndex]]        
                    allResponsesListsOfYArrays.append(listOfYArrays)
            else: #Else there is no mu_AP.
                for responseDimIndex in range(self.UserInput.num_response_dimensions):
                    listOfYArrays = [self.UserInput.responses['responses_observed'][responseDimIndex], self.mu_guess_SimulatedResponses[responseDimIndex], self.map_SimulatedResponses[responseDimIndex]]        
                    allResponsesListsOfYArrays.append(listOfYArrays)
        if plot_settings == {}: 
            plot_settings = self.UserInput.simulated_response_plot_settings
            if hasattr(self, 'mu_AP_parameter_set'): 
                plot_settings['legendLabels'] = ['experiments',  'mu_guess', 'MAP','mu_AP']
            else: #Else there is no mu_AP.
                plot_settings['legendLabels'] = ['experiments',  'mu_guess', 'MAP']
            #Other allowed settings are like this, but will be fed in as simulated_response_plot_settings keys rather than plot_settings keys.
            #plot_settings['x_label'] = 'T (K)'
            #plot_settings['y_label'] = r'$rate (s^{-1})$'
            #plot_settings['y_range'] = [0.00, 0.025] #optional.
            #plot_settings['figure_name'] = 'tprposterior'
        if 'figure_name' not in plot_settings:
            plot_settings['figurename'] = 'Posterior'
        import plotting_functions as plotting_functions #This is going to become import CheKIPEUQ.plotting_functions as plotting_functions
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
                figureObject = plotting_functions.createSimulatedResponsesPlot(allResponses_x_values[0], allResponsesListsOfYArrays[responseDimIndex], plot_settings, listOfYUncertainties=allResponsesListsOfYUncertainties[responseDimIndex])
                np.savetxt(individual_plot_settings['figure_name']+".csv", np.vstack((allResponses_x_values[0], allResponsesListsOfYArrays[responseDimIndex])).transpose(), delimiter=",", header='x_values, observed, sim_initial_guess, sim_MAP, sim_mu_AP', comments='')
            if np.shape(allResponses_x_values)[0] > 1: #This means a separate abscissa for each response.
                figureObject = plotting_functions.createSimulatedResponsesPlot(allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex], plot_settings, listOfYUncertainties=allResponsesListsOfYUncertainties[responseDimIndex])
                np.savetxt(individual_plot_settings['figure_name']+".csv", np.vstack((allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex])).transpose(), delimiter=",", header='x_values, observed, sim_initial_guess, sim_MAP, sim_mu_AP', comments='')
            allResponsesFigureObjectsList.append(figureObject)
        return allResponsesFigureObjectsList  #This is a list of matplotlib.pyplot as plt objects.

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
        for parameterIndex, parameterName in enumerate(self.UserInput.parameterNamesAndMathTypeExpressionsDict):
            individual_model_parameter_dictionary = {'parameter_number': parameterIndex, 'parameter_name': self.UserInput.parameterNamesAndMathTypeExpressionsDict[parameterName]} #we are actually putting the MathTypeExpression as the parameter name when feeding to mum_pce.
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
                firstParameter = int(self.UserInput.parameterNamesAndMathTypeExpressionsDict[pairIndex[0]])
                secondParameter = int(self.UserInput.parameterNamesAndMathTypeExpressionsDict[pairIndex[0]])
                pairs_of_parameter_indices[pairIndex] = [firstParameter, secondParameter]        
        figureObject_beta.mumpce_plots(model_parameter_info = self.UserInput.model_parameter_info, active_parameters = self.UserInput.active_parameters, pairs_of_parameter_indices = pairs_of_parameter_indices, posterior_mu_vector = posterior_mu_vector, posterior_cov_matrix = posterior_cov_matrix, prior_mu_vector = np.array(self.UserInput.mu_prior), prior_cov_matrix = self.UserInput.covmat_prior, contour_settings_custom = self.UserInput.contour_settings_custom)
        return figureObject_beta

    def createAllPlots(self):


        try:
            self.makeHistogramsForEachParameter()    
            self.makeSamplingScatterMatrixPlot()



            self.createMumpcePlots()
        except: #TODO: do something better than try & accept. Right now, this is because the above plots are designed for mcmc sampling and don't work if pure grid search or pure optimize is used.
            pass
        self.createSimulatedResponsesPlots()


class verbose_optimization_wrapper: #Modified slightly From https://stackoverflow.com/questions/16739065/how-to-display-progress-of-scipy-optimize-function
    def __init__(self, function):
        self.f = function # actual objective function
        self.num_calls = 0 # how many times f has been called
        self.callback_count = 0 # number of times callback has been called, also measures iteration count
        self.list_calls_inp = [] # input of all calls
        self.list_calls_res = [] # result of all calls
        self.decreasing_list_calls_inp = [] # input of calls that resulted in decrease
        self.decreasing_list_calls_res = [] # result of calls that resulted in decrease
        self.list_callback_inp = [] # only appends inputs on callback, as such they correspond to the iterations
        self.list_callback_res = [] # only appends results on callback, as such they correspond to the iterations
    
    def simulate(self, x):
        """Executes the actual simulation and returns the result, while
        updating the lists too. Pass to optimizer without arguments or
        parentheses."""
        result = self.f(x) # the actual evaluation of the function
        if not self.num_calls: # first call is stored in all lists
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
            self.list_callback_inp.append(x)
            self.list_callback_res.append(result)
        elif result < self.decreasing_list_calls_res[-1]:
            self.decreasing_list_calls_inp.append(x)
            self.decreasing_list_calls_res.append(result)
        self.list_calls_inp.append(x)
        self.list_calls_res.append(result)
        self.num_calls += 1
        return result
    
    def callback(self, xk, *_):
        """Callback function that can be used by optimizers of scipy.optimize.
        The third argument "*_" makes sure that it still works when the
        optimizer calls the callback function with more than one argument. Pass
        to optimizer without arguments or parentheses."""
        
        s1 = "{0:4d}  ".format(self.callback_count)
        xk = np.atleast_1d(xk)
        # search backwards in input list for input corresponding to xk
        for i, x in reversed(list(enumerate(self.list_calls_inp))):
            x = np.atleast_1d(x)
            if np.allclose(x, xk):
                break
    
        for comp in xk:
            s1 += f"{comp:10.5e}\t"
        s1 += f"{self.list_calls_res[i]:10.5e}"
        self.list_callback_inp.append(xk)
        self.list_callback_res.append(self.list_calls_res[i])
    
        if not self.callback_count:
            s0 = "Iter  "
            for j, _ in enumerate(xk):
                tmp = f"Par-{j+1}"
                s0 += f"{tmp:10s}\t"
            s0 += "ObjectiveF"
            print(s0)
        print(s1)
        self.callback_count += 1

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
        
if __name__ == "__main__":
    pass

