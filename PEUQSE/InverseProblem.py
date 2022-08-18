import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
#import scipy
from scipy.stats import multivariate_normal
#from scipy.integrate import odeint
import sys
import time
import copy
from collections.abc import Iterable 

try:
    import CiteSoft
    from CiteSoft import function_call_cite
except:
    import os #The below lines are to allow CiteSoftLocal to be called regardless of user's working directory.
    lenOfFileName = len(os.path.basename(__file__)) #This is the name of **this** file.
    absPathWithoutFileName = os.path.abspath(__file__)[0:-1*lenOfFileName]
    sys.path.append(absPathWithoutFileName)
    import CiteSoftLocal as CiteSoft
try:
    import UnitTesterSG.nestedObjectsFunctions as nestedObjectsFunctions
except:
    import PEUQSE.nestedObjectsFunctionsLocal as nestedObjectsFunctions

class parameter_estimation:
    #Inside this class, a UserInput namespace is provided. This has dictionaries of UserInput choices.
    #However, the code initally parses those choices and then puts processed versions in the SAME name space, but no longer in the dictionaries.
    #So functions in this class should (when possible) call the namespace variables that are not in dictionaries, unless the original userinput is desired.
    #'inverse problem'. Initialize chain with initial guess (prior if not provided) as starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  
    
    
    software_name = "PEUQSE Bayesian Parameter Estimation"
    software_version = "1.0.0"
    software_unique_id = "https://doi.org/10.1002/cctc.202000953"
    software_kwargs = {"version": software_version, "author": ["Aditya Savara", "Eric A. Walker"], "doi": "https://doi.org/10.1002/cctc.202000953", "cite": "Savara, A. and Walker, E.A. (2020), PEUQSE Intro 1: Bayesian Parameter Estimation Considering Uncertainty or Error from both Experiments and Theory. ChemCatChem. Accepted. doi:10.1002/cctc.202000953"} 
    @CiteSoft.after_call_compile_consolidated_log()
    @CiteSoft.function_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def __init__(self, UserInput = None):
        #TODO: settings that are supposed to be Booleans should get Boolean cast in here. Otherwise if they are strings they will cause problems in "or" statements (where strings can return true even if the string is 'False').        
        self.UserInput = UserInput #Note that this is a pointer, so the later lines are within this object.
        #Now will automatically populate some variables from UserInput
        #make subdirectories as needed.
        import os
        for directoryName in UserInput.directories:
            if not os.path.exists(directoryName):
                os.makedirs(directoryName)

        #Populate variables for bounds and reduced parameter space.
        self.parameterBoundsOn = bool(len(UserInput.model['InputParameterPriorValues_upperBounds']) + len(UserInput.model['InputParameterPriorValues_lowerBounds']))
        UserInput.InputParameterPriorValues_lowerBounds = UserInput.model['InputParameterPriorValues_lowerBounds']
        UserInput.InputParameterPriorValues_upperBounds = UserInput.model['InputParameterPriorValues_upperBounds']
        self.reducedParameterSpaceOn = bool(len(UserInput.model['reducedParameterSpace']))
        
        #Check for deprecated UserInput choices.
        if hasattr(UserInput.parameter_estimation_settings, 'multistart_gridsearchToSamples'):
            print("The UserInput feature parameter_estimation_settings['multistart_gridsearchToSamples'] has been renamed. Use parameter_estimation_settings['multistart_permutationsToSamples'].")
            if hasattr(UserInput.parameter_estimation_settings, 'multistart_permutationsToSamples') == False:    
                UserInput.parameter_estimation_settings['multistart_permutationsToSamples']= UserInput.parameter_estimation_settings['multistart_gridsearchToSamples']
        if hasattr(UserInput.parameter_estimation_settings, 'multistart_gridsearch_threshold_filter_coefficient'):
            print("The UserInput feature parameter_estimation_settings['multistart_gridsearch_threshold_filter_coefficient'] has been renamed. Use parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient'].")
            if hasattr(UserInput.parameter_estimation_settings, 'multistart_permutationsToSamples_threshold_filter_coefficient') == False:    
                UserInput.parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient']= UserInput.parameter_estimation_settings['multistart_gridsearch_threshold_filter_coefficient']
        #Check if EnsembleJumpSampling specification for multistart and change type to EnsembleModifiedMHSampling
        if UserInput.parameter_estimation_settings['multistart_searchType'] == 'doEnsembleJumpSampling':
            self.UserInput.parameter_estimation_settings['multistart_searchType'] = 'doEnsembleModifiedMHSampling'

        #Check if there are parameterNames provided. If not, we will make some.
        if len(UserInput.model['parameterNamesAndMathTypeExpressionsDict']) == 0:
            numParameters = len(UserInput.model['InputParameterPriorValues'])
            for parameterIndex in range(0,numParameters):
                UserInput.model['parameterNamesAndMathTypeExpressionsDict'][str(parameterIndex)]= 'ParInd_'+str(parameterIndex)
        elif type(UserInput.model['parameterNamesAndMathTypeExpressionsDict']) == type([1]): # if it's a list, make a dictionary.
                listToMakeDictionary = UserInput.model['parameterNamesAndMathTypeExpressionsDict']
                newDictionary = {}
                for paramName in listToMakeDictionary:
                    newDictionary[str(paramName)] = str(paramName)
                UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = newDictionary
        UserInput.parameterNamesList = list(UserInput.model['parameterNamesAndMathTypeExpressionsDict'].keys())
        UserInput.stringOfParameterNames = str(UserInput.parameterNamesList).replace("'","")[1:-1]
        UserInput.parameterNamesAndMathTypeExpressionsDict = UserInput.model['parameterNamesAndMathTypeExpressionsDict']
        if self.UserInput.parameter_estimation_settings['verbose']: 
            print("Parameter Estimation Object Initialized")
        
        
        #Set self.UserInput.user_requested_convergence_diagnostics to that from the UserInput. This may be changed to true or false during multi-start etc.
        self.UserInput.user_requested_convergence_diagnostics = self.UserInput.parameter_estimation_settings['convergence_diagnostics']
        
        if type(UserInput.parameter_estimation_settings['checkPointFrequency']) != type(None): #This is for backwards compatibility.
            UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] = UserInput.parameter_estimation_settings['checkPointFrequency']
            UserInput.parameter_estimation_settings['multistart_checkPointFrequency'] = UserInput.parameter_estimation_settings['checkPointFrequency']
        UserInput.request_mpi = False #Set as false as default.
        if ( \
            UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or \
            UserInput.parameter_estimation_settings['multistart_parallel_sampling'] or \
            UserInput.doe_settings['parallel_conditions_exploration'] or  \
            UserInput.doe_settings['parallel_parameter_modulation'] \
            ) \
            == True:
            UserInput.request_mpi = True
            if (UserInput.doe_settings['parallel_conditions_exploration'] or UserInput.doe_settings['parallel_parameter_modulation']) and (UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or UserInput.parameter_estimation_settings['multistart_parallel_sampling']):
                print("Warning: Parallelization of Design of experiments is not compatible with parallelization of either mcmc_parallel_sampling or multistart_parallel_sampling.  Those other features are being turned off.")
                UserInput.parameter_estimation_settings['multistart_parallel_sampling'] = False
                UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] = False
                
        if UserInput.request_mpi == True: #Rank zero needs to clear out the mpi_cached_files directory (unless we are continuing sampling), so check if we are using rank 0.
            import os; #import sys
            import PEUQSE.parallel_processing
            if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                if not os.path.exists('./mpi_cached_files'):
                    os.makedirs('./mpi_cached_files')
                if not os.path.exists(UserInput.directories['graphs']+"mpi_cached_files"):
                    os.makedirs(UserInput.directories['graphs']+"mpi_cached_files")                
                if not os.path.exists(UserInput.directories['logs_and_csvs']+"mpi_cached_files"):
                    os.makedirs(UserInput.directories['logs_and_csvs']+"mpi_cached_files")                
                if not os.path.exists(UserInput.directories['pickles']+"mpi_cached_files"):
                    os.makedirs(UserInput.directories['pickles']+"mpi_cached_files")                
                if ('mcmc_continueSampling' not in UserInput.parameter_estimation_settings) or (UserInput.parameter_estimation_settings['mcmc_continueSampling'] == False) or (UserInput.parameter_estimation_settings['mcmc_continueSampling'] == 'auto'):
                    os.chdir("./mpi_cached_files")                
                    deleteAllFilesInDirectory()
                    os.chdir("./..")
                    if UserInput.directories['graphs'] != "./":
                        os.chdir(UserInput.directories['graphs']+"mpi_cached_files")
                        deleteAllFilesInDirectory()
                        os.chdir("./../..")
                    if UserInput.directories['logs_and_csvs'] != "./":
                        os.chdir(UserInput.directories['graphs']+"mpi_cached_files")
                        deleteAllFilesInDirectory()
                        os.chdir("./../..")
                    if UserInput.directories['pickles'] != "./":
                        os.chdir(UserInput.directories['graphs']+"mpi_cached_files")
                        deleteAllFilesInDirectory()
                        os.chdir("./../..")
                
                #Now check the number of processor ranks to see if the person really is using parallel processing.
                if PEUQSE.parallel_processing.numProcessors > 1:    #This is the normal case.
                    sys.exit() #TODO: right now, processor zero just exits after making and emptying the directory. In the future, things will be more complex for the processor zero.
                elif PEUQSE.parallel_processing.numProcessors == 1: #This is the case where the person has only one process rank, so probably does not want code execution to stop just yet. (This is an intentional case for gridsearch for example, where running without mpi will print the number of grid Permutations).
                    print("Notice: you have requested parallel processing by MPI but have only 1 processor rank enabled or are not using mpi for this run. Parallel processing is being disabled for this run. If you are running to find the number of process ranks to use, another message will be printed out with the number of processor ranks to provide to mpi.")
                    UserInput.request_mpi = False
                    if UserInput.parameter_estimation_settings['mcmc_parallel_sampling']:
                        print("Your settings suggest that you are trying to use mcmc_parallel_sampling. Please use the mpi command from the prompt.  To do N parallel samplings requires N+1 process ranks. For example, if you wanted to have 4 parallel samplings, you would need 5 process ranks and would use: mpiexec -n 5 python runfile_for_your_analysis.py")
                        sys.exit()
                    
        
        #Setting this object so that we can make changes to it below without changing UserInput dictionaries.
        self.UserInput.mu_prior = np.array(UserInput.model['InputParameterPriorValues'], dtype='float')
        #Below code is mainly for allowing uniform distributions in priors.
        UserInput.InputParametersPriorValuesUncertainties = np.array(UserInput.model['InputParametersPriorValuesUncertainties'],dtype='float') #Doing this so that the -1.0 check below should work.
        if -1.0 in UserInput.InputParametersPriorValuesUncertainties: #This means that at least one of the uncertainties has been set to "-1" which means a uniform distribution. 
            UserInput.InputParametersPriorValuesUniformDistributionsIndices = [] #intializing.
            if len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) != 1:
                print("A value of '-1' in the uncertainties signifies a uniform distribution for PEUQSE. As of July 1st 2020, the uniform distribution feature is only compatible with a 1D of array for uncertainties and not compatible with providing a full covariance matrix. If you need such a feature, contact the developers because it could be implemented. Eventually, a more sophisiticated back end may be used which would allow such a feature.")
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
                    if self.parameterBoundsOn == False:
                        print("ERROR: An uncertaintyValue of -1.0 has been provided which indicates a Uniform distribution. Uniform distributions require both upper and lower bounds, but these have not been provided. CheKiPUEQ is exiting the program."); sys.exit()
                    std_prior_single_parameter = (UserInput.InputParameterPriorValues_upperBounds[parameterIndex] - UserInput.InputParameterPriorValues_lowerBounds[parameterIndex])/(12**0.5)
                    UserInput.InputParametersPriorValuesUncertainties[parameterIndex] = std_prior_single_parameter #Note that going forward the array InputParametersPriorValuesUncertainties cannot be checked to see if the parameter is from a uniform distribution. Instead, InputParametersPriorValuesUniformDistributionsKey must be checked. 
                    #We will also fill the model['InputParameterPriorValues'] to have the mean of the two bounds. This can matter for some of the scaling that occurs later.
                    self.UserInput.mu_prior[parameterIndex] = (UserInput.InputParameterPriorValues_upperBounds[parameterIndex] + UserInput.InputParameterPriorValues_lowerBounds[parameterIndex])/2
        
        #Now to make covmat. Leaving the original dictionary object intact, but making a new object to make covmat_prior.
        if len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) == 1 and (len(UserInput.InputParametersPriorValuesUncertainties) > 0): #If it's a 1D array/list that is filled, we'll diagonalize it.
            UserInput.std_prior = np.array(UserInput.InputParametersPriorValuesUncertainties, dtype='float') #using 32 since not everyone has 64.
            UserInput.var_prior = np.power(UserInput.InputParametersPriorValuesUncertainties,2)
            UserInput.covmat_prior = np.diagflat(self.UserInput.var_prior) 
        elif len(np.shape(UserInput.InputParametersPriorValuesUncertainties)) > 1: #If it's non-1D, we assume it's already a covariance matrix.
            UserInput.covmat_prior = np.array(UserInput.InputParametersPriorValuesUncertainties, dtype='float')
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
        
        #Below, we are generating samples of the prior for info gain purposes.  This requires considering random seeds.
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if type(self.UserInput.parameter_estimation_settings['mcmc_random_seed']) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        self.samples_of_prior = np.random.multivariate_normal(self.UserInput.mu_prior,UserInput.covmat_prior,UserInput.parameter_estimation_settings['mcmc_length'])
        
        #Now do some processing on the responses formatting and uncertainties.
        #Make them 2dNested if needed..

        UserInput.responses_observed = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses['responses_observed']))
        if UserInput.responses['num_responses']=='auto':
            self.UserInput.num_response_dimensions = np.shape(UserInput.responses_observed)[0]
        else:
            self.UserInput.num_response_dimensions = UserInput.responses['num_responses']
        if len(UserInput.responses['responses_abscissa']) == 0: #This means it has not been provided and we will make one.
            UserInput.responses_abscissa = [] #the one from input should already be a list, but we start a fresh one.
            for responseIndex in range(0,self.UserInput.num_response_dimensions):
                numPoints = len(UserInput.responses_observed[responseIndex])
                UserInput.responses_abscissa.append(np.linspace(0, numPoints,numPoints))
        else:
            UserInput.responses_abscissa = UserInput.responses['responses_abscissa']
        UserInput.responses_abscissa = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses_abscissa))        
        #Make sure all objects inside are arrays (if they are lists we convert them). This is needed to apply the heurestic.
        UserInput.responses_abscissa = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.responses_abscissa)
        UserInput.responses_observed = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.responses_observed)

        #Now to process responses_observed_uncertainties, there are several options so we need to process it according to the cases.
        #The normal case:
        if isinstance(self.UserInput.responses['responses_observed_uncertainties'], Iterable): #If it's an array or like one, we take it as is. The other options are a none object or a function.
            UserInput.responses_observed_uncertainties = UserInput.responses['responses_observed_uncertainties']
            #Processing of responses_observed_uncertainties for case that a blank list is received and not zeros.
            if len(UserInput.responses['responses_observed_uncertainties']) == 0:
                #if the response uncertainties is blank, we will use the heurestic of sigma = 5% of the observed value, and then add an orthogonal uncertainty of 2% of the maximum for that response. 
                #Note that we are actually checking in index[0], that is because as an atleast_2d array even a blank list / array in it will give a length of 1.
                UserInput.responses_observed_uncertainties = np.abs( UserInput.responses_observed) * 0.05
                for responseIndex in range(0,UserInput.num_response_dimensions): #need to cycle through to apply the "minimum" uncertainty of 0.02 times the max value.
                    maxResponseAbsValue = np.max(np.abs(UserInput.responses_observed[responseIndex]))
                    UserInput.responses_observed_uncertainties[responseIndex] = ( UserInput.responses_observed_uncertainties[responseIndex]**2 + (maxResponseAbsValue*0.02)**2 ) ** 0.5
                    #The below deprecated syntax is a bit hard to read, but it is similar to this: a[a==2] = 10 #replace all 2's with 10's                    #UserInput.responses_observed_uncertainties[responseIndex][UserInput.responses_observed_uncertainties[responseIndex] < maxResponseAbsValue * 0.02] = maxResponseAbsValue * 0.02
            elif nestedObjectsFunctions.sumNested(UserInput.responses['responses_observed_uncertainties']) == 0: #If a 0 (or list summing to 0) is provided, we will make the uncertainties zero.
                UserInput.responses_observed_uncertainties = UserInput.responses_observed * 0.0 #This will work because we've converted the internals to array already.
                #Below two lines not needed. Should be removed if everythig is working fine after Nov 2020.
                #for responseIndex in range(0,len(UserInput.responses_observed[0])):
                #    UserInput.responses_observed_uncertainties[0][responseIndex]= UserInput.responses_observed[0][responseIndex]*0.0
            UserInput.responses_observed_uncertainties = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses_observed_uncertainties))
            UserInput.responses_observed_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.responses_observed_uncertainties)
            #If the feature of self.UserInput.responses['responses_observed_weighting'] has been used, then we need to apply that weighting to the uncertainties.
            if len(self.UserInput.responses['responses_observed_weighting']) > 0:
                UserInput.responses_observed_weighting = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses['responses_observed_weighting']))
                UserInput.responses_observed_weighting = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.responses_observed_weighting)
                UserInput.responses_observed_weighting = UserInput.responses_observed_weighting.astype(np.float)
                UserInput.responses_observed_weight_coefficients = copy.deepcopy(UserInput.responses_observed_weighting).astype(np.float) #initialize the weight_coefficients
                #We'll apply it 1 response at a time.
                for responseIndex, responseWeightingArray in enumerate(UserInput.responses_observed_weighting):
                    if 0 in responseWeightingArray: #we can't have zeros in weights. So if we have any zeros, we will set the weighting of those to 1E6 times less than other values.
                        #Originally, used minNonZero/1E6. Now, use eps which is the smallest non-zero value allowed.
                        #minNonZero = np.min(UserInput.responses_observed_weighting[UserInput.responses_observed_weighting>0])
                        responseWeightingArray[responseWeightingArray==0] = np.finfo(float).eps #minNonZero/1E6  #set the 0 values to be 1E6 times smaller than minNonZero.
                        UserInput.responses_observed_weighting[responseIndex] = responseWeightingArray 
                #now calculate and apply the weight coefficients.
                for responseIndex in range(len(UserInput.responses_observed_weighting)):
                    UserInput.responses_observed_weight_coefficients[responseIndex] = (UserInput.responses_observed_weighting[responseIndex])**(-0.5) #this is analagous to the sigma of a variance weighted heuristic.
                UserInput.responses_observed_uncertainties = UserInput.responses_observed_uncertainties*UserInput.responses_observed_weight_coefficients
        else: #The other possibilities are a None object or a function. For either of thtose cases, we simply set UserInput.responses_observed_uncertainties equal to what the user provided.
            UserInput.responses_observed_uncertainties = copy.deepcopy(self.UserInput.responses['responses_observed_uncertainties'])

        #Now to process responses_simulation_uncertainties, there are several options so we need to process it according to the cases.
        #The normal case:
        if isinstance(self.UserInput.model['responses_simulation_uncertainties'], Iterable): #If it's an array or like one, we take it as is. The other options are a none object or a function.
            UserInput.responses_simulation_uncertainties = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.model['responses_simulation_uncertainties']))
            UserInput.responses_simulation_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.responses_simulation_uncertainties)
            #TODO: allow a length of zero 'responses_simulation_uncertainties' to mean that the heurestic function should be called after each simulation. That is not what it is doing right now (Dec 2020). Right now it is just using a static value from the heurestic applied to the observed_responses.
            #Processing of responses_simulation_uncertainties for case that a blank list is received and not zeros.
            if len(UserInput.responses_simulation_uncertainties[0]) == 0: 
                #if the response uncertainties is blank, we will use the heurestic of sigma = 5% of the observed value, with a floor of 2% of the maximum for that response. 
                #Note that we are actually checking in index[0], that is because as an atleast_2d array even a blank list / array in it will give a length of 1.
                UserInput.responses_simulation_uncertainties = np.abs(UserInput.responses_observed) * 0.05
                for responseIndex in range(0,UserInput.num_response_dimensions): #need to cycle through to apply the "minimum" uncertainty of 0.02 times the max value.
                    maxResponseAbsValue = np.max(np.abs(UserInput.responses_observed[responseIndex])) #Because of the "at_least2D" we actually need to use index 0.
                    #The below syntax is a bit hard to read, but it is similar to this: a[a==2] = 10 #replace all 2's with 10's
                    UserInput.responses_simulation_uncertainties[responseIndex][UserInput.responses_simulation_uncertainties[responseIndex] < maxResponseAbsValue * 0.02] = maxResponseAbsValue * 0.02
            elif nestedObjectsFunctions.sumNested(UserInput.responses_simulation_uncertainties) == 0: #If a 0 (or list summing to 0) is provided, we will make the uncertainties zero.
                UserInput.responses_simulation_uncertainties = UserInput.responses_observed * 0.0 #This will work because we've converted the internals to array already.
        else: #The other possibilities are a None object or a function. For either of thtose cases, we simply set UserInput.responses_simulation_uncertainties equal to what the user provided.
            UserInput.responses_simulation_uncertainties = copy.deepcopy(self.UserInput.model['responses_simulation_uncertainties'])

        #Now to process simulatedResponses_upperBounds and simulatedResponses_lowerBounds. Can be a blank list or a nested list.
        if len(UserInput.model['simulatedResponses_upperBounds']) > 0:
            UserInput.model['simulatedResponses_upperBounds'] = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.model['simulatedResponses_upperBounds']))
            UserInput.model['simulatedResponses_upperBounds'] = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.model['simulatedResponses_upperBounds'])
        if len(UserInput.model['simulatedResponses_lowerBounds']) > 0:
            UserInput.model['simulatedResponses_lowerBounds'] = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.model['simulatedResponses_lowerBounds']))
            UserInput.model['simulatedResponses_lowerBounds'] = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(UserInput.model['simulatedResponses_lowerBounds'])
            
        self.UserInput.num_data_points = len(nestedObjectsFunctions.flatten_2dNested(UserInput.responses_observed)) #This works if there is a single response series.
        #We need to figure out if the abscissa has length equal to the responses or not.
        if len(UserInput.responses_abscissa) == len(UserInput.responses_observed):
            self.separate_abscissa_per_response = True  #This means we **will** iterate across the abscissa when iterating across each response.
        else:
            self.separate_abscissa_per_response = False #This means we **won't** iterate across the abscissa when iterating across each response.

        self.staggeredResponses = nestedObjectsFunctions.checkIfStaggered_2dNested(UserInput.responses_observed)
        #TODO: This currently only is programmed for if the uncertainties are uncorrelated standard deviaions (so is not compatible with a directly fed cov_mat). Also, we need to figure out what do when they are not gaussian/symmetric.
        UserInput.responses_observed_transformed, UserInput.responses_observed_transformed_uncertainties  = self.transform_responses(UserInput.responses_observed, UserInput.responses_observed_uncertainties) #This creates transforms for any data that we might need it. The same transforms will also be applied during parameter estimation.
            
        #The below unusual code is because during doeParameterModulationPermutationsScanner, populate synthetic data calls init again.
        #So we will only call populateIndependentVariablesFunction if we're not in the middle of design of experiments.
        if not hasattr(self, 'middle_of_doe_flag'): #We check of the middle_of_doe_flag exists. #If the flag is not there and the populate function exists, we call it.
            if type(UserInput.model['populateIndependentVariablesFunction']) != type(None):
                UserInput.model['populateIndependentVariablesFunction'](UserInput.responses['independent_variables_values']) 
        if hasattr(self, 'middle_of_doe_flag'): #We check of the middle_of_doe_flag exists. If it's there, no problem.
            if self.middle_of_doe_flag == False: #If the flag is there, we only proceed to call the function if the flag is set to false.
                if type(UserInput.model['populateIndependentVariablesFunction']) != type(None):
                    UserInput.model['populateIndependentVariablesFunction'](UserInput.responses['independent_variables_values']) 
        if hasattr(self, 'middle_of_doe_flag'): #if the middle_of_doe_flag is there and true, then we also turn off the convergence diagnostics.
            if self.middle_of_doe_flag == True:
                self.UserInput.user_requested_convergence_diagnostics = False
                
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
            if np.shape(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.scaling_uncertainties))[0]==1:  #In this case, the uncertainties is not a covariance matrix.
                pass
            elif np.shape(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.scaling_uncertainties))[0]==np.shape(nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.scaling_uncertainties))[1]: #In his case, the uncertainties are a covariance matrix so we take the diagonal (which are variances) and the square root of them.
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
        #First, we have to make sure self.UserInput.responses_observed_transformed_uncertainties is an iterable. It could be a none-type or a function.
        if isinstance(self.UserInput.responses_observed_transformed_uncertainties, Iterable):
            self.observed_responses_covmat_transformed = returnShapedResponseCovMat(self.UserInput.num_response_dimensions, self.UserInput.responses_observed_transformed_uncertainties)
        else: #If responses_observed_transformed_uncertainties is a None type, then we don't need observed_responses_covmat_transformed.  If it is a function, then we have to create the object on the fly so can't create it now.
            pass

        #self.covmat_prior = UserInput.covmat_prior
        self.Q_mu = self.UserInput.mu_prior*0 # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.  
        self.Q_covmat = self.UserInput.covmat_prior # Take small steps. 
        #TODO: Make initial guess handle a string, then assume that it is a pickle file
        # check if InputParameterInitialGuess is a string
        # if so, read the pickle file in the same way as the sampling functions (look for file, then read)
        # After reading it, slice the matrix to only the first walker is saved.
        if isinstance(self.UserInput.model['InputParameterInitialGuess'], str):
            # try to read from current directory and pickle directory
            start_point_pkl_file_name = self.UserInput.model['InputParameterInitialGuess']
            if '.pkl' in start_point_pkl_file_name: # remove '.pkl' from string if it is there
                start_point_pkl_file_name = start_point_pkl_file_name.replace('.pkl', '')
            from os.path import exists
            # check if file exists in current working directory
            if exists(start_point_pkl_file_name + '.pkl'):
                initialGuessUnfiltered = unpickleAnObject(start_point_pkl_file_name)
            # check if file exists in pickle directory
            elif exists(self.UserInput.directories['pickles'] + start_point_pkl_file_name + '.pkl'):
                initialGuessUnfiltered = unpickleAnObject(self.UserInput.directories['pickles'] + start_point_pkl_file_name)
            else:
                print('The pickled object for initial guess points must exist in base directory or pickles directory. The pickled file should have extension of ".pkl"')
                sys.exit()
            # make sure the initial guess is a numpy array
            initialGuessUnfiltered = np.array(initialGuessUnfiltered)
            # check if the shape of the loaded initial guess is a single point or a list of walkers
            if len(initialGuessUnfiltered.shape) == 1:
                self.UserInput.InputParameterInitialGuess = initialGuessUnfiltered
            elif len(initialGuessUnfiltered.shape) == 2:
                # this is assumed to have the shape (num_walkers, num_parameters)
                self.UserInput.InputParameterInitialGuess = initialGuessUnfiltered[0] # take only the first walker as the initial points
            else:
                print('The shape of the initial guess pickled array should be (num_walkers, num_parameters). The current pickled array does not have this shape and has:', initialGuessUnfiltered.shape)
                sys.exit()
        else:
            #Getting initial guess of parameters and populating the internal variable for it.
            if ('InputParameterInitialGuess' not in self.UserInput.model) or (len(self.UserInput.model['InputParameterInitialGuess'])== 0): #if an initial guess is not provided, we use the prior.
                self.UserInput.model['InputParameterInitialGuess'] = np.array(self.UserInput.mu_prior, dtype='float')
            #From now, we switch to self.UserInput.InputParameterInitialGuess because this is needed in case we're going to do reducedParameterSpace or grid sampling.
            self.UserInput.InputParameterInitialGuess = np.array(self.UserInput.model['InputParameterInitialGuess'], dtype='float')
            # reassure the initial guess is a numpy array
            # check the shape to make sure it is a 1D array
            self.UserInput.InputParameterInitialGuess = np.array(self.UserInput.InputParameterInitialGuess)
            if len(self.UserInput.InputParameterInitialGuess.shape) > 1:
                if self.UserInput.InputParameterInitialGuess.shape[0] == 1: # make sure that the array was not double nested by accident
                    self.UserInput.InputParameterInitialGuess = np.ndarray.flatten(self.UserInput.InputParameterInitialGuess)
                else:
                    print('The Initial guess must be either a 1D array or pickle file string.')
                    sys.exit()
        #Now populate the simulation Functions. #NOTE: These will be changed if a reduced parameter space is used.
        self.UserInput.simulationFunction = self.UserInput.model['simulateByInputParametersOnlyFunction']
        self.UserInput.simulationOutputProcessingFunction = self.UserInput.model['simulationOutputProcessingFunction']

        #Check the shapes of the arrays for UserInput.responses_observed and UserInput.responses_observed_uncertainties by doing a simulation. Warn the user if the shapes don't match.
        initialGuessSimulatedResponses = self.getSimulatedResponses(self.UserInput.InputParameterInitialGuess)
        if np.shape(initialGuessSimulatedResponses) != np.shape(UserInput.responses_observed):
            print("PEUQSE Warning: the shape of the responses_observed is", np.shape(UserInput.responses_observed), ", but the shape using your provided simulation function is", np.shape(initialGuessSimulatedResponses), " .  PEUQSE is probably going to crash when trying to calculate the likelihood.")

        #Now reduce the parameter space if requested by the user. #Considered having this if statement as a function called outside of init.  However, using it in init is the best practice since it forces correct ordering of reduceParameterSpace and reduceResponseSpace
        if len(self.UserInput.model['reducedParameterSpace']) > 0:
            print("Notification: the UserInput.model['reducedParameterSpace'] is not blank. That means the only parameters allowed to change will be the ones in the indices inside 'reducedParameterSpace'.   All others will be held constant.  The values inside  'InputParameterInitialGuess will be used', and 'InputParameterPriorValues' if an initial guess was not provided.")
            self.reduceParameterSpace()
    
        #Now reduce the parameter space if requested by the user. #Considered having this if statement as a function called outside of init.  However, using it in init is the best practice since it forces correct ordering of reduceParameterSpace and reduceResponseSpace
        #This code must be **after** the reduceParameterSpace because this makes a wrapper for the simulationOutputProcessingFunction
        if len(self.UserInput.responses['reducedResponseSpace']) > 0:
            print("Important: the UserInput.model['reducedResponseSpace'] is not blank. That means the only responses examined will be the ones in the indices inside 'reducedReponseSpace'.   The values of all others will be discarded during each simulation.")
            self.reduceResponseSpace()
        
        #Check if we should plan to split the responses into separate likelihood terms. 
        self.prepareResponsesForSplitLikelihood = False #initialized as false, will change to True if needed.
        #First check if it is a single response, then if so a singlepoint with the initial values and see if it turns out okay.
        if len(self.UserInput.responses_observed_transformed)==1: #this means single response
            initialLogP = self.getLogP(self.UserInput.InputParameterInitialGuess)
            if (np.isnan(initialLogP) or initialLogP < -1E90):    
                if np.shape(UserInput.responses_observed) == np.shape(UserInput.responses_observed_uncertainties):
                    #this if statement only occurs if uncertainties are standard deviations and not covmat, 
                    #which means we can prepare the responses for split likelihood and should plan to.
                    self.prepareResponsesForSplitLikelihood = True
        
        self.permutation_and_doOptimizeNegLogP = False #just initializing this flag with its default.
        self.permutation_and_doOptimizeSSR = False #just initializing this flag with its default.
        self.permutation_and_doOptimizeLogP = False
    
    #fills the samples etc. By default, for multistart, it will also perform filtering and then create the post_burn_in_samples.
    def reload_samples(self, sampling_type='', filepath = ''):
        if (sampling_type != 'multistart') and (sampling_type != 'mcmc'):
            print("ERROR: reload_samples requires specifying either 'multistart' or 'mcmc' as the first argument"); sys.exit()
        if filepath == '':
            filepath = self.UserInput.directories['logs_and_csvs'] #take the default.
        #in both cases, multstart or mcmc, it's really an array of logP_and_parameter_values, just slightly different meanings.
        if sampling_type == 'multistart':  #load from the unfiltered values.
            self.permutations_MAP_logP_and_parameters_values = np.genfromtxt(filepath + "\\multistart_MAP_logP_and_parameters_values.csv", delimiter=",")
            multistart_permutationsToSamples_threshold_filter_coefficient = self.UserInput.parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient']
            permutations_to_samples_with_logP = convertPermutationsToSamples(self.permutations_MAP_logP_and_parameters_values, relativeFilteringThreshold = 10**(-1*multistart_permutationsToSamples_threshold_filter_coefficient))
            self.post_burn_in_samples = permutations_to_samples_with_logP[:,1:] #drop the first column which is logP.
            logP_and_parameter_values = np.array(permutations_to_samples_with_logP)
            self.post_burn_in_log_posteriors_un_normed_vec = permutations_to_samples_with_logP[:,0]
            self.post_burn_in_log_posteriors_un_normed_vec = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.post_burn_in_log_posteriors_un_normed_vec)).transpose()
            #need to populate post_burn_in_log_priors_vec this with an object, otherwise calculatePostBurnInStatistics will try to calculate all the priors.
            self.post_burn_in_log_priors_vec = None
            #Below is needed to avoid causing an error in the calculatePostBurnInStatistics since we don't have a real priors vec.
            self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_samples'] = False
            self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = False)    
        if sampling_type == 'mcmc':
            mcmc_logP_and_parameters_values = np.genfromtxt(filepath + "\\mcmc_logP_and_parameter_samples.csv", delimiter=",")
            logP_and_parameter_values = np.array(mcmc_logP_and_parameters_values)
            try:
                self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True)
            except:
                self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = False)
        # below are not needed because they occur in self.calculatePostBurnInStatistics
        # self.map_logP = max(np.array(logP_and_parameter_values[:,0]))
        # index_of_map_logP = np.where(logP_and_parameter_values[:,0] == self.map_logP)
        # self.map_parameter_set = logP_and_parameter_values[index_of_map_logP,1:]
        
        
    
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
        if self.parameterBoundsOn: #only need to reduce the iterables for the parameter bounds if they exist.
            UserInput.InputParameterPriorValues_lowerBounds = returnReducedIterable( UserInput.InputParameterPriorValues_lowerBounds    , reducedIndices )
            UserInput.InputParameterPriorValues_upperBounds = returnReducedIterable( UserInput.InputParameterPriorValues_upperBounds    , reducedIndices )
        #We have a list that keeps track of InputParameters which have Uniform Prior distributions. The logic for how to update that is a bit different.
        if hasattr(self.UserInput, 'InputParametersPriorValuesUniformDistributionsIndices') == True:
            NewInputParametersPriorValuesUniformDistributionsIndices = []
            for indexNumber, indexValue in enumerate(reducedIndices): #we need to know which position it has in the reducedIndices, and we'll append **that** position if the item is a uniform distribution.
                if indexValue in UserInput.InputParametersPriorValuesUniformDistributionsIndices:
                    NewInputParametersPriorValuesUniformDistributionsIndices.append(indexNumber)
                #else pass.
            UserInput.InputParametersPriorValuesUniformDistributionsIndices = NewInputParametersPriorValuesUniformDistributionsIndices
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
        if isinstance(self.UserInput.responses_simulation_uncertainties, Iterable): #If it's an array or like one, we take it as is. The other options are a non object or a function.
            responses_simulation_uncertainties = np.array(self.UserInput.responses_simulation_uncertainties)*1.0
        elif type(self.UserInput.responses_simulation_uncertainties) == type(None):
            responses_simulation_uncertainties = self.UserInput.responses_simulation_uncertainties
        else:  #Else we assume it's a function taking the discreteParameterVector.
            responses_simulation_uncertainties = self.UserInput.responses_simulation_uncertainties(discreteParameterVector) #This is passing an argument to a function.
            responses_simulation_uncertainties = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(responses_simulation_uncertainties))
            responses_simulation_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(responses_simulation_uncertainties)
        return responses_simulation_uncertainties #this normally a list of list where each element at the deepest level is 1 standard deviation of uncertainty (not variance).

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
        if type(simulationOutput) == type(None):
            return None #This is for the case that the simulation fails. User can have simulationOutput return a None type in case of failure. Perhaps should be made better in future. 
        if len(simulationOutput) == 1: #if a 'nan' is returned, then we'll treat the simulation like a failed simulation.
            if simulationOutput == np.float('nan'):
                return None 

            
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput #Is this the log of the rate? If so, Why?
        if type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        
        simulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(simulatedResponses)
        #This is not needed:
        #observedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.responses_observed)
        return simulatedResponses
    
    def transform_responses(self, nestedAllResponsesArray, nestedAllResponsesUncertainties = []):
        nestedAllResponsesArray_transformed = copy.deepcopy(nestedAllResponsesArray) #First make a copy to populate with transformed values.
        nestedAllResponsesUncertainties_transformed = copy.deepcopy(nestedAllResponsesUncertainties) #First make a copy to populate with transformed values. If blank, we won't populate it.        
        UserInput = self.UserInput
        #TODO: Make little function for interpolation in case it's necessary (see below).
#        def littleInterpolator():
#            abscissaRange = UserInput.responses_abscissa[responseIndex][-1] - UserInput.responses_abscissa[responseIndex][0] #Last value minus first value.
#            UserInput.responses_observed = nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses_observed)
#            UserInput.responses_observed_uncertainties = nestedObjectsFunctions.makeAtLeast_2dNested(UserInput.responses_observed_uncertainties)
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
                    if UserInput.responses['response_data_types'][responseIndex] == 'c':
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        #Perform the littleEuler twice.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray_transformed[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties_transformed[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'R':	#For reactant
                    if UserInput.responses['response_data_types'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        #TODO: use responses['points_if_transformed'] variable to interpolate the right number of points. This is for data that's not already evenly spaced.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'P':	 #For product
                    
                    if UserInput.responses['response_data_types'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        #TODO: use responses['points_if_transformed'] variable to interpolate the right number of points. This is for data that's not already evenly spaced.
                        t_values, nestedAllResponsesArray_transformed[responseIndex], dydt_values = littleEulerGivenArray(0, UserInput.responses_abscissa[abscissaIndex], nestedAllResponsesArray[responseIndex])
                        if len(nestedAllResponsesUncertainties) > 0:
                            nestedAllResponsesUncertainties_transformed[responseIndex] = littleEulerUncertaintyPropagation(nestedAllResponsesUncertainties[responseIndex], UserInput.responses_abscissa[abscissaIndex], np.mean(nestedAllResponsesUncertainties[responseIndex])/10) 
                if UserInput.responses['response_types'][responseIndex] == 'O': #O is for other.
                    if UserInput.responses['response_data_types'][responseIndex] == 'o': #other
                        pass
                    if UserInput.responses['response_data_types'][responseIndex] == 'c': #concentration
                        LittleEuler
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        LittleEulerTwice
        if UserInput.responses['data_overcategory'] == 'steady_state_kinetics': #TODO: so far, this does not do anything. It assumes that the abscissa is never time.
            for responseIndex, response in enumerate(UserInput.responses_observed):
                if UserInput.responses['response_types'][responseIndex] == 'T':	 #For abscissa of temperature dependence. Will probably do a log transform.
                    if UserInput.responses['response_data_types'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        pass
                if UserInput.responses['response_types'][responseIndex] == 'I' or UserInput.responses['response_types'][responseIndex] == 'P' or UserInput.responses['response_types'][responseIndex] == 'R': #For abscissa of concentration dependence.
                    if UserInput.responses['response_data_types'][responseIndex] == 'c':
                        pass
                    if UserInput.responses['response_data_types'][responseIndex] == 'r':
                        pass
        return nestedAllResponsesArray_transformed, nestedAllResponsesUncertainties_transformed  

    #Throughout this file, this function is called to generate initialStartPoint distributions or walkerInitialDistributions. 
    #    These are **not** identical variables and should not be messed up during editing. walkerInitialDistributions are a special case of initialStartPoints
    #    The distinction is hierarchical. Somebody could do a multiStart search with a uniformInitialDistributionType and have a walkerInitialDistribution started around each case within that.
    #    Effectively, the initialDistributionPoints are across parameter space, while the walkerInitialDistribution **could** be designed to find local modes using a smaller spread.
    
    def generateInitialPoints(self, numStartPoints=0, initialPointsDistributionType='uniform', relativeInitialDistributionSpread=1.0, numParameters = 0, centerPoint=None, gridsearchSamplingInterval = [], gridsearchSamplingRadii = []):
        #The initial points will be generated from a distribution based on the number of walkers and the distributions of the parameters.
        #The variable UserInput.std_prior has been populated with 1 sigma values, even for cases with uniform distributions.
        #The random generation at the front of the below expression is from the zeus example https://zeus-mcmc.readthedocs.io/en/latest/
        #The multiplication is based on the randn function using a sigma of one (which we then scale up) and then advising to add mu after: https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.random.randn.html
        #The actual numParameters cannot be 0. We just use 0 to mean not provided, in which case we pull it from the initial guess.
        #The arguments gridsearchSamplingInterval and gridsearchSamplingRadii are only for the distribution type 'grid', and correspond to the variables  gridsearchSamplingInterval = [], gridsearchSamplingRadii = [] inside getGridPermutations.
        if str(centerPoint).lower() == str(None).lower():
            centerPoint = np.array(self.UserInput.InputParameterInitialGuess)*1.0 #This may be a reduced parameter space.
        if initialPointsDistributionType.lower() not in ['grid', 'uniform', 'identical', 'gaussian', 'astroidal', 'sobol', 'shell']:
            print("Warning: initialPointsDistributionType must be from: 'grid', 'uniform', 'identical', 'gaussian', 'astroidal', 'sobol', 'shell'.  A different choice was received and is not understood.  initialPointsDistributionType is being set as 'sobol'.")
            initialPointsDistributionType = 'sobol'
        #For a multi-start with a grid, our algorithm is completely different than other cases.
        if initialPointsDistributionType.lower() =='grid':
            gridPermutations, numPermutations = self.getGridPermutations(centerPoint, gridsearchSamplingInterval=gridsearchSamplingInterval, gridsearchSamplingRadii=gridsearchSamplingRadii)
            initialPoints = gridPermutations
        #Below lines are for non-grid cases.
        if numParameters == 0:
            numParameters = len(centerPoint)
        if numStartPoints == 0: #This is a deprecated line. The function was originally designed for making mcmc walkers and then was generalized.
            numStartPoints = self.mcmc_nwalkers
        if initialPointsDistributionType.lower() =='uniform':
            initialPointsFirstTerm = np.random.uniform(-2,2, [numStartPoints,numParameters]) #<-- this is from me, trying to remove bias. This way we get sampling from a uniform distribution from -2 standard deviations to +2 standard deviations. That way the sampling is over 95% of the prior and is (according to the prior) likely to include the HPD region.
        elif initialPointsDistributionType.lower()  == 'identical':
            initialPointsFirstTerm = np.zeros((numStartPoints, numParameters)) #Make the first term all zeros.
        elif initialPointsDistributionType.lower() =='gaussian':
            initialPointsFirstTerm = np.random.randn(numStartPoints, numParameters) #<--- this was from the zeus example. TODO: change this to rng.standard_normal
        elif (initialPointsDistributionType.lower() == 'astroidal') or (initialPointsDistributionType.lower() == 'shell'):
            # The idea is to create a hypercube around the origin then apply a power law factor.
            # This factor is set as the numParameters to create an interesting distribution for Euclidean distance that starts as a uniform distribution then decays by a power law if the exponent is the number of dimensions. 
            from scipy.stats import qmc
            from warnings import catch_warnings, simplefilter #used to suppress warnings when sobol samples are not base2.
            # A sobol object has to be created to then extract points from the object.
            # The scramble (Owen Scramble) is always True. This option helps convergence and creates a more unbiased sampling.
            sobol_object = qmc.Sobol(d=numParameters, scramble=True)
            with catch_warnings():
                simplefilter("ignore")
                sobol_samples = sobol_object.random(numStartPoints)
            # now we must translate the sequence (from range(0,1) to range(-2,2)). This is analagous to the way we get sampling from a uniform distribution from -2 standard deviations to +2 standard deviations.
            initialPointsFirstTerm = -1 + 2*sobol_samples
            # This section assures that positive and negative values are generated.
            # create mapping scheme of negative values, then make matrix completely positive, apply negatives back later
            neg_map = np.ones((numStartPoints,numParameters), dtype=int)
            neg_map[initialPointsFirstTerm < 0] = -1
            initialPointsFirstTerm = np.abs(initialPointsFirstTerm)
            if initialPointsDistributionType.lower() == 'astroidal':
                initialPointsFirstTerm = initialPointsFirstTerm**numParameters
            elif initialPointsDistributionType.lower() == 'shell':
                initialPointsFirstTerm = initialPointsFirstTerm**(1/numParameters)
            initialPointsFirstTerm = neg_map*initialPointsFirstTerm
            # Apply a proportional factor of 2 to get bounds of 2 sigma. This is analagous to the way we get sampling from a uniform distribution from -2 standard deviations to +2 standard deviations.
            initialPointsFirstTerm *= 2
        elif initialPointsDistributionType.lower() == 'sobol':
            from scipy.stats import qmc
            from warnings import catch_warnings, simplefilter #used to suppress warnings when sobol samples are not base2.
            # A sobol object has to be created to then extract points from the object.
            # The scramble (Owen Scramble) is always True. This option helps convergence and creates a more unbiased sampling.
            sobol_object = qmc.Sobol(d=numParameters, scramble=True)
            with catch_warnings():
                simplefilter("ignore")
                sobol_samples = sobol_object.random(numStartPoints)
            # now we must translate the sequence (from range(0,1) to range(-2,2)). This is analagous to the way we get sampling from a uniform distribution from -2 standard deviations to +2 standard deviations.
            initialPointsFirstTerm = -2 + 4*sobol_samples
        if initialPointsDistributionType !='grid':
            #Now we add to centerPoint, usually self.UserInput.InputParameterInitialGuess. We don't use the UserInput initial guess directly because gridsearch and other things can change it -- so we need to use this one.
            initialPoints = relativeInitialDistributionSpread*initialPointsFirstTerm*self.UserInput.std_prior + centerPoint
        return initialPoints

    #This helper function has been made so that gridSearch and design of experiments can call it.
    #Although at first glance it may seem like it should be in the CombinationsGeneratorModule, that is a misconception. This is just a wrapper setting defaults for calling that module, such as using the prior for the grid interval when none is provided.
    #note that a blank list is okay for gridsearchSamplingInterval if doing a parameter grid, but not for other types of grids.
    def getGridPermutations(self, gridCenterVector, gridsearchSamplingInterval, gridsearchSamplingRadii, SpreadType="Addition",toFile=False):
        import PEUQSE.CombinationGeneratorModule as CombinationGeneratorModule
        numParameters = len(gridCenterVector)
        if len(gridsearchSamplingRadii) == 0:
            gridsearchSamplingRadii = np.ones(numParameters, dtype='int') #By default, will make ones.
            numPermutations = 3**numParameters
        else: 
            gridsearchSamplingRadii = np.array(gridsearchSamplingRadii, dtype='int')
            numPermutations = 1 #just initializing.
            for radius in gridsearchSamplingRadii:
                numPermutations=numPermutations*(2*radius+1)
        if len(gridsearchSamplingInterval) == 0:
            gridsearchSamplingInterval = self.UserInput.std_prior #By default, we use the standard deviations associated with the priors.
        else: gridsearchSamplingInterval = np.array(gridsearchSamplingInterval, dtype='float')
        gridPermutations = CombinationGeneratorModule.combinationGenerator(gridCenterVector, gridsearchSamplingInterval, gridsearchSamplingRadii, SpreadType=SpreadType,toFile=toFile)
        return gridPermutations, numPermutations  
        
    def doListOfPermutationsSearch(self, listOfPermutations, numPermutations = None, searchType='getLogP', exportLog = True, walkerInitialDistribution='UserChoice', passThroughArgs = {}, calculatePostBurnInStatistics=True,  keep_cumulative_post_burn_in_data = False, centerPoint=None, permutationsToSamples=False): #This is the 'engine' used by doGridSearch and  doMultiStartSearch
    #The listOfPermutations can also be another type of iterable.
        #Possible searchTypes are: 'getLogP', 'doEnsembleSliceSampling', 'doEnsembleJumpSampling', 'doMetropolisHastings', 'doOptimizeNegLogP', 'doOptimizeLogP' 'doOptimizeSSR'
        #permutationsToSamples should normally only be True if somebody is using gridsearch or uniform multistart with getLogP.
        self.listOfPermutations = listOfPermutations #This is being made into a class variable so that it can be used during parallelization
        if str(numPermutations).lower() == str(None).lower():
            numPermutations = len(self.listOfPermutations)
        if str(centerPoint).lower() == str(None).lower():
            centerPoint = self.UserInput.InputParameterInitialGuess*1.0
        if searchType == 'doGetLogP' or searchType == 'doSinglePoint': #Fixing a common input mistake.
            searchType = 'getLogP'
        if searchType == 'doEnsembleJumpSampling': #Fix the frontend name to the backend name. 
            searchType = 'doEnsembleModifiedMHSampling'
        self.permutation_searchType = searchType #This is mainly for consolidate_parallel_sampling_data
        verbose = self.UserInput.parameter_estimation_settings['verbose']
        if verbose:
            print("Starting multistart/permutations search.")
        if searchType in ['doEnsembleSliceSampling', 'doEnsembleJumpSampling', 'doMetropolisHastings']:
            if self.UserInput.user_requested_convergence_diagnostics == True:
                print("Notification: doMultistart / doListOfPermutationsSearch is being used with an MCMC searchType, which performs many MCMC runs. Convergence diagnostics are being turned off for the individual runs, and will be used only for the best run.")
                self.UserInput.user_requested_convergence_diagnostics = False
        

        file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix() #As of Nov 21st 2020, these should always be '' since multiStart_continueSampling is not intended to be used with parallel sampling.
        if (self.UserInput.parameter_estimation_settings['mcmc_continueSampling']  == 'auto') or (self.UserInput.parameter_estimation_settings['mcmc_continueSampling']  == False):
            mcmc_continueSampling = False #need to set this variable to false if it's an auto. The only time mcmc_continue sampling should be on for multistart is if someone is doing it intentionally, which would normally only be during an MPI case.
        if (self.UserInput.parameter_estimation_settings['mcmc_continueSampling']  == True):
            mcmc_continueSampling = True #The only time mcmc_continue sampling should be on for multistart is if someone is doing it intentionally, which would normally only be during an MPI case.
        #Check if we need to do multistart_continueSampling, and prepare for it if we need to.
        if ('multistart_continueSampling' not in self.UserInput.parameter_estimation_settings) or (self.UserInput.parameter_estimation_settings['multistart_continueSampling']  == 'auto'):
            if hasattr(self, 'multistart_MAP_logP_and_parameters_values'):
                multistart_continueSampling = True
            else:
                multistart_continueSampling = False
        else: multistart_continueSampling = self.UserInput.parameter_estimation_settings['multistart_continueSampling']
        if multistart_continueSampling == True:
            if hasattr(self, 'multistart_MAP_logP_and_parameters_values'): #if we are continuing from old results in the same instance
                self.last_permutations_MAP_logP_and_parameters_values = copy.deepcopy(self.permutations_MAP_logP_and_parameters_values)
            else: #Else we need to read from the file.
                self.last_permutations_MAP_logP_and_parameters_values_filename = file_name_prefix + "multistart_MAP_logP_and_parameters_values" + file_name_suffix
                self.last_permutations_MAP_logP_and_parameters_values = unpickleAnObject(self.UserInput.directories['pickles']+self.last_permutations_MAP_logP_and_parameters_values_filename)
            #extract he last_listOfPermutations from the array object.
            self.last_listOfPermutations =   np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_permutations_MAP_logP_and_parameters_values[:,1:])) #later columns are the permutations.
            if np.shape(self.last_listOfPermutations)[0] == 1: #In this case, need to transpose.
                self.last_listOfPermutations = self.last_listOfPermutations.transpose()
            #unlike in mcmc_continueSampling, we don't need the last_InputParameterInitialGuess information.
        #Initialize some things before permutations loop.
        allPermutationsResults = []
        self.permutations_MAP_logP_and_parameters_values = [] #Just initializing as fresh.
        if (type(self.UserInput.parameter_estimation_settings['multistart_checkPointFrequency']) != type(None)) or (verbose == True):
                timeAtPermutationSearchStart = time.time()
                timeAtLastPermutation = timeAtPermutationSearchStart #just initializing
        self.highest_logP = float('-inf') #just initializing
        highest_logP_parameter_set = np.ones(len(self.UserInput.InputParameterInitialGuess))*float('nan') #just initializing
        #bestResultSoFar has this form: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]   
        bestResultSoFar = [highest_logP_parameter_set, None, None, None, None, None, None] #just initializing. 
        highest_MAP_initial_point_index = None #just initializing
        highest_MAP_initial_point_parameters = None #just initializing
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
            self.permutations_unfiltered_map_simulated_outputs = []
        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doEnsembleModifiedMHSampling'):
            if str(self.UserInput.parameter_estimation_settings['mcmc_nwalkers']).lower() == 'auto':
                permutationSearch_mcmc_nwalkers = 2*len(centerPoint) #Lowest possible is 2 times num parameters for ESS.
            else:
                permutationSearch_mcmc_nwalkers = int(self.UserInput.parameter_estimation_settings['mcmc_nwalkers'])
        #Start grid search loop.
        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): #Choose the walker distribution type.
                if walkerInitialDistribution == 'UserChoice': #UserChoice comes from UserInput. It can still be auto.
                    walkerInitialDistribution = self.UserInput.parameter_estimation_settings['mcmc_walkerInitialDistribution']
                #The identical distribution is used by default because otherwise the walkers may be spread out too far and it could defeat the purpose of a gridsearch.
                if walkerInitialDistribution.lower() == 'auto':
                    walkerInitialDistribution = 'sobol'
        for permutationIndex,permutation in enumerate(self.listOfPermutations):
            #####Begin PEUQSE Parallel Processing During Loop Block####
            if (self.UserInput.parameter_estimation_settings['multistart_parallel_sampling'])== True:
                #We will only execute the sampling the permutationIndex matches the processor rank.
                #Additionally, if the rank is 0 and the simulation got here, it will be assumed the person is running this just to find the number of Permutations, so that will be spit out and the simulation ended.
                import PEUQSE.parallel_processing
                if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                    print("For the user input settings provided, the number of Permutations+1 will be",  numPermutations+1, ". Please use mpiexec or mpirun with this number for N. If you are not expecting to see this message, change your UserInput choices. You have chosen parallel processing for gridsearch and have run PEUQSE without mpi, which is a procedure to retrieve the number of processor ranks to use for parallelized gridsearch. A typical syntax now would be: mpiexec -n ",  numPermutations+1, " python runfile_for_your_analysis.py" )
                    sys.exit()
                elif PEUQSE.parallel_processing.currentProcessorNumber != permutationIndex+1:
                    continue #This means the permutation index does not match the processor rank so nothing should be executed.
                #elif PEUQSE.parallel_processing.currentProcessorNumber == permutationIndex+1:
                #    pass  #This is the "normal" case and is implied, so is commented out.
            #####End PEUQSE Parallel Processing During Loop Block####
            self.UserInput.InputParameterInitialGuess = permutation #We need to fill the variable InputParameterInitialGuess with the permutation being checked.
            if (searchType == 'getLogP'):
                self.map_logP = self.getLogP(permutation) #The getLogP function does not fill map_logP by itself.
                self.map_parameter_set = permutation
                thisResult = [self.map_parameter_set, None, None, None, None, None, None]
                #thisResultStr = [self.map_logP, str(self.map_parameter_set).replace(",","|").replace("[","").replace('(','').replace(')',''), 'None', 'None', 'None', 'None', 'None', 'None']
            if searchType == 'doMetropolisHastings':
                self.map_logP = np.float('-inf') #initializing as -inf to have a 'pure' mcmc sampling.
                thisResult = self.doMetropolisHastings(calculatePostBurnInStatistics=calculatePostBurnInStatistics, continueSampling=mcmc_continueSampling)
                #self.map_logP gets done by itself in doMetropolisHastings
                #Note that "thisResult" has the form: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
                if keep_cumulative_post_burn_in_data == True:
                    if permutationIndex == 0:
                        self.cumulative_post_burn_in_samples = self.post_burn_in_samples
                        self.cumulative_discrete_chains_post_burn_in_samples = self.discrete_chains_post_burn_in_samples
                        self.cumulative_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                    else: #This is basically elseif permutationIndex > 0:
                        self.cumulative_post_burn_in_samples = np.vstack((self.cumulative_post_burn_in_samples, self.post_burn_in_samples))
                        self.cumulative_discrete_chains_post_burn_in_samples = np.vstack((self.cumulative_discrete_chains_post_burn_in_samples, self.discrete_chains_post_burn_in_samples)) 
                        self.cumulative_post_burn_in_log_priors_vec = np.vstack((self.cumulative_post_burn_in_log_priors_vec, self.post_burn_in_log_priors_vec))
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = np.vstack((self.cumulative_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))
            if searchType == 'doEnsembleSliceSampling':
                self.map_logP = np.float('-inf') #initializing as -inf to have a 'pure' mcmc sampling.
                thisResult = self.doEnsembleSliceSampling(mcmc_nwalkers_direct_input=permutationSearch_mcmc_nwalkers, calculatePostBurnInStatistics=calculatePostBurnInStatistics, walkerInitialDistribution=walkerInitialDistribution, continueSampling=mcmc_continueSampling) 
                #Note that "thisResult" has the form: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
                #self.map_logP gets done by itself in doEnsembleSliceSampling
                if keep_cumulative_post_burn_in_data == True:
                    if permutationIndex == 0:
                        self.cumulative_post_burn_in_samples = self.post_burn_in_samples
                        self.cumulative_discrete_chains_post_burn_in_samples = self.discrete_chains_post_burn_in_samples
                        self.cumulative_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                    else: #This is basically elseif permutationIndex > 0:
                        self.cumulative_post_burn_in_samples = np.vstack((self.cumulative_post_burn_in_samples, self.post_burn_in_samples))
                        self.cumulative_discrete_chains_post_burn_in_samples = np.vstack((self.cumulative_discrete_chains_post_burn_in_samples, self.discrete_chains_post_burn_in_samples)) 
                        self.cumulative_post_burn_in_log_priors_vec = np.vstack((self.cumulative_post_burn_in_log_priors_vec, self.post_burn_in_log_priors_vec))
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = np.vstack((self.cumulative_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))
            if searchType == 'doEnsembleModifiedMHSampling':
                self.map_logP = np.float('-inf') #initializing as -inf to have a 'pure' mcmc sampling.
                thisResult = self.doEnsembleModifiedMHSampling(mcmc_nwalkers_direct_input=permutationSearch_mcmc_nwalkers, calculatePostBurnInStatistics=calculatePostBurnInStatistics, walkerInitialDistribution=walkerInitialDistribution, continueSampling=mcmc_continueSampling) 
                #Note that "thisResult" has the form: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
                #self.map_logP gets done by itself in doEnsembleJumpSampling
                if keep_cumulative_post_burn_in_data == True:
                    if permutationIndex == 0:
                        self.cumulative_post_burn_in_samples = self.post_burn_in_samples
                        self.cumulative_discrete_chains_post_burn_in_samples = self.discrete_chains_post_burn_in_samples
                        self.cumulative_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                    else: #This is basically elseif permutationIndex > 0:
                        self.cumulative_post_burn_in_samples = np.vstack((self.cumulative_post_burn_in_samples, self.post_burn_in_samples))
                        self.cumulative_discrete_chains_post_burn_in_samples = np.vstack((self.cumulative_discrete_chains_post_burn_in_samples, self.discrete_chains_post_burn_in_samples)) 
                        self.cumulative_post_burn_in_log_priors_vec = np.vstack((self.cumulative_post_burn_in_log_priors_vec, self.post_burn_in_log_priors_vec))
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = np.vstack((self.cumulative_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))
            if searchType == 'doOptimizeLogP':
                optimizationOutput = self.doOptimizeLogP(**passThroughArgs)
                self.map_logP = optimizationOutput[1] 
                self.map_parameter_set = optimizationOutput[0]
                thisResult = [self.map_parameter_set, None, None, None, None, None, None]
            if searchType == 'doOptimizeNegLogP':
                optimizationOutput = self.doOptimizeNegLogP(**passThroughArgs)
                self.map_logP = -1.0*optimizationOutput[1] #need to times by negative 1 to convert negLogP into P.
                self.map_parameter_set = optimizationOutput[0]
                thisResult = [self.map_parameter_set, None, None, None, None, None, None]
            if searchType == 'doOptimizeSSR':
                optimizationOutput = self.doOptimizeSSR(**passThroughArgs)
                self.map_logP = -1.0*optimizationOutput[1]  #The SSR is a minimizing objective function, so we multiply by -1 to make it analagous to a log_P.
                self.map_parameter_set = optimizationOutput[0]
                thisResult = [self.map_parameter_set, None, None, None, None, None, None]
            if (type(self.UserInput.parameter_estimation_settings['multistart_checkPointFrequency']) != type(None)) or (verbose == True):
                timeAtThisPermutation = time.time()
                timeOfThisPermutation = timeAtThisPermutation - timeAtLastPermutation
                averageTimePerPermutation = (timeAtThisPermutation - timeAtPermutationSearchStart)/(permutationIndex+1)
                numRemainingPermutations = numPermutations - permutationIndex+1
                timeAtLastPermutation = timeAtThisPermutation #Updating.
            if self.map_logP > self.highest_logP: #This is the grid point in space with the highest value found so far and will be kept.
                bestResultSoFar = thisResult #for mcmc: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
                self.highest_logP = np.copy(self.map_logP)
                highest_logP_parameter_set = np.copy(self.map_parameter_set)
                highest_MAP_initial_point_index = permutationIndex
                highest_MAP_initial_point_parameters = permutation
                if hasattr(self, 'discrete_chains_post_burn_in_samples'):
                    bestResultSoFar_discrete_chains_post_burn_in_samples = self.discrete_chains_post_burn_in_samples
            allPermutationsResults.append(thisResult)
            if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
                if (searchType == 'doEnsembleSliceSampling') or (searchType=='doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): #we need to run the map again, outside of mcmc, to populate 
                    self.map_logP = self.getLogP(self.map_parameter_set) #this has an implied return of self.lastSimulatedResponses.
                #else no extra work needs to be done since the last simulation was the map.
                self.permutations_unfiltered_map_simulated_outputs.append(np.array(self.lastSimulatedResponses).flatten())
            self.permutations_MAP_logP_and_parameters_values.append(np.hstack((self.map_logP, self.map_parameter_set)))    
            if verbose == True:
                print("Permutation", permutation, "number", permutationIndex+1, "out of", numPermutations, "timeOfThisPermutation", timeOfThisPermutation)
                print("Permutation", permutationIndex+1, "averageTimePerPermutation", "%.2f" % round(averageTimePerPermutation,2), "estimated time remaining", "%.2f" % round( numRemainingPermutations*averageTimePerPermutation,2), "s" )
                print("Permutation", permutationIndex+1, "current logP", self.map_logP, "highest logP", self.highest_logP, "highest logP Parameter Set", highest_logP_parameter_set)
            elif type(self.UserInput.parameter_estimation_settings['multistart_checkPointFrequency']) != type(None): #If verbose off but checkpoint frequency is on.
                if (permutationIndex ==0 or ((permutationIndex+1)/self.UserInput.parameter_estimation_settings['multistart_checkPointFrequency']).is_integer()):
                    print("Permutation", permutation, "number", permutationIndex+1, "out of", numPermutations, "timeOfThisPermutation", timeOfThisPermutation)
                    print("Permutation", permutationIndex+1, "averageTimePerPermutation", "%.2f" % round(averageTimePerPermutation,2), "estimated time remaining", "%.2f" % round( numRemainingPermutations*averageTimePerPermutation,2), "s" )
                    print("Permutation", permutationIndex+1, "current logP", self.map_logP, "highest logP", self.highest_logP)
        ####START BLOCK RELATED TO PARALLEL SAMPLING####
        if (self.UserInput.parameter_estimation_settings['multistart_parallel_sampling']) == True: #This is the parallel sampling mpi case. #Consider later adding self.UserInput.parameter_estimation_settings['permutation_parallel_sampling'])
            #We are going to export all of the relevant statistics for each permutation.
            self.exportPostPermutationStatistics(searchType = searchType) #this is needed for **each** permutation if parallel sampling is being done.
            self.checkIfAllParallelSimulationsDone("permutation"+"_map_logP_") #This checks if we are on the final process and also sets the global variable for it accordingly.
            if PEUQSE.parallel_processing.finalProcess == False:
                return self.map_logP #This is sortof like a sys.exit(), we are just ending the PermutationSearch function here if we are not on the finalProcess. 
            if PEUQSE.parallel_processing.finalProcess == True:
                self.UserInput.parameter_estimation_settings['multistart_parallel_sampling'] = False ##We are turning off the parallel sampling variable because the parallel sampling is over now. The export log will become export extra things if we keep this on for the next step.
                self.consolidate_parallel_sampling_data(parallelizationType="permutation", mpi_cached_files_prefix='permutation') #this parallelizationType means "keep only the best, don't average"
                
        ####END BLOCK RELATED TO PARALLEL SAMPLING####
        ####Doing some statistics across the full permutation set.  TODO: Consider merging this into exportPostPermutationStatistics and calling that same function again, which is what I think the mcmc parallel sampling does. But the filenames are different, so some care would be needed if that is going to be done.####
        #TODO: export the allPermutationsResults to file at end of search in a nicer format.        
        #set the initial guess back to the center of the grid.
        self.UserInput.InputParameterInitialGuess = centerPoint
        #populate the map etc. with those of the best result.
        self.map_logP = self.highest_logP 
        self.map_parameter_set = highest_logP_parameter_set 
        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
            #For MCMC, we can now calculate the post_burn_in statistics for the best sampling from the full samplings done. We don't want to lump all together because that would not be unbiased.
            #Note that "thisResult" and thus "bestResultSoFar" has the form: [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
            self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec = bestResultSoFar
            if calculatePostBurnInStatistics == True:
                #self.post_burn_in_samples = bestResultSoFar[5] #Setting the global variable will allow calculating the info gain and priors also.
                #self.post_burn_in_log_posteriors_un_normed_vec = bestResultSoFar[6]
                self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True)
                self.exportPostBurnInStatistics()
            if self.UserInput.parameter_estimation_settings['convergence_diagnostics'] == True: #at the end of the permutations, calculate convergence diagnostics for the last run.
                if keep_cumulative_post_burn_in_data == True:
                    self.getConvergenceDiagnostics(bestResultSoFar_discrete_chains_post_burn_in_samples)
                else:
                    print("Warning: could not get convergence diagnostics for multi-start because keep_cumulative_post_burn_in_data was set to False.")
            #One could call calculatePostBurnInStatistics() if one wanted the cumulative from all results. But we don't actually want that.
            #Below should not be used. These commented out lines are biased towards the center of the grid.
            #self.post_burn_in_samples = cumulative_post_burn_in_samples
            #self.post_burn_in_log_priors_vec = cumulative_post_burn_in_log_priors_vec
            #self.post_burn_in_log_posteriors_un_normed_vec = cumulative_post_burn_in_log_posteriors_un_normed_vec
            #implied return bestResultSoFar # [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec] 
        if (searchType == 'getLogP') or (searchType == 'doOptimizeNegLogP') or (searchType == 'doOptimizeSSR') or (searchType == 'doOptimizeLogP'):
            if (searchType == 'doOptimizeNegLogP'):
                self.permutation_and_doOptimizeNegLogP = True #turning on this flag for case of permutation_and_doOptimize. This is needed so that a warning can be put in the mcmc_log.
            if (searchType == 'doOptimizeSSR'): 
                self.permutation_and_doOptimizeSSR = True
            if (searchType == 'doOptimizeLogP'):
                self.permutation_and_doOptimizeLogP = True                
            #if it's getLogP gridsearch, we are going to convert it to samples if requested.
            if permutationsToSamples == True:
                self.permutations_MAP_logP_and_parameters_values = np.vstack( self.permutations_MAP_logP_and_parameters_values) #Note that vstack actually requires a tuple with multiple elements as an argument. So this list or array like structure is being converted to a tuple of many elements and then being stacked.
                #now stack with earlier results for multistart_continueSampling if needed.
                if multistart_continueSampling == True:
                        self.permutations_MAP_logP_and_parameters_values = np.vstack((self.last_permutations_MAP_logP_and_parameters_values,self.permutations_MAP_logP_and_parameters_values))                        
                        self.listOfPermutations = np.vstack((self.last_listOfPermutations, self.listOfPermutations))
                        highest_MAP_initial_point_index = "Not provided with continueSampling." #TODO: take self.map_parameter_set from after calculatePostBurnIn Statistics highest_MAP_initial_point_index and search for the right row in listOfPermutations.
                #First set the multistart_permutationsToSamples_threshold_filter_coefficient. We will take 10**-(thisnumber) later.
                if str(self.UserInput.parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient']).lower() == 'auto':
                    multistart_permutationsToSamples_threshold_filter_coefficient = 2.0
                else:
                    multistart_permutationsToSamples_threshold_filter_coefficient = self.UserInput.parameter_estimation_settings['multistart_permutationsToSamples_threshold_filter_coefficient']
                try:
                    logP_values_and_samples = convertPermutationsToSamples(self.permutations_MAP_logP_and_parameters_values, maxLogP=float(self.map_logP), relativeFilteringThreshold = 10**(-1*multistart_permutationsToSamples_threshold_filter_coefficient))
                    self.post_burn_in_log_posteriors_un_normed_vec = logP_values_and_samples[:,0]
                    self.post_burn_in_log_posteriors_un_normed_vec = np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.post_burn_in_log_posteriors_un_normed_vec)).transpose()
                    self.post_burn_in_samples = logP_values_and_samples[:,1:]
                    #need to populate post_burn_in_log_priors_vec this with an object, otherwise calculatePostBurnInStatistics will try to calculate all the priors.
                    self.post_burn_in_log_priors_vec = None
                    #Below is needed to avoid causing an error in the calculatePostBurnInStatistics since we don't have a real priors vec.
                    self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_samples'] = False
                    self.calculatePostBurnInStatistics()
                    # create discrete_chains_post_burn_in_samples
                    self.discrete_chains_post_burn_in_samples = np.expand_dims(self.post_burn_in_samples, axis=1)
                    if self.UserInput.parameter_estimation_settings['convergence_diagnostics']: #Run convergence diagnostics if UserInput defines it as True
                        self.getConvergenceDiagnostics(self.discrete_chains_post_burn_in_samples)
                except:
                    print("Could not convertPermutationsToSamples. This usually means there were no finite probability points sampled.")
                    permutationsToSamples = False #changing to false to prevent errors during exporting.
                
            #implied return bestResultSoFar# [self.map_parameter_set, self.map_logP]
        #This has to be below the later parts so that permutationsToSamples can occur first.
        if exportLog == True:
            pass #Later will do something with allPermutationsResults variable. It has one element for each result (that is, each permutation).
        with open(self.UserInput.directories['logs_and_csvs'] + "multistart_log_file.txt", 'w') as out_file:
                out_file.write("centerPoint: " + str(centerPoint) + "\n")
                if self.permutation_and_doOptimizeSSR == False:# In the normal case, we are not doing SSR.               
                    out_file.write("highest_MAP_logP: " + str(self.map_logP) + "\n")
                    out_file.write("highest_MAP_logP_parameter_set: " + str(bestResultSoFar[0])+ "\n")
                    out_file.write("highest_MAP_initial_point_index: " + str(highest_MAP_initial_point_index)+ "\n")
                    out_file.write("highest_MAP_initial_point_parameters: " + str( highest_MAP_initial_point_parameters)+ "\n")
                    if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling') or (permutationsToSamples == True):
                        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): 
                            caveat = ' (for the above initial point) '
                        elif permutationsToSamples == True:
                            caveat = ''
                        out_file.write("self.mu_AP_parameter_set : " + caveat + str( bestResultSoFar[1])+ "\n")
                        out_file.write("self.stdap_parameter_set : " + caveat  + str( bestResultSoFar[2])+ "\n")
                    if (self.permutation_and_doOptimizeNegLogP == True) or (self.permutation_and_doOptimizeLogP == True):
                        out_file.write("\n WARNING: It appears this run used a doOptimize with multi-start. In this case, the MAP_logP and map_parameter_set are the optimum from the run.  However, the mu_AP_parameter_set and stdap_parameter_set are not meaningful, since this was not an even weighted exploration of the posterior. \n")                        
                if self.permutation_and_doOptimizeSSR == True: #special case where we are doing SSR.
                    out_file.write("Below, negSSR means the SSR times -1. This is the optimum from the run. \n")
                    out_file.write("highest_negSSR: " + str(self.map_logP) + "\n")
                    out_file.write("highest_negSSR_parameter_set: " + str(self.map_parameter_set)+ "\n")
                    out_file.write("highest_negSSR_initial_point_index: " + str(highest_MAP_initial_point_index)+ "\n")
                    out_file.write("highest_negSSR_initial_point_parameters: " + str( highest_MAP_initial_point_parameters)+ "\n")
                    caveat = ' (actually just an analogue) '
                    out_file.write("self.mu_AP_parameter_set : " + caveat + str( self.mu_AP_parameter_set)+ "\n")
                    out_file.write("self.stdap_parameter_set : " + caveat  + str( self.stdap_parameter_set)+ "\n")
                    if self.permutation_and_doOptimizeSSR == True:
                        out_file.write("\n WARNING: It appears this run used a doOptimize with multi-start. In this case, the highest_negSSR and highest_negSSR_parameter_set are at the optimum from the run.  However, the mu_AP_parameter_set and stdap_parameter_set are not meaningful, since this was not an even weighted exploration of the posterior. However, the values have been reported that are analogous if the negSSR is treated as a logP, and the posterior graphs have been made accordingly. \n")                        

                    
        #do some exporting etc. This is at the end to avoid exporting every single time if parallelization is used.
        np.savetxt(self.UserInput.directories['logs_and_csvs']+'multistart_initial_points_parameters_values'+'.csv', self.listOfPermutations, delimiter=",")
        np.savetxt(self.UserInput.directories['logs_and_csvs']+'multistart_MAP_logP_and_parameters_values.csv',self.permutations_MAP_logP_and_parameters_values, delimiter=",")
        pickleAnObject(self.permutations_MAP_logP_and_parameters_values, self.UserInput.directories['pickles']+file_name_prefix+'multistart_MAP_logP_and_parameters_values'+file_name_suffix)
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
            np.savetxt(self.UserInput.directories['logs_and_csvs']+'multistart_unfiltered_map_simulated_outputs'+'.csv', self.permutations_unfiltered_map_simulated_outputs, delimiter=",")       
        print("Final map parameter results from PermutationSearch:", self.map_parameter_set,  " \nFinal map logP:", self.map_logP, "more details available in multistart_log_file.txt")        
        return bestResultSoFar# [self.map_parameter_set, self.map_logP, etc.]

    #@CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def doMultiStart(self, searchType='UserChoice', numStartPoints = 'UserChoice', relativeInitialDistributionSpread='UserChoice', exportLog = 'UserChoice', initialPointsDistributionType='UserChoice', passThroughArgs = 'UserChoice', calculatePostBurnInStatistics='UserChoice',  keep_cumulative_post_burn_in_data = 'UserChoice', walkerInitialDistribution='UserChoice', centerPoint = None, gridsearchSamplingInterval = 'UserChoice', gridsearchSamplingRadii = 'UserChoice'):
        #See doListOfPermutationsSearch for possible values of searchType variable
        #This function is basically a wrapper that creates a list of initial points and then runs a 'check each permutation' search on that list.
        #We set many of the arguments to have blank or zero values so that if they are not provided, the values will be taken from the UserInput choices.
        if str(searchType) == 'UserChoice': 
            searchType = self.UserInput.parameter_estimation_settings['multistart_searchType']
        if str(initialPointsDistributionType) == 'UserChoice': 
            initialPointsDistributionType = self.UserInput.parameter_estimation_settings['multistart_initialPointsDistributionType']
        if str(numStartPoints) =='UserChoice':
            numStartPoints = self.UserInput.parameter_estimation_settings['multistart_numStartPoints']
        if str(relativeInitialDistributionSpread) == 'UserChoice': 
            relativeInitialDistributionSpread = self.UserInput.parameter_estimation_settings['multistart_relativeInitialDistributionSpread']
        if str(gridsearchSamplingInterval) == 'UserChoice':
            gridsearchSamplingInterval = self.UserInput.parameter_estimation_settings['multistart_gridsearchSamplingInterval']
        if str(gridsearchSamplingRadii) == 'UserChoice':
            gridsearchSamplingRadii = self.UserInput.parameter_estimation_settings['multistart_gridsearchSamplingRadii']
        if str(exportLog) == 'UserChoice':
            exportLog = self.UserInput.parameter_estimation_settings['multistart_exportLog']
        if str(passThroughArgs) == 'UserChoice':
            passThroughArgs = self.UserInput.parameter_estimation_settings['multistart_passThroughArgs']
        if str(keep_cumulative_post_burn_in_data) == 'UserChoice':
            keep_cumulative_post_burn_in_data = self.UserInput.parameter_estimation_settings['multistart_keep_cumulative_post_burn_in_data']
        if str(calculatePostBurnInStatistics) == 'UserChoice':
            calculatePostBurnInStatistics = self.UserInput.parameter_estimation_settings['multistart_calculatePostBurnInStatistics']
        if numStartPoints == 0: #if it's still zero, we need to make it the default which is 3 times the number of active parameters.
            numStartPoints = len(self.UserInput.InputParameterInitialGuess)*3
        if relativeInitialDistributionSpread == 0: #if it's still zero, we need to make it the default which is 0.866
            relativeInitialDistributionSpread = 0.866 #This choice is to be helpful for uniform distribution cases, as described in the user input.
        if searchType == 'doGetLogP' or searchType == 'doSinglePoint': #Fixing a common input mistake.
            searchType = 'getLogP'
        #make the initial points list by mostly passing through arguments.
        multiStartInitialPointsList = self.generateInitialPoints(numStartPoints=numStartPoints, relativeInitialDistributionSpread=relativeInitialDistributionSpread, initialPointsDistributionType=initialPointsDistributionType, centerPoint = centerPoint, gridsearchSamplingInterval = gridsearchSamplingInterval, gridsearchSamplingRadii = gridsearchSamplingRadii)
        #we normally only turn on permutationsToSamples if grid or uniform and if getLogP or doOptimizeNegLogP.
        permutationsToSamples = False#initialize with default
        if self.UserInput.parameter_estimation_settings['multistart_permutationsToSamples'] == True:
            if (initialPointsDistributionType == 'grid') or (initialPointsDistributionType == 'uniform') or (initialPointsDistributionType == 'sobol') or (initialPointsDistributionType == 'astroidal') or (initialPointsDistributionType == 'shell'):
                if (searchType == 'getLogP') or (searchType=='doOptimizeNegLogP') or (searchType=='doOptimizeLogP') or (searchType=='doOptimizeSSR'):
                    permutationsToSamples = True
        #Look for the best result (highest map_logP) from among these permutations. Maybe later should add optional argument to allow searching for highest mu_AP to find HPD.
        bestResultSoFar = self.doListOfPermutationsSearch(listOfPermutations=multiStartInitialPointsList, searchType=searchType, exportLog=exportLog, walkerInitialDistribution=walkerInitialDistribution, passThroughArgs=passThroughArgs, calculatePostBurnInStatistics=calculatePostBurnInStatistics, keep_cumulative_post_burn_in_data=keep_cumulative_post_burn_in_data, centerPoint = centerPoint, permutationsToSamples=permutationsToSamples)
        return bestResultSoFar
  
    #@CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def doGridSearch(self, searchType='getLogP', exportLog = True, gridSamplingAbsoluteIntervalSize = [], gridSamplingNumOfIntervals = [], passThroughArgs = {}, calculatePostBurnInStatistics=True,  keep_cumulative_post_burn_in_data = False, walkerInitialDistribution='UserChoice'):
        print("Warning: You have called doGridSearch.  This function is deprecated and is only retained for old examples. Please use doMultiStart with multistart_initialPointsDistributionType = 'grid' ")
        # gridSamplingNumOfIntervals is the number of variations to check in units of variance for each parameter. Can be 0 if you don't want to vary a particular parameter in the grid search.
        #calculatePostBurnInStatistics will store all the individual runs in memory and will then provide the samples of the best one.
        #TODO: the upper part of the gridsearch may not be compatibile with reduced parameter space. Needs to be checked.
        gridCenter = self.UserInput.InputParameterInitialGuess*1.0 #This may be a reduced parameter space.    
        gridPermutations, numPermutations = self.getGridPermutations(gridCenter, gridSamplingAbsoluteIntervalSize, gridSamplingNumOfIntervals)
        bestResultSoFar = self.doListOfPermutationsSearch(gridPermutations, numPermutations = numPermutations, searchType=searchType, exportLog = exportLog, walkerInitialDistribution=walkerInitialDistribution, passThroughArgs=passThroughArgs, calculatePostBurnInStatistics=calculatePostBurnInStatistics,  keep_cumulative_post_burn_in_data = keep_cumulative_post_burn_in_data, centerPoint = gridCenter)
        return bestResultSoFar

    def checkIfAllParallelSimulationsDone(self, fileNameBase, fileNamePrefix='', fileNameSuffix=''):
        import PEUQSE.parallel_processing
        #PEUQSE.parallel_processing.currentProcessorNumber
        numSimulations = PEUQSE.parallel_processing.numSimulations
        import os
        os.chdir(self.UserInput.directories['pickles']+"mpi_cached_files")
        #now make a list of what we expect.
        simulationsKey = np.ones(numSimulations)
        working_dir=os.getcwd()
        filesInDirectory=os.listdir(working_dir)
        for simulationIndex in range(0,numSimulations): #For each simulation, we check if it's there and set the simulation key to 0 if it is done.
            simulationNumberString = str(simulationIndex+1)
            for name in filesInDirectory:
                if fileNamePrefix+fileNameBase+simulationNumberString+fileNameSuffix+".pkl" in name:
                    simulationsKey[simulationIndex] = 0
                    filesInDirectory.remove(name) #Removing so it won't be checked for again, to speed up next search.
        if self.UserInput.directories['pickles'] == "./":
            os.chdir("..") #change directory back regardless.
        else:
            os.chdir("../..") #change directory back regardless.
        if np.sum(simulationsKey) == 0:
            PEUQSE.parallel_processing.finalProcess = True
            return True
        else: #if simulationsKey is not zero, then we return False b/c not yet finsihed.
            PEUQSE.parallel_processing.finalProcess = False
            return False

    def consolidate_parallel_doe_data(self, parallelizationType='conditions'):
        import PEUQSE.parallel_processing
        #PEUQSE.parallel_processing.currentProcessorNumber
        numSimulations = PEUQSE.parallel_processing.numSimulations
        parModulationNumber = int(self.parModulationPermutationIndex + 1)
        #We will check **only** for this parModulationNumber. That way, it this processor is the last to finish this parModulation, it will do the infoGainMatrix stacking.
        if self.checkIfAllParallelSimulationsDone("conditionsPermutationAndInfoGain_mod"+str(parModulationNumber)+"_cond") == True:
            if parallelizationType.lower() == 'conditions':
                import os
                os.chdir(self.UserInput.directories['pickles'] + "mpi_cached_files")
                self.info_gain_matrix = [] #Initializing this as a blank list, it will be made into an array after the loop.
                for simulationIndex in range(0,numSimulations): #For each simulation, we need to grab the results.
                    simulationNumberString = str(simulationIndex+1)
                    #Getting the data out.    
                    current_conditionsPermutationAndInfoGain_filename = "conditionsPermutationAndInfoGain_mod"+str(parModulationNumber)+"_cond"+simulationNumberString
                    current_conditionsPermutationAndInfoGain_data = unpickleAnObject(current_conditionsPermutationAndInfoGain_filename)
                    #accumulating.
                    self.info_gain_matrix.append(current_conditionsPermutationAndInfoGain_data)                        
                #Now we'll make this info_gain_matrix into an array and pickle it. It will be an implied return.
                self.info_gain_matrix = np.array(self.info_gain_matrix)
                current_parModulationInfoGainMatrix_filename = "parModulationInfoGainMatrix_mod"+str(parModulationNumber)
                pickleAnObject(self.info_gain_matrix,current_parModulationInfoGainMatrix_filename)
                #Change back to the regular directory since we are done.
                if self.UserInput.directories['pickles'] == "./":
                    os.chdir("..")
                else:
                    os.chdir("../..")
                return True #so we know we're done.
        else:
            return False #this means we weren't done.
            
    def consolidate_parallel_doe_info_gain_matrices(self):
        import PEUQSE.parallel_processing
        numSimulations = PEUQSE.parallel_processing.numSimulations        
        import os
        os.chdir(self.UserInput.directories['pickles'] + "mpi_cached_files")
        info_gains_matrices_list = [] #Initializing this as a blank list, it will be made into an array after the loop.
        for parModulationIndex in range(0,self.numParModulationPermutations): #For each simulation, we need to grab the results.
            parModulationNumberString = str(parModulationIndex+1)
            #Getting the data out.    
            current_parModulationInfoGainMatrix_filename = "parModulationInfoGainMatrix_mod"+parModulationNumberString 
            current_parModulationInfoGainMatrix_data = unpickleAnObject(current_parModulationInfoGainMatrix_filename)
            #accumulating.
            info_gains_matrices_list.append(current_parModulationInfoGainMatrix_data)                        
        #nothing more needs to be done except making it into an array: self.info_gains_matrices_array is an implied return.
        self.info_gains_matrices_array=np.array(info_gains_matrices_list)
        if self.UserInput.directories['pickles'] == "./":
            os.chdir("..")
        else:
            os.chdir("../..")

 
    def consolidate_parallel_sampling_data(self, parallelizationType='equal', mpi_cached_files_prefix=''):
        #parallelizationType='equal' means everything will get averaged together. parallelizationType='permutation' will be treated differently, keeps only the best.
        #mpi_cached_files_prefix can be 'mcmc' or 'permutation' or '' and looks for a prefix before 'map_logP_6.pkl' where '6' would be the processor rank.
        import PEUQSE.parallel_processing
        #PEUQSE.parallel_processing.currentProcessorNumber
        numSimulations = PEUQSE.parallel_processing.numSimulations
        if self.checkIfAllParallelSimulationsDone(mpi_cached_files_prefix+"_map_logP_") == True: #FIXME: Need to make parallelization work even for non-mcmc
            if parallelizationType.lower() == 'permutation':
                searchType = self.permutation_searchType
                import os #All of the below happens in the pickles directory
                os.chdir(self.UserInput.directories['pickles']+"mpi_cached_files")
                self.listOfPermutations = [] #just initializing.
                self.permutations_MAP_logP_and_parameters_values = [] #just initializing.
                for simulationIndex in range(0,numSimulations): #For each simulation, we need to grab the results.
                    simulationNumberString = str(simulationIndex+1)
                    #Get the data out.    
                    
                    current_post_map_logP_filename = "permutation_map_logP_"+simulationNumberString
                    current_post_map_logP_data = unpickleAnObject(current_post_map_logP_filename)
                    self.map_logP = current_post_map_logP_data

                    current_post_initial_parameters_filename = "permutation_initial_point_parameters_"+simulationNumberString
                    current_post_initial_parameters_data = unpickleAnObject(current_post_initial_parameters_filename)
                    self.UserInput.InputParameterInitialGuess = current_post_initial_parameters_data

                    current_post_map_parameter_set_filename = "permutation_map_parameter_set_"+simulationNumberString
                    current_post_map_parameter_set_data = unpickleAnObject(current_post_map_parameter_set_filename)
                    self.map_parameter_set = current_post_map_parameter_set_data

                    if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
                        current_post_burn_in_statistics_filename = "permutation_post_burn_in_statistics_"+simulationNumberString
                        current_post_burn_in_statistics_data = unpickleAnObject(current_post_burn_in_statistics_filename)
                        [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec] = current_post_burn_in_statistics_data

                    #Still accumulating.
                    self.permutations_MAP_logP_and_parameters_values.append(np.hstack((self.map_logP, self.map_parameter_set)))
                    self.listOfPermutations.append(current_post_initial_parameters_data)
                    if simulationIndex == 0: #This is the first data set.
                        self.highest_logP = self.map_logP
                        self.highest_logP_parameter_set = self.map_parameter_set
                        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
                            self.highest_logP_post_burn_in_samples = self.post_burn_in_samples
                            self.highest_logP_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                            self.highest_logP_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                    else: #This is basically elseif permutationIndex > 0:
                        if self.highest_logP < self.map_logP:
                            self.highest_logP = self.map_logP
                            self.highest_logP_parameter_set = self.map_parameter_set
                            if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
                                self.highest_logP_post_burn_in_samples = self.post_burn_in_samples
                                self.highest_logP_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                                self.highest_logP_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                #After the loop is done, we want to keep the accumulated values and then do the regular final calculations.
                if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): #FIXME: These logic needs to be checked to make sure it is correct.
                    self.map_logP = max(self.post_burn_in_log_posteriors_un_normed_vec)
                    self.map_index = list(self.post_burn_in_log_posteriors_un_normed_vec).index(self.map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
                    self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the                 
                self.map_logP = self.highest_logP 
                self.map_parameter_set = self.highest_logP_parameter_set
                if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
                    self.post_burn_in_samples = self.highest_logP_post_burn_in_samples 
                    self.post_burn_in_log_priors_vec = self.highest_logP_post_burn_in_log_priors_vec 
                    self.post_burn_in_log_posteriors_un_normed_vec = self.highest_logP_post_burn_in_log_posteriors_un_normed_vec 
                #Now go back to the earlier directory since the consolidation is done.
                if self.UserInput.directories['pickles'] == "./":
                    os.chdir("..")
                else:
                    os.chdir("../..")
                if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
                    self.UserInput.request_mpi = False # we need to turn this off, because otherwise it will interfere with our attempts to calculate the post_burn_in statistics.
                    self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True) #The argument is provided because otherwise there can be some bad priors if ESS was used.
                    self.exportPostBurnInStatistics()
                    self.UserInput.request_mpi = True #Set this back to true so that consolidating plots etc. doesn't get messed up.
            elif parallelizationType.lower() == 'equal':
                import os #All of the below happens in the pickles directory
                os.chdir(self.UserInput.directories['pickles']+"./mpi_cached_files")
                #These pointers are initialized before the below loop. Mostly in case mpi never actually happened since then after the loop these would be empty.
                self.cumulative_post_burn_in_samples = self.post_burn_in_samples
                self.cumulative_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                self.cumulative_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                for simulationIndex in range(0,numSimulations): #For each simulation, we need to grab the results.
                    simulationNumberString = str(simulationIndex+1)
                    #Get the dat aout.    
                    current_post_burn_in_statistics_filename = "mcmc_post_burn_in_statistics_"+simulationNumberString
                    current_post_burn_in_statistics_data = unpickleAnObject(current_post_burn_in_statistics_filename)
                    #Populate the class variables.
                    [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec] = current_post_burn_in_statistics_data
                    #Still accumulating.
                    if simulationIndex == 0: #This is the first data set.
                        self.cumulative_post_burn_in_samples = self.post_burn_in_samples
                        self.cumulative_post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec
                    else: #This is basically elseif permutationIndex > 0:
                        self.cumulative_post_burn_in_samples = np.vstack((self.cumulative_post_burn_in_samples, self.post_burn_in_samples))
                        self.cumulative_post_burn_in_log_priors_vec = np.vstack((self.cumulative_post_burn_in_log_priors_vec, self.post_burn_in_log_priors_vec))
                        self.cumulative_post_burn_in_log_posteriors_un_normed_vec = np.vstack((self.cumulative_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))
                #After the loop is done, we want to keep the accumulated values and then do the regular final calculations.
                self.post_burn_in_samples = self.cumulative_post_burn_in_samples
                self.post_burn_in_log_priors_vec = self.cumulative_post_burn_in_log_priors_vec
                self.post_burn_in_log_posteriors_un_normed_vec = self.cumulative_post_burn_in_log_posteriors_un_normed_vec
                self.UserInput.request_mpi = False # we need to turn this off, because otherwise it will interfere with our attempts to calculate the post_burn_in statistics.
                self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] = False # we need to turn this off, because otherwise it will interfere with our attempts to calculate the post_burn_in statistics.
                if hasattr(self, "during_burn_in_samples"): #need to remove this so it doesn't get exported for the parallel case, since otherwise will export most recent one which is misleading.
                    delattr(self, "during_burn_in_samples")
                if self.UserInput.directories['pickles'] == "./":
                    os.chdir("..")
                else:
                    os.chdir("../..")
                self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True) #The argument is provided because otherwise there can be some bad priors if ESS was used.
                self.exportPostBurnInStatistics()
                self.UserInput.request_mpi = True #Set this back to true so that consolidating plots etc. doesn't get messed up.





    #The below function is a helper function that is used during doeInfoGainMatrix. However, it can certainly be used for other purposes.
    def populateResponsesWithSyntheticData(self, parModulationPermutation):
        #For each parameter Modulation Combination we are going to obtain a matrix of info_gains that is based on a grid of the independent_variables.
        #First we need to make some synthetic data using parModulationPermutation for the discreteParameterVector
        discreteParameterVector = parModulationPermutation
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
        simulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(simulatedResponses)
        #need to check if there are any 'responses_simulation_uncertainties'. #TODO: This isn't really implemented yet.
        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #if it's a None type, we keep it as a None type
            responses_simulation_uncertainties = None
        else:  #Else we get it based on the the discreteParameterVector
            responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(discreteParameterVector)
        
        synthetic_data  = simulatedResponses
        synthetic_data_uncertainties = responses_simulation_uncertainties
        #We need to populate the "observed" responses in UserInput with the synthetic data.
        self.UserInput.responses['responses_observed'] = simulatedResponses
        self.UserInput.responses['responses_observed_uncertainties'] = responses_simulation_uncertainties
        #Now need to do something unusual: Need to call the __init__ function again so that the arrays get reshaped as needed etc.
        self.__init__(self.UserInput)
    
    #This function requires first populating the doe_settings dictionary in UserInput in order to know which conditions to explore.
    software_name = "PEUQSE Bayesian Design of Experiments"
    software_version = "1.0.2"
    software_unique_id = "https://doi.org/10.1002/cctc.202000976"
    software_kwargs = {"version": software_version, "author": ["Eric A. Walker", "Kishore Ravisankar", "Aditya Savara"], "doi": "https://doi.org/10.1002/cctc.202000976", "cite": "Eric Alan Walker, Kishore Ravisankar, Aditya Savara. PEUQSE Intro 2: Harnessing Uncertainties from Data Sets, Bayesian Design of Experiments in Chemical Kinetics. ChemCatChem. Accepted. doi:10.1002/cctc.202000976"} 
    #@CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    @CiteSoft.function_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def doeGetInfoGainMatrix(self, parameterPermutation, searchType='doMetropolisHastings'):#Note: There is an implied argument of info_gains_matrices_array_format being 'xyz' or 'meshgrid'
        #At present, we *must* provide a parameterPermutation because right now the only way to get an InfoGainMatrix is with synthetic data assuming a particular parameterPermutation as the "real" or "actual" parameterPermutation.
        doe_settings = self.UserInput.doe_settings
        self.middle_of_doe_flag = True  #This is a work around that is needed because right now the synthetic data creation has an __init__ call which is going to try to modify the independent variables back to their original values if we don't do this.
        self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] = False #As of Oct 2020, mcmc_continueSampling is not compatible with design of experiments (doe) feature.
        self.info_gain_matrix = [] #Right now, if using KL_divergence, each item in here is a single array. It is a sum across all parameters. 
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
                conditionsGridPermutations, numPermutations = self.getGridPermutations(doe_settings['independent_variable_grid_center'], doe_settings['independent_variable_grid_interval_size'], doe_settings['independent_variable_grid_num_intervals'])
            #Here is the loop across conditions.                
            for conditionsPermutationIndex,conditionsPermutation in enumerate(conditionsGridPermutations):    
                #####Begin PEUQSE Parallel Processing During Loop Block####
                if (self.UserInput.doe_settings['parallel_conditions_exploration'])== True:
                    #We will only execute the sampling the permutationIndex matches the processor rank.
                    #Additionally, if the rank is 0 and the simulation got here, it will be assumed the person is running this just to find the number of Permutations, so that will be spit out and the simulation ended.
                    import PEUQSE.parallel_processing
                    if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                        print("For the user input settings provided, the number of Permutations+1 will be",  numPermutations+1, ". Please use mpiexec or mpirun with this number for N. If you are not expecting to see this message, change your UserInput choices. You have chosen parallel processing for gridsearch and have run PEUQSE without mpi, which is a procedure to retrieve the number of processor ranks to use for parallelized gridsearch. A typical syntax now would be: mpiexec -n ",  numPermutations+1, " python runfile_for_your_analysis.py" )
                        sys.exit()
                    elif PEUQSE.parallel_processing.currentProcessorNumber != conditionsPermutationIndex+1:
                        continue #This means the permutation index does not match the processor rank so nothing should be executed.
                    #elif PEUQSE.parallel_processing.currentProcessorNumber == permutationIndex+1:
                    #    pass  #This is the "normal" case and is implied, so is commented out.
                #####End PEUQSE Parallel Processing During Loop Block####
                #It is absolutely critical that we *do not* use syntax like self.UserInput.responses['independent_variables_values'] = xxxx
                #Because that would move where the pointer is going to. We need to instead populate the individual values in the simulation module's namespace.
                #This population Must occur here. It has to be after the indpendent variables have changed, before synthetic data is made, and before the MCMC is performed.
                self.UserInput.model['populateIndependentVariablesFunction'](conditionsPermutation)
                self.populateResponsesWithSyntheticData(parameterPermutation)
                if searchType=='doMetropolisHastings':
                    [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doMetropolisHastings()
                if searchType=='doEnsembleSliceSampling':
                    [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doEnsembleSliceSampling()
                if searchType == 'doEnsembleModifiedMHSampling':
                    [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doEnsembleModifiedMHSampling()
                conditionsPermutation = np.array(conditionsPermutation) #we're going to make this an array before adding to the info_gain matrix.
                conditionsPermutationAndInfoGain = np.hstack((conditionsPermutation, info_gain))
                self.info_gain_matrix.append(conditionsPermutationAndInfoGain)
                if (self.UserInput.doe_settings['parallel_conditions_exploration'])== True:
                    self.exportSingleConditionInfoGainMatrix(self.parameterPermutationNumber, conditionsPermutationAndInfoGain, conditionsPermutationIndex)
                if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above lines for the sum.
                    for parameterIndex in range(0,numParameters):#looping across number of parameters...
                        conditionsPermutationAndInfoGain = np.hstack((conditionsPermutation, np.array(self.info_gain_each_parameter[parameterIndex]))) #Need to pull the info gain matrix from the nested objected named info_gain_each_parameter
                        #Below mimics the line above which reads self.info_gain_matrix.append(conditionsPermutationAndInfoGain)
                        info_gain_matrices_each_parameter[parameterIndex].append(conditionsPermutationAndInfoGain)
            self.info_gain_matrix = np.array(self.info_gain_matrix) #this is an implied return in addition to the real return.
            if self.UserInput.doe_settings['parallel_conditions_exploration'] == True: #We will overwrite self.info_gain_matrix with a consolidated one if needed.
                self.consolidate_parallel_doe_data(parallelizationType='conditions')
            if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above line for the sum.
                for parameterIndex in range(0,numParameters):#looping across number of parameters...
                    self.info_gain_matrices_each_parameter[parameterIndex]= np.array(info_gain_matrices_each_parameter[parameterIndex])
            self.middle_of_doe_flag = False #Set this back to false once info gain matrix is ready.
            return np.array(self.info_gain_matrix)            
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
                doSimulation = True #This is a temporary (short-lived) variable being made for parallel processing purposes.
                conditionsPermutationIndex = 0
                #We will not be using the function "self.getGridPermutations" for the loops because the meshgrid needs a different loop format.
                for indValue2 in independentVariable2ValuesArray: #We know from experience that the outer loop should be over the YY variable.
                    for indValue1 in independentVariable1ValuesArray: #We know from experience that the inner loop should be over the XX variable.
                        #It is absolutely critical that we *do not* use syntax like self.UserInput.responses['independent_variables_values'] = xxxx
                        #Because that would move where the pointer is going to. We need to instead populate the individual values in the simulation module's namespace.
                        #This population Must occur here. It has to be after the indpendent variables have changed, before synthetic data is made, and before the MCMC is performed.
                        #####Begin PEUQSE Parallel Processing During Loop Block -- This block is custom for meshgrid since the loop is different.####
                        if (self.UserInput.doe_settings['parallel_conditions_exploration'])== True:
                            numPermutations = len(independentVariable2ValuesArray)*len(independentVariable1ValuesArray)
                            permutationIndex = conditionsPermutationIndex
                            #We will only execute the sampling the permutationIndex matches the processor rank.
                            #Additionally, if the rank is 0 and the simulation got here, it will be assumed the person is running this just to find the number of Permutations, so that will be spit out and the simulation ended.
                            import PEUQSE.parallel_processing
                            if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                                print("For the user input settings provided, the number of Permutations+1 will be",  numPermutations+1, ". Please use mpiexec or mpirun with this number for N. If you are not expecting to see this message, change your UserInput choices. You have chosen parallel processing for gridsearch and have run PEUQSE without mpi, which is a procedure to retrieve the number of processor ranks to use for parallelized gridsearch. A typical syntax now would be: mpiexec -n ",  numPermutations+1, " python runfile_for_your_analysis.py" )
                                sys.exit()
                            elif PEUQSE.parallel_processing.currentProcessorNumber != conditionsPermutationIndex+1:
                                doSimulation = False #This means the permutation index does not match the processor rank so nothing should be executed.
                            elif PEUQSE.parallel_processing.currentProcessorNumber == permutationIndex+1:
                                doSimulation = True  #This is the "normal" case.
                        #####End PEUQSE Parallel Processing During Loop Block####
                        if doSimulation == True:
                            self.UserInput.model['populateIndependentVariablesFunction']([indValue1,indValue2])
                            self.populateResponsesWithSyntheticData(parameterPermutation)
                            if searchType=='doMetropolisHastings':
                                [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doMetropolisHastings()
                            if searchType=='doEnsembleSliceSampling':
                                [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doEnsembleSliceSampling()
                            if searchType=='doEnsembleModifiedMHSampling':
                                [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = self.doEnsembleModifiedMHSampling()
                            conditionsPermutation = np.array([indValue1,indValue2])
                            conditionsPermutationAndInfoGain = np.hstack((conditionsPermutation, info_gain))
                            self.info_gain_matrix.append(conditionsPermutationAndInfoGain) #NOTE that the structure *includes* the Permutations.
                            if (self.UserInput.doe_settings['parallel_conditions_exploration'])== True:
                                self.exportSingleConditionInfoGainMatrix(self.parModulationPermutationIndex+1, conditionsPermutationAndInfoGain, conditionsPermutationIndex)
                            if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above lines for the sum.
                                for parameterIndex in range(0,numParameters):#looping across number of parameters...
                                    conditionsPermutationAndInfoGain = np.hstack((conditionsPermutation, np.array(self.info_gain_each_parameter[parameterIndex]))) #Need to pull the info gain matrix from the nested objected named info_gain_each_parameter
                                    #Below mimics the line above which reads self.info_gain_matrix.append(conditionsPermutationAndInfoGain)
                                    info_gain_matrices_each_parameter[parameterIndex].append(conditionsPermutationAndInfoGain)
                        conditionsPermutationIndex = conditionsPermutationIndex + 1 #This variable was added for and is used in parallelization.
                self.info_gain_matrix = np.array(self.info_gain_matrix) #this is an implied return in addition to the real return.
                if self.UserInput.doe_settings['parallel_conditions_exploration'] == True: #We will overwrite self.info_gain_matrix with a consolidated one if needed.
                    self.consolidate_parallel_doe_data(parallelizationType='conditions')
                if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #copy the above line for the sum.
                    for parameterIndex in range(0,numParameters):#looping across number of parameters...
                        self.info_gain_matrices_each_parameter[parameterIndex]= np.array(info_gain_matrices_each_parameter[parameterIndex])
                self.middle_of_doe_flag = False #Set this back to false once info gain matrix is ready.
                return np.array(self.info_gain_matrix)
    
    #This function requires population of the UserInput doe_settings dictionary. It automatically scans many parameter modulation Permutations.
    def doeParameterModulationPermutationsScanner(self, searchType='doMetropolisHastings'):
        print("Notification: doeParameterModulationPermutationsScanner is being used which performs many MCMC runs. Convergence diagnostics are being turned off for these runs.")
        self.UserInput.user_requested_convergence_diagnostics = False
        
        import PEUQSE.CombinationGeneratorModule as CombinationGeneratorModule
        doe_settings = self.UserInput.doe_settings 
        #For the parameters, we are able to use a default one standard deviation grid if gridSamplingAbsoluteIntervalSize is a blank list.
        #doe_settings['parameter_modulation_grid_center'] #We do NOT create such a variable in user input. The initial guess variable is used, which is the center of the prior if no guess has been provided.
        parModulationGridCenterVector = self.UserInput.InputParameterInitialGuess
        numParameters = len(parModulationGridCenterVector)
        parModulationGridIntervalSizeAbsolute = doe_settings['parameter_modulation_grid_interval_size']*self.UserInput.std_prior
        parModulationGridPermutations, numPermutations = self.getGridPermutations(parModulationGridCenterVector,parModulationGridIntervalSizeAbsolute, doe_settings['parameter_modulation_grid_num_intervals'])
        self.numParModulationPermutations = numPermutations
        parModulationGridPermutations= np.array(parModulationGridPermutations)
        
        if len(self.UserInput.parameterNamesList) == len(self.UserInput.InputParametersPriorValuesUncertainties): #then we assume variable names have been provided.
            headerString = self.UserInput.stringOfParameterNames #This variable is a string, no brackets.
        else: #else no variable names have been provided.
            headerString = ''
        np.savetxt(self.UserInput.directories['logs_and_csvs']+"Info_gain__parModulationGridPermutations.csv", parModulationGridPermutations, delimiter=",", encoding =None, header=headerString)
        #We will get a separate info gain matrix for each parModulationPermutation, we'll store that in this variable.
        info_gains_matrices_list = []
        if self.UserInput.doe_settings['info_gains_matrices_multiple_parameters'] == 'each': #just making analogous structure which exists for sum.
            info_gains_matrices_lists_one_for_each_parameter = [] #make a matrix ready to copy info_gains_matrices_list. 
            #need to make a list of lists (or similar) to fill it with the individual matrices necessary.
            numParameters = len(self.UserInput.InputParametersPriorValuesUncertainties)
            for parameterIndex in range(0,numParameters):#looping across number of parameters...
                info_gains_matrices_lists_one_for_each_parameter.append([]) #These are empty lists create to indices and initialize each parameter's info_gain_matrix. They will be appended to later.
        for parModulationPermutationIndex,parModulationPermutation in enumerate(parModulationGridPermutations):                
            #####Begin PEUQSE Parallel Processing During Loop Block####
            if (self.UserInput.doe_settings['parallel_parameter_modulation'])== True:
                #We will only execute the sampling the permutationIndex matches the processor rank.
                #Additionally, if the rank is 0 and the simulation got here, it will be assumed the person is running this just to find the number of Permutations, so that will be spit out and the simulation ended.
                import PEUQSE.parallel_processing
                if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                    print("For the user input settings provided, the number of Permutations+1 will be",  numPermutations+1, ". Please use mpiexec or mpirun with this number for N. If you are not expecting to see this message, change your UserInput choices. You have chosen parallel processing for gridsearch and have run PEUQSE without mpi, which is a procedure to retrieve the number of processor ranks to use for parallelized gridsearch. A typical syntax now would be: mpiexec -n ",  numPermutations+1, " python runfile_for_your_analysis.py" )
                    sys.exit()
                elif PEUQSE.parallel_processing.currentProcessorNumber != parModulationPermutationIndex+1:
                    continue #This means the permutation index does not match the processor rank so nothing should be executed.
                #elif PEUQSE.parallel_processing.currentProcessorNumber == permutationIndex+1:
                #    pass  #This is the "normal" case and is implied, so is commented out.
            #####End PEUQSE Parallel Processing During Loop Block####
            #We will get separate info gain matrix for each parameter modulation combination.
            self.parModulationPermutationIndex = parModulationPermutationIndex #This variable is being created for parallel processing of conditions.
            info_gain_matrix = self.doeGetInfoGainMatrix(parModulationPermutation, searchType=searchType)
            #Append the info gain matrix obtainend (unless doing a parallel_conditions_exploration).
            if self.UserInput.doe_settings['parallel_conditions_exploration'] == False:
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
        ####Start block for parallel_conditions_exploration #####
        if self.UserInput.doe_settings['parallel_conditions_exploration'] == True:
          #if we're doing a parallel_conditions_exploration, we need to check if we are on the last 
          #condition exploration of the last parameter modulation. [#things could be done differently, but this works.]          
          if self.parModulationPermutationIndex+1 != self.numParModulationPermutations:
            return #this means we do nothing because it's not the final parModulation.
          elif self.parModulationPermutationIndex+1 == self.numParModulationPermutations:
            import PEUQSE.parallel_processing #Even if it's the final parModulation, need to check if it's final combination.
            if PEUQSE.parallel_processing.finalProcess == False: #not final combination.
                return
            elif PEUQSE.parallel_processing.finalProcess == True: 
                #If final parModulation and final combination, we populate self.info_gain_matrices_array
                self.consolidate_parallel_doe_info_gain_matrices() #And now we continue on with the plotting.
        ####End block for parallel_conditions_exploration #####        
        return self.info_gains_matrices_array
    
    @CiteSoft.after_call_compile_consolidated_log(compile_checkpoints=True) #This is from the CiteSoft module.
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
        #self.info_gains_matrices_array is an implied argument that usually gets populated in doeParameterModulationPermutationsScanner (when that is used).
        #Right now, when using KL_divergence and design of experiments there is an option of UserInput.doe_settings['info_gains_matrices_multiple_parameters'] = 'each' or 'sum'
        #the default is sum. But when it is 'each', then it is possible to plot separate info_gains for each parameter.
        #Note: the below code *does not* add a suffix to inidicate when a parameter Index has been fed.
        #TODO: The variable "parameterInfoGainIndex" is made with the presumption that later we'll have to add another index when we have info_gains for each parameter. In that case it will become like this:
        #xValues = self.info_gains_matrices_array[modulationIndex][:,0] will become xValues = self.info_gains_matrices_array[modulationIndex][parameterInfoGainIndex][:,0]
        #self.meshGrid_independentVariable1ValuesArray will remain unchanged.      
        
        import PEUQSE.plotting_functions as plotting_functions
        setMatPlotLibAgg(self.UserInput.plotting_ouput_settings['setMatPlotLibAgg'])
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
            if hasattr(self, 'info_gains_matrices_arrays_one_for_each_parameter'): #this structure will only exist if doeParameterModulationPermutationsScanner has been called.
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
                    if self.UserInput.doe_settings['parallel_parameter_modulation'] == False: #This is the normal case.
                        plotting_functions.makeTrisurfacePlot(xValues, yValues, zValues, figure_name = "Info_gain_TrisurfacePlot_modulation_"+str(modulationIndex+1)+plot_suffix, directory = self.UserInput.directories['graphs'])
                    if self.UserInput.doe_settings['parallel_parameter_modulation'] == True: #This is the parallel case. In this case, the actual modulationIndex to attach to the filename is given by the processor rank.
                        import PEUQSE.parallel_processing
                        plotting_functions.makeTrisurfacePlot(xValues, yValues, zValues, figure_name = "Info_gain_TrisurfacePlot_modulation_"+str(PEUQSE.parallel_processing.currentProcessorNumber)+plot_suffix, directory = self.UserInput.directories['graphs'])
            if self.info_gains_matrices_array_format == 'meshgrid':        
                for modulationIndex in range(len(local_info_gains_matrices_array)):
                    #Now need to get things prepared for the meshgrid.
                    #NOTE: we do not pull XX and YY from local_info_gains_matrices_array because that is 1D and these are 2D arrays made a different way.
                    #xValues = local_info_gains_matrices_array[modulationIndex][:,0] #Still correct, but not being used.
                    #yValues = local_info_gains_matrices_array[modulationIndex][:,1] #Still correct, but not being used.
                    XX, YY = np.meshgrid(self.meshGrid_independentVariable1ValuesArray, self.meshGrid_independentVariable2ValuesArray)
                    zValues = local_info_gains_matrices_array[modulationIndex][:,2]
                    ZZ = zValues.reshape(XX.shape) #We know from experience to reshape this way.
                    if self.UserInput.doe_settings['parallel_parameter_modulation'] == False: #This is the normal case.
                        plotting_functions.makeMeshGridSurfacePlot(XX, YY, ZZ, figure_name = "Info_gain_Meshgrid_modulation_"+str(modulationIndex+1)+plot_suffix, directory = self.UserInput.directories['graphs'])
                    if self.UserInput.doe_settings['parallel_parameter_modulation'] == True: #This is the parallel case. In this case, the actual modulationIndex to attach to the filename is given by the processor rank.
                        import PEUQSE.parallel_processing
                        plotting_functions.makeMeshGridSurfacePlot(XX, YY, ZZ, figure_name = "Info_gain_Meshgrid_modulation_"+str(PEUQSE.parallel_processing.currentProcessorNumber)+plot_suffix, directory = self.UserInput.directories['graphs'])
        else:
            print("At present, createInfoGainPlots and createInfoGainModulationPlots only create plots when the length of  independent_variable_grid_center is 2. We don't currently support creation of other dimensional plots. The infogain data is being exported into the file _____.csv")
    def getLogP(self, proposal_sample): #The proposal sample is specific parameter vector.
        log_prior_proposal = self.getLogPrior(proposal_sample)
        if log_prior_proposal == float('-inf'): #due to the bounds check or other reasons, log_prior_proposal could be '-inf'.
            return float('-inf')
        #else continue as normal.
        log_likelihood_proposal, simulationOutput_proposal = self.getLogLikelihood(proposal_sample)
        logP = log_likelihood_proposal+log_prior_proposal #Of the Metropolis-Hastings accept/reject ratio
        return logP
        
    def getNegLogP(self, proposal_sample): #The proposal sample is specific parameter vector. We are using negative of log P because scipy optimize doesn't do maximizing. It's recommended minimize the negative in this situation.
        neg_log_postererior = -1*self.getLogP(proposal_sample)
        return neg_log_postererior

    #this is just a wrapper around doOptimizeNegLogP
    def doOptimizeLogP(self, simulationFunctionAdditionalArgs = (), method = None, optimizationAdditionalArgs = {}, printOptimum = True, verbose=True, maxiter=0):
        if printOptimum == True:
            print("doOptimizeLogP is a wrapper that calls doOptimizeNegLogP. The final results from doOptimizeNegLogP will be printed.")
        [self.map_parameter_set, negLogP] = self.doOptimizeNegLogP(simulationFunctionAdditionalArgs = simulationFunctionAdditionalArgs, method = method, optimizationAdditionalArgs = optimizationAdditionalArgs, printOptimum = printOptimum, verbose=verbose, maxiter=maxiter)
        self.map_logP = -1.0*negLogP
        return [self.map_parameter_set, self.map_logP]

    find_MAP = doOptimizeLogP # pointer to doOptimizeLogP to make our user options similar to existing popular BPE softwares (ie. pymc)

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
        negLogP = optimizeResult.fun
        self.map_logP = -1.0*negLogP #This is the map logP
        if printOptimum == True:
            print("Final results from doOptimizeNegLogP:", self.map_parameter_set, "final negLogP:", negLogP, "final logP:", self.map_logP)
        return [self.map_parameter_set, negLogP]


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
        if self.parameterBoundsOn:#Will force L-BFGS-B because it has bounds.
            if method != 'L-BFGS-B':
                print("Notification: Parameter bounds are on and doOptimizeSSR is being called. Forcing the optimization method to be L-BFGS-B because this it the only SSR method presently allowed for bounds with PEUQSE.")
                method = 'L-BFGS-B'
                zippedBounds = list(zip(self.UserInput.InputParameterPriorValues_lowerBounds, self.UserInput.InputParameterPriorValues_upperBounds))
                optimizationAdditionalArgs['bounds'] = zippedBounds
        
        initialGuess = self.UserInput.InputParameterInitialGuess
        import scipy.optimize
        if verbose == False:
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(self.getSSR, initialGuess, method = method, **optimizationAdditionalArgs)
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(self.getSSR, initialGuess, method = method, options={"maxiter": maxiter}, **optimizationAdditionalArgs)
        if verbose == True:
            verbose_simulator = verbose_optimization_wrapper(self.getSSR)
            if maxiter == 0:
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"disp": True}, **optimizationAdditionalArgs)
            if maxiter != 0:    
                optimizeResult = scipy.optimize.minimize(verbose_simulator.simulateAndStoreObjectiveFunction, initialGuess, method=method, callback=verbose_simulator.callback, options={"maxiter": maxiter}, **optimizationAdditionalArgs)
            #print(f"Number of calls to Simulator instance {verbose_simulator.num_calls}") <-- this is the same as the "Function evaluations" field that gets printed.
            
        self.opt_parameter_set = optimizeResult.x #This is the best fit parameter set.
        self.opt_SSR = optimizeResult.fun #This is the best fit SSR.
        if printOptimum == True:
            print("Final results from doOptimizeSSR:", self.opt_parameter_set, "final SSR:", self.opt_SSR, "final negSSR:", -1*self.opt_SSR)
        #TODO: Right now, the createAllPlots command will not work unless we populate the map parameter set, so that is what we are doing. But a better longterm solution needs to be made. In which the graph says "opt" rather than "MAP" and uses the appropriate variables.
        #TODO: Also need to add things like WSSR based on magnitude and variance weightings.
        self.map_parameter_set = self.opt_parameter_set
        return [self.opt_parameter_set, self.opt_SSR]

    
    #This function is meant to be called from the runfile when testing a new function etc. It allows a simulation plot to be created.
    #This is *not* recommended for use in other functions, where it is recommended that getLogP be called directly.
    def doSinglePoint(self, discreteParameterVector=None, objectiveFunction='logP'):
        #objectiveFunction can be 'logP' or 'SSR'
        if type(discreteParameterVector)!=type(None) and (self.reducedParameterSpaceOn): #if reduced parameter space is on, and the user is providing an discreteParameterVector, we that is supposed to be the full parameter vector so we will reduce it. 
            reducedIndices = self.UserInput.model['reducedParameterSpace']
            discreteParameterVector = returnReducedIterable(discreteParameterVector, reducedIndices)
        if type(discreteParameterVector)==type(None): #If somebody did not feed a specific vector, we take the initial guess. For the case of a reducedParameterSpace, this internal variable is already reduced.
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

    def calculateInfoGain(self):
        if self.UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] == 'KL_divergence':
            try: #If the info_gain_KL fails, we will return info_gain_log_ratio. 
                #Below is the KL_divergence info_gain calculation.
                from warnings import catch_warnings, simplefilter #Used to surpress warnings from 0 probability density in current_info_gain_KL calculation
                length, width = self.post_burn_in_samples.shape
                self.info_gain_KL = 0
                self.info_gain_KL_each_parameter  = []
                for param in range(width):
                    import matplotlib.pyplot as plt #FIXME: #TODO: this plotting needs to be moved into the plotting area and as optinoal.
                    (density0,bins0,pathces0)=plt.hist([self.samples_of_prior[:,param].flatten(),self.post_burn_in_samples[:,param].flatten()],bins=100,density=True)
                    plt.close()
                    # the following code handles surpressing the RuntimeWarning displayed when dealing with a 0 probability density in calculating current_info_gain_KL. It does not affect the calculation at all.
                    with catch_warnings():
                        simplefilter("ignore")
                        current_info_gain_KL = density0[1]*np.log(density0[1]/density0[0])
                    current_info_gain_KL = current_info_gain_KL[np.isfinite(current_info_gain_KL)]
                    current_info_gain_KL = np.sum(current_info_gain_KL)
                    self.info_gain_KL_each_parameter.append(current_info_gain_KL) #could make this optional, but normally shouldn't take much memory.
                    self.info_gain_KL = self.info_gain_KL + current_info_gain_KL
                self.info_gain_each_parameter = self.info_gain_KL_each_parameter #could make this optional, but normally shouldn't take much memory.
                self.info_gain = self.info_gain_KL
            except:
                print("unable to calculate KL_divergence info_gain. Calculating log_ratio info_gain.")
                self.UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] = 'log_ratio'
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
                #TODO: Below there are some "if verbose", but those should not be printed, they should be collected and exported to a file at the end.
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
            self.info_gain = self.info_gain_log_ratio
        return self.info_gain    


    def exportSingleConditionInfoGainMatrix(self, parameterPermutationNumber, conditionsPermutationAndInfoGain, conditionsPermutationIndex):
        #Note that parameterPermutationNumber is parameterPermutationIndex+1
        file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix() #Rather self explanatory.
        file_name_suffix=file_name_suffix[1:] #removing the '_' that comes by default, we will add the '_' back in later below.
        if int(conditionsPermutationIndex+1) != int(file_name_suffix):
            print("line 1199: There is a problem in the parallel processing of conditions info gain matrix calculation!", conditionsPermutationIndex+1, file_name_suffix)
        #I am commenting out the below line because the savetxt was causing amysterious "no such file or directory" error.
        np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'conditionsPermutationAndInfoGain_'+'mod'+str(int(parameterPermutationNumber))+'_cond'+str(conditionsPermutationIndex+1)+'.csv',conditionsPermutationAndInfoGain, delimiter=",")
        pickleAnObject(conditionsPermutationAndInfoGain, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'conditionsPermutationAndInfoGain_'+'mod'+str(int(parameterPermutationNumber))+'_cond'+str(conditionsPermutationIndex+1))
        
    #This function will calculate MAP and mu_AP, evidence, and related quantities.
    def calculatePostBurnInStatistics(self, calculate_post_burn_in_log_priors_vec = False):       
        #First need to create priors if not already there, because ESS does not store priors during the run (MH does).
        if not hasattr(self, 'post_burn_in_log_priors_vec'): 
            calculate_post_burn_in_log_priors_vec = True 
        if calculate_post_burn_in_log_priors_vec == True:
                #TODO: change below to use numpy vectorize. It will probably be faster then.
                #Below line following a line from https://github.com/threeML/threeML/blob/master/threeML/bayesian/zeus_sampler.py
                self.post_burn_in_log_priors_vec = np.array([self.getLogPrior(parameterPermutation) for parameterPermutation in self.post_burn_in_samples])
                self.post_burn_in_log_priors_vec = np.atleast_2d(self.post_burn_in_log_priors_vec).transpose()
        #Next need to apply filtering before getting statistics.
        filterSamples = bool(self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_samples'])
        filterCoeffient = self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_coefficient']
        if type(filterCoeffient) == type("string"):
            if filterCoeffient.lower() == "auto":
                filterCoeffient = 2.0
        if filterSamples == True:   
            #before filtering, we will keep an unfiltered version in case of ['exportAllSimulatedOutputs'] == True:
            if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True:
                try: #This try and except is primarily because as of Dec 6th 2020, the feature has been implemented for MH but not ESS. With ESS, post_burn_in_log_priors_vec_unfiltered is not gauranteed.
                    self.post_burn_in_samples_unfiltered = copy.deepcopy(self.post_burn_in_samples)
                    self.post_burn_in_log_posteriors_un_normed_vec_unfiltered = copy.deepcopy(self.post_burn_in_log_posteriors_un_normed_vec)
                    self.post_burn_in_log_priors_vec_unfiltered = copy.deepcopy(self.post_burn_in_log_priors_vec)
                except:
                    pass
            originalLength = np.shape(self.post_burn_in_log_posteriors_un_normed_vec)[0] 
            try:
                mergedArray = np.hstack( (self.post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_priors_vec, self.post_burn_in_samples) )
            except:
                print("Line 866: There has been an error, here are post_burn_in_log_posteriors_un_normed_vec, post_burn_in_samples, post_burn_in_log_priors_vec", np.shape(self.post_burn_in_log_posteriors_un_normed_vec), np.shape(self.post_burn_in_samples), np.shape(self.post_burn_in_log_priors_vec))
                # print(self.post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_samples, self.post_burn_in_log_priors_vec)
                sys.exit()
            #Now need to find cases where the probability is too low and filter them out.
            #Filtering Step 1: Find average and Stdev of log(-logP)
            logNegLogP = np.log(-1*self.post_burn_in_log_posteriors_un_normed_vec)
            meanLogNegLogP = np.mean(logNegLogP) 
            stdLogNegLogP = np.std(logNegLogP)
            maxLogP = np.max(self.post_burn_in_log_posteriors_un_normed_vec) #this is the MAP logP.
            #getting the mu_AP probability for filtering requires a couple of steps. We name it carefully because the mu_AP willl change after filtering.
            mu_AP_parameter_set_unfiltered_data = np.mean(self.post_burn_in_samples, axis=0) #now we are going to get the mu_AP of the unfiltered data.
            mu_AP_log_P_unfiltered_data = self.getLogP(mu_AP_parameter_set_unfiltered_data) #this is the mu_AP logP
            if self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_benchmark'] == 'mu_AP':
                #This benchmark is relative to the LogP values themselves
                filteringThreshold = mu_AP_log_P_unfiltered_data - filterCoeffient # filter values below threshold starting at the mu_AP or mu_AP proxy.
                removeValuesDirection = 'below'
                #Now, call the function I have made for filtering by deleting the rows above/below a certain value
                truncatedMergedArray = arrayThresholdFilter(mergedArray, filterKey=self.post_burn_in_log_posteriors_un_normed_vec, thresholdValue=filteringThreshold, removeValues = removeValuesDirection, transpose=False)
            elif self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_benchmark'] == 'MAP':
                #This benchmark is relative to the LogP values themselves
                filteringThreshold = maxLogP - filterCoeffient # filter values below threshold starting at the MAP
                removeValuesDirection = 'below'
                #Now, call the function I have made for filtering by deleting the rows above/below a certain value
                truncatedMergedArray = arrayThresholdFilter(mergedArray, filterKey=self.post_burn_in_log_posteriors_un_normed_vec, thresholdValue=filteringThreshold, removeValues = removeValuesDirection, transpose=False)
            elif self.UserInput.parameter_estimation_settings['mcmc_threshold_filter_benchmark'] == 'auto':
                #This benchmark is unusual be cause it it is related to logNegLogP
                filteringThreshold = meanLogNegLogP+filterCoeffient*stdLogNegLogP #Threshold. 
                removeValuesDirection = 'above'
                #Now, call the function I have made for filtering by deleting the rows above/below a certain value
                truncatedMergedArray = arrayThresholdFilter(mergedArray, filterKey=logNegLogP, thresholdValue=filteringThreshold, removeValues = removeValuesDirection, transpose=False)
            self.post_burn_in_log_posteriors_un_normed_vec = np.atleast_2d(truncatedMergedArray[:,0]).transpose()
            self.post_burn_in_log_priors_vec = np.atleast_2d(truncatedMergedArray[:,1]).transpose()
            self.post_burn_in_samples = truncatedMergedArray[:,2:]
        #Map calculation etc. is intentionally placed below the filtering so that the map_index is assigned correctly per the final values.
        self.mu_AP_parameter_set = np.mean(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.
        self.stdap_parameter_set = np.std(self.post_burn_in_samples, axis=0) #This is the mean of the posterior, and is the point with the highest expected value of the posterior (for most distributions). For the simplest cases, map and mu_AP will be the same.            
        map_logP = max(self.post_burn_in_log_posteriors_un_normed_vec)
        self.map_logP = map_logP
        self.map_index = list(self.post_burn_in_log_posteriors_un_normed_vec).index(map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
        self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.            
        #TODO: Probably should return the variance of each sample in the post_burn_in
        #posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        self.evidence = np.mean(np.exp(self.post_burn_in_log_posteriors_un_normed_vec))/np.linalg.norm(self.post_burn_in_samples)
        self.info_gain = self.calculateInfoGain()
        if self.UserInput.parameter_estimation_settings['verbose'] == True:
            print("map_parameter_set ", self.map_parameter_set)
            print("mu_AP_parameter_set ", self.mu_AP_parameter_set)
            print("stdap_parameter_set ",self.stdap_parameter_set)
        return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]

    #This function gets the prefix and suffix for saving files when doing ParallelProcessing with MPI.
    #Importantly, the **directory** of the parallel processing is included as part of the prefix.
    def getParallelProcessingPrefixAndSuffix(self):
        file_name_prefix = ''
        file_name_suffix = ''
        directory_name_suffix = ''
        if (self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or self.UserInput.parameter_estimation_settings['multistart_parallel_sampling'] or self.UserInput.doe_settings['parallel_conditions_exploration']) == True: 
            import PEUQSE.parallel_processing
            import os
            if PEUQSE.parallel_processing.currentProcessorNumber == 0:
                pass
            if PEUQSE.parallel_processing.currentProcessorNumber > 0:                
                file_name_suffix = "_"+str(PEUQSE.parallel_processing.currentProcessorNumber)
                file_name_prefix = ""  
                directory_name_suffix = "mpi_cached_files/" #TODO: FIX THIS, IT MAY NOT WORK ON EVERY OS. SHOULD USE 'os' MODULE TO FIND DIRECTION OF THE SLASH OR TO DO SOMETHING SIMILAR. CURRENTLY IT IS WORKING ON MY WINDOWS DESPITE BEING "/"
        return file_name_prefix, file_name_suffix, directory_name_suffix

    def exportPostBurnInStatistics(self):
        #TODO: Consider to Make header for mcmc_samples_array. Also make exporting the mcmc_samples_array optional. 
        file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix() #Rather self explanatory.
        mcmc_samples_array = np.hstack((self.post_burn_in_log_posteriors_un_normed_vec,self.post_burn_in_samples))
        np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_logP_and_parameter_samples'+file_name_suffix+'.csv',mcmc_samples_array, delimiter=",")
        pickleAnObject(mcmc_samples_array, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_logP_and_parameter_samples'+file_name_suffix)
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True: #By default, we should not keep this, it's a little too large with large sampling.
            try: #The main reason to use a try and except is because this feature has not been implemented for ESS. With ESS, the mcmc_unfiltered_post_burn_in_simulated_outputs would not be retained and the program would crash.
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_unfiltered_post_burn_in_simulated_outputs'+file_name_suffix+'.csv',self.post_burn_in_samples_simulatedOutputs, delimiter=",")         
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_unfiltered_post_burn_in_parameter_samples'+file_name_suffix+'.csv',self.post_burn_in_samples_unfiltered, delimiter=",")            
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_unfiltered_post_burn_in_log_priors_vec'+file_name_suffix+'.csv',self.post_burn_in_log_posteriors_un_normed_vec_unfiltered, delimiter=",")            
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_unfiltered_post_burn_in_log_posteriors_un_normed_vec'+file_name_suffix+'.csv',self.post_burn_in_log_priors_vec_unfiltered, delimiter=",")                        
            except:
                pass
        with open(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_log_file'+file_name_suffix+".txt", 'w') as out_file:
            out_file.write("self.initial_point_parameters:" + str( self.UserInput.InputParameterInitialGuess) + "\n")
            out_file.write("MAP_logP:" +  str(self.map_logP) + "\n")
            out_file.write("self.map_parameter_set:" + str( self.map_parameter_set) + "\n")
            out_file.write("self.map_index:" +  str(self.map_index) + "\n")
            out_file.write("self.mu_AP_parameter_set:" + str( self.mu_AP_parameter_set) + "\n")
            out_file.write("self.stdap_parameter_set:" + str( self.stdap_parameter_set) + "\n")
            out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
            out_file.write("evidence:" + str(self.evidence) + "\n")
            out_file.write("posterior_cov_matrix:" + "\n" + str(np.cov(self.post_burn_in_samples.T)) + "\n")
            if (self.permutation_and_doOptimizeNegLogP == True) or (self.permutation_and_doOptimizeLogP == True):
                out_file.write("\n WARNING: It appears this run used a doOptimize with multi-start. In this case, the MAP_logP and map_parameter_set are the optima.  However, the mu_AP_parameter_set and stdap_parameter_set are not meaningful, since this was not an even weighted exploration of the posterior.")
            relativeDifferenceArray = abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.std_prior)
            if np.array(relativeDifferenceArray > 0.10).any() == True:
                out_file.write("Convergence Status: Warning - The difference between the MAP parameter set and mu_AP parameter set is greater than 10% of the prior standard deviations in at least one parameter. This means that either your posterior distribution is asymmetric or that it is not yet converged. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632")
            else:
                out_file.write("Convergence Status: The difference between the MAP parameter set and mu_AP parameter set is less than 10% of the prior standard deviations for each parameter. There is no general method for knowing when the correct solution has been obtained, but this small difference is one sign of a converged solution.")
        postBurnInStatistics = [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
        if hasattr(self, 'during_burn_in_samples'):
            pickleAnObject(np.hstack((self.during_burn_in_log_posteriors_un_normed_vec, self.during_burn_in_samples)),self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_burn_in_logP_and_parameter_samples'+file_name_suffix)
            np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'mcmc_burn_in_logP_and_parameter_samples'+file_name_suffix+'.csv',np.hstack((self.during_burn_in_log_posteriors_un_normed_vec, self.during_burn_in_samples)), delimiter=",")
        pickleAnObject(postBurnInStatistics,self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_post_burn_in_statistics'+file_name_suffix)
        pickleAnObject(self.map_logP,self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_map_logP'+file_name_suffix)
        if self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] == True:
            pickleAnObject(self.UserInput.InputParameterInitialGuess,self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_continued_initial_point_parameters'+file_name_suffix)
        else:
            pickleAnObject(self.UserInput.InputParameterInitialGuess,self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_initial_point_parameters'+file_name_suffix)
        if hasattr(self, 'mcmc_last_point_sampled'):
            pickleAnObject(self.mcmc_last_point_sampled, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'mcmc_last_point_sampled'+file_name_suffix)

    #This function is modelled after exportPostBurnInStatistics. That is why it has the form that it does.
    def exportPostPermutationStatistics(self, searchType=''): #if it is an mcmc run, then we need to save the sampling as well.
        #TODO: Consider to Make header for mcmc_samples_array. Also make exporting the mcmc_samples_array optional. 
        file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix() #Rather self explanatory.
        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): #Note: this might be needed for parallel processing, not sure.
            mcmc_samples_array = np.hstack((self.post_burn_in_log_posteriors_un_normed_vec,self.post_burn_in_samples))
            np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_logP_and_parameter_samples'+file_name_suffix+'.csv',mcmc_samples_array, delimiter=",")
            pickleAnObject(mcmc_samples_array, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'permutation_logP_and_parameter_samples'+file_name_suffix)
        if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs'] == True: #By default, we should not keep this, it's a little too large with large sampling.
            try: #The main reason to use a try and except is because this feature has not been implemented for ESS. With ESS, the mcmc_unfiltered_post_burn_in_simulated_outputs would not be retained and the program would crash.
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_unfiltered_post_burn_in_simulated_outputs'+file_name_suffix+'.csv',self.post_burn_in_samples_simulatedOutputs, delimiter=",")            
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_unfiltered_post_burn_in_parameter_samples'+file_name_suffix+'.csv',self.post_burn_in_samples_unfiltered, delimiter=",")            
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_unfiltered_post_burn_in_log_priors_vec'+file_name_suffix+'.csv',self.post_burn_in_log_posteriors_un_normed_vec_unfiltered, delimiter=",")            
                np.savetxt(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_unfiltered_post_burn_in_log_posteriors_un_normed_vec'+file_name_suffix+'.csv',self.post_burn_in_log_priors_vec_unfiltered, delimiter=",")                        
            except:
                pass
        with open(self.UserInput.directories['logs_and_csvs']+directory_name_suffix+file_name_prefix+'permutation_log_file'+file_name_suffix+".txt", 'w') as out_file:
            out_file.write("self.initial_point_parameters:" + str( self.UserInput.InputParameterInitialGuess) + "\n")
            out_file.write("MAP_logP:" +  str(self.map_logP) + "\n")
            out_file.write("self.map_parameter_set:" + str( self.map_parameter_set) + "\n")
            if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'): #Below are only for mcmc_sampling
                out_file.write("self.map_index:" +  str(self.map_index) + "\n")
                out_file.write("self.mu_AP_parameter_set:" + str( self.mu_AP_parameter_set) + "\n")
                out_file.write("self.stdap_parameter_set:" + str( self.stdap_parameter_set) + "\n")
                out_file.write("self.info_gain:" +  str(self.info_gain) + "\n")
                out_file.write("evidence:" + str(self.evidence) + "\n")
                out_file.write("posterior_cov_matrix:" + "\n" + str(np.cov(self.post_burn_in_samples.T)) + "\n")
                relativeDifferenceArray = abs((self.map_parameter_set - self.mu_AP_parameter_set)/self.UserInput.std_prior)
                if np.array(relativeDifferenceArray > 0.10).any() == True:
                    out_file.write("Convergence Status: Warning - The difference between the MAP parameter set and mu_AP parameter set is greater than 10% of the prior standard deviations in at least one parameter. This means that either your posterior distribution is asymmetric or that it is not yet converged. This may mean that you need to increase your mcmc_length, increase or decrease your mcmc_relative_step_length, or change what is used for the model response.  There is no general method for knowing the right  value for mcmc_relative_step_length since it depends on the sharpness and smoothness of the response. See for example https://www.sciencedirect.com/science/article/pii/S0039602816300632")
                else:
                    out_file.write("Convergence Status: The difference between the MAP parameter set and mu_AP parameter set is less than 10% of the prior standard deviations for each parameter. There is no general method for knowing when the correct solution has been obtained, but this small difference is one sign of a converged solution.")
        if (searchType == 'doEnsembleSliceSampling') or (searchType == 'doMetropolisHastings') or (searchType == 'doEnsembleModifiedMHSampling'):
            postBurnInStatistics = [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]
            pickleAnObject(postBurnInStatistics, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'permutation_post_burn_in_statistics'+file_name_suffix)
        pickleAnObject(self.map_logP, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'permutation_map_logP'+file_name_suffix)
        pickleAnObject(self.map_parameter_set, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'permutation_map_parameter_set'+file_name_suffix)
        pickleAnObject(self.UserInput.InputParameterInitialGuess, self.UserInput.directories['pickles']+directory_name_suffix+file_name_prefix+'permutation_initial_point_parameters'+file_name_suffix)
        
    def getConvergenceDiagnostics(self, discrete_chains_post_burn_in_samples=[], showFigure = None):
        """
        Guides Integrated Autocorrelation Time and Geweke convergence analysis and makes plots.

        :param samplingFunctionstr (optional): String to define the sampler. (:type: str)
        :param discrete_chains_post_burn_in_samples (optional): Array that contains post burn in samples. Shape is (numSamples, numChains, numParams) (:type np.array)
        """
        if showFigure == None: showFigure = True
        if (len(discrete_chains_post_burn_in_samples)==0) and (hasattr(self, 'discrete_chains_post_burn_in_samples')):
            discrete_chains_post_burn_in_samples = self.discrete_chains_post_burn_in_samples
        # check if inputted array is a numpy array.
        if type(discrete_chains_post_burn_in_samples).__module__ != np.__name__:
            print('The input array needs to a numpy array. If a list was inputted, wrap the variable in np.array(var)')
            sys.exit()
        # check if array shape is right.
        if len(discrete_chains_post_burn_in_samples.shape) != 3:
            print('The input array needs to take the shape (numSamples, numChains, numParameters). The current inputs does not have the appropriate dimensions.')
            sys.exit()
       

        # convergence diagnostic function is called to run calculations and create plots
        # outputs are saved as a tuple 
        #TODO: make its own plot settings in different commit. 
        try:
            convergence_ouputs = calculateAndPlotConvergenceDiagnostics(discrete_chains_post_burn_in_samples, self.UserInput.model['parameterNamesAndMathTypeExpressionsDict'], self.UserInput.scatter_matrix_plots_settings, self.UserInput.directories['graphs'], showFigure = showFigure)
        except Exception as theError:
            print("Warning: Unable to calculate and plot convergence diagnostics. The error was:", theError)
            return None #this is to end the function, so it does not crash below.
            

        # initialize and populate convergence dictionary in class object
        self.convergence = {}
        self.convergence['AutoCorrTime'] = {}
        self.convergence['Geweke'] = {}
        # unpack tuple to save convergence information to self
        self.convergence['AutoCorrTime']['window_indices'] = convergence_ouputs[0]
        self.convergence['AutoCorrTime']['final_parameter_values'] = convergence_ouputs[1]
        self.convergence['AutoCorrTime']['parameter_act_for_each_window'] = convergence_ouputs[2]
        self.convergence['Geweke']['window_indices'] = convergence_ouputs[3]
        self.convergence['Geweke']['final_combined_parameter_values'] = convergence_ouputs[4]
        self.convergence['Geweke']['final_combined_parameter_percent_outlier'] = convergence_ouputs[5]

        
    #Our EnsembleModifiedMHSampling is done by the emcee back end. (pip install emcee)
    software_name = "emcee"
    software_version = "3.1.2"
    software_unique_id = "https://github.com/dfm/emcee"
    software_kwargs = {"version": software_version, "author": ['Foreman-Mackey, D.', 'Hogg, D.~W.', 'Lang, D.', 'Goodman, J.'], "cite": ["@article{emcee, author = {{Foreman-Mackey}, D. and {Hogg}, D.~W. and {Lang}, D. and {Goodman}, J.}, title = {emcee: The MCMC Hammer}, journal = {PASP}, year = 2013, volume = 125, pages = {306-312}, eprint = {1202.3665}, doi = {10.1086/670067}}"] }
    @CiteSoft.function_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def doEnsembleModifiedMHSampling(self, mcmc_nwalkers_direct_input = None, walkerInitialDistribution='UserChoice', walkerInitialDistributionSpread='UserChoice', calculatePostBurnInStatistics=True, mcmc_exportLog ='UserChoice', continueSampling='auto'):
        """
        TODO: make params and return definitions along with a short description.
        """
        # import the emcee backend module
        import emcee
        # handle the initial distribution spread with the UserInputs
        if walkerInitialDistribution == 'UserChoice':
            walkerInitialDistribution = self.UserInput.parameter_estimation_settings['mcmc_walkerInitialDistribution']
        if walkerInitialDistribution.lower() == 'auto':
            walkerInitialDistribution = 'sobol'
        if str(walkerInitialDistributionSpread) == 'UserChoice':
            walkerInitialDistributionSpread = self.UserInput.parameter_estimation_settings['mcmc_walkerInitialDistributionSpread']
        if str(walkerInitialDistributionSpread).lower() == 'auto':
            walkerInitialDistributionSpread = 0.866 #This choice is intended to be useful for uniform distribution priors, as described in the UserInput file.
        # First check whether the UserInput already has been used to run another MCMC run, and if so make sure it is the same algorithm type
        if hasattr(self.UserInput, 'last_MCMC_run_type'):
            if self.UserInput.last_MCMC_run_type == 'EJS':
                pass
            else:
                print("ERROR: Two different algorithms of MCMC runs were attempted to be run within the same python instance without reloading the UserInput. To run two MCMC runs of different algorithms in the same python instance is not supported, but can be accomplished by using imp.reload(UserInput) followed by imp.reload(PEUQSE) between changing MCMC algorithms, and which requires re-populating the UserInput with your choices. "); sys.exit()

        #Check if we need to continue sampling, and prepare for it if we need to.
        if continueSampling == 'auto':
            if ('mcmc_continueSampling' not in self.UserInput.parameter_estimation_settings) or self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] == 'auto': #check that UserInput does not overrule the auto.
                if hasattr(self, 'mcmc_last_point_sampled'): #if we have an existing mcmc_last_point_sampled in the object, we will assume more sampling is desired.
                    continueSampling = True
                else:
                    continueSampling = False
            else: continueSampling = self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] 
        if continueSampling == True:
            if hasattr(self, 'mcmc_last_point_sampled'): #If we are continuing from an old mcmc in this object.
                self.last_post_burn_in_log_posteriors_un_normed_vec = copy.deepcopy(self.post_burn_in_log_posteriors_un_normed_vec)
                self.last_post_burn_in_samples = copy.deepcopy(self.post_burn_in_samples)
            else: #Else we need to read from the file.                
                #First check if we are doing some kind of parallel sampling, because in that case we need to read from the file for our correct process rank. We put that info into the prefix and suffix.
                file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix()
                self.last_logP_and_parameter_samples_filename = file_name_prefix + "mcmc_logP_and_parameter_samples" + file_name_suffix
                self.last_logP_and_parameter_samples_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_logP_and_parameter_samples_filename)
                self.last_post_burn_in_log_posteriors_un_normed_vec =  np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,0]))  #First column is the logP
                if np.shape(self.last_post_burn_in_log_posteriors_un_normed_vec)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_log_posteriors_un_normed_vec = self.last_post_burn_in_log_posteriors_un_normed_vec.transpose()
                self.last_post_burn_in_samples =   np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,1:])) #later columns are the samples.
                if np.shape(self.last_post_burn_in_samples)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_samples = self.last_post_burn_in_samples.transpose()
                self.mcmc_last_point_sampled_filename = file_name_prefix + "mcmc_last_point_sampled" + file_name_suffix
                self.mcmc_last_point_sampled_data = unpickleAnObject(self.UserInput.directories['pickles']+self.mcmc_last_point_sampled_filename)
                self.mcmc_last_point_sampled = self.mcmc_last_point_sampled_data
                self.last_InputParameterInitialGuess_filename = file_name_prefix + "mcmc_initial_point_parameters" + file_name_suffix
                self.last_InputParameterInitialGuess_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_InputParameterInitialGuess_filename)
                self.UserInput.InputParameterInitialGuess = self.last_InputParameterInitialGuess_data #populating this because otherwise non-grid Multi-Start will get the wrong values exported. & Same for final plots.
        ####these variables need to be made part of UserInput####
        numParameters = len(self.UserInput.mu_prior) #This is the number of parameters.
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if isinstance(self.UserInput.parameter_estimation_settings['mcmc_random_seed'], int): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        if isinstance(mcmc_nwalkers_direct_input, type(None)): #This is the normal case.
            if 'mcmc_nwalkers' not in self.UserInput.parameter_estimation_settings: self.mcmc_nwalkers = 'auto'
            else: self.mcmc_nwalkers = self.UserInput.parameter_estimation_settings['mcmc_nwalkers']
            if isinstance(self.mcmc_nwalkers, str): 
                if self.mcmc_nwalkers.lower() == "auto":
                    self.mcmc_nwalkers = numParameters*16 # according to zeus paper (https://doi.org/10.1093/mnras/stab2867), 16*D is the optimal number of walkers for emcee
                else: #else it is an integer, or a string meant to be an integer.
                    self.mcmc_nwalkers =  int(self.mcmc_nwalkers)
        else: #this is mainly for PermutationSearch which will (by default) use the minimum number of walkers per point.
            self.mcmc_nwalkers = int(mcmc_nwalkers_direct_input)
        if (self.mcmc_nwalkers%2) != 0: #Check that it's even. If not, add one walker.
            print("The EnsembleJumpSampling requires an even number of Walkers. Adding one Walker.")
            self.mcmc_nwalkers = self.mcmc_nwalkers + 1
        requested_mcmc_steps = self.UserInput.parameter_estimation_settings['mcmc_length']
        nEnsembleSteps = int(requested_mcmc_steps/self.mcmc_nwalkers) #We calculate the calculate number of the Ensemble Steps from the total sampling steps requested divided by self.mcmc_nwalkers.
        if nEnsembleSteps == 0:
            nEnsembleSteps = 1
        if str(self.UserInput.parameter_estimation_settings['mcmc_burn_in']).lower() == 'auto': self.mcmc_burn_in_length = int(requested_mcmc_steps*0.1) # burn in is relative to entire ensemble
        else: self.mcmc_burn_in_length = self.UserInput.parameter_estimation_settings['mcmc_burn_in']
        if 'mcmc_maxiter' not in self.UserInput.parameter_estimation_settings: mcmc_maxiter = 1E6 #The default from zeus is 1E4, but I have found that is not always sufficient.
        else: mcmc_maxiter = self.UserInput.parameter_estimation_settings['mcmc_maxiter']
        ####end of user input variables####
        #now to do the mcmc
        # if the initial guess is a string, read its pickle file for starting points for walkers
        # if this read array is a single point, proceed as normal
        # if this read array is multiple points, make sure that it has the same shape as number of walkers
        # else, print error and exit
        loaded_initial_guess_flag = False # flag to switch allow normal generation of initial points to guided if the loaded points are a single point
        if isinstance(self.UserInput.model['InputParameterInitialGuess'], str): # if yes, then assume a pickle file is imported
            # try to read from current directory and pickle directory
            start_point_pkl_file_name = self.UserInput.model['InputParameterInitialGuess']
            if '.pkl' in start_point_pkl_file_name: # remove '.pkl' from string if it is there
                start_point_pkl_file_name = start_point_pkl_file_name.replace('.pkl', '')
            from os.path import exists
            # check if file exists in current working directory
            if exists(start_point_pkl_file_name + '.pkl'):
                walkerStartPoints = unpickleAnObject(start_point_pkl_file_name)
                loaded_initial_guess_flag = True # makes sure the pickled array is used for starting points
            # check if file exists in pickle directory
            elif exists(self.UserInput.directories['pickles'] + start_point_pkl_file_name + '.pkl'):
                walkerStartPoints = unpickleAnObject(self.UserInput.directories['pickles'] + start_point_pkl_file_name)
                loaded_initial_guess_flag = True # makes sure the pickled array is used for starting points
            else:
                print('The pickled object for initial guess points must exist in base directory or pickles directory.')
                sys.exit()
            # make sure the initial guess is a numpy array
            walkerStartPoints = np.array(walkerStartPoints)
            # check if start points have the same amount of walkers
            initial_guess_mcmc_nwalkers = walkerStartPoints.shape[0] # last points have shape (nwalkers, nparams)
            if len(walkerStartPoints.shape) == 1:
                loaded_initial_guess_flag = False # set this to false so initial points can be generated
            elif initial_guess_mcmc_nwalkers != self.mcmc_nwalkers:
                print(f'Initial guess walker number must be the same shape as mcmc_nwalkers: {self.mcmc_nwalkers} mcmc walkers but initial guess had {initial_guess_mcmc_nwalkers} walkers')
                sys.exit()
        if continueSampling == False:
            # check if initial guess is loaded
            if not loaded_initial_guess_flag:
                walkerStartPoints = self.generateInitialPoints(initialPointsDistributionType=walkerInitialDistribution, numStartPoints = self.mcmc_nwalkers,relativeInitialDistributionSpread=walkerInitialDistributionSpread) #making the first set of starting points.
        elif continueSampling == True:
            # make burn in samples 0 since burn in has already occurred
            self.mcmc_burn_in_length = 0
            # start points have the same amount of walkers
            previous_mcmc_nwalkers = self.mcmc_last_point_sampled.shape[0] # last points have shape (nwalkers, nparams)
            if previous_mcmc_nwalkers != self.mcmc_nwalkers:
                print(f'Setting walkers to be the same as the previous run for continue sampling: from {self.mcmc_nwalkers} to {previous_mcmc_nwalkers} walkers')
                self.mcmc_nwalkers = previous_mcmc_nwalkers
            # starts from the last set of chains
            walkerStartPoints = self.mcmc_last_point_sampled 
        emcee_sampler = emcee.EnsembleSampler(self.mcmc_nwalkers, numParameters, log_prob_fn=self.getLogP)    
        for trialN in range(0,1000):#Todo: This number of this range is hardcoded but should probably be a user selection.
            try:
                emcee_sampler.run_mcmc(walkerStartPoints, nEnsembleSteps,progress=True)
                break
            except Exception as exceptionObject:
                #TODO: come up with exception summaries like this.
                if "finite" in str(exceptionObject): #This means there is an error message from zeus saying " Invalid walker initial positions!  Initialise walkers from positions of finite log probability."
                    print("One of the starting points has a non-finite probability. Picking new starting points. If you see this message like an infinite loop, consider trying the doEnsembleSliceSampling optional argument of walkerInitialDistributionSpread. It has a default value of 1.0. Reducing this value to 0.25, for example, may work if your initial guess is near the maximum of the posterior distribution.")
                    #Need to make the sampler again, in this case, to throw away anything that has happened so far
                    walkerStartPoints = self.generateInitialPoints(initialPointsDistributionType=walkerInitialDistribution, numStartPoints = self.mcmc_nwalkers, relativeInitialDistributionSpread=walkerInitialDistributionSpread) 
                    emcee_sampler = emcee.EnsembleSampler(self.mcmc_nwalkers, numParameters, log_prob_fn=self.getLogP)    
                elif "maxiter" in str(exceptionObject): #This means there is an error message from zeus that the max iterations have been reached.
                    print("WARNING: One or more of the Ensemble Slice Sampling walkers encountered an error. The value of mcmc_maxiter is currently", mcmc_maxiter, "you should increase it, perhaps by a factor of 1E2.")
                else:
                    print(str(exceptionObject))
                    sys.exit()
        #Now to keep the results:
        #TODO: when implementing convergence diagnostics, extract more chain information from get_chain without flattening (i think)
        adjusted_mcmc_burn_in_length = int(self.mcmc_burn_in_length / self.mcmc_nwalkers)
        self.post_burn_in_samples = emcee_sampler.get_chain(flat=True, discard = adjusted_mcmc_burn_in_length )
        discrete_chains_post_burn_in_samples = emcee_sampler.get_chain(flat=False, discard = adjusted_mcmc_burn_in_length)
        if self.UserInput.parameter_estimation_settings['mcmc_store_samplingObject'] == True:
            self.discrete_chains_post_burn_in_samples = discrete_chains_post_burn_in_samples
        self.post_burn_in_log_posteriors_un_normed_vec = np.atleast_2d(emcee_sampler.get_log_prob(flat=True, discard=adjusted_mcmc_burn_in_length)).transpose() #Needed to make it 2D and transpose.
        self.mcmc_last_point_sampled = discrete_chains_post_burn_in_samples[-1] #gets the last sample of each chain.
        # Populate the last_MCMC_run_type in UserInput so that there is a record of which MCMC run type was last used with this UserInput instance.
        self.UserInput.last_MCMC_run_type = 'EJS'
        if continueSampling == True:
            self.post_burn_in_samples = np.vstack((self.last_post_burn_in_samples, self.post_burn_in_samples ))
            self.post_burn_in_log_posteriors_un_normed_vec = np.vstack( (self.last_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))        
        #####BELOW HERE SHOUD BE SAME FOR doMetropolisHastings and doEnsembleSliceSampling and doEnsembleModifiedMHSampling#####
        if (self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or self.UserInput.parameter_estimation_settings['multistart_parallel_sampling']) == True: #If we're using certain parallel processing, we need to make calculatePostBurnInStatistics into True.
            calculatePostBurnInStatistics = True
        if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling']: #mcmc_exportLog == True is needed for mcmc_parallel_sampling, but not for multistart_parallel_sampling
            mcmc_exportLog=True
        if self.UserInput.parameter_estimation_settings['mcmc_store_samplingObject']:
            self.samplerType = 'EnsembleModifiedMHSampling'
            self.samplingObject = emcee_sampler
        if calculatePostBurnInStatistics == True:
            self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True) #This function call will also filter the lowest probability samples out, when using default settings.
            if self.UserInput.user_requested_convergence_diagnostics: #Run convergence diagnostics if UserInput defines it as True
                self.getConvergenceDiagnostics(discrete_chains_post_burn_in_samples)
            if str(mcmc_exportLog) == 'UserChoice':
                mcmc_exportLog = bool(self.UserInput.parameter_estimation_settings['mcmc_exportLog'])
            if mcmc_exportLog == True:
                self.exportPostBurnInStatistics()
            if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] == True: #We don't call the below function at this time unless we are doing mcmc_parallel_sampling. For multistart_parallel_sampling the consolidation is done elsewhere and differently.
                self.consolidate_parallel_sampling_data(parallelizationType="equal", mpi_cached_files_prefix='mcmc')
            return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]   
        else: #In this case, we are probably doing a PermutationSearch or something like that and only want self.map_logP.
            self.map_logP = max(self.post_burn_in_log_posteriors_un_normed_vec)
            self.map_index = list(self.post_burn_in_log_posteriors_un_normed_vec).index(self.map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
            self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.            
            return self.map_logP
    
    # pointer to Ensemble Modified MH sampling
    # doEnsembleJumpSampling is provided for user convenience as a simpler choice of sampling
    doEnsembleJumpSampling = doEnsembleModifiedMHSampling


    #Our EnsembleSliceSampling is done by the Zeus back end. (pip install zeus-mcmc)
    software_name = "zeus"
    software_version = "2.0.0"
    software_unique_id = "https://github.com/minaskar/zeus"
    software_kwargs = {"version": software_version, "author": ["Minas Karamanis", "Florian Beutler"], "cite": ["Minas Karamanis and Florian Beutler. zeus: A Python Implementation of the Ensemble Slice Sampling method. 2020. ","https://arxiv.org/abs/2002.06212", "@article{ess,  title={Ensemble Slice Sampling}, author={Minas Karamanis and Florian Beutler}, year={2020}, eprint={2002.06212}, archivePrefix={arXiv}, primaryClass={stat.ML} }"] }
    #@CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    @CiteSoft.function_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
    def doEnsembleSliceSampling(self, mcmc_nwalkers_direct_input = None, walkerInitialDistribution='UserChoice', walkerInitialDistributionSpread='UserChoice', movesType='UserChoice', calculatePostBurnInStatistics=True, mcmc_exportLog ='UserChoice', continueSampling='auto'):
        #The distribution of walkers intial points can be uniform or gaussian or identical. As of Oct 2020, default is uniform spread around the intial guess.
        #The mcmc_nwalkers_direct_input is really meant for PermutationSearch to override the other settings, though of course people could also use it directly.  
        #The walkerInitialDistributionSpread is in relative units (relative to standard deviations). In the case of a uniform inital distribution the default level of spread is actually across two standard deviations, so the walkerInitialDistributionSpread is relative to that (that is, a value of 2 would give 2*2 = 4 for the full spread in each direction from the initial guess).
        import zeus
        if walkerInitialDistribution == 'UserChoice':
            walkerInitialDistribution = self.UserInput.parameter_estimation_settings['mcmc_walkerInitialDistribution']
        if walkerInitialDistribution.lower() == 'auto':
            walkerInitialDistribution = 'sobol'
        if str(walkerInitialDistributionSpread) == 'UserChoice':
            walkerInitialDistributionSpread = self.UserInput.parameter_estimation_settings['mcmc_walkerInitialDistributionSpread']
        if str(walkerInitialDistributionSpread).lower() == 'auto':
            walkerInitialDistributionSpread = 0.866  #This choice is intended to be useful for uniform distribution priors, as described in the UserInput file.
        # First check whether the UserInput already has been used to run another MCMC run, and if so make sure it is the same algorithm type
        if hasattr(self.UserInput, 'last_MCMC_run_type'):
            if self.UserInput.last_MCMC_run_type == 'ESS':
                pass
            else:
                print("ERROR: Two different algorithms of MCMC runs were attempted to be run within the same python instance without reloading the UserInput. To run two MCMC runs of different algorithms in the same python instance is not supported, but can be accomplished by using imp.reload(UserInput) followed by imp.reload(PEUQSE) between changing MCMC algorithms, and which requires re-populating the UserInput with your choices. "); sys.exit()

        #Check if we need to continue sampling, and prepare for it if we need to.
        if continueSampling == 'auto':
            if ('mcmc_continueSampling' not in self.UserInput.parameter_estimation_settings) or self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] == 'auto': #check that UserInput does not overrule the auto.
                if hasattr(self, 'mcmc_last_point_sampled'): #if we have an existing mcmc_last_point_sampled in the object, we will assume more sampling is desired.
                    continueSampling = True
                else:
                    continueSampling = False
            else: continueSampling = self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] 
        if continueSampling == True:
            if hasattr(self, 'mcmc_last_point_sampled'): #If we are continuing from an old mcmc in this object.
                self.last_post_burn_in_log_posteriors_un_normed_vec = copy.deepcopy(self.post_burn_in_log_posteriors_un_normed_vec)
                self.last_post_burn_in_samples = copy.deepcopy(self.post_burn_in_samples)
            else: #Else we need to read from the file.                
                #First check if we are doing some kind of parallel sampling, because in that case we need to read from the file for our correct process rank. We put that info into the prefix and suffix.
                file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix()
                self.last_logP_and_parameter_samples_filename = file_name_prefix + "mcmc_logP_and_parameter_samples" + file_name_suffix
                self.last_logP_and_parameter_samples_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_logP_and_parameter_samples_filename)
                self.last_post_burn_in_log_posteriors_un_normed_vec =  np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,0]))  #First column is the logP
                if np.shape(self.last_post_burn_in_log_posteriors_un_normed_vec)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_log_posteriors_un_normed_vec = self.last_post_burn_in_log_posteriors_un_normed_vec.transpose()
                self.last_post_burn_in_samples =   np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,1:])) #later columns are the samples.
                if np.shape(self.last_post_burn_in_samples)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_samples = self.last_post_burn_in_samples.transpose()
                self.mcmc_last_point_sampled_filename = file_name_prefix + "mcmc_last_point_sampled" + file_name_suffix
                self.mcmc_last_point_sampled_data = unpickleAnObject(self.UserInput.directories['pickles']+self.mcmc_last_point_sampled_filename)
                self.mcmc_last_point_sampled = self.mcmc_last_point_sampled_data
                self.last_InputParameterInitialGuess_filename = file_name_prefix + "mcmc_initial_point_parameters" + file_name_suffix
                self.last_InputParameterInitialGuess_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_InputParameterInitialGuess_filename)
                self.UserInput.InputParameterInitialGuess = self.last_InputParameterInitialGuess_data #populating this because otherwise non-grid Multi-Start will get the wrong values exported. & Same for final plots.
        ####these variables need to be made part of UserInput####
        numParameters = len(self.UserInput.mu_prior) #This is the number of parameters.
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if type(self.UserInput.parameter_estimation_settings['mcmc_random_seed']) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        if type(mcmc_nwalkers_direct_input) == type(None): #This is the normal case.
            if 'mcmc_nwalkers' not in self.UserInput.parameter_estimation_settings: self.mcmc_nwalkers = 'auto'
            else: self.mcmc_nwalkers = self.UserInput.parameter_estimation_settings['mcmc_nwalkers']
            if type(self.mcmc_nwalkers) == type("string"): 
                if self.mcmc_nwalkers.lower() == "auto":
                    self.mcmc_nwalkers = numParameters*4 # according to zeus paper (https://doi.org/10.1093/mnras/stab2867), 4*D is the optimal number of walkers for general problems. Multimodal problems benefit from additional walkers.
                else: #else it is an integer, or a string meant to be an integer.
                    self.mcmc_nwalkers =  int(self.mcmc_nwalkers)
        else: #this is mainly for PermutationSearch which will (by default) use the minimum number of walkers per point.
            self.mcmc_nwalkers = int(mcmc_nwalkers_direct_input)
        if (self.mcmc_nwalkers%2) != 0: #Check that it's even. If not, add one walker.
            print("The EnsembleSliceSampling requires an even number of Walkers. Adding one Walker.")
            self.mcmc_nwalkers = self.mcmc_nwalkers + 1
        requested_mcmc_steps = self.UserInput.parameter_estimation_settings['mcmc_length']
        nEnsembleSteps = int(requested_mcmc_steps/self.mcmc_nwalkers) #We calculate the calculate number of the Ensemble Steps from the total sampling steps requested divided by self.mcmc_nwalkers.
        if nEnsembleSteps == 0:
            nEnsembleSteps = 1
        if str(self.UserInput.parameter_estimation_settings['mcmc_burn_in']).lower() == 'auto': self.mcmc_burn_in_length = int(requested_mcmc_steps*0.1) # burn in is relative to entire ensemble
        else: self.mcmc_burn_in_length = self.UserInput.parameter_estimation_settings['mcmc_burn_in']
        if 'mcmc_maxiter' not in self.UserInput.parameter_estimation_settings: mcmc_maxiter = 1E6 #The default from zeus is 1E4, but I have found that is not always sufficient.
        else: mcmc_maxiter = self.UserInput.parameter_estimation_settings['mcmc_maxiter']
        # make the move objects from the UserInput or argument specification
        if movesType == 'UserChoice': movesType = self.UserInput.parameter_estimation_settings['mcmc_movesType']
        if movesType.lower() == 'auto': movesTypeObject = zeus.moves.DifferentialMove()
        elif movesType.lower() == 'differential': movesTypeObject = zeus.moves.DifferentialMove()
        elif movesType.lower() == 'global': movesTypeObject = zeus.moves.GlobalMove(n_components=4) #the n_components (believe it means number of modes to expect), has been hard coded to 4 from the default of 5. This was required to have Example00c2 to run without crashing.
        elif movesType.lower() == 'gaussian': movesTypeObject = zeus.moves.GaussianMove()
        elif movesType.lower() == 'kde': movesTypeObject = zeus.moves.KDEMove()
        else: print("The move type must be one of the following: ['auto', 'differential', 'global', 'gaussian', 'kde']")
        ####end of user input variables####
        #now to do the mcmc
        # if the initial guess is a string, read its pickle file for starting points for walkers
        # if this read array is a single point, proceed as normal
        # if this read array is multiple points, make sure that it has the same shape as number of walkers
        # else, print error and exit
        loaded_initial_guess_flag = False # flag to switch allow normal generation of initial points to guided if the loaded points are a single point
        if isinstance(self.UserInput.model['InputParameterInitialGuess'], str): # if yes, then assume a pickle file is imported
            # try to read from current directory and pickle directory
            start_point_pkl_file_name = self.UserInput.model['InputParameterInitialGuess']
            if '.pkl' in start_point_pkl_file_name: # remove '.pkl' from string if it is there
                start_point_pkl_file_name = start_point_pkl_file_name.replace('.pkl', '')
            from os.path import exists
            # check if file exists in current working directory
            if exists(start_point_pkl_file_name + '.pkl'):
                walkerStartPoints = unpickleAnObject(start_point_pkl_file_name)
                loaded_initial_guess_flag = True # makes sure the pickled array is used for starting points
            # check if file exists in pickle directory
            elif exists(self.UserInput.directories['pickles'] + start_point_pkl_file_name + '.pkl'):
                walkerStartPoints = unpickleAnObject(self.UserInput.directories['pickles'] + start_point_pkl_file_name)
                loaded_initial_guess_flag = True # makes sure the pickled array is used for starting points
            else:
                print('The pickled object for initial guess points must exist in base directory or pickles directory.')
                sys.exit()
            # make sure the initial guess is a numpy array
            walkerStartPoints = np.array(walkerStartPoints)
            # check if start points have the same amount of walkers
            initial_guess_mcmc_nwalkers = walkerStartPoints.shape[0] # last points have shape (nwalkers, nparams)
            if len(walkerStartPoints.shape) == 1:
                loaded_initial_guess_flag = False # set this to false so initial points can be generated
            elif initial_guess_mcmc_nwalkers != self.mcmc_nwalkers:
                print(f'Initial guess walker number must be the same shape as mcmc_nwalkers: {self.mcmc_nwalkers} mcmc walkers but initial guess had {initial_guess_mcmc_nwalkers} walkers')
                sys.exit()
        if continueSampling == False:
            # check if initial guess is loaded
            if not loaded_initial_guess_flag:
                walkerStartPoints = self.generateInitialPoints(initialPointsDistributionType=walkerInitialDistribution, numStartPoints = self.mcmc_nwalkers,relativeInitialDistributionSpread=walkerInitialDistributionSpread) #making the first set of starting points.
        elif continueSampling == True:
            # make burn in samples 0 since burn in has already occurred
            self.mcmc_burn_in_length = 0
            # start points have the same amount of walkers
            previous_mcmc_nwalkers = self.mcmc_last_point_sampled.shape[0] # last points have shape (nwalkers, nparams)
            if previous_mcmc_nwalkers != self.mcmc_nwalkers:
                print(f'Setting walkers to be the same as the previous run for continue sampling: from {self.mcmc_nwalkers} to {previous_mcmc_nwalkers} walkers')
                self.mcmc_nwalkers = previous_mcmc_nwalkers
            # start walkers at last walker points
            walkerStartPoints = self.mcmc_last_point_sampled
            # use global move if continue sampling occurs and movesType is on auto
            if movesType.lower() == 'auto': movesTypeObject = zeus.moves.GlobalMove(n_components=4); print('Using Global Move in continue sampling')
            # walkerStartPoints = self.map_parameter_set #used to be self.mcmc_last_point_sampled. However, ESS works best when sampling near the peak (if there is a monomodoal HPD).
        zeus_sampler = zeus.EnsembleSampler(self.mcmc_nwalkers, numParameters, logprob_fn=self.getLogP, maxiter=mcmc_maxiter, moves=movesTypeObject) #maxiter=1E4 is the typical number, but we may want to increase it based on some UserInput variable.        
        for trialN in range(0,1000):#Todo: This number of this range is hardcoded but should probably be a user selection.
            try:
                zeus_sampler.run_mcmc(walkerStartPoints, nEnsembleSteps)
                break
            except Exception as exceptionObject:
                if "finite" in str(exceptionObject): #This means there is an error message from zeus saying " Invalid walker initial positions!  Initialise walkers from positions of finite log probability."
                    print("One of the starting points has a non-finite probability. Picking new starting points. If you see this message like an infinite loop, consider trying the doEnsembleSliceSampling optional argument of walkerInitialDistributionSpread. It has a default value of 1.0. Reducing this value to 0.25, for example, may work if your initial guess is near the maximum of the posterior distribution.")
                    #Need to make the sampler again, in this case, to throw away anything that has happened so far
                    walkerStartPoints = self.generateInitialPoints(initialPointsDistributionType=walkerInitialDistribution, numStartPoints = self.mcmc_nwalkers, relativeInitialDistributionSpread=walkerInitialDistributionSpread) 
                    zeus_sampler = zeus.EnsembleSampler(self.mcmc_nwalkers, numParameters, logprob_fn=self.getLogP, maxiter=mcmc_maxiter, moves=movesTypeObject) #maxiter=1E4 is the typical number, but we may want to increase it based on some UserInput variable.        
                elif "maxiter" in str(exceptionObject): #This means there is an error message from zeus that the max iterations have been reached.
                    print("WARNING: One or more of the Ensemble Slice Sampling walkers encountered an error. The value of mcmc_maxiter is currently", mcmc_maxiter, "you should increase it, perhaps by a factor of 1E2.")
                else:
                    print(str(exceptionObject))
                    sys.exit()
        #Now to keep the results:
        adjusted_mcmc_burn_in_length = int(self.mcmc_burn_in_length / self.mcmc_nwalkers)
        self.post_burn_in_samples = zeus_sampler.samples.flatten(discard = adjusted_mcmc_burn_in_length )
        discrete_chains_post_burn_in_samples = zeus_sampler.get_chain(discard = adjusted_mcmc_burn_in_length)
        if self.UserInput.parameter_estimation_settings['mcmc_store_samplingObject'] == True:
            self.discrete_chains_post_burn_in_samples = discrete_chains_post_burn_in_samples
        self.post_burn_in_log_posteriors_un_normed_vec = np.atleast_2d(zeus_sampler.samples.flatten_logprob(discard=adjusted_mcmc_burn_in_length)).transpose() #Needed to make it 2D and transpose.
        self.mcmc_last_point_sampled = discrete_chains_post_burn_in_samples[-1] #Note that for **zeus** the last point sampled is actually an array of points equal to the number of walkers.        
        # Populate the last_MCMC_run_type in UserInput so that there is a record of which MCMC run type was last used with this UserInput instance.
        self.UserInput.last_MCMC_run_type = 'ESS'
        if continueSampling == True:
            self.post_burn_in_samples = np.vstack((self.last_post_burn_in_samples, self.post_burn_in_samples ))
            self.post_burn_in_log_posteriors_un_normed_vec = np.vstack( (self.last_post_burn_in_log_posteriors_un_normed_vec, self.post_burn_in_log_posteriors_un_normed_vec))        
        #####BELOW HERE SHOUD BE SAME FOR doMetropolisHastings and doEnsembleSliceSampling and doEnsembleModifiedMHSampling#####
        if (self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or self.UserInput.parameter_estimation_settings['multistart_parallel_sampling']) == True: #If we're using certain parallel processing, we need to make calculatePostBurnInStatistics into True.
            calculatePostBurnInStatistics = True;
        if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling']: #mcmc_exportLog == True is needed for mcmc_parallel_sampling, but not for multistart_parallel_sampling
            mcmc_exportLog=True
        if self.UserInput.parameter_estimation_settings['mcmc_store_samplingObject']:
            self.samplerType = 'EnsembleSliceSampling'
            self.samplingObject = zeus_sampler
        if calculatePostBurnInStatistics == True:
            self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True) #This function call will also filter the lowest probability samples out, when using default settings.
            if self.UserInput.user_requested_convergence_diagnostics: #Run convergence diagnostics if UserInput defines it as True
                self.getConvergenceDiagnostics(discrete_chains_post_burn_in_samples)
            if str(mcmc_exportLog) == 'UserChoice':
                mcmc_exportLog = bool(self.UserInput.parameter_estimation_settings['mcmc_exportLog'])
            if mcmc_exportLog == True:
                self.exportPostBurnInStatistics()
            if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] == True: #We don't call the below function at this time unless we are doing mcmc_parallel_sampling. For multistart_parallel_sampling the consolidation is done elsewhere and differently.
                self.consolidate_parallel_sampling_data(parallelizationType="equal", mpi_cached_files_prefix='mcmc')
            return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]   
        else: #In this case, we are probably doing a PermutationSearch or something like that and only want self.map_logP.
            self.map_logP = max(self.post_burn_in_log_posteriors_un_normed_vec)
            self.map_index = list(self.post_burn_in_log_posteriors_un_normed_vec).index(self.map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
            self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.            
            return self.map_logP
            
    def doEnsembleSliceSamplingGM(self):
        """A wrapper for Enseble Slice Sampler that uses Global Move"""
        self.doEnsembleSliceSampling(movesType='global')

    def doEnsembleSliceSamplingGlobal(self):
        """A wrapper for Ensemble Slice Sampler that uses Global Move"""
        self.doEnsembleSliceSampling(movesType='global')

    def doEnsembleSliceSamplingDifferential(self):
        """A wrapper for Ensemble Slice Sampler that uses Differential Move"""
        self.doEnsembleSliceSampling(movesType='differential')

    #main function to get samples #TODO: Maybe Should return map_log_P and mu_AP_log_P?
    #@CiteSoft.after_call_compile_consolidated_log() #This is from the CiteSoft module.
    def doMetropolisHastings(self, calculatePostBurnInStatistics = True, mcmc_exportLog='UserChoice', continueSampling = 'auto'):
        # First check whether the UserInput already has been used to run another MCMC run, and if so make sure it is the same algorithm type
        if hasattr(self.UserInput, 'last_MCMC_run_type'):
            if self.UserInput.last_MCMC_run_type == 'MH':
                pass
            else:
                print("ERROR: Two different algorithms of MCMC runs were attempted to be run within the same python instance without reloading the UserInput. To run two MCMC runs of different algorithms in the same python instance is not supported, but can be accomplished by using imp.reload(UserInput) followed by imp.reload(PEUQSE) between changing MCMC algorithms, and which requires re-populating the UserInput with your choices. "); sys.exit()

        #Check if we need to continue sampling, and prepare for it if we need to.
        if continueSampling == 'auto':
            if ('mcmc_continueSampling' not in self.UserInput.parameter_estimation_settings) or self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] == 'auto': #check that UserInput does not overrule the auto.
                if hasattr(self, 'mcmc_last_point_sampled'): #if we have an existing mcmc_last_point_sampled in the object, we will assume more sampling is desired.
                    continueSampling = True
                else:
                    continueSampling = False
            else: continueSampling = self.UserInput.parameter_estimation_settings['mcmc_continueSampling'] 
        if continueSampling == True:
            if hasattr(self, 'mcmc_last_point_sampled'): #if If we are continuing from an old
                self.last_post_burn_in_log_posteriors_un_normed_vec = copy.deepcopy(self.post_burn_in_log_posteriors_un_normed_vec)
                self.last_post_burn_in_samples = copy.deepcopy(self.post_burn_in_samples)
            else: #Else we need to read from the file.                
                #First check if we are doing some kind of parallel sampling, because in that case we need to read from the file for our correct process rank. We put that info into the prefix and suffix.
                file_name_prefix, file_name_suffix, directory_name_suffix = self.getParallelProcessingPrefixAndSuffix()
                self.last_logP_and_parameter_samples_filename = file_name_prefix + "mcmc_logP_and_parameter_samples" + file_name_suffix
                self.last_logP_and_parameter_samples_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_logP_and_parameter_samples_filename)
                self.last_post_burn_in_log_posteriors_un_normed_vec =  np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,0]))  #First column is the logP
                if np.shape(self.last_post_burn_in_log_posteriors_un_normed_vec)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_log_posteriors_un_normed_vec = self.last_post_burn_in_log_posteriors_un_normed_vec.transpose()
                self.last_post_burn_in_samples =   np.array(nestedObjectsFunctions.makeAtLeast_2dNested(self.last_logP_and_parameter_samples_data[:,1:])) #later columns are the samples.
                if np.shape(self.last_post_burn_in_samples)[0] == 1: #In this case, need to transpose.
                    self.last_post_burn_in_samples = self.last_post_burn_in_samples.transpose()
                self.mcmc_last_point_sampled = self.last_post_burn_in_samples[-1]        
                self.last_InputParameterInitialGuess_filename = file_name_prefix + "mcmc_initial_point_parameters" + file_name_suffix
                self.last_InputParameterInitialGuess_data = unpickleAnObject(self.UserInput.directories['pickles']+self.last_InputParameterInitialGuess_filename)
                self.UserInput.InputParameterInitialGuess = self.last_InputParameterInitialGuess_data #populating this because otherwise non-grid Multi-Start will get the wrong values exported. & Same for final plots.
        #Setting burn_in_length in below few lines (including case for continued sampling).
        if str(self.UserInput.parameter_estimation_settings['mcmc_burn_in']).lower() == 'auto': self.mcmc_burn_in_length = int(self.UserInput.parameter_estimation_settings['mcmc_length']*0.1)
        else: self.mcmc_burn_in_length = self.UserInput.parameter_estimation_settings['mcmc_burn_in']
        if continueSampling == True:
            self.mcmc_burn_in_length = 0
        if 'mcmc_random_seed' in self.UserInput.parameter_estimation_settings:
            if type(self.UserInput.parameter_estimation_settings['mcmc_random_seed']) == type(1): #if it's an integer, then it's not a "None" type or string, and we will use it.
                np.random.seed(self.UserInput.parameter_estimation_settings['mcmc_random_seed'])
        samples_simulatedOutputs = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],self.UserInput.num_data_points))
        samples = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],len(self.UserInput.mu_prior)))
        mcmc_step_modulation_history = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: Make this optional for efficiency. #This allows the steps to be larger or smaller. Make this same length as samples. In future, should probably be same in other dimension also, but that would require 2D sampling with each step.                                                                          
        if continueSampling == False:
            samples[0,:]=self.UserInput.InputParameterInitialGuess  # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean if an initial guess is not provided.. Do not use self.UserInput.model['InputParameterInitialGuess']  because that doesn't work with reduced parameter space feature.
        elif continueSampling == True:
            samples[0,:]= self.mcmc_last_point_sampled
        samples_drawn = samples*1.0 #this includes points that were rejected. #TODO: make this optional for efficiency.               
        log_likelihoods_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        log_posteriors_un_normed_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        log_postereriors_drawn = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'])) #TODO: make this optional for efficiency. We don't want this to be 2D, so we don't copy log_posteriors_un_normed_vec.
        log_priors_vec = np.zeros((self.UserInput.parameter_estimation_settings['mcmc_length'],1))
        #Code to initialize checkpoints.
        if type(self.UserInput.parameter_estimation_settings['mcmc_checkPointFrequency']) != type(None):
            print("Starting MCMC sampling.")
            timeOfFirstCheckpoint = time.time()
            timeCheckpoint = time.time() - timeOfFirstCheckpoint #First checkpoint at time 0.
            numCheckPoints = self.UserInput.parameter_estimation_settings['mcmc_length']/self.UserInput.parameter_estimation_settings['mcmc_checkPointFrequency']
        #Before sampling should fill in the first entry for the posterior vector we have created. #TODO: It would probably be better to start with i of 0 in below sampling loop. I believe that right now the "burn in" and "samples" arrays are actually off by an index of 1. But trying to change that alters their length relative to other arrays and causes problems. Since we always do many samples and this only affects the initial point being averaged in twice, it is not a major problem. It's also avoided if people use a burn in of at least 1.
        log_posteriors_un_normed_vec[0]= self.getLogP(samples[0])
        for i in range(1, self.UserInput.parameter_estimation_settings['mcmc_length']): #TODO: This should be checked... Don't we need to start with i of 0?
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
            [log_likelihood_current_location, simulationOutput_current_location] = self.getLogLikelihood(samples[i-1,:]) #TODO: the previous likelihood should be stored so that it doesn't need to be calculated again.
            #TODO: The below line (log_accept_probability calculation) should be based on the current and earlier logP (not adding the likelihoods and priors), that way we can rely on the getLogP function and can just store the current and previous logP.
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
                log_postereriors_drawn[i] = (log_likelihood_proposal+log_prior_proposal) #TODO: should be using getlogP for the proposed and prior location. log_accept_probability line needs to be changed.
                #samples_simulatedOutputs[i,:] = nestedObjectsFunctions.flatten_2dNested(simulationOutput_proposal)
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
                #samples_simulatedOutputs[i,:] = nestedObjectsFunctions.flatten_2dNested(simulationOutput_current_location)
                log_posteriors_un_normed_vec[i] = log_likelihood_current_location+log_prior_current_location
                log_likelihoods_vec[i] = log_likelihood_current_location
                log_priors_vec[i] = log_prior_current_location
            if type(self.UserInput.parameter_estimation_settings['mcmc_checkPointFrequency']) != type(None):
                if sampleNumber%self.UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] == 0: #The % is a modulus function.
                    timeSinceLastCheckPoint = (time.time() - timeOfFirstCheckpoint) -  timeCheckpoint
                    timeCheckpoint = time.time() - timeOfFirstCheckpoint
                    checkPointNumber = sampleNumber/self.UserInput.parameter_estimation_settings['mcmc_checkPointFrequency']
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
#                        np.savetxt(self.UserInput.directories['logs_and_csvs']+"recent_stacked.csv",recent_stacked, delimiter=',')
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
        if continueSampling == False: #Normally we enter this if statement and collect the burn_in samples. If continueSampling, either they already exist in the PE_object and will be exported again later, or they don't exist in the PE_object because are continuing from a previous python instance. If continuing from a previous python instance, during_burn_in_samples won't be created and also won't be exported again. (This logic is so we don't overwrite the old during_burn_in_samples).
            self.during_burn_in_samples = samples[0:self.mcmc_burn_in_length] 
            self.during_burn_in_log_posteriors_un_normed_vec = log_posteriors_un_normed_vec[0:self.mcmc_burn_in_length]
        self.post_burn_in_samples = samples[self.mcmc_burn_in_length:] 
        discrete_chains_post_burn_in_samples = np.expand_dims(self.post_burn_in_samples, axis=1) # MH only has one chain
        if self.UserInput.parameter_estimation_settings['mcmc_store_samplingObject'] == True:
            self.discrete_chains_post_burn_in_samples = discrete_chains_post_burn_in_samples
        #self.post_burn_in_samples_simulatedOutputs = copy.deepcopy(samples_simulatedOutputs)
        #self.post_burn_in_samples_simulatedOutputs[self.mcmc_burn_in_length:0] #Note: this feature is presently not compatible with continueSampling.
        self.post_burn_in_log_posteriors_un_normed_vec = log_posteriors_un_normed_vec[self.mcmc_burn_in_length:]
        self.mcmc_last_point_sampled = self.post_burn_in_samples[-1]
        self.post_burn_in_log_likelihoods_vec = log_likelihoods_vec[self.mcmc_burn_in_length:]
        self.post_burn_in_log_priors_vec = log_priors_vec[self.mcmc_burn_in_length:]        
        # Populate the last_MCMC_run_type in UserInput so that there is a record of which MCMC run type was last used with this UserInput instance.
        self.UserInput.last_MCMC_run_type = 'MH'
        #####BELOW HERE SHOUD BE SAME FOR doMetropolisHastings and doEnsembleSliceSampling and doEnsembleModifiedMHSampling#####
        if continueSampling == True:
            self.post_burn_in_samples = np.vstack((self.last_post_burn_in_samples, self.post_burn_in_samples ))
            self.post_burn_in_log_posteriors_un_normed_vec = np.vstack( (np.array(self.last_post_burn_in_log_posteriors_un_normed_vec), np.array(self.post_burn_in_log_posteriors_un_normed_vec)))
        if (self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] or self.UserInput.parameter_estimation_settings['multistart_parallel_sampling']) == True: #If we're using certain parallel processing, we need to make calculatePostBurnInStatistics into True.
            calculatePostBurnInStatistics = True;
        if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling']: #mcmc_exportLog == True is needed for mcmc_parallel_sampling, but not for multistart_parallel_sampling
            mcmc_exportLog=True
        if calculatePostBurnInStatistics == True:
            #TODO: I think Below, calculate_post_burn_in_log_priors_vec=True should be false unless we are using continue sampling. For now, will leave it since I am not sure why it is currently set to True/False.
            self.calculatePostBurnInStatistics(calculate_post_burn_in_log_priors_vec = True) #This function call will also filter the lowest probability samples out, when using default settings.
            if self.UserInput.user_requested_convergence_diagnostics: #Run convergence diagnostics if UserInput defines it as True
                self.getConvergenceDiagnostics(discrete_chains_post_burn_in_samples)
            if str(mcmc_exportLog) == 'UserChoice':
                mcmc_exportLog = bool(self.UserInput.parameter_estimation_settings['mcmc_exportLog'])
            if mcmc_exportLog == True:
                self.exportPostBurnInStatistics()
            if self.UserInput.parameter_estimation_settings['mcmc_parallel_sampling'] == True: #We don't call the below function at this time unless we are doing mcmc_parallel_sampling. For multistart_parallel_sampling the consolidation is done elsewhere and differently.
                self.consolidate_parallel_sampling_data(parallelizationType="equal", mpi_cached_files_prefix='mcmc')
            return [self.map_parameter_set, self.mu_AP_parameter_set, self.stdap_parameter_set, self.evidence, self.info_gain, self.post_burn_in_samples, self.post_burn_in_log_posteriors_un_normed_vec]   
        else: #In this case, we are probably doing a PermutationSearch or something like that and only want self.map_logP.
            self.map_logP = max(self.post_burn_in_log_posteriors_un_normed_vec)
            self.map_index = list(self.post_burn_in_log_posteriors_un_normed_vec).index(self.map_logP) #This does not have to be a unique answer, just one of them places which gives map_logP.
            self.map_parameter_set = self.post_burn_in_samples[self.map_index] #This  is the point with the highest probability in the posterior.            
            return self.map_logP
        
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
                print("WARNING: There is an error in your self.UserInput.scaling_uncertainties. This probably means that your uncertainties array does not have a size matching the number of parameters expected by your simulation function. If this is not the situation, contact the developers with a bug report. Send your input file and simulation function file.")
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
        if len(self.UserInput.InputParameterPriorValues_upperBounds) > 0:
            upperCheck = boundsCheck(discreteParameterVector, self.UserInput.InputParameterPriorValues_upperBounds, 'upper')
            if upperCheck == False:
                return False
        if len(self.UserInput.InputParameterPriorValues_lowerBounds) > 0:
            lowerCheck = boundsCheck(discreteParameterVector, self.UserInput.InputParameterPriorValues_lowerBounds, 'lower')
            if lowerCheck == False:
                return False
        return True #If the test has gotten here without failing any of the tests, we return true.

    #The responses bounds check is different from the parmeter bounds check because responses are nested.
    def doSimulatedResponsesBoundsChecks(self, simulatedResponses): #Bounds intended for the likelihood.
        if len(self.UserInput.InputParameterPriorValues_upperBounds) > 0:
            for responseIndex in self.UserInput.model['simulatedResponses_upperBounds']:
                upperCheck = boundsCheck(simulatedResponses[responseIndex], self.UserInput.model['simulatedResponses_upperBounds'][responseIndex], 'upper')
                if upperCheck == False:
                    return False
        if len(self.UserInput.InputParameterPriorValues_lowerBounds) > 0:
            for responseIndex in self.UserInput.model['simulatedResponses_lowerBounds']:
                lowerCheck = boundsCheck(simulatedResponses[responseIndex], self.UserInput.model['simulatedResponses_lowerBounds'][responseIndex], 'lower')
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
        try:#This warning will not always work if there are multiple responses. #TODO: make this a loop across the number of responses. For now, just making it a "try" and "except" statement.
            if np.array(simulationOutput).any()==float('nan'):
                print("WARNING: Your simulation output returned a 'nan' for parameter values " +str(discreteParameterVector) + ". 'nan' values cannot be processed by the PEUQSE software and this set of Parameter Values is being assigned a probability of 0.")
                return None #This is intended for the case that the simulation fails in some way without returning "None". 
        except:
            pass
        if type(simulationOutputProcessingFunction) == type(None):
            simulatedResponses = simulationOutput 
        elif type(simulationOutputProcessingFunction) != type(None):
            simulatedResponses = simulationOutputProcessingFunction(simulationOutput) 
        simulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(simulatedResponses)
        if self.doSimulatedResponsesBoundsChecks(simulatedResponses) == False:
            simulatedResponses = None
        #if self.UserInput.parameter_estimation_settings['exportAllSimulatedOutputs' == True: 
        #decided to always keep the lastSimulatedResponses in memory. Should be okay because only the most recent should be kept.
        #At least, that is my understanding after searching for "garbage" here and then reading: http://www.digi.com/wiki/developer/index.php/Python_Garbage_Collection
        self.lastSimulatedResponses = copy.deepcopy(simulatedResponses)
        return simulatedResponses
    
    def getLogLikelihood(self,discreteParameterVector): #The variable discreteParameterVector represents a vector of values for the parameters being sampled. So it represents a single point in the multidimensional parameter space.
        discreteParameterVectorTuple = tuple(np.ndarray.copy(np.array(discreteParameterVector))) #we use a tuple in case any of the functions try to change the parameters.
        
        #Check if user has provided a custom log likelihood function.
        if type(self.UserInput.model['custom_logLikelihood']) != type(None):
            logLikelihood, simulatedResponses = self.UserInput.model['custom_logLikelihood'](discreteParameterVectorTuple)
            simulatedResponses = np.array(simulatedResponses).flatten()
            return logLikelihood, simulatedResponses
        #else pass is implied.

        #Now get the simulated responses.
        simulatedResponses = self.getSimulatedResponses(discreteParameterVectorTuple)
        #Failure checks:
        if type(simulatedResponses) == type(None):
            print("Warning: There are likelihood points that have zero probability due to receiving a None type back for the simulated responses. If there are too many points like this during an MCMC or permutationsToSamples run, the MAP and mu_AP returned will not be meaningful. Parameters: " + str(discreteParameterVectorTuple))
            return float('-inf'), None #This is intended for the case that the simulation fails, indicated by receiving an 'nan' or None type from user's simulation function.
        #Check if there are any 'nan' in the simulations array, and treat that as a failure also.
        nans_in_array = np.isnan(simulatedResponses)        
        if True in nans_in_array:
            print("Warning: There are likelihood points that have zero probability due to not a number ('nan') values in the simulated responses. If there are too many points like this, the MAP and mu_AP returned will not be meaningful. Parameters: " + str(discreteParameterVectorTuple))
            return float('-inf'), None
        #need to check if there are any 'responses_simulation_uncertainties'.
        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #if it's a None type, we keep it as a None type
            responses_simulation_uncertainties = None
        else:  #Else we get it based on the the discreteParameterVectorTuple
            responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(discreteParameterVectorTuple)
        #Now need to do transforms. Transforms are only for calculating log likelihood. If responses_simulation_uncertainties is "None", then we need to have one less argument passed in and a blank list is returned along with the transformed simulated responses.
        if type(responses_simulation_uncertainties) == type(None):
            simulatedResponses_transformed, blank_list = self.transform_responses(simulatedResponses) #This creates transforms for any data that we might need it. The same transforms were also applied to the observed responses.
            responses_simulation_uncertainties_transformed = None
            simulated_responses_covmat_transformed = None
        else:
            simulatedResponses_transformed, responses_simulation_uncertainties_transformed = self.transform_responses(simulatedResponses, responses_simulation_uncertainties) #This creates transforms for any data that we might need it. The same transforms were also applied to the observed responses.
            simulated_responses_covmat_transformed = returnShapedResponseCovMat(self.UserInput.num_response_dimensions, responses_simulation_uncertainties_transformed)  #assume we got standard deviations back.
        log_probability_metric, simulatedResponses_transformed = self.getLogLikelihood_byResponses(simulatedResponses_transformed, simulated_responses_covmat_transformed)
        return log_probability_metric, simulatedResponses_transformed

    def getLogLikelihood_byResponses(self, simulatedResponses_transformed, simulated_responses_covmat_transformed=None, observedResponses_transformed=None,observed_responses_covmat_transformed=None):      
        #This function is meant for internal use only, but could be called by a user if desired.
        #simulatedResponses_transformed and observedResponses_transformed are just the simulatedResponses and observedResponses in their "final" form: lack of transform is completely normal. If a user has a need to call this function directly and is not planning on transforming their responses (such as log transform or integral) or if the user does not know what transform to use, then the user should simply provide simulatedResponses and observedResponses as the "already transformed" data.
        
        #This function expects the uncertainties received to be in covmat form already (1 covmat per response), but can instead be a vectors of variances per response (NOT standard deviations: if the user is calling this function directly and only has standard deviations, then the standard deviations must be squared before calling this funciton,such as by np.square).  The simulatedResponses_transformed must be the same shape as observedResponses_transformed.
        #The simulated_responses_covmat_transformed and observed_responses_covmat_transformed must be teh same shape is both are provided. 
        
        if type(observedResponses_transformed) == type(None):
            observedResponses_transformed = self.UserInput.responses_observed_transformed
        #If our likelihood is  “probability of Response given Theta”…  we have a continuous probability distribution for both the response and theta. That means the pdf  must use binning on both variables. Eric notes that the pdf returns a probability density, not a probability mass. So the pdf function here divides by the width of whatever small bin is being used and then returns the density accordingly. Because of this, our what we are calling likelihood is not actually probability (it’s not the actual likelihood) but is proportional to the likelihood.
        #Thus we call it a probability_metric and not a probability. #TODO: consider changing names of likelihood and get likelihood to "likelihoodMetric" and "getLikelihoodMetric"
        #Now we need to make the comprehensive_responses_covmat.
        #First we will check whether observed_responses_covmat_transformed is square or not. The multivariate_normal.pdf function requires a diagonal values vector to be 1D.
        if type(observed_responses_covmat_transformed) == type(None):
            observed_responses_covmat_transformed = self.observed_responses_covmat_transformed
        
        #For some cases, we should prepare to split the likelihood.
        if self.prepareResponsesForSplitLikelihood == True:
            simulatedResponses_transformed = nestedObjectsFunctions.convertToNested(simulatedResponses_transformed[0])
            if type(simulated_responses_covmat_transformed) != type(None):
                    simulated_responses_covmat_transformed = nestedObjectsFunctions.convertToNested(simulated_responses_covmat_transformed[0])
            observedResponses_transformed = nestedObjectsFunctions.convertToNested(observedResponses_transformed[0])
            observed_responses_covmat_transformed = nestedObjectsFunctions.convertToNested(observed_responses_covmat_transformed[0])
        #In general, the covmat could be a function of the responses magnitude and independent variables. So eventually, we will use non-linear regression or something to estimate it. However, for now we simply take the observed_responses_covmat_transformed which will work for most cases.
        #TODO: use Ashi's nonlinear regression code (which  he used in this paper https://www.sciencedirect.com/science/article/abs/pii/S0920586118310344).  Put in the response magnitudes and the independent variables.
        #in future it will be something like: if self.UserInput.covmat_regression== True: comprehensive_responses_covmat = nonLinearCovmatPrediction(self.UserInput['independent_variable_values'], observed_responses_covmat_transformed)
        #And that covmat_regression will be on by default.  We will need to have an additional argument for people to specify whether magnitude weighting and independent variable values should both be considered, or just one.
        #First, get the shape of the covmat.
        observed_responses_covmat_transformed_shape = np.shape(observed_responses_covmat_transformed)
        if type(simulated_responses_covmat_transformed) == type(None):
            comprehensive_responses_covmat = observed_responses_covmat_transformed
        else: #Else we add the uncertainties, assuming they are orthogonal. Note that these are already covmats so are already variances that can be added directly. 
            comprehensive_responses_covmat = observed_responses_covmat_transformed + simulated_responses_covmat_transformed #TODO: I think think this needs to be moved own into the responseIndex loop to correctly handle staggered uncertainties. [like one response having full covmatrix and others not]
        comprehensive_responses_covmat_shape = copy.deepcopy(observed_responses_covmat_transformed_shape) #no need to take the shape of the actual comprehensive_responses_covmat since they must be same. This is probably slightly less computation.
        if (len(comprehensive_responses_covmat_shape) == 1) and (comprehensive_responses_covmat_shape[0]==1): #Matrix is square because has only one value.
            log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed,x=observedResponses_transformed,cov=comprehensive_responses_covmat)
            if float(log_probability_metric) == float('-inf'):
                print("Warning: There are likelihood points that have zero probability. If there are too many points like this, the MAP and mu_AP returned will not be meaningful. Parameters: " + str(discreteParameterVectorTuple))
            return log_probability_metric, simulatedResponses_transformed #Return this rather than going through loop further.
        elif len(comprehensive_responses_covmat_shape) > 1 and (comprehensive_responses_covmat_shape[0] == comprehensive_responses_covmat_shape[1]):  #Else it is 2D, check if it's square.
            try:
                log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed,x=observedResponses_transformed,cov=comprehensive_responses_covmat)
                if float(log_probability_metric) == float('-inf'):
                    print("Warning: There are likelihood points that have zero probability. If there are too many points like this, the MAP and mu_AP returned will not be meaningful. Parameters: " + str(discreteParameterVectorTuple))
                return log_probability_metric, simulatedResponses_transformed #Return this rather than going through loop further.
            except:
                pass #If it failed, we assume it is not square. For example, it could be 2 responses of length 2 each, which is not actually square.
            #TODO: Put in near-diagonal solution described in github: https://github.com/AdityaSavara/PEUQSE/issues/3
        #If neither of the above return statements have occurred, we should go through the uncertainties per response.
        log_probability_metric = 0 #Initializing since we will be adding to it.
        for responseIndex in range(self.UserInput.num_response_dimensions):
            #We will check if the response has too many values. If has too many values, then the covmat will be too large and will evaluate each value separately (with only variance, no covariance) in order to achive a linear scaling.
            if len(simulatedResponses_transformed[responseIndex]) > self.UserInput.responses['responses_observed_max_covmat_size']:
                calculate_log_probability_metric_per_value = True
                response_log_probability_metric = 0 #initializing so that can check if it is a 'nan' or not a bit further down below.
            else:
                calculate_log_probability_metric_per_value = False
                #no need oto intialize response_log_probability_metric.
            #Now try to calculate response_log_probability_metric
            if calculate_log_probability_metric_per_value == False: #The normal case.            
                try: #try to evaluate, but switch to individual values if there is any problem.
                    response_log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed[responseIndex],x=observedResponses_transformed[responseIndex],cov=comprehensive_responses_covmat[responseIndex])  #comprehensive_responses_covmat has to be 2D or has to be 1D array/list of variances of length equal to x.
                except:
                    response_log_probability_metric = float('nan') #this keeps track of failure cases.
                    calculate_log_probability_metric_per_value = False
            if calculate_log_probability_metric_per_value == True:
                if response_log_probability_metric == float('nan'): # if a case failed...
                    response_log_probability_metric = -1E100 #Just initializing, then will add each probability separately. One for each **value** of this response dimension. The -1E100 is to penalize any cases responses that failed.
                else: 
                    response_log_probability_metric = 0 #No penalty if the 'per value' calculation is being done for non-failure reasons, like the number of values being too long to use a covmat directly.
                for responseValueIndex in range(len(simulatedResponses_transformed[responseIndex])):
                    try:
                        current_log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed[responseIndex][responseValueIndex],x=observedResponses_transformed[responseIndex][responseValueIndex],cov=comprehensive_responses_covmat[responseIndex][responseValueIndex])    
                    except: #The above is to catch cases when the multivariate_normal fails.
                        current_log_probability_metric = float('-inf')
                    #response_log_probability_metric = current_log_probability_metric + response_log_probability_metric
                    if float(current_log_probability_metric) == float('-inf'):
                        warnings.warn("Warning: There are cases of sampling where a response value has zero probability in the likelihood. If there are too many points like this, the MAP and mu_AP returned will not be meaningful. ResponseIndex: " + str( responseIndex))
                        current_log_probability_metric = -1E100 #Just choosing an arbitrarily very severe penalty. I know that I have seen 1E-48 to -303 from the multivariate pdf, and values inbetween like -171, -217, -272. I found that -1000 seems to be worse, but I don't have a systematic testing. I think -1000 was causing numerical errors.
                    response_log_probability_metric = current_log_probability_metric + response_log_probability_metric
            log_probability_metric = log_probability_metric + response_log_probability_metric
        return log_probability_metric, simulatedResponses_transformed


    def truncatePostBurnInSamples(self, post_burn_in_samples=[], parameterBoundsLower=None, parameterBoundsUpper=None):
        """
        Truncate the post_burn_in_samples variable along with other variables that are relative to it. 
        Apply a lower and upper bound to a parameter to truncate.
        Bounds are Inclusive. Set a None on a bound to not have it change. 
        Example of changing just lower bound for one parameter: parameterBounds = [(paramLowerBound, None)]

        :param post_burn_in_samples: Samples after mcmc run. (:type: np.array)
        :param parameterBoundsLower: List of parameter lower bounds. Use None to indicate a bound not being applied. (:type: list)
        :param parameterBoundsUpper: List of parameter upper bounds. Use None to indicate a bound not being applied. (:type: list)
        """
        if post_burn_in_samples == []:
            post_burn_in_samples = self.post_burn_in_samples
        # truncate the post burn in samples according to the parameterBounds
        truncated_post_burn_in_samples, truncated_mask = truncateSamples(post_burn_in_samples, parameterBoundsLower=parameterBoundsLower, parameterBoundsUpper=parameterBoundsUpper, returnMask=True)
        # reassign class variables to the truncated versions
        self.post_burn_in_samples = truncated_post_burn_in_samples
        self.post_burn_in_log_posteriors_un_normed_vec = self.post_burn_in_log_posteriors_un_normed_vec[truncated_mask, :]
        self.post_burn_in_log_priors_vec = self.post_burn_in_log_priors_vec[truncated_mask, :]
        # return post_burn_in_samples
        return truncated_post_burn_in_samples

    def makeHistogramsForEachParameter(self, showFigure=None):
        if showFigure == None: showFigure = True
        import PEUQSE.plotting_functions as plotting_functions 
        setMatPlotLibAgg(self.UserInput.plotting_ouput_settings['setMatPlotLibAgg'])
        parameterSamples = self.post_burn_in_samples
        parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        if hasattr(self.UserInput, 'histogram_plot_settings') == False: #put some defaults for backwards compatibility.
            self.UserInput.histogram_plot_settings={}
        plotting_functions.makeHistogramsForEachParameter(parameterSamples,parameterNamesAndMathTypeExpressionsDict, directory = self.UserInput.directories['graphs'], parameterInitialValue=self.UserInput.model['InputParameterPriorValues'], parameterMAPValue=self.map_parameter_set, parameterMuAPValue=self.mu_AP_parameter_set, histogram_plot_settings=self.UserInput.histogram_plot_settings, showFigure=showFigure)

    def makeSamplingScatterMatrixPlot(self, parameterSamples = [], parameterNamesAndMathTypeExpressionsDict={}, parameterNamesList =[], parameterMAPValue=[], parameterMuAPValue=[], parameterInitialValue = [], showFigure=None, plot_settings={'combined_plots':'auto'}):
        #showFigure will be set down below depending on whether combined plots is occurring or not.
        import pandas as pd #This is one of the only functions that use pandas.
        import matplotlib.pyplot as plt
        import PEUQSE.plotting_functions as plotting_functions
        if 'dpi' not in  plot_settings:  plot_settings['dpi'] = 220
        if 'figure_name' not in  plot_settings:  plot_settings['figure_name'] = 'scatter_matrix_posterior'
        if parameterSamples  ==[] : parameterSamples = self.post_burn_in_samples
        if parameterNamesAndMathTypeExpressionsDict == {}: parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        if parameterNamesList == []: parameterNamesList = self.UserInput.parameterNamesList #This is created when the parameter_estimation object is initialized.        
        if parameterMAPValue == []: parameterMAPValue = self.map_parameter_set
        if parameterMuAPValue == []: parameterMuAPValue = self.mu_AP_parameter_set
        if parameterInitialValue == []: parameterInitialValue = self.UserInput.model['InputParameterPriorValues']
        combined_plots = plot_settings['combined_plots']
        if 'individual_plots' in plot_settings: 
            individual_plots = plot_settings['individual_plots']
        else: #TODO: A. Savara added this code to prevent crashes for cases where individual_plots was not defined. However, it should be defined automatically before this point.
            individual_plots = True
        posterior_df = pd.DataFrame(parameterSamples,columns=[parameterNamesAndMathTypeExpressionsDict[x] for x in parameterNamesList])
        if combined_plots == 'auto': #by default, we will not make the scatter matrix when there are more than 5 parameters.
            if (len(parameterNamesList) > 5) or (len(parameterNamesAndMathTypeExpressionsDict) > 5):
                #For the case of 'auto' when the parameters is too large in number, we will turn off the combined plots.
                combined_plots = False 
                #For the case of 'auto' we will then turn on the individual plots.
                if self.UserInput.scatter_matrix_plots_settings['individual_plots'] == 'auto':
                    individual_plots = True
        if individual_plots == True: #This means we will return individual plots.
            if showFigure == None: showFigure = False #individual plots will not showFigure by default.
            #The below code was added by Troy Gustke and merged in to PEUQSE at end of June 2021.
            # create graph variable for plotting options
            graphs_directory = self.UserInput.directories['graphs']
            # combine all the solutions and meta data for each parameter posterior to the simulation.
            # Zip parameters contain parameter columns in dataframe, parameter names, MAP, muAP, and initial value
            finalParametersAndMetaData1 = zip(posterior_df.columns, parameterNamesAndMathTypeExpressionsDict.keys(), parameterMAPValue, parameterMuAPValue, parameterInitialValue)
            finalParametersAndMetaData2 = zip(posterior_df.columns, parameterNamesAndMathTypeExpressionsDict.keys(), parameterMAPValue, parameterMuAPValue, parameterInitialValue)
            # compare each parameter with only unique solutions
            # i and j represent the index of an abstract matrix created from comparing the parameter vectors.
            # The loop moves through the matrix and compares parameters by plotting but will only plot the bottom triangle of the matrix.
            for param_a_index, (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial) in enumerate(finalParametersAndMetaData1):
                for param_b_index, (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial) in enumerate(finalParametersAndMetaData2):
                    if param_a_index != param_b_index:
                        if self.UserInput.scatter_matrix_plots_settings['all_pair_permutations']:
                            plotting_functions.createScatterPlot(posterior_df[param_a_column], posterior_df[param_b_column], (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial),
                                            (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial), graphs_directory, plot_settings, showFigure=showFigure)
                        else:
                            if param_a_index<param_b_index: # only use the bottom triangle of the matrix and do not use the main diagonal
                                plotting_functions.createScatterPlot(posterior_df[param_a_column], posterior_df[param_b_column], (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial),
                                            (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial), graphs_directory, plot_settings, showFigure=showFigure)
        else:
            if showFigure == None: showFigure = True #combined plots will showFigure by default.
            pd.plotting.scatter_matrix(posterior_df)
            plt.savefig(self.UserInput.directories['graphs']+plot_settings['figure_name'],dpi=plot_settings['dpi'])
            if showFigure == False:
                plt.close()
        
    def makeScatterHeatMapPlots(self, parameterSamples = [], parameterNamesAndMathTypeExpressionsDict={}, parameterNamesList =[], parameterMAPValue=[], parameterMuAPValue=[], parameterInitialValue = [], showFigure=None, plot_settings={'combined_plots':'auto'}):
        if showFigure == None: showFigure = True
        import pandas as pd #This is one of the only functions that use pandas.
        import matplotlib.pyplot as plt
        import PEUQSE.plotting_functions as plotting_functions
        if 'dpi' not in  plot_settings:  plot_settings['dpi'] = 220
        if 'figure_name' not in  plot_settings:  plot_settings['figure_name'] = 'scatter_heatmap_posterior'
        if parameterSamples  ==[] : parameterSamples = self.post_burn_in_samples
        if parameterNamesAndMathTypeExpressionsDict == {}: parameterNamesAndMathTypeExpressionsDict = self.UserInput.parameterNamesAndMathTypeExpressionsDict
        if parameterNamesList == []: parameterNamesList = self.UserInput.parameterNamesList #This is created when the parameter_estimation object is initialized.        
        if parameterMAPValue == []: parameterMAPValue = self.map_parameter_set
        if parameterMuAPValue == []: parameterMuAPValue = self.mu_AP_parameter_set
        if parameterInitialValue == []: parameterInitialValue = self.UserInput.model['InputParameterPriorValues']
        # will always be seperate plots
        posterior_df = pd.DataFrame(parameterSamples,columns=[parameterNamesAndMathTypeExpressionsDict[x] for x in parameterNamesList])
        # create graph variable for plotting options
        graphs_directory = self.UserInput.directories['graphs']
        # combine all the solutions and meta data for each parameter posterior to the simulation.
        # Zip parameters contain parameter columns in dataframe, parameter names, MAP, muAP, and initial value
        # compare each parameter with only unique solutions
        # i and j represent the index of an abstract matrix created from comparing the parameter vectors.
        # The loop moves through the matrix and compares parameters by plotting but will only plot the bottom triangle of the matrix.
        finalParametersAndMetaData1 = zip(posterior_df.columns, parameterNamesAndMathTypeExpressionsDict.keys(), parameterMAPValue, parameterMuAPValue, parameterInitialValue)
        for param_a_index, (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial) in enumerate(finalParametersAndMetaData1):
            #The zipping must be done each loop, because the enumeration can only occur one time per zip.
            finalParametersAndMetaData2 = zip(posterior_df.columns, parameterNamesAndMathTypeExpressionsDict.keys(), parameterMAPValue, parameterMuAPValue, parameterInitialValue)
            for param_b_index, (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial) in enumerate(finalParametersAndMetaData2):
                if param_a_index != param_b_index:
                    if self.UserInput.scatter_heatmap_plots_settings['all_pair_permutations']:
                        plotting_functions.createScatterHeatMapPlot(posterior_df[param_a_column], posterior_df[param_b_column], (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial),
                                    (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial), graphs_directory, plot_settings) 
                    else:
                        if param_a_index<param_b_index: # only use the bottom triangle of the matrix and do not use the main diagonal
                            plotting_functions.createScatterHeatMapPlot(posterior_df[param_a_column], posterior_df[param_b_column], (param_a_column, param_a_name, param_a_MAP, param_a_mu_AP, param_a_initial),
                                        (param_b_column, param_b_name, param_b_MAP, param_b_mu_AP, param_b_initial), graphs_directory, plot_settings) 

    def createSimulatedResponsesPlots(self, allResponses_x_values=[], allResponsesListsOfYArrays =[], plot_settings={},allResponsesListsOfYUncertaintiesArrays=[], showFigure=None, flatten=False): 
        if showFigure == None: showFigure = True
        #allResponsesListsOfYArrays  is to have 3 layers of lists: Response > Responses Observed, mu_guess Simulated Responses, map_Simulated Responses, (mu_AP_simulatedResponses) > Values
        #flatten = True will convert the individual responses into a 'single response series'
        num_response_dimensions = self.UserInput.num_response_dimensions
        if flatten == True: #if we are flattening, we will have 1 response dimension at the end.
            num_response_dimensions = 1
        if allResponses_x_values == []: 
            if flatten == True:
                allResponses_x_values = np.array(self.UserInput.responses_abscissa).flatten()
                allResponses_x_values = nestedObjectsFunctions.makeAtLeast_2dNested(allResponses_x_values)
            else:
                allResponses_x_values = nestedObjectsFunctions.makeAtLeast_2dNested(self.UserInput.responses_abscissa)
        if allResponsesListsOfYArrays  ==[]: #In this case, we assume allResponsesListsOfYUncertaintiesArrays == [] also.
            allResponsesListsOfYArrays = [] #Need to make a new list in the case that there was one already, to avoid overwriting the default argument object.
            allResponsesListsOfYUncertaintiesArrays = [] #Set accompanying uncertainties list to a blank list in case it is not already one. Otherwise appending would mess up indexing.
            simulationFunction = self.UserInput.simulationFunction #Do NOT use self.UserInput.model['simulateByInputParametersOnlyFunction']  because that won't work with reduced parameter space requests.
            simulationOutputProcessingFunction = self.UserInput.simulationOutputProcessingFunction #Do NOT use self.UserInput.model['simulationOutputProcessingFunction'] because that won't work with reduced parameter space requests.
            
            #We already have self.UserInput.responses_observed, and will use that below. So now we get the simulated responses for the guess, MAP, mu_ap etc.
            
            #Get mu_guess simulated output and responses. 
            self.mu_guess_SimulatedOutput = simulationFunction( self.UserInput.InputParameterInitialGuess) #Do NOT use self.UserInput.model['InputParameterInitialGuess'] because that won't work with reduced parameter space requests.
            #Make in internal variable in case we need to flatten.
            mu_guess_SimulatedOutput = copy.deepcopy(self.mu_guess_SimulatedOutput)
            if type(simulationOutputProcessingFunction) == type(None):
                mu_guess_SimulatedResponses = mu_guess_SimulatedOutput
                if flatten == True:
                    mu_guess_SimulatedResponses = np.array(mu_guess_SimulatedResponses).flatten()
                mu_guess_SimulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(mu_guess_SimulatedResponses)
                mu_guess_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_guess_SimulatedResponses)
            if type(simulationOutputProcessingFunction) != type(None):
                mu_guess_SimulatedResponses = simulationOutputProcessingFunction(mu_guess_SimulatedOutput)
                if flatten == True:
                    mu_guess_SimulatedResponses = np.array(mu_guess_SimulatedResponses).flatten()
                mu_guess_SimulatedResponses =  nestedObjectsFunctions.makeAtLeast_2dNested(mu_guess_SimulatedResponses)     
                mu_guess_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_guess_SimulatedResponses)
            #Check if we have simulation uncertainties, and populate if so.
            if type(self.UserInput.responses_simulation_uncertainties) != type(None):
                #make an internal variable in case we need to flatten.
                mu_guess_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.UserInput.InputParameterInitialGuess)
                if flatten == True:
                    mu_guess_responses_simulation_uncertainties = np.array(mu_guess_responses_simulation_uncertainties).flatten()
                mu_guess_responses_simulation_uncertainties = nestedObjectsFunctions.makeAtLeast_2dNested(mu_guess_responses_simulation_uncertainties)
                mu_guess_responses_simulation_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_guess_responses_simulation_uncertainties)
            #Get map simiulated output and simulated responses.
            self.map_SimulatedOutput = simulationFunction(self.map_parameter_set)
            #Make an internal variable in case we need to flatten.
            map_SimulatedOutput = copy.deepcopy(self.map_SimulatedOutput)
            if type(simulationOutputProcessingFunction) == type(None):
                map_SimulatedResponses = map_SimulatedOutput
                if flatten == True:
                    map_SimulatedResponses = np.array(map_SimulatedResponses).flatten()
                map_SimulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(map_SimulatedResponses)
                map_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(map_SimulatedResponses)
            if type(simulationOutputProcessingFunction) != type(None):
                map_SimulatedResponses = simulationOutputProcessingFunction(map_SimulatedOutput)
                if flatten == True:
                    map_SimulatedResponses = np.array(map_SimulatedResponses).flatten()
                map_SimulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(map_SimulatedResponses)
                map_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(map_SimulatedResponses)
            #Check if we have simulation uncertainties, and populate if so.
            if type(self.UserInput.responses_simulation_uncertainties) != type(None):
                #make an internal variable in case we need to flatten.
                map_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.map_parameter_set)
                if flatten == True:
                    map_responses_simulation_uncertainties = np.array(map_responses_simulation_uncertainties).flatten()
                map_responses_simulation_uncertainties = nestedObjectsFunctions.makeAtLeast_2dNested(map_responses_simulation_uncertainties)
                map_responses_simulation_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(map_responses_simulation_uncertainties)
            if hasattr(self, 'mu_AP_parameter_set'): #Check if a mu_AP has been assigned. It is normally only assigned if mcmc was used.           
                #Get mu_AP simiulated output and simulated responses.
                self.mu_AP_SimulatedOutput = simulationFunction(self.mu_AP_parameter_set)
                #Make an internal variable in case we need to flatten.
                mu_AP_SimulatedOutput = copy.deepcopy(self.mu_AP_SimulatedOutput)
                if type(simulationOutputProcessingFunction) == type(None):
                    mu_AP_SimulatedResponses = mu_AP_SimulatedOutput
                    if flatten == True:
                        mu_AP_SimulatedResponses = np.array(mu_AP_SimulatedResponses).flatten()  
                    mu_AP_SimulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(mu_AP_SimulatedResponses)
                    mu_AP_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_AP_SimulatedResponses)
                if type(simulationOutputProcessingFunction) != type(None):
                    mu_AP_SimulatedResponses = simulationOutputProcessingFunction(mu_AP_SimulatedOutput)
                    if flatten == True:
                        mu_AP_SimulatedResponses = np.array(mu_AP_SimulatedResponses).flatten()  
                    mu_AP_SimulatedResponses = nestedObjectsFunctions.makeAtLeast_2dNested(mu_AP_SimulatedResponses)
                    mu_AP_SimulatedResponses = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_AP_SimulatedResponses)
                #Check if we have simulation uncertainties, and populate if so.
                if type(self.UserInput.responses_simulation_uncertainties) != type(None):
                    #make an internal variable in case we need to flatten.
                    mu_AP_responses_simulation_uncertainties = self.get_responses_simulation_uncertainties(self.mu_AP_parameter_set)
                    if flatten == True:
                        mu_AP_responses_simulation_uncertainties = np.array(mu_AP_responses_simulation_uncertainties).flatten()
                    mu_AP_responses_simulation_uncertainties = nestedObjectsFunctions.makeAtLeast_2dNested(mu_AP_responses_simulation_uncertainties)
                    mu_AP_responses_simulation_uncertainties = nestedObjectsFunctions.convertInternalToNumpyArray_2dNested(mu_AP_responses_simulation_uncertainties)
            #make internal variables for responses_observed and responses_observed_uncertainties in case we need to flatten them.
            responses_observed = copy.deepcopy(self.UserInput.responses_observed)
            responses_observed_uncertainties = copy.deepcopy(self.UserInput.responses_observed_uncertainties)
            if flatten == True: #flatten and then nest, as needed. ("lazy" way)
                responses_observed = [np.array(responses_observed).flatten()]
                responses_observed_uncertainties = [np.array(responses_observed_uncertainties).flatten()]
            #Now to populate the allResponsesListsOfYArrays and the allResponsesListsOfYUncertaintiesArrays
            #TODO: the below "if else" statements for "if num_response_dimensions == 1: " and "if num_response_dimensions > 1" are probably not necessary anymore. The code looks like it is now identical for both cases. If so, these if/else statements should be removed and the duplicate code removed.
            for responseDimIndex in range(num_response_dimensions):
                if not hasattr(self, 'mu_AP_parameter_set'): #Check if a mu_AP has been assigned. It is normally only assigned if mcmc was used.    
                    if num_response_dimensions == 1: 
                        listOfYArrays = [responses_observed[responseDimIndex], mu_guess_SimulatedResponses[responseDimIndex], map_SimulatedResponses[responseDimIndex]]        
                        allResponsesListsOfYArrays.append(listOfYArrays)
                        #Now to do uncertainties, there are two cases. First case is with only observed uncertainties and no simulation ones.
                        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                            allResponsesListsOfYUncertaintiesArrays.append( [responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                        else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                            allResponsesListsOfYUncertaintiesArrays.append([responses_observed_uncertainties[responseDimIndex],mu_guess_responses_simulation_uncertainties[responseDimIndex],map_responses_simulation_uncertainties[responseDimIndex]]) #We need to give a list for each response dimension.
                    elif num_response_dimensions > 1: 
                        listOfYArrays = [responses_observed[responseDimIndex], mu_guess_SimulatedResponses[responseDimIndex], map_SimulatedResponses[responseDimIndex]]        
                        allResponsesListsOfYArrays.append(listOfYArrays)
                        #Now to do uncertainties, there are two cases. First case is with only observed uncertainties and no simulation ones.
                        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                            allResponsesListsOfYUncertaintiesArrays.append( [responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                        else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                            allResponsesListsOfYUncertaintiesArrays.append([responses_observed_uncertainties[responseDimIndex],mu_guess_responses_simulation_uncertainties[responseDimIndex],map_responses_simulation_uncertainties[responseDimIndex]]) #We need to give a list for each response dimension.                    
                if hasattr(self, 'mu_AP_parameter_set'):
                    if num_response_dimensions == 1: 
                        listOfYArrays = [responses_observed[responseDimIndex], mu_guess_SimulatedResponses[responseDimIndex], map_SimulatedResponses[responseDimIndex], mu_AP_SimulatedResponses[responseDimIndex]]        
                        allResponsesListsOfYArrays.append(listOfYArrays)
                        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                            allResponsesListsOfYUncertaintiesArrays.append( [responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                        else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                            allResponsesListsOfYUncertaintiesArrays.append([responses_observed_uncertainties[responseDimIndex],mu_guess_responses_simulation_uncertainties[responseDimIndex],map_responses_simulation_uncertainties[responseDimIndex],mu_AP_responses_simulation_uncertainties[responseDimIndex]]) #We need to give a list for each response dimension.                                        
                    elif num_response_dimensions > 1: 
                        listOfYArrays = [responses_observed[responseDimIndex], mu_guess_SimulatedResponses[responseDimIndex], map_SimulatedResponses[responseDimIndex], mu_AP_SimulatedResponses[responseDimIndex]]        
                        allResponsesListsOfYArrays.append(listOfYArrays)
                        if type(self.UserInput.responses_simulation_uncertainties) == type(None): #This means there are no simulation uncertainties. So for each response dimension, there will be a list with only the observed uncertainties in that list.
                            allResponsesListsOfYUncertaintiesArrays.append( [responses_observed_uncertainties[responseDimIndex]] ) #Just creating nesting, we need to give a list for each response dimension.
                        else: #This case means that there are some responses_simulation_uncertainties to include, so allResponsesListsOfYUncertaintiesArrays will have more dimensions *within* its nested values.
                            allResponsesListsOfYUncertaintiesArrays.append([responses_observed_uncertainties[responseDimIndex],mu_guess_responses_simulation_uncertainties[responseDimIndex],map_responses_simulation_uncertainties[responseDimIndex],mu_AP_responses_simulation_uncertainties[responseDimIndex]]) #We need to give a list for each response dimension. 

        if plot_settings == {}: 
            plot_settings = self.UserInput.simulated_response_plot_settings
            if 'legendLabels' not in plot_settings: #The normal case:
                if hasattr(self, 'mu_AP_parameter_set'): 
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'MAP','mu_AP']
                else: #Else there is no mu_AP.
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'MAP']
                if hasattr(self, "opt_SSR"): #This means we are actually doing an optimization, and self.opt_SSR exists.
                    plot_settings['legendLabels'] = ['observed',  'mu_guess', 'CPE']
            #Other allowed settings are like this, but will be fed in as simulated_response_plot_settings keys rather than plot_settings keys.
            #plot_settings['x_label'] = 'T (K)'
            #plot_settings['y_label'] = r'$rate (s^{-1})$'
            #plot_settings['y_range'] = [0.00, 0.025] #optional.
            #plot_settings['figure_name'] = 'tprposterior'
        if 'figure_name' not in plot_settings:
            plot_settings['figurename'] = 'Posterior'
        import PEUQSE.plotting_functions as plotting_functions
        setMatPlotLibAgg(self.UserInput.plotting_ouput_settings['setMatPlotLibAgg'])
        allResponsesFigureObjectsList = []
        for responseDimIndex in range(num_response_dimensions): #TODO: Move the exporting out of the plot creation and/or rename the function and possibly have options about whether exporting graph, data, or both.
            #Some code for setting up individual plot settings in case there are multiple response dimensions.
            individual_plot_settings = copy.deepcopy(plot_settings) #we need to edit the plot settings slightly for each plot.
            if num_response_dimensions == 1:
                responseSuffix = '' #If there is only 1 dimension, we don't need to add a suffix to the files created. That would only confuse people, unless it is a 'flattened' case, in which case we will add a suffix of 'combined'.
                if flatten == True:
                    responseSuffix = '_combined'
            if num_response_dimensions > 1:
                responseSuffix = "_"+str(responseDimIndex)
            individual_plot_settings['figure_name'] = individual_plot_settings['figure_name']+responseSuffix
            if 'x_label' in plot_settings:
                if type(plot_settings['x_label']) == type(['list']) and len(plot_settings['x_label']) > 1: #the  label can be a single string, or a list of multiple response's labels. If it's a list of greater than 1 length, then we need to use the response index.
                    individual_plot_settings['x_label'] = plot_settings['x_label'][responseDimIndex]
            if 'y_label' in plot_settings:
                if type(plot_settings['y_label']) == type(['list']) and len(plot_settings['y_label']) > 1: #the  label can be a single string, or a list of multiple response's labels. If it's a list of greater than 1 length, then we need to use the response index.
                    individual_plot_settings['y_label'] = plot_settings['y_label'][responseDimIndex]                
            #TODO, low priority: we can check if x_range and y_range are nested, and thereby allow individual response dimension values for those.                               
            numberAbscissas = np.shape(allResponses_x_values)[0]
            #We have a separate abscissa for each response.              
            figureObject = plotting_functions.createSimulatedResponsesPlot(allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex], individual_plot_settings, listOfYUncertaintiesArrays=allResponsesListsOfYUncertaintiesArrays[responseDimIndex], directory = self.UserInput.directories['graphs'], showFigure=showFigure)
               # np.savetxt(self.UserInput.directories['logs_and_csvs']+individual_plot_settings['figure_name']+".csv", np.vstack((allResponses_x_values[responseDimIndex], allResponsesListsOfYArrays[responseDimIndex])).transpose(), delimiter=",", header='x_values, observed, sim_initial_guess, sim_MAP, sim_mu_AP', comments='')
            allResponsesFigureObjectsList.append(figureObject)
        return allResponsesFigureObjectsList  #This is a list of matplotlib.pyplot as plt objects.
    
    def createMumpcePlots(self, showFigure=None):
        if showFigure == None: showFigure = False
        import PEUQSE.plotting_functions as plotting_functions
        setMatPlotLibAgg(self.UserInput.plotting_ouput_settings['setMatPlotLibAgg'])
        from PEUQSE.plotting_functions import plotting_functions_class
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
        if len(self.UserInput.contour_plot_settings['active_parameters']) == 0:
            numParams = len(self.UserInput.model_parameter_info)
            active_parameters = np.linspace(0, numParams-1, numParams) #just a list of whole numbers.
            active_parameters = np.array(active_parameters, dtype='int')
        else:
            active_parameters = self.UserInput.contour_plot_settings['active_parameters']
        #TODO: reduce active_parameters by anything that has been set as a constant.
        pairs_of_parameter_indices = self.UserInput.contour_plot_settings['parameter_pairs']
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
        #Below we populate any custom fields as necessary. These go into a separate argument when making mumpce plots
        #Because these are basically arguments for a 'patch' on mumpce made by A. Savara and E. Walker.
        contour_settings_custom = {}
        # max_x/y_ticks and num_x/y_ticks are included in the following dict to keep backwards compatability
        contour_settings_custom_fields = {'figure_name','fontsize','max_num_y_ticks','max_num_x_ticks','num_y_ticks','num_x_ticks','colormap_posterior_customized','colormap_prior_customized','contours_normalized','colorbars','axis_limits','dpi', 'num_pts_per_axis','cmap_levels', 'space_between_subplots', 'zoom_std_devs', 'x_ticks', 'y_ticks', 'center_on'} #This is a set, not a dictionary.
        for custom_field in contour_settings_custom_fields:
            if custom_field in self.UserInput.contour_plot_settings:
                contour_settings_custom[custom_field] = self.UserInput.contour_plot_settings[custom_field]        
        #The colormap fields need to be removed if they are set to the default, because the default coloring is set in the mumpce class when they are not provided.
        if 'colormap_posterior_customized' in contour_settings_custom:
            if contour_settings_custom['colormap_posterior_customized'].lower() == 'default' or  contour_settings_custom['colormap_posterior_customized'].lower() == 'auto':
                del contour_settings_custom['colormap_posterior_customized']
        if 'colormap_prior_customized' in contour_settings_custom:
            if contour_settings_custom['colormap_prior_customized'].lower() == 'default' or contour_settings_custom['colormap_prior_customized'].lower() == 'auto':
                del contour_settings_custom['colormap_prior_customized']
        baseFigureName = contour_settings_custom['figure_name']
        #First make individual plots if requested.
        if self.UserInput.contour_plot_settings['individual_plots'] == 'auto':
            individual_plots = True
        else:
            individual_plots = self.UserInput.contour_plot_settings['individual_plots']
        if individual_plots == True:
            for pair in pairs_of_parameter_indices:
                contour_settings_custom['figure_name'] = self.UserInput.directories['graphs'] + baseFigureName + "__" + str(pair).replace('[','').replace(']','').replace(',','_').replace(' ','')
                figureObject_beta.mumpce_plots(model_parameter_info = self.UserInput.model_parameter_info, active_parameters = active_parameters, pairs_of_parameter_indices = [pair], posterior_mu_vector = posterior_mu_vector, posterior_cov_matrix = posterior_cov_matrix, prior_mu_vector = np.array(self.UserInput.mu_prior), prior_cov_matrix = self.UserInput.covmat_prior, contour_settings_custom = contour_settings_custom, showFigure=showFigure,)               
        #now make combined plots if requested.
        if self.UserInput.contour_plot_settings['combined_plots'] == 'auto':
            if len(pairs_of_parameter_indices) > 5:
                combined_plots = False
            else:
                combined_plots = True
        if combined_plots == True:
            contour_settings_custom['figure_name'] = self.UserInput.directories['graphs']+baseFigureName + "__combined"
            figureObject_beta.mumpce_plots(model_parameter_info = self.UserInput.model_parameter_info, active_parameters = active_parameters, pairs_of_parameter_indices = pairs_of_parameter_indices, posterior_mu_vector = posterior_mu_vector, posterior_cov_matrix = posterior_cov_matrix, prior_mu_vector = np.array(self.UserInput.mu_prior), prior_cov_matrix = self.UserInput.covmat_prior, contour_settings_custom = contour_settings_custom, showFigure=showFigure)
        return figureObject_beta

    createContourPlots = createMumpcePlots #This is just a pointer.

    @CiteSoft.after_call_compile_consolidated_log(compile_checkpoints=True) #This is from the CiteSoft module.
    def createAllPlots(self, verbose = False, showFigure=None):
        #if showFigure is none, then the default showFigure choices will occur for each of the plotting functions.
        print("Creating all plots for PEUQSE PE_object...")
        if self.UserInput.request_mpi == True: #need to check if UserInput.request_mpi is on, since if so we will only make plots after the final process.
            import PEUQSE.parallel_processing
            if PEUQSE.parallel_processing.finalProcess == True:
                pass#This will proceed as normal.
            elif PEUQSE.parallel_processing.finalProcess == False:
                return False #this will stop the plots creation.

        try:
            self.makeHistogramsForEachParameter(showFigure=showFigure)               
            if verbose: print("Finished with make histograms function call.")
        except:
            print("Unable to make histograms plots. This usually means your model is not returning simulated results for most of the sampled parameter possibilities.")


        try:
            self.makeSamplingScatterMatrixPlot(plot_settings=self.UserInput.scatter_matrix_plots_settings, showFigure=showFigure)             
            if verbose: print("Finished with makeSamplingScatterMatrixPlot function call.")
        except:
            print("Unable to make scatter matrix plot. This usually means your run is not an MCMC run, or that the sampling did not work well. If you are using Metropolis-Hastings, try EnsembleSliceSampling or try a uniform distribution multistart.")


        try:
            self.makeScatterHeatMapPlots(plot_settings=self.UserInput.scatter_heatmap_plots_settings, showFigure=showFigure)
            if verbose: print("Finished with make scatter heatmaps function call.")
        except:
            print("Unable to make scatter heatmap plots. This usually means your run is not an MCMC run, or that the sampling did not work well. If you are using Metropolis-Hastings, try one of the other samplers: EnsembleSliceSampling,  EnsembleJumpSampling,  astroidal distribution multistart, or uniform distribution multistart.")

        try:        
            self.createContourPlots(showFigure=showFigure)
            if verbose: print("Finished with create contour plots function call.")
        except:
            print("Unable to make contour plots. This usually means your run is not an MCMC run. However, it could mean that your prior and posterior are too far from each other for plotting.  You can change contour_plot_settings['colobars'] to false and can also change the contour_plot_settings['axis_limits'] if you know which region you wish to have plotted.")

        try:
            self.createSimulatedResponsesPlots(allResponses_x_values=[], allResponsesListsOfYArrays =[], plot_settings={},allResponsesListsOfYUncertaintiesArrays=[], showFigure=showFigure) #forcing the arguments to be blanks, because otherwise it might use some cached values.
            if verbose: print("Finished with create simulated responses plots function call.")
        except:
            print("Unable to make simulated response plots. This is unusual and typically means your observed values and simulated values are not the same array shape. If so, that needs to be fixed.")
            pass


        #Now we will call createSimulatedResponsesPlots again with flatten = True so that the series get plotted. This should only occur if all responses are scalars.
        try:
            self.createSimulatedResponsesPlots(allResponses_x_values=[], allResponsesListsOfYArrays =[], plot_settings={},allResponsesListsOfYUncertaintiesArrays=[], flatten = True) #forcing the arguments to be blanks, because otherwise it might use some cached values.
        except:
            print("Unable to make simulated response plots. This is unusual and typically means your observed values and simulated values are not the same array shape. If so, that needs to be fixed.")
            pass
            
        print("Finished creating all plots. Only some plots are shown on screen. The fulls set of plots are in:", self.UserInput.directories['graphs']) #TODO: take the graphs string, remove the '.' at the front if present, and print the full absolute path here.
            
    def save_to_dill(self, base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
        save_PE_object(self, base_file_name, file_name_prefix=file_name_prefix, file_name_suffix=file_name_suffix, file_name_extension=file_name_extension)
    def load_from_dill(self, base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
        theObject = load_PE_object(base_file_name, file_name_prefix=file_name_prefix, file_name_suffix=file_name_suffix, file_name_extension=file_name_extension)
        print("PE_object.load_from_dill executed. This function returns a new PE_object. To overwrite an existing PE_object, use PE_object = PE_object.load_from_dill(...)")
        return theObject
        
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


def convertPermutationsToSamples(permutations_MAP_logP_and_parameters_values, maxLogP=None, relativeFilteringThreshold=1E-2, priorsVector=None):
    #The relative filtering threshold removes anything which has a probability lower than that relative to maxLogP.
    #relativeFilteringThreshold should be a value between 0 and 1.
    #the permutations_MAP_logP_and_parameters_values should have the form logP, Parameter1, Parameter2, etc.
    #first get maxLogP if it's not provided.
    permutationsArray = permutations_MAP_logP_and_parameters_values
    if type(maxLogP) != type(None):
        maxLogP = maxLogP
    elif type(maxLogP) == type(None):
        maxLogP= -1*float('inf') #initializing.
        for element in permutationsArray:
            if element[0] > maxLogP:
                maxLogP = element[0]    
    #now calculate the absoluteFilteringThreshold:
    absoluteFilteringThreshold = maxLogP + np.log(relativeFilteringThreshold)    
    #Now make the samples repetitions based no the logP values.
    expandedArraysList = []
    for element in permutationsArray:
        if element[0] > absoluteFilteringThreshold:
            #If P2 is the smaller probability, here given by the absoluteFilteringThreshold, then...
            #it turns out we want P1/P2 = e^(logP1-logP2), where the logs here are all base e, which is our situation.
            #if it was base 10, we would want P1/P2 = 10^(log10(P1) -log10(P2)
            numberOfRepetitionsNeeded = np.exp(element[0]-absoluteFilteringThreshold)
            onesArray = np.ones((int(numberOfRepetitionsNeeded),len(element)))
            repeatedArray = onesArray * element
            expandedArraysList.append(repeatedArray)            
        elif element[0] < absoluteFilteringThreshold:
            pass        
    return np.vstack(expandedArraysList) #This stacks the expandedArraysList into a single array.


'''Below are a bunch of functions for Euler's Method.'''
#This takes an array of dydt values. #Note this is a local dydtArray, it is NOT a local deltaYArray.
software_name = "Integrated Production (Objective Function)"
software_version = "1.0.0"
software_unique_id = "https://doi.org/10.1016/j.susc.2016.07.001"
software_kwargs = {"version": software_version, "author": ["Aditya Savara"], "doi": "https://doi.org/10.1016/j.susc.2016.07.001", "cite": "Savara, Aditya. 'Simulation and fitting of complex reaction network TPR: The key is the objective function.' Surface Science 653 (2016): 169-180."} 
@CiteSoft.function_call_cite(unique_id=software_unique_id, software_name=software_name, **software_kwargs)
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
    lastIndex = len(y_values)-1
    delta_y_numerical = np.diff(np.insert(y_values,lastIndex,y_values[lastIndex])) #The diff command gives one less than what is fed in, so we insert the last value again. This gives a final value derivative of 0.
    delta_y_numerical[lastIndex] = last_point_derivative #now we set that last point to the optional argument.
    #It is ASSUMED that the t_values are evenly spaced.
    delta_t = t_values[1]-t_values[0]
    dydtNumerical = delta_y_numerical/delta_t
    return dydtNumerical
'''End of functions related to Euler's Method'''

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
    #The uncertainties, whether transformed or not, must be one of the folllowing: a) for a single dimension response can be a 1D array of standard deviations, b) for a single dimension response can be a covmat already (so already variances), c) for a multidimensional response we *only* support standard deviations at this time.
    if numResponseDimensions == 1:
        shapedUncertainties = np.array(uncertainties, dtype="float") #Initializing variable. 
        if np.shape(shapedUncertainties)[0] == (1): #This means it's just a list of standard deviations and needs to be squared to become variances.
            shapedUncertainties = np.square(shapedUncertainties) # Need to square standard deviations to make them into variances.
        else:
            shapedUncertainties = shapedUncertainties
    elif numResponseDimensions > 1:  #if the dimensionality of responses is greater than 1, we need to go through each one separately to check.
        for responseIndex in range(numResponseDimensions):
            shapedUncertainties = np.array(uncertainties, dtype="object") #Filling variable.   
            if np.shape(shapedUncertainties[responseIndex])[0] == (1): #This means it's just a list of standard deviations and needs to be squared to become variances.
                shapedUncertainties[responseIndex] = np.square(shapedUncertainties[responseIndex]) # Need to square standard deviations to make them into variances.
            else:
                shapedUncertainties[responseIndex] = shapedUncertainties[responseIndex]
    return shapedUncertainties

def boundsCheck(values, valuesBounds, boundsType):
    #Expects three arguments.
    #the first two are 1D array like arguments (values and a set of *either* upper bounds or lower bounds)
    #The third argumment is the type of bounds, either 'upper' or 'lower'
    #In practice, this means the function usually needs to be called twice.
    #A "None" type is expected for something that is not bounded in that direction. 
    
    #We first need to make arrays and remove anything that is None in the bounds.
    values = np.array(values).flatten()
    valuesBounds = np.array(valuesBounds).flatten()
    #to remove, we use brackets that pull out the indices where the comparison is not None. This is special numpy array syntax.
    parametersTruncated = values[valuesBounds !=  None].flatten() #flattening because becomes mysteriously nested.  On 6/28/22, removed the type call since python behavior changed. The line used to be: parametersTruncated = values[type(valuesBounds) != type(None)].flatten()
    parametersBoundsTruncated = valuesBounds[valuesBounds !=  None].flatten() #flattening because becomes mysteriously nested. On 6/28/22, removed the type call since python behavior changed. The line used to be: parametersBoundsTruncated = valuesBounds[type(valuesBounds) != type(None)].flatten()
    if boundsType.lower() == 'upper': #we make the input into lower case before proceeding.
        upperCheck = parametersTruncated <= parametersBoundsTruncated #Check if all are smaller.
        if False in upperCheck: #If any of them failed, we return False.
            return False
        else:
            pass #else we do the lower bounds check next.
    if boundsType.lower() == 'lower':
        lowerCheck = parametersTruncated >= parametersBoundsTruncated #Check if all are smaller.
        if False in lowerCheck: #If any of them failed, we return False.
            return False
        else:
            pass
    return True #If we have gotten down to here without returning False, both checks have passed and we return true.

def arrayThresholdFilter(inputArray, filterKey=[], thresholdValue=0, removeValues = 'below', transpose=False):
    #The thesholdFilter function takes an array and removes rows according to a filter key and thresholdValue.
    #The filterKey should be a 1D array and will be taken as the first column of the array if not provided.
    #The function finds where the filterKey is above or below the thresholdValue and then removes those rows from the original array.
    #removeValues can be "above" or "below".  
    if len(inputArray) == 0: #This should not happen for normal usage, but it has been observed in practice.
        return inputArray
    if transpose == True: #This is meant for 2D arrays.
        inputArray = np.array(inputArray).transpose()
    if len(np.shape(inputArray)) == 1:
        inputArray2D = np.atleast_2d(inputArray).transpose()
    if len(filterKey) == 0:
        filterKey == inputArray[0]
    #Now some masking type things to delete the rows above a certain value.
    filteringFailures = np.zeros(np.shape(filterKey))
    if removeValues.lower() == 'above':
        filteringFailures[filterKey>thresholdValue] = 1 #False and True, where True is beyond the filter
    if removeValues.lower() == 'below':
        filteringFailures[filterKey<thresholdValue] = 1 #False and True, where True is beyond the filter
    filteringFailuresStacked = filteringFailures*1.0 #initializing this, it is going to become a mask type shape we need.
    for dataVector in range(1, np.shape(inputArray)[1]): #This "shape()[1]" gives us the number of dataVectors we need to make the masked array.
        filteringFailuresStacked = np.hstack((filteringFailuresStacked, filteringFailures))
    inputArrayThresholdMarked = inputArray*1.0 #making a copy
    inputArrayThresholdMarked[filteringFailuresStacked==True] = float('nan') #setting any rows with large values to nan. The array filteringFailuresStacked has "True" across each of those rows.
    #Now need to remove what is masked. We need to actually remove it since we'll be doing more than mean and std after this.
    #https://stackoverflow.com/questions/22032668/numpy-drop-rows-with-all-nan-or-0-values
    #The first step is apply a "mask = np.all(..., axis=1)" line.
    mask = np.all(np.isnan(inputArrayThresholdMarked), axis=1) #This tells us which rows are all nan.
    filteredArray = inputArrayThresholdMarked[~mask] #This does the filtering where rows are deleted.
    return filteredArray

@CiteSoft.after_call_compile_consolidated_log(compile_checkpoints=True)
def exportCitations():
    pass

def setMatPlotLibAgg(matPlotLibAggSetting = 'auto'):
    #choice can be 'auto', True, or Fales.
    if matPlotLibAggSetting == 'auto':
        import platform #check if the system is Linux.
        if platform.system() == "Linux":
            matPlotLibAggSetting = True
        else:
            matPlotLibAggSetting = False
    if matPlotLibAggSetting == False:
        return #do nothing.
    if matPlotLibAggSetting == True:
        import PEUQSE.plotting_functions as plotting_functions
        plotting_functions.matplotlib.use('Agg') #added by A. Savara June 29th, 2021.
        #no return needed.

def pickleAnObject(objectToPickle, base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.pkl'):
    try:
        import pickle
        base_file_name = base_file_name.replace(file_name_extension, "") #remove the pkl extension if it’s already there, then we will add it back.
        data_filename = file_name_prefix + base_file_name + file_name_suffix + file_name_extension
        with open(data_filename, 'wb') as picklefile:
            pickle.dump(objectToPickle, picklefile)
    except Exception as theError:
        print("pickleAnObject was unable to perform pickling. The error was:", theError)

def unpickleAnObject(base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.pkl'):
    import pickle
    base_file_name = base_file_name.replace(file_name_extension, "") #remove the pkl extension if it’s already there, then we will add it back.
    data_filename = file_name_prefix + base_file_name + file_name_suffix + file_name_extension
    with open(data_filename, 'rb') as picklefile:
        theObject = pickle.load(picklefile)
    return theObject


def dillpickleAnObject(objectToPickle, base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
    #Can't use pickle. Need to use dill.
    try:
        import dill
    except:
        print("To use this feature requires dill. If you don't have it, open an anaconda prompt and type 'pip install dill' or use conda install. https://anaconda.org/anaconda/dill")
    base_file_name = base_file_name.replace(file_name_extension, "") #remove the pkl extension if it’s already there, then we will add it back.
    data_filename = file_name_prefix + base_file_name + file_name_suffix + file_name_extension
    try:
        with open(data_filename, 'wb') as picklefile:
            dill.dump(objectToPickle, picklefile)
    except Exception as theError:
        print("dillpickleAnObject was unable to perform dill pickling. The error was:", theError)
        

def unDillpickleAnObject(base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
    try:
        import dill
    except:
        print("To use this feature requires dill. If you don't have it, open an anaconda prompt and type 'pip install dill' or use conda install. https://anaconda.org/anaconda/dill")
    base_file_name = base_file_name.replace(file_name_extension, "") #remove the pkl extension if it’s already there, then we will add it back.
    data_filename = file_name_prefix + base_file_name + file_name_suffix + file_name_extension
    with open(data_filename, 'rb') as picklefile:
        theObject = dill.load(picklefile)
    return theObject

def save_PE_object(objectToPickle, base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
    dillpickleAnObject(objectToPickle, base_file_name, file_name_prefix=file_name_prefix, file_name_suffix=file_name_suffix, file_name_extension=file_name_extension)

def load_PE_object(base_file_name, file_name_prefix ='',  file_name_suffix='', file_name_extension='.dill'):
    theObject = unDillpickleAnObject(base_file_name, file_name_prefix=file_name_prefix, file_name_suffix=file_name_suffix, file_name_extension=file_name_extension)
    return theObject

def deleteAllFilesInDirectory(mydir=''):
    import os
    import copy
    if mydir == '':
        working_dir=os.getcwd()
        mydir = working_dir
    filelist = copy.deepcopy(os.listdir(mydir))
    for f in filelist:
        os.remove(os.path.join(mydir, f))

def getPointsNearExistingSample(numPointsToGet, existingSamples, logP_value = None, parameters_values = None, sortBy = 'relative_delta', pickleFileName='pointsNearExistingSample.pkl', unique_points=True):
    """
    This function retrieves a number of points near a point in existingSamples, either based on a logP_value or parameters_values.
    All of the rows in existingSamples must be of the form of LogP in the first column and parameters_values in the later columns.
    Either a logP_value can be provided to get the nearest points, or a parameters_values must be provided.
    if both are provided, the parameters_values vector will be how the 'reference' point will be chosen.

    :param numPointsToGet: Number of points to generate around point of interest. (:type: int)
    :param existingSamples: Array of existing points or string of CSV or PKL file that contains previous samples from another BPE run. It is recommended to import mcmc_logP_and_parameter_samples.csv. (:type: np.array or str)
    :param logP_value: LogP value to center points around. logP_value or parameters_value must be used as criteria. (:type: float)
    :param parameter_values: Parameter values to center points around. logP_value or parameters_value must be used as criteria. (:type: list or np.array)
    :param sortBy: Define the criteria to sort. Options are ['relative_delta', 'absolute_delta', 'absolute_distance', 'relative_distance']. (:type: str)
    :param pickleFileName: The name to export the pickle file of the nearest points. (:type: str)

    :return extracted_parameter_samples: The parameter sets that meet the defined criteria. Shape is (numPointsToGet, numParameters). (:type: np.array)
    :return extracted_logP_values: The logP values of the associated parameter sets. (:type: np.array)
    :return extracted_objective_values: The objective values to order the parameters on the defined criteria. (:type: np.array)
    """
    
    #The sorting can occur by 'relative_delta', 'absolute_delta', 'absolute_distance', 'relative_distance'. 

    #First check if existingSamples is a string. If it's a string, we assume it's a filename, including extension.
    if isinstance(existingSamples, str):
        if ".csv" in existingSamples:
            existingSamples = np.genfromtxt(existingSamples, delimiter=',', dtype=float)
        elif ".pkl" in existingSamples:
            existingSamples = unpickleAnObject(existingSamples)
    else: #else it is already an array like object.
       pass 

    #Now, we need to remove any rows where existingSamples has '-inf' in the first column, because we don't want to sample from points with 0 probability:
    existingSamples = existingSamples[ existingSamples[:,0] > float('-inf')] #this is numpy syntax for returning a filtered array with a particular condition, and here the condition is that for all rows the first value is > float('-inf').

    if not isinstance(parameters_values, type(None)):
        referencePoint = parameters_values
    else: 
        # find Index where existingSamples[:,0] = logP_value, takes the first value encountered that has the lowest value. This handles equal differences.
        indexToUse = np.abs(existingSamples[:,0] - logP_value).argmin() # finds the value that is closest to the input logP value
        parameters_values = existingSamples[indexToUse]
        referencePoint = parameters_values[1:]

    #now we will add a 0 in front to represent a logP_value. This will come in useful during subtractions.
    referencePointWithZeroInFront =  np.hstack(([0],referencePoint))

    #Now calculate the sortBy array.
    sortByCalculationsArray = copy.deepcopy(existingSamples)
    if sortBy == 'absolute_delta': #take the deltas relative to the reference point, then take their absolute values.
        sortByCalculationsArray = (existingSamples - referencePointWithZeroInFront)
        sortByCalculationsArray[:,1:] = np.abs(sortByCalculationsArray[:,1:])
    if sortBy == 'relative_delta': #take the deltas relative to the reference point, then take their absolute values, then divide by absolute magnitudes.
        sortByCalculationsArray = (existingSamples - referencePointWithZeroInFront)
        sortByCalculationsArray[:,1:] = np.abs(sortByCalculationsArray[:,1:])/ np.abs(referencePointWithZeroInFront[1:]) #divide by all of the parameter magnitudes.
    # Distance is calculated according to Sqrt(Sum((xi – xi0)^2)). This section only squares the residual.
    if sortBy == 'absolute_distance': # take distance relative to reference point
        sortByCalculationsArray = (existingSamples - referencePointWithZeroInFront)**2
        sortByCalculationsArray[:,1:] = np.abs(sortByCalculationsArray[:,1:])
    if sortBy == 'relative_distance': #take distance relative to reference point
        sortByCalculationsArray = (existingSamples - referencePointWithZeroInFront)**2
        sortByCalculationsArray[:,1:] = np.abs(sortByCalculationsArray[:,1:])/ np.abs(referencePointWithZeroInFront[1:]**2) #this converts into relative distance since (xi - xi0)**2/xi0**2 = ((xi - xi0)/xi)**2

    # we must sum over the rows only, so axis=1 is necessary. 
    # The axis is added to allow for hstack to work. Must be a 2d array.
    # Sum((xi – xi0)^2) is guided by:
    objectiveFunctionArray = np.sum(sortByCalculationsArray[:,1:], axis=1)[:, np.newaxis] # add an axis to make 2d array

    # check if criteria is by distance, if so then wrap the evaluation in a sqrt. 
    # Sqrt(Sum((xi – xi0)^2)) is guided by:
    if 'distance' in sortBy:
        objectiveFunctionArray = np.sqrt(objectiveFunctionArray)

    #now we can stack this in front of the existingSamples for sorting.
    arrayToSort = np.hstack((objectiveFunctionArray, existingSamples))
    
    #now sort it.
    sortedArray = arrayToSort[arrayToSort[:,0].argsort()] #by default this will be smallest to largest, which is what we want.
    
    # have only unique points if unique_points is True
    if unique_points:
        sortedArray = np.unique(sortedArray, axis=0)

    extractedSamples_with_objective_function_and_logP = sortedArray[0:numPointsToGet].T #extract the relevant rows and transform to appropriate shape of (numPointsToGet, numParameters)

    # seperate values into separate variables. 
    extracted_objective_values = extractedSamples_with_objective_function_and_logP[0]
    extracted_logP_values = extractedSamples_with_objective_function_and_logP[1]
    extracted_parameter_samples = extractedSamples_with_objective_function_and_logP[2:].T

    #pickle pointes to easily be used as initial guess
    pickleAnObject(extracted_parameter_samples, pickleFileName)
    # return parameter samples, logP values, and objective values separately
    return extracted_parameter_samples, extracted_logP_values, extracted_objective_values

def truncateSamples(samples, parameterBoundsLower=None, parameterBoundsUpper=None, returnMask=False):
    """
    Truncate samples by bounding the parameter space. Put None in bounds that are not truncated.
    Parameter bounds are inclusive.

    :param samples: Samples after mcmc run. (:type: np.array)
    :param parameterBoundsLower: List of parameter lower bounds. Use None to indicate a bound not being applied. (:type: list)
    :param parameterBoundsUpper: List of parameter upper bounds. Use None to indicate a bound not being applied. (:type: list)
    :param returnMask: Boolean value to return the mask along with the samples. (:type: bool)
    """
    parameterBounds = zip(parameterBoundsLower, parameterBoundsUpper)
    # check length of parameter bounds
    if samples.shape[1] != len(parameterBounds):
        print('The samples must have shape (numSamples, numParameters) and parameterBounds must be a list of tuples containing the lower and upper bounds for a parameter. The parameters bounds must be in the same order as the parameters ordered in the samples variable.')
        sys.exit()
    # apply a mask to truncate samples relative to parameter bounds
    for param_index, paramBounds in enumerate(parameterBounds):
        lowerBound, upperBound = paramBounds
        # check if Nones were input, this indicates to not truncate the parameter
        if (isinstance(lowerBound, type(None))) and (isinstance(upperBound, type(None))):
            continue # skips to the next iteration
        elif isinstance(lowerBound, type(None)): # create mask for lower bound
            truncatedMask = (samples[:, param_index] <= upperBound)
        elif isinstance(upperBound, type(None)): # create mask for upper bound
            truncatedMask = (samples[:, param_index] >= lowerBound)
        else: # create mask for lower and upper bound
            truncatedMask = ((samples[:, param_index] >= lowerBound) & (samples[:, param_index] <= upperBound))
        samples = samples[truncatedMask, :] # truncate number of samples from the mask
    if returnMask:
        return samples, truncatedMask
    else:
        return samples

def splitSamples(samples, parameter_indices, all_parameters_splitting_values):
    """
    Split samples at a given value(s).
    The samples will be returned in order. 
    Only supports splitting across one parameter for now. 

    :param samples: Parameter samples with shape (numSamples, numParams) (:type: np.array)
    :param parameter_indicies: Indices of parameters that are considered when splitting. Only supports one parameter. Ex: [3] (:type: list of ints)
    :param all_parameters_splitting_values: List of values to split at. Only supports splitting for one parameter. Ex: [[1,2,3]] (:type: list of lists)
    """
    # check if parameter indices is not in a list
    if isinstance(parameter_indices, int):
        parameter_indices = [parameter_indices]
    # make sure the splitting_values variable is a list of lists, or something similar
    all_parameters_splitting_values = nestedObjectsFunctions.makeAtLeast_2dNested(all_parameters_splitting_values)
    # check if parameter indices is not one, right now only one parameter can be handled
    if len(parameter_indices) != 1:
        print('This function (splitSamples) only supports splitting across one parameter for now.')
        sys.exit()
    #TODO: change code to make handle multiple parameters, currently it will not
    parameter_index = parameter_indices[0]
    splitting_values = all_parameters_splitting_values[0]
    # sort samples
    sorted_indices = np.argsort(samples[:, parameter_index], axis=0)[:, np.newaxis] # add newaxis to make samples and sorted_indices the same dimensions
    sorted_samples = np.take_along_axis(samples, sorted_indices, axis=0)
    # find splitting spots
    split_at = sorted_samples[:, parameter_index].searchsorted(splitting_values) 
    split_samples = np.split(sorted_samples, split_at) # split_at must be a list-like object
    # return list of arrays
    return split_samples

def calculateAndPlotConvergenceDiagnostics(discrete_chains_post_burn_in_samples, parameterNamesAndMathTypeExpressionsDict, plot_settings={}, graphs_directory='./', createPlots=True, showFigure = None):
    """
    Calls other convergence functions to do calculations and make plots.

    :param discrete_chains_post_burn_in_samples: Samples with specific arrays for chains. (:type: np.array)
    :param parameterNamesAndMathTypeExpressionsDict: Dictionary with parameter name and symbol (:type dict)
    :param plot_settings: Plotting settings from UserInput (:type: dict)
    :param graphs_directory: Path to save graphs. (:type: str)
    :param createPlots: Flag to create plots after convergence analysis. (:type: bool)
    """
    from warnings import catch_warnings, simplefilter
    if showFigure == None: showFigure = True
    # makes sure the plot settings is populated before plotting
    if len(plot_settings)==0:
        createPlots=False
    try:
        # use the zeus AutoCorrTime function for all calculations.
        from zeus.autocorr import AutoCorrTime
        # create window sizes that increase on a log scale. 
        window_indices_act = np.exp(np.linspace(0, np.log(discrete_chains_post_burn_in_samples.shape[0]), 21)).astype(int)[1:]
        # initialize array with shape (N_intervals, numParams)
        taus_zeus = np.empty((len(window_indices_act), discrete_chains_post_burn_in_samples.shape[2])) 
        # populate taus using zeus AutoCorrTime function where lag length is determined by Sokal 1989.
        # For more information on Integrated Autocorrelation time see https://emcee.readthedocs.io/en/stable/tutorials/autocorr/ 
        with catch_warnings(): # suppress warnings from ACT function
            simplefilter('ignore')
            for i, n in enumerate(window_indices_act): # loop through the window indices to get larger and larger windows
                # since size is (numSamples, numChains, numParameters), we pass in limited samples up to the window size index
                # while passing in every chain and parameter
                taus_zeus[i] = AutoCorrTime(discrete_chains_post_burn_in_samples[:n,:,:]) 
        # create plots using PEUQSE plotting functions file.
        from PEUQSE.plotting_functions import createAutoCorrTimePlot
        # create plots for each parameter. The parameter names and symbols are unpacked from the dictionary.
        # loop through each parameter act values to plot and assign to self convergence
        parameter_act_for_each_window = {}
        combined_parameter_act_for_each_window = np.ones((len(window_indices_act),))
        heuristic_exponent_value = 1 # for individual parameters, the heuristic line is 50tau
        for param_taus, (parameter_name, parameter_math_name) in zip(taus_zeus.T, parameterNamesAndMathTypeExpressionsDict.items()):
            # only plot if createPlots is True
            if createPlots:
                createAutoCorrTimePlot(window_indices_act, param_taus, parameter_name, parameter_math_name, heuristic_exponent_value, graphs_directory, showFigure=showFigure)
            parameter_act_for_each_window[parameter_name] = param_taus
            # combine parameters by adding log(ACT) or just by multiplying
            combined_parameter_act_for_each_window *= param_taus
        # create combined parameters plot for ACT
        heuristic_exponent_value = discrete_chains_post_burn_in_samples.shape[2] # reassign to the number of combined parameters, which is all parameters
        createAutoCorrTimePlot(window_indices_act, combined_parameter_act_for_each_window, 'Combined_Parameters', 'All Parameters', heuristic_exponent_value, graphs_directory, showFigure=showFigure)
    except Exception as theError:
        window_indices_act = None
        taus_zeus = None
        parameter_act_for_each_window = None
        print('The AutoCorrelation Time plots have failed to be created. The error was:', theError)
    try: # prevents crashing when running convergence diagnostics on short chains or weird models
        # We previously used ARVIZ version 0.11.0 with the below syntax and switched in July 2022.
        # from arviz import geweke
        # local_z_score = geweke(discrete_chains_post_burn_in_samples[:window, chain_num, param_num])
        from PEUQSE.plotting_functions import createGewekePlot
        # create a linearly space array for creating window sizes for Geweke percent diagnostic
        window_indices_geweke = np.linspace(0, discrete_chains_post_burn_in_samples.shape[0], 21).astype(int)[1:]
        # loop through each param, each chain, and each window size
        # Geweke function is called for each window size. The full window (last one) is saved for plotting.
        total_z_scores = [] # initialize list for combining parameters.
        for param_num, (parameter_name, parameter_math_name) in enumerate(parameterNamesAndMathTypeExpressionsDict.items()):
            z_scores_array_per_chain = [] # initialize list for number of windows
            for chain_num in range(discrete_chains_post_burn_in_samples.shape[1]):
                z_scores_array_per_window = [] # initialize the list
                for window in window_indices_geweke:
                    # calculate z scores for each window. Use default settings of first 10% and last 50% of points to compare.
                    with catch_warnings():
                        simplefilter('ignore')
                        local_z_score = geweke_diagnostic(discrete_chains_post_burn_in_samples[:window, chain_num, param_num])
                    # checks if it is the last window. If yes, save the indices. Save for plotting and all last windows are the same.
                    if window == window_indices_geweke[-1]:
                        z_scores_final_indices = local_z_score.T[0]
                    z_scores_array_per_window.append(local_z_score.T[1])
                z_scores_array_per_chain.append(z_scores_array_per_window)
            z_scores_array_per_chain = np.array(z_scores_array_per_chain)
            # save all chains for combining all parameters.
            total_z_scores.append(z_scores_array_per_chain)
            z_scores_array = np.mean(np.abs(z_scores_array_per_chain), axis=0)
            # use numpy function to count how many z values fall outside 1 std. Divide by total values to get percent (decimal)
            z_scores_percentage_outlier = np.count_nonzero(z_scores_array>1, axis=1) / z_scores_array.shape[1]
            # save last window for plotting.
            z_scores_final = z_scores_array[:,-1]
            z_scores_geweke_final_plot_inputs = [z_scores_final_indices, z_scores_final] # allows for easier plotting with unpacking.
            # now plot using PEUQSE.plotting function if createPlots is True
            if createPlots:
                createGewekePlot(z_scores_geweke_final_plot_inputs, window_indices_geweke, z_scores_percentage_outlier, parameter_name, parameter_math_name, graphs_directory, showFigure=showFigure)
        # get combined parameter Geweke plot
        total_z_scores = np.array(total_z_scores)
        # abs and average across the parameters.
        z_scores_sum_params = np.mean(np.abs(total_z_scores), axis=0)
        # average across each chain after params are averaged.
        z_scores_sum_params_and_chains = np.mean(z_scores_sum_params, axis=0)
        # save final window for plotting and self convergence.
        z_scores_sum_params_final = z_scores_sum_params_and_chains[:, -1]
        # use numpy function to count how many z values fall outside 1 std. Divide by total values to get percent (decimal)
        z_scores_sum_params_percentage_outlier = np.count_nonzero(z_scores_sum_params_and_chains>1, axis=1) / z_scores_sum_params_and_chains.shape[1]
        z_scores_sum_params_geweke_final_plot_inputs = [z_scores_final_indices, z_scores_sum_params_final] # allows for easier plotting with unpacking.
        # now plot using PEUQSE.plotting function if createPlots is True
        if createPlots:
            createGewekePlot(z_scores_sum_params_geweke_final_plot_inputs, window_indices_geweke, z_scores_sum_params_percentage_outlier, 'Combined_Parameters', 'All Parameters', graphs_directory, showFigure=showFigure)
    except Exception as theError:
        print('Could not calculated Geweke convergence analysis. The chain length may be too small, so more samples are recommended.')
        print('The Geweke diagnostic graphs failed to be created. The error was:', theError)
        window_indices_geweke = None
        z_scores_sum_params_final = None
        z_scores_sum_params_percentage_outlier = None
    # return both window_indicies, final ACT values each param, final ACT values each param and window, final z scores summed parameters, and final summed parameters percent outliers
    return (window_indices_act, taus_zeus[-1,:], parameter_act_for_each_window, window_indices_geweke, z_scores_sum_params_final, z_scores_sum_params_percentage_outlier)
        
def geweke_diagnostic(post_burn_in_samples, initial_window=0.1, comparison_window=0.5, intervals=20):
    """ Geweke diagnostic for convergence of MCMC sampling. 
    Z scores are compared from the initial window to a final window of the sampling. 

    :param post_burn_in_samples: Samples in MCMC sampling after burn in period. (:type: np.array)
    :param initial_window: Percent of initial samples to be compared. (:type: float)
    :param comparison_window: Percent of final samples to be compared to initial values. (:type: float)
    :param intervals: Number of intervals for windows to be compared. (:type: int)

    Geweke function based on logic from arviz module file diagnostics.py version=0.11.0
    """
    for interval in (initial_window, comparison_window):
        if interval <= 0 or interval >= 1:
            raise ValueError("Invalid intervals for Geweke convergence analysis:", (initial_window, comparison_window), "Must be between 0 and 1")
    if initial_window + comparison_window >= 1:
        raise ValueError("Invalid intervals for Geweke convergence analysis:", (initial_window, comparison_window), "first and last intervals should not overlap.")

    # Initialize list of z-scores
    z_scores_total = []

    # Last index value of the input samples
    last_index = len(post_burn_in_samples) - 1

    # Start intervals going up to the start of the comparison window of the chain
    last_start_index = (1 - comparison_window) * last_index

    # Calculate starting indices
    start_indices = np.linspace(0, last_start_index, num=intervals, endpoint=True, dtype=int)

    # Loop over start indices
    for start in start_indices:
        # Calculate slices
        first_slice = post_burn_in_samples[start : start + int(initial_window * (last_index - start))]
        last_slice = post_burn_in_samples[int(last_index - comparison_window * (last_index - start)) :]

        z_score = first_slice.mean() - last_slice.mean()
        z_score /= np.sqrt(first_slice.var() + last_slice.var())

        z_scores_total.append([start, z_score])

    return np.array(z_scores_total)



if __name__ == "__main__":
    pass

