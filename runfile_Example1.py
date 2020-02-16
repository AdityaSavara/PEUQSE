import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd
import copy
import CheKiPEUQ as CKPQ

if __name__ == "__main__":
    import UserInput_CKPQ_Example1 as UserInput
    UserInput.verbose = False    
    UserInput.mcmc_burn_in = 10000
    UserInput.mcmc_length = 40000
    UserInput.checkPointFrequency = 1000
    PE_object = CKPQ.parameter_estimation(UserInput)
    #[map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, samples_simulatedOutputs, logP] = PE_object.doMetropolisHastings()
    PE_object.doMetropolisHastings()
    
    PE_object.createAllPlots() #This function calls each of the below functions.
#    PE_object.makeHistogramsForEachParameter()    
#    PE_object.makeSamplingScatterMatrixPlot()
#    PE_object.createSimulatedResponsesPlot()
    
    #TODO: call the mum_pce plotting objects.