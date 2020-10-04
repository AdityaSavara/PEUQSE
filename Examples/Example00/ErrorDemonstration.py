# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 15:10:03 2020

@author: fvs
"""
from scipy.stats import multivariate_normal
import numpy as np

simulatedResponses_transformed_flattened = [ 160500,  810500, 1440500]
observedResponses_transformed_flattened = [ 360500 , 580500 , 1620500]
comprehensive_responses_covmat = [[200000, 200000, 200000]]

log_probability_metric = multivariate_normal.logpdf(mean=simulatedResponses_transformed_flattened,x=observedResponses_transformed_flattened,cov=comprehensive_responses_covmat[0])

print(log_probability_metric)



def returnShapedResponseCovMat(numResponseDimensions, uncertainties):
    #The uncertainties, whether transformed or not, must be one of the folllowing: a) for a single dimension response can be a 1D array of standard deviations, b) for as ingle dimension response can be a covmat already (so already variances), c) for a multidimensional response we *only* support standard deviations at this time.
    print("line 1459", uncertainties)
    if numResponseDimensions == 1:
        shapedUncertainties = np.array(uncertainties) #Initializing variable. 
        if np.shape(shapedUncertainties)[0] == (1): #This means it's just a list of standard deviations and needs to be squared to become variances.
            shapedUncertainties = shapedUncertainties**2 # Need to square standard deviations to make them into variances.
        else:
            shapedUncertainties = shapedUncertainties
    elif numResponseDimensions > 1:  #if the dimensionality of responses is greater than 1, we only support providing standard deviations. Will flatten and square.
        shapedUncertainties = np.array(uncertainties) #Filling variable.  
        shapedUncertainties = shapedUncertainties.flatten() 
        shapedUncertainties = shapedUncertainties**2 #Need to square standard deviations to make them into variances.
    print("line 1470", shapedUncertainties)
    return shapedUncertainties


shapedUncertainties = returnShapedResponseCovMat(1,  [[200000, 300000, 300000]])