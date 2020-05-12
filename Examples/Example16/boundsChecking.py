# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:29:45 2020

@author: fvs
"""

import numpy as np

upperArray = np.array([3,3,3])

middleArray = np.array([1,1,1])

lowerArray = np.array([0,0,0])

if middleArray.all() < upperArray.all():
    print(True)
else:
    print(False)
    print(middleArray)
    print(upperArray)
    
print(middleArray<upperArray)

print(middleArray.all()<upperArray.all())


comparisonOutput = middleArray<upperArray
if False in comparisonOutput:
    print("Failed")
if False not in comparisonOutput:
    print("Passed")    


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
            return True
    if boundsType.lower() == 'lower':
        lowerCheck = parametersTruncated > parametersBoundsTruncated #Check if all are smaller.
        if False in lowerCheck: #If any of them failed, we return False.
            return False
        else:
            return True

print('line 61')
print(boundsCheck(middleArray, upperArray, 'upper'))
print(boundsCheck(middleArray, upperArray, 'lower'))

middlingArray = np.array([1,6,1])
print(boundsCheck(middlingArray, upperArray, 'lower'))