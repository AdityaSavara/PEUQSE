# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:21:06 2020

@author: fvs
"""

import numpy as np
#The Ea will be returned in J/mol. T in K.
#k0 is the pre-exponential in s-1.
#Theta0 is the initial relative coverage
#beta_H is the heating rate in Kelvin per second.
#n is the reaction order.
#The variable "find" must be filled with either the string "Ea" or the string "Tp"
#Whether solving for Tp or Ea, a guess must be provided for the other. However, defaults are included so a guess does not actually need to be provided.
def redheadPeakMaximum_find_Ea_or_Tp(find="", k0=1E13, Theta0=1.0, beta_H=2.0, n=1, Ea=200000, Tp=500, Thetap=0.0):
    if find == "":
        print("A string containing either 'Ea' or 'Tp' must be provided as the first argument.")
        sys.exit()        
    R_constant = 8.314
    #We'll take Masel's Equation 7.59 from his 1996 book and use the left-hand side and right-hand side.
    if n==2:
        Thetap = 0.5*Theta0
    elif n==1:
        Thetap = 1.0*Theta0 #does not actually matter for n=1.
    else:
        print("only 1st order and 2nd order are currently supported, since Thetap will need to be solved for other cases.")
        sys.exit()

    def caseOfEaNotKnown(Ea): #Now in this little function Ea will be replaced by what's passed in.
        leftHandSide = Ea/(R_constant*Tp)
        rightHandSide = np.log(  (Thetap**(n-1))  *  k0*Tp*n/beta_H)-np.log(leftHandSide) #in numpy, "log" is the natural log.
        difference = leftHandSide-rightHandSide
        return difference #this is only the return of the little function, not of the outer function.
    def caseOfTpNotKnown(Tp): #Now in this little function Tp will be replaced by what's passed in.
        leftHandSide = Ea/(R_constant*Tp)
        rightHandSide = np.log(  (Thetap**(n-1))  *  k0*Tp*n/beta_H)-np.log(leftHandSide) #in numpy, "log" is the natural log.
        difference = leftHandSide-rightHandSide
        return difference #this is only the return of the little function, not of the outer function.

    import scipy.optimize
    if find=="Ea":
        solvedValue = scipy.optimize.root(caseOfEaNotKnown, Ea)
    if find=="Tp":
        solvedValue = scipy.optimize.root(caseOfTpNotKnown, Tp)   
    return float(solvedValue['x'])


# Ea_solved = redheadPeakMaximum_find_Ea_or_Tp("Ea", n=2, Theta0=0.25)
# print(Ea_solved)


# Tp_solved = redheadPeakMaximum_find_Ea_or_Tp("Tp", beta_H=10.0)
# print(Tp_solved)