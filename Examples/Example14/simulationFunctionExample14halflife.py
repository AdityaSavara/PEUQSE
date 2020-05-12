import numpy as np

global upperbound
global lowerbound
upperbound = 300 #seconds
lowerbound =  0.00231 #seconds

def likelihoodZeroIfOutsideKnownRange(discreteParametersArray):
    halflife = halflifeSimulator(discreteParametersArray)   #the * is a list expansion.
    if halflife > upperbound:
        return float('-inf'), 0 #This means zero probability because log of 0 approaches negative infinity. We return a nonphysical simulated halflife of zero along with that.
    elif halflife < lowerbound:
        return float('-inf'), 0 #This means zero probability because log of 0 approaches negative infinity.  We return a nonphysical simulated halflife of zero along with that.  
    else:# implies (halflife < upperbound) and (halflife > lowerbound):
        return 1, halflife  #If the halflife is in acceptable range, we return 1.
    
#This is *only* for first order halflives.
def halflifeSimulator(KineticParameters, T=100):
    log10A = KineticParameters[0]
    Ea = KineticParameters[1]
    A = 10**(log10A)
    k = A*np.exp(-1*Ea/(8.3145*T))
    halflife = 0.693/k
    return halflife
