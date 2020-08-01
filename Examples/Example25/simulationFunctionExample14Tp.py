from RedheadApproximation import redheadPeakMaximum_find_Ea_or_Tp

def getTpFromKineticParametersAndInitialCoverageWrapper(discreteParameterVector):
    Tp = getTpFromKineticParametersAndInitialCoverage(*discreteParameterVector)
    return Tp

def getTpFromKineticParametersAndInitialCoverage(log10A, Ea, Theta0=1.0, beta_H=2, n=1): #Thetap is only needed for second order.
    k0 = 10**log10A#k0 is the same as "A" 
    #In addition to changing log10A by "10**", this wrapper is also necessary to specify Tp seeking and to switch to named arguments.
    #Need to have a reasonable initial guess for Tp, otherwise the solver fails.
    Tp = redheadPeakMaximum_find_Ea_or_Tp(find="Tp", k0=k0, Ea=Ea, Theta0=Theta0, beta_H=beta_H, n=n, Tp=200)
    return Tp