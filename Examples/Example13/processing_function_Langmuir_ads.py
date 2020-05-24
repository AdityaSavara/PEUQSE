import numpy as np
global T
global pA
global pB
T = 298.15
pA = 0.1
pB = 0.1
def Langmuir_compete_ads(sample):
    sample=np.array(sample)
    kB = 8.61733035E-5 #eV/K
    KA = np.exp(-sample/(kB*T))
    KB = np.exp(-(-0.1)/(kB*T))
    theta_A = (KA*pA)/(1+KA*pA+KB*pB)
    return(theta_A)
