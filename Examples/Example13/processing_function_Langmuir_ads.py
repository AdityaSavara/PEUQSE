global connected_variables_values
#Initializing. Do *not* fill this before the function. This should remain blank *unless* it is called from outside of this module.
connected_variables_values = []


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


def Langmuir_compete_ads_connected(sample):
    if len(connected_variables_values) != 0:
        T = connected_variables_values[0]
        pA = connected_variables_values[1]
        pB = 0.1 #This example only has 2 connected variables so that the info gain is easy to plot as a 2D image.
    sample=np.array(sample)
    kB = 8.61733035E-5 #eV/K
    KA = np.exp(-sample/(kB*T))
    KB = np.exp(-(-0.1)/(kB*T))
    theta_A = (KA*pA)/(1+KA*pA+KB*pB)
    return(theta_A)



def Langmuir_replacement(sample):
    if len(connected_variables_values) != 0:
        T = connected_variables_values[0]
        pA = connected_variables_values[1]
        pB = 0.1 #This example only has 2 connected variables so that the info gain is easy to plot as a 2D image.
    sample=np.array(sample)
    kB = 8.61733035E-5 #eV/K
    K_R = np.exp(-sample/(kB*T))
    theta_A = (K_R*pA)/(1+K_R*pA)
    return(theta_A)