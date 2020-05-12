import numpy as np
global T
global pA
global pB
global S_CO
global G_H2O
S_CO = 0.00150
G_H2O = -1.32
T = 398.15
pA = 0.1
pB = 0.1
def Langmuir_compete_ads(sample):
    kB = 8.61733035E-5 #eV/K
    delta_G_CO = sample - T*S_CO # S_CO = 0.00150
    delta_G_H2O = G_H2O #G_H2O = -1.32
    KA = np.exp(-delta_G_CO/(kB*T))
    KB = np.exp(-delta_G_H2O/(kB*T))
    theta_A = (KA*pA)/(1+KA*pA+KB*pB)
    return(theta_A)
