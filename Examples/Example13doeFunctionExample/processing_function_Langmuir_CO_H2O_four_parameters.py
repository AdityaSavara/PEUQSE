import numpy as np
global T
global pA
global pB
#global S_CO
#global G_H2O
#S_CO = 0.00150
#G_H2O = -1.32
T = 398.15
pA = 0.1
pB = 0.1
global T_spread
T_spread = 25

global connected_variables_values
connected_variables_values = [] #Initializing. Do *not* fill this before the function. This should remain blank *unless* it is called from outside of this module.

def populate_pA_and_T(valuesToUse):
    global T #This global declaration is to emphasize that we are using the same value for this across the file.
    global pA #This global declaration is to emphasize that we are using the same value for this across the file.
    global connected_variables_values
    connected_variables_values = list(valuesToUse)
    #One could probably instead use T= valuesToUse[0] for example. Or perhaps T, pA = list(valuesToUse).

def Langmuir_replacement_three_temperatures_log(sample): # sample contains [DeltaH_rxn, DeltaS_rxn]
    global T #This global declaration is to emphasize that we are using the same value for this across the file.
    global pA #This global declaration is to emphasize that we are using the same value for this across the file.
    global connected_variables_values
    if len(connected_variables_values) != 0:
        T = connected_variables_values[0]
        pA = connected_variables_values[1]
    kB = 8.61733035E-5 #eV/K
    sample = np.array(sample) #Make sure it's an array and not a tuple.
    #
    theta_list = []
    #Temperature Middle:
    delta_G_rxn = sample[0] - T*sample[1] # DeltaH-T*DeltaS
    K_rxn = np.exp(-delta_G_rxn/(kB*T))
    theta_A = (K_rxn*pA)/(1+K_rxn*pA)
    theta_list.append(theta_A)

    #Temperature Low:
    delta_G_rxn = sample[0] - (T-T_spread)*sample[1] # DeltaH-T*DeltaS
    K_rxn = np.exp(-delta_G_rxn/(kB*(T-T_spread)))
    theta_A = (K_rxn*pA)/(1+K_rxn*pA)
    theta_list.append(theta_A)
    
    #Temperature Low:
    delta_G_rxn = sample[0] - (T+T_spread)*sample[1] # DeltaH-T*DeltaS
    K_rxn = np.exp(-delta_G_rxn/(kB*(T+T_spread)))
    theta_A = (K_rxn*pA)/(1+K_rxn*pA)
    theta_list.append(theta_A)
    return np.log10(theta_list)