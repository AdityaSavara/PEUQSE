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


def Langmuir_compete_ads(sample): # sample contains [H_CO, S_CO, H_H2O, S_H2O]
    kB = 8.61733035E-5 #eV/K
    delta_G_CO = sample[0] - T*sample[1] # S_CO = 0.00150
    delta_G_H2O = sample[2] - T*sample[3] #G_H2O = -1.32
    KA = np.exp(-delta_G_CO/(kB*T))
    KB = np.exp(-delta_G_H2O/(kB*T))
    theta_A = (KA*pA)/(1+KA*pA+KB*pB)
    return(np.log(theta_A))
    
    
    
def Langmuir_replacement(sample): # sample contains [H_CO, S_CO, H_H2O, S_H2O]
    kB = 8.61733035E-5 #eV/K
    delta_G_rxn = sample[0] - T*sample[1] # DeltaH-T*DeltaS
    K_rxn = np.exp(-delta_G_rxn/(kB*T))
    theta_A = (K_rxn*pA)/(1+K_rxn*pA)
    return(np.log(theta_A))



def Langmuir_replacement_three_temperatures(sample): # sample contains [DeltaH_rxn, DeltaS_rxn]
    kB = 8.61733035E-5 #eV/K
    
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
    return theta_list
    
    
def Langmuir_replacement_three_temperatures(sample): # sample contains [DeltaH_rxn, DeltaS_rxn]
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
    return theta_list
    
def Langmuir_replacement_three_temperatures_deltaH_plus1(sample): # sample contains [DeltaH_rxn, DeltaS_rxn]
    kB = 8.61733035E-5 #eV/K
    sample = np.array(sample) #Make sure it's an array and not a tuple.
    sample[0] = sample[0] - 1 #We had to add 1 to deltaH to prevent numerical errors during sampling. This removes that.
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
    return theta_list
    
def Langmuir_replacement_three_temperatures_log(sample): # sample contains [DeltaH_rxn, DeltaS_rxn]
    global T
    global pA
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