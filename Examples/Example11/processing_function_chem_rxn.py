import numpy as np
global T
T = 298.15
def delta_G_T_to_ln_CA(sample):
    kB = 8.61733035E-5 #eV/K
    K_eq = np.exp(-sample/(kB*T))
    ln_CA = np.log(1/(1+K_eq))
    return(ln_CA)
