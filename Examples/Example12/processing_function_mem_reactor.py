import numpy as np
global T
global volume
global F0
global P0
global R
T = 298.15
volume = 1000
F0 = np.array([5, 0])
P0 = 1
R = 82.057338

def cmr(theta,V,k_1,k_minus_1,k_B,T0):
    F_A, F_B = theta
    F_T = F_A + F_B
    C_T0 = P0 / (R * T0)
    C_A = C_T0 * F_A / F_T
    C_B = C_T0 * F_B / F_T
    r_A = k_1*C_A - k_minus_1*C_B
    R_B = k_B*C_B
    dFdV = [-r_A, r_A - R_B]
    return dFdV

def mem_reactor(sample):
    kB = 8.61733035E-5 #eV/K
    sol = odeint(cmr, F0, np.linspace(0,self.volume,2), args=(sample[0], sample[1], sample[2], T))
    conc_sol_last=sol[-1,:].T 
    return(conc_sol_last[0])

