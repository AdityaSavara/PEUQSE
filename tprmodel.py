
import numpy as np

# The below function must return a vector of rates. 
def tprequation(tpr_theta, t,Ea_1, Ea_2, log_A1, log_A2, gamma1, gamma2,beta_dTdt,start_T): #beta_dTdT is the heating rate. 
    if tpr_theta.ndim == 2: 
        theta_1 = tpr_theta[:,0] 
        theta_2 = tpr_theta[:,1] 
    else: 
        [theta_1, theta_2] = tpr_theta 
    T = start_T + beta_dTdt*t 
    kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1 
    rate_1 = theta_1*np.exp(-(Ea_1-kB*T*log_A1-gamma1*theta_1)/(kB*T)) 
    rate_2 = theta_2*np.exp(-(Ea_2-kB*T*log_A2-gamma2*theta_2)/(kB*T)) 
    return [-rate_1, -rate_2] 
    