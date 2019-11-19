import numpy as np
def tprequation(tpr_theta, t, Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean,dT,dt,start_T):
    if tpr_theta.ndim == 2:
        theta_1 = tpr_theta[:,0]
        theta_2 = tpr_theta[:,1]
    else:
        [theta_1, theta_2] = tpr_theta
    T = start_T + dT/dt*t
    kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1
    rate_1 = theta_1*np.exp(-(Ea1_mean-kB*T*log_A1_mean-gamma_1_mean*theta_1)/(kB*T))
    rate_2 = theta_2*np.exp(-(Ea2_mean-kB*T*log_A2_mean-gamma_2_mean*theta_2)/(kB*T))
    return [-rate_1, -rate_2]
