import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import scipy
from scipy.integrate import odeint
plt.rcParams.update({'font.size': 22})
log_A_mean = 1.07015244e+01
Ea_mean =  5.60317080e+04 * 1/1000 * 1/96
kB = 8.617333e-5

######################## Begin with the scalar mean only
def onecomponent(theta, t, log_A, Ea):
    T = 298.15 + t
    TOF = theta*np.exp(-(Ea-kB*T*log_A)/(kB*T))
    return -TOF

t = np.linspace(0, 300, 301) #101 points from 0 (s) to 100(s)
theta_t0 = 1.0
tpr_theta = odeint(onecomponent, theta_t0, t, args = (log_A_mean, Ea_mean))
############### We have the solution for theta, but we want to 
############### get out the solution for TOF, because the TOF
############### is what appears in the experimental measurements
T = 298.15 + t 
TOF = tpr_theta[:,0]*np.exp(-(Ea_mean-kB*T*log_A_mean)/(kB*T)) #tpr_theat
# is shape (100,1) and we need shape (100,) for the code to run.
plt.plot(t,TOF, 'r')
plt.xlabel('t (s) T=298.15 (K) + t')
plt.ylabel(r'$TOF (s^{-1})$')
plt.tight_layout()
plt.savefig('onecomponent.png', dpi=220)
####################### Now, a grid search will be conducted.
log_A_linspace = np.linspace(log_A_mean*0.9, log_A_mean*1.1, 5)

