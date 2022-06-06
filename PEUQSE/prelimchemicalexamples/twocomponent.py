# This script solves the prior mean and plots
# with the experimental data and errors.
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import pandas as pd
import scipy
from scipy.integrate import odeint
from scipy.stats import multivariate_normal
import sys
sys.path.append('..') # add one directory up to the path
from tprmodel import tprequation
plt.rcParams.update({'font.size': 22})
Ea1_mean = 41.5 # kJ mol^-1
Ea2_mean = 41.5
log_A1_mean = 13.0
log_A2_mean = 13.0
gamma_1_mean = 0.1
gamma_2_mean = 0.1
Ea1_std = 20.0
Ea2_std = 20.0
log_A1_std = 2.5
log_A2_std = 2.5
gamma_1_std = 0.1
gamma_2_std = 0.2
kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1

######################## Import experimental data
experiments_df = pd.read_csv('../ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedLargerErrors.csv')
dT = experiments_df['dT'][0] # assuming dT and dt are constant throughout
dt = experiments_df['dt'][0]
start_T = experiments_df['AcH - T'][0]
theta_1_t0 = 0.5 # half of the surface initially coverered with theta_1
theta_2_t0 = 0.5
tpr_theta = odeint(tprequation, [theta_1_t0, theta_2_t0], experiments_df['time'].to_numpy(), args = (Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean,dT,dt,start_T))

rate = tprequation(tpr_theta, experiments_df['time'].to_numpy(), Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean,dT,dt,start_T)
print(rate)
rate = np.asarray(rate)
print(rate.shape)
rate_tot = -np.sum(rate, axis=0)
fig1, ax1 = plt.subplots()
ax1.plot(experiments_df['AcH - T'].to_numpy(),rate_tot, 'r')
ax1.set_xlabel('T (K)')
ax1.set_ylabel(r'$rate (s^{-1})$')
ax1.legend(['model prior'])
fig1.tight_layout()
fig1.savefig('twocomponent.png', dpi=220)

fig2, ax2 = plt.subplots()
ax2.plot(experiments_df['AcH - T'].to_numpy(),experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000,'g')
ax2.set_xlabel('T (K)')
ax2.set_ylabel(r'$rate (s^{-1})/1000$')
ax2.legend(['experiment'])
fig2.tight_layout()
fig2.savefig('experiment.png', dpi=220)
