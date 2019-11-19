import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import scipy
from scipy.integrate import odeint
from scipy.stats import multivariate_normal
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
fig1, ax1 = plt.subplots()
ax1.plot(T,TOF, 'r')
ax1.set_xlabel('T (K)')
ax1.set_ylabel(r'$rate (s^{-1})$')
fig1.tight_layout()
fig1.savefig('onecomponent.png', dpi=220)
####################### Now, a grid search will be conducted.
log_A_linspace = np.linspace(log_A_mean*0.9, log_A_mean*1.1, 5)
Ea_linspace = np.linspace(Ea_mean*0.9, Ea_mean*1.1, 5)
############### First plot in the parameter space.
Ea_mesh, log_A_mesh = np.meshgrid(Ea_linspace, log_A_linspace)
pos = np.dstack((Ea_mesh, log_A_mesh))
cov_prior_param = np.array([[(Ea_mean*0.1)**2, 0], [0, (log_A_mean*0.1)**2]])
mu_prior_param = np.array([Ea_mean, log_A_mean])
prob_prior_param = multivariate_normal(mu_prior_param, cov_prior_param)
Z = prob_prior_param.pdf(pos)
levels = np.unique(Z.round(decimals=4))
normed_levels = (levels - levels[0])/(levels[-1] - levels[0])
colors = [cm.Reds(normed_level) for normed_level in normed_levels]
fig2, ax2 = plt.subplots()
cont_cm = ax2.contour(Ea_mesh, log_A_mesh, Z, levels = levels, colors=colors)
ax2.set_xlabel('Ea (eV)')
ax2.set_ylabel('log(A)')
fig2.colorbar(cont_cm)
ax2.set_xlim([Ea_linspace[0],Ea_linspace[-1]])
ax2.set_ylim([log_A_linspace[0], log_A_linspace[-1]])
fig2.tight_layout()
fig2.savefig('onecomponentpriorcontour.png', dpi=220)
############### Now map to the observable rate.
fig3, ax3 = plt.subplots()
for i, log_A in enumerate(log_A_linspace):
    for j, Ea in enumerate(Ea_linspace):
        tpr_theta = odeint(onecomponent, theta_t0, t, args = (log_A, Ea))
        prob_obs = Z[i, j]
        k = np.argmin((np.abs(levels-prob_obs))) # gets the closes level
        rate = tpr_theta[:,0]*np.exp(-(Ea-kB*T*log_A)/(kB*T))
        ax3.plot(T, rate, color = colors[k]) # colors and levels match indices.
ax3.set_xlabel('T (K)')
ax3.set_ylabel(r'$rate (s^{-1})$')
fig3.colorbar(cont_cm)
fig3.tight_layout()
fig3.savefig('onecomponentratecontours.png', dpi=220)
