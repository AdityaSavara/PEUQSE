import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.integrate import odeint
import processing_function_Langmuir_CO_H2O_four_parameters as fun
global T
global pA
global pB
fun.T = 598.15
fun.pA = 0.15
fun.pB = 0.1
mean = [-0.687, 1.50e-3, -0.858, 1.50e-3]
cov_diag = np.power([5.00e-2, 2.07e-4, 1.25e-1, 2.07e-4],2)
cov = np.zeros((4,4))
row,col = np.diag_indices(cov.shape[0])
cov[row,col] = cov_diag
prior = np.random.multivariate_normal(mean,cov,5000)
theta_A_prior = []
for i in range(len(prior)):
    theta_A_prior_sample = fun.Langmuir_compete_ads(prior[i,:])
    theta_A_prior.append(theta_A_prior_sample)
fig,ax = plt.subplots(figsize=(4,4))
ax.hist(theta_A_prior,density=True,bins=500)
ax.set_xlabel(r'$log(\theta_A)$')
ax.set_ylabel(r'$probability density$')
fig.tight_layout()
fig.savefig('four_parameter_theta_A_forward.png',dpi=220)

synthetic_theta_A = fun.Langmuir_compete_ads([-0.687+5.00e-2,1.50e-3, -0.858, 1.50e-3])
print('synthetic_theta_A',synthetic_theta_A)

print('forward standard deviation ',np.std(theta_A_prior))

theta_A_prior_H_CO = []
for j in range(len(prior)):
    theta_A_prior_sample = fun.Langmuir_compete_ads([prior[j,0],1.50e-3, -0.858, 1.50e-3])
    theta_A_prior_H_CO.append(theta_A_prior_sample)

print('forward standard deviation H_CO only uncertain ',np.std(theta_A_prior_H_CO))

fig1,ax1 = plt.subplots(figsize=(4,4))
ax1.hist(theta_A_prior_H_CO,density=True,bins=500)
ax1.set_xlabel(r'$log(\theta_A)$')
ax1.set_ylabel(r'$probability density$')
fig1.tight_layout()
fig1.savefig('four_parameter_theta_A_forward_H_CO_only_uncertain.png',dpi=220)
