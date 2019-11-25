import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy.stats import multivariate_normal
from scipy.integrate import odeint
import pandas as pd
from tprmodel import tprequation
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
verbose = True


class ip:
    #Ip for 'inverse problem'. Initialize prior chain starting point, chain burn-in length and total length, and Q (for proposal samples).  Initialize experimental data.  Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
    def __init__(self, UserInput = UserInput, verbose = True):
        if verbose: print("Bayes Model Initialized")
        self.mcmc_length = UserInput.mcmc_length
        self.mcmc_burn_in = UserInput.mcmc_burn_in # Number of samples trimmed off the beginning of the Markov chain.
        self.mu_prior = UserInput.mu_prior
        self.start_T = UserInput.T_0
        self.cov_prior = UserInput.cov_prior
        
    def import_experimental_settings(self):
        experiments_df = pd.read_csv('ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedLargerErrors.csv')
        # self.dT = experiments_df['dT'][0] # assuming dT and dt are constant throughout
        # self.dt = experiments_df['dt'][0]
        self.times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
        self.experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/1000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
        self.errors = np.array(experiments_df['Errors']) #.to_numpy()
        self.Q_mu = np.array([0.,0.,0.,0.,0.,0.]) # Q samples the next step at any point in the chain.  The next step may be accepted or rejected.  Q_mu is centered (0) around the current theta.
        self.Q_cov = self.cov_prior/20 # Take small steps.
    
    #main function to get samples
    def MetropolisHastings(self):
        self.import_experimental_settings()
        samples = np.zeros((self.mcmc_length,len(self.mu_prior)))
        samples[0,:] = self.mu_prior # Initialize the chain. Theta is initialized as the starting point of the chain.  It is placed at the prior mean.
        likelihoods_vec = np.zeros((self.mcmc_length,1))
        posteriors_un_normed_vec = np.zeros((self.mcmc_length,1))
        priors_vec = np.zeros((self.mcmc_length,1))
        for i in range(1,self.mcmc_length):
            if verbose: print(i)
            proposal_sample = samples[i-1,:] + np.random.multivariate_normal(self.Q_mu,self.Q_cov)
            prior_proposal = self.prior(proposal_sample)
            likelihood_proposal = self.likelihood(proposal_sample)
            prior_current_location = self.prior(samples[i-1,:])
            likelihood_current_location = self.likelihood(samples[i-1,:])
            accept_pro = (likelihood_proposal*prior_proposal)/(likelihood_current_location*prior_current_location) ###QUESTION: Is "pro" for probability of acceptance?
            if accept_pro> np.random.uniform():  #TODO: keep a log of the accept and reject. If the reject ratio is >90% or some other such number, warn the user.
                samples[i,:] = proposal_sample
                posteriors_un_normed_vec[i] = likelihood_proposal*prior_proposal
                likelihoods_vec[i] = likelihood_proposal
                priors_vec[i] = prior_proposal
            else:
                samples[i,:] = samples[i-1,:]
                posteriors_un_normed_vec[i] = likelihood_current_location*prior_current_location
                likelihoods_vec[i] = likelihood_current_location
                priors_vec[i] = prior_current_location
            ########################################
        samples = samples[self.mcmc_burn_in:]
        posteriors_un_normed_vec = posteriors_un_normed_vec[self.mcmc_burn_in:]
        likelihoods_vec = likelihoods_vec[self.mcmc_burn_in:]
        priors_vec = priors_vec[self.mcmc_burn_in:]
        # posterior probabilites are transformed to a standard normal (std=1) for obtaining the evidence:
        evidence = np.mean(posteriors_un_normed_vec)*np.sqrt(2*np.pi*np.std(samples)**2)
        posteriors_vec = posteriors_un_normed_vec/evidence
        log_ratios = np.log(posteriors_vec/priors_vec)
        log_ratios[np.isinf(log_ratios)] = 0
        log_ratios = np.nan_to_num(log_ratios)
        info_gain = np.mean(log_ratios)
        return [evidence, info_gain, samples]
    def prior(self,sample):
        probability = multivariate_normal.pdf(x=sample,mean=self.mu_prior,cov=self.cov_prior)
        return probability
    def likelihood(self,sample):
        tpr_theta = odeint(tprequation, [0.5, 0.5], self.times, args = (sample[0], sample[1], sample[2], sample[3], sample[4], sample[5],self.dT,self.dt,self.start_T)) # [0.5, 0.5] are the initial theta's.
        rate = tprequation(tpr_theta, self.times, sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], self.dT,self.dt,self.start_T)
        rate_tot = -np.sum(rate, axis=0)
        #intermediate_metric = np.mean(np.square(rate_tot - self.experiment) / np.square(self.errors ))
        probability_metric = multivariate_normal.pdf(x=rate_tot,mean=self.experiment,cov=self.errors)
        if verbose: print('likelihood probability',probability_metric)
        return probability_metric
        
    
ip_object = ip()
[evidence, info_gain, samples] = ip_object.MetropolisHastings()
############################################# The computation portion is contained above.
post_mean = np.mean(samples, axis=0)
experiments_df = pd.read_csv('ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedLargerErrors.csv')
dT = experiments_df['dT'][0] # assuming dT and dt are constant throughout
dt = experiments_df['dt'][0]
start_T = experiments_df['AcH - T'][0]
times = np.array(experiments_df['time']) #experiments_df['time'].to_numpy() #The to_numpy() syntax was not working for Ashi.
experiment = np.array(experiments_df['AcHBackgroundSubtracted'])/1000  #experiments_df['AcHBackgroundSubtracted'].to_numpy()/1000
errors = np.array(experiments_df['Errors']) #.to_numpy()
tpr_theta = odeint(tprequation, [0.5, 0.5], times, args = (post_mean[0], post_mean[1], post_mean[2], post_mean[3], post_mean[4], post_mean[5],dT,dt,start_T)) # [0.5, 0.5] are the initial theta's.
rate = tprequation(tpr_theta, times, post_mean[0], post_mean[1], post_mean[2], post_mean[3], post_mean[4], post_mean[5], dT,dt,start_T)
rate_tot = -np.sum(rate, axis=0)

fig1, ax1 = plt.subplots()


ax1.plot(np.array(experiments_df['AcH - T']),rate_tot, 'r')
ax1.plot(np.array(experiments_df['AcH - T']),np.array(experiments_df['AcHBackgroundSubtracted'])/1000,'g')
ax1.set_xlabel('T (K)')
ax1.set_ylabel(r'$rate (s^{-1})$')
ax1.legend(['model posterior', 'experiments'])
fig1.tight_layout()
fig1.savefig('tprposterior.png', dpi=220)

fig2, ax2 = plt.subplots()
ax2.hist(samples[:,0])
ax2.set_ylabel('frequency')
ax2.set_xlabel(r'$E_{a1}$')
fig2.tight_layout()
fig2.savefig('Ea1.png', dpi=220)

fig3, ax3 = plt.subplots()
ax3.hist(samples[:,1])
ax3.set_ylabel('frequency')
ax3.set_xlabel(r'$E_{a2}$')
fig3.tight_layout()
fig3.savefig('Ea2.png', dpi=220)

fig4, ax4 = plt.subplots()
ax4.hist(samples[:,2])
ax4.set_ylabel('frequency')
ax4.set_xlabel(r'$log(A_{1})$')
fig4.tight_layout()
fig4.savefig('logA1.png', dpi=220)

fig5, ax5 = plt.subplots()
ax5.hist(samples[:,3])
ax5.set_ylabel('frequency')
ax5.set_xlabel(r'$log(A_{2})$')
fig5.tight_layout()
fig5.savefig('logA2.png', dpi=220)

fig6, ax6 = plt.subplots()
ax6.hist(samples[:,4])
ax6.set_ylabel('frequency')
ax6.set_xlabel(r'$\gamma_{1}$')
fig6.tight_layout()
fig6.savefig('gamma1.png', dpi=220)

fig7, ax7 = plt.subplots()
ax7.hist(samples[:,5])
ax7.set_ylabel('frequency')
ax7.set_xlabel(r'$\gamma_{2}$')
fig7.tight_layout()
fig7.savefig('gamma2.png', dpi=220)
