import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.integrate import odeint
import sys; sys.path.insert(0, '../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

if __name__ == "__main__":
    
    observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    
    import processing_function_mem_reactor_km_only as fun

    UserInput.responses['responses_observed_uncertainties'] = 0.04
    
    UserInput.simulated_response_plot_settings['x_label'] = r'$reactor outlet$'
    UserInput.simulated_response_plot_settings['y_label'] = r'$ln(C_A)i$'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example12' #This creates the filename, also.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'delta_G':r'$\deltaG (eV)$'}
    UserInput.model['InputParameterPriorValues'] = [2e3] # eV
    UserInput.model['InputParametersPriorValuesUncertainties'] = [1e3] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D  This is a standard deviation!
    UserInput.model['InputParameterInitialGuess'] = [3e3] #This is where the mcmc chain will start.
    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    
    #InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
    UserInput.model['simulateByInputParametersOnlyFunction'] = fun.mem_reactor #This must simulate with *only* the parameters listed above, and no other arguments.
    #UserInput.model['simulationOutputProcessingFunction'] = processing_functions_tpd_odeint.no_log_wrapper_func #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "1e5"
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['mcmc_checkPointFrequency'] = 100
     
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10 #normally 1000 
    UserInput.parameter_estimation_settings['mcmc_length'] = 100 #normally 10000
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.05
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['mcmc_info_gain_cutoff'] = 1E-4
    
    UserInput.contour_plot_settings['parameter_pairs'] = [[0, 0]]
    UserInput.contour_plot_settings['contours_normalized'] = False
    UserInput.contour_plot_settings['figure_name'] = 'Mumpce_contour_plot_mem_reactor'
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    global T
    global volume
    global F0
    global P0
    global R
    global kB
    global k_1
    global k_minus_1
    kB = 8.61733035E-5 #eV/K\n
    k_1 = 1e2
    k_minus_1 = 1
    k_m = 3e3
    F0 = np.array([5, 0])
    temperatures = np.linspace(298.15,698.15,5)
    volumes = np.linspace(100,1100,5)

    flow_rates_a=[]
    for v in volumes:
        for t in temperatures:
            sol = odeint(fun.cmr, F0, np.linspace(0,v,50), args=(k_1,k_minus_1,k_m,t))
            conc_sol_last=sol[-1,0].T
            flow_rates_a.append(conc_sol_last)

    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(111, projection='3d')
    TT, V = np.meshgrid(temperatures, volumes)
    flow_rates_a=np.asarray(flow_rates_a)
    FRA = flow_rates_a.reshape(TT.shape)
    surf = ax.plot_surface(TT,V,FRA,cmap=matplotlib.cm.coolwarm)
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Volume')
    ax.set_zlabel('Flow Rate')
    ax.set_title('Surface Plot of F_A')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.savefig('synthetic_observables.png',dpi=220)

    fun.k_1 = 1e2
    fun.k_minus_1 = 1
    prior = np.random.normal(2e3,1e3,10000)
    info_gains=[]
    PE_object_list = []
    for v in volumes:
        for t in temperatures:
            sol = odeint(fun.cmr, F0, np.linspace(0,v,50), args=(k_1,k_minus_1,k_m,t))
            conc_sol_last=sol[-1,0].T
            print('conc_sol_last',conc_sol_last)
            UserInput.responses['responses_observed'] = conc_sol_last
            PE_object_list.append(PEUQSE.parameter_estimation(UserInput))
            fun.T = t
            fun.volume = v
            [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = PE_object_list[-1].doMetropolisHastings()
            info_gains.append(info_gain)
            fig, ax = plt.subplots()
            (density0,bins0,pathces0)=ax.hist([prior,PE_object_list[-1].post_burn_in_samples.flatten()],bins=100,label=['prior','posterior'],density=True)
            ax.legend()
            ax.set_ylabel('Probability density')
            ax.set_title('Prior and Posterior Density Plot at T = {} (K) volume = {} cm^3'.format(str(t),str(v)))
            fig.savefig('km_only_figures/prior_and_posterior_histogram_T_{}_V_{}.png'.format(str(t),str(v)), dpi=300)
    fig,ax = plt.subplots(figsize=(5,5))
    #ax = fig.add_subplot(111, projection='3d')
    T, V = np.meshgrid(temperatures, volumes)
    info_gains=np.asarray(info_gains)
    IG = info_gains.reshape(T.shape)
    surf = ax.pcolor(T,V,IG,cmap=matplotlib.cm.coolwarm)

    ax.set_xlabel('Temperature')
    ax.set_ylabel('Volume')
    #ax.set_zlabel('Information Gain')
    ax.set_xticks(temperatures)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_title('Information Gain Surface')
    fig.colorbar(surf, shrink=0.5, aspect=5) 
    fig.savefig('info_gain_surface_mem_reactor.png', dpi=220)
