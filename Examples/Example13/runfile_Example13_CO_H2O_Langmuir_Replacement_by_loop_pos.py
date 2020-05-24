import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.integrate import odeint
import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    
    #observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    
    import processing_function_Langmuir_ads as fun

    UserInput.responses['responses_observed_uncertainties'] = 0.05
    
    UserInput.simulated_response_plot_settings['x_label'] = r'$time$'
    UserInput.simulated_response_plot_settings['y_label'] = r'$\theta_A$'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example12' #This creates the filename, also.
    UserInput.responses['independent_variables_values'] = [400 , 0.1]
    UserInput.responses['independent_variables_names'] = ['T(K)', 'P(bar)']
    fun.connected_variables_values = UserInput.responses['independent_variables_values'] #It is important to push the list *into* the other module.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'delta_G':r'$\deltaG (eV)$'}
    UserInput.model['InputParameterPriorValues'] = [-0.2] # eV
    UserInput.model['InputParametersPriorValuesUncertainties'] = [0.05] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D  This is a standard deviation!
    UserInput.model['InputParameterInitialGuess'] = [-0.2] #This is where the mcmc chain will start.
    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    
    #InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
    UserInput.model['simulateByInputParametersOnlyFunction'] = fun.Langmuir_replacement #This must simulate with *only* the parameters listed above, and no other arguments.
    #UserInput.model['simulationOutputProcessingFunction'] = processing_functions_tpd_odeint.no_log_wrapper_func #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "off"
    
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 100
    UserInput.parameter_estimation_settings['mcmc'] = True 
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = None #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10
    UserInput.parameter_estimation_settings['mcmc_length'] = 500
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.3
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    
    UserInput.parameter_pairs_for_contour_plots = [[0, 0]]
    UserInput.contour_settings_custom['contours_normalized'] = False
    UserInput.contour_settings_custom['figure_name'] = 'Mumpce_contour_plot_Langmuir_compete_ads'
    #After making the UserInput, now we make a 'parameter_estimation' object from it.

    UserInput.doe_settings['independent_variable_grid_center'] = [400, 0.2]
    UserInput.doe_settings['independent_variable_grid_interval_size'] = [100, 0.1]
    UserInput.doe_settings['independent_variable_grid_num_intervals'] = [2,2] #This is the number in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the parameter to vary.
    
#    PE_object = CKPQ.parameter_estimation(UserInput)
#    PE_object.designOfExperiments()
#    print("line 58, runfile", len(PE_object.info_gains_matrices_array))
#    print(PE_object.info_gains_matrices_array[0])
#    
#    
#    import matplotlib.pyplot as plt
#    
#    from mpl_toolkits.mplot3d import Axes3D  
#    
#    fig, ax0 =  plt.subplots(1)
#    ax0 = plt.axes(projection ='3d')
#    xValues = PE_object.info_gains_matrices_array[0][:,0]
#    yValues = PE_object.info_gains_matrices_array[0][:,1]
#    zValues = PE_object.info_gains_matrices_array[0][:,2]
#    print(xValues)
#    print(yValues)
#    print(zValues)
#    
#    im = ax0.plot_trisurf(xValues,yValues, zValues)
#    plt.show()
#    
#    sys.exit()
    print("######################Below this line is the old way but plots more nicely.############################")
    print("######################For this example, we also need to use the loop because the uncertainty in DeltaG is a function of temperature.###########################)
    global T
    global pA
    global pB
    kB = 8.61733035E-5 #eV/K\n
    Rgas = 8.314 #J mol-1 K-1 
    
    
    
    temperatures = np.linspace(298.15,498.15,5)
    pressures = np.linspace(0.1,0.3,5)


    ###Here was some preliminary synthetic data creation  that is not used anymore.######
#    theta_A_list=[]
#    for p in pressures:
#        for t in temperatures:
#            fun.pA = p
#            fun.T = t
#            theta_A = fun.Langmuir_compete_ads(delta_G_1)
#            theta_A_list.append(theta_A)
#
#    fig = plt.figure(figsize = (5,5))
#    ax = fig.add_subplot(111, projection='3d')
#    TT, pp = np.meshgrid(temperatures, pressures)
#    theta_A_array=np.asarray(theta_A_list)
#    TA = theta_A_array.reshape(TT.shape)
#    surf = ax.plot_surface(TT,pp,TA,cmap=matplotlib.cm.coolwarm)
#    ax.set_xlabel('Temperature (K)')
#    ax.set_ylabel(r'$p_A$')
#    ax.set_zlabel(r'$\theta_A$')
#    ax.set_title(r'Surface Plot of $\theta_A$')
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    fig.savefig('synthetic_observables_theta_A.png',dpi=220)
    
    deltaHprior = 16.5 * 1000 #This is in J/mol
    deltaSprior = 5.41 #This is in J/(mol*K)
    
    prior = np.random.normal(-0.2,0.05,10000) #The first number is the mean, second is uncertainty, last is sampling number.
    info_gains=[]
    info_gains_KL=[]
    PE_object_list = []
    for p in pressures:
        for t in temperatures:
            fun.pA = p
            fun.T = t
            deltaG = 
            theta_A_obs = fun.Langmuir_compete_ads(delta_G_1) # Remember delta_G_1 is offset: -0.15 instead of the prior -0.2.
            UserInput.responses['responses_observed'] = theta_A_obs
            PE_object_list.append(CKPQ.parameter_estimation(UserInput))
            [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = PE_object_list[-1].doMetropolisHastings()
            info_gains.append(info_gain)
            fig, ax = plt.subplots()
            (density0,bins0,pathces0)=ax.hist([prior,PE_object_list[-1].post_burn_in_samples.flatten()],bins=100,label=['prior','posterior'],density=True)
            KL = density0[1]*np.log(density0[1]/density0[0])
            KL = KL[np.isfinite(KL)]
            print('KL divergence',np.sum(KL))
            info_gains_KL.append(np.sum(KL))
            ax.legend()
            ax.set_ylabel('Probability density')
            ax.set_title(r'Prior and Posterior Density Plot at T = {} (K) $p_A$ = {}'.format(str(t),str(p)))
            fig.savefig('Lang_prior_post_plots/prior_and_posterior_histogram_T_{}_p_A_{}.png'.format(str(t),str(p)), dpi=300)
    fig,ax = plt.subplots(figsize=(5,5))
    #ax = fig.add_subplot(111, projection='3d')
    info_gains=np.asarray(info_gains)
    IG = info_gains.reshape(TT.shape)
    surf = ax.pcolor(TT,pp,IG,cmap=matplotlib.cm.coolwarm)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'$p_A$')
    #ax.set_zlabel('Information Gain')
    ax.set_xticks(temperatures)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_title('Information Gain Surface')
    fig.colorbar(surf, shrink=0.5, aspect=5) 
    fig.savefig('info_gain_surface_Langmuir_compete_ads.png', dpi=220)

    fig1,ax1 = plt.subplots(figsize=(5,5))
    #ax = fig.add_subplot(111, projection='3d')
    info_gains_KL=np.asarray(info_gains_KL)
    IG_KL = info_gains_KL.reshape(TT.shape)
    print('info_gains_KL',info_gains_KL)
    surf = ax1.pcolor(TT,pp,IG_KL,cmap=matplotlib.cm.coolwarm)
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel(r'$p_A$')
    #ax.set_zlabel('Information Gain')
    ax1.set_xticks(temperatures)
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_title('Information Gain Surface from Kullback Leibler Divergence')
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    fig1.savefig('info_gain_surface_Langmuir_compete_ads_Kullback_Leibler.png', dpi=220)


