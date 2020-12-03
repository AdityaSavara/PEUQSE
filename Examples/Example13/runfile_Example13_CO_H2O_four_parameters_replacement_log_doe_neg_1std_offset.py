import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.integrate import odeint
#import dill
import sys; sys.path.append('../../');  import CheKiPEUQ as CKPQ
import CheKiPEUQ.UserInput as UserInput

if __name__ == "__main__":
    
    #observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
    
    import processing_function_Langmuir_CO_H2O_four_parameters as fun

    UserInput.responses['responses_abscissa'] = [fun.T-25,fun.T,fun.T+25]
    UserInput.responses['responses_observed'] = [0.00, 0.00, 0.00] #The initial values won't be used during the DOE.
    UserInput.responses['responses_observed_uncertainties'] = [np.log10(1.5),np.log10(1.5),np.log10(1.5)]
    
    UserInput.simulated_response_plot_settings['x_label'] = r'$Temperature (K)$'
    UserInput.simulated_response_plot_settings['y_label'] = r'$\theta_A$'
    #UserInput.simulated_response_plot_settings['y_range'] = [0.00, 0.025] #optional.
    UserInput.simulated_response_plot_settings['figure_name'] = 'Posterior_Example12' #This creates the filename, also.

    UserInput.model['parameterNamesAndMathTypeExpressionsDict'] = {'delta_H_rxn':'delta_H_rxn', 'delta_S_rxn':'delta_S_rxn'}
    UserInput.model['InputParameterPriorValues'] = [-0.687, 1.50e-3, -0.858, 1.50e-3] # [H_CO, S_CO, H_H2O, S_H2O] eV
    #Below uses numbers that Eric converted for eV units.
    delta_H_rxn = -0.687 - (-0.858)
    delta_S_rxn = -1.50e-3 - (-1.45E-3)
    delta_H_rxn_uncertainty = (5.00e-2**2 + 1.25e-1**2)**0.5
    delta_S_rxn_uncertainty = (2.07e-4**2 + 2.07e-4**2)**0.5
    UserInput.model['InputParameterPriorValues'] = [delta_H_rxn, delta_S_rxn]                                                                               
    UserInput.model['InputParametersPriorValuesUncertainties'] = [delta_H_rxn_uncertainty, delta_S_rxn_uncertainty] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D  This is a standard deviation!
#    UserInput.model['InputParameterInitialGuess'] = [-0.687, 1.50e-3, -0.858, 1.50e-3] #This is where the mcmc chain will start.
    #InputParameterInitialValues = [41.5, 41.5, 13.0, 13.0, 0.1, 0.1] # Ea1_mean, Ea2_mean, log_A1_mean, log_A2_mean, gamma_1_mean, gamma_2_mean 
    #InputParametersInitialValuesUncertainties = [200, 200, 13, 13, 0.1, 0.1] #If user wants to use a prior with covariance, then this must be a 2D array/ list. To assume no covariance, a 1D array can be used.
    UserInput.model['simulateByInputParametersOnlyFunction'] = fun.Langmuir_replacement_three_temperatures_log #This must simulate with *only* the parameters listed above, and no other arguments.
    #UserInput.model['simulationOutputProcessingFunction'] = processing_functions_tpd_odeint.no_log_wrapper_func #Optional: a function to process what comes out of the simulation Function and then return an observable vector.
    UserInput.parameter_estimation_settings['scaling_uncertainties_type'] = "std"
    
    UserInput.parameter_estimation_settings['verbose'] = False 
    UserInput.parameter_estimation_settings['checkPointFrequency'] = 100
     
    UserInput.parameter_estimation_settings['mcmc_mode'] = 'unbiased'
    UserInput.parameter_estimation_settings['mcmc_random_seed'] = 0 #Normally set to None so that mcmc is set to be random. To get the same results repeatedly, such as for testing purposes, set the random seed to 0 or another integer for testing purposes.
    UserInput.parameter_estimation_settings['mcmc_burn_in'] = 10000
    UserInput.parameter_estimation_settings['mcmc_length'] = 200000
    UserInput.parameter_estimation_settings['mcmc_relative_step_length'] = 0.05
    UserInput.parameter_estimation_settings['mcmc_modulate_accept_probability']  = 0 #Default value of 0. Changing this value sharpens or flattens the posterior. A value greater than 1 flattens the posterior by accepting low values more often. It can be useful when greater sampling is more important than accuracy. One way of using this feature is to try with a value of 0, then with the value equal to the number of priors for comparison, and then to gradually decrease this number as low as is useful (to minimize distortion of the result). A downside of changing changing this variable to greater than 1 is that it slows the the ascent to the maximum of the prior, so there is a balance in using it. In contrast, numbers increasingly less than one (such as 0.90 or 0.10) will speed up the ascent to the maximum of the posterior, but will also result in fewer points being retained.
    UserInput.parameter_estimation_settings['mcmc_info_gain_returned'] = 'log_ratio' #obtains the information gain using the Kullback-Leibler divergence    
   
    UserInput.parameter_pairs_for_contour_plots = [[0, 1]]
    UserInput.contour_plot_settings['contours_normalized'] = False
    UserInput.contour_settings_custom['figure_name'] = 'Mumpce_contour_plot_Langmuir_compete_ads'
    #After making the UserInput, now we make a 'parameter_estimation' object from it.
    
    #It's good to run a test before doing a design of experiments.
#    PE_object = CKPQ.parameter_estimation(UserInput)
#    PE_object.doMetropolisHastings()
    #PE_object.createAllPlots()
    
    temperatures = np.linspace(398.15,598.15,5)
    pressures = np.linspace(0.20,1.20, 5)

#   This block of code was made for some diagnostic in an earlier version.
#    theta_A_list=[]
#    for p in pressures:
#        for t in temperatures:
#            fun.pA = p
#            fun.T = t
#            theta_A = fun.Langmuir_compete_ads(parameter_set)
#            theta_A_list.append(theta_A)
#
#    fig = plt.figure(figsize = (5,5))
#    ax = fig.add_subplot(111, projection='3d')
#    TT, pp = np.meshgrid(temperatures, pressures)
#    theta_A_array=np.asarray(theta_A_list)
#    TA = theta_A_array.reshape(TT.shape)
#    surf = ax.plot_surface(TT,pp,TA,cmap=matplotlib.cm.coolwarm)
#    ax.set_xlabel('Temperature (K)')
#    ax.set_ylabel(r'$p_{CO}$')
#    ax.set_zlabel(r'$\theta_A$')
#    ax.set_title(r'Surface Plot of $\theta_A$')
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    fig.savefig('synthetic_observables_theta_A.png',dpi=220)

    pos_or_neg = 'neg' #Put 'pos' or 'neg' here.
    loadFromPickleName = pos_or_neg
    loadFromPickle = False
    
    prior = np.random.normal(delta_H_rxn, delta_H_rxn_uncertainty,2000000)                                                                                                                                        
    info_gains=[]
    info_gains_KL =[]
    PE_object_list = []
    if loadFromPickle == True:
        import dill
        PE_object_list = dill.load(open("PE_object_list_dill_log10_"+loadFromPickleName+".dill", 'rb'))
        info_gains = dill.load(open("info_gains_list_dill_log10_"+loadFromPickleName+".dill", 'rb'))
    counter = 0
    for p in pressures:
        for t in temperatures:
            fun.pA = p
            fun.T = t
            
            if loadFromPickle == False:
                #In this example, we use UserInput.model['InputParameterPriorValues'] with a positive offset for the deltaHrxn by 1 standard deviation.
                ParameterValues_synth = np.array(UserInput.model['InputParameterPriorValues'])*1.0 #Just initializing.
                if pos_or_neg =='neg':
                    ParameterValues_synth[0] = UserInput.model['InputParameterPriorValues'][0] - UserInput.model['InputParametersPriorValuesUncertainties'][0] #Taking deltaHrxn plus 1 standard deviation for it. 
                if pos_or_neg == 'pos':
                    ParameterValues_synth[0] = UserInput.model['InputParameterPriorValues'][0] + UserInput.model['InputParametersPriorValuesUncertainties'][0] #Taking deltaHrxn plus 1 standard deviation for it. 
                theta_A_obs_synth = fun.Langmuir_replacement_three_temperatures_log(ParameterValues_synth) 
                print('line 101', theta_A_obs_synth)
                UserInput.responses['responses_observed'] = theta_A_obs_synth
                #UserInput.responses['responses_observed_uncertainties'] = [0.05, 0.05, 0.05]
                PE_object_list.append(CKPQ.parameter_estimation(UserInput))
                [map_parameter_set, muap_parameter_set, stdap_parameter_set, evidence, info_gain, samples, logP] = PE_object_list[counter].doMetropolisHastings()
                info_gains.append(info_gain)
                
            fig, ax = plt.subplots()
            # for G_CO: (density0,bins0,pathces0)=ax.hist([prior,PE_object_list[-1].post_burn_in_samples[:,0].flatten() - t*PE_object_list[-1].post_burn_in_samples[:,1].flatten()],bins=100,label=['prior','posterior'],density=True)
            (density0,bins0,pathces0)=ax.hist([prior,PE_object_list[counter].post_burn_in_samples[:,0].flatten()],bins=100,label=['prior','posterior'],density=True)
            KL = density0[1]*np.log(density0[1]/density0[0])
            KL = KL[np.isfinite(KL)]
            print('KL divergence',np.sum(KL))
            info_gains_KL.append(np.sum(KL))

            ax.legend()
            ax.set_xlabel(r'$\Delta H_{rxn} (eV)$')
            ax.set_ylabel('Probability density')
            ax.set_title(r'Prior and Posterior Density Plot at T = {} (K) $p_A$ = {}'.format(str(t),str(p)))
            fig.tight_layout()
            fig.savefig('Lang_prior_post_plots/prior_and_posterior_histogram_T_{}_p_A_{}.png'.format(str(t),str(p)), dpi=300)
            counter = counter + 1         

    if loadFromPickle == False:
        try:
            import dill
            dillPickleFile = open("PE_object_list_dill_log10_"+loadFromPickleName+".dill", 'wb')
            dill.dump(PE_object_list, dillPickleFile)
            dillPickleFile.close() #It is important to remember to close! Otherwise there can be problems when trying to open later.
            dillPickleFile = open("info_gains_list_dill_log10_"+loadFromPickleName+".dill", 'wb')
            dill.dump(info_gains, dillPickleFile)
            dillPickleFile.close() #It is important to remember to close! Otherwise there can be problems when trying to open later.            
        except: #Things can crash during saving, so the except statement should try to open and close file again.
            dillPickleFile = open("PE_object_list_dill_log10_"+loadFromPickleName+".dill", 'wb')
            dillPickleFile.close()
            dillPickleFile = open("info_gains_list_dill_log10_"+loadFromPickleName+".dill", 'wb')
            dillPickleFile.close() #It is important to remember to close! Otherwise there can be problems when trying to open later.
    
    fig,ax = plt.subplots(figsize=(5,5))
    #ax = fig.add_subplot(111, projection='3d')
    info_gains=np.asarray(info_gains)
    TT, pp = np.meshgrid(temperatures, pressures)
    IG = info_gains.reshape(TT.shape)
    surf = ax.pcolor(TT,pp,IG,cmap=matplotlib.cm.coolwarm)
    ax.set_xlabel(r'$temperature (K)$')
    ax.set_ylabel(r'$p_{CO}$')
    #ax.set_zlabel('Information Gain')
    ax.set_xticks(temperatures)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_title('Information Gain Surface')
    fig.colorbar(surf, shrink=0.5, aspect=5) 
    fig.tight_layout()
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

