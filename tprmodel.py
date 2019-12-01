
import numpy as np

# The below function must return a vector of rates. 
def tprequation(tpr_theta,t,Ea_1, Ea_2, log_A1, log_A2, gamma1, gamma2,beta_dTdt,start_T): #beta_dTdT is the heating rate. 
    if tpr_theta.ndim == 1:  #for consistency, making tpr_theta a 2D array if it does not start as 2D. 
        tpr_theta2D = np.atleast_2d(tpr_theta)  
    if tpr_theta.ndim == 2: 
        tpr_theta2D = np.array(tpr_theta) 
    #Now find out how many species concentrations there are from the data: 
    num_of_concentrations = len(tpr_theta2D[0]) 
    Ea_Array = [Ea_1,Ea_2] 
    log_A_array = [log_A1, log_A2] 
    gamma_array = [gamma1, gamma2] 
    
    T = start_T + beta_dTdt*t  
    kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1  
    
    ratesList = [] 
    for rateIndex in range(num_of_concentrations): 
        rate = -tpr_theta2D[:,rateIndex]*np.exp(-(Ea_Array[rateIndex]-kB*T*log_A_array[rateIndex]-gamma_array[rateIndex]*tpr_theta2D[:,rateIndex])/(kB*T))  
    
        #Shortened below to one line (above) 
        # theta_i = tpr_theta2D[:,rateIndex] 
        # Ea_i = Ea_Array[rateIndex] 
        # log_A_i = log_A_array[rateIndex] 
        # gamma_i = gamma_array[rateIndex] 
        # rate = -theta_i*np.exp(-(Ea_i-kB*T*log_A_i-gamma_i*theta_i)/(kB*T))  
      
        #The above expression is the general form of this: rate_2 = -theta_2*np.exp(-(Ea_2-kB*T*log_A2-gamma2*theta_2)/(kB*T))       
        ratesList.append(rate) 
    
    if tpr_theta.ndim == 1: 
        ratesList = list(np.array(ratesList).flatten()) #for some reason, needs to be flattened for the MCMC. 
    return ratesList 
    