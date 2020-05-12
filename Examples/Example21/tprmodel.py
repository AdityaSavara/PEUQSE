import numpy as np
global exportThetaVsEaGlobal
exportThetaVsEaGlobal = False

# The below function must return a vector of rates. 
def tprequation(tpr_theta,t,Ea_1, Ea_2, log_A1, log_A2, gamma1, gamma2,beta_dTdt,start_T): #beta_dTdT is the heating rate. 
    if tpr_theta.ndim == 1:  #for consistency, making tpr_theta a 2D array if it does not start as 2D. 
        tpr_theta2D = np.atleast_2d(tpr_theta)  
    if tpr_theta.ndim == 2: 
        tpr_theta2D = np.array(tpr_theta) 
    #Now find out how many species concentrations there are from the data: 
    num_of_concentrations = len(tpr_theta2D[0]) 
    numTimes = len(tpr_theta2D)
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
    
    
# The below function must return a vector of rates. 
def tprequationPiecewiseWithOffset(tpr_theta,t,Ea_1, log_A1, gamma1, verticalOffset, beta_dTdt,start_T, *Gamma_offsets, exportThetaVsEa=False): #beta_dTdT is the heating rate. 
    global exportThetaVsEaGlobal #This was added later than the original writing, in order to provide a convenient way to export exportThetaVsEa
    if exportThetaVsEaGlobal == True:
        exportThetaVsEa = True
    if tpr_theta.ndim == 1:  #for consistency, making tpr_theta a 2D array if it does not start as 2D. 
        tpr_theta2D = np.atleast_2d(tpr_theta)  
    if tpr_theta.ndim == 2: 
        tpr_theta2D = np.array(tpr_theta) 
    #Now find out how many species concentrations there are from the data: 
    num_of_concentrations = len(tpr_theta2D[0]) 
    numTimes = len(tpr_theta2D)
    coverageIntervals = np.linspace(0,1,num=len(Gamma_offsets), endpoint=True) #This will be used below.
    cumulative_gamma_offsets = np.cumsum(Gamma_offsets)
#    print(coverageIntervals)
#    print(Eaoffsets)
#    print(len(coverageIntervals))
#    print(len(Eaoffsets)) 
#    sys.exit()
    Ea_Array = [Ea_1] 
    log_A_array = [log_A1] 
    gamma_array = [gamma1]    
    T = start_T + beta_dTdt*t  
    kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1  
    ratesList = [] 
    for rateIndex in range(num_of_concentrations): 
        if numTimes >= 1:
            ratesAtEachTime = []
            ThetAndEaForEachTime = []
            #print(tpr_theta2D)
            for timeIndex in range(numTimes):
                thisTheta = tpr_theta2D[timeIndex][rateIndex]
                baseEaForThisTheta = Ea_Array[rateIndex]
                gamma_dependance = gamma_array[rateIndex] + cumulative_gamma_offsets  
                this_gamma = np.interp(thisTheta, coverageIntervals, gamma_dependance)
                Ea_modified = baseEaForThisTheta-this_gamma*thisTheta*baseEaForThisTheta #+ Ea_1_Offest_for_this_theta
                #For the Temperature, either we have one temperature or many. Check if it is iterable.
                from typing import Iterable
                if isinstance(T, Iterable):
                    this_T = T[timeIndex]
                else:
                    this_T = T
                
                #print("line 83", np.shape(thisTheta), np.shape(Ea_modified), np.shape(this_T), np.shape(log_A_array[rateIndex]), np.shape()
                rate = -thisTheta*np.exp(-(Ea_modified-kB*this_T*log_A_array[rateIndex])/(kB*this_T))  - verticalOffset
                thisThetAndEa = [thisTheta,Ea_modified, this_gamma]
                ThetAndEaForEachTime.append(thisThetAndEa)
                ratesAtEachTime.append(np.array(rate))
                #this is appendin once for each time.
        
            #Shortened below to one line (above) 
            # theta_i = tpr_theta2D[:,rateIndex] 
            # Ea_i = Ea_Array[rateIndex] 
            # log_A_i = log_A_array[rateIndex] 
            # gamma_i = gamma_array[rateIndex] 
            # rate = -theta_i*np.exp(-(Ea_i-kB*T*log_A_i-gamma_i*theta_i)/(kB*T))  
          
            #The above expression is the general form of this: rate_2 = -theta_2*np.exp(-(Ea_2-kB*T*log_A2-gamma2*theta_2)/(kB*T))       
        ratesList.append(ratesAtEachTime) #This is appending once for each species..
        if exportThetaVsEa==True:
            np.savetxt("ThetaVersusEaCurve.csv", ThetAndEaForEachTime, delimiter=",")
    if tpr_theta.ndim == 1: 
        ratesList = list(np.array(ratesList).flatten()) #for some reason, needs to be flattened for the MCMC. 
    return ratesList 
        
# The below function must return a vector of rates. 
def tprequationPiecewise(tpr_theta,t,Ea_1, log_A1, gamma1, beta_dTdt,start_T, *Gamma_offsets, exportThetaVsEa=False): #beta_dTdT is the heating rate. 
    global exportThetaVsEaGlobal #This was added later than the original writing, in order to provide a convenient way to export exportThetaVsEa
    if exportThetaVsEaGlobal == True:
        exportThetaVsEa = True
    if tpr_theta.ndim == 1:  #for consistency, making tpr_theta a 2D array if it does not start as 2D. 
        tpr_theta2D = np.atleast_2d(tpr_theta)  
    if tpr_theta.ndim == 2: 
        tpr_theta2D = np.array(tpr_theta) 
    #Now find out how many species concentrations there are from the data: 
    num_of_concentrations = len(tpr_theta2D[0]) 
    numTimes = len(tpr_theta2D)
    coverageIntervals = np.linspace(0,1,num=len(Gamma_offsets), endpoint=True) #This will be used below.
    cumulative_gamma_offsets = np.cumsum(Gamma_offsets)
#    print(Eaoffsets)
#    print(len(coverageIntervals))
#    print(len(Eaoffsets)) 
#    sys.exit()
    Ea_Array = [Ea_1] 
    log_A_array = [log_A1] 
    gamma_array = [gamma1]    
    T = start_T + beta_dTdt*t  
    kB = 1.380649e-26*6.0221409e+23 #kJ mol^-1 K^-1  
    ratesList = [] 
    for rateIndex in range(num_of_concentrations): 
        if numTimes >= 1:
            ratesAtEachTime = []
            ThetAndEaForEachTime = []
            #print(tpr_theta2D)
            for timeIndex in range(numTimes):
                thisTheta = tpr_theta2D[timeIndex][rateIndex]
                baseEaForThisTheta = Ea_Array[rateIndex]
                gamma_dependance = gamma_array[rateIndex] + cumulative_gamma_offsets  #In future, this might become cumulative_gamma_offsets[rateIndex] but right now it is designed for a single case.  In essence, 
                this_gamma = np.interp(thisTheta, coverageIntervals, gamma_dependance)
                Ea_modified = baseEaForThisTheta-this_gamma*thisTheta*baseEaForThisTheta #+ Ea_1_Offest_for_this_theta
                #For the Temperature, either we have one temperature or many. Check if it is iterable.
                from typing import Iterable
                if isinstance(T, Iterable):
                    this_T = T[timeIndex]
                else:
                    this_T = T
                
                #print("line 83", np.shape(thisTheta), np.shape(Ea_modified), np.shape(this_T), np.shape(log_A_array[rateIndex]), np.shape()
                rate = -thisTheta*np.exp(-(Ea_modified-kB*this_T*log_A_array[rateIndex])/(kB*this_T)) 
                thisThetAndEa = [thisTheta,Ea_modified, this_gamma]
                ThetAndEaForEachTime.append(thisThetAndEa)
                ratesAtEachTime.append(np.array(rate))
                #this is appendin once for each time.
        
            #Shortened below to one line (above) 
            # theta_i = tpr_theta2D[:,rateIndex] 
            # Ea_i = Ea_Array[rateIndex] 
            # log_A_i = log_A_array[rateIndex] 
            # gamma_i = gamma_array[rateIndex] 
            # rate = -theta_i*np.exp(-(Ea_i-kB*T*log_A_i-gamma_i*theta_i)/(kB*T))  
          
            #The above expression is the general form of this: rate_2 = -theta_2*np.exp(-(Ea_2-kB*T*log_A2-gamma2*theta_2)/(kB*T))       
        ratesList.append(ratesAtEachTime) #This is appending once for each species..
        if exportThetaVsEa==True:
            np.savetxt("ThetaVersusEaCurve.csv", ThetAndEaForEachTime, delimiter=",")
    if tpr_theta.ndim == 1: 
        ratesList = list(np.array(ratesList).flatten()) #for some reason, needs to be flattened for the MCMC. 
    return ratesList 
        