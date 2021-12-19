import cantera as ct
import cantera.ck2cti as ck2cti
import numpy as np
from PEUQSE.simulationDriver import canteraSimulate
from PEUQSE.simulationDriver import canteraKineticsParametersParser
import copy

#First made a 'test' in main, then made a simulation wrapper to do the same thing.

#We will define some things *outside* of the function to load the model.
model_location = ".\\PEUQSE\\simulationDriver\\"
model_name = "ceO2"
import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a parameter optimization to reduce I/O.

cti_top_info_filename =model_location+ model_name + "_cti_top_info.cti"
with open(cti_top_info_filename, "r") as cti_top_info_file:
    cti_top_info_string = cti_top_info_file.read()     

initial_reactions_parameters_array = np.genfromtxt(model_location+ model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
#ow, looking at the paramters array file, we can see where the different parameters positions are that we can modify.
print(initial_reactions_parameters_array)

global observed_x_values #This can be set initially, then changed later if desired.
import processing_functions_tpd_odeint
observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
observed_x_values, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint.import_experimental_settings(observed_data_Filename)

global firstSimulation
firstSimulation = True


def cantera_simulation_wrapper_example2(parametersArray): #This takes in *only* the adjustable parameters.
    #Now we modify the 'initial' parameters array accordingly.
    #In this example, we're modifying reactions 2 and 4 which are indices 1 and 3.
    modified_reactions_parameters_array = copy.deepcopy(initial_reactions_parameters_array)
    
    #We now have... E_a1, E_a2, log_A1, log_A2, gamma_1, and gamma_2.
#    E_a1 = parametersArray[0] #We are receiving in kJ/mol
#    E_a2 =parametersArray[1] #We are receiving in kJ/mol
#    log_A1 =parametersArray[2]
#    log_A2 =parametersArray[3]
#    gamma_1 =parametersArray[4] <-- this becomes small e in cantera.
#    gamma_2 = parametersArray[5] <-- this becomes small e in cantera.
    modified_reactions_parameters_array[1][5] = parametersArray[0]*1000 #convert from kJ/mol to J/mol
    modified_reactions_parameters_array[3][5] = parametersArray[1]*1000 #convert from kJ/mol to J/mol
    modified_reactions_parameters_array[1][3] = 10**parametersArray[2]
    modified_reactions_parameters_array[3][3] = 10**parametersArray[3]    
    modified_reactions_parameters_array[1][8] = 10**parametersArray[4]
    modified_reactions_parameters_array[3][8] = 10**parametersArray[5]        
    
    #Now we do the simulation. The first time, we need to do a full simulation flow.  But second time and later we can just modify the cantera phases object and run again.
    global firstSimulation
    global canteraPhases #need to make it global so the object sticks around for next function call.
    if firstSimulation == True:
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.create_cti_and_SimulatePFRorTPRwithCantera(model_name, modified_reactions_parameters_array, ceO2_input_simulation_settings, cti_top_info_string = cti_top_info_string)
        firstSimulation = False
    elif firstSimulation == False: #This must be an elif, otherwise it will always be executed. In this function, the cantera phases object will be created the first time the function is called. Then will exist for later.
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.modify_reactions_and_SimulatePFRorTPRwithCantera(model_name, modified_reactions_parameters_array, ceO2_input_simulation_settings, canteraPhases=canteraPhases)        
    
    #Now we parse the output.
    times = cantera_phase_rates['surf'][:,1] #This is in seconds.
    temperatures = ceO2_input_simulation_settings.heating_rate*times+ceO2_input_simulation_settings.T_surf
    totalAbsoluteRate = -1*(cantera_phase_rates['surf'][:,3]+cantera_phase_rates['surf'][:,4]) 

    #Now we convert units.
    totalAbsoluteRate = totalAbsoluteRate * 1000 #going from kmol/m^2/s to mol/m^2/s
    totalAbsoluteRate = totalAbsoluteRate * (1.0/(100.0**2)) # 1/m^2 * (1^2 m^2)/(100^2 cm^2)  gives mol/cm^2/s which is the units we want to compare to site density
    totalAbsoluteRate = totalAbsoluteRate / 2.72E-9 #now should be moles produced per moles of site, which is molecules per site per second.
    
    #Now we need to interpolate to the x_values that we need.
    #we use a global called x_values and assume that our output is a continuous function.
    global observed_x_values
    interpolatedRate = np.interp(observed_x_values, times, totalAbsoluteRate)
    
    
    return interpolatedRate

totalAbsoluteRate = cantera_simulation_wrapper_example2([61.5, 41.5, 13.0, 13.0, 0.1, 0.1])

from matplotlib import pyplot as plt
import numpy as np
x = observed_x_values
y = totalAbsoluteRate
plt.plot(x,y)
plt.xlabel("time (s)")
plt.ylabel("total rate")
plt.title('Acetaldehyde Desorption Rate')
plt.show()    

totalAbsoluteRate = cantera_simulation_wrapper_example2([51.5, 51.5, 13.0, 13.0, 0.1, 0.1])

from matplotlib import pyplot as plt
import numpy as np
x = observed_x_values
y = totalAbsoluteRate
plt.plot(x,y)
plt.xlabel("time (s)")
plt.ylabel("total rate")
plt.title('Acetaldehyde Desorption Rate')
plt.show()    
if __name__ == "__main__":
    model_location = ".\\simulationDriver\\"
    model_name = "ceO2"
    
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
    ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a parameter optimization to reduce I/O.
    reactions_parameters_array = np.genfromtxt(model_location+ model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
#    print("Initial Reactions Parameters Array")
#    print(reactions_parameters_array)
#    
    
    cti_top_info_filename =model_location+ model_name + "_cti_top_info.cti"
    with open(cti_top_info_filename, "r") as cti_top_info_file:
        cti_top_info_string = cti_top_info_file.read()     
    
    modified_reactions_parameters_array = reactions_parameters_array
    concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
    canteraSimulate.create_cti_and_SimulatePFRorTPRwithCantera("ceO2", modified_reactions_parameters_array, ceO2_input_simulation_settings, cti_top_info_string = cti_top_info_string)
    
    #Playing with the variables we see we can get what we want from:
#    print(cantera_phase_rates['surf'])
    times = cantera_phase_rates['surf'][:,1] #This is in seconds.
    temperatures = ceO2_input_simulation_settings.heating_rate*times+ceO2_input_simulation_settings.T_surf
#    print(np.shape(cantera_phase_rates['surf']))
    totalAbsoluteRate = -1*(cantera_phase_rates['surf'][:,3]+cantera_phase_rates['surf'][:,4]) 
    
    
    #To convert totalAbsoluteCoverages to relative coverages, we recognize that cantera uses rates of kmol/m^2/s
    #See this link for more information on cantera units: https://github.com/Cantera/cantera/issues/760
    #Surface site density: 2.72E-9 mol/cm**2 , so now we need to convert...
    #Let's make some conversion factors for convenience....
    totalAbsoluteRate = totalAbsoluteRate * 1000 #going from kmol/m^2/s to mol/m^2/s
    totalAbsoluteRate = totalAbsoluteRate * (1.0/(100.0**2)) # 1/m^2 * (1^2 m^2)/(100^2 cm^2)  gives mol/cm^2/s which is the units we want to compare to site density
    totalAbsoluteRate = totalAbsoluteRate / 2.72E-9 #now should be moles produced per moles of site, which is molecules per site per second.
#    print(totalAbsoluteRate)
#    print(np.shape(totalAbsoluteRate))
    
    from matplotlib import pyplot as plt
    import numpy as np
    x = temperatures
    y = totalAbsoluteRate
    plt.plot(x,y)
    plt.xlabel("temperature")
    plt.ylabel("total rate")
    plt.title('Acetaldehyde Desorption Rate')
    plt.show()    
        
        
        
    #np.savetxt("ceO2_output_rates_all_"+"FullCTIsamplingCase.csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)              
        