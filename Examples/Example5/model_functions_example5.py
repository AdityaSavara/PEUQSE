import cantera as ct
import cantera.ck2cti as ck2cti
import numpy as np
from CheKiPEUQ.simulationDriver import canteraSimulate
from CheKiPEUQ.simulationDriver import canteraKineticsParametersParser
import copy

#First made a 'test' in main, then made a simulation wrapper to do the same thing.

#We will define some things *outside* of the function to load the model.
model_location = ".\\CheKiPEUQ\\simulationDriver\\"
model_name = "ceO2"
import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a parameter optimization to reduce I/O.
ceO2_input_simulation_settings.piecewise_coverage_dependence = True
reactions_parameters_array = np.genfromtxt(model_location+ model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
##ow, looking at the paramters array file, we can see where the different parameters positions are that we can modify.
print(reactions_parameters_array)

#This is before any offset modification:
ceO2_input_simulation_settings.original_reactions_parameters_array = reactions_parameters_array
piecewise_coverage_intervals = np.array([0,0.1,0.2,0.3,0.4,0.5,1.0])
species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
kineticParameterName = "A"
modifiers_A = ([0,0,0,0,0,0,0], [-1,1,-1,-1,-1,-1,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0])
canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
kineticParameterName = "E" 
modifiers_E = ( [60000,50000,40000,30000,20000,10000,0], [60000,50000,40000,30000,20000,10000,0], [0,0,0,0,0,0,0], [60000,50000,40000,30000,20000,10000,0])

#Now we call a helper function to make sure the activation energy is decreasing with coverage. It returns true if the final activation energy will decrease with coverage.
checkResult = canteraKineticsParametersParser.descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array=reactions_parameters_array, piecewise_coverage_intervals_all_reactions=piecewise_coverage_intervals, E_offsets_array_all_reactions=modifiers_E, verbose = False)    

if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
    pass #For now we will let all cases through.
    
cti_top_info_filename =model_location+ model_name + "_cti_top_info.cti"
with open(cti_top_info_filename, "r") as cti_top_info_file:
    cti_top_info_string = cti_top_info_file.read()     


global observed_x_values #This can be set initially, then changed later if desired.
import processing_functions_tpd_odeint
observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
observed_x_values, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint.import_experimental_settings(observed_data_Filename)

global firstSimulation
firstSimulation = True


global modifiers_A_original 
modifiers_A_original = modifiers_A
global modifiers_E_original 
modifiers_E_original = modifiers_E

def cantera_simulation_wrapper_example5(parametersArray): #This takes in *only* the adjustable parameters. The modifiers for A and the modifiers for E are part of that.
    #Now we modify the 'initial' parameters array accordingly.
    #In this example, we're modifying reactions 2 and 4 which are indices 1 and 3.
    modified_reactions_parameters_array = copy.deepcopy(reactions_parameters_array)
    
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
    
    numIntervals = len(piecewise_coverage_intervals)
    numReactions = 4 #This is hard coded right now.
    num_modifiers_A = numIntervals
    num_modifiers_E = numIntervals
#    modifiers_A_starting_index = 6
#    modifiers_A_ending_index = modifiers_A_starting_index+num_modifiers_A
    new_modifiers_A = list(parametersArray[6:6+num_modifiers_A]) #The first 5 parameters are the original kinetic parameters.
    new_modifiers_E = list(parametersArray[6+num_modifiers_A:6+num_modifiers_A+num_modifiers_E])
    
    global modifiers_A_original 
    global modifiers_E_original 
    modifiers_A_list = list(modifiers_A_original) #grabbing original values / last values
    modifiers_E_list = list(modifiers_E_original) #grabbing original values / last values
#    print("line 74", parametersArray)
#    print("line 74", parametersArray[modifiers_A_starting_index:modifiers_A_ending_index])
    modifiers_A_list[1]=new_modifiers_A #In this example, we're only modifying reaction at index 1.
    modifiers_E_list[1]=new_modifiers_E #In this example, we're only modifying reaction at index 1.
    
    
    
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "A"
    modifiers_A = tuple(modifiers_A_list)
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "E"    
    modifiers_E = tuple(modifiers_E_list)
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, modified_reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_E)
    
    
    
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

totalAbsoluteRate = cantera_simulation_wrapper_example5([61.5, 41.5, 13.0, 13.0, 0.1, 0.1, #Regular kinetic params
                                                         -1,-1,-1,-1,-1,-1,0 , #A modifiers.
                                                         60000,50000,40000,30000,20000,10000,0]) #E modifiers.

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
    ceO2_input_simulation_settings.piecewise_coverage_dependence = True
    reactions_parameters_array = np.genfromtxt(model_location+ model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    #This is before any offset modification:
    ceO2_input_simulation_settings.original_reactions_parameters_array = reactions_parameters_array
    piecewise_coverage_intervals = np.array([0,0.1,0.2,0.3,0.4,0.5,1.0])
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "A"
    modifiers_A = ([0,0,0,0,0,0,0], [-1,1,-1,-1,-1,-1,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0])
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
    kineticParameterName = "E" 
    modifiers_E = ( [60000,50000,40000,30000,20000,10000,0], [60000,50000,40000,30000,20000,10000,0], [0,0,0,0,0,0,0], [60000,50000,40000,30000,20000,10000,0])
    
    #Now we call a helper function to make sure the activation energy is decreasing with coverage. It returns true if the final activation energy will decrease with coverage.
    checkResult = canteraKineticsParametersParser.descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array=reactions_parameters_array, piecewise_coverage_intervals_all_reactions=piecewise_coverage_intervals, E_offsets_array_all_reactions=modifiers_E, verbose = False)    
    
#    print("Initial Reactions Parameters Array")
#    print(reactions_parameters_array)
#    
    if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
        canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_E)
        
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
        