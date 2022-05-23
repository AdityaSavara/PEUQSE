import sys; sys.path.append('../../'); #This is to get access to the simulationDriver directory below.
import numpy as np
from PEUQSE.simulationDriver import canteraSimulate
from PEUQSE.simulationDriver import canteraKineticsParametersParser
import copy

#First made a 'test' in main, then made a simulation wrapper to do the same thing.

#We will define some things *outside* of the function to load the model.
model_location = ".\\" #For example5, the model will be in the same directory.
model_name = "Example5" #this refers to the cantera model name, not PEUQSE.
import Example5_input_simulation_settings #The user may change settings in the python file with the same name.
Example5_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
Example5_input_simulation_settings.exportOutputs = True #This will be turned off during a parameter optimization to reduce I/O.
Example5_input_simulation_settings.piecewise_coverage_dependence = True
Example5_input_simulation_settings.surface_coverages = 'CeCation(S):0 Acetaldehyde1-Ce(S):1.00 Acetaldehyde2-Ce(S):0.0 OAnion(S):0.0'
reactions_parameters_array = np.genfromtxt(model_location+ model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
##ow, looking at the paramters array file, we can see where the different parameters positions are that we can modify.
print(reactions_parameters_array)

#This is before any offset modification:
Example5_input_simulation_settings.original_reactions_parameters_array = reactions_parameters_array

global piecewise_coverage_intervals
piecewise_coverage_intervals = np.array([0,0.3,0.50,0.70,0.90,.95,1.0])
species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
kineticParameterName = "A"
modifiers_A = ([0,0,0,0,0,0,0], [-1,1,-1,-1,-1,-1,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0])
canteraKineticsParametersParser.populatePiecewiseCoverageDependence(Example5_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
kineticParameterName = "E" 
modifiers_E = ( [60000,50000,40000,30000,20000,10000,0], [60000,50000,40000,30000,20000,10000,0], [0,0,0,0,0,0,0], [60000,50000,40000,30000,20000,10000,0])

#Now we call a helper function to make sure the activation energy is decreasing with coverage. It returns true if the final activation energy will decrease with coverage.
checkResult = canteraKineticsParametersParser.descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array=reactions_parameters_array, piecewise_coverage_intervals_all_reactions=piecewise_coverage_intervals, E_offsets_array_all_reactions=modifiers_E, verbose = False)    

if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
    pass #For now we will let all cases through.
    
yaml_top_info_filename =model_location+ model_name + "_yaml_top_info.yaml"
with open(yaml_top_info_filename, "r") as yaml_top_info_file:
    yaml_top_info_string = yaml_top_info_file.read()     


global observed_x_values #This can be set initially, then changed later if desired.
from processing_functions_tpd_odeint import import_experimental_settings
observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors.csv'
observed_x_values, responses_observed, observedResponses_uncertainties = import_experimental_settings(observed_data_Filename)

global firstSimulation
firstSimulation = True


global modifiers_A_original 
modifiers_A_original = modifiers_A
global modifiers_E_original 
modifiers_E_original = modifiers_E


def test_run():
#    
    if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
        canteraKineticsParametersParser.populatePiecewiseCoverageDependence(Example5_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_E)
        
        yaml_top_info_filename =model_location+ model_name + "_yaml_top_info.yaml"
        with open(yaml_top_info_filename, "r") as yaml_top_info_file:
            yaml_top_info_string = yaml_top_info_file.read()     
        
        modified_reactions_parameters_array = reactions_parameters_array
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.create_yaml_and_SimulatePFRorTPRwithCantera("Example5", modified_reactions_parameters_array, Example5_input_simulation_settings, yaml_top_info_string = yaml_top_info_string)
        
        #Playing with the variables we see we can get what we want from:
    #    print(cantera_phase_rates['surf'])
        times = cantera_phase_rates['surf'][:,1] #This is in seconds.
        temperatures = Example5_input_simulation_settings.heating_rate*times+Example5_input_simulation_settings.T_surf
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
            
        #np.savetxt("Example5_output_rates_all_"+"FullYAMLsamplingCase.csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)              
        


def cantera_simulation_wrapper_example5(parametersArray): #This takes in *only* the adjustable parameters. The modifiers for A and the modifiers for E can be part of that.
    Example5_input_simulation_settings.exportOutputs = False 
    #now we have to modify the cantera parameters that need to be modified.
    #In this example, we're modifying reaction 2, which is index 1 in the reactions list.
    #We are **only** going to modify the modifier_E values.
    modifiers_E = np.array(modifiers_E_original) #Note that by creating modifers_E here, we are making a *local* variable, not the global one.
    #What we are actually going to do is make the modifiers "stacked" to ensure a decrease is likely.    
    
    stackedValues = np.zeros(len(parametersArray))#just initializing.
    stackedValues[0] =-1*parametersArray[0] #We are subtracting from the activation energies.
    for valueIndex in range(len( parametersArray)):
        if valueIndex > 0:
            stackedValues[valueIndex] = stackedValues[valueIndex-1] - parametersArray[valueIndex] #become "more negative" each time.
    modifiers_E[1] = stackedValues  #replace E[1] with the new piecweise function based on the sampled parameters.

    #Now we call a helper function to make sure the activation energy is decreasing with coverage. It returns true if the final activation energy will decrease with coverage.
    checkResult = canteraKineticsParametersParser.descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array=reactions_parameters_array, piecewise_coverage_intervals_all_reactions=piecewise_coverage_intervals, E_offsets_array_all_reactions=modifiers_E, verbose = False)    
    if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
        pass #For now we will let all cases through.
    if checkResult == False:
        return None #we treat the simulation like a failed simulation if the Ea was not decreasing.
    
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "A"
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(Example5_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "E"    
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(Example5_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_E)
    
    #Now we do the simulation. The first time, we need to do a full simulation flow.  But second time and later we can just modify the cantera phases object and run again.
    global firstSimulation
    global canteraPhases #need to make it global so the object sticks around for next function call.
    if firstSimulation == True: #NOTE: Below we feed in reactions_parameters_array, but it will become overwritten later during simulation.
        try:
            concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
            canteraSimulate.create_yaml_and_SimulatePFRorTPRwithCantera(model_name, reactions_parameters_array, Example5_input_simulation_settings, yaml_top_info_string = yaml_top_info_string)
            firstSimulation = False #Toggle the flag.
        except:  #if a simulation fails we return a None type object.
            firstSimulation = False #Toggle the flag.
            return None 
    elif firstSimulation == False: #This must be an elif, otherwise it will always be executed. In this function, the cantera phases object will be created the first time the function is called. Then will exist for later.
        try: #NOTE: Below we feed in reactions_parameters_array, but it will become overwritten later during simulation.
            concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
            canteraSimulate.modify_reactions_and_SimulatePFRorTPRwithCantera(model_name, reactions_parameters_array, Example5_input_simulation_settings, canteraPhases=canteraPhases)        
        except: #if a simulation fails we return a None type object.
            return None
    
    #Now we parse the output.
    global times
    times = cantera_phase_rates['surf'][:,1] #This is in seconds.
    temperatures = Example5_input_simulation_settings.heating_rate*times+Example5_input_simulation_settings.T_surf
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


def integrated_cantera_simulation_wrapper_example5(discreteParameterVector): 
    simulationOutput = cantera_simulation_wrapper_example5(discreteParameterVector)
    if type(simulationOutput) == type(None): #return a None type object if that's what has been received.
        print("A cantera simulation has failed. If this happens for a large percentage of the samplings, the posterior will not be sufficiently sampled.")
        return None
    #Else is the normal case:
    rate = simulationOutput #This has already been interpolated to observed x values.

    global times
    global observed_x_values
    from PEUQSE import littleEulerGivenArray
    times, integrated_desorption, rate = littleEulerGivenArray(0, observed_x_values, rate)
    return integrated_desorption


if __name__ == "__main__":
    test_run()
    #Below are testing the wrapper. 
    #Below, #The x axis is not correct. We are simply looking for  "differences" with earlier desorption for wrapperOutput2, which we do see.
    #The first test case will be saved as wrapperOutput1...
    wrapperOutput1 = cantera_simulation_wrapper_example5([0,0,0,0,0,0,0])
    from matplotlib import pyplot as plt
    if type(wrapperOutput1) == type(None):
        print("WrappeOutput1 failed.")
    else:
        plt.plot(wrapperOutput1)
        plt.title('wrapperOutput1')
        plt.show()   
    #The second test case will be saved as wrapperOutput2...
    wrapperOutput2 = cantera_simulation_wrapper_example5([0,0,0,0,10000,10000, 10000])
    if type(wrapperOutput2) == type(None):
        print("WrappeOutput2 failed.")
    else:
        plt.plot(wrapperOutput2)
        plt.title('wrapperOutput2')
        plt.show()   
            