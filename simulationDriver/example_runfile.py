import cantera as ct
import cantera.ck2cti as ck2cti
import canteraSimulate
import numpy as np
import canteraKineticsParametersParser 


def main():
    #Here is an example usage. The things that must be defined in advance are...
    # a) a model name.
    # b) The simulation settings are set in a python file which then becomes imported and an argument. This can then be changed by a script.
    # c) The reactions_parameters_array must be fed as an array or as a file which contains such an array with a header. This includes the reaction equations and their Arrhenius initial values.
    # d) The cti_top_info_string. Right now, this file is made 'manually'.
    
    print("####EXAMPLE 1#####")
    #Here is an example that leaves the settings as they are in the settings file. 
    #For the cti_stop_info_string the setting "None" is used, but that just means the code will automatically look for it in a file. Though it could be provided instead.
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    ceO2_input_simulation_settings.print_frequency = 100 #This makes the simulation not print things out during checkpoints.
    canteraSimulate.create_cti_and_SimulatePFRorTPRwithCantera("ceO2", "ceO2_input_reactions_parameters.csv", ceO2_input_simulation_settings, cti_top_info_string = None, write_cti_to_file = True)
    #Note that one can change from static to PFR by changing a variable inside ceO2_simulation_settings. 
    #To get the output as python objects and not only as file exports, see example 2.

    print("####EXAMPLE 2#####")
    #Here is an example useage with a loop to scan different parameter settings.
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    model_name = "ceO2"
    reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    import copy
    modified_reactions_parameters_array = copy.deepcopy(reactions_parameters_array)
    #This time we are going to get the cti_top_info string so that we don't waste time reading it repeatedly.        
    cti_top_info_filename = model_name + "_cti_top_info.cti"
    with open(cti_top_info_filename, "r") as cti_top_info_file:
        cti_top_info_string = cti_top_info_file.read()     
    #Now we execute a loop with some changes...
    for multiplyingFactor in range(1,11):
        ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
        ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a paramter optimization.
        modified_reactions_parameters_array[1][5] = str(float(reactions_parameters_array[1][5])*(1+multiplyingFactor/10)) #This modifies the desorption rate constant. #note that the kinetic parameter starts as a string.
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.create_cti_and_SimulatePFRorTPRwithCantera("ceO2", modified_reactions_parameters_array, ceO2_input_simulation_settings, cti_top_info_string = cti_top_info_string)
        np.savetxt("ceO2_output_rates_all_"+str(multiplyingFactor)+".csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)
        print("round", multiplyingFactor, "Simulation Finished")
        
    print("####EXAMPLE 3#####")
    #This is like example 2, only now we use the modifyAllAdjustableReactionParametersArrayValues function rather than sticking values in by simple index.
    model_name = "ceO2"
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    
    #Let's get the modfiable parameter indices and the original values before any modifications...
    modifiableParameterIndices = canteraKineticsParametersParser.findModifiableParameterIndices(reactions_parameters_array)
    originalModifiableReactionParametersValues = canteraKineticsParametersParser.getAllModifiableReactionParametersValues(reactions_parameters_array, modifiableParameterIndices)
    print(modifiableParameterIndices)
    print(originalModifiableReactionParametersValues)
    #The original parameter values are:
    # ['1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '67400.0', '0.0', '0.0', '-6000.000000000001']
    
    #This time we are again going to get the cti_top_info string so that we don't waste time reading it repeatedly.        
    cti_top_info_filename = model_name + "_cti_top_info.cti"
    with open(cti_top_info_filename, "r") as cti_top_info_file:
        cti_top_info_string = cti_top_info_file.read()     
    #Now we execute a loop with some changes... but this time we're going to inject adjustedModifiableParameters
    adjustedModifiableParametersSamplings = [
    ['1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '67400.0', '0.0', '0.0', '-6000.000000000001', '1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '67400.0', '0.0', '0.0', '-6000.000000000001'], #First one is actually not modified.
    ['1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '87400.0', '0.0', '0.0', '-6000.000000000001', '1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '67400.0', '0.0', '0.0', '-6000.000000000001'],
    ['1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '107400.0', '0.0', '0.0', '-6000.000000000001', '1.0', '0.0', '2000.0', '10000000000000.0', '0.0', '67400.0', '0.0', '0.0', '-6000.000000000001']
    ] #IT IS IMPORTANT TO REALIZE THAT THESE ARE NOT THE FULL MOIDFIED PARAMETER ARRAYS, DON'T GET THE TWO CONFUSED!
    
    for caseIndex, adjustedModifiableParametersCase in enumerate(adjustedModifiableParametersSamplings):
        ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
        ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a paramter optimization.
        modified_reactions_parameters_array = canteraKineticsParametersParser.modifyAllAdjustableReactionParametersArrayValues(reactions_parameters_array, adjustedModifiableParametersCase)
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.create_cti_and_SimulatePFRorTPRwithCantera("ceO2", modified_reactions_parameters_array, ceO2_input_simulation_settings, cti_top_info_string = cti_top_info_string)
        np.savetxt("ceO2_output_rates_all_"+"FullCTIsamplingCase"+str(caseIndex)+".csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)              
        print("sampling case FullCTIsamplingCase", caseIndex, "Simulation Finished")

    print("####EXAMPLE 4#####")
    #This is like example 3, only now we use the modifyReactionsInOnePhase rather than full cti.
    model_name = "ceO2"
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    
    #This time we are again going to get the cti_top_info string so that we don't waste time reading it repeatedly.        
    cti_top_info_filename = model_name + "_cti_top_info.cti"
    with open(cti_top_info_filename, "r") as cti_top_info_file:
        cti_top_info_string = cti_top_info_file.read()     
    #Now we execute a loop with some changes... but this time we're going to inject entire already modified parameter arrays.
    adjustedModifiableParametersSamplings = [
            [
            ['0001','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde1-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0002','surface_reaction','Acetaldehyde1-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','67400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],
            ['0003','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde2-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0004','surface_reaction','Acetaldehyde2-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','67400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],
            ],#First one is actually not modified.
    
            [
            ['0001','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde1-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0002','surface_reaction','Acetaldehyde1-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','87400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],
            ['0003','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde2-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0004','surface_reaction','Acetaldehyde2-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','67400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],            
            ],
            
            
            [
            ['0001','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde1-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0002','surface_reaction','Acetaldehyde1-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','107400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],
            ['0003','surface_reaction','Acetaldehyde + CeCation(S) => Acetaldehyde2-Ce(S)','1.0','0.0','2000.0','nan','nan','nan','None','True'],
            ['0004','surface_reaction','Acetaldehyde2-Ce(S) => Acetaldehyde + CeCation(S)','10000000000000.0','0.0','67400.0','0.0','0.0','-6000.000000000001','Acetaldehyde(S)','False'],            
            ]
    ] 
    
    for caseIndex, adjustedModifiableParametersCase in enumerate(adjustedModifiableParametersSamplings):
        ceO2_input_simulation_settings.print_frequency = None #This makes the simulation not print things out during checkpoints.
        ceO2_input_simulation_settings.exportOutputs = True #This will be turned off during a paramter optimization.
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.modify_reactions_and_SimulatePFRorTPRwithCantera(model_name, adjustedModifiableParametersCase, ceO2_input_simulation_settings, canteraPhases=canteraPhases)        
        np.savetxt("ceO2_output_rates_all_"+"ModifyReactionssamplingCase"+str(caseIndex)+".csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)
        print("sampling case ModifyReactionssamplingCase", caseIndex, "Simulation Finished")
        
    print("####EXAMPLE 5b#####")
    #This is like example 4, only now we are going to use the applyPieceWiseCoverageDependence feature, which requires us to fill ceO2_input_simulation_settings.piecewise_coverage_dependences
    model_name = "ceO2"
    import ceO2_input_simulation_settings #The user may change settings in the python file with the same name.
    reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    #We are going to set the starting coverage of "Acetaldehyde1-Ce(S)"  to 1 and the rest to 0 so that it's effectively a 1 species PE-TPR experiment.
    ceO2_input_simulation_settings.surface_coverages = 'CeCation(S):0 Acetaldehyde1-Ce(S):1 Acetaldehyde2-Ce(S):0 OAnion(S):0'
        
    #This time we are again going to get the cti_top_info string so that we don't waste time reading it repeatedly.        
    cti_top_info_filename = model_name + "_cti_top_info.cti"
    with open(cti_top_info_filename, "r") as cti_top_info_file:
        cti_top_info_string = cti_top_info_file.read()     
        
    ceO2_input_simulation_settings.piecewise_coverage_dependence = True
    ceO2_input_simulation_settings.original_reactions_parameters_array = reactions_parameters_array
    piecewise_coverage_intervals = np.array([0,0.1,0.2,0.3,0.4,0.5,1.0])
    species_name = "Acetaldehyde1-Ce(S)" #In this example we make the parameters depend only on the coverage of species Acetaldehyde1-Ce(S).
    kineticParameterName = "A"
    modifiers_A = ( [0,0,0,0,0,0,0], [-1,-1,-1,-1,-1,-1,0], [0,0,0,0,0,0,0], [0,0,0,0,0,0,0])
    canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_A )
    kineticParameterName = "E" 
    modifiers_E = ( [60000,50000,40000,30000,20000,10000,0], [60000,50000,40000,30000,20000,10000,0], [0,0,0,0,0,0,0], [60000,50000,40000,30000,20000,10000,0])
    
    #Now we call a helper function to make sure the activation energy is decreasing with coverage. It returns true if the final activation energy will decrease with coverage.
    checkResult = canteraKineticsParametersParser.descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array=reactions_parameters_array, piecewise_coverage_intervals_all_reactions=piecewise_coverage_intervals, E_offsets_array_all_reactions=modifiers_E, verbose = False)
    if checkResult == True: #In this example, we are only proceeding if the coverage dependance is always decreasing with coverage.
        canteraKineticsParametersParser.populatePiecewiseCoverageDependence(ceO2_input_simulation_settings, reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_E)
        #Now simulate!
        concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
        canteraSimulate.modify_reactions_and_SimulatePFRorTPRwithCantera(model_name, adjustedModifiableParametersCase, ceO2_input_simulation_settings, canteraPhases=canteraPhases)
        np.savetxt("ceO2_output_rates_all_"+"Piecewise_coverage_Dependence"+".csv", rates_all_array, delimiter=",", comments='', header=rates_all_array_header)
        print("Piecewise_coverage_Dependence", "Simulation Finished")

if __name__ == '__main__':
    main()    