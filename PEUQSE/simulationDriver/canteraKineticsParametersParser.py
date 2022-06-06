# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:44:31 2019

#originally made for Cantera 2.4.  Not upgraded to work with later versions of Cantera, as of Oct 2021
"""

import xmltodict
import numpy as np
import cantera as ct 
import copy
#This module is written for reactions_parameters_arrays in the following format:
#headerString = "reactionID,canteraReactionType,reactionEquation,A,b,E,a,m,e,concentrationDependenceSpecies, is_sticking_coefficient"
#The csv files are comma separated, and have 1 header string as noted above.


def stackListsAndArrays(listOfItemsToStack):
    stackableItemsList = copy.deepcopy(listOfItemsToStack)
    for itemIndex in range(len(listOfItemsToStack)):
        #if it is a list, we turn it into a 2D array of the appropriate geometry.
        if str(type(listOfItemsToStack[itemIndex]))=="<class 'list'>":
            stackableItemsList[itemIndex] = np.atleast_2d(stackableItemsList[itemIndex]).transpose()
        else:
            pass
    return np.hstack(tuple(stackableItemsList)) #convert it to array and stack it. hstack requires a tuple as an argument.


def descendingLinearEWithPiecewiseOffsetCheckOneReactionAllReactions(reactions_parameters_array, piecewise_coverage_intervals_all_reactions, E_offsets_array_all_reactions, verbose = False):
    #Note that E_0 and g_slope are contained inside reactions_parameters_array
    #piecewise_coverage_intervals_all_reactions can either be 1D (defined the same for all reactions) or can be 2D (different intervals for each reaction)
    #piecewise_coverage_intervals_all_reactions and E_offsets_array_all_reactions can be staggered if they are 2D, but must have one row for each reaction.
    #NOTE: this calculates a single E_offsets_array_all_reactions, but actually there is going to be one for each species. The complication is that during simulations, there are many different coverage combinations possible
    #It would probably be better to sample a grid of possible species coverages that might be encountered, but for now this function will work for only one species.
    passedArray = np.zeros(len(reactions_parameters_array), dtype="bool") #This is an array that starts as all False, and we fill with true for each reaction that passes.
    for reactionIndex in range(len(reactions_parameters_array)):
        individualreactions_parameters_array = reactions_parameters_array[reactionIndex]  #TODO: it might be better to pass the individual reaction parameters array to the checking function, which would mean changing the checking function.
        if len(np.shape(reactions_parameters_array))==1: #Check if there is a single set of intervals to be used for all reactions.
            piecewise_coverage_intervals = piecewise_coverage_intervals_all_reactions
        else: #The else means different intervals for each reaction.
            piecewise_coverage_intervals = piecewise_coverage_intervals_all_reactions[reactionIndex]
        E_offsets_array = E_offsets_array_all_reactions[reactionIndex]
        #Populate the Arrhenius type parameters related to E and its slope...
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
#        A = str(individualreactions_parameters_array[3])
#        b = str(individualreactions_parameters_array[4])
        E = str(individualreactions_parameters_array[5])
#        a = str(individualreactions_parameters_array[6])
#        m = str(individualreactions_parameters_array[7])
        e = str(individualreactions_parameters_array[8])  #Note: in 2022, this is now "E" in cantera, but used to be "e".
        concentrationDependenceSpecies = str(individualreactions_parameters_array[9])

        #Get the numbers that we need out.
        E_0 = float(E)
        if np.isnan(float(e)): #Need to make e into 0.0 if it's an nan. Otherwise the descending check will fail.
            e = 0.0
        g_slope = -1.0*float(e) #Note that in 2022 in cantera "E" (formerly "e") is the negative of g_slope! https://cantera.org/science/kinetics.html  https://cantera.org/science/kinetics.html#sec-surface
        #print(reactionID, reactionEquation, e, g_slope)
        if '(S)' in reactionEquation: #This means the reaction is a surface reaction  and we will check it. Otherwise, we will not.
            passedArray[reactionIndex] = descendingLinearEWithPiecewiseOffsetCheckOneReaction(E_0, g_slope,piecewise_coverage_intervals, E_offsets_array)
            if verbose == True:
                if passedArray[reactionIndex] == False:
                    print("The reaction with ID of" + reactionID + "and equation of" + reactionEquation + "did not pass descendingLinearEWithPiecewiseOffsetCheckOneReaction" )
    if sum(passedArray )< len(passedArray): #For a boolean array, if the sum is less than the length, that means at least one thing remained false.
        return False
    elif sum(passedArray) == len(passedArray): #Else the sum must be equal and we return true.
        return True
    
    #This probably only works for the case of a single species, given how the g_slope is defined currently.
def descendingLinearEWithPiecewiseOffsetCheckOneReaction(E_0, g_slope,piecewise_coverage_intervals, E_offsets_array):
    E_array =E_0 + g_slope*piecewise_coverage_intervals+E_offsets_array
    delta_E_array = np.diff(E_array)
    if np.any(delta_E_array > 0):
        return False
    else:
        return True

def makeCanteraReactionObjectsListFromFile(FileName):
    #Creates a list of Reaction objects from all of the reactions defined in a YAML, CTI, or XML file.
    #The CTI and XML input formats are deprecated since Cantera 2.5 and will be removed in Cantera 3.0.
    #There is utility to convert old CTI files to the new YAML. https://cantera.org/tutorials/legacy2yaml.html
    #inside cti2yaml is a function called convert: convert(filename=None, output_name=None, text=None) so can be used to convert.
    canteraReactionObjectsList = ct.Reaction.listFromFile(FileName)
    return canteraReactionObjectsList


def extractReactionParametersFromFile(InputFileName, OutputFilename = ""): #The input fileName must point to a cti file or an xml file. An optional argument of OutputFilename will export the reactionParamtersArray
    canteraReactionObjectsList = makeCanteraReactionObjectsListFromFile(InputFileName)
    reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList = getReactionParametersFromCanteraReactionObjectsList(canteraReactionObjectsList)
    outputAsNumpyArray = stackListsAndArrays([reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList])
    if OutputFilename != "":
        headerString = "reactionID,canteraReactionType,reactionEquation,A,b,E,a,m,e,concentrationDependenceSpecies, is_sticking_coefficient"
        np.savetxt(OutputFilename,outputAsNumpyArray, fmt="%s", delimiter=",", comments='', header=headerString)
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList

def getReactionParametersFromCanteraReactionObjectsList(reactionObjectList):
    reactionIDsList = []
    reactionTypesList = []
    reactionEquationsList = []
    ArrheniusParametersList =[] #List of numpy arrays with A,b,E for each reaction.
    concentrationDependencesList =[]
    concentrationDependencesSpeciesList = []
    is_sticking_coefficientList =[]
    for reactionObject in reactionObjectList:
        reactionIDsList.append(reactionObject.ID)
        try:
            if int(reactionObject.reaction_type) == int(20):
                reactionTypesList.append("surface_reaction")
            if int(reactionObject.reaction_type) == int(1):
                reactionTypesList.append("reaction")
            if int(reactionObject.reaction_type) == int(2):
                reactionTypesList.append("three_body_reaction")
            if int(reactionObject.reaction_type) == int(4):
                reactionTypesList.append("falloff_reaction")
        except:
            reactionTypesList.append("None")
        #Now need to modify equation based on whether it's irreversible or not.
        reactionEquationsList.append( reactionObject.equation)
        A = float(reactionObject.rate.pre_exponential_factor)
        b = float(reactionObject.rate.temperature_exponent)
        E = float(reactionObject.rate.activation_energy)/1000 #Annoyingly, Cantera puts Energies out in J/kmol.
        singleReactionArrheniusParametersArray = np.array([A,b,E])
        ArrheniusParametersList.append(singleReactionArrheniusParametersArray)
        #Now check for coverage parameters.
        try:
            if len(reactionObject.concentration_deps) > 0:
                concentration_deps = True
            else:
                concentration_deps = False
        except:
                concentration_deps = False
        if concentration_deps == True:
            singleReactionconcentrationDependencesDict = reactionObject.concentration_deps
            temporaryconcentrationDependenceSpeciesList = list(reactionObject.concentration_deps.keys())
            if len(temporaryconcentrationDependenceSpeciesList) > 0:
                concentrationDependenceExists = True
                concentrationDependenceSpecies = temporaryconcentrationDependenceSpeciesList[0]
                a = reactionObject.concentration_deps[concentrationDependenceSpecies][0]
                m = reactionObject.concentration_deps[concentrationDependenceSpecies][1]
                e = float(reactionObject.concentration_deps[concentrationDependenceSpecies][2])/1000 #Annoyingly, Cantera puts Energies out in J/kmol.   
        else:#There is no coverage dependence, so will put 'NaN')
             a = float('nan')
             m = float('nan')
             e = float('nan')
             concentrationDependenceSpecies = 'None'
        concentrationDependencesList.append(np.array([a,m,e]))
        concentrationDependencesSpeciesList.append(concentrationDependenceSpecies)
        is_sticking_coefficientList.append(reactionObject.is_sticking_coefficient)
    reactionIDsList = reactionIDsList
    reactionTypesList = reactionTypesList
    ArrheniusParametersArray = np.array(ArrheniusParametersList)
    concentrationDependencesArray = np.array(concentrationDependencesList)
    concentrationDependencesSpeciesList = concentrationDependencesSpeciesList
    is_sticking_coefficientList = is_sticking_coefficientList
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList

def make_reaction_yaml_string(individualreactions_parameters_array, for_full_yaml = False, input_Ea_units="J/kmol"):
    #first input should be a an array of strings: 
    #reactionID	canteraReactionType	reactionEquation	A	b	E	a	m	e	concentrationDependenceSpecies	 is_sticking_coefficient
    # 1	20	H2 + 2 PT(S) => 2 H(S)	0.046	0	0	0	-1	0	PT(S)	TRUE
    # Actually looks like the below, in practice:
    #['4' 'surface_reaction' 'Acetaldehyde2-Ce(S) => Acetaldehyde + CeCation(S)' '10000000000000.0' '0.0' '67400.0' '0.0' '0.0' '-6000.000000000001' 'Acetaldehyde2-Ce(S)'  'FALSE']
    #canteraReactionType is typically "reaction". THe falloff and three_body reactions are not yet supported at this time.  For yaml as of Cantera >2.5, there is no "surface_reaction", instead, surface_reactions are now simply "reactions".
    ##
    #the "for_full_yaml" option can be True or False. It is important because when adding a reaction to a model, the yaml looks different compared to when it's added to a file.
    #For example, below is the string one would use to add a yaml reaction **to an existing cantera object**
    #'''
    #  equation: O2 + 2 PT(S) => 2 O(S)
    #  rate-constant: {A: 1.8900000000000004e+19, b: -0.5, Ea: 0.0}
    #'''
    #and below is how one would add the reaction to a yaml file (note the dash).
    #'''
    #- equation: O2 + 2 PT(S) => 2 O(S)
    #  rate-constant: {A: 1.8900000000000004e+19, b: -0.5, Ea: 0.0}
    #'''
    #In reality, we will use less linebreaks than the above two examples, but the main point is that the first character for each reaction is a dash when adding to file.
    #
    #Below is how coverage-dependencies are included (note that there is an extra indent).
    #'''
    #- equation: O2 + 2 PT(S) => 2 O(S)
    #  rate-constant: {A: 1.8900000000000004e+19, b: -0.5, Ea: 0.0}
    #  coverage-dependencies:
    #    PT(S): {a: 0.0, m: -1.0, E: 0.0}
    #'''    
    #NOTE: The default converage-dependencies will mean "no change" in the rate constant and thus would have values of a=0, m=0, and e=0 (e=0 is now E=0 in the yaml way of cantera).
    reaction_yaml_string = None #This is to help diagnose errors.
    reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
    individualReactionTypeReceived = str(individualreactions_parameters_array[1])
    reactionEquation = str(individualreactions_parameters_array[2])    
    if input_Ea_units.lower() == 'j/mol': #convert j/mol to j/kmol.
        individualreactions_parameters_array[5] = float(individualreactions_parameters_array[5])*1000.0 #This is the position of Ea in the array.
    A = str(individualreactions_parameters_array[3])
    b = str(individualreactions_parameters_array[4])
    E = str(individualreactions_parameters_array[5]) #Note: this is now "Ea" in cantera 2.6
    a = str(individualreactions_parameters_array[6])
    m = str(individualreactions_parameters_array[7])
    e = str(individualreactions_parameters_array[8]) #NOTE: This is now "E" in cantera 2.6.
    concentrationDependenceSpecies = str(individualreactions_parameters_array[9])  #TODO: it is now possible to have multiple concentration dependent species. My individualreactions_parameters_array syntax will need to be updated to accommodate this. It seems this will need to become something parseable, just as the reactionEquation is parseable. 
    is_sticking = str(individualreactions_parameters_array[10])
    if is_sticking.capitalize() == "False": #considered using distutils.util.strtobool(is_sticking) but decided that would slow the program down.
        ArrheniusString = "Arrhenius"
    if is_sticking.capitalize() == "True":
        ArrheniusString = "stick"
    #In the new yaml way of doing things, there is a field for the rate_constant type to distinguish sticking coefficients from other rate constants.
    if is_sticking.capitalize() == "True":
        rate_constant_type = "sticking-coefficient"
    if is_sticking.capitalize() == "False":
        rate_constant_type = "rate-constant"        
        
    if individualReactionTypeReceived == 'surface_reaction': #For Cantera yaml, there is no longer a special designation for surface_reactions.
        individualReactionTypeReceived = 'reaction' 
    #Now make the reaction_yaml strings...
    if individualReactionTypeReceived.lower() == "reaction":
        #The spacing will be "strange" here because the ''' does not ignore linebreaks, it keeps them, and they are not normal linebreaks.
        rxnStringTemplate = \
'''  equation: Eqn_string
  rate_constant_type: {A: A_value, b: b_value, Ea: Ea_value}'''
        reaction_yaml_string = rxnStringTemplate #just initializing.
        reaction_yaml_string = reaction_yaml_string.replace('Eqn_string', str(reactionEquation))
        reaction_yaml_string = reaction_yaml_string.replace('A_value', str(A))
        reaction_yaml_string = reaction_yaml_string.replace('b_value', str(b))
        reaction_yaml_string = reaction_yaml_string.replace('Ea_value', str(E))
        reaction_yaml_string = reaction_yaml_string.replace('rate_constant_type', str(rate_constant_type))
        if str(concentrationDependenceSpecies).lower() != 'none': #none means no species was provided. If a species was provided, then we will append the coverage dependance.
        #The spacing will be "strange" here because the ''' does not ignore linebreaks, it keeps them, and they are not normal linebreaks.
            coverageDependenceStringTemplate = \
'''
  coverage-dependencies:
    speciesName: {a: a_value, m: m_value, E: e_value}'''
            coverageDependenceString = coverageDependenceStringTemplate #just initializing
            coverageDependenceString = coverageDependenceString.replace('speciesName', str(concentrationDependenceSpecies))
            coverageDependenceString = coverageDependenceString.replace('a_value', str(a))
            coverageDependenceString = coverageDependenceString.replace('m_value', str(m))
            coverageDependenceString = coverageDependenceString.replace('e_value', str(e))
        else:
            coverageDependenceString = '' #make a blank string if there is no coverage dependence.
        #now add the two strings:
        reaction_yaml_string = reaction_yaml_string + coverageDependenceString

    #now modfiy the reaction in the actual model:
    if for_full_yaml == True:
        reaction_yaml_string = '-' + reaction_yaml_string[1:] #make the first character a dash if the string is going to be used for a yaml file.
    return reaction_yaml_string
    
def make_reaction_cti_string(individualreactions_parameters_array):
    #input should be a an array of strings: 
    #reactionID	canteraReactionType	reactionEquation	A	b	E	a	m	e	concentrationDependenceSpecies	 is_sticking_coefficient
    # 1	20	H2 + 2 PT(S) => 2 H(S)	0.046	0	0	0	-1	0	PT(S)	TRUE
    #canteraReactionType is typically "surface_reaction" or "reaction". THe falloff and three_body reactions are not yet supported at this time.
    reaction_cti_string = None #This is to help diagnose errors.
    reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
    individualReactionTypeReceived = str(individualreactions_parameters_array[1])
    reactionEquation = str(individualreactions_parameters_array[2])
    A = str(individualreactions_parameters_array[3])
    b = str(individualreactions_parameters_array[4])
    E = str(individualreactions_parameters_array[5])
    a = str(individualreactions_parameters_array[6])
    m = str(individualreactions_parameters_array[7])
    e = str(individualreactions_parameters_array[8])
    concentrationDependenceSpecies = str(individualreactions_parameters_array[9])
    is_sticking = str(individualreactions_parameters_array[10])
    if is_sticking.capitalize() == "False": #considered using distutils.util.strtobool(is_sticking) but decided that would slow the program down.
        ArrheniusString = "Arrhenius"
    if is_sticking.capitalize() == "True":
        ArrheniusString = "stick"

    #TODO: I need to make the coverage dependant defaults no coverage dependence if the person has entered "None".
    #NOTE: The default converage-dependencies will mean "no change" in the rate constant and thus would have values of a=0, m=0, and e=0 (e=0 is now E=0 in the yaml way of cantera).

    #Now make the reaction_cti strings...
    if individualReactionTypeReceived.lower() == "reaction":
        reaction_cti_string = '''reaction('{0}',
        {1}({2}, {3}, {4}])'''.format(reactionEquation, ArrheniusString, A,b,E)
    if individualReactionTypeReceived.lower() == "surface_reaction":
        if concentrationDependenceSpecies.capitalize() == "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
            reaction_cti_string= '''surface_reaction("{0}",
            {1}({2}, {3}, {4}))'''.format(reactionEquation,ArrheniusString,A,b,E) 
                #Like the below example.
                #        R3 = ct.InterfaceReaction.fromCti('''surface_reaction('O + H2 <=> H + OH',
                #        [3.87e1, 2.7, 2.619184e7])''')                      
        if concentrationDependenceSpecies.capitalize() != "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
            reaction_cti_string = '''surface_reaction("{0}",
            {1}({2}, {3}, {4},
            coverage = ['{5}', {6}, {7}, {8}]))'''.format(reactionEquation,ArrheniusString,A,b,E,concentrationDependenceSpecies,a,m,e)
                #Like the below example.
                #        R5 = ct.InterfaceReaction.fromCti('''surface_reaction( "CH4 + PT(S) + O(S) => CH3(S) + OH(S)",
                #                  Arrhenius(5.0e18, 0.7, 2000.0,
                #                            coverage = ['O(S)', 0.0, 0.0, 8000]))''')
    
    #now modfiy the reaction in the actual model:
    return reaction_cti_string


def create_full_cti(model_name, reactions_parameters_array = [], cti_top_info_string = None, write_cti_to_file = False):
    #It's convenient to use only the modelname. This then REQUIRES the reaction parameters and top info to
    #already be inside files with names of model_name+"_input_reactions_parameters.csv" and model_name+"_cti_top_info.cti"
    if cti_top_info_string == None:
       cti_top_info_filename = model_name + "_cti_top_info.cti"
       with open(cti_top_info_filename, "r") as cti_top_info_file:
           cti_top_info_string = cti_top_info_file.read() 
    #Normally, somebody will have to fill out the reaction parameters file ahead of time.
    #This time, it was made through the following line:
    #reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList = exportReactionParametersFromFile("ceO2_cti_full_existing.cti", "ceO2_input_reactions_parameters.csv")
    if len(reactions_parameters_array) == len([]):
       reaction_parameters_filename = model_name+"_input_reactions_parameters.csv"
       reactions_parameters_array = np.genfromtxt(reaction_parameters_filename, delimiter=",", skip_header=1, dtype="str")
       
    reactionsDataHeader = '\n\n\
#------------------------------------------------------------------------------- \n\
#  Reaction data \n\
#------------------------------------------------------------------------------- \n\
                           \n'
    #Now make the reactions string.
    cti_reactions_string = '\n' #Good to start with a new line.
    for element in reactions_parameters_array:
        cti_reactions_string += make_reaction_cti_string(element) + "\n" #need line breaks between each reaction_cti_string
    
    #Now generate the full cti.  
    cti_string = cti_top_info_string+reactionsDataHeader+cti_reactions_string
    if write_cti_to_file == True:
        cti_filename = model_name + "_cti_full.cti" 
        with open(cti_filename, 'w') as cti_file:
            cti_file.write(cti_string)
    return cti_string
    
    #This file takes the model_name, the reactions_parameters_array, and the yaml_top_info_string to create the **full** yaml string and/or yaml file.
def create_full_yaml(model_name, reactions_parameters_array = [], yaml_top_info_string = '', write_yaml_to_file = False):
    #It's convenient to use only the model_name. This then REQUIRES the reaction parameters and top info to
    #already be inside files with names of model_name+"_input_reactions_parameters.csv" and model_name+"_yaml_top_info.yaml"
    if yaml_top_info_string == '':
       yaml_top_info_filename = model_name + "_yaml_top_info.yaml"
       with open(yaml_top_info_filename, "r") as yaml_top_info_file:
           yaml_top_info_string = yaml_top_info_file.read() 
    #Normally, somebody will have to fill out the reaction parameters file ahead of time, either manually or with code.
    #This time, it was made through the following line:
    #reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList, is_sticking_coefficientList = exportReactionParametersFromFile("ceO2_cti_full_existing.cti", "ceO2_input_reactions_parameters.csv")
    #TODO: exportReactionParametersFromFile should be extended to work for YAML as well.
    if len(reactions_parameters_array) == len([]): #this means that the reactions_parameters_array needs to be loaded from file.
       reaction_parameters_filename = model_name+"_input_reactions_parameters.csv"
       reactions_parameters_array = np.genfromtxt(reaction_parameters_filename, delimiter=",", skip_header=1, dtype="str")

    reactionsDataHeader = 'reactions:\n'
    #Now make the reactions string.
    yaml_reactions_string = '' #Start with an empty string. Should not put extra line-breaks.
    for element in reactions_parameters_array:
        yaml_reactions_string += make_reaction_yaml_string(element, for_full_yaml=True) + "\n" #need line breaks between each reaction_yaml_string
    
    #Now generate the full yaml_string
    yaml_string = yaml_top_info_string+reactionsDataHeader+yaml_reactions_string
    
    if write_yaml_to_file == True:
        yaml_filename = model_name + "_yaml_full.yaml"
        with open(yaml_filename, 'w') as yaml_file:
            yaml_file.write(yaml_string)
           
           
           
    return yaml_string #this is the full yaml file string.


def findModifiableParameterIndices(reactions_parameters_array):
    reactions_parameters_array = np.atleast_2d(reactions_parameters_array) #If for some reason a person is using only 1 reaction, we still need to make it 2D.
    listOfModifiableParameterIndices = []
    listOfOriginalValues = []
    for reactionIndex, individualreactions_parameters_array in enumerate(reactions_parameters_array):
        #This is the format:
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
        A = str(individualreactions_parameters_array[3])
        b = str(individualreactions_parameters_array[4])
        E = str(individualreactions_parameters_array[5])
        a = str(individualreactions_parameters_array[6])
        m = str(individualreactions_parameters_array[7])
        e = str(individualreactions_parameters_array[8])
        concentrationDependenceSpecies = str(individualreactions_parameters_array[9])
        is_sticking = str(individualreactions_parameters_array[10])

        #For all reactions, we can modify the A,b,E parameters.
        listOfModifiableParameterIndices = listOfModifiableParameterIndices + [ [reactionIndex,3],[reactionIndex,4],[reactionIndex,5] ]
        #Now check for other parameters...
        if individualReactionTypeReceived == 'surface_reaction': #This can include adsorption, desorption, and surface reactions.
            if concentrationDependenceSpecies != 'None': #If this is not none, then that means the a, m, and e can also be modified.
                listOfModifiableParameterIndices = listOfModifiableParameterIndices + [ [reactionIndex,6],[reactionIndex,7],[reactionIndex,8] ]
                #Currently, we cannot modify the concentrationDependenceSpecies.
        #Other reaction types are not yet supported, but can be added in the future.                

    return listOfModifiableParameterIndices

def getAllModifiableReactionParametersValues(reactions_parameters_array, listOfModifiableParameterIndices=[]):
    if listOfModifiableParameterIndices == []: #If this has not been provided, we can get it.
        listOfModifiableParameterIndices = findModifiableParameterIndices(reactions_parameters_array)
    listOfOriginalValues = [] #Make the list that will be populated and returned.
    for modifiableParameterIndex in range(len(listOfModifiableParameterIndices)):
        reactionIndex = listOfModifiableParameterIndices[modifiableParameterIndex][0] #This is like the reaction ID, only it starts at 0.
        individualReactionParameterIndex = listOfModifiableParameterIndices[modifiableParameterIndex][1]#This says which column, basically whether it is A, b, E, etc.
        originalValue = reactions_parameters_array[reactionIndex][individualReactionParameterIndex] #Note that we use strings in the reaction parameters array.
        listOfOriginalValues.append(originalValue)
    return listOfOriginalValues

def modifyAllAdjustableReactionParametersArrayValues(reactions_parameters_array, newParametersList): #The first argument is a 2D numpy array. The second is a 1D list.
    reactions_parameters_array = np.atleast_2d(reactions_parameters_array) #If for some reason a person is using only 1 reaction, we still need to make it 2D.
    listOfModifiableParameterIndices = findModifiableParameterIndices(reactions_parameters_array)
    if len(listOfModifiableParameterIndices) != len(newParametersList):
        print("WARNING: modifyAllAdjustablereactions_parameters_arrayValues has been called with a newParametersList that is not of matching length.")
    #Now to make the modified parameters array...
    modified_reactions_parameters_array = copy.deepcopy(reactions_parameters_array)
    for newParameterIndex in range(len(newParametersList)):
        newValue = newParametersList[newParameterIndex]
        reactionIndex = listOfModifiableParameterIndices[newParameterIndex][0] #This is like the reaction ID, only it starts at 0.
        individualReactionParameterIndex = listOfModifiableParameterIndices[newParameterIndex][1]#This says which column, basically whether it is A, b, E, etc.
        modified_reactions_parameters_array[reactionIndex][individualReactionParameterIndex] = str(newValue) #Note that we use strings in the reaction parameters array.
    return modified_reactions_parameters_array

    
'''
#DEPRECATED BECAUSE NOW WE USE CANTERA CODE OR CANTERA OBJECTS
def makeDictfromXMLFile(FileName): #Filename is a string, can include path.
    with open("methane_pox_on_pt.xml") as xmlfile:
        xmldict = xmltodict.parse(xmlfile.read())
    return xmldict
#DEPRECATED BECAUSE NOW WE USE CANTERA CODE OR CANTERA OBJECTS
def makeReactionObjectsListFromCTMLasDict(CTMLasDict):
    reactionObjectsList = CTMLasDict["ctml"]["reactionData"]['reaction']
    return reactionObjectsList

#DEPRECATED BECAUSE NOW WE USE CANTERA CODE OR CANTERA OBJECTS  
    #NOte: does not have is_sticking_coefficientList and would need that to be useful.
def getReactionParametersFromReactionObjectsList(reactionObjectList):
    reactionIDsList = []
    reactionTypesList = []
    reactionEquationsList = []
    ArrheniusParametersList =[] #List of numpy arrays with A,b,E for each reaction.
    concentrationDependencesList =[]
    concentrationDependencesSpeciesList = []
    for reactionObject in reactionObjectList:
        reactionIDsList.append(reactionObject["@id"])
        try:
            reactionTypesList.append(reactionObject["@type"])
        except:
            reactionTypesList.append("None")
        #Now need to modify equation based on whether it's irreversible or not.
        reactionEquation = reactionObject["equation"]
        if reactionObject["@reversible"]=="no":
            reactionEquation = reactionEquation.replace("=]","=>")
        if reactionObject["@reversible"]=="yes":
            reactionEquation = reactionEquation.replace("=]","<=>")
        reactionEquationsList.append(reactionEquation)
        A = float(reactionObject['rateCoeff']['Arrhenius']['A'])
        b = float(reactionObject['rateCoeff']['Arrhenius']['b'])
        if reactionObject['rateCoeff']['Arrhenius']['E']['@units']=="J/mol":
            E = float(reactionObject['rateCoeff']['Arrhenius']['E']['#text']) #the xml is written in a way that this is necessary.
        else:
            E = 0.0
            print("Error: Currently, ony J/mol is supported for units of E. Fix parameters for this reaction: " + reactionObject['equation'])
        singleReactionArrheniusParametersArray = np.array([A,b,E])
        ArrheniusParametersList.append(singleReactionArrheniusParametersArray)
        #Now check for coverage parameters.
        if "coverage" in reactionObject['rateCoeff']['Arrhenius']:
            singleReactionconcentrationDependencesDict = reactionObject['rateCoeff']['Arrhenius']['coverage']
            speciesDependence = singleReactionconcentrationDependencesDict['@species']
            a = float(singleReactionconcentrationDependencesDict['a'])
            m = float(singleReactionconcentrationDependencesDict['m'])
            if singleReactionconcentrationDependencesDict['e']['@units']=="J/mol":
                e = float(singleReactionconcentrationDependencesDict['e']['#text']) #the xml is written in a way that this is necessary.
            else:
                e = 0.0
                print("Error: Currently, ony J/mol is supported for units of e. Fix coverage parameters for this reaction: " + reactionObject['equation'])
        else:#There is no coverage dependence, so will put 'NaN')
             a = float('nan')
             m = float('nan')
             e = float('nan')
             speciesDependence = 'None'
        concentrationDependencesList.append(np.array([a,m,e]))
        concentrationDependencesSpeciesList.append(speciesDependence)
    reactionIDsList = reactionIDsList
    reactionTypesList = reactionTypesList
    ArrheniusParametersArray = np.array(ArrheniusParametersList)
    concentrationDependencesArray = np.array(concentrationDependencesList)
    concentrationDependencesSpeciesList = concentrationDependencesSpeciesList
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, concentrationDependencesArray, concentrationDependencesSpeciesList
'''




def ArrheniusParameterAddedToInOnePhase(canteraModule, canteraPhaseObject, ArrheniusParametersOperandsArray, byProvidedReactionID = True):
    #ArrheniusParametersMultipliersArray should have the same form as the reactions_parameters_array. 
    reactions_parameters_array = ArrheniusParametersOperandsArray #just using same local variable for convenience.
    ct = canteraModule
    for individualreactions_parameters_array in reactions_parameters_array:      
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
        A_operand = float(individualreactions_parameters_array[3])
        b_operand = float(individualreactions_parameters_array[4])
        E_operand = float(individualreactions_parameters_array[5])
        a_operand = float(individualreactions_parameters_array[6])
        m_operand = float(individualreactions_parameters_array[7])
        e_operand = float(individualreactions_parameters_array[8])
        concentrationDependenceSpecies = str(individualreactions_parameters_array[9])

        if byProvidedReactionID == False: #optional.... 
                #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.

        #First get the existing object and its values so that they can be multiplied.
        existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
        existingA = float(existingReactionObject.rate.pre_exponential_factor)
        existingb = float(existingReactionObject.rate.temperature_exponent)
        existingE = float(existingReactionObject.rate.activation_energy)
        temporaryconcentrationDependenceSpeciesList = list(existingReactionObject.concentration_deps.keys())
        if len(temporaryconcentrationDependenceSpeciesList) > 0:
            concentrationDependenceExists = True
            concentrationDependenceSpecies = temporaryconcentrationDependenceSpeciesList[0]
            existinga = existingReactionObject.concentration_deps[concentrationDependenceSpecies][0]
            existingm = existingReactionObject.concentration_deps[concentrationDependenceSpecies][1]
            existinge = existingReactionObject.concentration_deps[concentrationDependenceSpecies][2]
        existingReactionObject.rate = ct.Arrhenius(float(existingA+A_operand), float(existingb+b_operand), float(existingE+E_operand))
        if concentrationDependenceExists == True:
            print("Warning: Coverage dependence modifiers not working yet")
            try:
                if a_operand+m_operand+e_operand != 0.0: #FIXME: Not working yet (because of Cantera side, not because of this script.)
                    tupleForconcentrationDependence = (existinga+a_operand ,existingm+m_operand , existinge+e_operand)
                    existingReactionObject.concentration_deps[concentrationDependenceSpecies]=tupleForconcentrationDependence
            except:
                pass
        canteraPhaseObject.modify_reaction(int(reactionID), existingReactionObject)           
               
def ArrheniusParametersMultiplierInOnePhase(canteraModule, canteraPhaseObject, ArrheniusParametersMultipliersArray, byProvidedReactionID = True):
    #ArrheniusParametersMultipliersArray should have the same form as the reactions_parameters_array. 
    reactions_parameters_array = ArrheniusParametersMultipliersArray #just using same local variable for convenience.
    ct = canteraModule
    for individualreactions_parameters_array in reactions_parameters_array:      
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
        A_multiplier = float(individualreactions_parameters_array[3])
        b_multiplier = float(individualreactions_parameters_array[4])
        E_multiplier = float(individualreactions_parameters_array[5])
        a_multiplier = float(individualreactions_parameters_array[6])
        m_multiplier = float(individualreactions_parameters_array[7])
        e_multiplier = float(individualreactions_parameters_array[8])
        concentrationDependenceSpecies = str(individualreactions_parameters_array[9])

        if byProvidedReactionID == False: #optional.... 
                #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.

        #First get the existing object and its values so that they can be multiplied.
        existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
        existingA = float(existingReactionObject.rate.pre_exponential_factor)
        existingb = float(existingReactionObject.rate.temperature_exponent)
        existingE = float(existingReactionObject.rate.activation_energy)
        temporaryconcentrationDependenceSpeciesList = list(existingReactionObject.concentration_deps.keys())
        if len(temporaryconcentrationDependenceSpeciesList) > 0:
            concentrationDependenceExists = True
            concentrationDependenceSpecies = temporaryconcentrationDependenceSpeciesList[0]
            existinga = existingReactionObject.concentration_deps[concentrationDependenceSpecies][0]
            existingm = existingReactionObject.concentration_deps[concentrationDependenceSpecies][1]
            existinge = existingReactionObject.concentration_deps[concentrationDependenceSpecies][2]
        existingReactionObject.rate = ct.Arrhenius(float(existingA*A_multiplier), float(existingb*b_multiplier), float(existingE*E_multiplier))
        if concentrationDependenceExists == True:
            print("Warning: Coverage dependence modifiers not working yet")
            try:
                if a_multiplier*m_multiplier*e_multiplier != 1.0: #FIXME: Not working yet (because of Cantera side, not because of this script.)
                    tupleForconcentrationDependence = (existinga*a_multiplier ,existingm*m_multiplier , existinge*e_multiplier)
                    existingReactionObject.concentration_deps[concentrationDependenceSpecies]=tupleForconcentrationDependence
            except:
                pass
        canteraPhaseObject.modify_reaction(int(reactionID), existingReactionObject)           

def populatePiecewiseCoverageDependence(simulation_settings_module, original_reactions_parameters_array, species_name, kineticParameterName, piecewise_coverage_intervals, modifiers_array):
    try:
        len(simulation_settings_module.piecewise_coverage_dependences)
    except: #if it has no length, then dictionary has not been created yet.
        simulation_settings_module.piecewise_coverage_dependences = {}
    try:
        len(simulation_settings_module.piecewise_coverage_dependences[species_name])
    except: #if it has no length, then dictionary has not been created yet.
        simulation_settings_module.piecewise_coverage_dependences[species_name] = {}
        simulation_settings_module.piecewise_coverage_dependences[species_name]['piecewise_coverage_intervals'] = piecewise_coverage_intervals
        simulation_settings_module.piecewise_coverage_dependences[species_name]['piecewise_kinetic_parameter_modifier_arrays']={}
    simulation_settings_module.original_reactions_parameters_array = original_reactions_parameters_array
    simulation_settings_module.piecewise_coverage_dependences[species_name]['piecewise_kinetic_parameter_modifier_arrays'][kineticParameterName] = np.array(modifiers_array)
    return simulation_settings_module.piecewise_coverage_dependences #Normally don't use the return. Normally just use the coverage dependences to populate the variables in the simulation_settings_module.

def getInterpolatedModifiersArray(existing_values,piecewise_coverage_dependences, species_name, species_coverage, kineticParameterKey, reaction_parameters_array_kinetic_parameter_index):
    all_coverages_modifiers_array = np.array(piecewise_coverage_dependences[species_name]['piecewise_kinetic_parameter_modifier_arrays'][kineticParameterKey],dtype="float") #This is for every coverage dependent value for this species for this kinetic parameter (i.e., coverages of 0 to 1.0 in bins)
    if len(all_coverages_modifiers_array) != len(existing_values):
        print("You have piecewise_coverage_dependence set to True in the simulation settings. However, the coverage dependence array for kinetic parameter " + kineticParameterKey + " has a function of the coverage of" + species_name  + "does not match the length of the kinetic parameters array. These lengths must match to use this feature.")
    all_piecewise_coverage_intervals = np.array(piecewise_coverage_dependences[species_name]["piecewise_coverage_intervals"], dtype="float") #This is the list of coverages to interpolate between to find out which index to extract from.
    index_of_relevant_interval_low = np.searchsorted(all_piecewise_coverage_intervals, species_coverage, side="right") - 1 #Prefer to find lower coverage, but can only get right side of interval by searchsorted. (side="left doesn't work usefully for us).
    index_of_relevant_interval_high = index_of_relevant_interval_low + 1 

     #If the coverage value is less than or equal to the lowest coverage defined, we use the lowest index's value.
    if species_coverage <= min(all_piecewise_coverage_intervals):
        indexOfMin = np.where(all_piecewise_coverage_intervals == min(all_piecewise_coverage_intervals) )#We don't assume the user gave an ascending list.
        interpolated_modifiers_array = np.array(all_coverages_modifiers_array[:,indexOfMin]).flatten() #What comes out is nested and so we flatten it.
     #If the coverage value is greater than or equal to the highest coverage defined, we use the highest index's value.
    elif species_coverage >= max(all_piecewise_coverage_intervals):
        indexOfMax = np.where(all_piecewise_coverage_intervals == max(all_piecewise_coverage_intervals) )#We don't assume the user gave an ascending list.
        interpolated_modifiers_array = np.array(all_coverages_modifiers_array[:,indexOfMax]).flatten()  #What comes out is nested and so we flatten it.
    else:
        #Normal cases... we will interpolate to find the right modifier value between the nearest coverage below and above (denoted "low" and "high").
        modifiers_array_low = np.array(all_coverages_modifiers_array[:,index_of_relevant_interval_low],dtype="float") # The ":" is for slicing across all reactions.
        modifiers_array_high = np.array(all_coverages_modifiers_array[:,index_of_relevant_interval_high],dtype="float") # The ":" is for slicing across all reactions.
        piecewise_interval_low = np.array(all_piecewise_coverage_intervals[index_of_relevant_interval_low],dtype="float")
        piecewise_interval_high = np.array(all_piecewise_coverage_intervals[index_of_relevant_interval_high],dtype="float")
        slopes_array = (modifiers_array_high-modifiers_array_low)/(piecewise_interval_high-piecewise_interval_low)
        intercepts_array = modifiers_array_high - 1.0*slopes_array*piecewise_interval_high #b = y - mx
        #now from interpolation y = mx+b means  modifiers_array = slopes_array*species_coverage+intercepts_array
        interpolated_modifiers_array = slopes_array*species_coverage+intercepts_array
    return interpolated_modifiers_array

def calculatePiecewiseCoverageDependentModifiedParametersArray(settingsModuleOrObject, species_names, species_coverages): #To use this feature, a module or object must be passed which has all of the arrays inside with the correct names. Typically, simulation_settings_module would be used.
    #The data structure here is important to understand know how this function works:
    #The settingsModuleOrObject is expected to have the original_reactions_parameters_array.
    #The settingsModuleOrObject is also expected to have the piecewise_coverage_dependences in it, which is nested. It has a structure of dictionaries inside such that arrays are pulled out like this:  
    #    piecewise_coverage_dependences[species_name][kineticParameterName] 
    original_reactions_parameters_array = settingsModuleOrObject.original_reactions_parameters_array #This feature requires the settingsModuleOrObject to already have the reactions_parameters_array assigned to it.
    modified_reactions_parameters_array = np.array(copy.deepcopy(original_reactions_parameters_array), dtype="str") #Normally the original shoudl already be a str based np array, but it could be a nested list. We need to make it a numpy array for slicing.
    piecewise_coverage_dependences = settingsModuleOrObject.piecewise_coverage_dependences
    #TODO: consider adding other parameters. For now, only A_values and E_values are supported, though of course others can be supported. In the case of A, the piecewise coverage dependence is a factor 10^modifier.
    for species_name in piecewise_coverage_dependences.keys(): #Looping across species for which there are modifications terms. #the piecewise_coverage_dependance object is a dictionary with species as keys. Within that, are dictionaries for each species' contributions to the coverage dependences. Note that this may not represent ALL species.
        species_index = species_names.index(species_name) #This is only the index of the species in the species_coverages list. The piecewise coverage dependences dictionary may have fewer species.
        species_coverage = float(species_coverages[species_index])
        for kineticParameterKey in piecewise_coverage_dependences[species_name]['piecewise_kinetic_parameter_modifier_arrays']: #looping across kinetic parameters for which there are coverage dependent modifier terms. #There could be all of the kinetic parameter keys present, or only some subset.
            if kineticParameterKey == "A":
                reaction_parameters_array_kinetic_parameter_index = 3 #This is the 'column' in the reaction parameters array where the kinetic parameter is kept.             

                #NOTE: we take the modified_reactions_parameters_array because this way modifications can stack.
                existing_parameter_values = np.array(modified_reactions_parameters_array[:,reaction_parameters_array_kinetic_parameter_index], dtype="float") #NOTE: we modify the current_reactions_parameters_array because this way modifications can stack.
                modifiers_array= getInterpolatedModifiersArray(modified_reactions_parameters_array,piecewise_coverage_dependences, species_name, species_coverage, kineticParameterKey, reaction_parameters_array_kinetic_parameter_index)
    
                #now to apply the modifiers...
                modified_values = existing_parameter_values*(10**modifiers_array) #In the case of A, we apply a coefficient of 10^modifiers_array_value.
                #updated the values in reaction parameters array, but make them into a string.
                modified_reactions_parameters_array[:,reaction_parameters_array_kinetic_parameter_index] = np.array(modified_values, dtype="str")
            if kineticParameterKey == "E":
                reaction_parameters_array_kinetic_parameter_index = 5  #This is the 'column' in the reaction parameters array where the kinetic parameter is kept.
                
                #NOTE: we take the modified_reactions_parameters_array because this way modifications can stack.
                existing_parameter_values = np.array(modified_reactions_parameters_array[:,reaction_parameters_array_kinetic_parameter_index], dtype="float") #NOTE: we modify the current_reactions_parameters_array because this way modifications can stack.
                modifiers_array= getInterpolatedModifiersArray(modified_reactions_parameters_array,piecewise_coverage_dependences, species_name, species_coverage, kineticParameterKey, reaction_parameters_array_kinetic_parameter_index)
                
                modified_values = existing_parameter_values+modifiers_array #For E, we add the modifiers.
                #updated the values in reaction parameters array, but make them into a string.
                modified_reactions_parameters_array[:,reaction_parameters_array_kinetic_parameter_index] = np.array(modified_values, dtype="str")    
    return modified_reactions_parameters_array 
        

def multiplyReactionsInOnePhase(canteraModule, canteraPhaseObject, reactions_parameters_array, rateMultipliersArray, byProvidedReactionID = True):
    ct = canteraModule
    for inputReactionIndex,individualreactions_parameters_array in enumerate(reactions_parameters_array):    
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        reactionEquation = str(individualreactions_parameters_array[2])
        if byProvidedReactionID == False: #optional.... 
            #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
            reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.
        canteraPhaseObject.set_multiplier(rateMultipliersArray[inputReactionIndex],int(reactionID)) #note that we are iterating across each individual reaction, but that the reactionID may not match the inputReactionIndex
    
def modifyReactionsInOnePhase(canteraPhaseObject, reactions_parameters_array, ArrheniusOnly = False, byProvidedReactionID = True, input_Ea_units = 'J/mol'):
    #if byProvidedReactionID is set to false, then the reaction ID is "looked up" from the reactions list in the canteraPhaseObject.
    #The "ArrheniusOnly = True" flag is because cantera previously did not support non-Arrhenius modify_reaction for surface reactions. (but it is now supported as of Cantera 2.6)
    #The "multiplierOnly = True" flag just modifies the existing pre-exponential by a specified factor.
    
    import cantera as ct
    # We needed to classify what type of phase we have, because surface reactions can only occur in interface objects, and pure gas phase reactions should happen in solution phase objects.
    if(str(type(canteraPhaseObject)))== "<class 'cantera.composite.Interface'>":
         phaseType = "Interface"
         inputReactionTypeNeeded = "InterfaceReaction" #this capitalization is for the python object type,  not the yaml_string capitalization.
         #outputReactionObjectType = "InterfaceReaction" 
    elif(str(type(canteraPhaseObject)))== "<class 'cantera.composite.Solution'>":
         phaseType = "Solution"
         inputReactionTypeNeeded = "Reaction"  #this capitalization is for the python object type,  not the yaml_string capitalization.
         #outputReactionObjectType = "Reaction" 
    else: #For now, we assume the only inputReaction type is a reaction.
         inputReactionTypeNeeded = "reaction"

    #Now add any reactions that match for this phase.
    for individualreactions_parameters_array in reactions_parameters_array:      
        individualreactions_parameters_array = np.array(individualreactions_parameters_array) #To make sure it is a mutable object.
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        #individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
        # The below comments are for informational purposes.
        # A = str(individualreactions_parameters_array[3])
        # b = str(individualreactions_parameters_array[4])
        # E = str(individualreactions_parameters_array[5]) #Note: this is now "Ea" in cantera 2.6
        # a = str(individualreactions_parameters_array[6])
        # m = str(individualreactions_parameters_array[7])
        # e = str(individualreactions_parameters_array[8]) #NOTE: This is now "E" in cantera 2.6.        

        #We need to check if the reaction has any surface phase species, to see if it is the right reaction type for this cantera phase.
        if '(S)' in reactionEquation:
            inputReactionType = "InterfaceReaction" #this capitalization is for the python object type,  not the yaml_string capitalization.
        else:
            inputReactionType = "Reaction" #this capitalization is for the python object type,  not the yaml_string capitalization.
        
        #The below block makes the actual modification to the cantera phase object.
        #It is necesary to make sure the inputReactionType matches for the phase being sought -- otherwise the matching reactionID won't exist, and also the attempt to modify a reaction in the wrong phase will crash.
        if inputReactionType == inputReactionTypeNeeded: #We will only make modifcations if the inputReactionType is the correct type for this phase.            
            if byProvidedReactionID == False: #optional.... 
                    #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                    reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.    
            #We need to change the units of the activation energy (index 5 in the array ) if the input_Ea_units are J/mol, because cantera has everything in j/kmol inside.
            if input_Ea_units.lower() == 'j/mol': #convert j/mol to j/kmol.
                individualreactions_parameters_array[5] = float(individualreactions_parameters_array[5])*1000.0
            if ArrheniusOnly == True: #If we are constrained to Arrhenius only, we will just modify individual reaction parameters. This case is unusual as of cantera 2.6
                A = str(individualreactions_parameters_array[3])
                b = str(individualreactions_parameters_array[4])            
                E = str(individualreactions_parameters_array[5])
                existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
                existingReactionObject.rate = ct.Arrhenius(float(A), float(b), float(E)) #If modifying cantera Arrhenius objects, need to convert from J/mol to J/kmol.
                canteraPhaseObject.modify_reaction(int(reactionID), existingReactionObject)     
            if ArrheniusOnly == False: #As of cantera 2.6, the This is the normal case. In this case, we modify the full reaction object                                                                                         
                #When we are not constrained to Arrhenius only, will make a new reaction object for each reaction using a yaml_string.
                #When we make the yaml_string for the reaction, we set "for_full_yaml" as false because we are modyfing an individual reaction.
                yaml_string = make_reaction_yaml_string(individualreactions_parameters_array, for_full_yaml = False, input_Ea_units="J/kmol") #put J/kmol because any conversion will have happened before now.
                from distutils.version import LooseVersion, StrictVersion
                if LooseVersion(ct.__version__) < LooseVersion('2.6'):
                    modifiedReactionObject = ct.Reaction.fromYaml(yaml_string , kinetics=canteraPhaseObject) #This is the correct way for cantera version 2.5        
                if LooseVersion(ct.__version__) >= LooseVersion('2.6'):
                    modifiedReactionObject = ct.Reaction.from_yaml(yaml_string , kinetics=canteraPhaseObject) #This is the correct way from version 2.6            
                canteraPhaseObject.modify_reaction(int(reactionID),modifiedReactionObject) 

def main():
    pass


if __name__ == '__main__':
    main()