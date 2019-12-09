# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:44:31 2019

@author: Yurik
"""

import xmltodict
import numpy as np
import cantera as ct 
import copy
#This module is written for reactions_parameters_arrays in the following format:
#headerString = "reactionID,canteraReactionType,reactionEquation,A,b,E,a,m,e,coverageDependenceSpecies, is_sticking_coefficient"
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
        if len(np.shape(reactions_parameters_array))==1:
            piecewise_coverage_intervals = piecewise_coverage_intervals_all_reactions
        else:
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
        e = str(individualreactions_parameters_array[8])
        coverageDependenceSpecies = str(individualreactions_parameters_array[9])

        #Get the numbers that we need out.
        E_0 = float(E)
        
        if np.isnan(float(e)): #Need to make e into 0.0 if it's an nan. Otherwise the descending check will fail.
            e = 0.0
        g_slope = -1.0*float(e) #Note that in cantera e is the negative of g_slope! https://cantera.org/science/reactions.html
        print(reactionID, reactionEquation, e, g_slope)
        if individualReactionTypeReceived != "surface_reaction":
            passedArray[reactionIndex] = True  #if it's not a surface reaction, we are not supposed to check it.
        if individualReactionTypeReceived == "surface_reaction":
            passedArray[reactionIndex] = descendingLinearEWithPiecewiseOffsetCheckOneReaction(E_0, g_slope,piecewise_coverage_intervals, E_offsets_array)
            if verbose == True:
                if passedArray[reactionIndex] == False:
                    print("The reaction with ID of" + reactionID + "and equation of" + reactionEquation + "did not pass descendingLinearEWithPiecewiseOffsetCheckOneReaction" )
    if sum(passedArray )< len(passedArray): #For a boolean array, if the sum is less than the length, that means at least one thing remained false.
        return False
    elif sum(passedArray) == len(passedArray): #Else the sum must be equal and we return true.
        return True
    
def descendingLinearEWithPiecewiseOffsetCheckOneReaction(E_0, g_slope,piecewise_coverage_intervals, E_offsets_array):
    E_array =E_0 + g_slope*piecewise_coverage_intervals+E_offsets_array
    delta_E_array = np.diff(E_array)
    if np.any(delta_E_array > 0):
        return False
    else:
        return True

def makeCanteraReactionObjectsListFromFile(FileName):
    canteraReactionObjectsList = ct.Reaction.listFromFile(FileName)
    return canteraReactionObjectsList


def extractReactionParametersFromFile(InputFileName, OutputFilename = ""): #The input fileName must point to a cti file or an xml file. An optional argument of OutputFilename will export the reactionParamtersArray
    canteraReactionObjectsList = makeCanteraReactionObjectsListFromFile(InputFileName)
    reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList, is_sticking_coefficientList = getReactionParametersFromCanteraReactionObjectsList(canteraReactionObjectsList)
    outputAsNumpyArray = stackListsAndArrays([reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList, is_sticking_coefficientList])
    if OutputFilename != "":
        headerString = "reactionID,canteraReactionType,reactionEquation,A,b,E,a,m,e,coverageDependenceSpecies, is_sticking_coefficient"
        np.savetxt(OutputFilename,outputAsNumpyArray, fmt="%s", delimiter=",", comments='', header=headerString)
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList, is_sticking_coefficientList

def getReactionParametersFromCanteraReactionObjectsList(reactionObjectList):
    reactionIDsList = []
    reactionTypesList = []
    reactionEquationsList = []
    ArrheniusParametersList =[] #List of numpy arrays with A,b,E for each reaction.
    coverageDependencesList =[]
    coverageDependencesSpeciesList = []
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
            if len(reactionObject.coverage_deps) > 0:
                coverage_deps = True
            else:
                coverage_deps = False
        except:
                coverage_deps = False
        if coverage_deps == True:
            singleReactionCoverageDependencesDict = reactionObject.coverage_deps
            temporaryCoverageDependenceSpeciesList = list(reactionObject.coverage_deps.keys())
            if len(temporaryCoverageDependenceSpeciesList) > 0:
                coverageDependenceExists = True
                coverageDependenceSpecies = temporaryCoverageDependenceSpeciesList[0]
                a = reactionObject.coverage_deps[coverageDependenceSpecies][0]
                m = reactionObject.coverage_deps[coverageDependenceSpecies][1]
                e = float(reactionObject.coverage_deps[coverageDependenceSpecies][2])/1000 #Annoyingly, Cantera puts Energies out in J/kmol.   
        else:#There is no coverage dependence, so will put 'NaN')
             a = float('nan')
             m = float('nan')
             e = float('nan')
             coverageDependenceSpecies = 'None'
        coverageDependencesList.append(np.array([a,m,e]))
        coverageDependencesSpeciesList.append(coverageDependenceSpecies)
        is_sticking_coefficientList.append(reactionObject.is_sticking_coefficient)
    reactionIDsList = reactionIDsList
    reactionTypesList = reactionTypesList
    ArrheniusParametersArray = np.array(ArrheniusParametersList)
    coverageDependencesArray = np.array(coverageDependencesList)
    coverageDependencesSpeciesList = coverageDependencesSpeciesList
    is_sticking_coefficientList = is_sticking_coefficientList
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList, is_sticking_coefficientList
    
def make_reaction_cti_string(individualreactions_parameters_array):
    #input should be a an array of strings: 
    #reactionID	canteraReactionType	reactionEquation	A	b	E	a	m	e	coverageDependenceSpecies	 is_sticking_coefficient
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
    coverageDependenceSpecies = str(individualreactions_parameters_array[9])
    is_sticking = str(individualreactions_parameters_array[10])
    if is_sticking.capitalize() == "False": #considered using distutils.util.strtobool(is_sticking) but decided that would slow the program down.
        ArrheniusString = "Arrhenius"
    if is_sticking.capitalize() == "True":
        ArrheniusString = "stick"
        
    #Now make the reaction_cti strings...
    if individualReactionTypeReceived.lower() == "reaction":
        reaction_cti_string = '''reaction('{0}',
        {1}({2}, {3}, {4}])'''.format(reactionEquation, ArrheniusString, A,b,E)
    if individualReactionTypeReceived.lower() == "surface_reaction":
        if coverageDependenceSpecies.capitalize() == "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
            reaction_cti_string= '''surface_reaction("{0}",
            {1}({2}, {3}, {4}))'''.format(reactionEquation,ArrheniusString,A,b,E) 
                #Like the below example.
                #        R3 = ct.InterfaceReaction.fromCti('''surface_reaction('O + H2 <=> H + OH',
                #        [3.87e1, 2.7, 2.619184e7])''')                      
        if coverageDependenceSpecies.capitalize() != "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
            reaction_cti_string = '''surface_reaction("{0}",
            {1}({2}, {3}, {4},
            coverage = ['{5}', {6}, {7}, {8}]))'''.format(reactionEquation,ArrheniusString,A,b,E,coverageDependenceSpecies,a,m,e)
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
    #reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList, is_sticking_coefficientList = exportReactionParametersFromFile("ceO2_cti_full_existing.cti", "ceO2_input_reactions_parameters.csv")
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
        coverageDependenceSpecies = str(individualreactions_parameters_array[9])
        is_sticking = str(individualreactions_parameters_array[10])

        #For all reactions, we can modify the A,b,E parameters.
        listOfModifiableParameterIndices = listOfModifiableParameterIndices + [ [reactionIndex,3],[reactionIndex,4],[reactionIndex,5] ]
        #Now check for other parameters...
        if individualReactionTypeReceived == 'surface_reaction': #This can include adsorption, desorption, and surface reactions.
            if coverageDependenceSpecies != 'None': #If this is not none, then that means the a, m, and e can also be modified.
                listOfModifiableParameterIndices = listOfModifiableParameterIndices + [ [reactionIndex,6],[reactionIndex,7],[reactionIndex,8] ]
                #Currently, we cannot modify the coverageDependenceSpecies.
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
    coverageDependencesList =[]
    coverageDependencesSpeciesList = []
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
            singleReactionCoverageDependencesDict = reactionObject['rateCoeff']['Arrhenius']['coverage']
            speciesDependence = singleReactionCoverageDependencesDict['@species']
            a = float(singleReactionCoverageDependencesDict['a'])
            m = float(singleReactionCoverageDependencesDict['m'])
            if singleReactionCoverageDependencesDict['e']['@units']=="J/mol":
                e = float(singleReactionCoverageDependencesDict['e']['#text']) #the xml is written in a way that this is necessary.
            else:
                e = 0.0
                print("Error: Currently, ony J/mol is supported for units of e. Fix coverage parameters for this reaction: " + reactionObject['equation'])
        else:#There is no coverage dependence, so will put 'NaN')
             a = float('nan')
             m = float('nan')
             e = float('nan')
             speciesDependence = 'None'
        coverageDependencesList.append(np.array([a,m,e]))
        coverageDependencesSpeciesList.append(speciesDependence)
    reactionIDsList = reactionIDsList
    reactionTypesList = reactionTypesList
    ArrheniusParametersArray = np.array(ArrheniusParametersList)
    coverageDependencesArray = np.array(coverageDependencesList)
    coverageDependencesSpeciesList = coverageDependencesSpeciesList
    return reactionIDsList, reactionTypesList, reactionEquationsList, ArrheniusParametersArray, coverageDependencesArray, coverageDependencesSpeciesList
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
        coverageDependenceSpecies = str(individualreactions_parameters_array[9])

        if byProvidedReactionID == False: #optional.... 
                #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.

        #First get the existing object and its values so that they can be multiplied.
        existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
        existingA = float(existingReactionObject.rate.pre_exponential_factor)
        existingb = float(existingReactionObject.rate.temperature_exponent)
        existingE = float(existingReactionObject.rate.activation_energy)
        temporaryCoverageDependenceSpeciesList = list(existingReactionObject.coverage_deps.keys())
        if len(temporaryCoverageDependenceSpeciesList) > 0:
            coverageDependenceExists = True
            coverageDependenceSpecies = temporaryCoverageDependenceSpeciesList[0]
            existinga = existingReactionObject.coverage_deps[coverageDependenceSpecies][0]
            existingm = existingReactionObject.coverage_deps[coverageDependenceSpecies][1]
            existinge = existingReactionObject.coverage_deps[coverageDependenceSpecies][2]
        existingReactionObject.rate = ct.Arrhenius(float(existingA+A_operand), float(existingb+b_operand), float(existingE+E_operand))
        if coverageDependenceExists == True:
            print("Warning: Coverage dependence modifiers not working yet")
            try:
                if a_operand+m_operand+e_operand != 0.0: #FIXME: Not working yet (because of Cantera side, not because of this script.)
                    tupleForCoverageDependence = (existinga+a_operand ,existingm+m_operand , existinge+e_operand)
                    existingReactionObject.coverage_deps[coverageDependenceSpecies]=tupleForCoverageDependence
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
        coverageDependenceSpecies = str(individualreactions_parameters_array[9])

        if byProvidedReactionID == False: #optional.... 
                #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.

        #First get the existing object and its values so that they can be multiplied.
        existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
        existingA = float(existingReactionObject.rate.pre_exponential_factor)
        existingb = float(existingReactionObject.rate.temperature_exponent)
        existingE = float(existingReactionObject.rate.activation_energy)
        temporaryCoverageDependenceSpeciesList = list(existingReactionObject.coverage_deps.keys())
        if len(temporaryCoverageDependenceSpeciesList) > 0:
            coverageDependenceExists = True
            coverageDependenceSpecies = temporaryCoverageDependenceSpeciesList[0]
            existinga = existingReactionObject.coverage_deps[coverageDependenceSpecies][0]
            existingm = existingReactionObject.coverage_deps[coverageDependenceSpecies][1]
            existinge = existingReactionObject.coverage_deps[coverageDependenceSpecies][2]
        existingReactionObject.rate = ct.Arrhenius(float(existingA*A_multiplier), float(existingb*b_multiplier), float(existingE*E_multiplier))
        if coverageDependenceExists == True:
            print("Warning: Coverage dependence modifiers not working yet")
            try:
                if a_multiplier*m_multiplier*e_multiplier != 1.0: #FIXME: Not working yet (because of Cantera side, not because of this script.)
                    tupleForCoverageDependence = (existinga*a_multiplier ,existingm*m_multiplier , existinge*e_multiplier)
                    existingReactionObject.coverage_deps[coverageDependenceSpecies]=tupleForCoverageDependence
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
    return simulation_settings_module.piecewise_coverage_dependences #Normally don't use the return. Normally just use it to populate the variables in the simulation_settings_module.

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
    
def modifyReactionsInOnePhase(canteraPhaseObject, reactions_parameters_array, ArrheniusOnly = True, byProvidedReactionID = True):
    #if byProvidedReactionID is set to false, then the reaction ID is "looked up" from the reactions list in the canteraPhaseObject.
    #The "ArrheniusOnly = True" flag is because cantera currently does not support non-Arrhenius modify_reaction for surface reactions.
    #The "multiplierOnly = True" flag just modifies the existing pre-exponential by a specified factor.
    
    import cantera as ct
    #First classify what type of phase we have.
    if(str(type(canteraPhaseObject)))== "<class 'cantera.composite.Interface'>":
        phaseType = "Interface"
        inputReactionTypeNeeded = "surface_reaction"
        #outputReactionObjectType = "InterfaceReaction" 
    if(str(type(canteraPhaseObject)))== "<class 'cantera.composite.Solution'>":
        phaseType = "Solution"
        inputReactionTypeNeeded = "reaction"
        #outputReactionObjectType = "Reaction" 

    #Now add any reactions that match for this phase.
    for individualreactions_parameters_array in reactions_parameters_array:      
        reactionID = str(int(individualreactions_parameters_array[0])-1) #we need to subtract one for how cantera expects the reaction IDs when modifying.
        individualReactionTypeReceived = str(individualreactions_parameters_array[1])
        reactionEquation = str(individualreactions_parameters_array[2])
        A = str(individualreactions_parameters_array[3])
        b = str(individualreactions_parameters_array[4])
        E = str(individualreactions_parameters_array[5])
        a = str(individualreactions_parameters_array[6])
        m = str(individualreactions_parameters_array[7])
        e = str(individualreactions_parameters_array[8])
        coverageDependenceSpecies = str(individualreactions_parameters_array[9])

        if byProvidedReactionID == False: #optional.... 
                #byProvidedReactionID is set to False, but this only works if the cantera model has the equation written exactly the same.
                reactionID = canteraPhaseObject.reaction_equations().index(reactionEquation) #don't add 1 because reaction index is not python array indexing.

        if ArrheniusOnly == True:
            if individualReactionTypeReceived == inputReactionTypeNeeded: #ignore any case that does not match what's needed for this phase.
                existingReactionObject = canteraPhaseObject.reactions()[int(reactionID)]
                existingReactionObject.rate = ct.Arrhenius(float(A), float(b), float(E)*1000) #If modifying cantera Arrhenius objects, need to convert from J/mol to J/kmol.
                canteraPhaseObject.modify_reaction(int(reactionID), existingReactionObject)     
        if ArrheniusOnly == False:
            if individualReactionTypeReceived == inputReactionTypeNeeded: #ignore any case that does not match what's needed for this phase.
                if phaseType == "Solution":
                    cti_string = '''reaction('{0}', [{1}, {2}, {3}])'''.format(reactionEquation, A,b,E)
                    modifiedReactionObject = ct.Reaction.fromCti(cti_string)
                if phaseType == "Interface":
                    if coverageDependenceSpecies == "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
                        cti_string= '''surface_reaction("{0}",
                            [{1}, {2}, {3}])'''.format(reactionEquation,A,b,E) 
                        modifiedReactionObject = ct.InterfaceReaction.fromCti(cti_string)
                            #Like the below example.
                            #        R3 = ct.InterfaceReaction.fromCti('''surface_reaction('O + H2 <=> H + OH',
                            #        [3.87e1, 2.7, 2.619184e7])''')                      
                    if coverageDependenceSpecies != "None": #FIXME: Not working yet (because of Cantera side, not because of this script.)
                        print("Warning: Coverage dependence modifiers not working yet")
                        cti_string = '''surface_reaction("{0}",Arrhenius({1}, {2}, {3},
                                            coverage = ['{4}', {5}, {6}, {7}]))'''.format(reactionEquation,A,b,E,coverageDependenceSpecies,a,m,e)
                        modifiedReactionObject = ct.InterfaceReaction.fromCti(cti_string)                    
                            #Like the below example.
                            #        R5 = ct.InterfaceReaction.fromCti('''surface_reaction( "CH4 + PT(S) + O(S) => CH3(S) + OH(S)",
                            #                  Arrhenius(5.0e18, 0.7, 2000.0,
                            #                            coverage = ['O(S)', 0.0, 0.0, 8000]))''')
               
                #now modfiy the reaction in the actual model:
                canteraPhaseObject.modify_reaction(int(reactionID),modifiedReactionObject) 


def main():
    pass


if __name__ == '__main__':
    main()