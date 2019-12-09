import cantera as ct
import cantera.ck2cti as ck2cti
import canteraSimulate
import numpy as np
import canteraKineticsParametersParser 




    

def main():
    model_name = "ceO2"
    reactions_parameters_array = np.genfromtxt(model_name + "_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    listOfModifiableParameterIndices = canteraKineticsParametersParser.findModifiableParameterIndices(reactions_parameters_array)
    print("For this example, we get 9 parameters:")
    print(listOfModifiableParameterIndices)
    #[[0, 3], [0, 4], [0, 5], [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [1, 8]]   
    print("Now we make a silly example, just to demonstrate.")
    print("Here is before modification")
    print(reactions_parameters_array)
    modified_reactions_parameters_array= canteraKineticsParametersParser.modifyAllAdjustableReactionParametersArrayValues(reactions_parameters_array, [1,2,3,4,5,6,7,8,9])
    print("Here is after modification")
    print(modified_reactions_parameters_array)


if __name__ == '__main__':
    main()    