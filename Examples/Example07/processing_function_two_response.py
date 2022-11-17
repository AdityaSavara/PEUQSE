import numpy as np
def direct_parameters_to_observations(inputArray):
    return inputArray

def split_to_list(inputArray):
    listAfterSplit = []
    for element in inputArray:
        listAfterSplit.append(element)
    return listAfterSplit

def split_to_separated_lists(inputArray):
    listAfterSplit = []
    for element in inputArray:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    return listAfterSplit

def split_to_separated_lists_None_and_Nan_filtered(inputArray):
    listAfterSplit = []
    if inputArray[0]*inputArray[1]> 2.5:
        return None
        #listAfterSplit.append([[float('nan')],[float('nan')]])
    elif inputArray[1] < 2.5:
        listAfterSplit.append([[float('nan')],[0]])
    else:
        for element in inputArray:
            listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    listAfterSplit = np.array(listAfterSplit, dtype="object")
    return listAfterSplit
        
def split_to_two_separated_lists(inputArray):
    listAfterSplit = [inputArray[0:2], inputArray[2:4]]
    return listAfterSplit


def split_to_separated_lists_plus_1h(inputArray):
    listAfterSplit = []
    for element in inputArray:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    listAfterSplitPlus1 = listAfterSplit #initializing new variable name.
    listAfterSplitPlus1[1].append(0)
    listAfterSplitPlus1 = np.array(listAfterSplitPlus1, dtype="object")
    return listAfterSplitPlus1


def split_to_separated_lists_plus_1i(inputArray):
    listAfterSplit = []
    for element in inputArray:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    listAfterSplitPlus1 = listAfterSplit #initializing new variable name.
    listAfterSplitPlus1[0].append(0)
    listAfterSplitPlus1 = np.array(listAfterSplitPlus1, dtype="object")
    return listAfterSplitPlus1
