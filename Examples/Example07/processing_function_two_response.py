def direct_parameters_to_observations(input):
    return input

def split_to_list(input):
    listAfterSplit = []
    for element in input:
        listAfterSplit.append(element)
    return listAfterSplit

def split_to_separated_lists(input):
    listAfterSplit = []
    for element in input:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    return listAfterSplit
