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


def split_to_separated_lists_plus_1h(input):
    listAfterSplit = []
    for element in input:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    listAfterSplitPlus1 = listAfterSplit #initializing new variable name.
    listAfterSplitPlus1[1].append(0)
    return listAfterSplitPlus1


def split_to_separated_lists_plus_1i(input):
    listAfterSplit = []
    for element in input:
        listAfterSplit.append([element]) #Note the extra bracket inside to create nesting.
    listAfterSplitPlus1 = listAfterSplit #initializing new variable name.
    listAfterSplitPlus1[0].append(0)
    return listAfterSplitPlus1
