import sys; sys.path.append('../../')
sys.path.append('../../../../');  import PEUQSE as PEUQSE                                                               
#import PEUQSE.UserInput as UserInput


import dill
PE_object = dill.load(open("PE_object_dill.pkl", 'rb'))
PE_object.createAllPlots()
