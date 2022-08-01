# -*- coding: utf-8 -*-
"""
There are two ways to use UnitTesterSG: either with a known expected response, or simply comparing to a stored response.
This file is an example and template for a case where we can't simply 'provide' a solution so must compare to an **existing from before** output.

TO USE THIS TEMPLATE/EXAMPLE, YOU MUST RUN THIS FILE ONE TIME TO INITIATE. THE FIRST TIME THIS FILE IS RUN, IT WILL SAY THAT THE CALCULATED RESULT AND EXPECTED RESULT DO NOT MATCH -- THAT IS BECAUSE THE EXPECTED RESULT DOES NOT YET EXIST. DURING RUNTIME, THE SECOND QUESTION THE PROGRAM ASKS WILL BE IF YOU WANT TO OVERRWRITE OR CREATE EXPECTED RESULTS FILES: CHOOSE "Y"  AFTER DOING SO, THE EXPECTED RESULT HAS BEEN SET.  RUN THIS FILE AGAIN. SINCE YOU ARE RUNNING ON THE SAME COMPUTER, THE RESULT WILL BE THE SAME AND THE TEST WILL PASS. THE RESULT HAS BEEN STORED, SO THE TEST DIRECTORY CAN NOW BE INCLUDED IN A REPOSITORY OR COPIED TO ANOTHER COMPUTER FOR TESTING.

IF YOU WISH TO RESET THE STORED EXPECTED RESULTS TO NOT BEING SET YET, DELETE THE FILES THAT BEGIN WITH THE WORD EXPECTED.

You may copy this file and modify it to make your own test. Just name your file test_N where N is an integer.
"""

#These "sys" lines are mainly because this are standard lines in our examples. Normally, you would not include these three lines.
import sys
sys.path.insert(1, ".\\lib")
sys.path.insert(1, "..")



#import the functions from UnitTesterSG
import UnitTesterSG as ut

#The below lines are typical code. There is no need to modify them.
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''


"""
#This file is an example/template for when we ***don't have an analytical result*** but we know our function is working.
#We know the function is working during template distribution because we are just using the test 12 example.
In this template, we ***will not*** use the "set_expected_result" command. So we are commenting out the below lines, and will go directly to using the function to create an actual output.
"""
import sys; sys.path.insert(0, '../../');  import PEUQSE as PEUQSE
import PEUQSE
import numpy as np
import runfile_for_unit_test_parallel_doe_control #This will run the file, given how it's structured.
expectedResult = PEUQSE.unpickleAnObject("runfile_for_unit_test_parallel_doe_control")

# # # #input for the unit that will be tested
# # # input = 4
#expectedFirstPart = runfile_for_unit_test_parallel_doe_a.PE_object.info_gains_matrices_array
# # # expectedSecondPart = [32,64]
# # # expectedResult = (expectedFirstPart,expectedSecondPart) #We are using a tuple, but it this could have been a list.

ut.set_expected_result(expectedResult,expected_result_str=str(expectedResult), prefix=prefix,suffix=suffix) #This is the typical syntax if you want to force an analytical result for your test.


"""
#Calculate our function outputs (actual results). We can use functions from another module in this section.
"""
import os
try:
    os.remove("runfile_for_unit_test_parallel_doe_conditions_exploration.pkl")
except:
    pass
try:
    os.system("mpiexec -n 10 python runfile_for_unit_test_parallel_doe_conditions_exploration.py")  #This will run the file, given how it's structured.
except:
    os.system("mpirun -n 10 python runfile_for_unit_test_parallel_doe_conditions_exploration.py")  #This will run the file, given how it's structured.
actualResult = PEUQSE.unpickleAnObject("runfile_for_unit_test_parallel_doe_conditions_exploration")

input = None
#outputs with the function being tested using the input

"""We put our actual result into the resultObj variable."""
#put this in the resultObject
resultObj = actualResult

#String must be provided provided. Make it '' if you do not want to use a result string.
resultStr = str(resultObj)


"""We set our tolerances. There can be some rounding when the tolerances get checked, so they are not exact."""
relativeTolerance = 5.0E-2
absoluteTolerance = 5.0E-2


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
"""#For any individual test, after finishing getting it working, set allowOverwrite to False in the line below calling doTest if you want to skip UnitTesterSG from stopping to notify user when results match but result strings don't. """        
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
