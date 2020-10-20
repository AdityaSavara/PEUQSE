# -*- coding: utf-8 -*-
"""
There are two ways to use UnitTesterSG: either with a known expected response, or simply comparing to a stored response.
Here, we will show a case with a known (e.g., analytical or otherwise calculated) expected response.
You may copy this file and modify it to make your own test. Just name your file test_N where N is an integer.

#NOTE: WHEN YOU RUN THIS FILE, IT WILL SAY THE **EXPECTED RESULT** MATCHES BUT THAT THE **EXPECTED RESULT STRING** DOES NOT MATCH.
#IT IS PERFECTLY FINE THAT THE RESULT STRING DOES NOT MATCH. THAT IS A TYPICAL SITUATION AND A FEATURE. IT DOES NOT MEAN THE TEST FAILED.
"""

#import the functions from UnitTesterSG
import UnitTesterSG as ut

#The below lines are typical code. There is no need to modify them.
#get the suffix argument for check_results
suffix = ut.returnDigitFromFilename(__file__)
#prefix. Make this '' if you do not want any prefix.
prefix = ''


"""
#For this test, we ***know*** what result to expect.
"""
expectedFirstPart = [1.00, 3.95]
expectedSecondPart = [0.90, 3.95]
expectedResult = (expectedFirstPart,expectedSecondPart) #We are using a tuple, but it this could have been a list.

ut.set_expected_result(expectedResult,expected_result_str=str(expectedResult), prefix=prefix,suffix=suffix) #This is the typical syntax if you want to force an analytical result for your test.


"""
#Calculate your actual outputs. 
"""
import runfile_Example07a
outputFirstPart = runfile_Example07a.PE_object.map_parameter_set
outputSecondPart = runfile_Example07a.PE_object.mu_AP_parameter_set
actualResult = (outputFirstPart, outputSecondPart)

"""Put your actual result into the resultObj variable."""
#put this in the resultObject
resultObj = actualResult

#String must be provided provided. Make it '' if you do not want to use a result string.
resultStr = str(resultObj)


"""Normally you will never touch the below lines except to change or make tolerances."""
relativeTolerance = 2.0E-2
absoluteTolerance = 2.0E-2


#this is so that pytest can do UnitTesterSG tests.
def test_pytest(): #note that it cannot have any required arguments for pytest to use it, and that it is using variables that are defined above in the module.
    ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance)
    
    
if __name__ == "__main__":
   #This is the normal way of using the UnitTesterSG module, and will be run by UnitTesterSG or by running this test file by itself.
   ut.doTest(resultObj, resultStr, prefix=prefix,suffix=suffix, allowOverwrite = False, relativeTolerance=relativeTolerance, absoluteTolerance=absoluteTolerance) 
