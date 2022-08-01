import UnitTesterSG

import sys

#We want runPytestDriver.py to give a pytest exitCode of 1 for travis CI purposes. That means we need to create an "actual" error when any test fails.
#The pytestDriver function of UnitTesterSG has an option to 'crash' with an error if "failWithError" is passed in.
#So the below lines of code are such that if the command line has 'python runPytestDriver.py failWithError"
if (len(sys.argv) <= 1):
    UnitTesterSG.pytestDriver(failWithError=False)
elif sys.argv[1].strip() == '':    
    UnitTesterSG.pytestDriver(failWithError=False)
elif sys.argv[1].strip().lower() == 'failwitherror':
    sys.argv.pop(1) #Need to remove the argument from sys.argv, because pytest also looks in the system arguments and if "failWithError" is still present, that will cause an undesired error.
    UnitTesterSG.pytestDriver(failWithError=True)