import numpy as np



#To use CheKiPEUQ, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(x,a,b): #here x is a scalar or an array and "a" and "b" are constants for the equation.
    x =np.array(x)
    y = (x-a)**2 + b
    #print("line 9", y)
    return y


#Now we will make a wrapper for the simulation function, since CheKiPEUQ needs that.
x_values_for_data = 0  #This is just initializing the global value to avoid confusion (see below)
def simulation_function_wrapper(parametersArray):#this has a and b in it.
    global x_values_for_data  #A good way of doing things is to create
    a_given = parametersArray[0]
    b_given = parametersArray[1]
    y = simulationFunction(x_values_for_data, a_given, b_given) 
    return y
    #an alternatie simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 