import numpy as np



#To use CheKiPEUQ, you can have a function, but you also need to make a function wrapper that takes *only* the parameters as a single vector.
def simulationFunction(x,a,b): #here x is a scalar or an array and "a" and "b" are constants for the equation.
    x =np.array(x)
    y = (x-a)**2 + b  #This is the same as d = (t-a)**2 + b
    #print("line 9", y)
    return y



#Now we will make a wrapper for the simulation function, since CheKiPEUQ needs that.
x_values_for_data = []  #This is just initializing the global value to avoid confusion (see below)
#Now to populate the global variable.
import observed_values_00 #I am using an import to show we can use the x_values associated with the observed data. Our wrapper should take only parameters and not x values.
x_values_for_data = observed_values_00.observed_data_x_values
def simulation_function_wrapper(parametersArray):#this has a and b in it.
    global x_values_for_data  #It is a good idea to create a global variable to pass in the x_values for your simulation.
    a_given = parametersArray[0] #a "given" just means this wrapper will simulate using whatever a value it receives.
    b_given = parametersArray[1] #b "given" just means this wrapper will simulate using whatever b value it receives.
    y = simulationFunction(x_values_for_data, a_given, b_given)  #an alternatie simpler syntax to unpack the parameters would be: simulationFunction(x_values_for_data, *parametersArray) 
    return y
    