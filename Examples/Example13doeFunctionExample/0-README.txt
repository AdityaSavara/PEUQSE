In this example, two files are provided.

runfile_Example13_doe_function_meshgrid.py

and

runfile_Example13_doe_function_xyz.py

They each use the design of experiments function which is a class function for the parameter estimation (PE_object) class.



Below explains the variables etc. Check the UserInput file and other documentation as this might be out of date.
    
    UserInput.doe_settings['info_gains_matrices_array_format'] = 'meshgrid'  #<-- this can be 'meshgrid' or 'xyz'. It *must* be defined before calling the doe function.
    UserInput.doe_settings['independent_variable_grid_center'] = [500, 0.5]  #<-- this is the center of the conditions space.
    UserInput.doe_settings['independent_variable_grid_interval_size'] = [100, 0.1]  #<-- this how big of a step to take in each direction.
    UserInput.doe_settings['independent_variable_grid_num_intervals'] = [2,2] #This is the number of steps in each direction outward from center. So a 2 here gives 5 evaluations. A zero means we don't allow the condition to vary.
    
    UserInput.doe_settings['parameter_modulation_grid_interval_size'] = [1,1] #This defines the paramater modulation step sizes *relative* to their prior standard deviations. Use a non-zero value even for parameters that you will not vary.
    UserInput.doe_settings['parameter_modulation_grid_num_intervals'] = [1,1] #make the number of intervals (steps) in each dicrection. zero for a parameter that you don't want to vary
    
    #Make your PE_object as usual:
    PE_object = CKPQ.parameter_estimation(UserInput)
    
    #Call the design of experiments function, which will either do 'xyz' or 'meshgrid'.
    PE_object.doeParameterModulationCombinationsScanner()

    #The below larger object is returned.  the first index is the first combination of parameter modulations.
    PE_object.info_gains_matrices_array[0]


    #To get a single info gain plot, you can use the below, but you still have to fill out the doe_settings for the independent_variable_grid.
    PE_object.doeGetInfoGainMatrix()  