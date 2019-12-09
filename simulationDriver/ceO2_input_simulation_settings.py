#First add in a couple of variables just as unit converters... users can add more.
cm = 0.01 #This will convert things to m as required by cantera.
minute = 60.0 #This will convert things to seconds, as required by cantera.
atm = 101325.0 #This will convert atm to Pascal, as requred by cantera

'''START OF settings that users should change.'''
flow_type = "Static" #can be "PFR" or "Static". It does not only affect the next line, it affects the simulation.
T_gas_feed = 300 
T_surf = 100 #Initial temperature in K.
velocity = 30 * cm/minute #This must be in m/s. User can ignore this if using flow_type = "Static". It will be set to 0 below.
reactor_cross_area = 1.0 * cm**2  # this is our own variable. We are converting to meters.
cat_area_per_vol = 1000.0 / cm  # Catalyst particle surface area per unit volume
porosity = 0.3 #Catalyst bed porosity. <-- this will matter for getting our surface area.
length = 0.3 * cm  # Catalyst bed length or reactive area bed length.
P_gas = 1E-20* atm #We cannot have a pressure of 0. So for UHV, we just use something very small.
gas_composition = 'Acetaldehyde:1' #The composition cannot add up to 0, but does not have to add up to 1. Syntax like 'CH4:1, O2:1.5, AR:0.1'
t_step_size = 1.0 #This is the time resolution (s) to export. #This variable only matters for the "Static" case. Ignore it if one is using PFR to steady state. For a TPR experiment, a value of 1.0 is typically sufficient.
t_final = 500 #This is the final time to reach.  #This variable only matters for the "Static" case. Ignore it if one is using PFR to steady state.
NReactors = 201 # The PFR will be simulated by a chain of 'NReactors' stirred reactors. (CSTRs). Set it to 1 higher than what you need.
print_frequency = 10 #This is how often to print to screen, either in reactors (PFR) or in time (Static). Can be set to "None" or an integer. 
heating_rate = 2.0 #K/s. This variable is currently ignored if flow_type of PFR is being used, since PFR with TPR is not yet implemented.
surface_coverages = 'CeCation(S):0 Acetaldehyde1-Ce(S):0.25 Acetaldehyde2-Ce(S):0.25 OAnion(S):0.5' #This must be a list or string that is either like this: 'CeCation(S):0 Acetaldehyde-Ce(S):0 OAnion(S):2' or like this [0,1,2]. The sum will become normalized to 1, unless a species takes up two sites.
piecewise_coverage_dependence = False
rtol = 1.0e-9 #numerical relative tolerance.
atol = 1.0e-21 #numerical absolute tolerance.
exportOutputs = True #This exports results to files.
'''END OF settings that users should change.'''


