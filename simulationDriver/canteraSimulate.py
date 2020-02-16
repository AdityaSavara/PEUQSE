#These lines are to make sure that the module has access to canteraKineticsParametersParser
import os
import sys
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)

import cantera as ct
import numpy as np
import csv
from math import ceil
import canteraKineticsParametersParser 


def create_cti_and_SimulatePFRorTPRwithCantera(model_name, reactions_parameters_array, simulation_settings_module, cti_top_info_string = None, write_cti_to_file = False ):
    #The things that must be defined in advance are...
    # a) a model name.
    # b) The simulation settings are set in a python file which then becomes imported and an argument. This can then be changed by a script.
    # c) The reactions_parameters_array must be fed as an array or as a filename that points to a file that contains such an array with a header.
    # d) The cti_top_info_string. 
    cti_string, canteraPhases  = create_cti_and_cantera_phases(model_name=model_name, simulation_settings_module=simulation_settings_module, reactions_parameters_array=reactions_parameters_array, cti_top_info_string = cti_top_info_string, write_cti_to_file = write_cti_to_file)

    #Now import the simulation settings before running the simulation...
    concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
    simulatePFRorTPRwithCantera(model_name, canteraPhases['gas'], canteraPhases['surf'], simulation_settings_module)  
    return concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject
    
def create_cti_and_cantera_phases(model_name, reactions_parameters_array, simulation_settings_module, cti_top_info_string = None, write_cti_to_file = False ):
    #The things that must be defined in advance are...
    # a) a model name.
    # b) The simulation settings are set in a python file which then becomes imported and an argument. This can then be changed by a script.
    # c) The reactions_parameters_array must be fed as an array or as a filename that points to a file that contains such an array with a header.
    # d) The cti_top_info_string. 
    if type(reactions_parameters_array) == type("string"):
        reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    #In our example, the cti_top_info file is already made, and the functions below know where to find that file.
    
    #The below function makes a cti_file and then returns the filename.
    cti_string = canteraKineticsParametersParser.create_full_cti(model_name=model_name, reactions_parameters_array=reactions_parameters_array, cti_top_info_string = cti_top_info_string, write_cti_to_file = write_cti_to_file)
    
    canteraPhases = {}
    #NOTE: the reaction parameters are in the reaction objects. Currently, we cannot update the coverage dependence after cti file creation. That is why they must be created after the cti_file is made.
    # import the gas model and surface model.
    from distutils.version import LooseVersion, StrictVersion
    if LooseVersion(ct.__version__) < LooseVersion('2.5'):
        canteraPhases['gas'] = ct.Solution(source=cti_string, phaseid='gas')
        #canteraPhases['gas'].reactionsParametersArray = reactions_parameters_array
        # import the surface model
        canteraPhases['surf']  = ct.Interface(source=cti_string,phaseid='surf', phases=[canteraPhases['gas']]) #The word gas here is passing that same gas phase object in that was created above.
        #canteraPhases['surf'].reactionsParametersArray = reactions_parameters_array
    if LooseVersion(ct.__version__) >= LooseVersion('2.5'):
        canteraPhases['gas'] = ct.Solution(source=cti_string, name='gas')
        #canteraPhases['gas'].reactionsParametersArray = reactions_parameters_array
        # import the surface model
        canteraPhases['surf']  = ct.Interface(source=cti_string, name='surf', adjacent=[canteraPhases['gas']]) #The word gas here is passing that same gas phase object in that was created above.
        #canteraPhases['surf'].reactionsParametersArray = reactions_parameters_array
    #NOTE: we intentionally use "surf" rather than surface because in principle a person can make a model with more than 1 surface.
    #The choice of "surf" will make the variables easier to keep track of when it is time to make such a change.
    return cti_string, canteraPhases
    
#TODO: This can probably be generalized to run any simulation rather than only simulatePFRorTPRwithCantera.
def modify_reactions_and_SimulatePFRorTPRwithCantera(model_name, reactions_parameters_array, simulation_settings_module, canteraPhases, ArrheniusOnly = True, byProvidedReactionID = True):
    #The things that must be defined in advance are...
    # a) a model name.
    # b) The simulation settings are set in a python file which then becomes imported and an argument. This can then be changed by a script.
    # c) The reactions_parameters_array must be fed as an array or as a filename that points to a file that contains such an array with a header.
    # d) The cti_top_info_string. 
    # Note: for this function, what should be passed in is a modified version of eactions_parameters_array. Otherwise there is no point in using this function. 
    # If the original parameter array is going to be used, then simulatePFRorTPRwithCantera could be used directly.
    if type(reactions_parameters_array) == type("string"):
        reactions_parameters_array = np.genfromtxt(model_name + "_input_reactions_parameters.csv", delimiter=",", dtype="str", skip_header=1)
    #In our example, the cti_top_info file is already made, and the functions below know where to find that file.
    
    #The below function modifies the reactions in the gas phase.
    canteraKineticsParametersParser.modifyReactionsInOnePhase(canteraPhases['gas'], reactions_parameters_array, ArrheniusOnly = ArrheniusOnly, byProvidedReactionID = byProvidedReactionID)   
        #The below function modifies the reactions in the surface phase.
    canteraKineticsParametersParser.modifyReactionsInOnePhase(canteraPhases['surf'], reactions_parameters_array, ArrheniusOnly = ArrheniusOnly, byProvidedReactionID = byProvidedReactionID)
    
    #Now run the simulation...
    concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject = \
    simulatePFRorTPRwithCantera(model_name, canteraPhases['gas'], canteraPhases['surf'], simulation_settings_module)  
    return concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject
    
def simulatePFRorTPRwithCantera(model_name, canteraGasPhaseObject, canteraSurfacePhaseObject, simulation_settings_module):
    #simulation_settings_module_name must be a string if it is provided. The module itself should not be passed as an argument. This is intentional.
    canteraPhases ={}
    gas = canteraGasPhaseObject
    surf = canteraSurfacePhaseObject
    canteraPhases['gas'] = gas
    canteraPhases['surf'] = surf
       
    
    '''Now the code that handles *either* isothermal PFR or surface TPR. In future, will allow PFR TPR as well'''
    #We are going to model things as a flow reactor made of CSTRs, with no flow for the surface TPR case, following this example: https://cantera.org/examples/python/reactors/surf_pfr.py.html
    #it could probably also have been done as a flow reactor with no flow:  https://cantera.org/examples/python/surface_chemistry/catalytic_combustion.py.html
      

    
    '''START OF settings that users should change.'''
    flow_type = simulation_settings_module.flow_type
    T_gas_feed = simulation_settings_module.T_gas_feed
    T_surf = simulation_settings_module.T_surf
    velocity = simulation_settings_module.velocity
    reactor_cross_area = simulation_settings_module.reactor_cross_area
    cat_area_per_vol = simulation_settings_module.cat_area_per_vol
    porosity = simulation_settings_module.porosity
    length = simulation_settings_module.length
    P_gas = simulation_settings_module.P_gas
    gas_composition = simulation_settings_module.gas_composition
    t_step_size = simulation_settings_module.t_step_size
    t_final = simulation_settings_module.t_final
    NReactors = simulation_settings_module.NReactors
    print_frequency = simulation_settings_module.print_frequency
    heating_rate = simulation_settings_module.heating_rate
    surface_coverages =  simulation_settings_module.surface_coverages
    rtol = simulation_settings_module.rtol
    atol = simulation_settings_module.atol
    exportOutputs = simulation_settings_module.exportOutputs
    '''END OF settings that users should change.'''
    


    #Initiate concentrations output file and headers.
    if exportOutputs == True:
        concentrations_output_filename = model_name + "_output_concentrations.csv"
        outfile = open(concentrations_output_filename,'w')
        writer = csv.writer(outfile)
    concentrationsArrayHeaderList = ['Distance (m)', 'time(s)',  'T_gas (K)', 'T_surf (K)', 'P (atm)'] + \
                    gas.species_names + surf.species_names
    concentrationsArrayHeader = str(concentrationsArrayHeaderList)[1:-1] #The species names were imported when "surf" and "gas" objects were created.
    if exportOutputs == True:    
        writer.writerow(concentrationsArrayHeaderList)
    
    
    if flow_type == "Static":
        velocity = 0
        NReactors = 2 #In the surf_pfr example. For static, we only need 1 reactor, but the rlen formula below has a minimum value of 2. #FIXME
    
    num_t_steps = ceil(t_final / t_step_size) #rounds up.
    
    rlen = length/(NReactors-1) #This is each individual CSTR's length. #FIXME: Why is this not just Nreactors? Is it because of reservoirs?...
    #rvol is the reactor volume. We're modeling this as as a single CSTR with no flow.
    rvol = reactor_cross_area * rlen * porosity #Each individual CSTR gas volume reduced by porosity.
    # catalyst area in one reactor
    cat_area = cat_area_per_vol * rvol
    
    
    #Set the initial conditions for the gas and the surface.
    gas.TPX = T_gas_feed, P_gas, gas_composition
    surf.TP = T_surf, P_gas
    surf.X = surface_coverages
    
    mass_flow_rate = velocity * reactor_cross_area * gas.density  #linear velocity times area is volume/s, times kg/vol becomes kg/s. I think that for this example we neglect effects of surface adsorption regarding how much mass is in the gas phase?
    
    TDY = gas.TDY #Get/Set temperature [K] and density [kg/m^3 or kmol/m^3], and mass fractions. From: https://cantera.org/documentation/docs-2.4/sphinx/html/cython/thermo.html
    cov = surf.coverages #This is normalized coverages, built in from InerfacePhase class: https://cantera.github.io/docs/sphinx/html/cython/thermo.html#cantera.InterfacePhase
    #It would also be possible to type surface.site_density, giving [kmol/m^2] for surface phases. Also possible to use surf.set_unnormalized_coverages(self, cov) for cases when don't want to use normalized coverages.
    
    # create a new reactor
    gas.TDY = TDY #<-- If TDY was changing, and if original gas wanted to be kept, this could have been gas = copy.deepcopy(gas)?
    reactor = ct.IdealGasReactor(gas, energy='off')
    reactor.volume = rvol
    
    # create a reservoir to represent the reactor immediately upstream. Note
    # that the gas object is set already to the state of the upstream reactor
    upstream_of_CSTR = ct.Reservoir(gas, name='upstream') #A cantera reservoir never changes no matter what goes in/out: https://cantera.org/documentation/docs-2.4/sphinx/html/cython/zerodim.html#cantera.Reservoir
    
    # create a reservoir for the reactor to exhaust into. The composition of
    # this reservoir is irrelevant.
    downstream_of_CSTR = ct.Reservoir(gas, name='downstream') #A cantera reservoir never changes no matter what goes in/out: https://cantera.org/documentation/docs-2.4/sphinx/html/cython/zerodim.html#cantera.Reservoir
    
    #Note: these are reservoirs for individual CSTRs. It's used in a clever way in this example, as will be seen later.
    #I would prefer to call these reservoirs "feed" and "exhaust", and would consider to fill the downstream/exhaust with inert to make it clear that it is different. Or maybe make deep copies of gas called "feed_gas" and "exhaust_gas".
    
    # Add the reacting surface to the reactor. The area is set to the desired
    # catalyst area in the reactor.
    rsurf = ct.ReactorSurface(surf, reactor, A=cat_area)  #Here is where the catalyst site density gets used.  Note that cat_area is actually in meters now, even though site_density is defined as mol/cm^2 in the cti file.
    
    # The mass flow rate into the reactor will be fixed by using a
    # MassFlowController object.
    f_mfr = ct.MassFlowController(upstream_of_CSTR, reactor, mdot=mass_flow_rate) #mass flow controller makes flow that goes from the first argument to the second argument with mdot providing "The mass flow rate [kg/s] through this device at time t [s]."s
    
    # We need an outlet to the downstream reservoir. This will determine the
    # pressure in the reactor. The value of K will only affect the transient
    # pressure difference.
    e_mfr = ct.PressureController(reactor, downstream_of_CSTR, master=f_mfr, K=1e-5)  #This makes flow that goes from first argument to second argument with the mass flow controller "f_mfr" controlling the pressure here. K has units of kg/s/Pa times the pressure difference. this "v" that comes out is in same units as the mass flow rate, it's in kg/s.   https://cantera.org/documentation/docs-2.4/sphinx/html/cython/zerodim.html#flowdevice e_mfr is the exhaust mass flow rate.
    
    sim = ct.ReactorNet([reactor]) #This is normally a list of reactors that are coupled.
    sim.max_err_test_fails = 12 #Even after looking at docs, I'm not sure what this does. I think this is how many relative and absolute tolerance failures there can be in a single time step before the resolution of the ODE integrator does that step again with finer resolution? https://cantera.org/documentation/docs-2.4/sphinx/html/cython/zerodim.html#reactor-networks
    
    #Static versus PFR flag
    
    
    
    # set relative and absolute tolerances on the simulation
    sim.rtol = rtol
    sim.atol = atol
    gas_rates = [] #NOTE: These are just the rates from the surface phase. These are *not* including any homogeneous rates.
    surface_rates = [] #NOTE: These are just the rates from the surface phase. These are *not* including any homogeneous rates.
    sim_times = [] #NOTE: This is less useful if one uses advance_to_steady_state
    sim_dist = [] #This is the distance in the reactor. If one is using flow and advance_to_steady_state, then this is representative of the kinetics.
    concentrationsArray = []
    #Print some things out for before the simulation, these are basically headers.
    if print_frequency != None:
        print(concentrationsArrayHeader)
    
    if flow_type == "PFR":
        for n in range(NReactors): #iterate across the CSTRs.
            # Set the state of the reservoir to match that of the previous reactor
            gas.TDY = reactor.thermo.TDY #<-- setting gas to *most recent* state of the reactor.
            upstream_of_CSTR.syncState() #<-- this is an interesting trick. Once one CSTR has been simulated, we set the current values (the one from the just simulated CSTR) to be upstream value. This works because we're doing implicit ODE down a reactor length, and probably only gives the right answer for steady state kinetics. It also only works because we don't care what's happening downstream. We have a fixed K for the e_mfr, and also f_mfr.
            sim.reinitialize()
            sim.advance_to_steady_state() #<-- we advance the current CSTR to steady state before going to the next one.  For a tpr, we'd probably use advance(self, double t) instead. The syntax would be something like sim.advance(10.0) for 10 seconds. https://cantera.github.io/docs/sphinx/html/cython/examples/reactors_reactor1.html#py-example-reactor1-py
        
            dist = n * rlen # distance in m <-- this could have been defined above, and I would have... but it's probably good to have it down here to make it clear that only 'n" is being used rather than this value, to simulate the reactor.
            sim_dist.append(dist)
            sim_times.append(sim.time)
            gas_rates.append(surf.get_net_production_rates('gas'))
            surface_rates.append(surf.get_net_production_rates('surf'))
            # write the gas mole fractions and surface coverages vs. distance
            rowListOfStrings = [dist, sim.time, gas.T, surf.T, reactor.thermo.P/ct.one_atm] + \
                        list(gas.concentrations) + list(surf.coverages)
            if exportOutputs == True:
                writer.writerow(rowListOfStrings)
                concentrationsArray.append(np.array(rowListOfStrings))
            if print_frequency != None: 
                if not n % print_frequency: #This only prints every specified number of steps.
                    try:
                        with np.printoptions(precision=3):
                            print(np.array(rowListOfStrings))
                    except:
                            print(np.array(rowListOfStrings))
    
    if flow_type == "Static":
        T_surf_0 = T_surf
        #sim.set_max_time_step(t_step_size) #Cantera's sim.step normally advances in a variable way.
        for i_step in range(num_t_steps): #iterate across the CSTRs.
            time = i_step*t_step_size
            T_surf = T_surf_0 + time*heating_rate
            surf.TP = T_surf, P_gas
            if simulation_settings_module.piecewise_coverage_dependence == True:
                modified_reactions_parameters_array = canteraKineticsParametersParser.calculatePiecewiseCoverageDependentModifiedParametersArray(simulation_settings_module, surf.species_names, surf.coverages) #This feature requires the piecewise coverage dependence settings AND the reactions_parameters_array to already be inside the surf object **in advance**
                canteraKineticsParametersParser.modifyReactionsInOnePhase(surf, modified_reactions_parameters_array, ArrheniusOnly=True) #TODO: Right now this is constrainted to Arrhenius only because cantera does not yet allow modifyReactionsInOnePhase to do more (but it's planned to change in developers roadmap)                
            surf.advance_coverages(t_step_size)  #sim.advance(time) would not work. Changing T with time is not officially supported but happens to work with surf.advance_coverages. Supported way to change temperature during simulation for arbitrary reactors is to use custom integrator: https://cantera.org/examples/python/reactors/custom.py.html
            dist = 0.0
            sim_dist.append(dist)
            sim_times.append(time)
            gas_rates.append(surf.get_net_production_rates('gas'))
            surface_rates.append(surf.get_net_production_rates('surf'))        
            rowListOfStrings = [dist, time, gas.T, surf.T, reactor.thermo.P/ct.one_atm] + \
                    list(gas.X) + list(surf.coverages)
            if exportOutputs == True:
                writer.writerow(rowListOfStrings)
                concentrationsArray.append(np.array(rowListOfStrings))
            if print_frequency != None:
                if not i_step % print_frequency: #This only prints every specified number of steps.
                    try:
                        with np.printoptions(precision=3):
                            print(np.array(rowListOfStrings))
                    except:
                            print(np.array(rowListOfStrings))

           
    #Need to get things into a stackable state. Figured out from numpy shape that I needed to do at_least2D and transpose.
    sim_dist =  np.atleast_2d(sim_dist).transpose()
    sim_times =  np.atleast_2d(sim_times).transpose()
    gas_rates =  np.array(gas_rates)
    surface_rates =  np.array(surface_rates)    
    gasRatesArray = np.hstack((sim_dist,sim_times, gas_rates))
    surfaceRatesArray = np.hstack((sim_dist,sim_times, surface_rates))
    rates_all_array = np.hstack((sim_dist,sim_times, gas_rates, surface_rates))
    concentrationsArray = np.array(concentrationsArray)
    concentrationsArrayHeader = concentrationsArrayHeader
    gasRatesArrayHeader = 'dist(m), time(s),'+str(gas.species_names).replace("'","")[1:-1]
    surfaceRatesArrayHeader = 'dist(m),time(s),'+str(surf.species_names).replace("'","")[1:-1]
    rates_all_array_header = 'dist(m),time(s),'+str(gas.species_names).replace("'","")[1:-1]+"," + str(surf.species_names).replace("'","")[1:-1]
    canteraSimulationsObject = sim
    cantera_phase_rates = {"gas":gasRatesArray, "surf":surfaceRatesArray}
    cantera_phase_rates_headers = {"gas":gasRatesArrayHeader, "surf":surfaceRatesArrayHeader}
    
    if exportOutputs == True:
        np.savetxt(model_name + "_output_rates_all.csv", rates_all_array, delimiter=",", comments = '', header = rates_all_array_header)
        np.savetxt(model_name + "_output_rates_gas.csv", gasRatesArray, delimiter=",", comments = '', header = gasRatesArrayHeader )
        np.savetxt(model_name + "_output_rates_surf.csv", surfaceRatesArray, delimiter=",", comments = '', header = surfaceRatesArrayHeader)     
    if exportOutputs == True:    
        outfile.close()
    return concentrationsArray, concentrationsArrayHeader, rates_all_array, rates_all_array_header, cantera_phase_rates, canteraPhases, cantera_phase_rates_headers, canteraSimulationsObject