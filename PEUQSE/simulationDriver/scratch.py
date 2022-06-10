# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import cantera
#cantera.cti2yaml.convert(text='surface_reaction("Acetaldehyde + CeCation(S) => Acetaldehyde1-Ce(S)", stick(1, 0, 2000))', output_name="Testing.yaml")
import ceO2_input_simulation_settings
import canteraKineticsParametersParser
import canteraSimulate

reactions_parameters_array = "ceO2_input_reactions_parameters.csv"
yaml_string, canteraPhases     = canteraSimulate.create_yaml_and_cantera_phases('ceO2', reactions_parameters_array, ceO2_input_simulation_settings)

print(canteraPhases)


