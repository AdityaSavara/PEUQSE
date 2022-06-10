# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:02:29 2020

@author: fvs
"""

import sys; sys.path.append('../../');  import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput

import processing_functions_tpd_odeint_two_site_NineParameters

observed_data_Filename = 'ExperimentalDataAcetaldehydeTPDCeO2111MullinsTruncatedConstantErrors005.csv'
times, responses_observed, observedResponses_uncertainties = processing_functions_tpd_odeint_two_site_NineParameters.import_integrals_settings(observed_data_Filename)

paramsToTest = [1.00528604e+00, 3.35978420e-04, 9.84034693e-01, 8.58890662e+01,
 3.25868789e+01, 1.26661158e+01 ,1.96633345e+01 ,3.48339352e-01,
 3.42164921e-01]
integrated_Desorption = processing_functions_tpd_odeint_two_site_NineParameters.TPR_integerated_simulationFunctionWrapperNineParameters(paramsToTest)
print(integrated_Desorption)
