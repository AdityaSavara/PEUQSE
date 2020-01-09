# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 16:10:05 2019

@author: fvs
"""
import sys
sys.path.insert(0, '/mumpce/')
import mumpce.Project as mumpceProject
import mumpce.solution as mumpceSolution
import numpy as np
import matplotlib
from matplotlib import cm #EAW 2020/01/07
matplotlib.use('Agg') #EAW 2020/01/07

#TODO: Make 2D parameter response surfaces like in the perspective figures. Can make it for each variable pair. Should have it as an option in the UserInput as True, False, or a list of pairs for which ones to make.
#Below is pseudo code to begin doing that. See more info at https://github.com/AdityaSavara/ODE-KIN-BAYES-SG-EW/issues/9
#active_params= ['a']
mumpceProjectObject = mumpceProject.Project()
mumpceProjectObject.model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': 'Parameter 0', 'parameter_value': 1.0},
 {'parameter_number': 1, 'parameter_name': 'Parameter 1', 'parameter_value': 1.0}, #The parameter values are not used in making these contour plots.
 {'parameter_number': 2, 'parameter_name': 'Parameter 2', 'parameter_value': 1.0},
 {'parameter_number': 3, 'parameter_name': 'Parameter 3', 'parameter_value': 1.0},
 {'parameter_number': 4, 'parameter_name': 'Parameter 4', 'parameter_value': 1.0},
 {'parameter_number': 5, 'parameter_name': 'Parameter 5', 'parameter_value': 1.0},
 {'parameter_number': 6, 'parameter_name': 'Parameter 6', 'parameter_value': 1.0}]) #This must be made into a numpy array.
mumpceProjectObject.active_parameters = np.array([0, 1, 2, 4, 6]) #this must be made into a numpy array.
#mumpceProjectObject.set_active_parameters = [0, 1, 2, 4, 6]
Posterior_mu_vector = np.array([-0.58888733,1.1200355, 0.00704044, -1.62385888,0.80439847]) #this must be made into a numpy array. #This will become solution.x
Posterior_cov_vector = np.array([[ 0.0148872,-0.01894579, -0.01047339,0.01325883,0.04734254],
 [-0.01894579,0.04284732, -0.00131389, -0.04801795, -0.04545703],
 [-0.01047339, -0.00131389,0.02343653,0.01588293, -0.05618226],
 [ 0.01325883, -0.04801795,0.01588293,0.08171972,0.00875017],
 [ 0.04734254, -0.04545703, -0.05618226,0.00875017,0.20669273]]) #This will become solution.cov. It does not need to be a numpy array, but we make it one for consistency.

mumpceSolutionsObject = mumpceSolution.Solution(Posterior_mu_vector, Posterior_cov_vector)
mumpceProjectObject.solution = mumpceSolutionsObject
mumpceProjectObject.pairsOfParameterIndices = [[0, 1], [1, 2],[3, 4]]
mumpceProjectObject.solution.mu_prior = np.array([-0.98888733,0.8200355, 0.01204044, -1.02385888,0.40439847])
mumpceProjectObject.solution.cov_prior = 10*Posterior_cov_vector
#np.array([[ 0.0148872,-0.01894579, -0.01047339,0.01325883,0.04734254],
# [-0.01894579,0.04284732, -0.00131389, -0.04801795, -0.04545703],
# [-0.01047339, -0.00131389,0.02343653,0.01588293, -0.05618226],
# [ 0.01325883, -0.04801795,0.01588293,0.08171972,0.00875017],
# [ 0.04734254, -0.04545703, -0.05618226,0.00875017,0.20669273]])

#This makes the figures as originally programmed. It assumes/requries things be normalized to 1.
mumpceProjectObject.solution.figure_name='mumpce_plots_alpha.png'
mumpceProjectObject.solution.num_pts = 500

mumpceProjectObject.plot_pdfs(mumpceProjectObject.pairsOfParameterIndices, contour_settings_custom = False)
cmap = cm.Reds


#I have expanded the code to allow more versatility an optional argument called contour_settings_custom.
#  It does not assume/require that things be normalized to 1, but requires the below variables to be populated.
mumpceProjectObject.solution.contour_axis_range = np.array([[-3.0,1.5],[-1.5,3.0],[-1.5,1.5],[-1.5,1.5],[-1,1],[-1,1],[-1,1]]) #this must be made into a numpy array.
mumpceProjectObject.solution.contour_resolution = np.array([0.01,0.01,0.01,0.01,0.01,0.01])
mumpceProjectObject.solution.axis_tick_spacing = np.array([1,1.5,1,1,1,1])
mumpceProjectObject.solution.contour_levels = np.exp((np.arange(-2,0,0.5) ** 2) * -1)
mumpceProjectObject.solution.figure_name='mumpce_plots_beta.png'
mumpceProjectObject.plot_pdfs(mumpceProjectObject.pairsOfParameterIndices, contour_settings_custom = True)

#Here is the graph after Eric has made the needed changes to have contour_settings_auto....
#Right now it produces fairly similar looking graphs. Later, it will change things more!
mumpceProjectObject.solution.figure_name='mumpce_plots_gamma.png'
mumpceProjectObject.plot_pdfs(mumpceProjectObject.pairsOfParameterIndices, contour_settings_auto= True)
