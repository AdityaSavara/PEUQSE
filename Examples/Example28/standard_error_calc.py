# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:54:04 2022

@author: NWT

This program estimates the uncertainty associated with the diffusion 
coefficient for Cr from weight loss data.
"""
from scipy.stats import sem
import math


####################################Weight Loss################################

# Weight loss (mg/cm^2) values estimated from data provided in Dr. Zheng's paper
weight_loss = [0.45625, 0.5499]

# Standard error calculated
std_err = sem(weight_loss) 

# Mean value of weight loss
DW = sum(weight_loss)/len(weight_loss)

# Surface area calculation
l = 13E1 # cm
w = 7E1 # cm
h = 1E1 # cm
So = 2*l*w + 2*w*h + 2*l*h
print("Surface area: {So} cm^2".format(So = So))

# Starting concentration
Co = 16.825

# time
t = 3000

# Conversion to absolute weight
DWa = DW * So * 10**-6
print("Weight loss: {DWa} g".format(DWa = DW))

# Diffusion coefficient
Deff = (DWa**2*math.pi)/(4*So**2*Co**2*t)
print(Deff)