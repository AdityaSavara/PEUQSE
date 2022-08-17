# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 07:37:49 2022

@author: NWT
This file graphs the value of Deff against the uncertainty provided to PEUQSE
It graphs MAP and MU_AP values, as well as standard deviations
"""
import matplotlib.pyplot as plt
import math

map_values = [-18.37662834,-18.36896044,-18.0852596,-17.95542376] # ,-17.95128962]
mu_ap_values = [-18.37667027,-18.36823056,-18.08384853,-17.9557197] # ,-14.26723996]
std_devs = [0.00093476,0.00948315,0.05376846,0.06246781] # ,5.86600188]
uncertainties = [-3.0,-2.0,-1.0,0.0] # ,1.0]
literature_values = []
for i in range(4): # 5):
    literature_values.append(math.log10(4.2*10**-19))

# Create a graph of MAP values against uncertainty
plt.figure(0)
plt.plot(uncertainties, map_values, 'ko-')
plt.plot(uncertainties, literature_values, 'g--')
plt.xlim([-3.1, 1.1])
plt.xlabel("Log base10 of uncertainty")
plt.ylabel("Log base10 of diffusion coefficient")
plt.legend(["PEUQSE values","Experimental value"])
plt.title("MAP values vs log base10 of uncertainty")
plt.show()

plt.figure(1)
plt.plot(uncertainties, mu_ap_values, 'ko-')
plt.plot(uncertainties, literature_values, 'g--')
plt.xlim([-3.1, 1.1])
plt.xlabel("Log base10 of uncertainty")
plt.ylabel("Log base10 of diffusion coefficient")
plt.legend(["PEUQSE values","Experimental value"])
plt.title("MU_AP values vs log base10 of uncertainty")
plt.show()

plt.figure(2)
plt.plot(uncertainties, std_devs, 'ko-')
plt.xlim([-3.1, 1.1])
plt.xlabel("Log base10 of uncertainty")
plt.ylabel("Standard deviation")
plt.legend(["PEUQSE values","Experimental value"])
plt.title("Standard deviation vs log base10 of uncertainty")
plt.show()