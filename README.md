[![PyPI version](https://badge.fury.io/py/PEUQSE.svg)](https://badge.fury.io/py/PEUQSE)
[![DOI:<10.11578/dc.20200908.7>](http://img.shields.io/badge/DOI-10.11578/dc.20200908.7-green.svg)](<https://doi.org/10.11578/dc.20200908.7>) [![Build Status](https://app.travis-ci.com/AdityaSavara/PEUQSE.svg?branch=master)](https://app.travis-ci.com/AdityaSavara/PEUQSE)  [![UnitTests:<UnitTests>](http://img.shields.io/badge/UnitTests-Passing-blue.svg)](<Passing-blue>) 
[![Downloads](https://pepy.tech/badge/peuqse/month)](https://pepy.tech/project/peuqse)
[![Open Source Love](https://badges.frapsoft.com/os/v3/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)



# Parameter Estimation and Uncertainty Quantification for Science and Engineering (PEUQSE)

Parameter estimation for complex physical problems often suffers from finding ‘solutions’ that are not physically realistic. The PEUQSE software (formerly named CheKiPEUQ) provides tools for finding physically realistic parameter estimates, graphs of the parameter estimate positions within parameter space, and plots of the final simulation results.

The MANUAL is located in the PEUQSE directory at https://github.com/AdityaSavara/PEUQSE/tree/master/PEUQSE
The recommended installation is to get Anaconda, then to install by conda or pip inside an anaconda prompt:

`pip install PEUQSE[COMPLETE]` (includes all optional dependencies) or `conda install conda-forge::peuqse`

For those who are not familiar with python, longer instructions are provided in the How to Get Started section below.

The primary documentation is in the Examples directory in a file named  ExamplesAndTutorialAndGettingStarted.rtf, but new users should read the 'How to Get Started' section, below.

# A Quick Introduction

Why should you use PEUQSE to get parameters from observed data? Because a few lines of code will give you realistic parameter estimates, possibly their uncertainties, and also some graphs.

If you don't see images in this quick introduction, then read online at: https://github.com/AdityaSavara/PEUQSE 

Consider a situation where we have three observed experimental data points with uncertainties:

![Image 001](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image001.png)

Their values, including uncertainties, are:<br>
160500 +/- 200000 <br>
810500 +/- 300000 <br>
1440500 +/- 200000 <br>

Consider it given that this situation is known to be described the following equation:<br>
<b>  y=(x-a)^2 + b  </b>

Assume we know that the physically realistic values of "a" and "b" are: <br>
a is expected to be 200 +/- 100   (this is the 1 sigma confidence interval) <br>
b is expected to be 500 +/- 200   (this is the 1 sigma confidence interval) <br>

If one tries to do a regular sum of squares fitting (conventional parameter estimation, CPE), we will not get realistic values for "a" and "b"  We get <b> a = 255, b = 139153 </b>.  The value for "a" is fine, but the value for "b" is not realistic.

<b> However, if we do a Bayesian Parameter Estimation (BPE) -- what PEUQSE is designed for -- </b> then we get the following answers: <b> a = 166 +/- 57, b= 509 +/- 198 </b>.  Where these errors are the 1 sigma credible intervals. Notice that now both of the parameters have physically realistic values.  We even have error bars that took into account the uncertainty! Additionally, the covariance matrix for the parameter estimates -- which includes the correlated uncertainties of estimated parameters -- is automatically exported to the directory (and is also accessible from python shells/scripts).

How good is the match in this example? <br>
The fitting (CPE) gives the red line below. This plot automatically generated by PEUQSE: 

![Image 003](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image003.png)

The Bayesian Parameter Estimation gives the black line below (and the red, not explained here). This plot automatically generated by PEUQSE: 

![Image 005](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image005.png)


We see that for this example, the CPE result from fitting and the BPE results do not look very different from each other. Both parameter estimation methods manage to stay in the error bars, yet the BPE result has a far more physically realistic pair of parameters!  This is the main purpose using PEUQSE in order to do BPE: it will tend to give more realistic parameter estimates, and can even give a type of uncertainty (called credible intervals) on the final estimates. More info can be found in the first PEUQSE paper: https://doi.org/10.1002/cctc.202000953

Here is the short and simple code that was required after making the model equation:

```python
import PEUQSE as PEUQSE
import PEUQSE.UserInput as UserInput 
UserInput.model['InputParameterPriorValues'] = [200, 500] #prior expected values for a and b
UserInput.model['InputParametersPriorValuesUncertainties'] = [100, 200] #these are 1 sigma uncertainties (1 standard deviation). The uncertainties are not correlated in this example, but a covariance matrix can be used instead.
UserInput.model['simulateByInputParametersOnlyFunction'] = simulation_model_00.simulation_function_wrapper #This just points to the User created model equation.
PE_object = PEUQSE.parameter_estimation(UserInput) #This creates a PEUQSE object of the parameter_estimation class from the UserInput.
PE_object.doMetropolisHastings() #This does the actual parameter_estimation by sampling.
PE_object.createAllPlots() #This creates all of the plots below!
```
PEUQSE has four types of posterior distribution sampling options as of Oct 2020. Two Markov Chain Monte Carlo sampling options: EnsembleSliceSampling and MetropolisHastings; and two even distribution (unbiased and unguided) sampling options: multistart grid sampling and multistart uniform distribution sampling. These are compared a little bit in Example00 in the examples. The latter two sampling types are rarely included in other packages but can be very important for scientific applications and for rough response surfaces. It is often important to know whether the solution found is a unique solution, and  the unbiased and unguided sampling can be used for verification of the uniqueness of a solution.

There is a logfile generated called mcmc_log_file.txt or permutations_log_file.txt (along with other files in the directory).
You will also get the following plots, some of which can be further customized, such as removing the bars from the contour plots. These plots automatically generated by PEUQSE: 


![Image 007](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image007.png)
![Image 009](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image009.png)
![Image 011](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image011.png)
![Image 013](https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/readmeImages/image013.png)


We can see that in this example the position and uncertainty in "a" narrowed more than that of "b".

# How To Get Started

1) First get anaconda (https://www.anaconda.com/products/individual). Then open an anaconda terminal (for windows, press the windows key and type "Anaconda Prompt"). Inside the anaconda terminal, type `pip install PEUQSE[COMPLETE]`

2) Then download the zipfile which has the Examples directory: https://github.com/AdityaSavara/PEUQSE/archive/master.zip

3) Open the examples directory, open the file ExamplesAndTutorialAndGettingStarted.rtf  (any program like Word can open this)

4) Start running the examples' run files (it is recommended to run them in the program called spyder, if you are not an advanced user. You can open it by typing "spyder" in the anaconda terminal).

5) Join the Google Group to ask for further help: https://groups.google.com/g/PEUQSE-users and add to existing issues discussions github https://github.com/AdityaSavara/PEUQSE/issues

6) To download new updates of PEUQSE, type 'pip install PEUQSE[COMPLETE] --upgrade'

* * *
# More info about the File Structure

The file structure is such that the file `./PEUQSE/__init__.py` is where the modules functions are loaded from. The main functionalities are inside `InverseProblem.py` , and various dependencies can be traced from those two files.

At present, the Examples directory contains various examples and a file named  ExamplesAndTutorialAndGettingStarted which provides some explanation.

In the relatively near future, the Examples may be separated into a complementary repository, and/or may be tucked away into the installation directory.

* * *
# Authors and Contributors

Authors: Aditya "Ashi" Savara, Eric A. Walker, Troy Gustke

The LICENSE is a BSD-3-Clause LICENSE.