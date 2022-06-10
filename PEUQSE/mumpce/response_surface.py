import copy
import math
import numpy as np

class ResponseSurface(object):
    """A top level class describing a polynomial response surface.
    
    This class contains the methods for a response surface object. It largely parallels :func:`model` in that it has an :func:`evaluate` method which will return a value and a :func:`sensitivity` method which will return a value and a vector of sensitivities.
    
    The response surface is assumed here to be a second order polynomial :math:`y = z + a^{\\text{T}}x + x^{\\text{T}}bx`
    
    :param zero_term: The zero-order term of the response surface, :math:`z`
    :param a_terms: The first order terms of the response surface, :math:`a`
    :param b_terms: The second order terms of the response surface, :math:`b`
    :param c_terms: The third order terms of the response surface (not implemented)
    :param d_terms: The estimated error in the second order terms of the response surface
    :param active_parameters: If the list of active paremeters differs from measurement to measurement, it would be specified here.
    
    :type zero_term: float
    :type a_terms: ndarray(float), len(active_parameters)
    :type b_terms: ndarray(float), len(active_parameters)xlen(active_parameters)
    :type c_terms: ndarray(float)
    :type d_terms: ndarray(float), len(active_parameters)xlen(active_parameters)
    :type active_parameters: ndarray(int)
    
    
    :todo: implement active_parameters, implement c_terms
    
    
    """
    def __init__(self,
                 zero_term=None,
                 a_terms=None,b_terms=None,c_terms=None,
                 d_terms=None,
                 active_parameters=None):
        
        self.z = zero_term #: :math:`z`        
        self.a = a_terms #: :math:`a` 
        self.b = b_terms #: :math:`b` 
        
        self.c = c_terms
        self.d = d_terms
        self.active_parameters = active_parameters
        return
    
    def evaluate(self,x,cov_x=None):#x_full,exp_value,weight):
        """ Evaluates the response surface
        
        Computes the response value :math:`y = z + a^{\\text{T}}x + x^{\\text{T}}bx`
        If cov_x is supplied, returns the uncertainty :math:`\\sigma^2 = a^{\\text{T}}\\Sigma a + 2\\text{tr}((b\\Sigma)^2)``
        
        :param x: The parameter vector at which the response value is to be calculated
        :param cov_x: The covariance matrix among the parameters
        :type x: ndarray(float), len(active_parameters)
        :type cov_x: ndarray(float), len(active_parameters)xlen(active_parameters)
        :returns: response_value: :math:`y = z + a^{\\text{T}}x + x^{\\text{T}}bx` and  computed_unc: :math:`\\sigma^2 = a^{\\text{T}}\\Sigma a + 2\\text{tr}((b\\Sigma)^2)` if cov_x is supplied
        :rtype: tuple
        
        """
        #Select only the active parameters from the full parameter list
        #x = x[active_parameters]
        
        
        #Evaluate the response surface
        #Zero and first order terms
        response_value  = copy.deepcopy(self.z)
        response_value += np.dot(self.a,x)
        
        if cov_x is not None:
            a_times_cov = np.dot(self.a,cov_x)
            variance = np.dot(self.a,a_times_cov.T)
        
        #Second order terms (might not exist)
        if self.b is not None:
            b_times_x = np.dot(self.b,x)
            response_value += np.dot(b_times_x.T,x)
            if cov_x is not None:
                b_times_cov = np.dot(self.b,cov_x)
                variance += 2*np.trace(np.dot(b_times_cov,b_times_cov))
              
        #Third order terms not implemented
        
        if cov_x is not None:
            computed_unc = math.sqrt(variance)
            return response_value,computed_unc
        
        
        return response_value
    
    def sensitivity(self,x):
        """ Evaluates the response surface and response surface gradient
        
        Computes the response value :math:`y = z + a^{\\text{T}}x + x^{\\text{T}}bx` and the response surface gradient :math:`\\frac{dy}{dx} = a + 2bx`
        
        :param x: The parameter vector at which the response value is to be calculated
        :type x: ndarray(float), len(active_parameters)
        :returns: response_value: :math:`y = z + a^{\\text{T}}x + x^{\\text{T}}bx` and response_gradient: :math:`\\frac{dy}{dx} = a + 2bx`
        :rtype: tuple
        
        """
        #Select only the active parameters from the full parameter list
        #x = x[active_parameters]
        
        #Evaluate the response surface
        #Zero and first order terms
        response_value  = copy.deepcopy(self.z)
        response_value += np.dot(self.a,x)
        
        response_grad = copy.deepcopy(self.a)
        
        #Second order terms (might not exist)
        if not(self.b is None):
            b_times_x = np.dot(self.b,x)
            response_value += np.dot(b_times_x.T,x)
            response_grad += 2*b_times_x
              
        #Third order terms not implemented

        return response_value,response_grad