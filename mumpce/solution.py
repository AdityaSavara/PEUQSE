import numpy as np

class Solution(object):
    """A top level class for a constrained model with uncertainty
    
    :param solution_x: The solution vector
    :param covariance_x: The covariance matrix among the elements of the solution vector
    :param second_order_x: A structure describing the second order variation in the elements of the solution vector
    :type solution_x: ndarray,float
    :type covariance_x: ndarray,float
    :type second_order_x:
    
    """
    def __init__(self,
                solution_x,covariance_x=None,second_order_x=None,initial_x=None,initial_covariance=None):
        self.x     = solution_x #: The solution vector, :math:`x_{opt}`
        self.cov   = covariance_x #: The covariance matrix, :math:`\Sigma`
        self.alpha = np.linalg.cholesky(covariance_x) #: The lower triangular decomposition :math:`\alpha` where :math:`\Sigma = \alpha \alpha^{\text{T}}` 
        self.beta  = second_order_x
        
        self.x_i   = initial_x
        #"""The initial guess vector :math:`x_{init}` used to start the optimization, if not zero"""
        self.cov_i = initial_covariance
        #"""The initial covariance matrix :math:`\Sigma_{init}` used to start the optimization, if not zero"""
        self.alpha_i = None
        if initial_covariance is not None:
            self.alpha_i = np.linalg.cholesky(initial_covariance)
        
        
        return
    
    def update(self,new_x=None,new_cov=None):
        """Updates the solution and covariance in the Solution
        """
        if new_x is not None:
            self.x = new_x
        if new_cov is not None:
            self.cov = new_cov
            self.alpha = np.linalg.cholesky(new_cov)
        return