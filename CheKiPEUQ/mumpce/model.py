from abc import ABCMeta, abstractmethod

class Model(object):
    """This is the top-level class for a model object. Methods are defined as abstract methods, which must be 
    implemented as a subclass of model.
    
    The following methods are abstract methods within this class and must all be defined in order 
       
       * :func:`evaluate`: Returns the model value :math:`y`
       * :func:`sensitivity`: Returns the sensitivity coefficients :math:`S_{ij} = \\frac{d\ln y_i}{d\ln x_j}`
       * :func:`get_parameter`: Takes a parameter ID and returns the value of the corresponding model parameter
       * :func:`perturb_parameter`: Takes a parameter ID and replaces the corrsponding value with a new value
       * :func:`reset model`: Resets all model parameter values to their default values.
       * :func:`get_model_parameter_info`: Returns a list of dicts containing, at least, the parameter's name and possibly additional information
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def __str__(self):
        """Return some interesting information about the model
        """
        pass
    
    @abstractmethod
    def evaluate(self):
        """Runs the model once and return a single value
        
        :returns: model_value
        :rtype: float
        """
        pass
    @abstractmethod
    def sensitivity(self,perturbation,parameter_list,logfile):
        """Evaluates the sensitivity of the model value with respect to the model parameters
        
        :param perturbation: The amount to perturb each parameter during the sensitivity analysis
        :param parameter_list: The list of parameters to perturb. Usually this will be a list of parameter identifiers, which are usually ints or strs.
        :param logfile: The logging file that will contain the sensitivity calculation output.
        :type perturbation: float
        :type parameter_list: array_like
        :type logfile: str
        :returns: model_value,sensitivity_vector
        :rtype: tuple of float and 1d array
        """
        pass
    @abstractmethod
    def get_parameter(self,parameter_id):
        """Retrieves a model parameter's value
        
        :param parameter_id: The parameter identifier. This might be an int or a str or possibly another type, which will depend on the model.
        :returns: parameter_value
        :rtype: float
        """
        pass
    @abstractmethod
    def perturb_parameter(self,parameter_id,new_value):
        """Replaces a model parameter's value by a new value.
        
        :param parameter_id: The parameter identifier. This might be an int or a str or possibly another type, which will depend on the model.
        :param new_value: The amount to change the parameters value.
        :type new_value: float
        """
        pass
    @abstractmethod
    def reset_model(self):
        """Resets all model parameters to their original values"""
        pass
    @abstractmethod
    def get_model_parameter_info(self,number_parameters):
        """Gets information about the parameters, which will go up to the hosting measurement. This is called during instantiation of the model and normally would not be called at any other time.
       
        :returns: model_parameter_info
        :rtype: list
        """
        pass
    #@abstractmethod
    def prepare_for_save(self):
        pass