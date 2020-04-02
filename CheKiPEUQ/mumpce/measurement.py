import numpy as np
import pickle
from response_surface import ResponseSurface

def idfunc(*arg,**kwargs):
    if len(arg) == 1:
        return arg[0]
    return arg

#print('loading')

try:
    import tqdm
    tqfunc = tqdm.tqdm_notebook
    #print('tqdm found')
except ImportError:
    #is not available
    tqfunc = idfunc

class Measurement(object):
    """A top level class for a measurement object
    
    This class is intended to contain metadata for the measurement it defines as well as a simulation model and a response 
    surface for optimization.
    
    :param name: The name of the measurement
    :param comment: An arbitrary text string describing the measurement
    :param model: The simulation model
    :param value: The measured value from the experiment
    :param uncertainty: The uncertainty in the measured value
    :param active_parameters: The list of active parameters for this simulation. Not normally defined at object creation
    :param parameter_uncertainties: The list of uncertainties in the active parameters. Not normally defined at object creation
    :param response: The response surface for this simulation. Not normally defined at object creation
    :param response_type: Whether the response is to be linear ('linear', :math:`y = a^{\\text{T}}x + z`) or logarithmic ('log', :math:`\ln y = a^{\\text{T}}x + z`)
    :param response_perturbation: The perturbation in the parameters for response surface generation
    :param model_value: The computed value from the simulation model for this measurement
    :param sensitivity_list: The list of sensitivities of the model value to each parameter
    :type name: str
    :type comment: str
    :type model: model object
    :type value: float
    :type uncertainty: float
    :type active_parameters: int list
    :type parameter_uncertainties: float list
    :type response: response surface object
    :type response_type: str
    :type response_perturbation: float
    :type model_value: float
    :type sensitivity_list: float list
    
    """
    def __init__(self,
                 name=None,
                 model=None,
                 value=None,
                 uncertainty=None,
                 active_parameters=None,
                 parameter_uncertainties=None,
                 response=None,
                 response_type='linear',
                 response_perturbation=0.3,
                 model_value=None,
                 sensitivity_list=None,
                 comment=None
                ):
        
        self.tqfunc = tqfunc
        
        #Define the name
        self.name = name
        self.comment = comment
        
        #Define the computer model that will be used to simulate this measurement
        self.model = model
        
        
        #Define the experimental value
        self.value = value
        self.uncertainty = uncertainty
        
        #Define the status of the experiment (used by Project for reporting purposes)
        self._status = 'Active'
        
        #Define the response surface and active parameters for this measurement
        self.active_parameters = active_parameters
        self.response = response
        self.parameter_uncertainties = parameter_uncertainties
        
        #Define the perturbations for sensitivity analysis
        self.response_perturbation = response_perturbation
        self.response_sensitivity = 1.0e-3
        self.response_type = response_type
        
        #Take the sensitivity list from the initialization
        self.model_value = model_value
        self.model_uncertainty = None
        self.sensitivity_list=sensitivity_list
        
        #Create attributes for optimized values and uncertainties
        self.optimized_value = None
        self.optimized_uncertainty = None
        self.consistency = None
        self.weighted_consistency = None
        self.entropy = None
        
        return
    
    def __str__(self):
        """Returns Name (Status): str(self.model)
        """
        modelstr = str(self.model)
        
        meas_str = "{} ({}): {}".format(self.name,self._status,modelstr,)
        
        if self.comment:
            meas_str = meas_str + "\n     {}".format(self.comment)
        
        return meas_str
        
    def set_value(self,value,uncertainty):
        """Sets the experimental value and uncertainty for this measurement
        
        :param value: The experiment's measured value
        :param uncertainty: The measurement uncertainty in the value
        :type value: float
        :type uncertainty: float
        
        """
        self.value = value
        self.uncertainty = uncertainty
        return
    
    def make_response(self): #(self,zero_term,perterbations,sensitivities):
        """Generates a sensitivity_analysis_based response surface for this measurement
        """
        #zero_term = self.evaluate
        self.model.reset_model()
        logfile_name = self.name + '_resp_log.out'
        response_logfile = open(logfile_name,'w')
        
        sensitivity_args = (self.response_sensitivity,
                           self.active_parameters,
                           response_logfile)
        sensitivity_kw = dict(tq=True)
        
        zero_term, sens_zero = self.model.sensitivity(*sensitivity_args)
        
        if self.response_type == 'log':
            zero_term = np.log(zero_term)
        
        #print zero_term
        #print sens_zero
        
        number_params = len(self.active_parameters)
        
        #Calculate the multipliers that will be used for the SAB sensitivity calculations
        multipliers = self.parameter_uncertainties ** self.response_perturbation
        #print multipliers
        
        perturbations = np.zeros((number_params,2))
        sens_positive = np.zeros((number_params,number_params))
        sens_negative = np.zeros_like(sens_positive)
        
        #Changed this so that tqdm will be used if it is available, but otherwise not
        #for (parameter_number,parameter) in tqdm.tqdm(enumerate(self.active_parameters)):
        for (parameter_number,parameter) in enumerate(self.tqfunc(self.active_parameters)):
            self.model.reset_model()
            base_value = self.model.get_parameter(parameter)
            param_name = self.model.model_parameter_info[parameter]['parameter_name']
            
            #print 'Parameter = ', parameter
            
            
            positive_perturbation = multipliers[parameter_number]
            negative_perturbation = 1/positive_perturbation
            
            #print positive_perturbation
            #print negative_perturbation
            
            #Positive perturbation
            response_logfile.write('\nParameter number = {: 4d} {:30s}\n'.format(parameter,param_name))
            response_logfile.write('Positive perturbation = {: 10.5e}\n'.format(positive_perturbation))
            self.model.perturb_parameter(parameter,positive_perturbation*base_value)
            value_pos, sens_pos = self.model.sensitivity(*sensitivity_args,**sensitivity_kw)
            
            #Negative perturbation
            response_logfile.write('Negative perturbation = {: 10.5e}\n'.format(negative_perturbation))
            self.model.perturb_parameter(parameter,negative_perturbation*base_value)
            value_neg, sens_neg = self.model.sensitivity(*sensitivity_args,**sensitivity_kw)
            
            perturbations[parameter_number,:] = [value_pos, value_neg]
                
            sens_positive[parameter_number,:] = sens_pos
            sens_negative[parameter_number,:] = sens_neg
            
            #print "Positive sensitivity:", sens_pos
            #print "Negative sensitivity:", sens_neg
            #print ''
        
        self.model.reset_model()
        #First order terms of response surface
        if self.response_type == 'log':
            perturbations = np.log(perturbations)
        a_terms = (perturbations[:,0] - perturbations[:,1]) / (2 * self.response_perturbation)
        #Second order terms of response surface
        b_terms_first = (sens_positive - sens_negative) * np.log(self.parameter_uncertainties) / (4 * self.response_perturbation)
        
        if self.response_type == 'linear':
            #a_terms = a_terms * zero_term
            b_terms_first = b_terms_first * zero_term
        
        b_terms = (b_terms_first + b_terms_first.T)/2# - np.diag(np.diag(b_terms_first))
        d_terms = (b_terms_first - b_terms_first.T)/2
        
        
        self.response = ResponseSurface(zero_term=zero_term,
                                         a_terms=a_terms,
                                         b_terms=b_terms,
                                         d_terms=d_terms,
                                         active_parameters=self.active_parameters)
        return
    
    def evaluate(self):
        """Evaluates the model once and sets self.model_value to the returned value.
        
        :returns: model_value
        """
        self.model_value = self.model.evaluate()
        if self.response_type == 'log':
            self.model_value = np.log(self.model_value)
        return self.model_value
    
    def evaluate_sensitivity(self,perturbation=0.05):
        """Conducts a sensitivity analysis on the model and storee the nominal value in self.model_value and the sensitivity in self.sensitivity_list
        
        :param perturbation: The amount to perturb each parameter when conducting the sensitivity analysis
        :type perturbation: float
        """
        all_parameters = np.arange(self.model.number_parameters,dtype=int)
        logfile_name = self.name + '_sen.out'
        with open(logfile_name,'w') as logfile:
                self.model_value,self.sensitivity_list = self.model.sensitivity(perturbation=perturbation,
                                                                                parameter_list=all_parameters,
                                                                                logfile=logfile
                                                                               )
        if self.response_type == 'log':
            self.model_value = np.log(self.model_value)
        return
    
    def print_sorted_sensitivity(self,sensitivity=None,max_number=None):
        """Sorts the parameters by sensitivity coefficient and prints them.
        
        :param sensitivity: The list of sensitivity coefficients. If None, will use self.sensitivity_list
        :param max_number: If not None, will only print the top max_number parameters
        :type sensitivity: array_like or None
        :type max_number: int or None
        """
        
        if sensitivity is None:
            sensitivity = self.sensitivity_list
        
        sorted_param_nums = np.argsort(abs(sensitivity))
        
        if max_number is None:
            max_number = self.model.number_parameters
        
        print_params = sorted_param_nums[-1:-1*max_number:-1]
        
        for print_param in print_params:
            param_info = self.model.model_parameter_info[print_param]
            param_name = param_info['parameter_name']
            print('{: 4d} {: 10.4e}  {}'.format(print_param,
                                                sensitivity[print_param],
                                                param_name
                                               )
                )
            
           
        return
    
    def evaluate_response(self,x):
        """Evaluates the response surface for this measurement.
        
        :param x: The set of parameter values at which the surface must be evaluated
        :type x: 1d array
        :returns: response
        :rtype: float
        """
        response = self.response.evaluate(x)
        return response
    
    def sensitivity_response(self,x):
        """Evaluates the response surface and the response surface gradient for this measurement.
        
        :param x: The set of parameter values at which the surface must be evaluated
        :type x: 1d array
        :returns: response,response_gradient
        :rtype: tuple
        """
        response,response_grad = self.response.sensitivity(x)
        return response,response_grad
    
    def evaluate_uncertainty(self,x,cov):
        """Evaluates the response surface for this measurement and compute its uncertainty
        
        :param x: The set of parameter values at which the surface must be evaluated
        :param cov: The covariance matrix among the model parameters
        :type x: 1d array
        :type cov: 2d array
        :returns: response,response_uncertainty
        :rtype: tuple
        """
        response,response_uncertainty = self.response.evaluate(x,cov)
        return response,response_uncertainty
    
    def save(self):
        """Saves a pickled representation of the measurement
        """
        filename = self.name + '.save'
        
        with open(filename,'wb') as f:
            pickle.dump(self,f)
        return
    
    def prepare_for_save(self):
        self.model.prepare_for_save()
        return
    
    def _save(self):
        """Saves the model value, sensitivity list, and response surface to disk. By default, they are saved to the file 
        'self.name'.npz
        """
        
        #The response surface may not exist, so create some dummy variables to ensure that the measurement can be saved
        localz = None
        locala = None
        localb = None
        locald = None
        
        print('Checking for response surface ...')
        if self.response is not None:
            print('Response surface exists, saving response surface ...')
            localz = self.response.z
            locala = self.response.a
            localb = self.response.b
            locald = self.response.d
        
        outputfilename = self.name + '.npz'
        print ('Saving to output file: {}'.format(outputfilename))
        with open(outputfilename, 'w+') as outputfile:
            outputfile.seek(0)
            np.savez(outputfile,
                     model_value = self.model_value,
                     sensitivity_list = self.sensitivity_list,
                     z = localz,
                     a = locala,
                     b = localb,
                     d = locald,
                     active_parameters = self.active_parameters,
                     active_parameter_uncertainties = self.parameter_uncertainties,
                     value=self.value,
                     uncertainty=self.uncertainty
                    )
        return
    
    def load(self):
        """Loads the model value, sensitivity list, and response surface from disk. By default, they are loaded from the 
        file 'self.name'.npz
        """
        inputfilename = self.name + '.npz'
        print ('Loading from output file: {}'.format(inputfilename))
        with open(inputfilename, 'r') as inputfile:
            response_data = np.load(inputfile)
            
            #Check to see if the file has response surface data in it
            localz = response_data['z']
            print('Checking response surface...')
            #Do not load the response surface if there is no response surface data
            if localz is not None:
                print('Response surface exists, loading response surface ...')
                self.response = ResponseSurface(zero_term=response_data['z'],
                                                 a_terms=response_data['a'],
                                                 b_terms=response_data['b'],
                                                 d_terms=response_data['d'],
                                                 active_parameters=response_data['active_parameters']
                                                )
            self.model_value = response_data['model_value']
            if self.model_value: print('Model value loaded')
                
            self.sensitivity_list = response_data['sensitivity_list']
            if self.sensitivity_list.any(): print('Sensitivity vector loaded')
                
            self.value = response_data['value']
            if self.value: print('Experimental value loaded')
                
            self.uncertainty = response_data['uncertainty']
            if self.uncertainty: print('Experimental uncertainty loaded')
                
            self.active_parameters = response_data['active_parameters']
            if self.active_parameters.any(): print('Active parameters loaded')
                
            self.parameter_uncertainties = response_data['active_parameter_uncertainties']
            if self.parameter_uncertainties.any(): print('Parameter uncertainties loaded')
                
        return
    
    def print_model_values(self):
        parameter_info = self.model.model_parameter_info
        
        headname = 'Parameter name'
        headv = 'Value'
        headf = 'Uncert'
        headx = 'FactVal'
        heads = 'FactUnc'
        headfs = 'NewVal'
        headsig = 'Uncert'
        
        head_args = (headname,headv,headf)
        
        head_format = '{:40s}   {:8s} {:7s}'
        
        #print(head_format.format(*head_args))
        
        header = head_format.format(*head_args)
        
        carriage_return = '\n'
        output = ''
        
        output = carriage_return.join((output,header))
        
        for active_num,param in enumerate(self.active_parameters):
            
            param_name = parameter_info[param]['parameter_name']
            
            value = parameter_info[param]['parameter_value']#self.model.get_parameter(param)
            this_unc = self.parameter_uncertainties[active_num]
            
            print_args = (param_name[:40],value,this_unc)
            
            print_string = '{:40s} : {: 8.2g} {: 7.2f}'
            #if abs(value) > 5000:
            #    print_string = '{:30s} : {: 7.2e} {: 7.2f}'
            
            #print(print_string.format(*print_args))
            line = print_string.format(*print_args)
            output = carriage_return.join((output,line))
        
        return output
    
    def get_active_names(self):
        names = []
        for param in self.active_parameters:
            param_name = self.model.model_parameter_info[param]['parameter_name']
            names += [param_name]
        return names
    
    def interpret_model(self,x,cov):
        parameter_info = self.model.model_parameter_info
        
        headname = 'Parameter name'
        headv = 'Value'
        headf = 'Uncert'
        headx = 'FactVal'
        heads = 'FactUnc'
        headfs = 'NewVal'
        headsig = 'Uncert'
        
        head_args = (headname,headv,headf,headx,heads,headfs,headsig)
        
        head_format = '{:40s}   {:8s} {:7s} {:7s} {:7s} {:8s} {:7s}'
        
        #print(head_format.format(*head_args))
        
        carriage_return = '\n'
        output = ''
        
        header = head_format.format(*head_args)
        output = carriage_return.join((output,header))
        
        for active_num,param in enumerate(self.active_parameters):
            
            param_name = parameter_info[param]['parameter_name']
            
            value = parameter_info[param]['parameter_value']#self.model.get_parameter(param)
            
            this_x = x[active_num]
            this_std = 2*np.sqrt(cov[active_num,active_num])
            this_unc = self.parameter_uncertainties[active_num]
            
            multiplier = this_unc ** this_x
            new_value = value*multiplier
            new_uncertainty = this_unc ** (this_std)
            
            print_args = (param_name[:40],value,this_unc,this_x,this_std,new_value,new_uncertainty)
            
            print_string = '{:40s} : {: 8.2g} {: 7.2f} {: 7.2f} {: 7.2f} {: 8.2g} {: 7.2f}'
            #if abs(value) > 5000:
            #    print_string = '{:30s} : {: 7.2e} {: 7.2f} {: 7.2f} {: 7.2f} {: 7.2e} {: 7.2f}'
            line = print_string.format(*print_args)
            output = carriage_return.join((output,line))
            #print(print_string.format(*print_args))
        
        return output
    
    def modify_model(self,x):
        for active_num,param in enumerate(self.active_parameters):
            #value = self.model.get_parameter(param)
            this_x = x[active_num]
            this_unc = self.parameter_uncertainties[active_num]
            
            multiplier = this_unc ** this_x
            new_value = multiplier#*value
            
            self.model.perturb_parameter(param,new_value)
        return
    
    def get_model_values(self):
        values = []
        uncertainties = []
        parameter_info = self.model.model_parameter_info
        for active_num,param in enumerate(self.active_parameters):
            value = parameter_info[param]['parameter_value']#self.model.get_parameter(param)
            this_unc = self.parameter_uncertainties[active_num]
            values += [value]
            uncertainties += [this_unc]
        values = np.array(values)
        uncertainties = np.array(uncertainties)
        return values,uncertainties
    
    def get_opt_values(self,x,cov):
        new_values = []
        new_uncertainties = []
        parameter_info = self.model.model_parameter_info
        for active_num,param in enumerate(self.active_parameters):
            value = parameter_info[param]['parameter_value']#self.model.get_parameter(param)
            this_x = x[active_num]
            this_std = 2*np.sqrt(cov[active_num,active_num])
            this_unc = self.parameter_uncertainties[active_num]

            multiplier = this_unc ** this_x
            new_value = value*multiplier
            new_uncertainty = this_unc ** (this_std)

            new_values += [new_value]
            new_uncertainties += [new_uncertainty]
        new_values = np.array(new_values)
        new_uncertainties = np.array(new_uncertainties)

        return new_values,new_uncertainties
            