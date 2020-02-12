import sys
sys.path.insert(0, '/mumpce/')
import mumpce.Project as mumpceProject
import mumpce.solution as mumpceSolution
import numpy as np
import matplotlib
from matplotlib import cm #EAW 2020/01/07
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput

class plotting_functions():
    def __init__(self, UserInput = UserInput, samples = False): # The plots require samples.  Other plot settings are probably plotting-package specific.
        self.UserInput = UserInput
        if not samples:
            print("Warning: Pass in the 'samples' keyword argument containing a numpy array of samples to plot.")
    
    def mu_and_cov_from_samples(self):
        mu = np.mean(self.samples, axis = 0)
        cov = np.cov(self.samples,rowvar=False)
        return mu, cov

    def mumpce_plots(self, model_parameter_info = {}, active_parameters = [], pairs_of_parameter_indices = [], posterior_mu_vector = [], posterior_cov_matrix = [], prior_mu_vector = [], prior_cov_matrix = []): # Pass empty keyword arguments for important parameters.  That way, warnings may be issued if they are not set.  There is not really a good default for these keyword arguments.  They depend entirely on the nature of the data being plotted.
        contour_settings_custom = {'figure_name','fontsize','num_y_ticks','num_x_ticks','colormap_posterior_customized','colormap_prior_customized','contours_normalized','colorbars'} # Here are the dictionary options, unassigned yet
        mumpceProjectObject = mumpceProject.Project() # A mumpce project object must be created.
        if len(model_parameter_info) == 0:
            print("Pass the 'model_parameter_info' argument to the mumpce_plots function.")
            model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': 'Parameter 0'},{'parameter_number': 1, 'parameter_name': 'Parameter 1'}])
        if len(active_parameters) == 0:
            print("Pass the 'active_parameters' argument to the mumpce_plots function.")
            active_parameters = np.array([0, 1]) 
        mumpceProjectObject.active_parameters = active_parameters
        if len(pairs_of_parameter_indices) == 0:
            print("Pass the 'pairs_of_parameter_indices' argument to the mumpce_plots function.")
            mumpceProjectObject.pairsOfParameterIndices = [[0, 1]]
        else:
            mumpceProjectObject.pairsOfParameterIndices = pairs_of_parameter_indices
        if len(posterior_mu_vector) == 0:
            print("Pass the 'posterior_mu_vector' argument to the mumpce_plots function.")
            posterior_mu_vector = np.array([-0.58888733,1.1200355])
        if len(posterior_cov_matrix) == 0:
            print("Pass the 'posterior_cov_matrix' argument to the mumpce_plots function.")
            posterior_cov_matrix = np.array([[ 0.0148872,-0.01894579],[-0.01894579,0.04284732]])
        if len(prior_mu_vector) == 0:
            print("Pass the 'prior_mu_vector' argument to the mumpce_plots functions.")
            prior_cov_matrix = np.array([-0.98888733,0.8200355])
        if len(prior_cov_matrix) == 0:
            prior_cov_matrix = 10*posterior_cov_matrix
            print("Pass the 'prior_cov_matrix' argument to the mumpce_plots functions.")
        mumpceProjectObject.model_parameter_info = model_parameter_info
        mumpceSolutionsObject = mumpceSolution.Solution(posterior_mu_vector, posterior_cov_matrix, initial_x=prior_mu_vector, initial_covariance=prior_cov_matrix)
        mumpceProjectObject.solution = mumpceSolutionsObject
        contour_settings_custom = {}
        if hasattr(UserInput,'figure_name'):
            contour_settings_custom['figure_name']=UserInput.figure_name
        else:
            contour_settings_custom['figure_name']='mumpce_plots_04'
        if hasattr(UserInput,'fontsize'):
            contour_settings_custom['fontsize'] = UserInput.fontsize        
        else:
            contour_settings_custom['fontsize'] = 'auto'
        if hasattr(UserInput,'num_y_ticks'):
            contour_settings_custom['num_y_ticks'] = UserInput.num_y_ticks
        else:
            contour_settings_custom['num_y_ticks'] = 'auto'
        if hasattr(UserInput,'num_x_ticks'):
            contour_settings_custom['num_x_ticks'] = UserInput.num_x_ticks
        else:
            contour_settings_custom['num_x_ticks'] = 'auto'
        if hasattr(UserInput,'colormap_posterior_customized'):
            contour_settings_custom['colormap_posterior_customized'] = UserInput.colormap_posterior_customized
        else:
            contour_settings_custom['colormap_posterior_customized'] = "Oranges"
        if hasattr(UserInput,'colormap_prior_customized'):
            contour_settings_custom['colormap_prior_customized'] = UserInput.colormap_prior_customized
        else:
            contour_settings_custom['colormap_prior_customized'] = "Greens"
        if hasattr(UserInput,'contours_normalized'):
            contour_settings_custom['contours_normalized'] = UserInput.contours_normalized
        else:
            contour_settings_custom['contours_normalized'] = False
        if hasattr(UserInput,'center_on'):
            contour_settings_custom['center_on'] = UserInput.center_on
        else:
            contour_settings_custom['center_on'] = 'prior'
        if hasattr(UserInput,'colorbars'):
            contour_settings_custom['colorbars'] = UserInput.colorbars
        else:
            contour_settings_custom['colorbars'] = True
        mumpceProjectObject.plot_pdfs(mumpceProjectObject.pairsOfParameterIndices, contour_settings_custom = contour_settings_custom)


    def seaborn_scatterplot_matrix(self):
        #fig0, ax0 = plt.subplots()
        #if UserInput.verbose:
        #    print(np.mean(rate_tot_array,axis = 0))
        #ax0.plot(np.array(experiments_df['AcH - T']),np.mean(rate_tot_array,axis = 0), 'r')
        #ax0.plot(np.array(experiments_df['AcH - T']),np.array(experiments_df['AcHBackgroundSubtracted'])/2000,'g')
        #ax0.set_ylim([0.00, 0.025])
        #ax0.set_xlabel('T (K)')
        #ax0.set_ylabel(r'$rate (s^{-1})$')
        #ax0.legend(['model posterior', 'experiments'])
        #fig0.tight_layout()
        #fig0.savefig('tprposterior.png', dpi=220)
        #posterior_df = pd.DataFrame(samples,columns=[UserInput.parameterNamesAndMathTypeExpressionsDict[x] for x in UserInput.parameterNamesList])
        #pd.plotting.scatter_matrix(posterior_df)
        #plt.savefig('scatter_matrix_posterior.png',dpi=220)
        return

    def rate_tot_plot(self):
        return