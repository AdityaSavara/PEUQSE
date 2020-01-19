# Type 'pytest' in the current directory to run tests.
# EAW 2020/01/17
import plotting_functions
from plotting_functions import plotting_functions
import UserInput_ODE_KIN_BAYES_SG_EW as UserInput
def test_mumpce_plots():
    plot_object = plotting_functions()
    assert plot_object.mumpce_plots(model_parameter_info = UserInput.model_parameter_info, active_parameters = UserInput.active_parameters, pairs_of_parameter_indices = UserInput.pairs_of_parameter_indices, posterior_mu_vector = UserInput.posterior_mu_vector, posterior_cov_matrix = UserInput.posterior_cov_matrix, prior_mu_vector = UserInput.prior_mu_vector, prior_cov_matrix = UserInput.prior_cov_matrix) == None
