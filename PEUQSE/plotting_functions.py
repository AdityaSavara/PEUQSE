import sys
#sys.path.insert(0, '/mumpce/')
import PEUQSE.mumpce.Project as mumpceProject
import PEUQSE.mumpce.solution as mumpceSolution
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm #EAW 2020/01/07
import platform
if platform.system() == 'posix':
    matplotlib.use('Agg') #added by A. Savara June 29th, 2021.
import copy

class plotting_functions_class():
    def __init__(self, UserInput, samples = False): # The plots require samples.  Other plot settings are probably plotting-package specific.
        self.UserInput = UserInput
#        if not samples: #Deprecated.
#            print("Warning: Pass in the 'samples' keyword argument containing a numpy array of samples to plot.")
    
    def mu_and_cov_from_samples(self):
        mu = np.mean(self.samples, axis = 0)
        cov = np.cov(self.samples,rowvar=False)
        return mu, cov

    def mumpce_plots(self, model_parameter_info = {}, active_parameters = [], pairs_of_parameter_indices = [], posterior_mu_vector = 0, posterior_cov_matrix = 0, prior_mu_vector = 0, prior_cov_matrix = 0, contour_settings_custom = {'figure_name','fontsize','max_num_y_ticks','max_num_x_ticks','colormap_posterior_customized','colormap_prior_customized','contours_normalized','colorbars','axis_limits'}): # Pass empty keyword arguments for important parameters.  That way, warnings may be issued if they are not set.  There is not really a good default for these keyword arguments.  They depend entirely on the nature of the data being plotted.
        mumpceProjectObject = mumpceProject.Project() # A mumpce project object must be created.
        if len(model_parameter_info) == 0:
            print("Pass the 'model_parameter_info' argument to the mumpce_plots function.")
            model_parameter_info = np.array([{'parameter_number': 0, 'parameter_name': 'Parameter 0'},{'parameter_number': 1, 'parameter_name': 'Parameter 1'}])
        if len(active_parameters) == 0:
            print("Pass the 'active_parameters' argument to the mumpce_plots function.")
            active_parameters = np.array([0, 1]) 
        mumpceProjectObject.active_parameters = active_parameters
        #if len(pairs_of_parameter_indices) == 0:
        #    print("Pass the 'pairs_of_parameter_indices' argument to the mumpce_plots function.")
        #    mumpceProjectObject.pairsOfParameterIndices = [[0, 1]]
        #else:
        mumpceProjectObject.pairsOfParameterIndices = pairs_of_parameter_indices
        #if not posterior_mu_vector == 0:
        #    print("Pass the 'posterior_mu_vector' argument to the mumpce_plots function.")
        #    posterior_mu_vector = np.array([-0.58888733,1.1200355])
        #if not posterior_cov_matrix == 0:
        #    print("Pass the 'posterior_cov_matrix' argument to the mumpce_plots function.")
        #    posterior_cov_matrix = np.array([[ 0.0148872,-0.01894579],[-0.01894579,0.04284732]])
        #if not prior_mu_vector == 0:
        #    print("Pass the 'prior_mu_vector' argument to the mumpce_plots functions.")
        #    prior_cov_matrix = np.array([-0.98888733,0.8200355])
        #if not prior_cov_matrix == 0:
        #    prior_cov_matrix = 10*posterior_cov_matrix
        #    print("Pass the 'prior_cov_matrix' argument to the mumpce_plots functions.")
        mumpceProjectObject.model_parameter_info = model_parameter_info
        mumpceSolutionsObject = mumpceSolution.Solution(posterior_mu_vector, posterior_cov_matrix, initial_x=prior_mu_vector, initial_covariance=prior_cov_matrix)
        mumpceProjectObject.solution = mumpceSolutionsObject
        #if hasattr(UserInput,'figure_name'):
        #    contour_settings_custom['figure_name']=UserInput.figure_name
        #else:
        #    contour_settings_custom['figure_name']='mumpce_plots_04'
        #if hasattr(UserInput,'fontsize'):
        #    contour_settings_custom['fontsize'] = UserInput.fontsize        
        #else:
        #    contour_settings_custom['fontsize'] = 'auto'
        #if hasattr(UserInput,'max_num_y_ticks'):
        #    contour_settings_custom['max_num_y_ticks'] = UserInput.max_num_y_ticks
        #else:
        #    contour_settings_custom['max_num_y_ticks'] = 'auto'
        #if hasattr(UserInput,'max_num_x_ticks'):
        #    contour_settings_custom['max_num_x_ticks'] = UserInput.max_num_x_ticks
        #else:
        #    contour_settings_custom['max_num_x_ticks'] = 'auto'
        #if hasattr(UserInput,'colormap_posterior_customized'):
        #    contour_settings_custom['colormap_posterior_customized'] = UserInput.colormap_posterior_customized
        #else:
        #    contour_settings_custom['colormap_posterior_customized'] = "Oranges"
        #if hasattr(UserInput,'colormap_prior_customized'):
        #    contour_settings_custom['colormap_prior_customized'] = UserInput.colormap_prior_customized
        #else:
        #    contour_settings_custom['colormap_prior_customized'] = "Greens"
        #if hasattr(UserInput,'contours_normalized'):
        #    contour_settings_custom['contours_normalized'] = UserInput.contours_normalized
        #else:
        #    contour_settings_custom['contours_normalized'] = False
        #if hasattr(UserInput,'center_on'):
        #    contour_settings_custom['center_on'] = UserInput.center_on
        #else:
        #    contour_settings_custom['center_on'] = 'prior'
        #if hasattr(UserInput,'colorbars'):
        #    contour_settings_custom['colorbars'] = UserInput.colorbars
        #else:
        #    contour_settings_custom['colorbars'] = True
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

def sampledParameterHistogramMaker(parameterSamples, parameterName, parameterNamesAndMathTypeExpressionsDict, sampledParameterFiguresDictionary, sampledParameterAxesDictionary, directory='', parameterInitialValue=None, parameterMAPValue=None, parameterMuAPValue=None, histogram_plot_settings={}):

        if len(histogram_plot_settings) == 0: #this means histogram_plot_settings was not provided, or was blank, in which case we will populate variables with defaults.
            histogram_plot_settings = copy.deepcopy(histogram_plot_settings) #first make a fresh copy so the original dictionary is not changed.
            histogram_plot_settings['histograms_as_density'] = False
            histogram_plot_settings['y_label'] = ''
            histogram_plot_settings['show_initial_value'] = True
            histogram_plot_settings['show_MAP_value'] = True
            histogram_plot_settings['show_muAP_value'] = True
            histogram_plot_settings['x_label_size'] = 16
            histogram_plot_settings['y_label_size'] = 16
            histogram_plot_settings['axis_font_size'] = 16
            histogram_plot_settings['dpi'] = 220
            histogram_plot_settings['vertical_linewidth'] = 1.5


        #The axis font size argument needs to be parsed into another form. #Code was made following answer by "binaryfunt" here https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
        axis_font = {'size':str(histogram_plot_settings['axis_font_size'])}
        
        #Now to get some more strings needed and then make the plot object.
        parameterIndex = list(parameterNamesAndMathTypeExpressionsDict).index(parameterName)
        sampledParameterFiguresDictionary[parameterName], sampledParameterAxesDictionary[parameterName] = plt.subplots()   #making plt objects    
        sampledParameterAxesDictionary[parameterName].hist(parameterSamples[:,parameterIndex], density=histogram_plot_settings['histograms_as_density']) #filling the object with data
        #setting the labels etc. and then exporting.

        if histogram_plot_settings['y_label'] == '': #will use defaults if ''        
            if histogram_plot_settings['histograms_as_density'] == False:
                sampledParameterAxesDictionary[parameterName].set_ylabel('Frequency', **axis_font)
            if histogram_plot_settings['histograms_as_density'] == True:
                sampledParameterAxesDictionary[parameterName].set_ylabel('Probability Density', **axis_font)
        sampledParameterAxesDictionary[parameterName].set_xlabel(parameterNamesAndMathTypeExpressionsDict[parameterName], **axis_font)

        # Add vertical lines at initial value, MAP, and, muAP if determined by UserInput
        # Make sure that the value passed has a value
        if histogram_plot_settings['show_initial_value'] and type(parameterInitialValue) != type(None):
            sampledParameterAxesDictionary[parameterName].axvline(x=parameterInitialValue, color='#00A5DF', linestyle='--', linewidth=histogram_plot_settings['vertical_linewidth'])
        if histogram_plot_settings['show_MAP_value'] and type(parameterMAPValue) != type(None):
            sampledParameterAxesDictionary[parameterName].axvline(x=parameterMAPValue, color='r', linestyle='--', linewidth=histogram_plot_settings['vertical_linewidth'])
        if histogram_plot_settings['show_muAP_value'] and type(parameterMuAPValue) != type(None):
            sampledParameterAxesDictionary[parameterName].axvline(x=parameterMuAPValue, color='k', linestyle='--', linewidth=histogram_plot_settings['vertical_linewidth'])


        sampledParameterAxesDictionary[parameterName].tick_params(axis='x', labelsize=histogram_plot_settings['x_label_size']) #TODO: make these labels sizes a setting that can be changed.
        sampledParameterAxesDictionary[parameterName].tick_params(axis='y', labelsize=histogram_plot_settings['y_label_size'])
        sampledParameterFiguresDictionary[parameterName].tight_layout()
        sampledParameterFiguresDictionary[parameterName].savefig(directory+'Histogram_sampling_'+str(parameterIndex)+'_'+parameterName+'.png', dpi=histogram_plot_settings['dpi'])

        
        
        #The above block makes code kind of like this in a dynamic fashion. Since we know how many we will need, a dictionary is used to avoid the need for 'exec' statements when making new parameters.
        # fig2, ax2 = plt.subplots()
        # ax2.hist(samples[:,1])
        # ax2.set_ylabel('frequency')
        # ax2.set_xlabel(r'$E_{a2}$')
        # fig2.tight_layout()
        # fig2.savefig('Ea2.png', dpi=220)

    #Make histograms for each parameter. Need to make some dictionaries where relevant objects will be stored.
def makeHistogramsForEachParameter(parameterSamples,parameterNamesAndMathTypeExpressionsDict, directory='', parameterInitialValue=None, parameterMAPValue=None, parameterMuAPValue=None, histogram_plot_settings={}):
    sampledParameterFiguresDictionary = copy.deepcopy(parameterNamesAndMathTypeExpressionsDict) #This must be a deep copy to perserve original.
    sampledParameterAxesDictionary = copy.deepcopy(parameterNamesAndMathTypeExpressionsDict) #This must be a deep copy to preserve original.
    #The below code was originally added by Eric Walker, then modified by Troy Gustke to add in the initialValue, MAPvalue, and mu_apvALUE in June 2021, and merged in by A. Savara.    
    for key, iv, mp, mup in zip(parameterNamesAndMathTypeExpressionsDict,parameterInitialValue,parameterMAPValue,parameterMuAPValue):
        parameterName = key
        initialValue = iv
        MAPValue = mp
        Mu_APValue = mup
        sampledParameterHistogramMaker(parameterSamples, parameterName, parameterNamesAndMathTypeExpressionsDict, sampledParameterFiguresDictionary, sampledParameterAxesDictionary, directory=directory, parameterInitialValue=initialValue, parameterMAPValue=MAPValue, parameterMuAPValue=Mu_APValue, histogram_plot_settings=histogram_plot_settings)        

def createSimulatedResponsesPlot(x_values, listOfYArrays, plot_settings={}, listOfYUncertaintiesArrays=[], showFigure=True, directory=''):
    exportFigure = True #This variable should be moved to an argument or something in plot_settings.
    #First put some defaults in if not already defined.
    x_values = np.array(x_values)
    if len(np.shape(x_values)) == 0: #This 2d line breaks regular arrays, but works when a 'zero length' array comes in (I don't understand how a zero length array can occur, but it has happened.)
        x_values = np.atleast_2d(x_values)
    if 'x_label' not in plot_settings: plot_settings['x_label'] = ''
    if 'y_label' not in plot_settings: plot_settings['y_label'] = ''
    if 'legendLabels' not in plot_settings: plot_settings['legendLabels'] = ''
    if 'figure_name' not in plot_settings: plot_settings['figure_name'] = 'simulatedResponsesPlot'
    if 'dpi' not in plot_settings: plot_settings['dpi']=220
    fig0, ax0 = plt.subplots()
    if 'fontdict' in plot_settings: 
        #There are various things that could be added to this fontdict. #https://www.tutorialexample.com/understand-matplotlib-fontdict-a-beginner-guide-matplotlib-tutorial/
        fontdict = plot_settings['fontdict']
        if 'size' in fontdict:
            ax0.tick_params(axis='x', labelsize=fontdict['size'])
            ax0.tick_params(axis='y', labelsize=fontdict['size'])
    else:
        fontdict = None #initializing with the matplotlib default
    ax0.set_xlabel(plot_settings['x_label'], fontdict=fontdict)
    ax0.set_ylabel(plot_settings['y_label'], fontdict=fontdict) #TODO: THis is not yet generalized (will be a function)
    #The error linewidth is going to be set to a thick value if we have a small number of points.
    if 'error_linewidth' in plot_settings: error_linewidth = plot_settings['error_linewidth'] # making a convenient local variable.
    if str(error_linewidth).lower() == 'auto':
        if len(x_values) == 1:
            error_linewidth = 20
        elif len(x_values) > 10:
            error_linewidth = 1
        elif len(x_values) <= 10:
            error_linewidth = 10
    if str(error_linewidth).lower() == 'none': error_linewidth = 0 #This will hide the rror bars if they are not desired.
    if 'y_range' in plot_settings: ax0.set_ylim(plot_settings['y_range'] )
    if len(listOfYArrays) == 3: #This generally means observed, mu_guess, map, in that order.
        if len(x_values) > 1: #This means there are enough data to make lines.        
            for seriesIndex in range(len(listOfYArrays)):           
                ax0.plot(x_values,listOfYArrays[0],'g')
                if len(listOfYUncertaintiesArrays) >= 1: #If length is >=1, uncertainties for first data set.
                    ax0.errorbar(x_values, listOfYArrays[0], yerr=np.array(listOfYUncertaintiesArrays[0]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for experiments data series, just a line..
                ax0.plot(x_values,listOfYArrays[1], '#00A5DF')
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[1], yerr=np.array(listOfYUncertaintiesArrays[1]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.
                ax0.plot(x_values,listOfYArrays[2], 'r') 
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[2], yerr=np.array(listOfYUncertaintiesArrays[2]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.
                    
        if len(x_values) == 1: #This means there are single points, and we need to make symbols, by adding an "o".
            for seriesIndex in range(len(listOfYArrays)):           
                ax0.plot(x_values,listOfYArrays[0],'go')
                if len(listOfYUncertaintiesArrays) >= 1: #If length is >=1, uncertainties for first data set.
                    ax0.errorbar(x_values, listOfYArrays[0], yerr=np.array(listOfYUncertaintiesArrays[0]).flatten(), fmt='o', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for experiments data series, just a line.
                ax0.plot(x_values,listOfYArrays[1], 'co')
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[1], yerr=np.array(listOfYUncertaintiesArrays[1]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
                ax0.plot(x_values,listOfYArrays[2], 'ro') 
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[2], yerr=np.array(listOfYUncertaintiesArrays[2]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.
    elif len(listOfYArrays) == 4: #This generally means observed, mu_guess, map, mu_app
        if len(x_values) > 1: #This means there are enough data to make lines.        
            for seriesIndex in range(len(listOfYArrays)):
                ax0.plot(x_values,listOfYArrays[0],'g')
                if len(listOfYUncertaintiesArrays) >= 1: #If length is >=1, uncertainties for first data set.
                    ax0.errorbar(x_values, listOfYArrays[0], yerr=np.array(listOfYUncertaintiesArrays[0]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for experiments data series, just a line.
                ax0.plot(x_values,listOfYArrays[1], '#00A5DF')
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[1], yerr=np.array(listOfYUncertaintiesArrays[1]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                                             
                ax0.plot(x_values,listOfYArrays[2], 'r') 
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[2], yerr=np.array(listOfYUncertaintiesArrays[2]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
                ax0.plot(x_values,listOfYArrays[3], 'k')  #k is black.
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[3], yerr=np.array(listOfYUncertaintiesArrays[3]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
        if len(x_values) == 1: #This means there are single points, and we need to make symbols, by adding an "o".
                ax0.plot(x_values,listOfYArrays[0],'go')
                if len(listOfYUncertaintiesArrays) >= 1: #If length is >=1, uncertainties for first data set.
                    ax0.errorbar(x_values, listOfYArrays[0], yerr=np.array(listOfYUncertaintiesArrays[0]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for experiments data series, just a line.
                ax0.plot(x_values,listOfYArrays[1], 'co')
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[1], yerr=np.array(listOfYUncertaintiesArrays[1]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
                ax0.plot(x_values,listOfYArrays[2], 'ro') 
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[2], yerr=np.array(listOfYUncertaintiesArrays[2]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
                ax0.plot(x_values,listOfYArrays[3], 'ko')  #k is black. https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.plot.html#matplotlib.pyplot.plot
                if len(listOfYUncertaintiesArrays) > 1: #If length is >1, uncertainties for all data sets
                    ax0.errorbar(x_values, listOfYArrays[2], yerr=np.array(listOfYUncertaintiesArrays[3]).flatten(), fmt='.', barsabove=False, markersize=0, linewidth=error_linewidth, color="gray", ecolor="lightgray") #markersize=0 because we want no marker for this.                    
    else:
        if len(x_values) > 1: #This means there are enough data to make lines.        
            for seriesIndex in range(len(listOfYArrays)):
                ax0.plot(x_values[0],listOfYArrays[seriesIndex])
        if len(x_values) == 1: #This means there are single points, and we need to make symbols.
            for seriesIndex in range(len(listOfYArrays)):
                ax0.plot(x_values[0],listOfYArrays[seriesIndex], 'o')            
    if plot_settings['legend'] == True:
        ax0.legend(plot_settings['legendLabels']) #legends must be after plots are made.
    fig0.tight_layout()
    if exportFigure==True:
        fig0.savefig(directory + plot_settings['figure_name'] + '.png', dpi=plot_settings['dpi'])
    if showFigure==False:
        plt.close(fig0)
    return fig0

def makeTrisurfacePlot(xValues, yValues, zValues, exportFigure = True, figure_name="TrisurfacePlot", showFigure=True, directory=''):
    from mpl_toolkits.mplot3d import Axes3D #Although it does not look like this is called here, the InfoGain plots will fail without this line.
    fig1, ax1 =  plt.subplots(1)
    ax1 = plt.axes(projection ='3d')
    image = ax1.plot_trisurf(xValues,yValues, zValues)
    if exportFigure == True:
        fig1.savefig(directory + figure_name + '.png')
    if showFigure==False:
        plt.close(fig1)
    return fig1, ax1, image

def makeMeshGridSurfacePlot(XX, YY, ZZ,  plot_settings = {}, exportFigure = True, figure_name="MeshGridSurfacePlot", showFigure=True, directory=''):
    #TODO: plot_settings should be used for axis labels etc, like above.
    #TODO: create a UserInput variable named info_gain_plot_settings (like what the other cases have).
    from mpl_toolkits.mplot3d import Axes3D #I am not sure if this line is needed here, it might b.
    fig1,ax1 = plt.subplots(figsize=(5,5))
    #ax = fig.add_subplot(111, projection='3d')
    surf = ax1.pcolor(XX,YY,ZZ,cmap=matplotlib.cm.coolwarm)
    # ax1.set_xlabel('Temperature (K)')
    # ax1.set_ylabel(r'$p_A$')
    #ax.set_zlabel('Information Gain')
    #ax1.set_xticks(temperatures)
    #ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.set_title('Information Gain Surface')
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    if exportFigure==True:
        fig1.savefig(directory + figure_name + '.png')
    if showFigure==False:
        plt.close(fig1)
    return fig1, ax1#, image
    

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    Sampled from Guillaume post on https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib/53865762#53865762
    """
    from matplotlib import cm
    from matplotlib.colors import Normalize 
    from scipy.interpolate import interpn
    if ax is None:
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Density')
    return fig, ax                   


def createScatterPlot(data_a, data_b, a_tuple, b_tuple, graphs_directory, plot_settings):
    """Generates and saves a scatter matrix plot.

    :param data_a: parameter points in first for loop (:type: pd.Series)
    :param data_b: parameter points in second for loop (:type: pd.Series)
    :param a_tuple: set of parameter names, MAP, muAP, and initial poitns (:type: tuple)
    :param b_tuple: set of parameter names, MAP, muAP, and initial poitns (:type: tuple)
    :param graphs_directory: path to graphs directory (:type: str)
    :param plot_settings: plot settings from User Input (:type: dict) https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/UserInput.py
    """
    import matplotlib.pyplot as plt
    point_plot_settings = (plot_settings['sampled_point_sizes'], plot_settings['sampled_point_transparency'])
    cross_plot_settings = (plot_settings['cross_marker_size'], plot_settings['cross_marker_transparency'])
    fig = plt.figure()
    # create a scatter plot of the posterior data between two parameters
    plt.scatter(data_a, data_b, s=point_plot_settings[0], alpha=point_plot_settings[1])
    # Allows the user to have no crosses if size is set to 0
    if cross_plot_settings[0] != 0:
        cross_size = cross_plot_settings[0]
        cross_transparency = cross_plot_settings[1]
        # creates the extra points for the MAP, muAP, and initial point denoted by crosses
        plt.scatter(a_tuple[2], b_tuple[2], s=cross_size, alpha=cross_transparency, c='r', marker='x') 
        plt.scatter(a_tuple[3], b_tuple[3], s=cross_size, alpha=cross_transparency, c='k', marker='x') 
        plt.scatter(a_tuple[4], b_tuple[4], s=cross_size, alpha=cross_transparency, c='#00A5DF', marker='x')
    # create labels and save the image to the graphs directory
    plt.xlabel(a_tuple[0], fontsize=plot_settings['fontsize'])
    plt.ylabel(b_tuple[0], fontsize=plot_settings['fontsize'])
    if plot_settings['max_num_x_ticks'] != 'auto' and isinstance(plot_settings['max_num_x_ticks'], int):
        plt.locator_params(axis='x', nbins=plot_settings['max_num_x_ticks'])
    if plot_settings['max_num_y_ticks'] != 'auto' and isinstance(plot_settings['max_num_y_ticks'], int):
        plt.locator_params(axis='y', nbins=plot_settings['max_num_y_ticks'])
    fig.savefig(graphs_directory+f'Scatter_{a_tuple[1]}_{b_tuple[1]}',dpi=plot_settings['dpi'])
    plt.close(fig)

def createScatterHeatMapPlot(data_a, data_b, a_tuple, b_tuple, graphs_directory, plot_settings):
    """Generates and saves a scatter heat map plot.

    :param data_a: parameter points in first for loop (:type: pd.Series)
    :param data_b: parameter points in second for loop (:type: pd.Series)
    :param a_tuple: set of parameter names, MAP, muAP, and initial poitns (:type: tuple)
    :param b_tuple: set of parameter names, MAP, muAP, and initial poitns (:type: tuple)
    :param graphs_directory: path to graphs directory (:type: str)
    :param plot_settings: plot settings from User Input (:type: dict) https://github.com/AdityaSavara/PEUQSE/blob/master/PEUQSE/UserInput.py
    """
    import matplotlib.pyplot as plt
    point_plot_settings = (plot_settings['sampled_point_sizes'], plot_settings['sampled_point_transparency'])
    cross_plot_settings = (plot_settings['cross_marker_size'], plot_settings['cross_marker_transparency'])
    # create a scatter plot of the posterior data between two parameters
    fig, ax = density_scatter(data_a, data_b, s=point_plot_settings[0], alpha=point_plot_settings[1])
    # Allows the user to have no crosses if size is set to 0
    if cross_plot_settings[0] != 0:
        cross_size = cross_plot_settings[0]
        cross_transparency = cross_plot_settings[1]
        # creates the extra points for the MAP, muAP, and initial point denoted by crosses
        ax.scatter(a_tuple[2], b_tuple[2], s=cross_size, alpha=cross_transparency, c='r', marker='x') 
        ax.scatter(a_tuple[3], b_tuple[3], s=cross_size, alpha=cross_transparency, c='k', marker='x') 
        ax.scatter(a_tuple[4], b_tuple[4], s=cross_size, alpha=cross_transparency, c='#00A5DF', marker='x')
    # create labels and save the image to the graphs directory
    ax.set_xlabel(a_tuple[0], fontsize=plot_settings['fontsize']) # , fontsize=plot_settings['fontsize']
    ax.set_ylabel(b_tuple[0], fontsize=plot_settings['fontsize']) # , fontsize=plot_settings['fontsize']
    if plot_settings['max_num_x_ticks'] != 'auto' and isinstance(plot_settings['max_num_x_ticks'], int):
        plt.locator_params(axis='x', nbins=plot_settings['max_num_x_ticks'])
    if plot_settings['max_num_y_ticks'] != 'auto' and isinstance(plot_settings['max_num_y_ticks'], int):
        plt.locator_params(axis='y', nbins=plot_settings['max_num_y_ticks'])
    fig.savefig(graphs_directory+f'Heat_Scatter_{a_tuple[1]}_{b_tuple[1]}',dpi=plot_settings['dpi'])
    plt.close(fig)

def createAutoCorrPlot(N, taus, param_name, param_symbol, graphs_directory):
    """
    Creates Integrated Autocorrelated Time Plots to show MCMC convergence.
    Convergence can be inferred when the AutoCorrelatedTime converges.

    :param N: Indices for the window size. (:type: np.array)
    :param taus: Autocorrelated Time values. (:type: np.array)
    :param param_name: Parameter name. (:type: str)
    :param param_symbol: Parameter symbol. Usually in unicode raw string. (:type: str)
    :param graphs_directory: Directory for storing graphs. (:type: str)
    """
    # create loglog plot
    fig = plt.figure()
    plt.loglog(N, taus, "o-", label='AutoCorr')
    # get current y boundaries
    ylim = plt.gca().get_ylim()
    # plot 50*tau heuristical line for possible convergence
    plt.plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
    # reset the y boundary
    plt.ylim(ylim)
    # create labels and legend
    plt.xlabel("number of samples, $N$")
    plt.ylabel(r"$\tau$ estimates")
    plt.title(f'Integrate Autocorrelation Time for {param_symbol}')
    plt.legend(fontsize=14)
    # save and close figure
    plt.savefig(graphs_directory+'AutoCorrelationPlot_'+param_name)
    plt.close(fig)

def createGewekePlot(z_scores, N, z_percents, param_name, param_symbol, graphs_directory):
    """
    Creates Gewekes indicies plot of entire post burn in samples
    and a plot for percentage of points that fall outside 1 std when compared to the 
    first 10% of the relative window size to the last 50%.

    :param z_scores: z scores for overall Geweke plot. (:type: np.array)
    :param N: Indices for the window size. (:type: np.array)
    :param z_percents: Percent of z scores outside 1 std for each window size. (:type: np.array)
    :param param_name: Parameter name. (:type: str)
    :param param_symbol: Parameter symbol. Usually in unicode raw string. (:type: str)
    :param graphs_directory: Directory for storing graphs. (:type: str)
    """
    # create 2 subplots for one figure
    fig, (ax1, ax2) = plt.subplots(1, 2, tight_layout=True, squeeze=True)
    plt.suptitle(f'Geweke Diagnostics for {param_symbol}')
    # plot Geweke plot on first graph
    ax1.scatter(*z_scores)
    # plot 1 std line. since this is abs of values only positive is needed.
    ax1.hlines([1], 0, round(z_scores[0][-1]), linestyles='dotted')
    # set labels and x boundaries
    ax1.set_xlabel('Last 50% Samples')
    ax1.set_ylabel('Z scores')
    ax1.set_xlim(0, round(z_scores[0][-1]))
    # plot Geweke percents
    ax2.plot(N, z_percents, "o-")
    # set labels, save figure, and close
    ax2.set_xlabel('Sample Indices')
    ax2.set_ylabel('Percent Outside '+u"\u00B1"+"1"+u"\u03C3")
    fig.savefig(graphs_directory+'GewekeDiagnostic_'+param_name)
    plt.close(fig)