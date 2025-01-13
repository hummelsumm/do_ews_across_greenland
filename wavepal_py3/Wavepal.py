import numpy as np

#__path__ = ["/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/wavepal_test/wavepal_py3"]
class Wavepal:

	def __init__(self,t,mydata,t_axis_label="",mydata_axis_label="",t_units=None,mydata_units=None):
		
		""" Constructor of Wavepal class. It Initializes all the variables of Wavepal class with an access outside Wavepal (the user can access them). The list of available variables is given in this function.
			Required Inputs:
			- t [1-dim numpy array of floats]: the times of the time series, distinct and in ascending order.
			- mydata [1-dim numpy array of floats - size=time.size]: the data at the times given by 't'.
			Optional Inputs:
			- t_axis_label="": label for the time axis in figures
			- mydata_axis_label="": label for the data axis in figures
			- t_units=None: units of 't' (string type)
			- mydata_units=None: units of 'mydata' (string type).
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		# check inputs
		try:
			assert np.issubsctype(t,float)
		except:
			print("Error at input 't': must be a numpy array of 'float' type")
			return
		try:
			assert np.issubsctype(mydata,float)
		except:
			print("Error at input 'mydata': must be a numpy array of 'float' type")
			return
		try:
			assert type(t_axis_label) is str
		except AssertionError:
			print("Error at input 't_axis_label': must be of 'str' type")
			return
		try:
			assert type(mydata_axis_label) is str
		except AssertionError:
			print("Error at input 'mydata_axis_label': must be of 'str' type")
			return
		try:
			assert (t_units is None) or (type(t_units) is str)
		except AssertionError:
			print("Error at input 't_units': must be None or of 'str' type")
			return
		try:
			assert (mydata_units is None) or (type(mydata_units) is str)
		except AssertionError:
			print("Error at input 'mydata_units': must be None or of 'str' type")
			return
		# Initializes all the variables the user may have access when using Wavepal
		self.t=t
		self.mydata=mydata
		self.t_axis_label=t_axis_label
		self.mydata_axis_label=mydata_axis_label
		self.t_units=t_units
		self.mydata_units=mydata_units
		if t_units is None:
			self.t_label=""
			self.freq_label=""
		else:
			self.t_label=" ("+t_units+")"
			self.freq_label=" ("+t_units+"${}^{-1}$)"
		if mydata_units is None:
			self.mydata_label=""
			self.power_label=""
			self.varpow_label=""
		else:
			self.mydata_label=" ("+mydata_units+")"
			self.power_label=" ("+mydata_units+"${}^2$)"
			self.varpow_label=" ("+mydata_units+"${}^4$)"
		# Defined in function 'check_data'
		self.nt=None
		self.run_check_data=False
		# Defined in function 'plot_timestep'
		self.dt=None
		# Defined in function 'choose_trend_degree'
		self.pol_degree=None
		self.trend=None
		self.run_choose_trend_degree=False
		# Defined in function 'carma_params'
		self.p=None
		self.q=None
		self.signif_level_type=""
		self.nmcmc=None
		self.mylength=None
		self.beta_gam=None
		self.sigwn_unique=None
		self.alpha_unique=None
		self.beta_unique=None
		self.ARMA_mat_unique=None
		self.myn=None					# contains the carma(p,q) MCMC time series. It remains 'None' if p=q=0
		self.run_carma_params=False
		# Defined in function 'trend_vectors'
		self.myprojvec=None
		self.Vmat=None
		self.run_trend_vectors=False
		# Defined in function 'freq_analysis'
		self.freq=None
		self.tau=None
		self.myind_time=None
		self.myind_freq=None
		self.myind_Q=None
		self.D=None
		self.nsmooth_vec=None
		self.nsmooth=None
		self.tapwindow=None
		self.weighted_WOSA=None
		self.computes_amplitude=None
		self.n_moments=None
		self.percentile=None
		self.periodogram=None
		self.periodogram_cl_mcmc=None
		self.periodogram_cl_anal=None
		self.f_periodogram=None
		self.f_periodogram_cl=None
		self.amplitude=None
		self.amplitude_cos=None
		self.amplitude_sin=None
		self.pseudo_spectrum_mcmc=None
		self.pseudo_spectrum_anal=None
		self.variance_anal=None
		self.periodogram_cl_anal_check_percentile=None
		self.run_freq_analysis=False
		# Defined in function 'freq_filtering'
		self.freq_filtered_signal=None
		self.freq_filtered_signal_bounds=None
		self.run_freq_filtering=False
		# Defined in function 'timefreq_analysis'
		self.w0=None
		self.theta=None
		self.period_cwt=None
		self.period_ampl=None
		self.smoothing_coeff=None
		self.weighted_CWT=None
		self.shannonnyquistexclusionzone=None
		self.coi1=None
		self.coi2=None
		self.coi1_smooth=None
		self.coi2_smooth=None
		self.perlim1_smooth_cwt=None
		self.perlim1_smooth_ampl=None
		self.perlim2_smooth_scal=None
		self.perlim2_smooth_ampl=None
		self.computes_cwtamplitude=None
		self.n_moments_cwt=None
		self.percentile_cwt=None
		self.scalogram=None
		self.scalogram_cl_mcmc=None
		self.scalogram_cl_anal=None
		self.cwtamplitude=None
		self.cwtamplitude_cos=None
		self.cwtamplitude_sin=None
		self.pseudo_cwtspectrum_mcmc=None
		self.pseudo_cwtspectrum_anal=None
		self.cwt_variance_anal=None
		self.scalogram_cl_anal_check_convergence=None
		self.computes_global_scalogram=None
		self.global_scalogram=None
		self.global_amplitude=None
		self.pseudo_global_spectrum_mcmc=None
		self.global_scalogram_cl_mcmc=None
		self.pseudo_global_spectrum_anal=None
		self.global_scalogram_variance_anal=None
		self.global_scalogram_cl_anal=None
		self.global_scalogram_cl_anal_check_convergence=None
		self.minscal=None
		self.maxscal=None
		self.minampl=None
		self.maxampl=None
		self.minampl_sq=None
		self.maxampl_sq=None
		self.min_pseudo_cwtspectrum_anal=None
		self.max_pseudo_cwtspectrum_anal=None
		self.min_pseudo_cwtspectrum_mcmc=None
		self.max_pseudo_cwtspectrum_mcmc=None
		self.min_cwt_variance_anal=None
		self.max_cwt_variance_anal=None
		self.n_outside_scalelim1=None
		self.weight_cwt=None
		self.r=None
		self.run_timefreq_analysis=False
		# Defined in function 'timefreq_ridges_filtering'
		self.skeleton=None
		self.run_timefreq_ridges_filtering=False
		# Defined in function 'timefreq_band_filtering'
		self.timefreq_band_filtered_signal=None
		self.timefreq_band_filtered_signal_bounds=None
		self.run_timefreq_band_filtering=False
		


	def check_data(self):
		
		""" check_data checks the data, and modifies them if needed. 
				-> check for NaN or Inf values. Returns an error if NaN or Inf found.
				-> Put the time axis in ascending order and check that all the times are distinct. Correction if needed.
				-> Check that there are at least 50 data points. Returns an error if this condition is not fulfilled.
			Inputs:
			/
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		from .dt_min import dt_min
		from .distinct_ages import distinct_ages

		# Check that the data are not Nan or Inf
		for k in range(self.t.size):
			if np.isnan(self.t[k]) or np.isnan(self.mydata[k]) or np.isinf(self.t[k]) or np.isinf(self.mydata[k]):
				print("Error: Nan or Inf value at time/age ", self.t[k])
				return
		# Put the time axis in ascending order and check that all the times are distinct - Correction if needed
		tind=np.argsort(self.t)
		self.t=self.t[tind]
		self.mydata=self.mydata[tind]
		try:
			assert dt_min(self.t)>0.0
		except AssertionError:
			print("WARNING: The times/ages of the time series must be distinct")
			print("The program will automatically select the first occurrence of the time/age (and skip the others) where they are the same (times/ages are previously set in ascending order)")
			self.t, self.mydata=distinct_ages(self.t,self.mydata)
		# Check that we have at least 50 data points
		self.nt=self.t.size
		try:
			assert self.nt>=50 # useful for the ACF of residual noise and its confidence levels (in CARMA pack)
		except AssertionError:
			print("Error: Not enough data points - please provide at least 50 data points")
			return
		self.run_check_data=True



	def plot_timestep(self,hist=True,nbins=10,log_yaxis=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,reverse_xaxis=False):
	
		""" plot_timestep computes and generates the figure of the time step in function of time.
			Optional Inputs:
			- hist=True: if True, draws the histogram of the distribution of the time steps.
			- nbins=10: number of bins for the histogram.
			- log_yaxis=False: If True, the vertical axis is in log scale.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- reverse_xaxis=False: Reverse the horizontal axis if True
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert type(hist) is bool
		except AssertionError:
			print("Error at input 'hist': must be True or False")
			return
		try:
			assert (type(nbins) is int) and (nbins>0)
		except AssertionError:
			print("Error at input 'nbins': must be of type 'int' and >0")
			return
		try:
			assert type(log_yaxis) is bool
		except AssertionError:
			print("Error at input 'log_yaxis': must be True or False")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print("Error: Must have run function 'check_data'")
			return
		self.dt=np.zeros(self.nt-1)
		for k in range(1,self.nt):
			self.dt[k-1]=self.t[k]-self.t[k-1]
		if log_yaxis is True:
			plt.semilogy(self.t[1:],self.dt,"k.",zorder=1)
		else:
			plt.plot(self.t[1:],self.dt,"k.",zorder=1)
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.xlim(self.t[1], self.t[-1])
		plt.ylabel(self.t_axis_label+" step"+self.t_label,fontsize=fontsize_axes)
		plt.title(self.t_axis_label+" step in function of "+self.t_axis_label,fontsize=fontsize_title)
		if hist is True:
			pdf, bin_edges = np.histogram(self.dt, bins=nbins)
			bin_edges = bin_edges[0:pdf.size]
			# Stretch the PDF so that it is readable on the residual plot when plotted horizontally
			pdf = pdf / float(pdf.max()) * 0.4 * (self.t[-1]-self.t[1])
			# Add the histogram to the plot
			plt.barh(bin_edges, pdf, height=bin_edges[1] - bin_edges[0],left=self.t[1],alpha=0.5,zorder=2)
		plt.tick_params(labelsize=fontsize_ticks)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		return plt



	def plot_trend(self,pol_degree,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='small',linewidth_data=1.0,linewidth_trend=2.0,reverse_xaxis=False,reverse_yaxis=False):
	
		""" plot_trend computes and generates the figure the polynomial trend of the time series.
			Required Inputs:
			- pol_degree [int]: the degree of the polynomial trend. pol_degree=-1 means no trend (trend is imposed to zero). pol_degree=0 means trend=constant=data average. pol_degree=1 means linear detrending (trend=a*t+b). etc.
			Optional Inputs:
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='small': fontsize for the figure legend.
			- linewidth_data=1.0: linewidth for the data
			- linewidth_trend=2.0: linewidth for the trend
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt
		from .detrending import detrending
		
		# check inputs
		try:
			assert (type(pol_degree) is int) and pol_degree>=-1
		except AssertionError:
			print("Error at input 'pol_degree': must be an integer >= -1")
			return
		try:
			assert (type(linewidth_data) is int) or (type(linewidth_data) is float)
		except AssertionError:
			print("Error at input 'linewidth_data': must be an integer or float")
			return
		try:
			assert (type(linewidth_trend) is int) or (type(linewidth_trend) is float)
		except AssertionError:
			print("Error at input 'linewidth_trend': must be an integer or float")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print("Error: Must have run function 'check_data'")
			return
		trend=detrending(self.t,self.mydata,pol_degree)
		# Figure with: raw data + trend
		plt.plot(self.t,self.mydata,"b",label="Data",linewidth=linewidth_data)
		plt.plot(self.t,trend,"r",label="Trend",linewidth=linewidth_trend)
		plt.legend(fancybox=True,fontsize=fontsize_legend)
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel(self.mydata_axis_label+self.mydata_label,fontsize=fontsize_axes)
		plt.title("Data and Trend",fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		return plt
	
	
	
	def choose_trend_degree(self,pol_degree):

		""" choose_trend_degree records the user choice for the degree of the polynomial trend.
			Required Inputs:
			- pol_degree [int]: the degree of the polynomial trend. pol_degree=-1 means no trend (trend is imposed to zero). pol_degree=0 means trend=constant=data average. pol_degree=1 means linear detrending (trend=a*t+b). etc.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		from .detrending import detrending
		
		# check inputs
		try:
			assert (type(pol_degree) is int) and pol_degree>=-1
		except AssertionError:
			print("Error at input 'pol_degree': must be an integer >= -1")
			return
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print("Error: Must have run function 'check_data'")
			return
		self.pol_degree=pol_degree
		self.trend=detrending(self.t,self.mydata,self.pol_degree)
		self.run_choose_trend_degree=True



	def trend_vectors(self):
		
		""" trend_vectors computes some arrays concerning the trend of the time series.
			Inputs:
			/
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import numpy.linalg as la
		import copy
		
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print("Error: Must have run function 'check_data'")
			return
		try:
			assert self.run_choose_trend_degree is True
		except AssertionError:
			print("Error: Must have run function 'choose_trend_degree'")
			return

		tbis=self.t/self.t[-1]    # for numerical stability with the powers of the time, for the polynomial trend
		myprojvec=np.zeros((self.nt,self.pol_degree+3))
		Vmat=np.zeros((self.nt,self.pol_degree+3))
		for o in range(self.pol_degree+1):
			myprojvec_o=tbis**o
			Vmat[:,o]=copy.copy(myprojvec_o)
			for k in range(o):	# Gram-Schmidt
				h=myprojvec[:,k]
				myprojvec_o-=np.dot(h,myprojvec_o)*h
			myprojvec[:,o]=myprojvec_o/la.norm(myprojvec_o)
		self.myprojvec=myprojvec
		self.Vmat=Vmat
		self.run_trend_vectors=True



	def carma_params(self,p=1,q=0,signif_level_type="",min_autocorrelation=0.2,nmcmc=None,nmcmc_carma_max=1000000,make_carma_fig=False,path_to_figure_folder="",nbins=10,maxlag=50,golden_fact=100.,dpi=None,reverse_xaxis=False):
		
		""" carma_params builds all the material related to the CARMA(p,q) background noise. It uses 'carma pack' python package to estimate the parameters of the carma process and to generate MCMC samples. 'carma pack' relies on the following paper:
			B. C. Kelly, A. C. Becker, M. Sobolewska, A. Siemiginowska, and P. Uttley. Flexible and scalable methods for quantifying stochastic variability in the era of massive time-domain astronomical data sets. The Astrophysical Journal, 788(1):33, 2014.
			Note that 'carma pack' is not used when dealing with a white noise (p=q=0).
			Optional Inputs:
			- p=1 and q=0: the orders of the CARMA(p,q) process. Other options are
			-> p=q=0, or
			-> p>=1 and p>q
			- signif_level_type="a". Other options are:
			-> signif_level_type="n"
			-> signif_level_type="an"
			-> signif_level_type=""
			If "a" is in signif_level_type, analytical confidence levels will be computed in a later step, and carma_params builds matrix K. Matrix K is defined in
			'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
			If "n" is in signif_level_type, MCMC confidence levels will be computed in a later step, and carma_params generates MCMC samples of a carma(p,q) time series.
			If signif_level_type="", confidence levels won't be computed in a later step.
			- min_autocorrelation=0.2 [Not used if p=q=0]: value of the autocorrelation for the variables of the posterior carma distribution for which the lag is returned. That lag is the 'decorrelation length' of the posterior distribution, and is used ta skim off the distribution to get decorrelated samples. Typical values are 0.1 or 0.2. For the same final number of samples, the smallest min_autocorrelation is, the longest is the computing time, because a small value for min_autocorrelation means a lot of raw samples to be generated.
			- nmcmc=None: Number of MCMC samples. 
                If "a" is in signif_level_type:
                    If p>0, that's the approximate number of samples for the parameters, that are used to compute the median values. If p=q=0, nmcmc is not used, because this case is fully analytical.
                If "n" is in signif_level_type: 
                    If p>0, that's the approximate number of generated carma(p,q) time series. Note that, due to the skimming off explained above, a bigger number of MCMC samples is first generated. That number may be limited to save computing time (see nmcmc_carma_max below) but in that case, nmcmc may be automatically adpated toward a smaller value. The final value of nmcmc is, indicated at the terminal. If p=q=0, that's the exact number of MCMC generated white noise time series.
                Default initial values: nmcmc=1000 if signif_level_type="a", or nmcmc=10000 if signif_level_type="an" or "n".
			- nmcmc_carma_max=1000000 [Not used if p=q=0]: Maximum number of MCMC samples in 'carma pack'.
			- make_carma_fig=False: generates the figures of posterior distributions and the quality of the fit, most of them coming from 'carma pack'. Change to True to activate it.
			- path_to_figure_folder="": path to the folder where the figures are to be saved. Used if make_carma_fig is True.
			- nbins=10: Number of bins for the histogram of the standradized residuals. Used if make_carma_fig is True.
			- maxlag=50: Displayed max lag for the figure of the ACF of the (squared) residuals. Used if make_carma_fig is True.
			- golden_fact=100. [Only used if p=q=0.]: parameter controlling the interval in which the search is performed in order to find the max. of the posterior distribution of the white noise variance. Increase golden_fact to widen the interval. Typical values are 10., 100., 1000., ...
			- dpi=None: Used for the figures of 'carma pack'. It imposes the figure resolution in dots per inch (dpi). If None it will default to the value savefig.dpi in the matplotlibrc file.
			- reverse_xaxis=False: Used for the figures of 'carma pack'. When the horizontal axis represents the times of the time series, axis is reversed if reverse_xaxis is True.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		
		# check inputs
		try:
			assert (type(signif_level_type) is str) and (signif_level_type.lower()=="")
		except AssertionError:
			print("Error at input 'signif_level_type': must be '' because carmapack doesn't work")
			return
		signif_level_type=signif_level_type.lower()
	
		# check that some functions were previously run
		try:
			assert self.run_check_data is True
		except AssertionError:
			print("Error: Must have run function 'check_data'")
			return
		try:
			assert self.run_choose_trend_degree is True
		except AssertionError:
			print("Error: Must have run function 'choose_trend_degree'")
			return
		self.p=p
		self.q=q
		self.signif_level_type=signif_level_type
		
		self.run_carma_params=True
	


	def freq_analysis(self,freqmin=None,freqmax=None,freqstep=None,dt_GCD=None,freq_min_bound=True,freq_max_bound=True,mywindow=1,D=None,betafact=0.75,coverage=90.,WOSA_segments=None,percentile=None,weighted_WOSA=True,n_moments=10,MaxFunEvals=100000,algo_moments='gamma-polynomial',computes_amplitude=False):

		""" freq_analysis computes the WOSA periodogram and its confidence levels and the amplitude periodogram.
			Optional Inputs:
			- freqmin=None: minimal frequency. Default value is freqmin=1.0/(t[-1]-t[0]) (where t is the time vector)
			- freqmax=None: maximal frequency. Default value is freqmax=1.0/2.0/dt_GCD
			- freqstep=None: frequency step. Default value is freqstep=(freqmax-freqmin)/t.size (where t is the time vector)
			- dt_GCD=None: the greatest common divisor of the time steps of the time vector. Default value is the smallest time step (which is not exactly dt_GCD; it is an upper bound on the estimation of dt_GCD; but it is very simple to compute). 
			- freq_min_bound=True: limit (if True), or not (if False), the lower bound of the frequency range, for each WOSA segment. More details in
			'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
			- freq_max_bound=True: limit (if True), or not (if False), the upper bound of the frequency range, for each WOSA segment. More details in the above cited article.
			- mywindow=1: window choice for the windowing of the WOSA segments.
				-> 1: Square window
				-> 2: Triangular window
				-> 3: sin window
				-> 4: sin**2 (Hanning) window
				-> 5: sin**3 window
				-> 6: sin**4 window
				-> 7: Hamming window, defined as 0.54-0.46*np.cos(2.0*np.pi*time/D) [D is defined below]
				-> 8: 4-term Blackman-Harris window, with a0=0.35875 and a1=0.48829 and a2=0.14128 and a3=0.01168
				-> 9: Kaiser-Bessel window, with parameter alpha=2.5
				-> 10: Gaussian window, with standard dev. sigma=D/6.0 [D is defined below]
				Terminology and formulas come from:
				F. Harris. On the use of windows for harmonic analysis with the discrete fourier transform. Proceedings of the IEEE, 66(1):51-83, January 1978.
			- D=None: the temporal length of the WOSA segments. Default value is D=(t[-1]-t[0])/10.0 (where t is the time vector).
			- betafact=0.75: overlapping factor for the WOSA segments. Must take a value in [0.,1.[
			- coverage=90.: minimal coverage (in percent) of the data points along the segment length. Below this value, a WOSA segment is not considered. Must take a value between 0. and 100.
			- WOSA_segments=None: Choose the minimal number of WOSA segments to be present at each frequency to take it into account, thus defining the frequency range for the analysis.
				-> WOSA_segments='all': No restrictions on the number of segments per frequency.
				-> WOSA_segments='max': Consider only the frequencies for which the number of WOSA segments is maximal. This is the most restrictive case.
				-> WOSA_segments=None: Consider only the frequencies for which the number of WOSA segments is at least 10, or maximal if there are less than 10 segments.
				-> WOSA_segments=n (n is an integer): Consider only the frequencies for which the number of WOSA segments is at least n.
			- percentile=None: The x^th percentiles for the confidence levels. Must be a 1-dim numpy array. Default is the 95^th percentile (i.e. the 95% confidence level):percentile=np.zeros(1); percentile[0]=95.0
			- weighted_WOSA=True: True to get the weighted periodogram, or False for the classical periodogram.
			- n_moments=10: number of conserved moments for the analytical confidence levels. Must be >= 2. Used if "a" is in signif_level_type (see function 'carma_params')
			- MaxFunEvals=100000: "max_nfev" option for "least_squares" - see python help of "scipy.optimize.least_squares". Used if algo_moments='generalized-gamma-polynomial'.
			- algo_moments='gamma-polynomial': Choice for the algorithm for the computation of analytical confidence levels. Used if "a" is in signif_level_type (see function 'carma_params'). algo_moments='gamma-polynomial' or algo='generalized-gamma-polynomial'. 
			- computes_amplitude=False: computes the amplitude periodogram (if True) or not (if False).
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import numpy.linalg as la
		from .dt_min import dt_min
		from tqdm import trange
		from .freq_analysis_prelims import freq_analysis_prelims
		from .LS_WOSA import LS_WOSA
		from .LS_WOSA_and_Ampl import LS_WOSA_and_Ampl
		from scipy.stats import f as fdistr

		# check inputs
		try:
			assert (freqmin is None) or ((type(freqmin) is float) and freqmin>0.)
		except AssertionError:
			print("Error at input 'freqmin': must be None or of 'float' type and >0.")
			return
		try:
			assert (freqmax is None) or ((type(freqmax) is float) and freqmax>0.)
		except AssertionError:
			print("Error at input 'freqmax': must be None or of 'float' type and >0.")
			return
		if (freqmin is not None) and (freqmax is not None):
			try:
				assert freqmax>=freqmin
			except AssertionError:
				print("Error at input 'freqmin' and 'freqmax': must have freqmax>=freqmin")
				return
		try:
			assert (freqstep is None) or ((type(freqstep) is float) and freqstep>0.)
		except AssertionError:
			print("Error at input 'freqstep': must be None or of 'float' type and >0.")
			return
		try:
			assert (dt_GCD is None) or ((type(dt_GCD) is float) and dt_GCD>0.)
		except AssertionError:
			print("Error at input 'dt_GCD': must be None or of 'float' type and >0.")
			return
		try:
			assert type(freq_min_bound) is bool
		except AssertionError:
			print("Error at input 'freq_min_bound': must be True or False")
			return
		try:
			assert type(freq_max_bound) is bool
		except AssertionError:
			print("Error at input 'freq_max_bound': must be True or False")
			return
		try:
			assert (type(mywindow) is int) and (mywindow>=1)
		except AssertionError:
			print("Error at input 'mywindow': must be of 'int' type and >=1")
			return
		try:
			assert (D is None) or ((type(D) is float) and (D>0.))
		except AssertionError:
			print("Error at input 'D': must be None or of 'float' type and >0.")
			return
		try:
			assert (type(betafact) is float) and (betafact>=0.) and (betafact<1.)
		except AssertionError:
			print("Error at input 'betafact': must be of 'float' type and must take a value in [0.,1.[")
			return
		try:
			assert (type(coverage) is float) and (coverage>=0.) and (coverage<100.)
		except AssertionError:
			print("Error at input 'coverage': must be of 'float' type and must take a value between 0. and 100.")
			return
		try:
			assert (WOSA_segments is None) or ((type(WOSA_segments) is str) and WOSA_segments.lower()=="all") or ((type(WOSA_segments) is str) and WOSA_segments.lower()=="max") or (type(WOSA_segments) is int and WOSA_segments>=0)
		except AssertionError:
			print("Error at input 'WOSA_segments': must be None, 'all', 'max' or an integer >0")
			return
		try:
			assert (percentile is None) or np.issubsctype(percentile,float)
		except:
			print("Error at input 'percentile': must be None or a numpy array of 'float' type")
			return
		try:
			assert type(weighted_WOSA) is bool
		except AssertionError:
			print("Error at input 'weighted_WOSA': must be True or False")
			return
		try:
			assert (type(n_moments) is int) and (n_moments>=2)
		except AssertionError:
			print("Error at input 'n_moments': must be of 'int' type and >=2")
			return
		try:
			assert (type(MaxFunEvals) is int) and (MaxFunEvals>=10)
		except AssertionError:
			print("Error at input 'MaxFunEvals': must be of 'int' type and >=10")
			return
		try:
			assert (type(algo_moments) is str) and ((algo_moments.lower()=="gamma-polynomial") or (algo_moments.lower()=="generalized-gamma-polynomial"))
		except AssertionError:
			print("Error at input 'algo_moments': must be 'gamma-polynomial' or 'generalized-gamma-polynomial'")
			return
		algo_moments=algo_moments.lower()
		try:
			assert type(computes_amplitude) is bool
		except AssertionError:
			print("Error at input 'computes_amplitude': must be True or False")
			return
		# check that some functions were previously run
		if ('a' in self.signif_level_type) or ('n' in self.signif_level_type):
			try:
				assert self.run_carma_params is True
			except AssertionError:
				print("Error: Must have run function 'carma_params'")
				return
		else:
			try:
				assert self.run_check_data is True
			except AssertionError:
				print("Error: Must have run function 'check_data'")
				return
			try:
				assert self.run_choose_trend_degree is True
			except AssertionError:
				print("Error: Must have run function 'choose_trend_degree'")
				return
		try:
			assert self.run_trend_vectors is True
		except AssertionError:
			print("Error: Must have run function 'trend_vectors'")
			return
		# Set Default values for input arguments
		if dt_GCD is None:
			dt_GCD=dt_min(self.t)
		if freqmin is None:
			freqmin=1.0/(self.t[-1]-self.t[0])
		if freqmax is None:
			freqmax=1.0/2.0/dt_GCD
		if freqstep is None:
			freqstep=(freqmax-freqmin)/self.t.size
		if D is None:
			D=(self.t[-1]-self.t[0])/10.0
		if percentile is None:
			percentile=np.zeros(1)
			percentile[0]=95.0
		self.percentile=percentile
		# Adjust freqmin and freqmax and builds the frequency vector
		freqmin=np.maximum(1.0/(self.t[-1]-self.t[0]),freqmin)
		freqmax=np.minimum(1.0/2.0/dt_GCD,freqmax)
		freq=np.linspace(freqmin,freqmax,int((freqmax-freqmin)/freqstep)+1)   # way better than np.arange !!!
		# Check that freq.size is at least 6
		try:
			assert freq.size>5 # in order to be able to compute the convergence of the percentiles at 6 frequencies
		except AssertionError:
			print("Error: Not enough frequency points - please provide at least 6 frequency points")
			return
		# Build the WOSA Lomb-Scargle components
		self.freq,self.tau,self.myind_time,self.myind_freq,self.myind_Q,self.D,self.nsmooth_vec,self.nsmooth,weight_WOSA=freq_analysis_prelims(self.t,freq,D,betafact,mywindow,coverage,freq_min_bound,freq_max_bound,self.pol_degree,WOSA_segments,weighted_WOSA)
		self.tapwindow=mywindow
		self.weighted_WOSA=weighted_WOSA
		J=self.freq.size
		# WOSA periodogram: significance levels, data periodogram, amplitude - This is the main loop, over the frequencies
		npercentile=percentile.size
		self.periodogram=np.zeros(J)
		if self.p==0 and self.q==0 and self.nsmooth==1 and self.tapwindow==2 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1:
			compl_spectrum_data=np.zeros(J)
			norm_sq_mydata=(la.norm(self.mydata)**2)
			if weighted_WOSA is True:
				corr_WOSA_weight=np.sqrt(float(self.nt)/2.)
			else:
				corr_WOSA_weight=1.
		if computes_amplitude is True:
			self.amplitude=np.zeros(J)
			if (self.nsmooth==1 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1):
				self.amplitude_cos=np.zeros(J)
				self.amplitude_sin=np.zeros(J)
		self.computes_amplitude=computes_amplitude
		print("Main loop, over the frequencies:")
		for k in trange(J):
			if computes_amplitude is True:
				Amplitude, M2=LS_WOSA_and_Ampl(self.t,self.mydata,self.freq[k],k,self.myprojvec,self.Vmat,self.D,self.tau,self.nsmooth,self.nsmooth_vec[k],self.myind_time,self.myind_freq,self.tapwindow,self.pol_degree,weight_WOSA)
			else:
				M2=LS_WOSA(self.t,self.freq[k],k,self.myprojvec,self.D,self.tau,self.nsmooth,self.nsmooth_vec[k],self.myind_time,self.myind_freq,self.tapwindow,self.pol_degree,weight_WOSA)
			# Data WOSA periodogram
			self.periodogram[k]=la.norm(np.dot(np.transpose(self.mydata),M2))**2
			# Data F-periodogram
			if self.p==0 and self.q==0 and self.nsmooth==1 and self.tapwindow==2 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1:
				self.myprojvec[:,self.pol_degree+1]=M2[:,0]*corr_WOSA_weight
				self.myprojvec[:,self.pol_degree+2]=M2[:,1]*corr_WOSA_weight
				compl_spectrum_data[k]=(norm_sq_mydata-la.norm(np.dot(np.transpose(self.mydata),self.myprojvec))**2)/corr_WOSA_weight**2
			if computes_amplitude is True:
				# Data Amplitude
				self.amplitude[k]=np.sum(la.norm(Amplitude,axis=0)**2)
				if (self.nsmooth==1 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1):
					self.amplitude_cos[k]=Amplitude[0,0]
					self.amplitude_sin[k]=Amplitude[1,0]
				self.amplitude[k]/=float(self.nsmooth_vec[k])
			self.periodogram[k]/=float(self.nsmooth_vec[k])
		if self.p==0 and self.q==0 and self.nsmooth==1 and self.tapwindow==2 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1:
			self.f_periodogram=float(self.nt-self.pol_degree-3)*self.periodogram/compl_spectrum_data/2.0
			# signif levels
			self.f_periodogram_cl=np.zeros((J,npercentile))
			for k in range(npercentile):
				self.f_periodogram_cl[:,k]=fdistr.ppf(percentile[k]/100.0,2,self.nt-self.pol_degree-3)
		if computes_amplitude is True:
			self.amplitude=np.sqrt(self.amplitude)
		# Linear combination of chi squares: approx at the n_moment order
		# and all the moments are computed at 6 particular frequencies, in order to visualize the convergence
		self.run_freq_analysis=True



	def freq_filtering(self,freq_bounds):

		""" freq_filtering filters the data in a given frequency band. The amplitude periodogram must have been computed and without smoothing (i.e. only with 1 WOSA segment).
			Required Inputs:
			- freq_bounds: a list of tuples. Each tuple contains the two frequencies delimitating the band. Example: [(1./23.,1./19.), (1./43.,1./37.), (1./105.,1.95.)] will filter in three bands. First one: between periods 19. and 23., second one: between periods 37. and 43. and third one: between periods 95. and 105.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		# check inputs
		try:
			assert type(freq_bounds) is list
		except AssertionError:
			print("Error at input 'freq_bounds': must be of type 'list'")
			return
		for k in range(len(freq_bounds)):
			try:
				assert type(freq_bounds[k]) is tuple
			except AssertionError:
				print("Error at input 'freq_bounds': must be a list containing tuples")
				return
			try:
				assert (type(freq_bounds[k][0]) is float) and (type(freq_bounds[k][1]) is float)
			except:
				print("Error at input 'freq_bounds': must be a list containing tuples, with each tuple containing 2 floats")
				return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		try:
			assert self.computes_amplitude is True
		except AssertionError:
			print("Error: In function 'freq_analysis', amplitude must have been computed")
			return
		try:
			assert (self.nsmooth==1 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1)
		except AssertionError:
			print("Error: In function 'freq_analysis', WOSA smoothing is not compatible with filtering")
			return

		self.freq_filtered_signal=np.zeros((self.nt,len(freq_bounds)))
		self.freq_filtered_signal_bounds=np.zeros((len(freq_bounds),2))
		for l in range(len(freq_bounds)):
			ind_low=np.argmin(np.absolute(self.freq-freq_bounds[l][0]))
			ind_high=np.argmin(np.absolute(self.freq-freq_bounds[l][1]))
			try:
				assert ind_high>=ind_low
			except AssertionError:
				print("Error in 'freq_bounds' input parameters: must have frequency_high>=frequency_low")
				return
			for k in range(ind_low,(ind_high+1)):
				self.freq_filtered_signal[:,l]+=self.amplitude_cos[k]*np.cos(2.*np.pi*self.freq[k]*self.t)+self.amplitude_sin[k]*np.sin(2.*np.pi*self.freq[k]*self.t)
			self.freq_filtered_signal_bounds[l,0]=self.freq[ind_low]
			self.freq_filtered_signal_bounds[l,1]=self.freq[ind_high]

		self.run_freq_filtering=True



	def plot_WOSA_segmentation(self,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='small',linewidth_data=1.0,linewidth_segment=1.0,reverse_xaxis=False,reverse_yaxis=False):
		
		""" plot_WOSA_segmentation generates the figure of the WOSA segmentation, superposed to the time series.
			Optional Inputs:
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='small': fontsize for the figure legend.
			- linewidth_data=1.0: linewidth for the data
			- linewidth_segment=1.0: linewidth for the segments
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert (type(linewidth_data) is int) or (type(linewidth_data) is float)
		except AssertionError:
			print("Error at input 'linewidth_data': must be an integer or float")
			return
		try:
			assert (type(linewidth_segment) is int) or (type(linewidth_segment) is float)
		except AssertionError:
			print("Error at input 'linewidth_segment': must be an integer or float")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		mindata=np.amin(self.mydata)
		maxdata=np.amax(self.mydata)
		stepdata=(maxdata-mindata)/100.0
		plt.plot(self.t,self.mydata,"b--.",label="Data",linewidth=linewidth_data)
		counts=0
		countd=0
		for l in range(self.nsmooth):
			myt=[self.tau[l], self.tau[l]+self.D]
			myt_height=mindata-(float(self.myind_Q[l])+1.0)*stepdata
			y=[myt_height, myt_height]
			if (self.myind_freq[l].size==self.freq.size):
				counts+=1
				if (counts==1):
					plt.plot(myt,y,"r",label="All the frequencies are present",linewidth=linewidth_segment)    # solid line if all the frequencies are covered
				else:
					plt.plot(myt,y,"r",linewidth=linewidth_segment)    # solid line if all the frequencies are covered
			else:
				countd+=1
				if (countd==1):
					plt.plot(myt,y,"r:",label="Not all the frequencies are present",linewidth=linewidth_segment)    # dashed line if not all the frequencies are covered
				else:
					plt.plot(myt,y,"r:",linewidth=linewidth_segment)      # dashed line if not all the frequencies are covered
		plt.legend(fancybox=True,fontsize=fontsize_legend)
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel(self.mydata_axis_label+self.mydata_label,fontsize=fontsize_axes)
		plt.title("WOSA Segmentation", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		return plt



	def plot_number_WOSA_segments(self,xaxis="frequency",fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,linewidth=1.0):
	
		""" plot_number_WOSA_segments generates the figure of the number of WOSA segments per frequency.
			Optional Inputs:
			- xaxis="frequency": choice for the horizontal figure axis: "frequency" or "period" (=1/frequency).
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- linewidth=1.0: linewidth.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert (type(linewidth) is int) or (type(linewidth) is float)
		except AssertionError:
			print("Error at input 'linewidth': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		if xaxis=="period":
			mypow=-1
		elif xaxis=="frequency":
			mypow=1
		plt.plot(self.freq**mypow,self.nsmooth_vec,linewidth=linewidth)
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("# segments",fontsize=fontsize_axes)
		if xaxis=="period":
			plt.title("Number of WOSA segments per period", fontsize=fontsize_title)
		elif xaxis=="frequency":
			plt.title("Number of WOSA segments per frequency", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt



	def plot_periodogram(self,xaxis="frequency",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='xx-small',linewidth_per=1.0,linewidth_cl=1.0):
		
		""" plot_periodogram generates the figure of the peridogram and its confidence levels.
			Optional Inputs:
			- xaxis="frequency": choice for the horizontal figure axis: "frequency" or "period" (=1/frequency).
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='small': fontsize for the figure legend.
			- linewidth_per=1.0: linewidth for the periodogram.
			- linewidth_cl=1.0: linewidth for the confidence levels.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth_per) is int) or (type(linewidth_per) is float)
		except AssertionError:
			print("Error at input 'linewidth_per': must be an integer or float")
			return
		try:
			assert (type(linewidth_cl) is int) or (type(linewidth_cl) is float)
		except AssertionError:
			print("Error at input 'linewidth_cl': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		if xaxis=="period":
			mypow=-1
		elif xaxis=="frequency":
			mypow=1
		if loglog is False:
			plt.plot(self.freq**mypow,self.periodogram,"k",label="Data periodogram",linewidth=linewidth_per)
		elif loglog is True:
			plt.loglog(self.freq**mypow,self.periodogram,"k",label="Data periodogram",linewidth=linewidth_per)
		
		plt.legend(fancybox=True,fontsize=fontsize_legend,bbox_to_anchor=(1.1, 1.05))
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		if self.signif_level_type!="" and self.percentile.size>0:
			plt.suptitle("WOSA periodogram and Confidence levels", fontsize=fontsize_title)
		else:
			plt.suptitle("WOSA periodogram", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt
			


	def plot_f_periodogram(self,xaxis="frequency",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='small',linewidth_per=1.0,linewidth_cl=1.0):

		""" plot_f_periodogram generates the figure of the f-peridogram and its confidence levels. Only available if p=q=0.
			Optional Inputs:
			- xaxis="frequency": choice for the horizontal figure axis: "frequency" or "period" (=1/frequency).
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='small': fontsize for the figure legend.
			- linewidth_per=1.0: linewidth for the f-periododgram.
			- linewidth_cl=1.0: linewidth for the confidence levels.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			Note: If you just want the confidence levels for the F-periodogram, in the method 'carma_params' choose signif_level_type="" and p=q=0. 
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth_per) is int) or (type(linewidth_per) is float)
		except AssertionError:
			print("Error at input 'linewidth_per': must be an integer or float")
			return
		try:
			assert (type(linewidth_cl) is int) or (type(linewidth_cl) is float)
		except AssertionError:
			print("Error at input 'linewidth_cl': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		try:
			assert self.p==0 and self.q==0 and self.nsmooth==1 and self.tapwindow==2 and self.myind_time[0,0]==0 and self.myind_time[0,1]==self.nt-1
		except AssertionError:
			print("Error: plot_f_periodogram cannot be applied")
			return
		if xaxis=="period":
			mypow=-1
		elif xaxis=="frequency":
			mypow=1
		if loglog is False:
			plt.plot(self.freq**mypow,self.f_periodogram,"k",label="Data F-periodogram",linewidth=linewidth_per)
		elif loglog is True:
			plt.loglog(self.freq**mypow,self.f_periodogram,"k",label="Data F-periodogram",linewidth=linewidth_per)
		for k in range(self.percentile.size):
			if loglog is False:
				plt.plot(self.freq**mypow,self.f_periodogram_cl[:,k],label="CL at "+str(self.percentile[k])+"%",linewidth=linewidth_cl)
			elif loglog is True:
				plt.loglog(self.freq**mypow,self.f_periodogram_cl[:,k],label="CL at "+str(self.percentile[k])+"%",linewidth=linewidth_cl)
		plt.legend(fancybox=True,fontsize=fontsize_legend)
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Relative Power",fontsize=fontsize_axes)
		if self.percentile.size>0:
			plt.title("F-periodogram and Confidence levels", fontsize=fontsize_title)
		else:
			plt.title("F-periodogram", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt
			




	def plot_amplitude(self,xaxis="frequency",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,plot_band_filtering=False,linewidth=1.0):
		
		""" plot_amplitude generates the figure of the amplitude periodogram. Only available if the option computes_amplitude is True in 'freq_analysis' function.
			Optional Inputs:
			- xaxis="frequency": choice for the horizontal figure axis: "frequency" or "period" (=1/frequency).
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- plot_band_filtering=False: draws the bands where filtering was performed. Must have previously run 'freq_filtering' function. Must be True or False.
			- linewidth=1.0: linewidth.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert type(plot_band_filtering) is bool
		except AssertionError:
			print("Error at input 'plot_band_filtering': Must be True or False")
			return
		try:
			assert (type(linewidth) is int) or (type(linewidth) is float)
		except AssertionError:
			print("Error at input 'linewidth': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		try:
			assert self.computes_amplitude is True
		except AssertionError:
			print("Error: plot_amplitude cannot be applied")
			return
		if xaxis=="period":
			mypow=-1
		elif xaxis=="frequency":
			mypow=1
		if loglog is False:
			plt.plot(self.freq**mypow,self.amplitude,linewidth=linewidth)
		elif loglog is True:
			plt.loglog(self.freq**mypow,self.amplitude,linewidth=linewidth)
		if (plot_band_filtering is True and self.run_freq_filtering is True):
			y=np.zeros(2)
			y[0]=0.
			y[1]=np.amax(self.amplitude)
			for k in range(self.freq_filtered_signal_bounds.shape[0]):
				x1=np.zeros(2)
				x1[0]=self.freq_filtered_signal_bounds[k,0]
				x1[1]=self.freq_filtered_signal_bounds[k,0]
				x2=np.zeros(2)
				x2[0]=self.freq_filtered_signal_bounds[k,1]
				x2[1]=self.freq_filtered_signal_bounds[k,1]
				plt.fill_betweenx(y,x1,x2,edgecolors=None,facecolor='black',alpha=0.5)
		else:
			print("WARNING: function 'freq_filtering' was not run => unable to draw the bands")
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Amplitude"+self.mydata_label,fontsize=fontsize_axes)
		plt.title("Amplitude", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt
		
		
		
	def plot_amplitude_vs_periodogram(self,xaxis="frequency",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='xx-small',linewidth=1.0):
	
		""" plot_amplitude_vs_periodogram generates the figures of the squared amplitude periodogram and of the weighted periodogram. Both are usually very similar, as expected from the theory. Only available if the options computes_amplitude is True and weighted_WOSA is True in 'freq_analysis' function.
			Optional Inputs:
			- xaxis="frequency": choice for the horizontal figure axis: "frequency" or "period" (=1/frequency).
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='small': fontsize for the figure legend.
			- linewidth=1.0: linewidth.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth) is int) or (type(linewidth) is float)
		except AssertionError:
			print("Error at input 'linewidth': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_freq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'freq_analysis'")
			return
		try:
			assert self.computes_amplitude is True and self.weighted_WOSA is True
		except AssertionError:
			print("Error: plot_amplitude_vs_periodogram cannot be applied")
			return
		if xaxis=="period":
			mypow=-1
		elif xaxis=="frequency":
			mypow=1
		if loglog is False:
			plt.plot(self.freq**mypow,self.amplitude**2,label="Squared data amplitude",linewidth=linewidth)
			plt.plot(self.freq**mypow,self.periodogram,label="Data weighted periodogram",linewidth=linewidth)
		elif loglog is True:
			plt.loglog(self.freq**mypow,self.amplitude**2,label="Squared data amplitude",linewidth=linewidth)
			plt.loglog(self.freq**mypow,self.periodogram,label="Data weighted periodogram",linewidth=linewidth)
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		plt.suptitle("Squared Amplitude vs Weighted Periodogram",fontsize=fontsize_title)
		plt.legend(fancybox=True,fontsize=fontsize_legend,bbox_to_anchor=(1.1, 1.05))
		return plt



	def timefreq_analysis(self,theta=None,permin=None,permax=None,deltaj=0.05,w0=5.5,gauss_spread=3.0,eps=0.000001,dt_GCD=None,shannonnyquistexclusionzone=True,weighted_CWT=True,smoothing_coeff=0.5,smoothing_type="fixed",percentile=None,n_moments=2,MaxFunEvals=100000,algo_moments='gamma-polynomial',computes_amplitude=False,computes_global_scalogram=True,activate_perlim2=True):
		
		""" timefreq_analysis computes the (smoothed) scalogram and its confidence levels and the amplitude scalogram.
			Optional Inputs:
			- theta=None: the times at which the CWT is computed. Default is theta=np.linspace(t[0],t[-1],t.size), where t is the time vector of the data. N.B.: theta may be irregularly sampled, although we recommend a regular grid.
			- permin=None: minimal period. Default value is approximately permin=2.*dt_GCD (see below for the definition of dt_GCD).
			- permax=None: maximal period. Default value is approximately permax=np.pi*(t[-1]-t[0])/w0/gauss_spread, where t is the time vector of the data. w0 and gauss_spread are defined below.
			- deltaj=0.05: parameter controlling the density of periods/scales. Scale vector is defined as scale_min*2.0**(float(j)*deltaj) for j=0, 1, 2, ... A smaller deltaj implies denser scales and a better precision. We have the following approximate relation between periods and scales: period=2*np.pi*scale.
			- w0=5.5: the usual parameter for the Morlet wavelet controlling the time-frequency resolution. Minimal allowed value is w0=5.5. Increasing w0 offers a better scale/period resolution but a worse time resolution. There is always a trade-off due to the so-called 'Heisenberg uncertainty principle'.
			- gauss_spread=3.0: parameter for the spread of gaussian. 2*gauss_spread*std (where std is the standard dev. of the gaussian) is the approximate SUPPORT (i.e. where the function is not zero) of a gaussian. Typical values are gauss_spread=3.0 (conservative choice) or sqrt(2.0) (value taken in Torrence and Compo, 1998, and some subsequent papers). This is used for the computation of the cone of influence and for the max. allowed scale (Gaussian in fct of time), and for the refining of the SNEZ (Gaussian in fct of scale).
			- eps=1.0e-06: parameter controlling the flexibility on the value of the border of the SNEZ.
			- dt_GCD=None: the greatest common divisor of the time steps of the data time vector. Default value is the smallest time step (which is not exactly dt_GCD; it is an upper bound on the estimation of dt_GCD; but it is very simple to compute).
			- shannonnyquistexclusionzone=False: activate (True) or not (False) the Shannon-Nyquist exclusion zone (SNEZ).
			- weighted_CWT=False: True to get weighted scalogram, or False for the classical scalogram.
			- smoothing_coeff=0.5: smoothing the CWT is performed this way: at a given (theta,scale), the CWT_smoothed is the average of the CWT along neighbouring values of theta (at the same scale). The interval of those neighbouring values of theta is delimited by: theta-smoothing_coeff*std and theta+smoothing_coeff*std, where std is the standard deviation of the gaussian. Note that:
				-> the std depends on the scale.
				-> Near the edges of the time series or close the the SNEZ, the full interval over theta-smoothing_coeff*std to theta+smoothing_coeff*std cannot be considered. In such case, the CWT_smoothed at (theta,scale) is either ignored (if smoothing_type='fixed') or the interval around theta is shortened (if smoothing_type='variable').
			- smoothing_type='fixed': See above the explanations for 'smoothing_coeff'.
			- percentile=None: The x^th percentiles for the confidence levels. Must be a 1-dim numpy array. Default is the 95^th percentile (i.e. the 95% confidence level):percentile=np.zeros(1); percentile[0]=95.0
			- n_moments=2: number of conserved moments for the analytical confidence levels. Used if "a" is in signif_level_type (see function 'carma_params').
			- MaxFunEvals=100000: "max_nfev" option for "least_squares" - see python help of "scipy.optimize.least_squares". Used if algo_moments='generalized-gamma-polynomial'.
			- algo_moments='gamma-polynomial': Choice for the algorithm for the computation of analytical confidence levels. Used if "a" is in signif_level_type (see function 'carma_params'). algo_moments='gamma-polynomial' or algo='generalized-gamma-polynomial'.
			- computes_amplitude=False: computes the amplitude scalogram (if True) or not (if False).
			- computes_global_scalogram=True: computes the global scalogram (if True) or not (if False). It also computes the global amplitude (if computes_amplitude=True) and the global values of other time-frequency variables.
			- activate_perlim2=True: computes (and draws, for the functions making the figures) the refinement of the SNEZ. True or False.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import numpy.linalg as la
		from .dt_min import dt_min
		from tqdm import trange
		from .CWT import CWT
		from .CWT_and_Ampl import CWT_and_Ampl
		# from percentile_n_moments import percentile_n_moments
		# from white_noise_mcmc import white_noise_mcmc
		from .timefreq_analysis_prelims import timefreq_analysis_prelims
		import copy
		from scipy.optimize import brenth
		
		# check inputs
		try:
			assert (theta is None) or np.issubsctype(theta,float)
		except:
			print("Error at input 'theta': must be None or a numpy array of 'float' type")
			return
		try:
			assert (permin is None) or ((type(permin) is float) and (permin>0.))
		except AssertionError:
			print("Error at input 'permin': must be of 'float' type and >0.")
			return
		try:
			assert (permax is None) or ((type(permax) is float) and (permax>0.))
		except AssertionError:
			print("Error at input 'permax': must be of 'float' type and >0.")
			return
		if (permin is not None) and (permax is not None):
			try:
				assert permax>=permin
			except AssertionError:
				print("Error at input 'permin' and 'permax': must have permax>=permin")
				return
		try:
			assert (type(deltaj) is float) and (deltaj>0.)
		except AssertionError:
			print("Error at input 'deltaj': must be of 'float' type and >0.")
			return
		try:
			assert (type(w0) is float) and (w0>=5.5)
		except AssertionError:
			print("Error at input 'w0': must be of 'float' type and >=5.5")
			return
		try:
			assert (type(gauss_spread) is float) and (gauss_spread>0.)
		except AssertionError:
			print("Error at input 'gauss_spread': must be of 'float' type and >0.")
			return
		try:
			assert (type(eps) is float) and (eps>=0.)
		except AssertionError:
			print("Error at input 'eps': must be of 'float' type and >=0.")
			return
		try:
			assert (dt_GCD is None) or ((type(dt_GCD) is float) and (dt_GCD>0.))
		except AssertionError:
			print("Error at input 'dt_GCD': must be None or of 'float' type and >0.")
			return
		try:
			assert type(shannonnyquistexclusionzone) is bool
		except AssertionError:
			print("Error at input 'shannonnyquistexclusionzone': Must be True or False")
			return
		try:
			assert type(weighted_CWT) is bool
		except AssertionError:
			print("Error at input 'weighted_CWT': Must be True or False")
			return
		try:
			assert (type(smoothing_coeff) is float) and (smoothing_coeff>=0.)
		except AssertionError:
			print("Error at input 'smoothing_coeff': must be of 'float' type and >=0.")
			return
		try:
			assert (type(smoothing_type) is str) and ((smoothing_type.lower()=="fixed") or (smoothing_type.lower()=="variable"))
		except AssertionError:
			print("Error at input 'smoothing_type': Must be 'fixed' or 'variable'")
			return
		smoothing_type=smoothing_type.lower()
		try:
			assert (percentile is None) or np.issubsctype(percentile,float)
		except:
			print("Error at input 'percentile': must be None or a numpy array of 'float' type")
			return
		try:
			assert (type(n_moments) is int) and (n_moments>=2)
		except AssertionError:
			print("Error at input 'n_moments': must be of 'int' type and >=2")
			return
		try:
			assert (type(MaxFunEvals) is int) and (MaxFunEvals>=10)
		except AssertionError:
			print("Error at input 'MaxFunEvals': must be of 'int' type and >=10")
			return
		try:
			assert (type(algo_moments) is str) and ((algo_moments.lower()=="gamma-polynomial") or (algo_moments.lower()=="generalized-gamma-polynomial"))
		except AssertionError:
			print("Error at input 'algo_moments': must be 'gamma-polynomial' or 'generalized-gamma-polynomial'")
			return
		algo_moments=algo_moments.lower()
		try:
			assert type(computes_amplitude) is bool
		except AssertionError:
			print("Error at input 'computes_amplitude': must be True or False")
			return
		try:
			assert type(computes_global_scalogram) is bool
		except AssertionError:
			print("Error at input 'computes_global_scalogram': must be True or False")
			return
		try:
			assert type(activate_perlim2) is bool
		except AssertionError:
			print("Error at input 'activate_perlim2': must be True or False")
			return
		# check that some functions were previously run
		if ('a' in self.signif_level_type) or ('n' in self.signif_level_type):
			try:
				assert self.run_carma_params is True
			except AssertionError:
				print("Error: Must have run function 'carma_params'")
				return
		else:
			try:
				assert self.run_check_data is True
			except AssertionError:
				print("Error: Must have run function 'check_data'")
				return
			try:
				assert self.run_choose_trend_degree is True
			except AssertionError:
				print("Error: Must have run function 'choose_trend_degree'")
				return
		try:
			assert self.run_trend_vectors is True
		except AssertionError:
			print("Error: Must have run function 'trend_vectors'")
			return
		# Period -> Scale conversion for 'permin' and 'permax' if they are provided by the user
		if permin is None:
			scalemin=None
		else:
			if weighted_CWT is True:
				scalemin=permin/(2.*np.pi)
			elif weighted_CWT is False:
				scalemin=permin*(w0+np.sqrt(2+w0**2))/(4.*np.pi*w0)
		if permax is None:
			scalemax=None
		else:
			if weighted_CWT is True:
				scalemax=permax/(2.*np.pi)
			elif weighted_CWT is False:
				scalemax=permax*(w0+np.sqrt(2+w0**2))/(4.*np.pi*w0)
		# Set Default values for input arguments
		if dt_GCD is None:
			dt_GCD=dt_min(self.t)
		if scalemin is None:
			scalemin=dt_GCD/np.pi*(1.+eps)
		if scalemax is None:
			scalemax=(self.t[-1]-self.t[0])/2./w0/(gauss_spread+smoothing_coeff)
		if theta is None:
			self.theta=np.linspace(self.t[0],self.t[-1],self.t.size)
		else:
			self.theta=theta
		if percentile is None:
			percentile=np.zeros(1)
			percentile[0]=95.0
		self.percentile_cwt=percentile
		self.smoothing_coeff=smoothing_coeff
		self.n_moments_cwt=n_moments
		self.w0=w0
		# theta is put in ascending order and check whether theta range is included in self.t range
		theta_ind_sort=np.argsort(self.theta)
		self.theta=self.theta[theta_ind_sort]
		try:
			assert self.theta[0]>=self.t[0]
		except AssertionError:
			print("Error: theta[0] must be >= than the smallest time of the time series")
			return
		try:
			assert self.theta[-1]<=self.t[-1]
		except AssertionError:
			print("Error: theta[-1] must be <= than the biggest time of the time series")
			return
		# Adjust scalemin and scalemax and builds the scale vector
		scalemin=np.maximum(dt_GCD/np.pi*(1.+eps),scalemin)
		scalemax=np.minimum((self.t[-1]-self.t[0])/2./w0/(gauss_spread+smoothing_coeff),scalemax)
		J=int(np.floor(np.log2(scalemax/scalemin)/deltaj))
		scale=np.zeros(J+1)
		for k in range(J+1):
			scale[k]=scalemin*2.**(float(k)*deltaj)
		# Build the CWT components
		scale,self.coi1,self.coi2,self.coi1_smooth,self.coi2_smooth,coi_smooth_ind,weight_cwt,scalelim1_ind,scalelim1_smooth,scalelim1_ind_smooth,Qmax,n_outside_scalelim1=timefreq_analysis_prelims(self.t,self.theta,scale,w0,gauss_spread,eps,dt_GCD,shannonnyquistexclusionzone,weighted_CWT,smoothing_coeff,smoothing_type)
		self.weighted_CWT=weighted_CWT
		J=scale.size
		Q=self.theta.size
		self.shannonnyquistexclusionzone=shannonnyquistexclusionzone
		# Scale -> Period conversion
		if weighted_CWT is True:
			self.period_cwt=scale*2.*np.pi
			self.perlim1_smooth_cwt=scalelim1_smooth*2.*np.pi
		elif weighted_CWT is False:
			self.period_cwt=scale*4.*np.pi*w0/(w0+np.sqrt(2+w0**2))
			self.perlim1_smooth_cwt=scalelim1_smooth*4.*np.pi*w0/(w0+np.sqrt(2+w0**2))
		if computes_amplitude is True:
			self.period_ampl=scale*2.*np.pi
			self.perlim1_smooth_ampl=scalelim1_smooth*2.*np.pi
		print("Re-estimated period range: from ", self.period_cwt[0]," to ", self.period_cwt[-1])
		# CWT: significance levels, data scalogram, amplitude - This is the main loop, over the scales
		npercentile=percentile.size
		self.scalogram=np.zeros((Q,J))
		if computes_amplitude is True:
			self.cwtamplitude=np.zeros((Q,J))
			if smoothing_coeff==0.:
				self.cwtamplitude_cos=np.zeros((Q,J))
				self.cwtamplitude_sin=np.zeros((Q,J))
				self.n_outside_scalelim1=n_outside_scalelim1
				self.weight_cwt=weight_cwt
		self.computes_cwtamplitude=computes_amplitude
		self.computes_global_scalogram=computes_global_scalogram
		if computes_global_scalogram is True:
			self.global_scalogram=np.zeros(J)
			if computes_amplitude is True:
				self.global_amplitude=np.zeros(J)
		mydata_transp=np.transpose(self.mydata)
		self.r=np.zeros((Q,J),dtype=int)
		print("Main loop, over the time-frequency plane:")
		for l in trange(J):
			scale_l=scale[l]
			if computes_amplitude is True:
				Amplitude, M2, mytheta, mytheta_ind=CWT_and_Ampl(self.t,self.mydata,self.theta,l,scale_l,self.myprojvec,self.Vmat,self.pol_degree,weight_cwt[:,l],scalelim1_ind,n_outside_scalelim1[l],w0)
			else:
				M2, mytheta, mytheta_ind=CWT(self.t,self.theta,l,scale_l,self.myprojvec,self.pol_degree,weight_cwt[:,l],scalelim1_ind,n_outside_scalelim1[l],w0)
			# Data scalogram - intermediate calculus
			scalogram_int=np.dot(mydata_transp,M2)
			#  Store the amplitude relative to cos/sin if no smoothing (to be used for band and ridges filtering)
			if computes_amplitude is True:
				if smoothing_coeff==0.:
					self.cwtamplitude_cos[mytheta_ind,l]=copy.copy(Amplitude[0,:])
					self.cwtamplitude_sin[mytheta_ind,l]=copy.copy(Amplitude[1,:])
			count=0
			if computes_global_scalogram is True:
				if self.p==0:
					mymat_full=np.dot(np.transpose(M2),M2)
				weight_gs=np.zeros(n_outside_scalelim1[l])
				for k in range(n_outside_scalelim1[l]):
					if (smoothing_type=="fixed" and (l<scalelim1_ind_smooth[mytheta_ind[k]] or l>coi_smooth_ind[mytheta_ind[k]])):
						continue
					else:
						count+=1
						mytheta_k=mytheta[k]
						ind_left=np.argmin(np.absolute(mytheta-(mytheta_k-smoothing_coeff*w0*scale_l)))
						ind_right=np.argmin(np.absolute(mytheta-(mytheta_k+smoothing_coeff*w0*scale_l)))
						mylength=ind_right+1-ind_left
						weight_gs[ind_left:(ind_right+1)]+=1./float(mylength)
						# Smoothing factor
						self.r[mytheta_ind[k],l]=mylength
						# Data scalogram
						self.scalogram[mytheta_ind[k],l]=la.norm(scalogram_int[2*ind_left:2*(ind_right+1)])**2/float(mylength)
						# Data Amplitude
						if computes_amplitude is True:
							self.cwtamplitude[mytheta_ind[k],l]=np.sqrt(np.sum(la.norm(Amplitude[:,ind_left:(ind_right+1)],axis=0)**2)/float(mylength))
						# Old stuff for next loop:
						ind_right_old=ind_right
						ind_left_old=ind_left
				if count>0:
					weight_gs/=float(count)
					# Global Scalogram
					self.global_scalogram[l]=np.sum(self.scalogram[:,l])/float(count)
					# Global Amplitude
					if computes_amplitude is True:
						self.global_amplitude[l]=np.sqrt(np.sum(self.cwtamplitude[:,l]**2)/float(count))

			else:
				mymat_full=np.dot(np.transpose(M2),np.dot(self.ARMA_mat_unique,M2))
				for k in range(n_outside_scalelim1[l]):
					if (smoothing_type=="fixed" and (l<scalelim1_ind_smooth[mytheta_ind[k]] or l>coi_smooth_ind[mytheta_ind[k]])):
						continue
					else:
						count+=1
						mytheta_k=mytheta[k]
						ind_left=np.argmin(np.absolute(mytheta-(mytheta_k-smoothing_coeff*w0*scale_l)))
						ind_right=np.argmin(np.absolute(mytheta-(mytheta_k+smoothing_coeff*w0*scale_l)))
						mylength=ind_right+1-ind_left
						# Smoothing factor
						self.r[mytheta_ind[k],l]=mylength
						# Data scalogram
						self.scalogram[mytheta_ind[k],l]=la.norm(scalogram_int[2*ind_left:2*(ind_right+1)])**2/float(mylength)
						# Data Amplitude
						if computes_amplitude is True:
							self.cwtamplitude[mytheta_ind[k],l]=np.sqrt(np.sum(la.norm(Amplitude[:,ind_left:(ind_right+1)],axis=0)**2)/float(mylength))
						# Old stuff for next loop:
						ind_right_old=ind_right
						ind_left_old=ind_left
		coi_ind=np.zeros(Q,dtype=int)
		for k in range(Q):
			theta_k=self.theta[k]
			for l in range(J-1,scalelim1_ind_smooth[k]-1,-1):
				if theta_k>max(self.coi1[l],self.coi1_smooth[l]) and theta_k<min(self.coi2[l],self.coi2_smooth[l]):
					coi_ind[k]=l
					break
		if shannonnyquistexclusionzone is False:
			if activate_perlim2 is True:
				print("shannonnyquistexclusionzone is False => activate_perlim2 is changed to False")
				activate_perlim2=False
		if activate_perlim2 is True:
			scalelim2_smooth_scal=copy.copy(scalelim1_smooth)
			# computes the spread in scale at scalelim1_smooth
			if weighted_CWT is True:
				mystd_scalelim1_smooth_scal=gauss_spread*self.perlim1_smooth_cwt/np.sqrt(2.)/2./np.pi/w0
				scalelim2_smooth_scal=scalelim1_smooth+mystd_scalelim1_smooth_scal
			elif weighted_CWT is False:
				mystd_scalelim1_smooth_scal=np.zeros(Q)
				for k in range(Q):
					scale_max_k=self.perlim1_smooth_cwt[k]*(w0+np.sqrt(w0**2+2.))/4./np.pi/w0
					myfun=lambda xx: xx*np.exp(-((2.*np.pi/self.perlim1_smooth_cwt[k]*xx-1.)**2)*w0**2)-scale_max_k*np.exp(-((2.*np.pi/self.perlim1_smooth_cwt[k]*scale_max_k-1.)**2)*w0**2)*np.exp(-gauss_spread**2/2.)
					scalelim2_smooth_scal[k], rootresults=brenth(myfun,scalelim1_smooth[k],scale[-1],maxiter=1000,full_output=True,disp=False) # brenth is a scipy root finding algorithm
					if rootresults.converged==False:
						scalelim2_smooth_scal[k]=scale[-1]
			if computes_amplitude is True:
				mystd_scalelim1_smooth_ampl=gauss_spread*self.perlim1_smooth_ampl/2./np.pi/w0
				scalelim2_smooth_ampl=scalelim1_smooth+mystd_scalelim1_smooth_ampl
		else:
			my_scalelim1_smooth=np.zeros(Q)
			for k in range(Q):
				my_scalelim1_smooth[k]=scale[max(0,scalelim1_ind_smooth[k]-1)]
			scalelim2_smooth_scal=copy.copy(my_scalelim1_smooth)
			scalelim2_smooth_ampl=copy.copy(my_scalelim1_smooth)
		# computes the min/max for the color scale
		minscal=float("inf")
		maxscal=-float("inf")
		minampl=float("inf")
		maxampl=-float("inf")
		minampl_sq=float("inf")
		maxampl_sq=-float("inf")
		min_pseudo_cwtspectrum_anal=float("inf")
		max_pseudo_cwtspectrum_anal=-float("inf")
		min_pseudo_cwtspectrum_mcmc=float("inf")
		max_pseudo_cwtspectrum_mcmc=-float("inf")
		min_cwt_variance_anal=float("inf")
		max_cwt_variance_anal=-float("inf")
		for k in range(Q):
			scalelim2_smooth_scal_ind=np.argmin(np.absolute(scalelim2_smooth_scal[k]-scale))+1
			if (coi_ind[k]>scalelim2_smooth_scal_ind):
				minscal=min(minscal,np.amin(self.scalogram[k,scalelim2_smooth_scal_ind:coi_ind[k]+1]))
				maxscal=max(maxscal,np.amax(self.scalogram[k,scalelim2_smooth_scal_ind:coi_ind[k]+1]))
		self.minscal=minscal
		self.maxscal=maxscal
		
		if computes_amplitude is True:
			for k in range(Q):
				scalelim2_smooth_ampl_ind=np.argmin(np.absolute(scalelim2_smooth_ampl[k]-scale))+1
				scalelim2_smooth_scal_ind=np.argmin(np.absolute(scalelim2_smooth_scal[k]-scale))+1
				if (coi_ind[k]>scalelim2_smooth_ampl_ind):
					minampl=min(minampl,np.amin(self.cwtamplitude[k,scalelim2_smooth_ampl_ind:coi_ind[k]+1]))
					maxampl=max(maxampl,np.amax(self.cwtamplitude[k,scalelim2_smooth_ampl_ind:coi_ind[k]+1]))
					minampl_sq=min(minampl_sq,np.amin(self.cwtamplitude[k,scalelim2_smooth_scal_ind:coi_ind[k]+1]**2))
					maxampl_sq=max(maxampl_sq,np.amax(self.cwtamplitude[k,scalelim2_smooth_scal_ind:coi_ind[k]+1]**2))
			self.minampl=minampl
			self.maxampl=maxampl
			self.minampl_sq=minampl_sq
			self.maxampl_sq=maxampl_sq
		# Scale -> Period conversion
		if weighted_CWT is True:
			self.perlim2_smooth_scal=scalelim2_smooth_scal*2.*np.pi
		elif weighted_CWT is False:
			self.perlim2_smooth_scal=scalelim2_smooth_scal*4.*np.pi*w0/(w0+np.sqrt(2+w0**2))
		if computes_amplitude is True:
			self.perlim2_smooth_ampl=scalelim2_smooth_ampl*2.*np.pi

		self.run_timefreq_analysis=True



	def timefreq_ridges_filtering(self,N=1.5,maxampl_point_coeff=0.5,maxampl_ridge_coeff=0.8,chaining_ridges_coeff=0.5):
		
		""" timefreq_ridges_filtering computes the ridges of the wavelet transform of the data. The amplitude scalogram must have been computed and without smoothing.
			Optional Inputs:
			- N=1.5: Removes all ridges of less than N periods in length
			- maxampl_point_coeff=0.5: Removes all small amplitude ridge points having ampl<maxampl_point_coeff*maxampl (where maxampl is the max of the amplitude over the time-frequency plane).
			- maxampl_ridge_coeff=0.8: Selects all the ridges for which at least one ridge point is >maxampl_ridge_coeff*maxampl (where maxampl is the max of the amplitude over the time-frequency plane).
			- chaining_ridges_coeff=0.5: Controls agressiveness of chaining ridge points across scales. Increase the value to chain ridges more agressively across scales, or decrease it to supress chaining across scales.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import sys
		import ridges as ridges
		import copy
		
		# check inputs
		try:
			assert (type(N) is float) and (N>0.)
		except AssertionError:
			print("Error at input 'N': must be of 'float' type and >0.")
			return
		try:
			assert (type(maxampl_point_coeff) is float) and (maxampl_point_coeff>0.)
		except AssertionError:
			print("Error at input 'maxampl_point_coeff': must be of 'float' type and >0.")
			return
		try:
			assert (type(maxampl_ridge_coeff) is float) and (maxampl_ridge_coeff>=0.) and (maxampl_ridge_coeff<=1.)
		except AssertionError:
			print("Error at input 'maxampl_ridge_coeff': must be of 'float' type and between 0. and 1.")
			return
		try:
			assert (type(chaining_ridges_coeff) is float) and (chaining_ridges_coeff>=0.)
		except AssertionError:
			print("Error at input 'chaining_ridges_coeff': must be of 'float' type and >=0.")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True
		except AssertionError:
			print("Error: In function 'timefreq_analysis', amplitude must have been computed")
			return
		try:
			assert self.smoothing_coeff==0.
		except AssertionError:
			print("Error: In function 'timefreq_analysis', variable 'smoothing_coeff' must be equal to zero")
			return

		omega=2.*np.pi/self.period_ampl
		IR,JR,XR,_=ridges.ridgewalk(self.cwtamplitude,omega,N,maxampl_point_coeff*self.maxampl,chaining_ridges_coeff)

		if IR.size>0:
			NR=0
			ENDR=np.zeros(IR.size,dtype=int)
			STARTR=np.zeros(IR.size,dtype=int)
			STARTR[0]=0
			for k in range(IR.size):
				if IR[k]==-999:
					NR+=1
					ENDR[NR-1]=k-1
					STARTR[NR]=k+1
			skeleton=NR*[None]     # Store ridges informations
			for k in range(NR):
				skeleton[k]=5*[None]
				skeleton[k][0]=self.theta[IR[STARTR[k]:ENDR[k]+1]]
				skeleton[k][1]=self.period_ampl[JR[STARTR[k]:ENDR[k]+1]]
				length_skel=len(range(STARTR[k],ENDR[k]+1))
				skeleton[k][2]=np.zeros(length_skel)
				skeleton[k][3]=np.zeros(length_skel)
				skeleton[k][4]=np.zeros(length_skel)
				for l in range(length_skel):
					skeleton[k][2][l]=self.cwtamplitude_cos[IR[STARTR[k]+l],JR[STARTR[k]+l]]
					skeleton[k][3][l]=self.cwtamplitude_sin[IR[STARTR[k]+l],JR[STARTR[k]+l]]
					skeleton[k][4][l]=self.cwtamplitude[IR[STARTR[k]+l],JR[STARTR[k]+l]]
			maxampl_ridge=maxampl_ridge_coeff*self.maxampl
			signifskel=np.zeros(0,dtype=int)
			for k in range(NR):
				mybool=skeleton[k][4]>maxampl_ridge				#<maxampl_ridge;
				if np.sum(mybool)>0:							#np.prod(bool)==0
					signifskel=np.append(signifskel,k)
			nsk=signifskel.size
			print("Number of ridges is "+str(nsk))
			skeletonbis=copy.copy(skeleton)
			self.skeleton=[]
			for k in range(nsk):
				self.skeleton.append(skeletonbis[signifskel[k]])
			del skeletonbis
		else:
			nsk=0
			self.skeleton=0
		
		self.run_timefreq_ridges_filtering=True



	def timefreq_band_filtering(self,period_bounds):
		
		""" timefreq_band_filtering filters the data in a given period band. The amplitude scalogram must have been computed and without smoothing.
			Required Inputs:
			- period_bounds: a list of tuples. Each tuple contains the two periods delimitating the band. Example: [(19.,23.), (37.,43.), (95.,105.)] will filter in three bands. First one: between periods 19. and 23., second one: between periods 37. and 43. and third one: between periods 95. and 105.
			Outputs:
			/
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
			
		# check inputs
		try:
			assert type(period_bounds) is list
		except AssertionError:
			print("Error at input 'period_bounds': must be of type 'list'")
			return
		for k in range(len(period_bounds)):
			try:
				assert type(period_bounds[k]) is tuple
			except AssertionError:
				print("Error at input 'period_bounds': must be a list containing tuples")
				return
			try:
				assert (type(period_bounds[k][0]) is float) and (type(period_bounds[k][1]) is float)
			except:
				print("Error at input 'period_bounds': must be a list containing tuples, with each tuple containing 2 floats")
				return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True
		except AssertionError:
			print("Error: In function 'timefreq_analysis', amplitude must have been computed")
			return
		try:
			assert self.smoothing_coeff==0.
		except AssertionError:
			print("Error: In function 'timefreq_analysis', variable 'smoothing_coeff' must be equal to zero")
			return
		
		period_min=self.period_ampl[0]
		deltaj=np.log2(self.period_ampl[1]/self.period_ampl[0])
		self.timefreq_band_filtered_signal=np.zeros((self.theta.size,len(period_bounds)))
		self.timefreq_band_filtered_signal_bounds=np.zeros((len(period_bounds),2))
		for l in range(len(period_bounds)):
			ind_low=int(np.floor(np.log2(period_bounds[l][0]/period_min)/deltaj))
			ind_high=int(np.ceil(np.log2(period_bounds[l][1]/period_min)/deltaj))
			try:
				assert ind_high>=ind_low
			except AssertionError:
				print("Error in 'period_bounds' input parameters: must have period_high>=period_low")
				return
			count=np.zeros(self.theta.size,dtype=int)
			for k in range(ind_low,(ind_high+1)):
				self.timefreq_band_filtered_signal[:,l]+=self.cwtamplitude_cos[:,k]*np.cos(2.*np.pi/self.period_ampl[k]*self.theta)+self.cwtamplitude_sin[:,k]*np.sin(2.*np.pi/self.period_ampl[k]*self.theta)
				for m in range(self.theta.size):
					if self.weight_cwt[m,k]>-1.:
						count[m]+=1
			count[count==0]=1
			for k in range(self.theta.size):
				self.timefreq_band_filtered_signal[k,l]/=float(count[k])
			self.timefreq_band_filtered_signal_bounds[l,0]=self.period_ampl[ind_low]
			self.timefreq_band_filtered_signal_bounds[l,1]=self.period_ampl[ind_high]

		self.run_timefreq_band_filtering=True



	def plot_scalogram(self,with_global_scalogram=True,time_string=None,period_string=None,power_string=None,dashed_periods=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,left_padding=0.,right_padding=0.,middle_padding=0.,color_cl_anal=None,color_cl_mcmc=None,linewidth_cl=2,global_scal_xlabel="top",global_scal_xlabel_ticks="top",linewidth_gscal=1.0,linewidth_gcl=1.0,cmap="jet",nlevels=50,plot_coi="fill",linewidth_coi=1.0,plot_perlim2="fill",linewidth_perlim2=1.0,reverse_xaxis=False,reverse_yaxis=False,alpha_SNEZ=0.5,decimals=3,minscal=None,maxscal=None,permin_ext=None):
		
		""" plot_scalogram generates the figure of the scalogram and its confidence levels. It also generates the figure of the global scalogram and its confidence levels.
			Optional Inputs:
			- with_global_scalogram=True: with the global scalogram, or False.
			- time_string=None: list of floats containing the location of the ticks for the time axis.
			- period_string=None: list of floats containing the location of the ticks for the period axis.
			- power_string=None: list of floats containing the location of the ticks for the power axis of the global scalogram.
			- dashed_periods=None: list of floats containing the periods for which a dashed line is drawn.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- left_padding=0.: padding on the left of the scalogram figure.
			- right_padding=0.: padding on the right of the global scalogram.
			- middle_padding=0.: padding between the scalogram and the global scalogram.
			N.B.: paddings are only active if with_global_scalogram is True.
			WARNING with the paddings: If there is an overlap between figures, python can choose to not display one.
			- color_cl_anal=None: colors for the analytical confidence levels/contours. Ex.: ['m','g','r'] will draw in magenta, green and red the 3 levels specified in 'percentile' input of function 'timefreq_analysis'. Default is a magenta contour for all the confidence levels.
			- color_cl_mcmc=None: colors for the mcmc confidence levels/contours. Ex.: ['m','g','r'] will draw in magenta, green and red the 3 levels specified in 'percentile' input of function 'timefreq_analysis'. Default is a green contour for all the confidence levels.
			- linewidth_cl=2: linewidth of the confidence contours in the scalogram
			- global_scal_xlabel="top": location of the xlabel for the global scalogram: "top" or "bottom".
			- global_scal_xlabel_ticks="top": location of the ticks of the xlabel for the global scalogram: "top" or "bottom".
			- linewidth_gscal=1.0: linewidth for the global scalogram.
			- linewidth_gcl=1.0: linewidth for the confidence levels of the global scalogram.
			- cmap="jet": colormap for the scalogram. Other choices on http://matplotlib.org/users/colormaps.html
			- nlevels=50: number of automatically-chosen color levels. 
			- plot_coi="fill": plotting-type for the cone of influence: "fill" or "line".
			- linewidth_coi=1.0: linewidth for the coi (if plot_coi="line").
			- plot_perlim2="fill": plotting-type for the refinement of the SNEZ: "fill" or "line".
			- linewidth_perlim2=1.0: linewidth for perlim2 (if plot_perlim2="line").
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			- alpha_SNEZ=0.5: Transparency for the SNEZ. It must take a value between 0 (completely transparent) and 1 (completely opaque). Only used if shannonnyquistexclusionzone=False in the method 'timefreq_analysis'.
			- decimals=3: Numbers of decimals for the colorbar scale ticks under the scientific notation (e.g. 1.567e-10).
			- minscal=None: Fix the lower bound of the color scale. Default (None) takes the min value of the scalogram over the non-shaded area.
			- maxscal=None: Fix the higher bound of the color scale. Default (None) takes the max value of the scalogram over the non-shaded area.
			- permin_ext=None: If a float value is given, the SNEZ is broaden towards this given period. Available only if shannonnyquistexclusionzone=True in the method 'timefreq_analysis'.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		from copy import copy
		
		# check inputs
		try:
			assert type(with_global_scalogram) is bool
		except AssertionError:
			print("Error at input 'with_global_scalogram': must be True or False")
			return
		try:
			assert (time_string is None) or (type(time_string) is list)
		except AssertionError:
			print("Error at input 'time_string': must be None or of type 'list'")
			return
		if type(time_string) is list:
			for k in range(len(time_string)):
				try:
					assert type(time_string[k]) is float
				except AssertionError:
					print("Error at input 'time_string': must be a list containing floats")
					return
		try:
			assert (period_string is None) or (type(period_string) is list)
		except AssertionError:
			print("Error at input 'period_string': must be None or of type 'list'")
			return
		if type(period_string) is list:
			for k in range(len(period_string)):
				try:
					assert type(period_string[k]) is float
				except AssertionError:
					print("Error at input 'period_string': must be a list containing floats")
					return
		try:
			assert (power_string is None) or (type(power_string) is list)
		except AssertionError:
			print("Error at input 'power_string': must be None or of type 'list'")
			return
		if type(power_string) is list:
			for k in range(len(power_string)):
				try:
					assert type(power_string[k]) is float
				except AssertionError:
					print("Error at input 'power_string': must be a list containing floats")
					return
		try:
			assert (dashed_periods is None) or (type(dashed_periods) is list)
		except AssertionError:
			print("Error at input 'dashed_periods': must be None or of type 'list'")
			return
		if type(dashed_periods) is list:
			for k in range(len(dashed_periods)):
				try:
					assert type(dashed_periods[k]) is float
				except AssertionError:
					print("Error at input 'dashed_periods': must be a list containing floats")
					return
		try:
			assert type(left_padding) is float
		except AssertionError:
			print("Error at input 'left_padding': must be of type 'float'")
			return
		try:
			assert type(right_padding) is float
		except AssertionError:
			print("Error at input 'right_padding': must be of type 'float'")
			return
		try:
			assert type(middle_padding) is float
		except AssertionError:
			print("Error at input 'middle_padding': must be of type 'float'")
			return
		try:
			assert (color_cl_anal is None) or (type(color_cl_anal) is list)
		except AssertionError:
			print("Error at input 'color_cl_anal': must be None or of type 'list'")
			return
		try:
			assert (color_cl_mcmc is None) or (type(color_cl_mcmc) is list)
		except AssertionError:
			print("Error at input 'color_cl_mcmc': must be None or of type 'list'")
			return
		try:
			assert (type(linewidth_cl) is int) or (type(linewidth_cl) is float)
		except AssertionError:
			print("Error at input 'linewidth_cl': must be an integer or float")
			return
		try:
			assert (type(global_scal_xlabel) is str) and ((global_scal_xlabel.lower()=="bottom") or (global_scal_xlabel.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_scal_xlabel': must be 'bottom' or 'top'")
			return
		global_scal_xlabel=global_scal_xlabel.lower()
		try:
			assert (type(global_scal_xlabel_ticks) is str) and ((global_scal_xlabel_ticks.lower()=="bottom") or (global_scal_xlabel_ticks.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_scal_xlabel_ticks': must be 'bottom' or 'top'")
			return
		global_scal_xlabel_ticks=global_scal_xlabel_ticks.lower()
		try:
			assert (type(linewidth_gscal) is int) or (type(linewidth_gscal) is float)
		except AssertionError:
			print("Error at input 'linewidth_gscal': must be an integer or float")
			return
		try:
			assert (type(linewidth_gcl) is int) or (type(linewidth_gcl) is float)
		except AssertionError:
			print("Error at input 'linewidth_gcl': must be an integer or float")
			return
		try:
			assert type(cmap) is str
		except AssertionError:
			print("Error at input 'cmap': must be of type 'str'")
			return
		try:
			assert (type(nlevels) is int) and nlevels>0
		except AssertionError:
			print("Error at input 'nlevels': must be an integer >0")
			return
		try:
			assert (type(plot_coi) is str) and ((plot_coi.lower()=="fill") or (plot_coi.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_coi': must be 'fill' or 'line'")
			return
		plot_coi=plot_coi.lower()
		try:
			assert (type(linewidth_coi) is int) or (type(linewidth_coi) is float)
		except AssertionError:
			print("Error at input 'linewidth_coi': must be an integer or float")
			return
		try:
			assert (type(plot_perlim2) is str) and ((plot_perlim2.lower()=="fill") or (plot_perlim2.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_perlim2': must be 'fill' or 'line'")
			return
		plot_perlim2=plot_perlim2.lower()
		try:
			assert (type(linewidth_perlim2) is int) or (type(linewidth_perlim2) is float)
		except AssertionError:
			print("Error at input 'linewidth_perlim2': must be an integer or float")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		try:
			assert ((type(alpha_SNEZ) is float) or (type(alpha_SNEZ) is int)) and (alpha_SNEZ>=0 and alpha_SNEZ<=1)
		except AssertionError:
			print("Error at input 'alpha_SNEZ': must of type float or int and must take a value in [0,1]")
			return
		try:
			assert (type(decimals) is int) and decimals>=0
		except AssertionError:
			print("Error at input 'decimals': must be an integer >=0")
			return
		try:
			assert (minscal is None) or (type(minscal) is float and minscal>=0.)
		except AssertionError:
			print("Error at input 'minscal': must be a float >=0")
			return
		try:
			assert (maxscal is None) or (type(maxscal) is float and maxscal>=0.)
		except AssertionError:
			print("Error at input 'maxscal': must be a float, >=0, and >= 'minscal'")
			return
		try:
			assert (permin_ext is None) or (type(permin_ext) is float and permin_ext>=0.)
		except AssertionError:
			print("Error at input 'permin_ext': must be a float >=0")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		if with_global_scalogram is True:
			try:
				assert self.computes_global_scalogram is True
			except AssertionError:
				print("Error: Global scalogram was not computed")
				print("=> drawing the figure without the global scalogram")
				with_global_scalogram=False
		if (self.shannonnyquistexclusionzone is False) and (permin_ext is not None):
			print("parameter 'permin_ext' changed to 'None', since 'shannonnyquistexclusionzone' is 'False'")
			permin_ext=None
		if period_string is None:
			if (permin_ext is not None) and (permin_ext<self.period_cwt[0]):
				myperiod=np.ceil(np.log2(permin_ext))
			else:
				myperiod=np.ceil(np.log2(self.period_cwt[0]))
			myperiod=2.**myperiod
			period_string=[myperiod]
			while True:
				myperiod*=2.
				if myperiod>self.period_cwt[-1]:
					break
				period_string.append(myperiod)
		if minscal is None:
			minscal=self.minscal
		if maxscal is None:
			maxscal=self.maxscal
		try:
			assert maxscal>=minscal
		except AssertionError:
			print("Error: variable 'maxscal' must be greater than or equal to 'minscal'")
			return
		if with_global_scalogram is True:
			gs1 = gridspec.GridSpec(1,5)
			gs1.update(left=0.05+left_padding)
			plt.subplot(gs1[0,:-1])
		if color_cl_anal is None:
			color_cl_anal=['m']*self.percentile_cwt.size
		if color_cl_mcmc is None:
			color_cl_mcmc=['g']*self.percentile_cwt.size
		minminscal=np.amin(np.amin(self.scalogram))
		maxmaxscal=np.amax(np.amax(self.scalogram))
		myscal=copy(self.scalogram)
		if maxscal<maxmaxscal:
			myscal[myscal>maxscal]=maxscal
			myvmax=maxscal
		else:
			myvmax=maxscal
		if minscal>minminscal:
			myscal[myscal<minscal]=minscal
			myvmin=minscal
		else:
			myvmin=minscal
		mycontourf=plt.contourf(self.theta,np.log2(self.period_cwt),np.transpose(myscal),nlevels,vmin=myvmin,vmax=myvmax,cmap=cmap)
		if plot_coi=="fill":
			plt.fill_betweenx(np.log2(self.period_cwt),self.t[0],self.coi1,edgecolors=None,facecolor='black',alpha=0.5)
			plt.fill_betweenx(np.log2(self.period_cwt),self.coi2,self.t[-1],edgecolors=None,facecolor='black',alpha=0.5)
		elif plot_coi=="line":
			plt.plot(self.coi1,np.log2(self.period_cwt),'k')
			plt.plot(self.coi2,np.log2(self.period_cwt),'k')
		plt.fill_betweenx(np.log2(self.period_cwt),self.theta[0],self.coi1_smooth,edgecolors=None,facecolor='black')
		plt.fill_betweenx(np.log2(self.period_cwt),self.coi2_smooth,self.theta[-1],edgecolors=None,facecolor='black')
		if self.shannonnyquistexclusionzone is True:
			if (permin_ext is not None) and (permin_ext<self.period_cwt[0]):
				plt.fill_between(self.theta,np.log2(permin_ext)*np.ones(self.theta.size),np.log2(self.perlim1_smooth_cwt),edgecolors=None,facecolor='black')
			else:
				plt.fill_between(self.theta,np.log2(self.period_cwt[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_cwt),edgecolors=None,facecolor='black')
			if plot_perlim2=="fill":
				plt.fill_between(self.theta,np.log2(self.period_cwt[0])*np.ones(self.theta.size),np.log2(self.perlim2_smooth_scal),edgecolors=None,facecolor='black',alpha=0.5)
			elif plot_perlim2=="line":
				plt.plot(self.theta,np.log2(self.perlim2_smooth_scal),'k')
		elif self.shannonnyquistexclusionzone is False:
			plt.fill_between(self.theta,np.log2(self.period_cwt[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_cwt),edgecolors=None,facecolor='black',alpha=alpha_SNEZ)
		ax = plt.gca()
		ax.tick_params(length=5, width=1, color='w')
		if dashed_periods is not None:
			for k in range(len(dashed_periods)):
				plt.plot(self.theta,np.log2(dashed_periods[k])*np.ones(self.theta.size),'w--')
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel("Period"+self.t_label,fontsize=fontsize_axes)
		mytitle="Wavelet Scalogram"
		if time_string is not None:
			plt.xticks(time_string, time_string)
		plt.xticks(fontsize=fontsize_ticks)
		plt.yticks(np.log2(period_string), period_string, fontsize=fontsize_ticks)
		plt.xlim([self.theta[0], self.theta[-1]])
		if (permin_ext is not None) and (permin_ext<self.period_cwt[0]):
			plt.ylim([np.log2(permin_ext), np.log2(self.period_cwt[-1])])
		else:
			plt.ylim([np.log2(self.period_cwt[0]), np.log2(self.period_cwt[-1])])
		plt.suptitle(mytitle, fontsize=fontsize_title)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		if with_global_scalogram is True:
			gs2 = gridspec.GridSpec(1,5)
			gs2.update(left=-0.3+middle_padding,right=0.95-right_padding)
			plt.subplot(gs2[0,-1])
			glob_scal_min=np.amin(self.global_scalogram)
			glob_scal_max=np.amax(self.global_scalogram)
			plt.plot(self.global_scalogram,np.log2(self.period_cwt),"b",linewidth=linewidth_gscal)
			myrange=np.arange(glob_scal_min,glob_scal_max,(glob_scal_max-glob_scal_min)/1000.)
			if dashed_periods is not None:
				for k in range(len(dashed_periods)):
					plt.plot(myrange,np.log2(dashed_periods[k])*np.ones(len(myrange)),'k--')
			plt.yticks(np.log2(period_string),[])
			if (permin_ext is not None) and (permin_ext<self.period_cwt[0]):
				plt.ylim([np.log2(permin_ext), np.log2(self.period_cwt[-1])])
			else:
				plt.ylim([np.log2(self.period_cwt[0]), np.log2(self.period_cwt[-1])])
			plt.xlim([glob_scal_min, glob_scal_max])
			plt.xlabel("Power"+self.power_label,fontsize=fontsize_axes)
			ax=plt.gca()
			ax.xaxis.set_label_position(global_scal_xlabel)
			if power_string is not None:
				plt.xticks(power_string, power_string)
			if global_scal_xlabel_ticks=="top":
				plt.tick_params(axis='x', labelbottom='off', labeltop='on', labelsize=fontsize_ticks)
			elif global_scal_xlabel_ticks=="bottom":
				plt.tick_params(axis='x', labelbottom='on', labeltop='off', labelsize=fontsize_ticks)
			if reverse_yaxis is True:
				plt.gca().invert_yaxis()
		# Colorbar and its rescaling: min and max of the levels of color are defined over the non-shaded regions
		cbar=plt.colorbar(mycontourf)
		cbar.ax.set_ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		cbar.ax.tick_params(labelsize=fontsize_ticks)
		my_color_ticks=np.linspace(max(minscal,minminscal),min(maxscal,maxmaxscal),10)
		cbar.set_ticks(my_color_ticks)
		my_color_ticklabels=copy(my_color_ticks)
		my_power_color_ticklabels=float(10**(int(np.floor(np.log10(my_color_ticklabels[-1])))))
		my_color_ticklabels=[el/my_power_color_ticklabels for el in my_color_ticklabels]
		my_color_ticklabels=np.around(my_color_ticklabels,decimals=decimals)
		mystring=' '.join(str(e) for e in my_color_ticklabels).split(' ')
		if minscal>minminscal and (('1' in mystring[0]) or ('2' in mystring[0]) or ('3' in mystring[0]) or ('4' in mystring[0]) or ('5' in mystring[0]) or ('6' in mystring[0]) or ('7' in mystring[0]) or ('8' in mystring[0]) or ('9' in mystring[0])):
			mystring[0]="<="+mystring[0]
		if maxscal<maxmaxscal:
			mystring[-1]=">="+mystring[-1]
		cbar.set_ticklabels(mystring)
		if int(np.log10(my_power_color_ticklabels))!=0:
			cbar.ax.set_title('1e'+str(int(np.log10(my_power_color_ticklabels))),fontsize=fontsize_ticks)
		return plt

	
	
	def plot_smoothing_factor(self,time_string=None,period_string=None,power_string=None,dashed_periods=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,cmap="jet",reverse_xaxis=False,reverse_yaxis=False,alpha_SNEZ=0.5):
		
		""" plot_smoothing_factor draws the half dimension of matrix M2 in the time-frequency plane. This gives the number of points over which smoothing is performed. No smoothing results in a value of 1 everywhere.
			Optional Inputs:
			- time_string=None: list of floats containing the location of the ticks for the time axis.
			- period_string=None: list of floats containing the location of the ticks for the period axis.
			- power_string=None: list of floats containing the location of the ticks for the power axis of the global scalogram.
			- dashed_periods=None: list of floats containing the periods for which a dashed line is drawn.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- cmap="jet": colormap for the scalogram. Other choices on http://matplotlib.org/users/colormaps.html
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			- alpha_SNEZ=0.5: Transparency for the SNEZ. It must take a value between 0 (completely transparent) and 1 (completely opaque). Only used if shannonnyquistexclusionzone=False in the method 'timefreq_analysis'.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
			-> plt.show(): to draw the figure
			-> plt.savefig(figname.pdf): to save a figure
			etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2017 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		from copy import copy
		
		# check inputs
		try:
			assert (time_string is None) or (type(time_string) is list)
		except AssertionError:
			print("Error at input 'time_string': must be None or of type 'list'")
			return
		if type(time_string) is list:
			for k in range(len(time_string)):
				try:
					assert type(time_string[k]) is float
				except AssertionError:
					print("Error at input 'time_string': must be a list containing floats")
					return
		try:
			assert (period_string is None) or (type(period_string) is list)
		except AssertionError:
			print("Error at input 'period_string': must be None or of type 'list'")
			return
		if type(period_string) is list:
			for k in range(len(period_string)):
				try:
					assert type(period_string[k]) is float
				except AssertionError:
					print("Error at input 'period_string': must be a list containing floats")
					return
		try:
			assert (power_string is None) or (type(power_string) is list)
		except AssertionError:
			print("Error at input 'power_string': must be None or of type 'list'")
			return
		if type(power_string) is list:
			for k in range(len(power_string)):
				try:
					assert type(power_string[k]) is float
				except AssertionError:
					print("Error at input 'power_string': must be a list containing floats")
					return
		try:
			assert (dashed_periods is None) or (type(dashed_periods) is list)
		except AssertionError:
			print("Error at input 'dashed_periods': must be None or of type 'list'")
			return
		if type(dashed_periods) is list:
			for k in range(len(dashed_periods)):
				try:
					assert type(dashed_periods[k]) is float
				except AssertionError:
					print("Error at input 'dashed_periods': must be a list containing floats")
					return
		try:
			assert type(cmap) is str
		except AssertionError:
			print("Error at input 'cmap': must be of type 'str'")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		try:
			assert ((type(alpha_SNEZ) is float) or (type(alpha_SNEZ) is int)) and (alpha_SNEZ>=0 and alpha_SNEZ<=1)
		except AssertionError:
			print("Error at input 'alpha_SNEZ': must of type float or int and must take a value in [0,1]")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		if period_string is None:
			myperiod=np.ceil(np.log2(self.period_cwt[0]))
			myperiod=2.**myperiod
			period_string=[myperiod]
			while True:
				myperiod*=2.
				if myperiod>self.period_cwt[-1]:
					break
				period_string.append(myperiod)
		mycontourf=plt.contourf(self.theta,np.log2(self.period_cwt),np.transpose(self.r),levels=range(self.r.max()+1),cmap=cmap)
		plt.fill_betweenx(np.log2(self.period_cwt),self.theta[0],self.coi1_smooth,edgecolors=None,facecolor='black')
		plt.fill_betweenx(np.log2(self.period_cwt),self.coi2_smooth,self.theta[-1],edgecolors=None,facecolor='black')
		if self.shannonnyquistexclusionzone is True:
			plt.fill_between(self.theta,np.log2(self.period_cwt[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_cwt),edgecolors=None,facecolor='black')
		elif self.shannonnyquistexclusionzone is False:
			plt.fill_between(self.theta,np.log2(self.period_cwt[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_cwt),edgecolors=None,facecolor='black',alpha=alpha_SNEZ)
		ax = plt.gca()
		ax.tick_params(length=5, width=1, color='w')
		if dashed_periods is not None:
			for k in range(len(dashed_periods)):
				plt.plot(self.theta,np.log2(dashed_periods[k])*np.ones(self.theta.size),'w--')
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel("Period"+self.t_label,fontsize=fontsize_axes)
		mytitle="Smoothing Factor"
		if time_string is not None:
			plt.xticks(time_string, time_string)
		plt.xticks(fontsize=fontsize_ticks)
		plt.yticks(np.log2(period_string), period_string, fontsize=fontsize_ticks)
		plt.xlim([self.theta[0], self.theta[-1]])
		plt.ylim([np.log2(self.period_cwt[0]), np.log2(self.period_cwt[-1])])
		plt.suptitle(mytitle, fontsize=fontsize_title)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		# Colorbar and its rescaling: min and max of the levels of color are defined over the non-shaded regions
		cbar=plt.colorbar(mycontourf)
		cbar.ax.set_ylabel("Smoothing Factor",fontsize=fontsize_axes)
		cbar.ax.tick_params(labelsize=fontsize_ticks)
		return plt
	
					 
						 
	def plot_cwtamplitude(self,with_global_amplitude=True,time_string=None,period_string=None,power_string=None,dashed_periods=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,left_padding=0.,right_padding=0.,middle_padding=0.,global_amplitude_xlabel="top",global_amplitude_xlabel_ticks="top",cmap="jet",nlevels=50,plot_coi="fill",linewidth_coi=1.0,plot_perlim2="fill",linewidth_perlim2=1.0,plot_ridges=False,k_skeleton=[],plot_band_filtering=False,linewidth_gampl=1.0,reverse_xaxis=False,reverse_yaxis=False,alpha_SNEZ=0.5,decimals=3,linewidth_ridges=1.0):
		
		""" plot_cwtamplitude generates the figure of the amplitude scalogram. It also generates the figure of the global amplitude scalogram. Only available if computes_amplitude is True (see 'timefreq_analysis').
			Optional Inputs:
			- with_global_amplitude=True: with the global amplitude or False.
			- time_string=None: list of floats containing the location of the ticks for the time axis.
			- period_string=None: list of floats containing the location of the ticks for the period axis.
			- power_string=None: list of floats containing the location of the ticks for the power axis of the global amplitude.
			- dashed_periods=None: list of floats containing the periods for which a dashed line is drawn.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- left_padding=0.: padding on the left of the amplitude figure.
			- right_padding=0.: padding on the right of the global amplitude.
			- middle_padding=0.: padding between the amplitude and the global amplitude.
			N.B.: paddings are only active if with_global_amplitude is True.
			WARNING with the paddings: If there is an overlap between subfigures, python can choose to not display one.
			- global_amplitude_xlabel="top": location of the xlabel for the global amplitude: "top" or "bottom".
			- global_amplitude_xlabel_ticks="top": location of the ticks of the xlabel for the global amplitude: "top" or "bottom".
			- cmap="jet": colormap for the amplitude. Other choices on http://matplotlib.org/users/colormaps.html
			- nlevels=50: number of automatically-chosen color levels. 
			- plot_coi="fill": plotting-type for the cone of influence: "fill" or "line".
			- linewidth_coi=1.0: linewidth for the coi (if plot_coi="line").
			- plot_perlim2="fill": plotting-type for the refinement of the SNEZ: "fill" or "line".
			- linewidth_perlim2=1.0: linewidth for perlim2 (if plot_perlim2="line").
			- plot_ridges=False: adds the ridges on the figure: True or False. There are drawn in white. Ridges are computed with function 'timefreq_ridges_filtering'.
			- k_skeleton=[]: list of integers or empty list, containing the indices of the ridges to be drawn in black. The smallest index corresponds to the ridge starting at the bottom-left coner of the time-period panel. The highest index corresponds to the ridge starting at the top-right corner.
			- plot_band_filtering=False: adds a shaded zone corresponding to where band filtering was performed, according to the function 'timefreq_band_filtering'. Must be True or False.
			- linewidth_gampl=1.0: linewidth for the global amplitude scalogram.
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			- alpha_SNEZ=0.5: Transparency for the SNEZ. It must take a value between 0 (completely transparent) and 1 (completely opaque). Only used if shannonnyquistexclusionzone=False in the method 'timefreq_analysis'.
			- decimals=3: Numbers of decimals for the colorbar scale ticks under the scientific notation (e.g. 1.567e-10).
			- linewidth_ridges=1.0: Linewidth for the ridges
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		
		# check inputs
		try:
			assert type(with_global_amplitude) is bool
		except AssertionError:
			print("Error at input 'with_global_amplitude': must be True or False")
			return
		try:
			assert (time_string is None) or (type(time_string) is list)
		except AssertionError:
			print("Error at input 'time_string': must be None or of type 'list'")
			return
		if type(time_string) is list:
			for k in range(len(time_string)):
				try:
					assert type(time_string[k]) is float
				except AssertionError:
					print("Error at input 'time_string': must be a list containing floats")
					return
		try:
			assert (period_string is None) or (type(period_string) is list)
		except AssertionError:
			print("Error at input 'period_string': must be None or of type 'list'")
			return
		if type(period_string) is list:
			for k in range(len(period_string)):
				try:
					assert type(period_string[k]) is float
				except AssertionError:
					print("Error at input 'period_string': must be a list containing floats")
					return
		try:
			assert (power_string is None) or (type(power_string) is list)
		except AssertionError:
			print("Error at input 'power_string': must be None or of type 'list'")
			return
		if type(power_string) is list:
			for k in range(len(power_string)):
				try:
					assert type(power_string[k]) is float
				except AssertionError:
					print("Error at input 'power_string': must be a list containing floats")
					return
		try:
			assert type(left_padding) is float
		except AssertionError:
			print("Error at input 'left_padding': must be of type 'float'")
			return
		try:
			assert type(right_padding) is float
		except AssertionError:
			print("Error at input 'right_padding': must be of type 'float'")
			return
		try:
			assert type(middle_padding) is float
		except AssertionError:
			print("Error at input 'middle_padding': must be of type 'float'")
			return
		try:
			assert (type(global_amplitude_xlabel) is str) and ((global_amplitude_xlabel.lower()=="bottom") or (global_amplitude_xlabel.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_amplitude_xlabel': must be 'bottom' or 'top'")
			return
		global_amplitude_xlabel=global_amplitude_xlabel.lower()
		try:
			assert (type(global_amplitude_xlabel_ticks) is str) and ((global_amplitude_xlabel_ticks.lower()=="bottom") or (global_amplitude_xlabel_ticks.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_amplitude_xlabel_ticks': must be 'bottom' or 'top'")
			return
		global_amplitude_xlabel_ticks=global_amplitude_xlabel_ticks.lower()
		try:
			assert type(cmap) is str
		except AssertionError:
			print("Error at input 'cmap': must be of type 'str'")
			return
		try:
			assert (type(nlevels) is int) and nlevels>0
		except AssertionError:
			print("Error at input 'nlevels': must be an integer >0")
			return
		try:
			assert (type(plot_coi) is str) and ((plot_coi.lower()=="fill") or (plot_coi.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_coi': must be 'fill' or 'line'")
			return
		plot_coi=plot_coi.lower()
		try:
			assert (type(linewidth_coi) is int) or (type(linewidth_coi) is float)
		except AssertionError:
			print("Error at input 'linewidth_coi': must be an integer or float")
			return
		try:
			assert (type(plot_perlim2) is str) and ((plot_perlim2.lower()=="fill") or (plot_perlim2.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_perlim2': must be 'fill' or 'line'")
			return
		plot_perlim2=plot_perlim2.lower()
		try:
			assert (type(linewidth_perlim2) is int) or (type(linewidth_perlim2) is float)
		except AssertionError:
			print("Error at input 'linewidth_perlim2': must be an integer or float")
			return
		try:
			assert type(plot_ridges) is bool
		except AssertionError:
			print("Error at input 'plot_ridges': must be True or False")
			return
		try:
			assert type(k_skeleton) is list
		except AssertionError:
			print("Error at input 'k_skeleton': must be of type 'list'")
			return
		for k in range(len(k_skeleton)):
			try:
				assert type(k_skeleton[k]) is int
			except AssertionError:
				print("Error at input 'k_skeleton': must be a list containing integers")
				return
		try:
			assert type(plot_band_filtering) is bool
		except AssertionError:
			print("Error at input 'plot_band_filtering': must be True or False")
			return
		try:
			assert (type(linewidth_gampl) is int) or (type(linewidth_gampl) is float)
		except AssertionError:
			print("Error at input 'linewidth_gampl': must be an integer or float")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		try:
			assert ((type(alpha_SNEZ) is float) or (type(alpha_SNEZ) is int)) and (alpha_SNEZ>=0 and alpha_SNEZ<=1)
		except AssertionError:
			print("Error at input 'alpha_SNEZ': must of type float or int and must take a value in [0,1]")
			return
		try:
			assert (type(decimals) is int) and decimals>=0
		except AssertionError:
			print("Error at input 'decimals': must be an integer >=0")
			return
		try:
			assert (type(linewidth_ridges) is int) or (type(linewidth_ridges) is float)
		except AssertionError:
			print("Error at input 'linewidth_ridges': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True
		except AssertionError:
			print("Error: plot_cwtamplitude cannot be applied")
			return
		if with_global_amplitude is True:
			try:
				assert self.computes_global_scalogram is True
			except AssertionError:
				print("Error: Global scalogram was not computed")
				print("=> drawing the figure without the global amplitude")
				with_global_amplitude=False
		if period_string is None:
			myperiod=np.ceil(np.log2(self.period_ampl[0]))
			myperiod=2.**myperiod
			period_string=[myperiod]
			while True:
				myperiod*=2.
				if myperiod>self.period_ampl[-1]:
					break
				period_string.append(myperiod)
		if with_global_amplitude is True:
			gs1 = gridspec.GridSpec(1,5)
			gs1.update(left=0.05+left_padding)
			plt.subplot(gs1[0,:-1])
		minminampl=np.amin(np.amin(self.cwtamplitude))
		maxmaxampl=np.amax(np.amax(self.cwtamplitude))
		mynlevels=int(float(nlevels)/(self.maxampl-self.minampl)*(maxmaxampl-minminampl))
		mycontourf=plt.contourf(self.theta,np.log2(self.period_ampl),np.transpose(self.cwtamplitude),mynlevels,vmin=self.minampl,vmax=self.maxampl,cmap=cmap)
		if (self.run_timefreq_ridges_filtering is True) and (plot_ridges is True):
			nsk=len(self.skeleton)
			for k in range(nsk):
				plt.plot(self.skeleton[k][0],np.log2(self.skeleton[k][1]),'w',linewidth=linewidth_ridges)
			for k in range(len(k_skeleton)):
				if k>=nsk:
					print( "WARNING: Element "+str(k)+" of 'k_skeleton' is not a correct value")
				else:
					plt.plot(self.skeleton[k_skeleton[k]][0],np.log2(self.skeleton[k_skeleton[k]][1]),'k',linewidth=linewidth_ridges)
		elif (self.run_timefreq_ridges_filtering is not True) and (plot_ridges is True):
			print("WARNING: function 'timefreq_ridges_filtering' was not run => unable to draw the ridges")
		if (self.run_timefreq_band_filtering is True) and (plot_band_filtering is True):
			for k in range(self.timefreq_band_filtered_signal_bounds.shape[0]):
				plt.fill_between(self.theta,np.log2(self.timefreq_band_filtered_signal_bounds[k,0]),np.log2(self.timefreq_band_filtered_signal_bounds[k,1]),edgecolors=None,facecolor='black',alpha=0.5)
		elif (self.run_timefreq_band_filtering is not True) and (plot_band_filtering is True):
			print("WARNING: function 'timefreq_band_filtering' was not run => unable to draw the bands")
		if plot_coi=="fill":
			plt.fill_betweenx(np.log2(self.period_ampl),self.t[0],self.coi1,edgecolors=None,facecolor='black',alpha=0.5)
			plt.fill_betweenx(np.log2(self.period_ampl),self.coi2,self.t[-1],edgecolors=None,facecolor='black',alpha=0.5)
		elif plot_coi=="line":
			plt.plot(self.coi1,np.log2(self.period_ampl),'k')
			plt.plot(self.coi2,np.log2(self.period_ampl),'k')
		plt.fill_betweenx(np.log2(self.period_ampl),self.theta[0],self.coi1_smooth,edgecolors=None,facecolor='black')
		plt.fill_betweenx(np.log2(self.period_ampl),self.coi2_smooth,self.theta[-1],edgecolors=None,facecolor='black')
		if self.shannonnyquistexclusionzone is True:
			plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_ampl),edgecolors=None,facecolor='black')
			if plot_perlim2=="fill":
				plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim2_smooth_ampl),edgecolors=None,facecolor='black',alpha=0.5)
			elif plot_perlim2=="line":
				plt.plot(self.theta,np.log2(self.perlim2_smooth_ampl),'k')
		elif self.shannonnyquistexclusionzone is False:
			plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_ampl),edgecolors=None,facecolor='black',alpha=alpha_SNEZ)
		ax = plt.gca()
		ax.tick_params(length=5, width=1, color='w')
		if dashed_periods is not None:
			for k in range(len(dashed_periods)):
				plt.plot(self.theta,np.log2(dashed_periods[k])*np.ones(self.theta.size),'w--')
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel("Period"+self.t_label,fontsize=fontsize_axes)
		if time_string is not None:
			plt.xticks(time_string, time_string)
		plt.xticks(fontsize=fontsize_ticks)
		plt.yticks(np.log2(period_string), period_string, fontsize=fontsize_ticks)
		plt.xlim([self.theta[0], self.theta[-1]])
		plt.ylim([np.log2(self.period_ampl[0]), np.log2(self.period_ampl[-1])])
		plt.suptitle("Wavelet Amplitude", fontsize=fontsize_title)
		if reverse_xaxis is True:
			plt.gca().invert_xaxis()
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		if with_global_amplitude is True:
			gs2 = gridspec.GridSpec(1,5)
			gs2.update(left=-0.3+middle_padding,right=0.95-right_padding)
			plt.subplot(gs2[0,-1])
			glob_ampl_min=np.amin(self.global_amplitude)
			glob_ampl_max=np.amax(self.global_amplitude)
			plt.plot(self.global_amplitude,np.log2(self.period_ampl),"b",linewidth=linewidth_gampl)
			myrange=np.arange(glob_ampl_min,glob_ampl_max,(glob_ampl_max-glob_ampl_min)/1000.)
			if dashed_periods is not None:
				for k in range(len(dashed_periods)):
					plt.plot(myrange,np.log2(dashed_periods[k])*np.ones(len(myrange)),'k--')
			plt.yticks(np.log2(period_string),[])
			plt.ylim([np.log2(self.period_ampl[0]), np.log2(self.period_ampl[-1])])
			plt.xlim([glob_ampl_min, glob_ampl_max])
			plt.xlabel("Power"+self.power_label,fontsize=fontsize_axes)
			ax=plt.gca()
			ax.xaxis.set_label_position(global_amplitude_xlabel)
			if power_string is not None:
				plt.xticks(power_string, power_string)
			if global_amplitude_xlabel_ticks=="top":
				plt.tick_params(axis='x', labelbottom='off', labeltop='on', labelsize=fontsize_ticks)
			elif global_amplitude_xlabel_ticks=="bottom":
				plt.tick_params(axis='x', labelbottom='on', labeltop='off', labelsize=fontsize_ticks)
			if reverse_yaxis is True:
				plt.gca().invert_yaxis()
		# Colorbar and its rescaling: min and max of the levels of color are defined over the non-shaded regions
		cbar=plt.colorbar(mycontourf)
		cbar.ax.set_ylabel("Amplitude"+self.mydata_label,fontsize=fontsize_axes)
		cbar.ax.tick_params(labelsize=fontsize_ticks)
		#cbar.set_clim(minminampl,maxmaxampl)
		mycontourf.set_clim(minminampl,maxmaxampl)
		my_color_ticks=np.linspace(minminampl,maxmaxampl,10)
		cbar.set_ticks(my_color_ticks)
		my_color_ticklabels=np.linspace(self.minampl,self.maxampl,10)
		my_power_color_ticklabels=float(10**(int(np.floor(np.log10(my_color_ticklabels[-1])))))
		my_color_ticklabels=[el/my_power_color_ticklabels for el in my_color_ticklabels]
		my_color_ticklabels=np.around(my_color_ticklabels,decimals=decimals)
		mystring=' '.join(str(e) for e in my_color_ticklabels).split(' ')
		if ('1' in mystring[0]) or ('2' in mystring[0]) or ('3' in mystring[0]) or ('4' in mystring[0]) or ('5' in mystring[0]) or ('6' in mystring[0]) or ('7' in mystring[0]) or ('8' in mystring[0]) or ('9' in mystring[0]):
			mystring[0]="<="+mystring[0]
		mystring[-1]=">="+mystring[-1]
		cbar.set_ticklabels(mystring)
		if int(np.log10(my_power_color_ticklabels))!=0:
			cbar.ax.set_title('1e'+str(int(np.log10(my_power_color_ticklabels))),fontsize=fontsize_ticks)
		return plt
						 
						 
						 
	def plot_cwtamplitude_squared(self,with_global_amplitude=True,time_string=None,period_string=None,power_string=None,dashed_periods=None,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,left_padding=0.,right_padding=0.,middle_padding=0.,global_amplitude_xlabel="top",global_amplitude_xlabel_ticks="top",cmap="jet",nlevels=50,plot_coi="fill",linewidth_coi=1.0,plot_perlim2="fill",linewidth_perlim2=1.0,plot_ridges=False,k_skeleton=[],plot_band_filtering=False,linewidth_gampl=1.0,reverse_xaxis=False,reverse_yaxis=False,alpha_SNEZ=0.5,decimals=3,linewidth_ridges=1.0,minampl_sq=None,maxampl_sq=None,permin_ext=None):
		
		""" plot_cwtamplitude_squared generates the figure of the squared amplitude scalogram. It also generates the figure of the global squared amplitude scalogram. Only available if computes_amplitude is True (see 'timefreq_analysis').
			Optional Inputs:
			- with_global_amplitude=True: with the global squared amplitude or False.
			- time_string=None: list of floats containing the location of the ticks for the time axis.
			- period_string=None: list of floats containing the location of the ticks for the period axis.
			- power_string=None: list of floats containing the location of the ticks for the power axis of the global squared amplitude.
			- dashed_periods=None: list of floats containing the periods for which a dashed line is drawn.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- left_padding=0.: padding on the left of the squared amplitude figure.
			- right_padding=0.: padding on the right of the global squared amplitude.
			- middle_padding=0.: padding between the squared amplitude and the global squared amplitude.
			N.B.: paddings are only active if with_global_amplitude is True.
			WARNING with the paddings: If there is an overlap between subfigures, python can choose to not display one.
			- global_amplitude_xlabel="top": location of the xlabel for the global squared amplitude: "top" or "bottom".
			- global_amplitude_xlabel_ticks="top": location of the ticks of the xlabel for the global squared amplitude: "top" or "bottom".
			- cmap="jet": colormap for the squared amplitude. Other choices on http://matplotlib.org/users/colormaps.html
			- nlevels=50: number of automatically-chosen color levels. 
			- plot_coi="fill": plotting-type for the cone of influence: "fill" or "line".
			- linewidth_coi=1.0: linewidth for the coi (if plot_coi="line").
			- plot_perlim2="fill": plotting-type for the refinement of the SNEZ: "fill" or "line".
			- linewidth_perlim2=1.0: linewidth for perlim2 (if plot_perlim2="line").
			- plot_ridges=False: adds the ridges on the figure: True or False. There are drawn in white. Ridges are computed with function 'timefreq_ridges_filtering'.
			- k_skeleton=[]: list of integers or empty list, containing the indices of the ridges to be drawn in black. The smallest index corresponds to the ridge starting at the bottom-left coner of the time-period panel. The highest index corresponds to the ridge starting at the top-right corner.
			- plot_band_filtering=False: adds a shaded zone corresponding to where band filtering was performed, according to the function 'timefreq_band_filtering'. Must be True or False.
			- linewidth_gampl=1.0: linewidth for the global squared amplitude scalogram.
			- reverse_xaxis=False: Reverse the horizontal axis if True
			- reverse_yaxis=False: Reverse the vertical axis if True
			- alpha_SNEZ=0.5: Transparency for the SNEZ. It must take a value between 0 (completely transparent) and 1 (completely opaque). Only used if shannonnyquistexclusionzone=False in the method 'timefreq_analysis'.
			- decimals=3: Numbers of decimals for the colorbar scale ticks under the scientific notation (e.g. 1.567e-10).
			- linewidth_ridges=1.0: Linewidth for the ridges
			- minampl_sq=None: Fix the lower bound of the color scale. Default (None) takes the min value of the squared amplitude scalogram over the non-shaded area.
			- maxampl_sq=None: Fix the higher bound of the color scale. Default (None) takes the max value of the squared amplitude scalogram over the non-shaded area.
			- permin_ext=None: If a float value is given, the SNEZ is broaden towards this given period. Available only if shannonnyquistexclusionzone=True in the method 'timefreq_analysis'.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
		
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gridspec
		from copy import copy
		
		# check inputs
		try:
			assert type(with_global_amplitude) is bool
		except AssertionError:
			print("Error at input 'with_global_amplitude': must be True or False")
			return
		try:
			assert (time_string is None) or (type(time_string) is list)
		except AssertionError:
			print("Error at input 'time_string': must be None or of type 'list'")
			return
		if type(time_string) is list:
			for k in range(len(time_string)):
				try:
					assert type(time_string[k]) is float
				except AssertionError:
					print("Error at input 'time_string': must be a list containing floats")
					return
		try:
			assert (period_string is None) or (type(period_string) is list)
		except AssertionError:
			print("Error at input 'period_string': must be None or of type 'list'")
			return
		if type(period_string) is list:
			for k in range(len(period_string)):
				try:
					assert type(period_string[k]) is float
				except AssertionError:
					print("Error at input 'period_string': must be a list containing floats")
					return
		try:
			assert (power_string is None) or (type(power_string) is list)
		except AssertionError:
			print("Error at input 'power_string': must be None or of type 'list'")
			return
		if type(power_string) is list:
			for k in range(len(power_string)):
				try:
					assert type(power_string[k]) is float
				except AssertionError:
					print("Error at input 'power_string': must be a list containing floats")
					return
		try:
			assert type(left_padding) is float
		except AssertionError:
			print("Error at input 'left_padding': must be of type 'float'")
			return
		try:
			assert type(right_padding) is float
		except AssertionError:
			print("Error at input 'right_padding': must be of type 'float'")
			return
		try:
			assert type(middle_padding) is float
		except AssertionError:
			print("Error at input 'middle_padding': must be of type 'float'")
			return
		try:
			assert (type(global_amplitude_xlabel) is str) and ((global_amplitude_xlabel.lower()=="bottom") or (global_amplitude_xlabel.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_amplitude_xlabel': must be 'bottom' or 'top'")
			return
		global_amplitude_xlabel=global_amplitude_xlabel.lower()
		try:
			assert (type(global_amplitude_xlabel_ticks) is str) and ((global_amplitude_xlabel_ticks.lower()=="bottom") or (global_amplitude_xlabel_ticks.lower()=="top"))
		except AssertionError:
			print("Error at input 'global_amplitude_xlabel_ticks': must be 'bottom' or 'top'")
			return
		global_amplitude_xlabel_ticks=global_amplitude_xlabel_ticks.lower()
		try:
			assert type(cmap) is str
		except AssertionError:
			print("Error at input 'cmap': must be of type 'str'")
			return
		try:
			assert (type(nlevels) is int) and nlevels>0
		except AssertionError:
			print("Error at input 'nlevels': must be an integer >0")
			return
		try:
			assert (type(plot_coi) is str) and ((plot_coi.lower()=="fill") or (plot_coi.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_coi': must be 'fill' or 'line'")
			return
		plot_coi=plot_coi.lower()
		try:
			assert (type(linewidth_coi) is int) or (type(linewidth_coi) is float)
		except AssertionError:
			print("Error at input 'linewidth_coi': must be an integer or float")
			return
		try:
			assert (type(plot_perlim2) is str) and ((plot_perlim2.lower()=="fill") or (plot_perlim2.lower()=="line"))
		except AssertionError:
			print("Error at input 'plot_perlim2': must be 'fill' or 'line'")
			return
		plot_perlim2=plot_perlim2.lower()
		try:
			assert (type(linewidth_perlim2) is int) or (type(linewidth_perlim2) is float)
		except AssertionError:
			print("Error at input 'linewidth_perlim2': must be an integer or float")
			return
		try:
			assert type(plot_ridges) is bool
		except AssertionError:
			print("Error at input 'plot_ridges': must be True or False")
			return
		try:
			assert type(k_skeleton) is list
		except AssertionError:
			print("Error at input 'k_skeleton': must be of type 'list'")
			return
		for k in range(len(k_skeleton)):
			try:
				assert type(k_skeleton[k]) is int
			except AssertionError:
				print("Error at input 'k_skeleton': must be a list containing integers")
				return
		try:
			assert type(plot_band_filtering) is bool
		except AssertionError:
			print("Error at input 'plot_band_filtering': must be True or False")
			return
		try:
			assert (type(linewidth_gampl) is int) or (type(linewidth_gampl) is float)
		except AssertionError:
			print("Error at input 'linewidth_gampl': must be an integer or float")
			return
		try:
			assert type(reverse_xaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_xaxis': must be True or False")
			return
		try:
			assert type(reverse_yaxis) is bool
		except AssertionError:
			print("Error at input 'reverse_yaxis': must be True or False")
			return
		try:
			assert ((type(alpha_SNEZ) is float) or (type(alpha_SNEZ) is int)) and (alpha_SNEZ>=0 and alpha_SNEZ<=1)
		except AssertionError:
			print("Error at input 'alpha_SNEZ': must of type float or int and must take a value in [0,1]")
			return
		try:
			assert (type(decimals) is int) and decimals>=0
		except AssertionError:
			print("Error at input 'decimals': must be an integer >=0")
			return
		try:
			assert (type(linewidth_ridges) is int) or (type(linewidth_ridges) is float)
		except AssertionError:
			print("Error at input 'linewidth_ridges': must be an integer or float")
			return
		try:
			assert (minampl_sq is None) or (type(minampl_sq) is float and minampl_sq>=0.)
		except AssertionError:
			print("Error at input 'minampl_sq': must be a float >=0")
			return
		try:
			assert (maxampl_sq is None) or (type(maxampl_sq) is float and maxampl_sq>=0.)
		except AssertionError:
			print("Error at input 'maxampl_sq': must be a float, >=0, and >= 'minampl_sq'")
			return
		try:
			assert (permin_ext is None) or (type(permin_ext) is float and permin_ext>=0.)
		except AssertionError:
			print("Error at input 'permin_ext': must be a float >=0")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True
		except AssertionError:
			print("Error: plot_cwtamplitude_squared cannot be applied")
			return
		if with_global_amplitude is True:
			try:
				assert self.computes_global_scalogram is True
			except AssertionError:
				print("Error: Global scalogram was not computed")
				print("=> drawing the figure without the global amplitude")
				with_global_amplitude=False
		if (self.shannonnyquistexclusionzone is False) and (permin_ext is not None):
			print("parameter 'permin_ext' changed to 'None', since 'shannonnyquistexclusionzone' is 'False'")
			permin_ext=None
		if period_string is None:
			if (permin_ext is not None) and (permin_ext<self.period_ampl[0]):
				myperiod=np.ceil(np.log2(permin_ext))
			else:
				myperiod=np.ceil(np.log2(self.period_ampl[0]))
			myperiod=2.**myperiod
			period_string=[myperiod]
			while True:
				myperiod*=2.
				if myperiod>self.period_ampl[-1]:
					break
				period_string.append(myperiod)
		if with_global_amplitude is True:
			gs1 = gridspec.GridSpec(1,5)
			gs1.update(left=0.05+left_padding)
			plt.subplot(gs1[0,:-1])
		if minampl_sq is None:
			minampl_sq=self.minampl_sq
		if maxampl_sq is None:
			maxampl_sq=self.maxampl_sq
		try:
			assert maxampl_sq>=minampl_sq
		except AssertionError:
			print("Error: variable 'maxampl_sq' must be greater than or equal to 'minampl_sq'")
			return
		minminampl_sq=np.amin(np.amin(self.cwtamplitude**2))
		maxmaxampl_sq=np.amax(np.amax(self.cwtamplitude**2))
		mycwtamplitude=copy(self.cwtamplitude**2)
		if maxampl_sq<maxmaxampl_sq:
			mycwtamplitude[mycwtamplitude>maxampl_sq]=maxampl_sq
			myvmax=maxampl_sq
		else:
			myvmax=maxampl_sq
		if minampl_sq>minminampl_sq:
			mycwtamplitude[mycwtamplitude<minampl_sq]=minampl_sq
			myvmin=minampl_sq
		else:
			myvmin=minampl_sq
		mycontourf=plt.contourf(self.theta,np.log2(self.period_ampl),np.transpose(mycwtamplitude),nlevels,vmin=myvmin,vmax=myvmax,cmap=cmap)
		if (self.run_timefreq_ridges_filtering is True) and (plot_ridges is True):
			nsk=len(self.skeleton)
			for k in range(nsk):
				plt.plot(self.skeleton[k][0],np.log2(self.skeleton[k][1]),'w',linewidth=linewidth_ridges)
			for k in range(len(k_skeleton)):
			#for k in range(len(kskel)):
				if k>=nsk:
					print("WARNING: Element "+str(k)+" of 'k_skeleton' is not a correct value")
				else:
					plt.plot(self.skeleton[k_skeleton[k]][0],np.log2(self.skeleton[k_skeleton[k]][1]),'k',linewidth=linewidth_ridges)
		elif (self.run_timefreq_ridges_filtering is not True) and (plot_ridges is True):
			print("WARNING: function 'timefreq_ridges_filtering' was not run => unable to draw the ridges")
		if (self.run_timefreq_band_filtering is True) and (plot_band_filtering is True):
			for k in range(self.timefreq_band_filtered_signal_bounds.shape[0]):
				plt.fill_between(self.theta,np.log2(self.timefreq_band_filtered_signal_bounds[k,0]),np.log2(self.timefreq_band_filtered_signal_bounds[k,1]),edgecolors=None,facecolor='black',alpha=0.5)
		elif (self.run_timefreq_band_filtering is not True) and (plot_band_filtering is True):
			print("WARNING: function 'timefreq_band_filtering' was not run => unable to draw the bands")
		if plot_coi=="fill":
			plt.fill_betweenx(np.log2(self.period_ampl),self.t[0],self.coi1,edgecolors=None,facecolor='black',alpha=0.5)
			plt.fill_betweenx(np.log2(self.period_ampl),self.coi2,self.t[-1],edgecolors=None,facecolor='black',alpha=0.5)
		elif plot_coi=="line":
			plt.plot(self.coi1,np.log2(self.period_ampl),'k')
			plt.plot(self.coi2,np.log2(self.period_ampl),'k')
		plt.fill_betweenx(np.log2(self.period_ampl),self.theta[0],self.coi1_smooth,edgecolors=None,facecolor='black')
		plt.fill_betweenx(np.log2(self.period_ampl),self.coi2_smooth,self.theta[-1],edgecolors=None,facecolor='black')
		if self.shannonnyquistexclusionzone is True:
			if (permin_ext is not None) and (permin_ext<self.period_ampl[0]):
				plt.fill_between(self.theta,np.log2(permin_ext)*np.ones(self.theta.size),np.log2(self.perlim1_smooth_ampl),edgecolors=None,facecolor='black')
			else:
				plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_ampl),edgecolors=None,facecolor='black')
			if plot_perlim2=="fill":
				plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim2_smooth_scal),edgecolors=None,facecolor='black',alpha=0.5)
			elif plot_perlim2=="line":
				plt.plot(self.theta,np.log2(self.perlim2_smooth_scal),'k')
		elif self.shannonnyquistexclusionzone is False:
			plt.fill_between(self.theta,np.log2(self.period_ampl[0])*np.ones(self.theta.size),np.log2(self.perlim1_smooth_ampl),edgecolors=None,facecolor='black',alpha=alpha_SNEZ)
		ax = plt.gca()
		ax.tick_params(length=5, width=1, color='w')
		if dashed_periods is not None:
			for k in range(len(dashed_periods)):
				plt.plot(self.theta,np.log2(dashed_periods[k])*np.ones(self.theta.size),'w--')
		plt.xlabel(self.t_axis_label+self.t_label,fontsize=fontsize_axes)
		plt.ylabel("Period"+self.t_label,fontsize=fontsize_axes)
		if time_string is not None:
			plt.xticks(time_string, time_string)
		plt.xticks(fontsize=fontsize_ticks)
		plt.yticks(np.log2(period_string), period_string, fontsize=fontsize_ticks)
		plt.xlim([self.theta[0], self.theta[-1]])
		if (permin_ext is not None) and (permin_ext<self.period_ampl[0]):
			plt.ylim([np.log2(permin_ext), np.log2(self.period_ampl[-1])])
		else:
			plt.ylim([np.log2(self.period_ampl[0]), np.log2(self.period_ampl[-1])])
		plt.suptitle("Wavelet Squared Amplitude", fontsize=fontsize_title)
		if reverse_yaxis is True:
			plt.gca().invert_yaxis()
		if with_global_amplitude is True:
			gs2 = gridspec.GridSpec(1,5)
			gs2.update(left=-0.3+middle_padding,right=0.95-right_padding)
			plt.subplot(gs2[0,-1])
			glob_ampl_min=np.amin(self.global_amplitude**2)
			glob_ampl_max=np.amax(self.global_amplitude**2)
			plt.plot(self.global_amplitude**2,np.log2(self.period_ampl),"b",linewidth=linewidth_gampl)
			myrange=np.arange(glob_ampl_min,glob_ampl_max,(glob_ampl_max-glob_ampl_min)/1000.)
			if dashed_periods is not None:
				for k in range(len(dashed_periods)):
					plt.plot(myrange,np.log2(dashed_periods[k])*np.ones(len(myrange)),'k--')
			plt.yticks(np.log2(period_string),[])
			if (permin_ext is not None) and (permin_ext<self.period_ampl[0]):
				plt.ylim([np.log2(permin_ext), np.log2(self.period_ampl[-1])])
			else:
				plt.ylim([np.log2(self.period_ampl[0]), np.log2(self.period_ampl[-1])])
			plt.xlim([glob_ampl_min, glob_ampl_max])
			plt.xlabel("Power"+self.power_label,fontsize=fontsize_axes)
			ax=plt.gca()
			ax.xaxis.set_label_position(global_amplitude_xlabel)
			if power_string is not None:
				plt.xticks(power_string, power_string)
			if global_amplitude_xlabel_ticks=="top":
				plt.tick_params(axis='x', labelbottom='off', labeltop='on', labelsize=fontsize_ticks)
			elif global_amplitude_xlabel_ticks=="bottom":
				plt.tick_params(axis='x', labelbottom='on', labeltop='off', labelsize=fontsize_ticks)
			if reverse_yaxis is True:
				plt.gca().invert_yaxis()
		# Colorbar and its rescaling: min and max of the levels of color are defined over the non-shaded regions
		cbar=plt.colorbar(mycontourf)
		cbar.ax.set_ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		cbar.ax.tick_params(labelsize=fontsize_ticks)
		my_color_ticks=np.linspace(max(minampl_sq,minminampl_sq),min(maxampl_sq,maxmaxampl_sq),10)
		cbar.set_ticks(my_color_ticks)
		my_color_ticklabels=copy(my_color_ticks)
		my_power_color_ticklabels=float(10**(int(np.floor(np.log10(my_color_ticklabels[-1])))))
		my_color_ticklabels=[el/my_power_color_ticklabels for el in my_color_ticklabels]
		my_color_ticklabels=np.around(my_color_ticklabels,decimals=decimals)
		mystring=' '.join(str(e) for e in my_color_ticklabels).split(' ')
		if minampl_sq>minminampl_sq and (('1' in mystring[0]) or ('2' in mystring[0]) or ('3' in mystring[0]) or ('4' in mystring[0]) or ('5' in mystring[0]) or ('6' in mystring[0]) or ('7' in mystring[0]) or ('8' in mystring[0]) or ('9' in mystring[0])):
			mystring[0]="<="+mystring[0]
		if maxampl_sq<maxmaxampl_sq:
			mystring[-1]=">="+mystring[-1]
		cbar.set_ticklabels(mystring)
		if int(np.log10(my_power_color_ticklabels))!=0:
			cbar.ax.set_title('1e'+str(int(np.log10(my_power_color_ticklabels))),fontsize=fontsize_ticks)
		return plt
						 
						 
						 
	def plot_global_scalogram(self,xaxis="period",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='xx-small',linewidth_gscal=1.0,linewidth_gcl=1.0):
		
		""" plot_global_scalogram generates the figure of the global scalogram and its confidence levels.
			Optional Inputs:
			- xaxis="period": choice for the horizontal figure axis: "frequency" or "period".
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='xx-small': fontsize for the figure legend.
			- linewidth_gscal=1.0: linewidth for the global scalogram. 
			- linewidth_gcl=1.0: linewidth for the confidence levels.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt
		
		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth_gscal) is int) or (type(linewidth_gscal) is float)
		except AssertionError:
			print("Error at input 'linewidth_gscal': must be an integer or float")
			return
		try:
			assert (type(linewidth_gcl) is int) or (type(linewidth_gcl) is float)
		except AssertionError:
			print("Error at input 'linewidth_gcl': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_global_scalogram is True
		except AssertionError:
			print("Error: plot_global_scalogram cannot be applied")
			return
		if xaxis=="period":
			mypow=1
		elif xaxis=="frequency":
			mypow=-1
		if loglog is False:
			plt.plot(self.period_cwt**mypow,self.global_scalogram,"k",label="Data global scal.",linewidth=linewidth_gscal)
		elif loglog is True:
			plt.loglog(self.period_cwt**mypow,self.global_scalogram,"k",label="Data global scal.",linewidth=linewidth_gscal)
		
		plt.legend(fancybox=True,fontsize=fontsize_legend,bbox_to_anchor=(1.1, 1.05))
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		if self.signif_level_type!="" and self.percentile_cwt.size>0:
			plt.suptitle("Global Scalogram and Confidence levels", fontsize=fontsize_title)
		else:
			plt.suptitle("Global Scalogram", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt
						 
						 


									 
						 
						 
	def plot_global_amplitude(self,xaxis="period",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,linewidth=1.0):

		""" plot_global_amplitude generates the figure of the global amplitude scalogram. Only available if computes_amplitude is True (see 'timefreq_analysis').
			Optional Inputs:
			- xaxis="period": choice for the horizontal figure axis: "frequency" or "period".
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- linewidth=1.0: linewidth.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""

		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth) is int) or (type(linewidth) is float)
		except AssertionError:
			print("Error at input 'linewidth': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True
		except AssertionError:
			print("Error: plot_global_amplitude cannot be applied")
			return
		try:
			assert self.computes_global_scalogram is True
		except AssertionError:
			print("Error: plot_global_amplitude cannot be applied")
			return
		if xaxis=="period":
			mypow=1
		elif xaxis=="frequency":
			mypow=-1
		if loglog is False:
			plt.plot(self.period_ampl**mypow,self.global_amplitude,linewidth=linewidth)
		elif loglog is True:
			plt.loglog(self.period_ampl**mypow,self.global_amplitude,linewidth=linewidth)
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Amplitude"+self.mydata_label,fontsize=fontsize_axes)
		plt.title("Global Amplitude", fontsize=fontsize_title)
		plt.tick_params(labelsize=fontsize_ticks)
		return plt

						 
						 
	def plot_global_amplitude_vs_global_scalogram(self,xaxis="period",loglog=False,fontsize_title=14,fontsize_axes=12,fontsize_ticks=12,fontsize_legend='xx-small',linewidth=1.0):
	
		""" plot_global_amplitude_vs_global_scalogram generates the figure comparing of the global amplitude scalogram and the global weighted scalogram. Only available if computes_amplitude is True and weighted_CWT is True (see 'timefreq_analysis').
			Optional Inputs:
			- xaxis="period": choice for the horizontal figure axis: "frequency" or "period".
			- loglog=False: draw axis in log-log if True.
			- fontsize_title=14: fontsize for the figure title.
			- fontsize_axes=12: fontsize for the figure axes.
			- fontsize_ticks=12: fontsize for the figure ticks.
			- fontsize_legend='xx-small': fontsize for the figure legend.
			- linewidth=1.0: linewidth.
			Outputs:
			- plt: matplotlib.pyplot object that gives the user an access to the figure.
				-> plt.show(): to draw the figure
				-> plt.savefig(figname.pdf): to save a figure
				etc. See matplotlib documentation.
			-----------------------------
			This is part of WAVEPAL
			(C) 2016 G. Lenoir"""
	
		import matplotlib.pyplot as plt

		# check inputs
		try:
			assert (type(xaxis) is str) and ((xaxis.lower()=="frequency") or (xaxis.lower()=="period"))
		except AssertionError:
			print("Error at input 'xaxis': Must be 'frequency' or 'period'")
			return
		xaxis=xaxis.lower()
		try:
			assert type(loglog) is bool
		except AssertionError:
			print("Error at input 'loglog': Must be True or False")
			return
		try:
			assert (type(linewidth) is int) or (type(linewidth) is float)
		except AssertionError:
			print("Error at input 'linewidth': must be an integer or float")
			return
		# check that some functions were previously run
		try:
			assert self.run_timefreq_analysis is True
		except AssertionError:
			print("Error: Must have run function 'timefreq_analysis'")
			return
		try:
			assert self.computes_cwtamplitude is True and self.weighted_CWT is True
		except AssertionError:
			print("Error: plot_global_amplitude_vs_global_scalogram cannot be applied")
			return
		try:
			assert self.computes_global_scalogram is True
		except AssertionError:
			print("Error: plot_global_amplitude_vs_global_scalogram cannot be applied")
			return
		if xaxis=="period":
			mypow=1
		elif xaxis=="frequency":
			mypow=-1
		if loglog is False:
			plt.plot(self.period_cwt**mypow,self.global_amplitude**2,label="Squared global amplitude",linewidth=linewidth)
			plt.plot(self.period_cwt**mypow,self.global_scalogram,label="Weighted global scalogram",linewidth=linewidth)
		elif loglog is True:
			plt.loglog(self.period_cwt**mypow,self.global_amplitude**2,label="Squared global amplitude",linewidth=linewidth)
			plt.loglog(self.period_cwt**mypow,self.global_scalogram,label="Weighted global scalogram",linewidth=linewidth)
		if xaxis=="period":
			plt.xlabel("Period"+self.t_label,fontsize=fontsize_axes)
		elif xaxis=="frequency":
			plt.xlabel("Frequency"+self.freq_label,fontsize=fontsize_axes)
		plt.ylabel("Power"+self.power_label,fontsize=fontsize_axes)
		plt.suptitle("Squared Global Amplitude vs Weighted Global Scalogram",fontsize=fontsize_title)
		plt.legend(fancybox=True,fontsize=fontsize_legend,bbox_to_anchor=(1.1, 1.05))
		return plt


