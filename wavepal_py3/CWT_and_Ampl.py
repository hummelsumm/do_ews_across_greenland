import numpy as np
import numpy.linalg as la

def CWT_and_Ampl(time,data,theta,scale_ind,myscale,myprojvec,Vmat,pol_degree,my_weight_cwt,scalelim1_ind,n_outside_scalelim1,w0):
	
	""" CWT_and_Ampl computes the basic bricks for all the computations related to the scalogram and amplitude scalogram in 'Wavepal' class. It is called once for each scale.
		Inputs:
		- time [1-dim numpy array of floats]: the times of the time series
		- data [1-dim numpy array of floats - size=time.size]: the data at the times given by 'time'.
		- theta [1-dim numpy array of floats]: the times at which the CWT is to be computed
		- scale_ind [int]: the index of the scale at which the CWT is to be computed
		- myscale [float]: the scale at which the CWT is to be computed
		- myprojvec [numpy array of floats - dimension=(time.size,pol_degree+3)]: array with content related to the trend. This is an output from 'trend_vectors', in 'Wavepal' class.
		- Vmat [numpy array of floats - dimension=(time.size,pol_degree+3)]: array with content related to the trend. This is an output from 'trend_vectors', in 'Wavepal' class.
		- pol_degree [int]: degree of the polynomial trend. pol_degree=-1 means no trend.
		- my_weight_cwt [1-dim numpy array of floats - size=theta.size]: the weights, for each value of theta, to compute the weighted scalogram.
		- scalelim1_ind [1-dim numpy array of floats - size=theta.size]: indices of the scales defining the border of the Shannon-Nyquist exclusion zone. 
		- n_outside_scalelim1 [int]: number of theta's outside the Shannon-Nyquist exclusion zone at the scale given by myscale.
		- w0 [float]: the usual parameter for the Morlet wavelet
		Outputs:
		- Ampl [numpy array of floats - dimension=(2,n_outside_scalelim1)]: Amplitude scalogram. Ampl[0,:] gives the amplitudes of the 'cos' part and Ampl[1,:] gives the amplitudes of the 'sin' part, for all the theta's outside the Shannon-Nyquist exclusion zone, at the scale given by 'myscale'. 
		- M2 [numpy array of floats - dimension=(time.size,2*n_outside_scalelim1)]: array containing the vectors on which we perform the orthogonal projection, in order to compute the scalogram. See:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. II. Extension to time-frequency Analysis', G. Lenoir and M. Crucifix
		- mytheta [1-dim numpy array of floats - size=n_outside_scalelim1]: values of 'theta' outside the Shannon-Nyquist exclusion zone.
		- mytheta_ind [1-dim numpy array of ints - size=n_outside_scalelim1]: indices of 'theta' corresponding to 'mytheta' values.
		-------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	c0=1./2./w0**2
	N=time.size
	Q=theta.size
	M2=np.zeros((N,2*n_outside_scalelim1))
	Ampl=np.zeros((2,n_outside_scalelim1))
	mytheta=np.zeros(n_outside_scalelim1)
	mytheta_ind=np.zeros(n_outside_scalelim1,dtype=int)
	count=-1
	for k in range(Q):
		if scale_ind<scalelim1_ind[k]:
			continue
		else:
			count+=1
			theta_k=theta[k]
			mytheta[count]=theta_k
			mytheta_ind[count]=k
			gvec=np.exp(-c0/myscale**2*(time-theta_k)**2)
			mycos=np.cos(time/myscale)
			Vmat[:,pol_degree+1]=mycos[:]
			mycos[:]=gvec*mycos
			for p in range(pol_degree+1):     # Gram-Schmidt
				h=myprojvec[:,p]
				mycos[:]-=np.dot(h,mycos)*h
			mycos[:]/=la.norm(mycos)
			myprojvec[:,pol_degree+1]=mycos[:]
			mysin=np.sin(time/myscale)
			Vmat[:,pol_degree+2]=mysin[:]
			mysin[:]=gvec*mysin
			for p in range(pol_degree+2):     # Gram-Schmidt
				h=myprojvec[:,p]
				mysin[:]-=np.dot(h,mysin)*h
			mysin[:]/=la.norm(mysin)
			myprojvec[:,pol_degree+2]=mysin[:]
			M2[:,2*count]=mycos[:]*my_weight_cwt[k]
			M2[:,2*count+1]=mysin[:]*my_weight_cwt[k]
			# Amplitude
			myprojvec_transp=np.transpose(myprojvec)
			amat=np.dot(myprojvec_transp,Vmat)
			bmat=np.dot(myprojvec_transp,data)
			mysol=la.solve(amat,bmat)
			# Check that the solution is correct
			try:
				assert np.allclose(np.dot(amat,mysol),bmat)==True
			except AssertionError:
				print("WARNING: Error when computing signal Amplitude: Matrix is singular")
			finally:   # this is executed even if there is an AssertionError
				Ampl[:,count]=mysol[-2:]

	return Ampl, M2, mytheta, mytheta_ind

