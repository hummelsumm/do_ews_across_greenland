# using Pkg
# Pkg.activate(".")

using Distributions, StatsBase
using LinearRegressionKit, DataFrames, StatsModels
using ImageFiltering
using HypothesisTests
using JLD2, FileIO
using DelimitedFiles
using Dierckx #for Interpolations
using FFTW, StatsFuns, Statistics, SpecialFunctions
using Polynomials

function cheby_highpass(cutoff, fs, order, rp)
    nyq = 0.5*fs
    normal_cutoff = cutoff / nyq
    responsetype = Highpass(normal_cutoff)
    design = Chebyshev1(order, rp)
    return digitalfilter(responsetype, design)
end

function cheby_highpass_filter(x, cutoff, fs, order, rp)
    che= cheby_highpass(cutoff, fs, order, rp)
    y = filtfilt(che, x)
    return y
end

function cheby_lowpass(cutoff, fs, order, rp)
    nyq = 0.5*fs
    normal_cutoff = cutoff / nyq
    responsetype = Lowpass(normal_cutoff)
    design = Chebyshev1(order, rp)
    return digitalfilter(responsetype, design)
end

function cheby_lowpass_filter(x, cutoff, fs, order, rp)
    che= cheby_lowpass(cutoff, fs, order, rp)
    y = filtfilt(che, x)
    return y
end

function runmean(x, w; only_full = true)
    n = length(x)
    xs = Array{Union{Missing,Float64}}(missing,n)
    half_w = div(w, 2)
    for i in 1:half_w
        if !only_full && i<=n
            xs[i] = mean(x[1:min(i + half_w, n)])
        end
    end
    for i in (n - half_w + 1):n
        if !only_full && 1<=i<=n
            xs[i] = mean(x[max(1,i - half_w):end])
        end
    end
    for i in (half_w + 1):(n - half_w)
        if 1<=i<=n
            xs[i] = mean(x[(i - half_w):min(i + half_w, n)])
        end
    end
    return xs
end

function runvar(x, w; only_full = true)
    n = length(x)
    xs = Array{Union{Missing,Float64}}(missing,n)
    half_w = div(w, 2)
    for i in 1:half_w
        if !only_full && i<=n
            xs[i] = var(x[1:min(i + half_w, n)])
        end
    end
    for i in (n - half_w + 1):n
        if !only_full && 1<=i<=n
            xs[i] = var(x[max(1,i - half_w):end])
        end
    end
    for i in (half_w + 1):(n - half_w)
        if 1<=i<=n
            xs[i] = var(x[max(1,i - half_w):min(i + half_w,n)])
        end
    end
    return xs
end

function sample_autocor1(x,k)
    return only(autocor(x,[k], demean=true))
end

function sample_autocor2(x,k)
    mx = mean(x)
    u = mean((x[1+k:end] .- mx).*(x[1:end-k] .- mx))
    return u/var(x)
end

function sample_autocor3(x,k)
    y = x .- mean(x)
    dot(y[1:end-k], y[1+k:end])/dot(y,y)
end

function sample_autocor4(x,k)
    return only(autocor(x,[k], demean=false))
end

function sample_autocor(x,k)
    n= length(x)
    mx = mean(x)
    uu=0
    ll = 0
    for i = 1:(n-k)
        uu += (x[i+k] - mx)*(x[i] - mx)
    end
    for i = 1:n
        ll += (x[i] - mx)^2
    end
    return uu/ll
end


function runac(x, w; k=1, v=1, only_full = true)
    if v==0
        ac = sample_autocor
    elseif v==1
        ac = sample_autocor1
    elseif v==2
        ac = sample_autocor2
    elseif v==3
        ac = sample_autocor3
    elseif v==4
        ac = sample_autocor4
    end
    n = length(x)
    xs = Array{Union{Missing,Float64}}(missing,n)
    half_w = div(w, 2)
    for i in 1:half_w
        if !only_full && i<=n
            xs[i] = ac(x[1:min(i+half_w, n)],k) #aa
        end
    end
    for i in (n - half_w + 1):n
        if !only_full && 1<=i<=n
            xs[i] = ac(x[max(i - half_w,1):end],k)
        end
    end
    for i in (half_w + 1):(n - half_w)
        if 1<=i<=n
            xs[i] = ac(x[max(i - half_w,1):min(i+half_w, n)],k)
        end
    end
    return xs
end

function get_slopes(timeh, dat)
    yh = Array{Float64}(dat[findall(!ismissing, dat)])
    xh = timeh[findall(!ismissing, dat)]
    if length(yh) > 2 #need more than 2 data points
        slope = Polynomials.fit(xh, yh, 1)[1]
    else
        slope = nothing
    end
    slope
end

################################################
############## Wavelet functions ###############
################################################

### The following is an adapted version of the WaveletAnalysis code by Bernd Blasius (bernd.blasius@gmail.com) 
### on https://github.com/berndblasius/WaveletAnalysis, last visited on 06.12.2024:

# This code is based on a straightforward Julia translation of the
# code for 1d Wavelet transform from Torrence and Compo
#
# There is no intention to have a fancy, performant, or elegant "translation"
# it just should be able to "do stuff"
# The copyright remains by Torrence and Compo (see below)
# This code is provided without any warranties
# by Bernd Blasius (bernd.blasius@gmail.com)
#
#  HERE GOES THE ORIGINAL COPYRIGHT
#----------------------------------------------------------------------------
#   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made. This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------


heavi(x) = x > 0. ? 1. : 0.    # Heaviside function

function wave_bases(mother,k,scale,param)
    #   Computes the wavelet function as a function of Fourier frequency,
    #   used for the wavelet transform in Fourier space.
    #   (This program is called automatically by WAVELET)
    #
    # INPUTS:
    #    mother = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
    #    k = a vector, the Fourier frequencies at which to calculate the wavelet
    #    scale = a number, the wavelet scale
    #    param = the nondimensional parameter for the wavelet function
    #
    # OUTPUTS:
    #    daughter = a vector, the wavelet function
    #    fourier_factor = the ratio of Fourier period to scale
    #    coi = a number, the cone-of-influence size at the scale
    #    dofmin = a number, degrees of freedom for each point in the wavelet power
    #             (either 2 for Morlet and Paul, or 1 for the DOG)
    #             (either 2 for Morlet and Paul, or 1 for the DOG)
    #
    #    Copyright statement of the original Matlab version
    #----------------------------------------------------------------------------
    #   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
    #   University of Colorado, Program in Atmospheric and Oceanic Sciences.
    #   This software may be used, copied, or redistributed as long as it is not
    #   sold and this copyright notice is reproduced on each copy made.  This
    #   routine is provided as is without any express or implied warranties
    #   whatsoever.
    #----------------------------------------------------------------------------
        mother = uppercase(mother)
        n = length(k)
        if mother == "MORLET"  #-----------------------------------  Morlet
           if param == -1; param = 6.; end
           k0 = param
           expnt = -(scale*k .- k0).^2 / 2 .* heavi.(k)
           norm = sqrt(scale*k[2])*(pi^(-0.25))*sqrt(n)      # total energy=N   [Eqn(7)]
           daughter = norm*exp.(expnt)
           daughter = daughter .* heavi.(k)                  # Heaviside step function
           fourier_factor = (4.0*pi)/(k0 + sqrt(2.0 + k0^2)) # Scale-->Fourier [Sec.3h]
           coi = fourier_factor/sqrt(2.0)                    # Cone-of-influence [Sec.3g]
           dofmin = 2                                        # Degrees of freedom
         elseif mother == "PAUL"  #-----------------------------------  Paul
           if param == -1; param = 4. ; end
           m = param
           expnt = -(scale.*k).* heavi.(k)
           norm = sqrt(scale*k[2])*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n)
           daughter = norm*((scale.*k).^m).*exp.(expnt)
           daughter = daughter.* heavi.(k)        # Heaviside step function
           fourier_factor = 4*pi/(2.0*m+1)
           coi = fourier_factor*sqrt(2.0)
           dofmin = 2
         elseif mother == "DOG"  #-----------------------------------  DOG
           if param == -1; param = 2.; end
           m = param
           expnt = -(scale.*k).^2 ./ 2.0
           norm = sqrt(scale*k[2]/gamma(m+0.5))*sqrt(n)
           daughter = -norm*(im^m)*((scale.*k).^m).*exp.(expnt)
           fourier_factor = 2*pi*sqrt(2 ./ (2*m+1))
           coi = fourier_factor/sqrt(2)
           dofmin = 1
         else
           error("Mother must be one of MORLET,PAUL,DOG")
         end
         daughter, fourier_factor, coi, dofmin
end
    
    
# **********************************************************************************
# **********************************************************************************

#   1D contunous Wavelet transform
#
#   Computes the wavelet transform of the vector Y (length N),
#   with sampling rate DT.
#
#   By default, the Morlet wavelet (k0=6) is used.
#   The wavelet basis is normalized to have total energy=1 at all scales.
#
#
# INPUTS:
#
#    y = the time series of length N.
#    dt = amount of time between each Y value, i.e. the sampling time.
#
# OUTPUTS:
#
#    wave is the WAVELET transform of y. This is a complex array
#    of dimensions (n,j1+1). real(wave) gives the WAVELET amplitude,
#    atan(imag(wave),real(wave))  OR angle(wave) gives the WAVELET phase.
#    The WAVELET power spectrum is abs2(wave).
#    Its units are sigma^2 (the time series variance).
#
#
# OPTIONAL INPUTS:
#
# *** Note *** setting any of the following to -1 will cause the default
#               value to be used.
#
#    pad = if set to 1 (default is 0), pad time series with enough zeroes to get
#         N up to the next higher power of 2. This prevents wraparound
#         from the end of the time series to the beginning, and also
#         speeds up the FFT's used to do the wavelet transform.
#         This will not eliminate all edge effects (see COI below).
#
#    dj = the spacing between discrete scales. Default is 0.25.
#         A smaller # will give better scale resolution, but be slower to plot.
#
#    s0 = the smallest scale of the wavelet.  Default is 2*DT.
#
#    j1 = the # of scales minus one. Scales range from s0 up to s0*2^(j1*dj),
#        to give a total of (j1+1) scales. Default is j1 = (log2(n dt/s0))/dj.
#
#    mother = the mother wavelet function.
#             The choices are 'MORLET', 'PAUL', or 'DOG'
#
#    param = the mother wavelet parameter.
#            For 'MORLET' this is k0 (wavenumber), default is 6.
#            For 'PAUL' this is m (order), default is 4.
#            For 'DOG' this is m (m-th derivative), default is 2.
#
#
# OUTPUTS:
#
#    period = the vector of "Fourier" periods (in time units) that corresponds
#           to the scaels.
#
#    scale = the vector of scale indices, given by s0*2^(j*dj), j=0...j1
#            where j1+1 is the total # of scales.
#
#    coi = if specified, then return the Cone-of-Influence, which is a vector
#        of N points that contains the maximum period of useful information
#        at that particular time.
#        Periods greater than this are subject to edge effects.
#        This can be used to plot coi lines on a contour plot by doing:
#
#             IN MATLAB:
#              contour(time,log(period),log(power))
#              plot(time,log(coi),'k')
#
#----------------------------------------------------------------------------
#   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
#
#   This software may be used, copied, or redistributed as long as it is not
#   sold and this copyright notice is reproduced on each copy made. This
#   routine is provided as is without any express or implied warranties
#   whatsoever.
#
# Notice: Please acknowledge the use of the above software in any publications:
#    ``Wavelet software was provided by C. Torrence and G. Compo,
#      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
#
# Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
#            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
#
# Please send a copy of such publications to either C. Torrence or G. Compo:
#  Dr. Christopher Torrence               Dr. Gilbert P. Compo
#  Research Systems, Inc.                 Climate Diagnostics Center
#  4990 Pearl East Circle                 325 Broadway R/CDC1
#  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
#  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
#----------------------------------------------------------------------------

function wavelet(y,dt;pad=0,dj=-1,s0=-1,j1=-1,mother=-1,param=-1)

    n1 = length(y)
    if s0 == -1; s0=2*dt; end
    if dj == -1; dj = 1. / 4.; end
    if j1 == -1; j1=trunc(Int,(log(n1*dt/s0)/log(2))/dj); end
    if mother == -1; mother = "MORLET"; end

    #....construct time series to analyze, pad if necessary
    x = y .- mean(y)
    if pad == 1
        base2 = trunc(Int,log(n1)/log(2) + 0.4999)   # power of 2 nearest to N
        x = [x; zeros(2^(base2+1)-n1)]
    end
    n = length(x)

    #....construct wavenumber array used in transform [Eqn(5)]
    k = 1:trunc(Int,n/2)  
    k = k .* ((2. * pi)/(n*dt))
    k = [0.; k; -k[trunc(Int,(n-1)/2):-1:1]]
    #
    #....compute FFT of the (padded) time series
    f = fft(x)    # [Eqn(3)]

    #....construct SCALE array & empty PERIOD & WAVE arrays
    scale = s0*2 .^ ((0:j1)*dj)
    period = scale
    wave = zeros(Complex,j1+1,n)   # define the wavelet array

    # loop through all scales and compute transform
    fourier_factor = 0.0  # julia seems to need a default
    for a1 = 1:j1+1
        daughter,fourier_factor,coi,dofmin = wave_bases(mother,k,scale[a1],param)
        wave[a1,:] = ifft(f.*daughter)  # wavelet transform[Eqn(4)]
    end

    period = fourier_factor*scale
    #coi = coi*dt*[1e-5; collect(1:((n1+1)/2-1)); collect((n1/2-1):-1:1); 1e-5]  # COI [Sec.3g]
    coi = coi*dt*[1e-5; (1:((n1+1)/2-1)); ((n1/2-1):-1:1); 1e-5]  # COI [Sec.3g]
    wave = wave[:,1:n1]  # get rid of padding before returning

    wave, period, scale, coi
end

# continuoous wavelet transform (cwt)
function cwt(x,dt; 
    pad=1,        # pad the time series with zeros
    nvoices = 32, # total number of voices
    noctave = 4,  # total number of octaves
    s0 = 2*dt,    # starting scale
    mother = "MORLET",
    param = 6,
    dj = 1 ./ nvoices )
    # x : vector of time series
    # dt: sampling time interval
    # this function assumes that the time series is sampled equally
    # no checks are made to ensure this!!

    #dj = 1 ./ nvoices
    j1 = noctave * nvoices

    # normalize time series, but store the variance
    variance = var(x)
    x = x .- mean(x)    
    x = x ./ std(x)

    wave, period, scale, coi = wavelet(x,dt,pad,dj,s0,j1,mother,param)

    # wavelet power spetrum
    power = variance * abs2.(wave)   

    # Global wavelet spectrum 
    global_ws = variance * mean(power,dims=2)   # time-average over all times

    return wave, period, power, global_ws, variance, scale, coi
end

################################################
################################################

