# this is python2 (v 2.7)
import numpy as np
from scipy.interpolate import UnivariateSpline
import pandas as pd
from scipy.signal import filtfilt, cheby1
import os

def cheby_lowpass(cutoff, fs, order, rp):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = cheby1(order, rp, normal_cutoff, btype='low', analog=False)
    return b, a

def cheby_lowpass_filter(x, cutoff, fs, order, rp):
    b, a = cheby_lowpass(cutoff, fs, order, rp)
    y = filtfilt(b, a, x)
    return y

def cheby_highpass(cutoff, fs, order, rp):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = cheby1(order, rp, normal_cutoff, btype='high', analog=False)
    return b, a

def cheby_highpass_filter(x, cutoff, fs, order, rp):
    b, a = cheby_highpass(cutoff, fs, order, rp)
    y = filtfilt(b, a, x)
    return y

# set parameters for chebychef filter:
order = 8
# the maximum ripple allowed below unity gain in the passband (in dB)
rp = .05
# sampling frequency, should be equal to 1 because we have 1-yr steps
fss = 1.
# the sampling width to which we want to interpolate the data
samp = 5
# MTM time half-bandwidth
NW = 4

sec_factor = .95
sw = 100

# the corresponding cutoff frequency
cutoff = .5 * 1. / samp
# offset at beginning and end:
offset = 19


#cpath = "/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/ice_core_data/in_use/NGRIP_d18O_and_dust_5cm.xlsx"
cpath = "ice_core_data/NGRIP_d18O_and_dust_5cm.xlsx"
cdat = pd.read_excel(cpath, sheet_name=2)

cdelta = cdat["Delta O18 (permil)"]
cdelta = cdelta.to_numpy()
ctime = cdat["GICC05 age (yr b2k)"]
ctime = ctime.to_numpy()
#print(cdelta.shape)

#ndat = np.genfromtxt('/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/niklas_code/ngrip_age_depth_dust_d18o_unc_10k_60k.csv', delimiter = ';', filling_values = np.nan)
ndat = np.genfromtxt('ice_core_data/ngrip_age_depth_dust_d18o_unc_10k_60k.csv', delimiter = ';', filling_values = np.nan)
ntime = ndat[:,1]
ndelta = ndat[:,2]

#print("len(ntime), len(ctime)",len(ntime), len(ctime))

new_ctime = ctime[len(ctime)- len(ntime):]
new_cdelta = cdelta[len(cdelta)- len(ndelta):]

#print("len(ntime), len(new_ctime)", len(ntime), len(new_ctime))
#print("np.all(ndelta == new_cdelta)", np.all(ndelta == new_cdelta))

#check that the difference is really just rounding! 
for i,k in enumerate(new_ctime):
    if np.round(k,1)!= ntime[i]:
        #to make sure that we round up
        k2 = k + 0.01
        if np.round(k2,1)!= ntime[i]:
            print(k, k2, np.round(k2,1), ntime[i])


#savedir = '/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/new_surrogate_files/ice_cores_py/'
savedir = 'new_surrogate_files/ice_cores_py/'
names = ["N_py2","C_py2"]

for age,ngrip,saveto in zip([ntime,new_ctime],[ndelta,new_cdelta],names):
    age_int = np.array(age, dtype = 'int')
    end_age = 10277
    start = np.argmin(np.abs(age - end_age))
    end = -1
    age = age[start : end]
    age_int = age_int[start : end]
    age_ipd = np.arange(age_int.min(), age_int.max())

    time = np.arange(age.min(), age.max(), samp)[::-1]
    d18o = ngrip[start : end]
    
    #print(len(age), len(d18o))

    f_d18o = UnivariateSpline(age, d18o, s = 0.)
    d18o_ipd_cub_temp = f_d18o(age_ipd)
    d18o_ipd_cub = cheby_lowpass_filter(d18o_ipd_cub_temp, sec_factor * cutoff, fss, order, rp)
    #print(offset)
    dat = d18o_ipd_cub[offset : -offset][::samp][::-1]
    
    dat_no_lowpass = d18o_ipd_cub_temp[offset : -offset][::samp][::-1]
    
    time = age_ipd[offset : -offset][::samp][::-1]

    #print(time.shape)
    #print(dat.shape)

    mdat = np.mean(dat)
    sdat = np.std(dat, ddof = 1)
    dat_normed = (dat - mdat) / sdat
    
    #print(np.mean(dat_no_lowpass))
    #print(np.std(dat_no_lowpass, ddof = 1))
    dat_normed_no_lowpass = (dat_no_lowpass - np.mean(dat_no_lowpass)) / np.std(dat_no_lowpass, ddof = 1)

    filt = 100
    scavf = 200 / samp

    lpf = 800

    filt_dat = cheby_highpass_filter(dat_normed, .95 * 1. / filt, 1. / samp, 8, .05)
    
    filt_dat_no_lowpass = cheby_highpass_filter(dat_normed_no_lowpass, .95 * 1. / filt, 1. / samp, 8, .05)
        
    print(savedir + saveto)
    

    np.savez(savedir + saveto,
        time = time,
        dat = dat,
        dat_normed = dat_normed,
        filt_dat = filt_dat,
        dat_no_lowpass = dat_no_lowpass,
        dat_normed_no_lowpass = dat_normed_no_lowpass,
        filt_dat_no_lowpass = filt_dat_no_lowpass)