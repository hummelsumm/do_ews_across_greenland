#this is in python 3
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import wavepal_py3.Wavepal as wv
from scipy import optimize
import os.path
import os
import warnings
warnings.simplefilter('ignore', np.RankWarning)


# cpath = "/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/ice_core_data/in_use/NGRIP_d18O_and_dust_5cm.xlsx"
cpath = "ice_core_data/NGRIP_d18O_and_dust_5cm.xlsx"
cdat = pd.read_excel(cpath, sheet_name=2)

cdelta = cdat["Delta O18 (permil)"]
cdelta = cdelta.to_numpy()
ctime = cdat["GICC05 age (yr b2k)"]
ctime = ctime.to_numpy()
#print(cdelta.shape)

# ndat = np.genfromtxt('/Users/chu017/Library/CloudStorage/OneDrive-UiTOffice365/Documents/do_stuff/niklas_code/ngrip_age_depth_dust_d18o_unc_10k_60k.csv', delimiter = ';', filling_values = np.nan)
ndat = np.genfromtxt('ice_core_data/ngrip_age_depth_dust_d18o_unc_10k_60k.csv', delimiter = ';', filling_values = np.nan)
ntime = ndat[:,1]
ndelta = ndat[:,2]

## ages are rounded to first decimal in N's data!!!!!!!!
##Checkin that:
# for i,t in enumerate(ntime):
#     if t != ctime[ctime>=ntime[0]][i]:
#         ct = np.round(ctime[ctime>=ntime[0]][i], decimals = 1)
#         if t != ct:
#             ct +=0.1
#             if t != np.round(ct, decimals = 1):
#                 print(t, ct, ctime[ctime>=ntime[0]][i])


gs_end = [11705, 14690, 23375, 27790, 28910, 32520, 33735, 35505, 38235, 40165, 41480, 43365, 46860, 49315, 54235, 55815, 58280]
gs_start = [12900, 23105, 27460, 28510, 32025, 33390, 34725, 36590, 39935, 40790, 42100, 44285, 49105, 51660, 54745, 56440, 58515]

cts = []
nts = []
cxs = []
nxs = []

#tsp = plt.plot(ntime, ndelta, alpha = 0.5, c= "grey")
for st,en in zip(gs_end, gs_start):
    cidx = (st <= ctime) & (ctime <= en)
    nidx = (st <= ntime) & (ntime <= en)
    ctt = ctime[cidx]
    cxx = cdelta[cidx]
    ntt = ntime[nidx]
    nxx = ndelta[nidx]
    #print(type(xx))
    cts.append(ctt)
    cxs.append(cxx)
    nts.append(ntt)
    nxs.append(nxx)
    #plt.plot(tt,xx, c="darkred")
#plt.show()




def tau_estimation(y, t):
    ''' Estimates the  temporal decay scale of an (un)evenly spaced time series.

    Esimtates the temporal decay scale of an (un)evenly spaced time series. 
    Uses `scipy.optimize.minimize_scalar <https://docs.scipy.org/doc/scipy/reference/optimize.minimize_scalar-bounded.html>`_.

    Parameters
    ----------

    y : array
        A time series
    t : array
        Time axis of the time series

    Returns
    -------

    tau_est : float
        The estimated persistence

    References
    ----------

    Mudelsee, M. TAUEST: A Computer Program for Estimating Persistence in Unevenly Spaced Weather/Climate Time Series.
        Comput. Geosci. 28, 69–72 (2002).

    Code taken from
    ----------
    
    https://pyleoclim-util.readthedocs.io/en/latest/_modules/pyleoclim/utils/tsmodel.html#tau_estimation

    '''
    dt = np.diff(t)

    def ar1_fun(a):
        return np.sum((y[1:] - y[:-1]*a**dt)**2)
    a_est = optimize.minimize_scalar(ar1_fun, bounds=[0, 1], method='bounded').x
    #tau_est = -1 / np.log(a_est)

    return a_est#, tau_est

def runmean(t,x,w=200.0, only_full = True):
    t_start=t[0]
    t_end=t[-1]
    m=w/2.0
    n=x.size
    if only_full:
        mov_av_x=np.ones(n)*np.nan
    else:
        mov_av_x=np.zeros(n)
    ind_left=0
    ind_right=0

    for k in range(n):
        t_left=np.max([t_start,t[k]-m])
        ind_left=ind_left+np.argmin(np.absolute(t[ind_left:]-t_left))
        t_right=np.min([t_end,t[k]+m])
        ind_right=ind_right+np.argmin(np.absolute(t[ind_right:]-t_right))
        if only_full:
            if (t[k]-m >=t_start) & (t[k]+m <= t_end):
                mov_av_x[k] = np.mean(x[ind_left:ind_right+1])
        else:
            mov_av_x[k] = np.mean(x[ind_left:ind_right+1])
    return mov_av_x

def runstd(t,x,w=200.0, only_full = True):
    t_start=t[0]
    t_end=t[-1]
    m=w/2.0
    n=x.size
    if only_full:
        mov_av_x=np.ones(n)*np.nan
    else:
        mov_av_x=np.zeros(n)
    ind_left=0
    ind_right=0

    for k in range(n):
        t_left=np.max([t_start,t[k]-m])
        ind_left=ind_left+np.argmin(np.absolute(t[ind_left:]-t_left))
        t_right=np.min([t_end,t[k]+m])
        ind_right=ind_right+np.argmin(np.absolute(t[ind_right:]-t_right))
        if only_full:
            if (t[k]-m >=t_start) & (t[k]+m <= t_end):
                mov_av_x[k] = np.std(x[ind_left:ind_right+1])
        else:
            mov_av_x[k] = np.std(x[ind_left:ind_right+1])
    return mov_av_x

def runac(t,x,w=200.0, only_full = True):
    t_start=t[0]
    t_end=t[-1]
    m=w/2.0
    n=x.size
    if only_full:
        mov_av_x=np.ones(n)*np.nan
    else:
        mov_av_x=np.zeros(n)
    ind_left=0
    ind_right=0

    for k in range(n):
        t_left=np.max([t_start,t[k]-m])
        ind_left=ind_left+np.argmin(np.absolute(t[ind_left:]-t_left))
        t_right=np.min([t_end,t[k]+m])
        ind_right=ind_right+np.argmin(np.absolute(t[ind_right:]-t_right))
        if only_full:
            if (t[k]-m >=t_start) & (t[k]+m <= t_end):
                d  = np.mean(np.diff(t[ind_left:ind_right+1]))
                a_est= tau_estimation(x[ind_left:ind_right+1] - np.mean(x[ind_left:ind_right+1]), t[ind_left:ind_right+1])
                mov_av_x[k] = a_est**d         
        else:
            d  = np.mean(np.diff(t[ind_left:ind_right+1]))
            a_est= tau_estimation(x[ind_left:ind_right+1] - np.mean(x[ind_left:ind_right+1]), t[ind_left:ind_right+1])
            mov_av_x[k] = a_est**d
    return mov_av_x



def get_indicators(tt,xx, trend_deg = 0,scale_bounds = [(10,50)], 
                   theta_type = ["time", "regular", "default"][1], 
                   theta_dt = 5, w=200.0, filt = 100.0):
    x=wv.Wavepal(tt, xx, "Age", "$\\delta^{18}O$", t_units="yr (b2k)", mydata_units="permil")
    x.check_data()
    x.choose_trend_degree(trend_deg)
    x.trend_vectors()
    x.carma_params(signif_level_type="")

    w0=6.0
    scale_to_per = (4*np.pi)/(w0 + np.sqrt(2+w0**2))
    
    theta_type = theta_type.lower()
    if theta_type == "time":
        theta = x.t
    elif theta_type == "regular":
        theta_dt = int(theta_dt)
        if np.round(x.t)[0]%theta_dt==0:
            lb5 = np.round(x.t)[0]
        else:
            lb5 = np.round(x.t)[0] + theta_dt-np.round(x.t)[0]%theta_dt
        if np.round(x.t)[-1]%theta_dt==0:
            ub5 = np.round(x.t)[-1]
        else:
            ub5 = np.round(x.t)[-1] - np.round(x.t)[-1]%theta_dt
        theta = np.arange(lb5,ub5, step = theta_dt, dtype = float)
    else:
        theta = None

    pmin = float(scale_to_per*5)
    pmax = float(scale_to_per*300) #float(scale_to_per2*300)
    dj =0.1# 0.05# 0.1

    x.timefreq_analysis(theta=theta,
                        w0=w0,
                        permin=pmin,
                        permax=pmax, 
                        deltaj = dj, 
                        smoothing_coeff=0.0, #??
                        weighted_CWT = False, #??
                        computes_amplitude=True,
                        shannonnyquistexclusionzone=True) 
    scale = x.period_cwt /scale_to_per
    C_delta = 0.776
    scavs = np.zeros((len(scale_bounds), len(x.theta)))
    hursts = np.zeros((len(scale_bounds), len(x.theta)))
    
    for j,(s1,s2) in enumerate(scale_bounds):
        pow_in_ranges = np.ones_like(x.scalogram)*np.nan
        
        for i, sca in enumerate(scale):
            if (sca>= s1) & (sca <= s2):
                outside_coi = (x.theta > x.coi1[i]) & (x.theta < x.coi2[i])
                scalogr = x.scalogram[:,i][outside_coi]
                pow_in_ranges[:,i][outside_coi] = scalogr /sca
                scavs[j,outside_coi] += scalogr /sca
                scavs[j,~outside_coi] += np.nan
        avg = (scale>= s1) & (scale <= s2)
        if len(scale[avg])>0:
            pa = np.polyfit(np.log(scale[avg]), np.log(pow_in_ranges[:,avg]).T, 1)
            hursts[j,:] = (pa[0] + 1) / 2.
            ### to exclude times where at least 1 point lies in COI
            hursts[j,np.isnan(scavs[j,:])] = np.nan
        else:
            hursts[j,:] = np.ones(len(x.theta))*np.nan
   

    scavs *= (np.mean(np.diff(x.theta)) * dj)/C_delta
    
    l2 = np.min([filt/scale_to_per, x.period_cwt[-1]])
    x.timefreq_band_filtering([(float(x.period_ampl[0]),float(l2))])
    filt_dat = x.timefreq_band_filtered_signal[:,0]

    std = runstd(x.theta, filt_dat, w=w)
    ac = runac(x.theta, filt_dat,w=w)
    return x.theta, scavs, hursts, std, ac, filt_dat, x
    

    
def get_tfts(x, eps = 0.05, ns=100):
    #first shuffling, then keeping low freqs, as done in timeseriessurrogates.jl 
    tfts = np.zeros((ns,len(x.mydata)))
    n = len(x.freq)
    n_preserve = int(np.ceil(n*eps))
    
    for js in range(ns):
        freq_shuffled = np.random.permutation(x.freq)
        freq_shuffled[:n_preserve] = x.freq[:n_preserve]
        for (i,f) in enumerate(freq_shuffled):
            tfts[js,:] += x.amplitude_cos[i]*np.cos(2*np.pi*f*x.t) + x.amplitude_sin[i]*np.sin(2*np.pi*f*x.t)
    return tfts

def get_shuffled_phases(x,ns=100):
    shuffled = np.zeros((ns,len(x.mydata)))
    
    for js in range(ns):
        freq_shuffled = np.random.permutation(x.freq)
        for (i,f) in enumerate(freq_shuffled):
            shuffled[js,:] += x.amplitude_cos[i]*np.cos(2*np.pi*f*x.t) + x.amplitude_sin[i]*np.sin(2*np.pi*f*x.t)
    return shuffled

#### make proper arrays / data structures and save ####

n_surrs = 100
s1s = [10,20,30,40,50,60,70,80,90,100]
s2s = [20,30,40,50,60,70,80,90,100,110]
range_step = int(np.mean(np.diff(s1s)))
scale_bounds = []
for i in s1s:
    for j in s2s:
        if i<j:
            scale_bounds.append((i,j))

data_type = "c" #"n"

if data_type =="c":
    ts = cts
    xs = cxs
elif data_type == "n":
    ts = nts
    xs = nxs

method = "TFTS"#"FOURIER"
theta_type = "time" #"regular" "default"


savedir = "wavepal_test/indicators_and_surrs/"



def save_slopes_and_surrs(n_surrs, s1s,s2s,data_type, 
                          method = ["FOURIER", "TFTS"][0], 
                          trend_deg = 0,
                          theta_type = ["time", "regular", "default"][0], theta_dt = 5, 
                          w=200.0, filt = 100.0, 
                          savedir = "wavepal_test/indicators_and_surrs/"):
    if len(s1s)>1:
        range_step = int(np.mean(np.diff(s1s)))
    else:
        range_step = 0

    scale_bounds = []
    for i in s1s:
        for j in s2s:
            if i<j:
                #nscales +=1
                scale_bounds.append((i,j))

    if data_type =="c":
        ts = cts
        xs = cxs
    elif data_type == "n":
        ts = nts
        xs = nxs
    
    method = method.lower()

    saveto_temp =  f"TEMP_indicators_NGRIP_{theta_dt}y_{data_type}_wavepal_theta_{theta_type}_sranges_{s1s[0]}:{range_step}:{s2s[-1]}_{method.lower()}_ns_{n_surrs}.npz"
    saveto = f"indicators_NGRIP_{theta_dt}y_{data_type}_wavepal_theta_{theta_type}_sranges_{s1s[0]}:{range_step}:{s2s[-1]}_{method.lower()}_ns_{n_surrs}.npz"

    if not os.path.exists(savedir + saveto):
        print("################################################")
        print("################################################")
        print("################################################")
        print("calculating", saveto)
        print("################################################")
        print("################################################")
        print("################################################")
        print("")

        if not os.path.exists(savedir + saveto_temp):
            slopes_w = np.zeros((len(xs), n_surrs+1, len(scale_bounds)))
            slopes_h = np.zeros((len(xs), n_surrs+1, len(scale_bounds)))
            slopes_sigma = np.zeros((len(xs), n_surrs+1))
            slopes_alpha = np.zeros((len(xs), n_surrs+1))


            thetas = [[] for i in range(len(xs))]
            scavs = [[] for i in range(len(xs))]
            hursts = [[] for i in range(len(xs))]
            stds = [[] for i in range(len(xs))]
            acs = [[] for i in range(len(xs))]
            filts = [[] for i in range(len(xs))]
            gs_counter = -1
            surr_counter = -1
        
        else:
            temp_file =  np.load(savedir + saveto_temp)
            
            gs_counter = temp_file["gs_counter"]
            surr_counter = temp_file["surr_counter"]
            #print("surr_counter", surr_counter)

            slopes_w = temp_file["slopes_w"]
            slopes_h = temp_file["slopes_h"]
            slopes_sigma = temp_file["slopes_sigma"]
            slopes_alpha = temp_file["slopes_alpha"]

            thetas = list(temp_file["thetas"])
            scavs = list(temp_file["scavs"])
            hursts = list(temp_file["hursts"])
            stds = list(temp_file["stds"])
            acs = list(temp_file["acs"])
            filts = list(temp_file["filts"])

            for v in range(len(thetas)):
                thetas[v] = list(thetas[v])
                scavs[v] = list(scavs[v])
                hursts[v] = list(hursts[v])
                stds[v] = list(stds[v])
                acs[v] = list(acs[v])
                filts[v] = list(filts[v])

            thetas.extend([[]]*(16-gs_counter))
            scavs.extend([[]]*(16-gs_counter))
            hursts.extend([[]]*(16-gs_counter))
            stds.extend([[]]*(16-gs_counter))
            acs.extend([[]]*(16-gs_counter))
            filts.extend([[]]*(16-gs_counter))    

        for i,(tt,xx) in enumerate(zip(ts, xs)):
            print("GS", i)
            if (i>gs_counter) or ((i==gs_counter) and (surr_counter < n_surrs -1)):# or (surr_counter < n_surrs-1):
                
                thet, scav, hurst, std, ac, filtdat, yr= get_indicators(tt,xx, 
                                                trend_deg = trend_deg,scale_bounds = scale_bounds, 
                                                theta_type = theta_type, 
                                                theta_dt = theta_dt, w=w, filt = filt)
                if i>gs_counter:
                    thetas[i].extend(thet)
                    scavs[i].append(scav) #scav has len len(scale_bounds)
                    hursts[i].append(hurst) #hurst has len len(scale_bounds)
                    stds[i].append(std)
                    acs[i].append(ac)
                    filts[i].extend(filtdat)
                    surr_counter = -1
                
                D = yr.t[-1]-yr.t[0]
                yr.freq_analysis(computes_amplitude=True, D=float(D))
                if method == "tfts":
                    surs = get_tfts(yr,ns=n_surrs)
                else:
                    surs = get_shuffled_phases(yr, ns=n_surrs)
                
                for j,s in enumerate(surs):
                    if j%10 ==0:
                        print("GS", i, "surrogate", j)
                    if j>surr_counter:

                        _,ww,h,sigma,alpha,_,_ = get_indicators(tt,s, 
                                                    trend_deg = trend_deg,scale_bounds = scale_bounds, 
                                                    theta_type = theta_type, 
                                                    theta_dt = theta_dt, w=w, filt = filt)
                        scavs[i].append(ww)
                        hursts[i].append(h)
                        stds[i].append(sigma)
                        acs[i].append(alpha)

                        

                        if  (j+1)%10 == 0:
                            #saving every 10th surrogate
                            print("saving at GS", i, "surrogate", j)
                                      
                            start_to_change = max(0,j-9)
                            for ks in range(start_to_change,j+1):
                                sigma = stds[i][ks][:len(thet)]
                                alpha = acs[i][ks][:len(thet)]
                                #print(i,k, slopes_sigma.shape)
                                if len(sigma[~np.isnan(sigma)])>0:
                                    slopes_sigma[i,ks] = -np.polyfit(thet[~np.isnan(sigma)],sigma[~np.isnan(sigma)],1)[0]
                                else:
                                    slopes_sigma[i,ks] = np.nan
                                if len(alpha[~np.isnan(alpha)])>0:
                                    slopes_alpha[i,ks] = -np.polyfit(thet[~np.isnan(alpha)],acs[i][ks][~np.isnan(alpha)],1)[0]
                                else:
                                    slopes_alpha[i,ks] = np.nan

                                for l, www in enumerate(scavs[i][ks]): #the single bands
                                    ww = www[:len(thet)]
                                    h = hursts[i][ks][l][:len(thet)]
                                    if len(ww[~np.isnan(ww)])>0:
                                        slopes_w[i,ks,l] = -np.polyfit(thet[~np.isnan(ww)],ww[~np.isnan(ww)],1)[0]
                                    else:
                                        slopes_w[i,ks,l] = np.nan
                                    if len(h[~np.isnan(h)])>0:
                                        slopes_h[i,ks,l] = -np.polyfit(thet[~np.isnan(h)],h[~np.isnan(h)],1)[0]
                                    else:
                                        slopes_h[i,ks,l] = np.nan
                                
                            print("transforming lists to matrices to save")

                            lens_thetas = [len(l) for l in thetas[:i+1]]
                            max_len_thetas = max(lens_thetas)
                            thetas2 = np.empty((len(thetas[:i+1]), max_len_thetas))*np.nan
                            filts2 = np.empty((len(filts[:i+1]), max_len_thetas))*np.nan
                            mask = np.arange(max_len_thetas) < np.array(lens_thetas)[:,None]
                            thetas2[mask] = np.concatenate(thetas[:i+1])
                            filts2[mask] = np.concatenate(filts[:i+1])


                            #std and ac
                            #array size: 17 with n_surrs +1 subarrays with lengths theta subarrays
                            if i>0:
                                dim2 = n_surrs + 1
                            else:
                                dim2 = j+2
                            stds2 = np.empty((len(stds[:i+1]), dim2, max_len_thetas))*np.nan
                            acs2 = np.empty((len(acs[:i+1]), dim2, max_len_thetas))*np.nan
                            mask2 = np.tile(mask[:,None], (1,dim2,1))
                            #######
                            mask2[i,j+2:,:] = False
                            stds2[mask2] = np.concatenate([np.concatenate(l) for l in stds[:i+1]])
                            acs2[mask2] = np.concatenate([np.concatenate(l) for l in acs[:i+1]])

                            #hurst and scas
                            #array size: 17 with n_surrs+1 subarrays with n_scales subarrays with lens theta subarrays
                            hursts2 = np.empty((len(hursts[:i+1]), dim2, len(scale_bounds) ,max_len_thetas))*np.nan
                            scavs2 = np.empty((len(scavs[:i+1]),dim2, len(scale_bounds) ,max_len_thetas))*np.nan
                            mask3 = np.tile(mask[:,None, None],(1,dim2,len(scale_bounds),1))
                            ########
                            mask3[i,j+2:,:,:] = False
                            hursts2[mask3] = np.concatenate([np.concatenate([np.concatenate(l) for l in k]) for k in hursts[:i+1]])
                            scavs2[mask3] = np.concatenate([np.concatenate([np.concatenate(l) for l in k]) for k in scavs[:i+1]])

                            surr_counter +=10
                            gs_counter = i
                            
                            print(f"saving indicators and slopes at gs {gs_counter} surrogate {surr_counter} to", saveto_temp)

                            np.savez(savedir + saveto_temp, 
                                thetas = thetas2, 
                                scavs = scavs2, 
                                hursts = hursts2,
                                stds = stds2,
                                acs = acs2,
                                filts = filts2,
                                slopes_w = slopes_w,
                                slopes_h = slopes_h,
                                slopes_sigma = slopes_sigma,
                                slopes_alpha = slopes_alpha,
                                scale_bounds = scale_bounds,
                                gs_counter = gs_counter,
                                surr_counter = surr_counter
                                )
                print("")    
            else:
                print(saveto_temp, "already exists :)")
                print("moving on to next GS")

        print("all done, changing temporary file to", saveto)
        os.rename(savedir + saveto_temp, savedir + saveto)
 
    else:
        print(saveto, "already exists :)")


#choose n_surrs to be divisible by 10 !!


#just "the good one"
save_slopes_and_surrs(1000, [10],[50],"c", 
                          method = ["FOURIER", "TFTS"][1], 
                          trend_deg = 0,
                          theta_type = ["time", "regular", "default"][0], theta_dt = 5, 
                          w=200.0, filt = 100.0, 
                        #   savedir = "wavepal_test/indicators_and_surrs/"
                            # savedir = "indicators_and_surrs/"
                            savedir = "new_surrogate_files/NGRIP_irreg/"
                          )
