import numpy as np
import scipy.stats as st
import statsmodels.api as sm
#import matplotlib.pyplot as plt


def tfts_surrogates_shuffled(ts2,eps=0.05, ns=100):
    m = np.mean(ts2)
    ts = ts2 - m
    tfts = np.zeros((ns, ts.shape[0]))
    ts_fourier  = np.fft.rfft(ts) #forward
    n_preserve = int(np.ceil((ts.shape[0] // 2 + 1)*eps))
    #random_phases = np.exp(np.random.uniform(0, 2 * np.pi, (ns, ts.shape[0] // 2 + 1)) * 1.0j)
    #old_phases = np.exp(np.linspace(0, 2 * np.pi, ts.shape[0] // 2 + 1, endpoint=False) * 1.0j)
    old_phases = np.exp(np.linspace(0, 2 * np.pi, ts.shape[0] // 2 + 1) * 1.0j)
    for js in range(ns):
        shuffled_phases = np.random.permutation(old_phases)
        shuffled_phases[:n_preserve] = old_phases[:n_preserve]
        ts_fourier_new = ts_fourier * shuffled_phases
        tfts[js,:] = np.real(np.fft.irfft(ts_fourier_new,n= ts.shape[0])) + m
    return tfts

def run_fit_a_ar1(x, w, adj=False):
    n = x.shape[0]
    xs = np.zeros_like(x)
    
    for i in range(w // 2):
        xs[i] = np.nan
    
    for i in range(n - w // 2, n):
        xs[i] = np.nan
    
    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]
        xw = xw - xw.mean()
        
        #p0, p1 = np.polyfit(np.arange(xw.shape[0]), xw, 1)
        lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
        p0 = lg[0]
        p1 = lg[1]
        
        xw = xw - p0 * np.arange(xw.shape[0]) - p1
        

        # dxw = xw[2:] - xw[:-2] / 2
        dxw = xw[1:] - xw[:-1]
        
        xw = sm.add_constant(xw)
        model = sm.GLSAR(dxw, xw[:-1], rho=1)
        results = model.iterative_fit(maxiter=10)
        # results = model.iterative_fit()
        
        a = results.params[1]
        
        if adj:
            xs[i] = np.log(a + 1)
        else:
            xs[i] = a
    return xs


xx = np.load("ngrip_5.npz")
ns = 10_000
w=200//5

max_len = 1684
times = np.empty((17,max_len))*np.nan#[[] for i in range(17)]
lambdas = np.empty((17,max_len))*np.nan#[[] for i in range(17)]
slopes = np.empty(17)
p_ones = np.empty(17)
surr_slopes = np.empty((17,ns))


for gsi,gs in enumerate(np.arange(1,18)): 
    print("calculating GS", gs)
    x = xx[f"{gs}"][:,1]
    n = len(x)
    times[gsi,:n] = xx[f"{gs}"][:,0]
    
    l = run_fit_a_ar1(x,w, adj=False)
    lambdas[gsi,:n] = l
    
    lmask = ~np.isnan(l)
    tm = np.arange(1,l.shape[0]+1)[lmask]
    lm = l[lmask]
    slope = st.linregress(tm, lm)[0]
    slopes[gsi] = slope

    tfts = tfts_surrogates_shuffled(x, eps=0.05,ns=ns)
    
    #surr_slopes = np.zeros(ns)
    for i in range(ns):
        ls = run_fit_a_ar1(tfts[i,:],w, adj=False)
        lsmask = ~np.isnan(ls)
        tsm = np.arange(1,ls.shape[0]+1)[lsmask]
        lsm = ls[lsmask]
        surr_slopes[gsi,i] = st.linregress(tsm, lsm)[0]
    
    p_ones[gsi] = np.count_nonzero(surr_slopes >= slope)/ns

print("saving")
np.savez(f"lambda_ngrip5_w_{w}_{ns}_TFTS.npz",
         times = times,
         vals= lambdas,
         slopes = slopes,
         surr_slopes = surr_slopes,
         p_one = p_ones)

print("Done :)")