# do_ews_across_greenland
Code for the analysis of EWS for DO events across Greenland ice core records related to

Hummel, C., Boers, N., and Rypdal, M. (2025): Inconclusive early warning signals for Dansgaard-Oeschger events across Greenland ice cores, [preprint], https://doi.org/10.5194/egusphere-2024-3567

### $\delta^{18}\text{O}$ records from Greenland ice cores
- NGRIP with irregular, 5-, 10- and 20-year temporal resolution 
- NEEM with 10-year temporal resolution
- GRIP with 20-year temporal resolution
- GISP2 with 20-year temporal resolution
  
### Early warning indicators
- **based on critical slowing down**: Variance $V$, lag-one autocorrelation parameter $\alpha_1$, and restoring rate $\lambda$
  - analyses of $V$ and $\alpha_1$ of regularly sampled records are based on and adapted from Niklas Boers' code related to
    Boers, N. (2018): Early-warning signals for Dansgaard-Oeschger events in a high-resolution ice core record, https://doi.org/10.1038/s41467-018-04881-7
  - analyses of $\lambda$ of NGRIP in 5-year resolution are based on Niklas Boers' code related to
    Boers, N. (2021): Observation-based early-warning signals for a collapse of the Atlantic Meridional Overturning Circulation, https://doi.org/10.1038/s41558-021-01097-4
  - analyses of $\alpha_1$ of the irregularly sampled NGRIP record are based on
    Mudelsee, M. (2002): TAUEST: a computer program for estimating persistence in unevenly spaced weather/climate time series, https://doi.org/10.1016/S0098-3004(01)00041-3,
    as well as Lenoir, G. and Crucifix, M. (2018): A general theory on frequency and time–frequency analysis of irregularly sampled time series based on projection methods – Part 2: Extension to time–frequency analysis, https://doi.org/10.5194/npg-25-175-2018, 2018.
    The corresponding python2 code is provided in the WAVEPAL package (https://github.com/guillaumelenoir/WAVEPAL), which has been adapted to python3.
  
  
- **based on wavelet analysis**: scale-averaged wavelet coefficient $\hat{w}^2$, and local Hurst exponent $\hat{H}^{\text{loc}}$
  - analyses of regularly sampled records are based on
    Torrence, C. and Compo, G. P. (1998): A Practical Guide to Wavelet Analysis, https://doi.org/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2.
    Wavelet software was provided by C. Torrence and G. Compo, and is available at URL: http://paos.colorado.edu/research/wavelets/''.
    We use code adapted from Bernd Blasius' implementation of his WaveletAnalysis code (https://github.com/berndblasius/WaveletAnalysis)
  - analyses of irregularly sampled records are based on
    Lenoir, G. and Crucifix, M. (2018): A general theory on frequency and time–frequency analysis of irregularly sampled time series based on projection methods – Part 2: Extension to time–frequency analysis, https://doi.org/10.5194/npg-25-175-2018.
    The corresponding python2 code is provided in the WAVEPAL package (https://github.com/guillaumelenoir/WAVEPAL), which has been adapted to python3.


### Data sources:
- NGRIP:
  NGRIP members (2004): High-resolution record of Northern Hemisphere climate extending into the last interglacial period, https://doi.org/10.1038/nature02805, and
Gkinis et al. (2014): Water isotope diffusion rates from the North-GRIP ice core for the last 16,000 years – Glaciological and paleoclimatic implications, https://doi.org/10.1016/j.epsl.2014.08.022,
provided on https://www.iceandclimate.nbi.ku.dk/data/
- NEEM:
  Gkinis, V., et al. (2020): NEEM ice core High Resolution (0.05 m) Water Isotope Ratios (18O/16O, 2H/1H) covering 8-129 ky b2k, https://doi.org/10.1594/PANGAEA.925552
- NGRIP, GRIP, GISP2 (20-year resampling):
  Seierstad et al. (2014): Consistently dated records from the Greenland GRIP, GISP2 and NGRIP ice cores for the past 104 ka reveal regional millennial-scale δ18O gradients with possible Heinrich event imprint, https://doi.org/10.1016/j.quascirev.2014.10.032, 
  and
Rasmussen et al. (2014): A stratigraphic framework for abrupt climatic changes during the Last Glacial period based on three synchronized Greenland ice-core records: refining and extending the INTIMATE event stratigraphy, https://doi.org/10.1016/j.quascirev.2014.09.007, 
  provided on https://www.iceandclimate.nbi.ku.dk/data/
- GS/GI onsets: as given in Table S1 in the supplement of
  Boers, N. (2018): Early-warning signals for Dansgaard-Oeschger events in a high-resolution ice core record, https://doi.org/10.1038/s41467-018-04881-7
