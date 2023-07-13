correct\_ts\_from\_gacos.py
============
Correct InSAR Time Series data from Gacos atmospheric models (data to be download and cited on: ceg-research.ncl.ac.uk/v2/gacos/). 1) Convert .ztd files to .tif format, 2) crop, re-project and re-resample atmospheric models to data geometry 3) correct time series data.

```
Usage: correct_ts_from_gacos.py -h | --help
```

Reference:
* [Dini, B., Daout, S., Manconi, A., & Loew, S. (2019). Classification of slope processes based on multitemporal DInSAR analyses in the Himalaya of NW Bhutan. Remote Sensing of Environment, 233, 111408.] (https://www.sciencedirect.com/science/article/abs/pii/S0034425719304274)

invert\_ramp\_topo\_unw.py
============
Removes atmospheric phase/elevation correlations or/and azimuthal and range ramps polynomial coefficeints on unwrapped interferograms (Compatible with NSBAS/ GAMMA/ GTiff formats). Possibility to include cross-dependent polynomial functions of range, azimuth and topography, use masks based on additional file or amplide of the RMG interferogram. Reconstruction of the empirical phase correction by time series inversion and impose time series closure of all coefficients.

```
Usage: invert_ramp_topo_unw.py -h | --help 
```

The program has a lot of options that are detailed with the help command. For instance, you might want to correct your list of interferograms (interf_pair_success.txt) from a quadratic ramp in range, a linear ramp in azimuth (flat=5), from a quadratic relationship (nfit=1) with the topography and refer all your interferograms to zero between lines 3000 and 3400 and colunmns 0 and 1000 (ref_zone). In order to eliminate outliers for the estimations, you also might want to mask pixels with LOS above the 95% percentile (perc=95), pixels with slope (DEM gradient) bellow their 90% percentile (perc_slope=90) and with coherence bellow 0.2 (cohpixel=0.2). To do so, run:

```
Exemple: invert_ramp_topo_unw.py --ref_zone=3000,3400,0,1000 --int_list=../interf_pair_success.txt --prefix =filt_col_ --suffix=_sd --rlook=4 --cohpixel=yes --threshold_coh=0.2 --flat=5 --perc=95 -- perc_slope=90 --topofile=../20170222/radar_4rlks.hgt --nfit=1 --tsinv=yes
```

The program produces several output figures.
– The first one corresponds to the phase-elevation plot (Fig. 20141016-20141109phase-topo.png) 
– The second figure correspond is a comparison between the LOS data, the model and the residual LOS in map view
– The third figure is very similar to the second one with the difference that the model is the one obtained after time series inversion of all independent parameters of the empirical function (time series inversion performed if tsinv=yes). This step is very important to impose the consistency of the correction within the network. Stronger weights are given for short baselines interferograms that are less likely to contains deformation (Fig. 20141016-20141109corrections.png)


Reference:
* [Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017692)


gacos2rdr.py / atmos_corr_histo.py / correct_ifg_from_gacos.py  
============
Tools to correct interferograms from ERA5 or GACOS atmospheric models

Reference:
* Dodds, N., Daout, S., Walker, R. T., Begenjev, G., Bezmenov, Y., Mirzin, R., & Parsons, B. (2022). Interseismic deformation and strain-partitioning along the Main Köpetdag Fault, Turkmenistan, with Sentinel-1 InSAR time-series. Geophysical Journal International, 230(3), 1612-1629.


fetch\_t2m.py / fetch\_era5.py / extract\_param\_from\_datanc.py / plot_histo_datanc.py 
============
Fetch any ERAI parameters from the ECMWF server defined in a input text file. For more information about the API format, visit: https://confluence.ecmwf.int/display/CKB/Global+data%3A+Download+data+from+ECMWF+for+a+particular+area+and+resolution 

Fetch any ERA5 parameters from the ECMWF server (hard coding)

Open and plot parameter in fonction of the time from a \*.nc file containing multiple bands for sevral time steps or differents pressure levels doanload for ECMWF server. 

Reference:
* [Daout, Simon, et al. "Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau." Geophysical Research Letters 44.2 (2017): 901-909](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2016GL070781) 


