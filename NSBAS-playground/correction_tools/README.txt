nsb_unw_ERAcorr_stats.py / nsb_unw_GACOScorr_stats.py
============
Tools to correct interferograms from ERA5 or GACOS atmospheric models

Reference:
* Dodds, N., Daout, S., Walker, R. T., Begenjev, G., Bezmenov, Y., Mirzin, R., & Parsons, B. (2022). Interseismic deformation and strain-partitioning along the Main Köpetdag Fault, Turkmenistan, with Sentinel-1 InSAR time-series. Geophysical Journal International, 230(3), 1612-1629.


nsb_filtflatunw.py
============
Equivalent to nsb_filtflatunw.pl in NSBAS

Reference:
* [Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017692)

Contrary to other processing chains, NSBAS provides various empirical corrections on the wrapped phase in order to reduce the variability of phase and thus enable unwrapping over large areas (Doin et al. (2015); Daout et al. (2017, 2018)). These corrections include (1) ramp corrections in range or azimuth to remove residual orbital errors, ionospheric delays, and effects of the satellite clock drift (Fig.), (2) empirical estimation of the phase delays due to differences in the stratified troposphere with a possible mask on deforming areas, (3) empirical removal of a fringe pattern (e.g co-seismic fringe pattern, seasonal ground deformation making the unwrapping difficult) (Fig. 23d). All those corrections are estimated by searching on small moving sub-windows the proportionality between (1) phase and range, (2) phase and elevation, or (3) phase and the template that maximises the local phase coherence within the subwindows (Doin et al., 2011; Doin et al., 2015), respectively.

![My Image](processing.png)

The available jobs are check_look, ecmwf, replace_amp, filterSW, filterROI, flatr, flat_az, flat_atmo, flat_model, colin, look_int, unwrapping, add_flatr_back, add_flata_back, add_model_back, add_atmo_back. Bellow, some details about each of this processing step.


make_fit_az.sh / correct_rgaz_unw.py / correct_ramp_unw.py
============
Old NSBAS Tools

