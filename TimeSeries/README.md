To use this package
=============
* cube file: time series displacement maps in a binary format of size N\*Ncol\*Nlines, where N is the number of acquisitions dates, Ncol is the number of columns per date and Nlines, the number of lines per dates. To create a cube, use create\_cube.py
* list\_images file: text file of list of images made of 5 columns containing for each images 1) number 2) date in YYYYMMDD format 3) what you want (not use) 4) numerical date 5) perpendicular baseline.
To convert date in YYYYMMDD format to numerical format, use date2dec.py.
* lect.in: text file containing the number of columns and lines (eg. 1420 4000). (Depricated: it has to be replace by cube.hdr)
* aps: RMSpixel.txt file from invers_pixel output: 3 columns txt file: 1) nb, 2) date, 3) RMS 

invers\_disp2coef.py
============
Spatial and temporal inversions of the time series delay maps based on an iteration procedure. Requiered at least a cube of images (BIP format) and a list of images with temporal and perpendicular baselines. 
At each iteration, (1) estimation of spatial ramps, (2) linear decomposition in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal...) or a discretised vector file.  
(3) estimation of RMS of each maps that will be then used as weight for the next iteration. Possibility to also to correct for a term proportional to the topography.

```
invers_disp2coef.py -h | --help
```

For example, to decompose all your time series pixels into a linear and seasonal function with a simple linear ramp inversion in range and azimuth on your cumulative maps prior to the decomposition, run:

```
invers_disp2coef.py --seasonal=yes --flat=3
```

Possibility, the program can also correct for a term proportional to the topography in the first spatial step given a topographic file as input (topofile) and defining the order of the polynomial fit (nfit, ivar). Remember to use the crop (ibeg, iend, jbeg, jend) or mask (mask) options to mask deforming areas for the spatial estimations.

```
invers_disp2coef.py --seasonal=yes --flat=3 --topofile=DEM.r4 --mask=mask.tif --threshold_mask=-3 --scale_mask=-1 --nfit=1 --ivar=0
```

invers\_disp\_pixel.py
============
Temporal decomposition of the time series delay maps of selected pixels. 

```
Usage: invers_disp_pixel.py  -h | --help 
invers_disp_pixel.py --cube=depl_cumule_clean --cols=917,704,1076 --ligns=1108,1169,1131 --aps= inaps.txt
```

This program also offers the possibility to decompose the time series into various sources. For example if you want to decompose your time series into a linear and a seasonal function with an additional DEM error term proportional to the perpendicular baseline, run:

```
invers_disp_pixel.py --cube=depl_cumule_clean --cols=917,704,1076 --ligns=1108,1169,1131 --aps= inaps.txt --seasonal=yes --dem=yes
```

![Alt text](exemple_ts_plot.png)

lect\_disp\_pixel.py
=============

Plot time series results obtained with invers\_disp2coef.py for given pixels.

```
Usage: lect_cube_pixel.py   -h | --help
```

invers\_disp\_gps.py
============

Temporal decomposition of gps time series and projection into a LOS vector. 

```
Usage: invers_disp_gps.py  -h | --help 
``` 

invers\_pca.py / invers\_ica.py / reconstruct\_ica.py
============

PCA/ICA time series decomposition 

Reference:
* [Maubant, L., Pathier, E., Daout, S., Radiguet, M., Doin, M. P., Kazachkina, E., ... & Walpersdorf, A. (2020). Independent component analysis and parametric approach for source separation in InSAR time series at regional scale: application to the 2017–2018 Slow Slip Event in Guerrero (Mexico). Journal of Geophysical Research: Solid Earth, 125(3), e2019JB018187.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB018187)


clean\_phi-amp.py / compute\_modseas.py
============

Seasonal analysis of the output file of invers\_disp2coef.py with seasonal=yes as done in:

* [Daout, S., Dini, B., Haeberli, W., Doin, M. P., & Parsons, B. (2020). Ice loss in the Northeastern Tibetan Plateau permafrost as seen by 16 yr of ESA SAR missions. Earth and Planetary Science Letters, 545, 116404.] (https://www.sciencedirect.com/science/article/abs/pii/S0012821X20303484)

For instance to clean all maps with low sesaonl amplitude displacements anc compute statistics from it:

```
clean_phi-amp.py --slopefile=Slope.r4 --threshold_amp=.8 --plotcorr=yes
```

where Slope.r4 is the DEM slope file. 

You can then compute the mod of seasonal displacements in a crop area with:

```
compute_modseas.py --crop=500,1200,1000,2000 --minamp=3.5 --maxamp=40 --name=valley_flank 
```

![Alt text](clean_phi-amp_era_histo.png)

clean\_ts.py
============
Clean a time series file (cube in binary format) given an other real4 file (mask) and a threshold on this mask

```
Usage: clean_ts.py -h | --help
```


plot\_ts.py / plot\_annimate\_ts.py
============

Plot Time series or create a .mp4 animation of your time series cube. 
 

geocode\_cube.py
============
Geocode cube of cumulative deplacements: create date.unw, geo\_date.unw, and geo\_date.tiff for each dates. !!! Need geocode.pl from ROI\_PAC

```
Usage: geocode_cube.py -h | --help
```

- to plot the time series, use the command line

```
plot_geots.py --vmax=30 --vmin=-30 --geocrop=37,38,95.5,97
```

where vmin, vmax defined the min and max of the colorscale and geocrop defines the crop boundaries. You can also create an .mp4 animation of your time series with the command line:

```
plot_annimate_geots.py --cube=depl_cumule
```

![Alt text](geots_zom.png)

date2dec.py
============
convert dates in YYYYMMDD format to numerical format


 References
============

* [Daout, Simon, et al. "Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau." Geophysical Research Letters 44.2 (2017): 901-909](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2016GL070781)

* [Daout, S., Sudhaus, H., Kausch, T., Steinberg, A., & Dini, B. (2019). Interseismic and postseismic shallow creep of the North Qaidam Thrust faults detected with a multitemporal InSAR analysis. Journal of Geophysical Research: Solid Earth, 124(7), 7259-7279.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JB017692)



