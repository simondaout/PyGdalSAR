# ReversibleDeformation: example GACOS corrections, empirical phase/topo corrections and extraction of reversible derformations with seasonal functions.

* Convert DEM to ASPECT file in the LOS geometry given average Lattitute and LOS angle
```
dem2slope.py --help
dem2slope.py --infile=DEM_UTM.tif --outfile=DEM_slope_UTM.tif --los=34.5 --heading=-14 --lat=27.5
```

* Plot Tiff file
```
plot_tiff.py --help
plot_tiff.py --geotiff=DEM_UTM.tif 
```

* Correct Data from GACOS atmospheric models: convert .ztd to tif, crop (xmin,xmax,ymin,ymax in UTM output geometry), resample (take width and lenght data), re-project (UTM, EPSG:32645) and correct data
```
correct_ts_from_gacos.py --help
correct_ts_from_gacos.py --cube=depl_cumule --crop=714887.000,3049145.000,719887.000,3056645.000  --proj=32645 --gacos2data=-824.126
```

* Temporal decomposition of two selected pixels into linear and seasonal fonctions 
```
invers_disp_pixel.py --help
invers_disp_pixel.py --cube=depl_cumule --list_images=list_images.txt --lectfile=lect_ts.in --interseismic=yes --seasonal=yes --windowsize=2 --cols=156,100 --ligns=204,80 --name=TS --rad2mm=1 --dem=no --plot=yes
```

* Temporal decomposition of all acquisitions corrected from GACOS. First perform a ramp (function of latitude and longitude) and empirical corrections (quadratic pahse/elevation relationship) on data masked from mask\_landslide\_utm.tif (absolute value of cumulative displacements) and negative aspects, then, perform the decomposition of all pixels into linear and seasonal fonctions
```
invers_disp2coef.py --help
invers_disp2coef.py --cube=depl_cumule_gacos --list_images=list_images.txt --interseismic=yes --seasonal=yes --flat=3 --topofile=DEM_UTM.tif --niter=2 --mask=mask_landslide_utm.tif --threshold_mask=-3 --scale_mask=-1 --nfit=1 --geotiff=DEM_UTM.tif --aspect=DEM_slope_UTM.tif --plot=yes
``` 

* display time series output for selected pixels with uncertainties equal to the RMS of the first iteration (aps\_0.txt)
```
lect_disp_pixel.py --cols=156,100 --ligns=204,80 --list_images=list_images.txt --cos=coswt_coeff.tif --sin=sinwt_coeff.tif --slope=lin_coeff.tif --aps=aps_0.txt --rad2mm=1
```

* Visualize amplitude of the seasonal deformation in map view
```
plot_tiff.py --geotiff=ampwt_coeff.tif --vmax=10 --vmin=0
```
