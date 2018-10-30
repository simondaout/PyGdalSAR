# ReversibleDeformation: example GACOS corrections, empirical phase/topo corrections and extraction of reversible derformations with seasonal functions.

* Convert DEM to ASPECT file in the LOS geometry given average Lattitute and LOS angle
```
dem2slope.py --infile=DEM_UTM.tif --outfile=DEM_slope_UTM.tif --los=34.5 --heading=-14 --lat=27.5
```

* Correct Data from GACOS atmospheric models: convert .ztd to tif, crop, resample, re-project and correct data
```
correct_ts_from_gacos.py --crop=714887.000,3049145.000,719887.000,3056645.000  --proj=32645
```

* Temporal decomposition of two selected pixels into linear and seasonal fonctions 
```
invers_disp_pixel.py --cube=depl_cumule --list_images=list_images.txt --lectfile=lect_ts.in --interseismic=yes --seasonal=yes --windowsize=2 --cols=156,100 --ligns=204,80 --name=TS --rad2mm=1 --dem=no --plot=yes
``

* Temporal decomposition of all acquisitions corrected from GACOS
* First perform a ramp and empirical corrections on data masked from mask_landslide_utm.tif and negative aspects and then perform the decomposition into linear and seasonal fonctions
```
invers_disp2coef.py --cube=depl_cumule_gacos --list_images=list_images.txt --interseismic=yes --seasonal=yes --flat=3 --topofile=DEM_UTM.tif --niter=2 --mask=mask_landslide_utm.tif --threshold_mask=-3 --scale_mask=-1 --nfit=1 --geotiff=DEM_UTM.tif --aspect=DEM_slope_UTM.tif --plot=yes
``` 
