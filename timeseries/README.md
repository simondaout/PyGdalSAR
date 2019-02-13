To use this package
=============
* cube file: time series displacement maps in a binary format of size N\*Ncol\*Nlines, where N is the number of acquisitions dates, Ncol is the number of columns per date and Nlines, the number of lines per dates. To create a cube, use create\_cube.py
* list\_images file: text file of list of images made of 5 columns containing for each images 1) number 2) date in YYYYMMDD format 3) what you want (not use) 4) numerical date 5) perpendicular baseline.
To convert date in YYYYMMDD format to numerical format, use date2dec.py.
* lect.in: text file containing the number of columns and lines (eg. 1420 4000). (Depricated: it has to be replace by cube.hdr) 


invers\_disp2coef.py
============
Spatial and temporal inversions of the time series delay maps based on an iteration procedure. Requiered at least a cube of images (BIP format) and a list of images with temporal and perpendicular baselines. 
At each iteration, (1) estimation of spatial ramps, (2) linear decomposition in time based on a library of temporal functions (linear, heaviside, logarithm, seasonal...) or a discretised vector file.  
(3) estimation of RMS of each maps that will be then used as weight for the next iteration. Possibility to also to correct for a term proportional to the topography.

```
invers_disp2coef.py -h | --help
```

correct\_ts\_from\_gacos.py
============
Correct InSAR Time Series data from Gacos atmospheric models (data to be download and cited on: ceg-research.ncl.ac.uk/v2/gacos/). 1) Convert .ztd files to .tif format, 2) crop, re-project and re-resample atmospheric models to data geometry 3) correct time series data.

```
correct_ts_from_gacos.py -h | --help
```


invers\_disp\_pixel.py
============
Temporal decomposition of the time series delay maps of selected pixels. 

```
Usage: invers_disp_pixel.py  -h | --help 
```

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

clean\_ts.py
============
Clean a time series file (cube in binary format) given an other real4 file (mask) and a threshold on this mask

```
Usage: clean_ts.py -h | --help
```

geocode\_cube.py
============
Geocode cube of cumulative deplacements: create date.unw, geo\_date.unw, and geo\_date.tiff for each dates. !!! Need geocode.pl from ROIi\_PAC

```
Usage: geocode_cube.py -h | --help
```



