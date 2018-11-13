Tools for correction from the GACOS atmospheric models (ceg-research.ncl.ac.uk/v2/gacos/) or empirical correction on wappred interferograms.

correct\_ts\_from\_gacos.py
============
Correct InSAR Time Series data from Gacos atmospheric models (data to be download and cited on: ceg-research.ncl.ac.uk/v2/gacos/). 1) Convert .ztd files to .tif format, 2) crop, re-project and re-resample atmospheric models to data geometry 3) correct time series data.

```
Usage: correct_ts_from_gacos.py -h | --help
```

invert\_ramp\_topo\_unw.py
============
Removes atmospheric phase/elevation correlations or/and azimuthal and range ramps polynomial coefficeints on unwrapped interferograms (2 bands BIL format). Reconstruction of the empirical phase correction by time series inversion.

```
Usage: invert_ramp_topo_unw.py
```

fetch\_t2m.py
============
Fetch any ERAI parameters from the ECMWF server defined in a input text file. For more information about the API format, visit: https://confluence.ecmwf.int/display/CKB/Global+data%3A+Download+data+from+ECMWF+for+a+particular+area+and+resolution 

```
Usage:
  fetch_erai.py [-v] <proc_file>
  fetch_erai.py -h | --help
```

fetch\_era5.py
============
Fetch any ERA5 parameters from the ECMWF server (hard coding)

extract\_param\_from\_datanc.py
============
Open and plot parameter in fonction of the time from a \*.nc file containing multiple bands for sevral time steps or differents pressure levels doanload for ECMWF server. 

```
Usage:
  extract\_param\_from\_datanc.py.py [-v] <data.nc>
  extract\_param\_from\_datanc.py -h | --help
```



