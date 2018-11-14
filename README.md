# PyGdalSAR
Post-processing InSAR tool package written in the Python and Gdal programming language. 
It can be utilizes for predictives and empirical atmospheric corrections on the unwrapped interferograms, or time series correction from the GACOS atmospheric models (ceg-research.ncl.ac.uk/v2/gacos/). Additional package for iterative spatial and temporal decompositions of geodetic Time Series (InSAR, GPS or Pixel Offsets), plotting, cleaning of interferograms or time series, radar to geographic conversions, etc... 


To install the package
=============
```git clone https://github.com/simondaout/TSdecomp.git```

To update the software: 
```git pull```

In case of fire:
```
git commit
git push
leave building
```

Organisation
=============
This project contains the following folders:
 * timeseries: programms for time series decomposition (see documantation bellow)
 * atmocorr: programms for atmospehric corrections.
 * plots: plotting programms for .tiff, .r4, time series or .unw files
 * utils: additional tools for cleaning .r4 or .unw formats, georeferenced formats (default: .tif), geocoding, conversion dates to decimal dates.
 * src: fortrtan programms, in progress
 * tutorial: 
 	 * ReversibleDeformation: example GACOS corrections, empirical phase/topo corrections and extraction of reversible derformations with seasonal functions.
 	 * Vector: in progress

Requirements
=============
This project needs the following external components:
 * Python (2.7)
 * py27-gdal
 * py27-numPy
 * py27-sciPy
 * Gdal (https://www.gdal.org)
 * docopt (Available on https://github.com/docopt/docopt or in ./python/).
 * py27-matplotlib
 * py27-datetime
 * NSBAS (optional): http://efidir.poleterresolide.fr/index.php/effidir-tools/nsbas

Developpers & Contact
=============
```
Simon Daout: simon.daout@earth.ox.ac.uk
Postdoctoral Research Fellow, Department of Earth Sciences, University of Oxford
```

```
Collaborators:  
Louise Maubant: louise.maubant@univ-grenoble-alpes.fr
Benedetta Dini: benedetta.dini@erdw.ethz.ch
Marie-Pierre Doin: marie-pierre.doin@univ-grenoble-alpes.fr
Mathieu Volat: matthieu.volat@univ-lyon1.fr
```
 References
============

* [Daout, Simon, et al. "Strain Partitioning and Present‐Day Fault Kinematics in NW Tibet From Envisat SAR Interferometry." Journal of Geophysical Research: Solid Earth 123.3 (2018): 2462-2483](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/2017JB015020)

* [Daout, Simon, et al. "Large‐scale InSAR monitoring of permafrost freeze‐thaw cycles on the Tibetan Plateau." Geophysical Research Letters 44.2 (2017): 901-909](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2016GL070781)
