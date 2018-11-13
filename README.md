# PyGdalSAR
Post-processing InSAR tool package written in the Python and Gdal programming language.
Predictives and empirical atmospheric corrections for unwrapped interferograms. Iterative spatial and temporal decompositions of Geodetic Time Series (InSAR, GPS or Pixel Offsets). Additional tools for time series correction from the GACOS atmospheric models (ceg-research.ncl.ac.uk/v2/gacos/), plotting, cleaning, radar to geographic conversion... 


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

