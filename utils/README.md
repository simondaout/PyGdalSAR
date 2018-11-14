add\_r4.py
============
Add, substract, compute the amplitude (```infile1**2,infile2**2```) ar the phase (```arctan(infile1/infile2)) from two real4 files.

```
Usage: add_r4.py -h | --help
```

add\_rmg.py
============
Add, substract, and plot RMG (BIL) unwrapped files.

```
Usage: add_rmg.py -h | --help
```

clean\_r4.py
============
Clean a real4 file given an other real4 file (mask) and a threshold on this mask

```
Usage: clean_r4.py -h | --help
```


cut\_r4.py
============
Crop r4 file

```
Usage: cut_r4.py -h | --help
```

cut\_rmg.py
============
Crop RMG (BIL) file

```
Usage: cut_rmg.py -h | --help
```

date2dec.py
============
convert dates in YYYYMMDD format to numerical format

```
Usage: date2dec.py -h | --help
```

dem2slope.py
============
Convert dem to slope file in the LOS geeometry knowing mean LOS heading and latitudes.

```
Usage: dem2slope.py -h | --help
```

extend\_r4.py
============
Extend real4 file to new number of lines and columns

```
Usage: extend_r4.py -h | --help
```

geo2r4.py
============
Convert geotiff or any raster made of 1 band into real4 format. 

```
Usage: geo2r4.py -h | --help
```

histo\_r4.py
============
Compute histogram distribution of a real4 file

```
Usage: histo_r4.py -h | --help
```

invert\_phi.py
============
Time series inversion of given parameters: solve phase/noise/spectrum for each acquisition dates using phase/noise/spectrum of all interferograms such as phi(AB) = phi(B) - phi(A) or sigma(AB)^2 = sigma(A)^2 + sigma(B)^2

```
Usage: invert_phi.py -h | --help
```

ll\_to\_radar.py
============
Transform a list of geographic coordinates into radar corrdinates

```
Usage: ll_to_radar.py -h | --help
```

mask\_unw.py
============
Mask unwrapped interferogram (BIL format) given a mask file in real4 format. 

```
Usage: mask_unw.py -h | --help
```

