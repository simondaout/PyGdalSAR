#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################
#
# PyGdalSAR: An InSAR post-processing package 
# written in Python-Gdal
#
############################################
# Author        : Simon DAOUT (Oxford)
############################################


"""\
radar_to_ll.py
-------------
transform a list of radar corrdinates into geographic coordinates

Usage: ll_to_radar.py --cols=<values> --ligns=<values> [--latlonfile=<path>] 

Options:
-h --help           Show this screen.
--cols VALUE        Pixel columss 
--ligns VALUE        Pixel lines
--latlonfile PATH   Path to lat_erai_4rlks.r4  [default: ./latlon_4rlks.trans]  
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided
from decimal import Decimal
from osgeo import gdal

import docopt
arguments = docopt.docopt(__doc__)
if arguments["--latlonfile"] ==  None:
   infile = 'latlon_4rlks.trans'
else:
   infile = arguments["--latlonfile"]

ds = gdal.Open(infile, gdal.GA_ReadOnly)
lat_band = ds.GetRasterBand(1)
lat = lat_band.ReadAsArray(0, 0,
         ds.RasterXSize, ds.RasterYSize,
         ds.RasterXSize, ds.RasterYSize)
lon_band = ds.GetRasterBand(2)
lon = lon_band.ReadAsArray(0, 0,
         ds.RasterXSize, ds.RasterYSize,
         ds.RasterXSize, ds.RasterYSize)
nlines, ncol = ds.RasterYSize, ds.RasterXSize

list_cols = list(map(int,arguments["--cols"].replace(',',' ').split()))
list_lines = list(map(int,arguments["--ligns"].replace(',',' ').split()))
if len(list_cols) != len(list_lines):
   raise Exception("cols and ligns lists are not the same size")

index = np.array(list_lines),np.array(list_cols)

list_lats = lat[index]
list_lons = lon[index]

print('lats:',list_lats)
print('lons:',list_lons)
