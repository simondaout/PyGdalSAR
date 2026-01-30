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
r42geo.py
-------------
Convert r4 file to geotif.

Usage: r42geo.py --infile=<path> --ref=<path> [--lectfile=<path>] 
geo2r4.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       r4 file 
--ref PATH          ref tiff file
--lectfile= PATH     Path of the lect.in file for r4 format
"""

from osgeo import gdal
import numpy as np 
import docopt
import os, sys
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
if arguments["--lectfile"] ==  None:
    lecfile = "lect.in"
else:
    lecfile = arguments["--lectfile"]

# open infile
ncols, nlines = map(int, open(lecfile).readline().split(None, 2)[0:2])
phi = np.fromfile(infile,dtype=np.float32)[:nlines*ncols].reshape((nlines,ncols))
print("> Driver:   REAL4  band file")
print("> Size:     ", ncols,'x',nlines,'x')
print("> Datatype: FLOAT32")

# open ref
geotiff = arguments["--ref"]
georef = gdal.Open(geotiff)
gt = georef.GetGeoTransform()
proj = georef.GetProjection()
driver = gdal.GetDriverByName('GTiff')

# save output
outfile = os.path.splitext(infile)[0] + '.tiff'
dst_ds = driver.Create(outfile, ncols, nlines, 1, gdal.GDT_Float32)
dst_band2 = dst_ds.GetRasterBand(1)
dst_band2.WriteArray(phi,0,0)
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)
dst_band2.FlushCache()
del dst_ds

