#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
correl2displacement.py
==========================

A very simple script that convert the correlation measurement between two
images to a physical displacement and remove the median term from an image.

The script also copy georeferencing infos

Usage:
  tio_correl2displacement.py [--direction=<arg>] <img1_file> <correl_file> <output_file>

Options:
  --direction=<arg>  Direction of displacement (x or y) [default: x] 
  -h --help          Show this screen

"""

import sys
import numpy as np
from osgeo import gdal
import docopt

# Read command line & extract arguments
arguments = docopt.docopt(__doc__)

# Open inputs
img1_ds = gdal.OpenEx(arguments["<img1_file>"], gdal.OF_VERBOSE_ERROR)
band = img1_ds.GetRasterBand(1)
data = img1_ds.ReadAsArray(0, 0,
                          band.XSize, band.YSize)
if img1_ds == None:
    exit(1)
correl_ds = gdal.OpenEx(arguments["<correl_file>"], gdal.OF_VERBOSE_ERROR)
if correl_ds == None:
    exit(1)

# Get georeferncing from img1, this will be the reference
gt = img1_ds.GetGeoTransform()
proj = img1_ds.GetProjection()

# Create output file
drv = gdal.GetDriverByName("GTiff")
dst_ds = drv.Create(arguments["<output_file>"],
                    correl_ds.RasterXSize, correl_ds.RasterYSize, 1,
                    gdal.GDT_Float32)

# Set georeferencing from input
dst_ds.SetGeoTransform(gt)
dst_ds.SetProjection(proj)

median = np.nanmedian(data)

# Now, compute displacement and write it to output (done line by line
# to save memory)
correl_band = correl_ds.GetRasterBand(1)
dst_band = dst_ds.GetRasterBand(1)
resol = gt[1] if arguments["--direction"] == "x" else gt[5]
for y in range(correl_ds.RasterYSize):
    line = correl_band.ReadAsArray(0, y, correl_ds.RasterXSize, 1).astype(float)
    #line = (line - median) * 0.1 * resol
    line = (line - median) * resol
    dst_band.WriteArray(line, 0, y)

# Close output dataset
del dst_ds
