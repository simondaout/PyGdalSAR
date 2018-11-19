#!/usr/bin/env python2
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
mask_unw.py
-------------
Mask unwrapped interferogram (BIL format) given a mask file in real4 format.

Usage: mask_unw.py --infile=<path> --maskfile=<path> --outfile=<path> --threshold=<value> --plot=<yes/no>

Options:
-h --help           Show this screen.
--infile PATH       File to be masked
--maskfile PATH     File to mask int
--outfile PATH      Masked int
--threshold         Threshold for the mask (mask pixels < threshold)
--plot VALUE        Display interferograms
"""


# lib
from __future__ import print_function
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

import os, sys

import gdal
import time

# docopt (command line parser)
import docopt

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
maskfile = arguments["--maskfile"] 
outfile = arguments["--outfile"]
seuil = float(arguments["--threshold"])
plot = arguments["--plot"]

gdal.UseExceptions()
# Open dataset (image)
ds = gdal.Open(infile, gdal.GA_ReadOnly)
# Get the band that have the data we want
ds_band1 = ds.GetRasterBand(1)
ds_band2 = ds.GetRasterBand(2)
# Attributes
print("> Driver:   ", ds.GetDriver().ShortName)
print("> Size:     ", ds.RasterXSize,'x',ds.RasterYSize,'x',ds.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds_band2.DataType))
amp = ds_band1.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
data = ds_band2.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
nlign,ncol = ds.RasterYSize, ds.RasterXSize
del ds

# Open mask
ds2 = gdal.Open(maskfile, gdal.GA_ReadOnly)
ds2_band1 = ds2.GetRasterBand(1)
ds2_band2 = ds2.GetRasterBand(2)
print("> Driver:   ", ds2.GetDriver().ShortName)
print("> Size:     ", ds2.RasterXSize,'x',ds2.RasterYSize,'x',ds2.RasterCount)
print("> Datatype: ", gdal.GetDataTypeName(ds2_band2.DataType))
maskamp = ds2_band1.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
maskdata = ds2_band2.ReadAsArray(0, 0, ds2.RasterXSize, ds2.RasterYSize)
del ds2

# clean
outamp = np.copy(amp)
outdata = np.copy(data)
kk = np.nonzero(maskdata[:nlign,:ncol]>seuil)
outdata[kk], outamp[kk] = 0, 0

# create new GDAL image
drv = gdal.GetDriverByName("roi_pac")
dst_ds = drv.Create(outfile, ncol, nlign, 2, gdal.GDT_Float32)
dst_band1 = dst_ds.GetRasterBand(1)
dst_band2 = dst_ds.GetRasterBand(2)
dst_band1.WriteArray(outamp,0,0)
dst_band2.WriteArray(outdata,0,0)
del dst_ds

if plot=="yes":
      vmax = np.percentile(data, 98)
      # Plot
      fig = plt.figure(0,figsize=(8,9))
      ax = fig.add_subplot(1,2,1)
      cax = ax.imshow(data,cmap=cm.gist_rainbow,vmax=vmax)
      ax.set_title(infile)
      setp( ax.get_xticklabels(), visible=False)

      ax = fig.add_subplot(1,2,2)
      cax = ax.imshow(outdata,cmap=cm.gist_rainbow,vmax=vmax)
      ax.set_title(outfile)
      setp( ax.get_xticklabels(), visible=False)

      cbar = fig.colorbar(cax, orientation='vertical',aspect=5)

      # Display the data
      fig.canvas.set_window_title(sys.argv[1])
      plt.show()
