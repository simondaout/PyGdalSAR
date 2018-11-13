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
geo2r4.py
-------------

Usage: geo2r4.py --infile=<path> 
geo2r4.py -h | --help

Options:
-h --help           Show this screen
--infile PATH       tiff file 

"""

import gdal
import numpy as np 
import docopt
import os, sys
from matplotlib import pyplot as plt
import matplotlib.cm as cm

# read arguments
arguments = docopt.docopt(__doc__)
infile = arguments["--infile"]
ds = gdal.Open(infile, gdal.GA_ReadOnly)

ds_extension = os.path.splitext(infile)[1]
ds_name = os.path.splitext(infile)[0]

if ds_extension == ".unw":
	band = ds.GetRasterBand(2)
	m = band.ReadAsArray(0, 0,
           ds.RasterXSize, ds.RasterYSize,
           ds.RasterXSize, ds.RasterYSize)
if ds_extension == ".tif":
	band = ds.GetRasterBand(1)
	m = band.ReadAsArray()

# save output
outfile = ds_name + '.r4'
fid = open(outfile,'wb')
m.flatten().astype('float32').tofile(fid)

# initiate figure depl
fig = plt.figure(1,figsize=(6,5))
# plot topo
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(m,cmap=cm.jet)
ax.set_title(outfile,fontsize=6)
fig.colorbar(cax, orientation='vertical',aspect=10)

plt.show()



