#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
################################################################################
# 
# NSBAS - New Small Baseline Chain
# 
################################################################################
# Author        : Simon DAOUT (ISTerre)
################################################################################

"""\
plot_tiff.py
-------------

Usage: plot_tiff.py --geotiff=<path> [--vmin=<value>] [--vmax=<value>] [--wrap=<yes/no>] [--crop=<values>]

Options:
-h --help           Show this screen.
--geotiff PATH      Tiff file to be displayed
--vmax 				Max colorscale [default: 90th percentile]
--vmin 				Min colorscale [default: 10th percentile]
--wrap 				Wrapped phase [default: no]
--crop VALUE        Define a region of interest [default: 0,nlign,0,ncol]
"""

# numpy
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from pylab import *

import gdal

# docopt (command line parser)
from nsbas import docopt

# read arguments
arguments = docopt.docopt(__doc__)
geotiff = arguments["--geotiff"]

ds = gdal.Open(geotiff)
nlign,ncol = ds.RasterYSize, ds.RasterXSize

if arguments["--crop"] ==  None:
    crop = [0,nlign,0,ncol]
else:
    crop = map(float,arguments["--crop"].replace(',',' ').split())
ibeg,iend,jbeg,jend = int(crop[0]),int(crop[1]),int(crop[2]),int(crop[3])

band = ds.GetRasterBand(1)
m = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize)
m = m[ibeg:iend,jbeg:jend]
kk = np.nonzero(np.logical_or(np.logical_or(~np.isnan(m), np.abs(m)<999.),m==0.0))
mprim = m[kk]

if arguments["--vmax"] ==  None:
	vmax = np.nanpercentile(mprim, 90)
else:
	vmax = np.float(arguments["--vmax"])

if arguments["--vmin"] ==  None:
	vmin = np.nanpercentile(mprim, 10)
else:
	vmin = np.float(arguments["--vmin"])

if arguments["--wrap"] ==  'yes':	
	m = np.mod(m,2*np.pi)-np.pi 
	vmax=np.pi
	vmin=-np.pi

# Plot
fig = plt.figure(0,figsize=(6,4))
ax = fig.add_subplot(1,1,1)
cax = ax.imshow(m,cmap=cm.jet,vmax=vmax,vmin=vmin)
ax.set_title(geotiff)
setp( ax.get_xticklabels(), visible=False)
cbar = fig.colorbar(cax, orientation='vertical',aspect=5)

fig.savefig('{}.eps'.format(geotiff), format='EPS',dpi=150)
plt.show()

